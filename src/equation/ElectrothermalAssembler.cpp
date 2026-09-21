#include "vela/equation/ElectrothermalAssembler.h"
#include "vela/discretization/ThermalSgCurrent.h"
#include "vela/equation/LatticeBandEdgeWork.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include <cmath>
#include <stdexcept>
#include <limits>

namespace vela {
ElectrothermalAssembler::ElectrothermalAssembler(const DeviceMesh& mesh,
    const DopingModel& doping, ElectrothermalGeometry geometry, LatticeHeatAssembler heat,
    MobilityModelConfig mobility, SiliconThermalPhysics physics, Real muE, Real muH,bool reusePreparation,bool reuseIalScreening,bool reuseNeutralRoots,bool safeguardedNeutralNewton)
    :mesh_(mesh),doping_(doping),geometry_(std::move(geometry)),heat_(std::move(heat)),
     mobility_(std::move(mobility)),physics_(std::move(physics)),muE_(muE),muH_(muH),reusePreparation_(reusePreparation),reuseIalScreening_(reuseIalScreening),reuseNeutralRoots_(reuseNeutralRoots),safeguardedNeutralNewton_(safeguardedNeutralNewton) {
    const auto n=mesh.numNodes(),e=mesh.numEdges();
    if(geometry_.recombinationArea_m2.size()==0)geometry_.recombinationArea_m2=geometry_.siliconArea_m2;
    if(doping.numNodes()!=n || geometry_.siliconArea_m2.size()!=n || geometry_.recombinationArea_m2.size()!=n ||
       geometry_.fixedCharge_C_per_m.size()!=n || geometry_.poissonEdge_F_per_m.size()!=e ||
       geometry_.transportWeight.size()!=e || heat_.nodalAreas_m2().size()!=n)
        throw std::invalid_argument("Electrothermal geometry size mismatch");
    for(const auto* v:{&geometry_.siliconArea_m2,&geometry_.recombinationArea_m2,&geometry_.fixedCharge_C_per_m,
                      &geometry_.poissonEdge_F_per_m,&geometry_.transportWeight})
        if(!v->allFinite())throw std::invalid_argument("Nonfinite electrothermal geometry");
    if((geometry_.siliconArea_m2.array()<0.).any() || (geometry_.recombinationArea_m2.array()<0.).any() || (geometry_.transportWeight.array()<0.).any() ||
       !std::isfinite(muE_) || !std::isfinite(muH_) || muE_<=0. || muH_<=0.)
        throw std::invalid_argument("Invalid electrothermal area/weight/mobility");
    if(mobility_.model!="constant" && mobility_.model!="ialmob")
        throw std::invalid_argument("Electrothermal operator supports explicit constant or IALMob mobility");
    if(mobility_.internalConcentrationToM3!=1. || mobility_.internalMobilityToM2PerVS!=1. ||
       mobility_.internalFieldToVPerM!=1.)
        throw std::invalid_argument("Electrothermal operator requires SI material units");
    for(const auto& edge:mesh.edges())
        if(geometry_.transportWeight[edge.id]>0. &&
           (geometry_.siliconArea_m2[edge.n0]<=0. || geometry_.siliconArea_m2[edge.n1]<=0.))
            throw std::invalid_argument("Transport edge touches a non-silicon node");
    prepareIalTransportGeometry(mobility_,mesh_);
    if(reusePreparation_){dopingPreparation_.resize(n);temperaturePreparation_.resize(n);}
    if(reuseNeutralRoots_)neutralRoots_.resize(n);
}

const SiliconThermalPhysics::TemperaturePreparation& ElectrothermalAssembler::preparedAt(Index node,Real temperature) const {
    auto& dop=dopingPreparation_.at(node);auto& state=temperaturePreparation_.at(node);
    const Real nd=doping_.donors(node),na=doping_.acceptors(node);
    if(!dop || !dop->matches(nd,na)){
        dop=physics_.prepareDoping(nd,na);state.reset();++preparationCounts_[0];
    }
    if(!state || !state->matches(temperature,nd,na)){
        state=physics_.prepareTemperature(temperature,*dop);++preparationCounts_[1];
    }else ++preparationCounts_[2];
    return *state;
}
std::pair<Real,Real> ElectrothermalAssembler::neutralPotential(Index node,Real bias,Real t) const {
    if(node>=mesh_.numNodes() || geometry_.siliconArea_m2[node]<=0. || !std::isfinite(bias))
        throw std::invalid_argument("Neutral boundary requires a silicon node and finite bias");
    SiliconThermalState s{bias,bias,bias,t,doping_.donors(node),doping_.acceptors(node)};
    if(reuseNeutralRoots_)if(const auto& cached=neutralRoots_[node];cached &&
       cached->bias==bias && std::signbit(cached->bias)==std::signbit(bias) &&
       cached->temperature==t && cached->donors==s.donors_m3 && cached->acceptors==s.acceptors_m3) {
        ++neutralRootCounts_[1];return cached->value;
    }
    ++neutralRootCounts_[0];
    const auto* prepared=reusePreparation_?&preparedAt(node,t):nullptr;
    const auto densities=[&](){
        ++neutralRootIterationCounts_[0];
        if(prepared)return physics_.carrierDensities(s,*prepared);
        const auto p=physics_.evaluate(s);return std::array<ThermalQuantity,2>{p.electrons_m3,p.holes_m3};
    };
    const Real net=s.donors_m3-s.acceptors_m3;
    Real lo=bias-1.5,hi=bias+1.5;
    s.potential_V=lo;const auto lower=densities();
    s.potential_V=hi;const auto upper=densities();
    if(lower[0].value-lower[1].value>net || upper[0].value-upper[1].value<net)
        throw std::invalid_argument("Neutral potential outside the audited silicon bracket");
    bool resolved=false;
    if(safeguardedNeutralNewton_) {
        Real next=(lo+hi)/2.;bool lastNeighbor=false;
        for(int k=0;k<70;++k) {
            s.potential_V=next;const auto p=densities();
            const Real charge=p[0].value-p[1].value;
            // Preserve the original strict comparison, including flat rounded roots.
            if(charge>net)hi=next;else lo=next;
            const Real middle=(lo+hi)/2.;
            if(middle==lo || middle==hi){resolved=true;break;}
            const Real slope=p[0].derivative[0]-p[1].derivative[0];
            const Real proposal=std::isfinite(slope) && slope>0.?next-(charge-net)/slope:
                std::numeric_limits<Real>::quiet_NaN();
            if(std::isfinite(slope) && slope>0. && std::isfinite(proposal) &&
               proposal>lo && proposal<hi && std::abs(proposal-next)<=.5*(hi-lo)) {
                next=proposal;lastNeighbor=false;++neutralRootIterationCounts_[1];
            } else if(proposal==next && !lastNeighbor) {
                // Probe the other side of a rounded Newton root once; repeated
                // adjacent-float walking across a flat interval is not useful.
                const Real neighbor=std::nextafter(next,charge>net?lo:hi);
                if(neighbor>lo && neighbor<hi) {
                    next=neighbor;lastNeighbor=true;++neutralRootIterationCounts_[3];
                } else {next=middle;lastNeighbor=false;++neutralRootIterationCounts_[2];}
            } else {next=middle;lastNeighbor=false;++neutralRootIterationCounts_[2];}
        }
        // Close to zero, 70 original bisections need not reach adjacent floats.
        // Preserve that finite-budget convention instead of silently returning
        // a different rounded contact target or temperature derivative.
        const Real root=(lo+hi)/2.;
        const Real spacing=std::abs(std::nextafter(root,std::numeric_limits<Real>::infinity())-root);
        if(spacing<std::ldexp((bias+1.5)-(bias-1.5),-70))resolved=false;
        if(!resolved){lo=bias-1.5;hi=bias+1.5;++neutralRootIterationCounts_[4];}
    }
    if(!resolved)for(int k=0;k<70;++k){s.potential_V=(lo+hi)/2.;const auto p=densities();
        ++neutralRootIterationCounts_[2];
        if(p[0].value-p[1].value>net)hi=s.potential_V;else lo=s.potential_V;}
    s.potential_V=(lo+hi)/2.;const auto p=densities();
    const Real slope=p[0].derivative[0]-p[1].derivative[0];
    const std::pair<Real,Real> result{s.potential_V,-(p[0].derivative[3]-p[1].derivative[3])/slope};
    if(reuseNeutralRoots_)neutralRoots_[node]=NeutralRoot{bias,t,s.donors_m3,s.acceptors_m3,result};
    return result;
}

ElectrothermalAssembly ElectrothermalAssembler::assemble(const VectorXd& x,
    const ElectrothermalBoundary& bc,const VectorXd& eReference,const VectorXd& hReference,bool buildJacobian,bool skipEquilibriumTransport,
    bool diagnosticFreezeMobilityDerivatives,bool diagnosticFreezeRecombinationDerivatives,
    std::vector<ElectrothermalHoleRowAudit>* holeRowAudit) const {
    const Index n=mesh_.numNodes();
    std::map<Index,ElectrothermalHoleRowAudit*> audited;
    if(holeRowAudit)for(auto& row:*holeRowAudit){
        if(row.node>=n || geometry_.siliconArea_m2[row.node]<=0. || bc.holeQf_V.contains(row.node) ||
           (bc.neutralContactBias_V.contains(row.node) && !bc.holeRecombination.contains(row.node)) || !audited.emplace(row.node,&row).second)
            throw std::invalid_argument("Hole row audit requires distinct active silicon hole rows");
        row.edges.clear();row.finiteContact=false;
    }
    if(x.size()!=4*n || !x.allFinite())throw std::invalid_argument("Invalid electrothermal state");
    for(const auto* r:{&eReference,&hReference})
        if(r->size()!=0 && (r->size()!=n || !r->allFinite()))throw std::invalid_argument("Invalid electrothermal QF reference");
    const auto ref=[&](Index i,bool e){const auto& r=e?eReference:hReference;return r.size()?r[i]:0.;};
    ElectrothermalAssembly out;out.residual=VectorXd::Zero(4*n);
    out.electronOutflow_A_per_m=VectorXd::Zero(n);out.holeOutflow_A_per_m=VectorXd::Zero(n);
    out.electronFluxAbs_A_per_m=VectorXd::Zero(n);out.holeFluxAbs_A_per_m=VectorXd::Zero(n);out.recombination_A_per_m=VectorXd::Zero(n);
    std::vector<Eigen::Triplet<Real>> entries;
    std::vector<bool> constrained(4*n,false);
    const std::array<const std::map<Index,Real>*,4> fixed{&bc.potential_V,&bc.electronQf_V,&bc.holeQf_V,&bc.temperature_K};
    for(int k=0;k<4;++k)for(const auto& [i,v]:*fixed[k]){
        if(i>=n || !std::isfinite(v))throw std::invalid_argument("Invalid electrothermal boundary");
        constrained[4*i+k]=true;
    }
    for(const auto& [i,v]:bc.neutralContactBias_V){
        if(i>=n || geometry_.siliconArea_m2[i]<=0. || !std::isfinite(v))
            throw std::invalid_argument("Invalid neutral contact");
        for(int k=0;k<3;++k){
            if(constrained[4*i+k])throw std::invalid_argument("Conflicting neutral contact constraint");
            constrained[4*i+k]=!(k==2 && bc.holeRecombination.contains(i));
        }
    }
    for(const auto& [i,c]:bc.holeRecombination)
        if(!bc.neutralContactBias_V.contains(i) || !std::isfinite(c.velocity_m_per_s) ||
           c.velocity_m_per_s<0. || !std::isfinite(c.boundaryLength_m) || c.boundaryLength_m<=0.)
            throw std::invalid_argument("Finite hole contact requires neutral bias, nonnegative SI velocity and positive boundary length");
    for(Index i=0;i<n;++i)if(geometry_.siliconArea_m2[i]==0.){
        constrained[4*i+1]=true;constrained[4*i+2]=true;
    }
    const auto add=[&](Index row,Index column,Real value){
        if(buildJacobian && !constrained[row] && value!=0.)entries.emplace_back(row,column,value);
    };
    std::vector<SiliconThermalState> states(n);std::vector<SiliconThermalResult> properties(n);
    VectorXd t(n),psi(n),fn(n),fp(n),ne=VectorXd::Zero(n),nh=ne,dne=ne,dnh=ne,dneT=ne,dnhT=ne;
    for(Index i=0;i<n;++i){
        psi[i]=x[4*i];fn[i]=ref(i,true)+x[4*i+1];fp[i]=ref(i,false)+x[4*i+2];t[i]=x[4*i+3];
        if(t[i]<50.)throw std::invalid_argument("Electrothermal temperature below 50 K");
        if(geometry_.siliconArea_m2[i]==0.)continue;
        states[i]={psi[i],x[4*i+1],x[4*i+2],t[i],doping_.donors(i),doping_.acceptors(i),ref(i,true),ref(i,false)};
        auto& p=properties[i];p=reusePreparation_?physics_.evaluate(states[i],preparedAt(i,t[i])):physics_.evaluate(states[i]);
        ne[i]=p.electrons_m3.value;nh[i]=p.holes_m3.value;
        dne[i]=p.electrons_m3.derivative[0];dnh[i]=p.holes_m3.derivative[2];
        dneT[i]=p.electrons_m3.derivative[3];dnhT[i]=p.holes_m3.derivative[3];
        const Real qArea=constants::q*geometry_.siliconArea_m2[i];
        const Real sourceArea=constants::q*geometry_.recombinationArea_m2[i];
        const Real source=sourceArea*(p.srhRate_m3_per_s.value+p.augerRate_m3_per_s.value);
        if(auto it=audited.find(i);it!=audited.end()){
            auto& row=*it->second;row.state=states[i];row.properties=p;row.sourceAreaCharge=sourceArea;row.source=source;
        }
        out.residual[4*i]-=qArea*(nh[i]-ne[i]+doping_.netDoping(i));
        out.recombination_A_per_m[i]=source;
        out.residual[4*i+1]-=source;out.residual[4*i+2]+=source;
        for(int k=0;k<4;++k){
            add(4*i,4*i+k,-qArea*(p.holes_m3.derivative[k]-p.electrons_m3.derivative[k]));
            const Real d=sourceArea*(p.srhRate_m3_per_s.derivative[k]+p.augerRate_m3_per_s.derivative[k]);
            if(!diagnosticFreezeRecombinationDerivatives){add(4*i+1,4*i+k,-d);add(4*i+2,4*i+k,d);}
        }
    }
    // Only fully constrained, isothermal, flat-QF Poisson prebias has exactly
    // zero edge current/work and no active transport Jacobian rows.
    // Compare references/increments separately; their rounded sums can hide
    // physical sub-ULP QF differences and nonzero edge currents.
    bool identicalQfRepresentation=skipEquilibriumTransport && n>0;
    for(Index i=0;i<n && identicalQfRepresentation;++i)
        identicalQfRepresentation=ref(i,true)==ref(0,true) && ref(i,false)==ref(0,true) &&
            x[4*i+1]==x[1] && x[4*i+2]==x[1];
    const bool skipTransport=skipEquilibriumTransport && identicalQfRepresentation &&
        bc.electronQf_V.size()==n && bc.holeQf_V.size()==n && bc.temperature_K.size()==n &&
        bc.neutralContactBias_V.empty() && bc.holeRecombination.empty() &&
        (t.array()==t[0]).all();
    if(!skipTransport)updateIalTransportState(mobility_,mesh_,doping_,psi,ne,nh,fn,fp,dne,dnh,t,dneT,dnhT,reuseIalScreening_,buildJacobian);
    for(Index i=0;i<n;++i)out.residual[4*i]-=geometry_.fixedCharge_C_per_m[i];
    for(const auto& edge:mesh_.edges()){
        const Index a=edge.n0,b=edge.n1;const Real eps=geometry_.poissonEdge_F_per_m[edge.id];
        const Real flux=eps*(psi[a]-psi[b]);out.residual[4*a]+=flux;out.residual[4*b]-=flux;
        add(4*a,4*a,eps);add(4*a,4*b,-eps);add(4*b,4*a,-eps);add(4*b,4*b,eps);
        const Real weight=geometry_.transportWeight[edge.id];if(weight==0. || skipTransport)continue;
        const auto& pa=properties[a];const auto& pb=properties[b];
        Real currents[2]{};std::array<std::map<Index,Real>,2> derivatives;
        for(int carrier=0;carrier<2;++carrier){
            const bool electron=carrier==0;
            const Real mu=mobility_.model=="ialmob"?ialEdgeMobility(mobility_,edge.id,
                electron?CarrierType::Electron:CarrierType::Hole):(electron?muE_:muH_);
            const auto current=thermalSgCurrent(states[a],pa,states[b],pb,mu,weight,electron);
            if(!electron && !audited.empty())for(Index node:{a,b})if(auto it=audited.find(node);it!=audited.end())
                it->second->edges.push_back({edge.id,a,b,states[a],states[b],pa,pb,mu,weight,current.current_A_per_m});
            currents[carrier]=current.current_A_per_m;auto& d=derivatives[carrier];
            if(buildJacobian)for(int k=0;k<4;++k){d[4*a+k]+=current.derivative[k];d[4*b+k]+=current.derivative[4+k];}
            if(buildJacobian && !diagnosticFreezeMobilityDerivatives && mobility_.model=="ialmob")for(const auto& support:mobility_.ialmobGeometry->edges[edge.id]){
                const auto& r=mobility_.ialmobState->cells[support.support];
                const auto& cell=mesh_.getCell(mobility_.ialmobGeometry->cells[support.support].cellId);
                const auto& md=electron?r.electron:r.hole;
                const auto& td=electron?r.electronTemperatureDerivative:r.holeTemperatureDerivative;
                const Real factor=current.mobilityDerivative*support.weight;
                for(int j=0;j<3;++j){
                    for(int k=0;k<3;++k)d[4*cell.node_ids[j]+k]+=factor*md.derivative[3*j+k];
                    d[4*cell.node_ids[j]+3]+=factor*td[j];
                }
            }
            const Index offset=electron?1:2;
            out.residual[4*a+offset]+=current.current_A_per_m;out.residual[4*b+offset]-=current.current_A_per_m;
            auto& terminal=electron?out.electronOutflow_A_per_m:out.holeOutflow_A_per_m;
            auto& fluxAbs=electron?out.electronFluxAbs_A_per_m:out.holeFluxAbs_A_per_m;
            fluxAbs[a]+=std::abs(current.current_A_per_m);fluxAbs[b]+=std::abs(current.current_A_per_m);
            terminal[a]+=current.current_A_per_m;terminal[b]-=current.current_A_per_m;
            for(const auto& [column,value]:d){add(4*a+offset,column,value);add(4*b+offset,column,-value);}
        }
        const auto power=latticeBandEdgeWork(currents[0],currents[1],pa.conductionBand_eV.value,
            pb.conductionBand_eV.value,pa.valenceBand_eV.value,pb.valenceBand_eV.value);
        out.latticeSource_W_per_m+=power.power_W_per_m;
        // Conservative endpoint deposition of this very same edge work.
        out.residual[4*a+3]-=.5*power.power_W_per_m;out.residual[4*b+3]-=.5*power.power_W_per_m;
        std::map<Index,Real> dp;
        for(int c=0;c<2;++c)for(const auto& [column,value]:derivatives[c])dp[column]+=power.derivative[c]*value;
        if(buildJacobian)for(int k=0;k<4;++k){
            dp[4*a+k]+=power.derivative[2]*pa.conductionBand_eV.derivative[k]+power.derivative[4]*pa.valenceBand_eV.derivative[k];
            dp[4*b+k]+=power.derivative[3]*pb.conductionBand_eV.derivative[k]+power.derivative[5]*pb.valenceBand_eV.derivative[k];
        }
        for(const auto& [column,value]:dp){add(4*a+3,column,-.5*value);add(4*b+3,column,-.5*value);}
    }
    const auto heat=heat_.assemble(t,VectorXd::Zero(mesh_.numCells()));
    out.boundaryHeat_W_per_m=heat.outward_boundary_heat_W_per_m;
    for(Index i=0;i<n;++i)out.residual[4*i+3]+=heat.residual_W_per_m[i];
    for(int k=0;k<heat.jacobian_W_per_m_K.outerSize();++k)
        for(SparseMatrixd::InnerIterator it(heat.jacobian_W_per_m_K,k);it;++it)
            add(4*it.row()+3,4*it.col()+3,it.value());
    for(Index i=0;i<n;++i)if(geometry_.siliconArea_m2[i]==0.)for(int k=1;k<=2;++k){
        if(!fixed[k]->contains(i)){out.residual[4*i+k]=x[4*i+k]+ref(i,k==1);entries.emplace_back(4*i+k,4*i+k,1.);}
    }
    for(int k=0;k<4;++k)for(const auto& [i,value]:*fixed[k]){
        out.residual[4*i+k]=x[4*i+k]-(value-(k==1?ref(i,true):k==2?ref(i,false):0.));entries.emplace_back(4*i+k,4*i+k,1.);
    }
    for(const auto& [i,bias]:bc.neutralContactBias_V){
        const auto [potential,derivative]=neutralPotential(i,bias,t[i]);
        out.residual[4*i]=psi[i]-potential;entries.emplace_back(4*i,4*i,1.);entries.emplace_back(4*i,4*i+3,-derivative);
        for(int k=1;k<=2;++k){
            if(k==2 && bc.holeRecombination.contains(i))continue;
            out.residual[4*i+k]=x[4*i+k]-(bias-ref(i,k==1));entries.emplace_back(4*i+k,4*i+k,1.);
        }
        if(const auto found=bc.holeRecombination.find(i);found!=bc.holeRecombination.end()){
            const auto& contact=found->second;
            const SiliconThermalState equilibriumState{potential,0.,0.,t[i],doping_.donors(i),doping_.acceptors(i),bias,bias};
            const auto equilibrium=reusePreparation_?physics_.evaluate(equilibriumState,preparedAt(i,t[i])):physics_.evaluate(equilibriumState);
            const auto& p=properties[i].holes_m3;
            const auto& p0=equilibrium.holes_m3;
            const Real vt=constants::kb*t[i]/constants::q;
            // Preserve sub-ULP QF increments and avoid subtracting two nearly
            // equal densities. The small-eta expansion has relative O(delta^2)
            // error below 1e-12 and exactly vanishes at equilibrium.
            const Real delta=static_cast<Real>((static_cast<long double>(ref(i,false))-bias+
                x[4*i+2]-(static_cast<long double>(psi[i])-potential))/vt);
            const Real eta=equilibrium.holeEta.value;
            const Real excess=std::abs(delta)<1e-6 ? equilibrium.Nv_m3.value*delta*
                (fermiDiracHalfDerivative(eta)+.5*delta*fermiDiracHalfSecondDerivative(eta)) : p.value-p0.value;
            const Real coefficient=constants::q*contact.velocity_m_per_s*contact.boundaryLength_m;
            // Positive hole loss is OUTWARD from the device; stored terminal
            // currents and UG Eq.103 use the contact-to-semiconductor direction.
            const Real outward=coefficient*excess;
            if(auto it=audited.find(i);it!=audited.end()){
                auto& row=*it->second;row.finiteContact=true;row.contactBias=bias;row.neutralPotential=potential;
                row.neutralTemperatureDerivative=derivative;row.contactVt=vt;row.contactDelta=delta;
                row.contactCoefficient=coefficient;row.contactOutward=outward;row.equilibriumNv=equilibrium.Nv_m3.value;
                row.equilibriumEta=eta;row.equilibriumDensity=p0.value;
                row.fermiDerivative=fermiDiracHalfDerivative(eta);row.fermiSecondDerivative=fermiDiracHalfSecondDerivative(eta);
            }
            out.residual[4*i+2]+=outward;
            out.holeFluxAbs_A_per_m[i]+=std::abs(outward);
            out.holeOutflow_A_per_m[i]=-outward;
            for(int k=0;k<4;++k){
                const Real equilibriumDerivative=k==3?p0.derivative[3]+p0.derivative[0]*derivative:0.;
                add(4*i+2,4*i+k,coefficient*(p.derivative[k]-equilibriumDerivative));
            }
        }
    }
    if(buildJacobian){out.jacobian.resize(4*n,4*n);out.jacobian.setFromTriplets(entries.begin(),entries.end());}
    for(const auto& [node,row]:audited){row->residual=out.residual[4*node+2];row->fluxAbs=out.holeFluxAbs_A_per_m[node];}
    return out;
}
} // namespace vela
