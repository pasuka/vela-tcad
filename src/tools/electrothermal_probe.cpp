// Experimental four-equation state audit; qualification is recorded externally.
#include "vela/equation/ElectrothermalAssembler.h"
#include "ElectrothermalIterationControl.h"
#include "vela/io/MeshReader.h"
#include <Eigen/SparseLU>
#include "vela/core/PhysicsCallCounters.h"
#include <chrono>
#include "vela/solver/NewtonSolver.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace vela;
using json=nlohmann::json;
namespace {
VectorXd values(const json& a) {
    const auto data=a.get<std::vector<Real>>();
    return Eigen::Map<const VectorXd>(data.data(),data.size());
}
std::vector<Real> list(const VectorXd& a) {return {a.data(),a.data()+a.size()};}
}
int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::invalid_argument("Usage: electrothermal_probe INPUT.json OUTPUT.json");
        if (std::filesystem::exists(argv[2])) throw std::invalid_argument("Output already exists");
        std::ifstream input(argv[1]);if (!input) throw std::invalid_argument("Cannot open input");
        json cfg;input>>cfg;
        JsonMeshReader reader;
        auto mesh=reader.read(cfg.at("mesh_file").get<std::string>());
        std::map<Index,LatticeConductivity> laws;
        for (const auto& r:cfg.at("region_conductivity")) {
            LatticeConductivity law;
            const auto model=r.at("model").get<std::string>();
            if (model=="constant") law.constant_W_per_m_K=r.at("value_W_per_m_K");
            else if (model=="inverse_quadratic") {
                law.model=LatticeConductivity::Model::InverseQuadratic;
                law.numerator=r.at("numerator");law.denominator=r.at("denominator").get<std::array<Real,3>>();
            } else throw std::invalid_argument("Unknown conductivity model");
            if (!laws.emplace(r.at("region_id").get<Index>(),law).second)
                throw std::invalid_argument("Duplicate region conductivity");
        }
        std::vector<LatticeThermodeEdge> edges;
        for (const auto& b:cfg.at("thermodes"))
            edges.push_back({b.at("nodes").get<std::array<Index,2>>(),b.at("ambient_K"),b.at("conductance_W_per_m2_K")});

        const Real lengthFactor=cfg.at("coordinate_to_metres");
        ElectrothermalGeometry geometry;
        geometry.siliconArea_m2=values(cfg.at("silicon_area_m2"));
        if(cfg.contains("recombination_area_m2"))geometry.recombinationArea_m2=values(cfg.at("recombination_area_m2"));
        geometry.fixedCharge_C_per_m=values(cfg.at("fixed_charge_C_per_m"));
        geometry.poissonEdge_F_per_m=VectorXd::Zero(mesh.numEdges());
        geometry.transportWeight=VectorXd::Zero(mesh.numEdges());
        std::map<std::pair<Index,Index>,json> inputEdges;
        for(const auto& e:cfg.at("edge_geometry")){
            const auto pair=e.at("nodes").get<std::array<Index,2>>();
            if(!inputEdges.emplace(std::minmax(pair[0],pair[1]),e).second)
                throw std::invalid_argument("Duplicate edge geometry");
        }
        if(inputEdges.size()!=mesh.numEdges())throw std::invalid_argument("Incomplete edge geometry");
        for(const auto& edge:mesh.edges()){
            const auto& e=inputEdges.at(std::minmax(edge.n0,edge.n1));
            geometry.poissonEdge_F_per_m[edge.id]=e.at("poisson_F_per_m");
            geometry.transportWeight[edge.id]=e.at("transport_weight");
            mesh.setTransportCouple(edge.id,geometry.transportWeight[edge.id]*edge.length);
        }
        VectorXd nd=values(cfg.at("donors_m3")),na=values(cfg.at("acceptors_m3"));
        if(nd.size()!=mesh.numNodes() || na.size()!=mesh.numNodes())throw std::invalid_argument("Doping size mismatch");
        DopingModel doping(mesh.numNodes());for(Index i=0;i<mesh.numNodes();++i)doping.setNodeDoping(i,nd[i],na[i]);
        auto mobility=mobilityModelConfigFromJson(cfg.at("mobility_SI"));
        mobility.internalLengthToM=lengthFactor;
        LatticeHeatAssembler heat(mesh,lengthFactor,laws,edges);
        SiliconThermalParameters siliconParameters;
        siliconParameters.augerWithGeneration=cfg.value("auger_with_generation",false);
        ElectrothermalAssembler assembler(mesh,doping,std::move(geometry),heat,mobility,SiliconThermalPhysics(siliconParameters),.1,.04,cfg.value("reuse_physics_preparation",false),cfg.value("reuse_ialmob_screening",false));
        ElectrothermalBoundary bc;
        for(const auto& b:cfg.at("boundaries")){
            const Index node=b.at("node");const std::string kind=b.at("kind");const Real value=b.at("value");
            auto* target=kind=="neutral_contact"?&bc.neutralContactBias_V:
                kind=="psi"?&bc.potential_V:kind=="fn"?&bc.electronQf_V:kind=="fp"?&bc.holeQf_V:kind=="temperature"?&bc.temperature_K:nullptr;
            if(!target || !target->emplace(node,value).second)throw std::invalid_argument("Invalid/duplicate boundary");
            if(b.contains("hole_recombination_velocity_m_per_s")){
                if(kind!="neutral_contact")throw std::invalid_argument("Hole recombination requires a neutral contact");
                bc.holeRecombination.emplace(node,ElectrothermalBoundary::HoleRecombination{
                    b.at("hole_recombination_velocity_m_per_s"),b.at("boundary_length_m")});
            }else if(b.contains("boundary_length_m"))throw std::invalid_argument("Boundary length requires a recombination velocity");
        }
        VectorXd x=values(cfg.at("state_interleaved"));
        if(x.size()!=4*mesh.numNodes() || !x.allFinite())
            throw std::invalid_argument("Invalid electrothermal state size or values");



        VectorXd eReference=VectorXd::Zero(mesh.numNodes()),hReference=eReference;
        if(cfg.value("use_qf_references",true)){
            if(cfg.contains("referenced_state_interleaved")){
                x=values(cfg.at("referenced_state_interleaved"));
                eReference=values(cfg.at("electron_qf_reference_V"));hReference=values(cfg.at("hole_qf_reference_V"));
                if(x.size()!=4*mesh.numNodes() || eReference.size()!=mesh.numNodes() || hReference.size()!=mesh.numNodes())
                    throw std::invalid_argument("Referenced input size mismatch");
            }else for(Index i=0;i<mesh.numNodes();++i){
                eReference[i]=x[4*i+1];hReference[i]=x[4*i+2];x[4*i+1]=0.;x[4*i+2]=0.;
            }
        }

        // When enabling finite exchange on a previously qualified ideal-contact
        // seed, impose its prescribed potential exactly before Newton. Otherwise
        // an admissible old Dirichlet residual can create a spurious p-p0 flux
        // even at zero bias. This does not project the now-free hole QF.
        if(!bc.holeRecombination.empty()){
            if(!x.allFinite() || !eReference.allFinite() || !hReference.allFinite())
                throw std::invalid_argument("Invalid referenced state values");
            for(const auto& [i,contact]:bc.holeRecombination){
                if(i>=mesh.numNodes())throw std::invalid_argument("Invalid finite contact node");
                x[4*i]=assembler.neutralPotential(i,bc.neutralContactBias_V.at(i),x[4*i+3]).first;
            }
        }
        const auto recenter=[&](){
            bool changed=false;
            if(!cfg.value("use_qf_references",true))return changed;
            for(Index i=0;i<mesh.numNodes();++i)for(int k=1;k<=2;++k){
                Real& reference=k==1?eReference[i]:hReference[i];Real& increment=x[4*i+k];
                if(std::abs(increment)<=1e-3)continue;
                const Real next=static_cast<Real>(static_cast<long double>(reference)+increment);
                increment=static_cast<Real>((static_cast<long double>(reference)-next)+increment);
                reference=next;changed=true;
            }
            return changed;
        };
        recenter();
        const VectorXd siliconArea=values(cfg.at("silicon_area_m2"));
        const Real currentScale=cfg.value("electrical_current_scale_A_per_m",1.);
        const auto gateConfig=newtonConfigFromJson(cfg.value("electrical_gate_solver",json::object()));
        const auto rowGate=[&](const ElectrothermalAssembly& result){
            std::vector<CoupledDDCarrierTermDiagnostic> rows(mesh.numNodes());
            for(Index i=0;i<mesh.numNodes();++i){auto& row=rows[i];row.nodeId=i;
                row.electronContinuityActive=siliconArea[i]>0. && !bc.neutralContactBias_V.contains(i) && !bc.electronQf_V.contains(i);
                row.holeContinuityActive=siliconArea[i]>0. && (!bc.neutralContactBias_V.contains(i) || bc.holeRecombination.contains(i)) && !bc.holeQf_V.contains(i);
                row.electronFluxAbsSum=result.electronFluxAbs_A_per_m[i]/currentScale;
                row.holeFluxAbsSum=result.holeFluxAbs_A_per_m[i]/currentScale;
                row.electronRecombination=row.holeRecombination=result.recombination_A_per_m[i]/currentScale;
                row.electronResidual=result.residual[4*i+1]/currentScale;row.holeResidual=result.residual[4*i+2]/currentScale;
            }
            return evaluateCarrierRowConvergence(rows,gateConfig.carrierRowConvergence);
        };
        const auto blockGates=[&](const ElectrothermalAssembly& result){
        std::array<Real,3> rawNorm{},weightedNorm{};
        const Real vt=constants::Vt_300,poissonScale=11.7*constants::eps0*vt;
        for(Index i=0;i<mesh.numNodes();++i)for(int k=0;k<3;++k){
            const bool constrained=(bc.neutralContactBias_V.contains(i) && !(k==2 && bc.holeRecombination.contains(i))) || (k==0?bc.potential_V.contains(i):
                k==1?bc.electronQf_V.contains(i)||siliconArea[i]==0.:bc.holeQf_V.contains(i)||siliconArea[i]==0.);
            const Real r=result.residual[4*i+k]/(constrained?vt:k==0?poissonScale:currentScale);
            Real weight=1.;const auto& scaling=gateConfig.continuityRowScaling;
            if(k>0 && !constrained && scaling.enabled){
                const Real source=std::abs(result.recombination_A_per_m[i])/currentScale;
                const Real flux=(k==1?result.electronFluxAbs_A_per_m[i]:result.holeFluxAbs_A_per_m[i])/currentScale;
                if(source>=scaling.minSourceScale && source>=scaling.fluxFraction*flux)
                    weight=std::clamp(1./std::max(source,scaling.scaleFloor),scaling.minWeight,scaling.maxWeight);
            }
            rawNorm[k]+=r*r;weightedNorm[k]+=r*r*weight*weight;
        }
        json blockGate=json::array();
        const std::array<Real,3> ceilings{gateConfig.blockAbsoluteConvergence.psiResidualCeiling,
            gateConfig.blockAbsoluteConvergence.electronResidualCeiling,gateConfig.blockAbsoluteConvergence.holeResidualCeiling};
        for(int k=0;k<3;++k)blockGate.push_back({{"raw_l2",std::sqrt(rawNorm[k])},{"weighted_l2",std::sqrt(weightedNorm[k])},
            {"limit",ceilings[k]},{"satisfied",std::sqrt(weightedNorm[k])<=ceilings[k]}});
            return blockGate;
        };
        const auto blocksSatisfied=[&](const ElectrothermalAssembly& result){const auto b=blockGates(result);return std::all_of(b.begin(),b.end(),[](const json& v){return v.at("satisfied").get<bool>();});};
        const bool profiling=cfg.value("performance_profiling",false);
        using Clock=std::chrono::steady_clock;
        const auto seconds=[](Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();};
        Real assemblySeconds=0.,factorizationSeconds=0.,solveSeconds=0.;
        unsigned assemblyCalls=0,factorizations=0;
        const auto countersBefore=physicsCallCounters;
        const auto assemble=[&](const VectorXd& state){
            const auto start=Clock::now();++assemblyCalls;
            try{auto result=assembler.assemble(state,bc,eReference,hReference);
                assemblySeconds+=seconds(start);return result;}
            catch(...){assemblySeconds+=seconds(start);throw;}
        };
        auto a=assemble(x);

        const unsigned maximum=cfg.value("diagnostic_newton_max_iterations",0u);
        const bool reuseSymbolic=cfg.value("reuse_sparselu_symbolic",false);
        if(cfg.contains("diagnostic_stagnation_window") && !cfg.at("diagnostic_stagnation_window").is_number_unsigned())
            throw std::invalid_argument("Stagnation window must be a nonnegative integer");
        experimental::ElectrothermalSparseLU lu;
        experimental::ElectrothermalStagnationWatch stagnation(cfg.value("diagnostic_stagnation_window",0u));
        json history=json::array();std::string stop="frozen_state";
        VectorXd columns(x.size());for(int i=0;i<x.size();++i)columns[i]=i%4==3?10.:.025;
        for(unsigned iteration=0;iteration<maximum;++iteration){
            if(recenter())a=assemble(x);
            VectorXd rows=VectorXd::Zero(x.size());
            for(int k=0;k<a.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(a.jacobian,k);it;++it)
                rows[it.row()]+=std::abs(it.value()*columns[it.col()]);
            rows=rows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
            const auto merit=[&](const ElectrothermalAssembly& r){return (r.residual.array()*rows.array()).matrix().norm();};
            const Real before=merit(a);if(before<1e-9 && rowGate(a).satisfied && blocksSatisfied(a)){stop="diagnostic_scaled_residual";break;}
            SparseMatrixd matrix=a.jacobian;
            for(int k=0;k<matrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(matrix,k);it;++it)
                it.valueRef()*=rows[it.row()]*columns[it.col()];
            const auto factorStart=Clock::now();
            lu.compute(matrix,reuseSymbolic);
            factorizationSeconds+=seconds(factorStart);++factorizations;
            if(lu.info()!=Eigen::Success){stop="factorization_failed";break;}
            const VectorXd rhs=-(a.residual.array()*rows.array()).matrix();
            const auto solveStart=Clock::now();
            VectorXd direction=lu.solve(rhs);direction.array()*=columns.array();
            solveSeconds+=seconds(solveStart);
            if(lu.info()!=Eigen::Success || !direction.allFinite()){stop="linear_solve_failed";break;}
            Real alpha=1.;int limiter=-1;
            std::array<Real,4> maxDirection{};
            for(int i=0;i<x.size();++i)if(direction[i]!=0.){
                const Real cap=(i%4==3?30.:.2)/std::abs(direction[i]);
                if(cap<alpha)limiter=i;
                alpha=std::min(alpha,cap);
                maxDirection[i%4]=std::max(maxDirection[i%4],std::abs(direction[i]));
            }
            const Real initialAlpha=alpha;
            unsigned trials=0;
            bool accepted=false;
            for(int trial=0;trial<24;++trial){
                ++trials;
                VectorXd candidate=x+alpha*direction;
                try{auto next=assemble(candidate);
                    if(merit(next)<before || (before<1e-9 && merit(next)<1e-9 && rowGate(next).maxRatio<rowGate(a).maxRatio)){x=candidate;a=std::move(next);accepted=true;break;}}
                catch(const std::exception&){}
                alpha*=.5;
            }
            history.push_back({{"iteration",iteration+1},{"scaled_l2_before",before},{"scaled_l2_after",merit(a)},{"alpha",alpha},{"accepted",accepted}});
            if(profiling){
                history.back()["initial_alpha"]=initialAlpha;
                history.back()["line_search_trials"]=trials;
                history.back()["limiter_node"]=limiter<0?-1:limiter/4;
                history.back()["limiter_component"]=limiter<0?-1:limiter%4;
                history.back()["limiter_direction"]=limiter<0?0.:direction[limiter];
                history.back()["max_abs_direction_by_block"]=maxDirection;
            }
            std::cout<<history.back().dump()<<std::endl;
            if(!accepted){stop="line_search_failed";break;}
            if(stagnation.update(before,merit(a),alpha,accepted)){stop="diagnostic_stagnation_reject";break;}
            stop="diagnostic_iteration_limit";
        }

        const auto rows=rowGate(a);
        json violations=json::array();for(const auto& row:rows.violations)if(violations.size()<20)violations.push_back({{"node",row.nodeId},{"carrier",row.carrier},{"ratio",row.ratio},{"residual_scaled",row.residual},{"scale",row.scale}});

        const auto blockGate=blockGates(a);
        json blocks=json::array(),contacts=json::array();

        for(int k=0;k<4;++k){Real maximum=0.,squared=0.;
            for(Index i=0;i<mesh.numNodes();++i){const Real r=a.residual[4*i+k];maximum=std::max(maximum,std::abs(r));squared+=r*r;}
            blocks.push_back({{"inf",maximum},{"l2",std::sqrt(squared)}});
        }
        for(const auto& c:mesh.contacts()){
            if(c.name=="th_lat")continue;
            Real ne=0.,nh=0.;for(Index i:c.node_ids){ne+=a.electronOutflow_A_per_m[i];nh+=a.holeOutflow_A_per_m[i];}
            contacts.push_back({{"contact",c.name},{"electron_outflow_A_per_m",ne},{"hole_outflow_A_per_m",nh},{"total_outflow_A_per_m",ne+nh}});
        }
        VectorXd t(mesh.numNodes());for(Index i=0;i<mesh.numNodes();++i)t[i]=x[4*i+3];
        VectorXd physicalState=x;
        for(Index i=0;i<mesh.numNodes();++i){physicalState[4*i+1]+=eReference[i];physicalState[4*i+2]+=hReference[i];}
        const auto appliedUpdates=std::count_if(history.begin(),history.end(),[](const json& v){return v.at("accepted").get<bool>();});
        json result={
            {"schema","vela.experimental_electrothermal_probe.v1"},
            {"potential_origin_V",cfg.value("potential_origin_V",0.)},
            {"auger_with_generation",siliconParameters.augerWithGeneration},
            {"scope","Experimental four-equation diagnostic; scaled Newton criterion is not the LDMOS acceptance gate"},
            {"newton_updates",appliedUpdates},{"newton_attempts",history.size()},{"diagnostic_stop",stop},{"history",history},{"state_interleaved",list(physicalState)},{"referenced_state_interleaved",list(x)},{"electron_qf_reference_V",list(eReference)},{"hole_qf_reference_V",list(hReference)},{"contacts",contacts},{"residual_blocks",blocks},
            {"electrical_block_gates",blockGate},
            {"carrier_row_gate",{{"satisfied",rows.satisfied},{"max_ratio",rows.maxRatio},{"eps_row",rows.epsRow},{"qualified_rows",rows.qualifiedRowCount},{"violations",violations}}},
            {"jacobian_nonzeros",a.jacobian.nonZeros()},
            {"lattice_source_W_per_m",a.latticeSource_W_per_m},{"boundary_heat_W_per_m",a.boundaryHeat_W_per_m},
            {"temperature_K",list(t)},{"nodal_area_m2",list(heat.nodalAreas_m2())},
            {"residual",list(a.residual)}};
        if(profiling)result["performance"]={
            {"scope","Newton phase wall times; factorization includes symbolic analysis; counters do not change decisions"},
            {"assembly_seconds",assemblySeconds},{"assembly_calls",assemblyCalls},
            {"factorization_seconds",factorizationSeconds},{"factorizations",factorizations},
            {"linear_solve_seconds",solveSeconds},
            {"fermi_half_calls",physicsCallCounters.fermiDiracHalf-countersBefore.fermiDiracHalf},
            {"fermi_half_derivative_calls",physicsCallCounters.fermiDiracHalfDerivative-countersBefore.fermiDiracHalfDerivative},
            {"inverse_fermi_half_calls",physicsCallCounters.inverseFermiDiracHalf-countersBefore.inverseFermiDiracHalf}};
        if(profiling)result["performance"]["preparation_counts"]=assembler.preparationCounts();
        if(profiling)result["performance"]["symbolic_analyses"]=lu.analyses();
        if(profiling)result["performance"]["ialmob_screening_minimum_solves"]=physicsCallCounters.ialScreeningMinimumSolves-countersBefore.ialScreeningMinimumSolves;
        if(cfg.contains("diagnostic_jacobian_direction")){
            const VectorXd direction=values(cfg.at("diagnostic_jacobian_direction"));
            const Real step=cfg.value("diagnostic_jacobian_step",1e-6);
            if(direction.size()!=x.size() || !direction.allFinite() || !std::isfinite(step) || step<=0.)
                throw std::invalid_argument("Invalid Jacobian audit direction/step");
            const auto plus=assembler.assemble(x+step*direction,bc,eReference,hReference);
            const auto minus=assembler.assemble(x-step*direction,bc,eReference,hReference);
            result["jacobian_audit"]={{"step",step},{"action",list(a.jacobian*direction)},
                {"central_difference",list((plus.residual-minus.residual)/(2.*step))}};
        }
        std::ofstream output(argv[2]);if(!output)throw std::runtime_error("Cannot open output");
        output<<result.dump(2)<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
