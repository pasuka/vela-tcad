// Experimental four-equation state audit; qualification is recorded externally.
#include "vela/equation/ElectrothermalAssembler.h"
#include "vela/io/MeshReader.h"
#include <Eigen/SparseLU>
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
        ElectrothermalAssembler assembler(mesh,doping,std::move(geometry),heat,mobility);
        ElectrothermalBoundary bc;
        for(const auto& b:cfg.at("boundaries")){
            const Index node=b.at("node");const std::string kind=b.at("kind");const Real value=b.at("value");
            auto* target=kind=="neutral_contact"?&bc.neutralContactBias_V:
                kind=="psi"?&bc.potential_V:kind=="fn"?&bc.electronQf_V:kind=="fp"?&bc.holeQf_V:kind=="temperature"?&bc.temperature_K:nullptr;
            if(!target || !target->emplace(node,value).second)throw std::invalid_argument("Invalid/duplicate boundary");
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
                row.holeContinuityActive=siliconArea[i]>0. && !bc.neutralContactBias_V.contains(i) && !bc.holeQf_V.contains(i);
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
            const bool constrained=bc.neutralContactBias_V.contains(i) || (k==0?bc.potential_V.contains(i):
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
        auto a=assembler.assemble(x,bc,eReference,hReference);

        const unsigned maximum=cfg.value("diagnostic_newton_max_iterations",0u);
        json history=json::array();std::string stop="frozen_state";
        VectorXd columns(x.size());for(int i=0;i<x.size();++i)columns[i]=i%4==3?10.:.025;
        for(unsigned iteration=0;iteration<maximum;++iteration){
            if(recenter())a=assembler.assemble(x,bc,eReference,hReference);
            VectorXd rows=VectorXd::Zero(x.size());
            for(int k=0;k<a.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(a.jacobian,k);it;++it)
                rows[it.row()]+=std::abs(it.value()*columns[it.col()]);
            rows=rows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
            const auto merit=[&](const ElectrothermalAssembly& r){return (r.residual.array()*rows.array()).matrix().norm();};
            const Real before=merit(a);if(before<1e-9 && rowGate(a).satisfied && blocksSatisfied(a)){stop="diagnostic_scaled_residual";break;}
            SparseMatrixd matrix=a.jacobian;
            for(int k=0;k<matrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(matrix,k);it;++it)
                it.valueRef()*=rows[it.row()]*columns[it.col()];
            Eigen::SparseLU<SparseMatrixd> lu;lu.compute(matrix);
            if(lu.info()!=Eigen::Success){stop="factorization_failed";break;}
            const VectorXd rhs=-(a.residual.array()*rows.array()).matrix();
            VectorXd direction=lu.solve(rhs);direction.array()*=columns.array();
            if(lu.info()!=Eigen::Success || !direction.allFinite()){stop="linear_solve_failed";break;}
            Real alpha=1.;
            for(int i=0;i<x.size();++i)if(direction[i]!=0.)alpha=std::min(alpha,(i%4==3?30.:.2)/std::abs(direction[i]));
            bool accepted=false;
            for(int trial=0;trial<24;++trial){
                VectorXd candidate=x+alpha*direction;
                try{auto next=assembler.assemble(candidate,bc,eReference,hReference);
                    if(merit(next)<before || (before<1e-9 && merit(next)<1e-9 && rowGate(next).maxRatio<rowGate(a).maxRatio)){x=candidate;a=std::move(next);accepted=true;break;}}
                catch(const std::exception&){}
                alpha*=.5;
            }
            history.push_back({{"iteration",iteration+1},{"scaled_l2_before",before},{"scaled_l2_after",merit(a)},{"alpha",alpha},{"accepted",accepted}});
            std::cout<<history.back().dump()<<std::endl;
            if(!accepted){stop="line_search_failed";break;}
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
            {"scope","Experimental four-equation diagnostic; scaled Newton criterion is not the LDMOS acceptance gate"},
            {"newton_updates",appliedUpdates},{"newton_attempts",history.size()},{"diagnostic_stop",stop},{"history",history},{"state_interleaved",list(physicalState)},{"referenced_state_interleaved",list(x)},{"electron_qf_reference_V",list(eReference)},{"hole_qf_reference_V",list(hReference)},{"contacts",contacts},{"residual_blocks",blocks},
            {"electrical_block_gates",blockGate},
            {"carrier_row_gate",{{"satisfied",rows.satisfied},{"max_ratio",rows.maxRatio},{"eps_row",rows.epsRow},{"qualified_rows",rows.qualifiedRowCount},{"violations",violations}}},
            {"jacobian_nonzeros",a.jacobian.nonZeros()},
            {"lattice_source_W_per_m",a.latticeSource_W_per_m},{"boundary_heat_W_per_m",a.boundaryHeat_W_per_m},
            {"temperature_K",list(t)},{"nodal_area_m2",list(heat.nodalAreas_m2())},
            {"residual",list(a.residual)}};
        std::ofstream output(argv[2]);if(!output)throw std::runtime_error("Cannot open output");
        output<<result.dump(2)<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
