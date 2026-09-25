// Experimental four-equation state audit; qualification is recorded externally.
#include "vela/equation/ElectrothermalAssembler.h"
#include "vela/solver/ElectrothermalIterationControl.h"
#include "vela/solver/ElectrothermalDensityUpdate.h"
#include "vela/solver/ElectrothermalNearSteady.h"
#include "vela/solver/ElectrothermalNaturalDamping.h"
#include "vela/solver/ElectrothermalPseudoTransient.h"
#include "vela/solver/ElectrothermalPredictorQuality.h"
#include "vela/solver/ElectrothermalTangent.h"
#include "vela/solver/ElectrothermalResidualMixing.h"
#include "vela/simulation/ElectrothermalSimulation.h"
#include "vela/io/MeshReader.h"
#include <Eigen/SparseLU>
#include "vela/core/PhysicsCallCounters.h"
#include "vela/core/IalKernelProfiling.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/IalMobilityJson.h"
#include <chrono>
#include "vela/solver/NewtonSolver.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
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
namespace {
using PreparationClock=std::chrono::steady_clock;
std::string sourceBytes(const std::string& path) {
    std::ifstream file(path,std::ios::binary);
    if(!file)throw std::runtime_error("Cannot read preparation dependency: "+path);
    std::string result((std::istreambuf_iterator<char>(file)),{});
    if(file.bad())throw std::runtime_error("Failed reading preparation dependency: "+path);
    return result;
}
std::string preparationIdentity(const json& cfg) {
    json key=json::object();
    for(const auto* name:{"mesh_file","coordinate_to_metres","region_conductivity","thermodes",
        "silicon_area_m2","recombination_area_m2","fixed_charge_C_per_m","edge_geometry",
        "donors_m3","acceptors_m3","mobility_SI","reuse_ialmob_local_preparation",
        "residual_ialmob_values_only","reuse_ialmob_thermal_high_field",
        "diagnostic_ialmob_explicit_high_field","diagnostic_ialmob_generated_low_field",
        "diagnostic_ialmob_screening_method"})if(cfg.contains(name))key[name]=cfg.at(name);
    // Compare complete bytes, not timestamps or a hash with possible collisions.
    // Mesh bytes include contacts; mobility JSON includes crystal axes and units.
    key["mesh_source_bytes"]=sourceBytes(cfg.at("mesh_file"));
    const auto& mobility=cfg.at("mobility_SI");
    if(mobility.contains("ialmob"))
        key["ialmob_source_bytes"]=sourceBytes(mobility.at("ialmob").at("geometry_file"));
    return key.dump();
}
}
struct vela::ElectrothermalPreparationContext::Impl {
    using Clock=PreparationClock;
    static double elapsed(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
    std::string identity;
    DeviceMesh mesh;
    DopingModel doping{0};
    VectorXd nd,na;
    ElectrothermalGeometry geometry;
    MobilityModelConfig mobility;
    std::unique_ptr<LatticeHeatAssembler> heat;
    std::shared_ptr<ElectrothermalAssembler::StructureCache> structure;
    double meshSeconds=0.,inputSeconds=0.,geometrySeconds=0.;
    explicit Impl(const json& cfg) {
        auto stage=Clock::now();
        JsonMeshReader reader;
        mesh=reader.read(cfg.at("mesh_file").get<std::string>());
        meshSeconds=elapsed(stage);stage=Clock::now();
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
        nd=values(cfg.at("donors_m3"));na=values(cfg.at("acceptors_m3"));
        if(nd.size()!=mesh.numNodes() || na.size()!=mesh.numNodes())throw std::invalid_argument("Doping size mismatch");
        doping=DopingModel(mesh.numNodes());for(Index i=0;i<mesh.numNodes();++i)doping.setNodeDoping(i,nd[i],na[i]);
        auto mobilityInput=cfg.at("mobility_SI");
        const auto screeningMethod=ial_json::electrothermalScreeningMethod(cfg);
        if (mobilityInput.is_object() && mobilityInput.contains("ialmob"))
            mobilityInput["ialmob"]["screening_method"]=ialScreeningMethodName(screeningMethod);
        mobility=mobilityModelConfigFromJson(mobilityInput);
        if(mobility.ialmob){
            auto options=std::make_shared<IalTransportOptions>(*mobility.ialmob);
            options->element.reuseLocalPreparation=cfg.value("reuse_ialmob_local_preparation",false);
            options->element.residualValuesOnly=cfg.value("residual_ialmob_values_only",false);
            options->element.reuseThermalHighField=cfg.value("reuse_ialmob_thermal_high_field",false);
            options->element.explicitHighFieldPartials=cfg.value("diagnostic_ialmob_explicit_high_field",false);
            options->element.generatedLowFieldPartials=cfg.value("diagnostic_ialmob_generated_low_field",false);
            mobility.ialmob=std::move(options);
        }
        mobility.internalLengthToM=lengthFactor;
        heat=std::make_unique<LatticeHeatAssembler>(mesh,lengthFactor,laws,edges);
        inputSeconds=elapsed(stage);stage=Clock::now();
        prepareIalTransportGeometry(mobility,mesh);
        geometrySeconds=elapsed(stage);
    }
};
nlohmann::json vela::solveElectrothermalPoint(const nlohmann::json& cfg, std::ostream& progress,
    ElectrothermalPreparationContext* context) {
        // Failed or throwing point requests must not leave reusable linear state.
        struct LinearRequestGuard {
            ElectrothermalPreparationContext* context;
            bool accepted=false;
            ~LinearRequestGuard(){if(context && !accepted)context->clearLinearContext();}
        } linearGuard{context};
        const bool reuseLinear=cfg.value("reuse_linear_analysis",true);
        if(context && !reuseLinear)context->clearLinearContext();
        IalKernelProfilingScope ialProfileScope(cfg.value("diagnostic_ialmob_kernel_timing",false));
        ElectrothermalCostScope costScope(cfg.value("diagnostic_electrothermal_cost",std::string("off")));
        const auto preparationStart=PreparationClock::now();
        using Preparation=ElectrothermalPreparationContext::Impl;
        const auto keyStart=PreparationClock::now();
        std::string identity=context?preparationIdentity(cfg):std::string{};
        const double preparationKeySeconds=Preparation::elapsed(keyStart);
        const bool preparationHit=context && context->prepared_ && context->prepared_->identity==identity;
        auto prepared=preparationHit?context->prepared_:std::make_shared<Preparation>(cfg);
        if(context && !preparationHit) {
            context->clearLinearContext();
            if(preparationIdentity(cfg)!=identity)throw std::runtime_error("Preparation inputs changed during construction");
            prepared->identity=std::move(identity);context->prepared_=prepared;
        }
        const auto& mesh=prepared->mesh;const auto& doping=prepared->doping;
        const auto& nd=prepared->nd;const auto& na=prepared->na;
        const auto& heat=*prepared->heat;
        const auto assemblerStart=PreparationClock::now();
        SiliconThermalParameters siliconParameters;
        siliconParameters.augerWithGeneration=cfg.value("auger_with_generation",false);
        ElectrothermalAssembler assembler(mesh,doping,prepared->geometry,heat,prepared->mobility,SiliconThermalPhysics(siliconParameters),.1,.04,cfg.value("reuse_physics_preparation",false),cfg.value("reuse_ialmob_screening",false),cfg.value("reuse_neutral_contact_roots",false),cfg.value("diagnostic_neutral_root_newton",false));
        const bool reuseStructure=cfg.value("reuse_jacobian_structure",true);
        if(reuseStructure){
            if(!prepared->structure)prepared->structure=std::make_shared<ElectrothermalAssembler::StructureCache>();
            assembler.setStructureCache(prepared->structure);
        }else prepared->structure.reset();
        const auto structureBefore=reuseStructure?std::array<std::size_t,4>{prepared->structure->coupled->builds,
            prepared->structure->coupled->hits,prepared->structure->heat->builds,prepared->structure->heat->hits}:std::array<std::size_t,4>{};
        const double assemblerPreparationSeconds=Preparation::elapsed(assemblerStart);
        const double pointPreparationSeconds=Preparation::elapsed(preparationStart);
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
        const std::string initialization=cfg.value("initialization",std::string("provided_state"));
        const std::string solveMode=cfg.value("solve_mode",std::string("coupled"));
        const bool frozenRowAudit=cfg.contains("diagnostic_hole_row_audit_nodes");
        if(frozenRowAudit && (cfg.value("diagnostic_newton_max_iterations",0u)!=0 ||
           solveMode!="coupled" || initialization!="provided_state"))
            throw std::invalid_argument("Hole row audit requires a zero-update coupled provided state");
        if(solveMode!="coupled" && solveMode!="poisson")throw std::invalid_argument("Unknown electrothermal solve_mode");
        VectorXd x;
        if(initialization=="neutral_300K") {
            if(cfg.contains("state_interleaved")||cfg.contains("referenced_state_interleaved"))
                throw std::invalid_argument("Neutral initialization cannot also supply a state");
            const Real bias=-cfg.value("potential_origin_V",0.);
            const auto areas=values(cfg.at("silicon_area_m2"));
            x=VectorXd::Zero(4*mesh.numNodes());
            for(Index i=0;i<mesh.numNodes();++i) {
                x[4*i]=areas[i]>0.?assembler.neutralPotential(i,bias,300.).first:bias;
                x[4*i+1]=x[4*i+2]=bias;x[4*i+3]=300.;
            }
        } else if(initialization=="provided_state") x=values(cfg.at("state_interleaved"));
        else throw std::invalid_argument("Unknown electrothermal initialization");
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
        if(!bc.holeRecombination.empty() && !frozenRowAudit){
            if(!x.allFinite() || !eReference.allFinite() || !hReference.allFinite())
                throw std::invalid_argument("Invalid referenced state values");
            for(const auto& [i,contact]:bc.holeRecombination){
                if(i>=mesh.numNodes())throw std::invalid_argument("Invalid finite contact node");
                x[4*i]=assembler.neutralPotential(i,bc.neutralContactBias_V.at(i),x[4*i+3]).first;
            }
        }
        // Poisson-only initialization/prebias keeps both quasi-Fermi potentials
        // and temperature fixed, matching the original inactive carrier solves.
        if(solveMode=="poisson") {
            for(const auto& [i,bias]:bc.neutralContactBias_V)
                bc.potential_V[i]=assembler.neutralPotential(i,bias,x[4*i+3]).first;
            bc.neutralContactBias_V.clear();bc.holeRecombination.clear();
            for(Index i=0;i<mesh.numNodes();++i) {
                bc.electronQf_V[i]=eReference[i]+x[4*i+1];
                bc.holeQf_V[i]=hReference[i]+x[4*i+2];
                bc.temperature_K[i]=x[4*i+3];
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
        // Frozen arithmetic observations must preserve even sub-ULP boundary
        // defects and the supplied reference/increment representation.
        if(!frozenRowAudit)recenter();
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
        Real fullAssemblySeconds=0.,residualAssemblySeconds=0.;
        unsigned assemblyCalls=0,factorizations=0,residualOnlyCalls=0;
        const auto countersBefore=physicsCallCounters;
        const auto assemble=[&](const VectorXd& state,bool buildJacobian=true){
            const auto start=Clock::now();++assemblyCalls;if(!buildJacobian)++residualOnlyCalls;
            try{auto result=assembler.assemble(state,bc,eReference,hReference,buildJacobian,
                    solveMode=="poisson" && cfg.value("skip_equilibrium_poisson_transport",false));
                const Real elapsed=seconds(start);assemblySeconds+=elapsed;
                (buildJacobian?fullAssemblySeconds:residualAssemblySeconds)+=elapsed;return result;}
            catch(...){const Real elapsed=seconds(start);assemblySeconds+=elapsed;
                (buildJacobian?fullAssemblySeconds:residualAssemblySeconds)+=elapsed;throw;}
        };
        auto a=assemble(x);
        json predictorCandidates=cfg.value("diagnostic_predictor_candidates",json::array());
        json tangentPreparation;unsigned tangentAnalyses=0;
        if(cfg.contains("diagnostic_tangent_predictor")) {
            if(solveMode!="coupled" || initialization!="provided_state" || !cfg.value("use_qf_references",true))
                throw std::invalid_argument("Tangent requires a referenced coupled provided state");
            const auto start=Clock::now();tangentPreparation={{"prepared",false}};
            try {
                const auto& source=cfg.at("diagnostic_tangent_predictor");
                const Real delta=source.at("delta_bias_V"),bias=source.at("source_bias_V");
                if(!(delta>0. && std::isfinite(delta) && std::isfinite(bias)) ||
                   source.value("potential_origin_V",0.)!=cfg.value("potential_origin_V",0.))
                    throw std::invalid_argument("Invalid tangent bias or origin");
                const auto nodes=source.at("moving_nodes").get<std::set<Index>>();
                auto sourceBoundary=bc;
                const VectorXd derivative=experimental::electrothermalContactBiasDerivative(mesh.numNodes(),bc,nodes);
                for(Index i:nodes) {
                    if(std::abs(bc.neutralContactBias_V.at(i)-bias-delta)>1e-12)
                        throw std::invalid_argument("Tangent source and target bias mismatch");
                    sourceBoundary.neutralContactBias_V[i]=bias;
                }
                VectorXd seed=values(source.at("referenced_state_interleaved"));
                const VectorXd er=values(source.at("electron_qf_reference_V")),hr=values(source.at("hole_qf_reference_V"));
                if(seed.size()!=x.size() || er.size()!=eReference.size() || hr.size()!=hReference.size() ||
                   !seed.allFinite() || !er.allFinite() || !hr.allFinite())
                    throw std::invalid_argument("Invalid tangent source state");
                for(Index i=0;i<mesh.numNodes();++i) {
                    seed[4*i+1]=static_cast<Real>((static_cast<long double>(er[i])-eReference[i])+seed[4*i+1]);
                    seed[4*i+2]=static_cast<Real>((static_cast<long double>(hr[i])-hReference[i])+seed[4*i+2]);
                }
                // Assemble at the accepted source bias/state, before any target
                // contact projection. Its final Jacobian is freshly factorized.
                const auto assemblyStart=Clock::now();++assemblyCalls;
                ElectrothermalAssembly sourceAssembly;
                try {sourceAssembly=assembler.assemble(seed,sourceBoundary,eReference,hReference);}
                catch(...) {assemblySeconds+=seconds(assemblyStart);throw;}
                assemblySeconds+=seconds(assemblyStart);
                VectorXd columns(x.size()),rows=VectorXd::Zero(x.size());
                for(int i=0;i<x.size();++i)columns[i]=i%4==3?10.:.025;
                auto matrix=sourceAssembly.jacobian;
                for(int k=0;k<matrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(matrix,k);it;++it)
                    rows[it.row()]+=std::abs(it.value()*columns[it.col()]);
                rows=rows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
                for(int k=0;k<matrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(matrix,k);it;++it)
                    it.valueRef()*=rows[it.row()]*columns[it.col()];
                experimental::ElectrothermalDirectSolver tangentSolver(cfg.value("electrothermal_linear_solver",experimental::ElectrothermalDirectSolver::defaultBackend()));
                const auto factorStart=Clock::now();++factorizations;
                tangentSolver.compute(matrix,false);factorizationSeconds+=seconds(factorStart);tangentAnalyses=tangentSolver.analyses();
                if(tangentSolver.info()!=Eigen::Success)throw std::runtime_error("Tangent factorization failed");
                const auto solveStart=Clock::now();const VectorXd rhs=-(derivative.array()*rows.array()).matrix();
                VectorXd direction=tangentSolver.solve(rhs);direction.array()*=columns.array();solveSeconds+=seconds(solveStart);
                if(tangentSolver.info()!=Eigen::Success || !direction.allFinite())throw std::runtime_error("Tangent linear solve failed");
                tangentPreparation["scaled_linear_residual"]=(matrix*(direction.array()/columns.array()).matrix()-rhs).norm();
                seed+=delta*direction;
                predictorCandidates.push_back({{"label","tangent"},{"potential_origin_V",cfg.value("potential_origin_V",0.)},
                    {"referenced_state_interleaved",list(seed)},{"electron_qf_reference_V",list(eReference)},{"hole_qf_reference_V",list(hReference)}});
                tangentPreparation["prepared"]=true;
            } catch(const std::exception& error) {tangentPreparation["fallback_reason"]=error.what();}
            tangentPreparation["wall_seconds"]=seconds(start);
        }
        // Opt-in predictor screening. Compare every candidate under the SAME
        // primary Jacobian row scale and physical QF reference representation.
        json predictorSelection;
        if(!predictorCandidates.empty()) {
            if(solveMode!="coupled" || initialization!="provided_state" || !cfg.value("use_qf_references",true))
                throw std::invalid_argument("Predictor screening requires a referenced coupled provided state");
            const auto& candidates=predictorCandidates;
            if(!candidates.is_array() || candidates.size()>4)
                throw std::invalid_argument("Predictor screening accepts at most four candidates");
            const auto screenStart=Clock::now();
            VectorXd scale=VectorXd::Zero(x.size());
            for(int k=0;k<a.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(a.jacobian,k);it;++it)
                scale[it.row()]+=std::abs(it.value()*(it.col()%4==3?10.:.025));
            scale=scale.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
            const auto quality=[&](const ElectrothermalAssembly& r) {
                const auto blocks=blockGates(r);
                std::array<Real,5> q{(r.residual.array()*scale.array()).matrix().norm(),rowGate(r).maxRatio,0.,0.,0.};
                for(int k=0;k<3;++k)q[k+2]=blocks[k].at("weighted_l2");
                return q;
            };
            const std::array<Real,5> floor{1e-9,1.,gateConfig.blockAbsoluteConvergence.psiResidualCeiling,
                gateConfig.blockAbsoluteConvergence.electronResidualCeiling,gateConfig.blockAbsoluteConvergence.holeResidualCeiling};
            auto best=quality(a);VectorXd selected=x;
            predictorSelection={{"selected","primary"},{"primary_quality",best},{"quality_order",{"common_scaled_l2","carrier_row_ratio","psi_weighted_l2","electron_weighted_l2","hole_weighted_l2"}},
                {"candidates",json::array()}};
            for(const auto& candidate:candidates) {
                json record={{"label",candidate.value("label",std::string("candidate"))},{"improved_best",false}};
                try {
                    if(record.at("label")=="primary")throw std::invalid_argument("Predictor candidate label primary is reserved");
                    VectorXd trial=values(candidate.at("referenced_state_interleaved"));
                    const VectorXd er=values(candidate.at("electron_qf_reference_V")),hr=values(candidate.at("hole_qf_reference_V"));
                    if(trial.size()!=x.size() || er.size()!=eReference.size() || hr.size()!=hReference.size() ||
                       !trial.allFinite() || !er.allFinite() || !hr.allFinite() ||
                       candidate.value("potential_origin_V",cfg.value("potential_origin_V",0.))!=cfg.value("potential_origin_V",0.))
                        throw std::invalid_argument("Invalid predictor state/reference");
                    for(Index i=0;i<mesh.numNodes();++i) {
                        if(!(trial[4*i+3]>50. && trial[4*i+3]<5000.))throw std::invalid_argument("Invalid predictor temperature");
                        trial[4*i+1]=static_cast<Real>((static_cast<long double>(er[i])-eReference[i])+trial[4*i+1]);
                        trial[4*i+2]=static_cast<Real>((static_cast<long double>(hr[i])-hReference[i])+trial[4*i+2]);
                    }
                    // Apply exactly the same finite-contact projection as the primary seed.
                    for(const auto& [i,contact]:bc.holeRecombination)
                        trial[4*i]=assembler.neutralPotential(i,bc.neutralContactBias_V.at(i),trial[4*i+3]).first;
                    const auto trialAssembly=assemble(trial,false);const auto q=quality(trialAssembly);
                    const bool better=experimental::improvesElectrothermalPredictor(q,best,floor);
                    record["quality"]=q;record["improved_best"]=better;
                    if(better){best=q;selected=std::move(trial);predictorSelection["selected"]=record["label"];}
                } catch(const std::exception& error){record["rejected_reason"]=error.what();}
                predictorSelection["candidates"].push_back(record);
            }
            if(predictorSelection.at("selected")!="primary") {x=std::move(selected);recenter();a=assemble(x);}
            predictorSelection["selected_quality"]=best;
            predictorSelection["wall_seconds"]=seconds(screenStart);
        }
        const bool deferJacobian=cfg.value("defer_recentered_candidate_jacobian",false);
        const auto willRecenter=[&](const VectorXd& state){
            if(!cfg.value("use_qf_references",true))return false;
            for(Index i=0;i<mesh.numNodes();++i)for(int k=1;k<=2;++k)
                if(std::abs(state[4*i+k])>1e-3)return true;
            return false;
        };

        const unsigned maximum=cfg.value("diagnostic_newton_max_iterations",0u);
        if(cfg.contains("diagnostic_density_update_iterations")) {
            const auto& count=cfg.at("diagnostic_density_update_iterations");
            if(!count.is_number_integer() || count<0 || count>std::numeric_limits<unsigned>::max())
                throw std::invalid_argument("Density update iterations must be a nonnegative unsigned-range integer");
        }
        unsigned densityIterations=cfg.value("diagnostic_density_update_iterations",0u);
        if(densityIterations && cfg.value("diagnostic_density_requires_primary_prediction",false) &&
           !predictorSelection.is_null() && predictorSelection.at("selected")!="primary") {
            densityIterations=0;
            predictorSelection["density_update_disabled_after_fallback"]=true;
        }
        const std::string projectionMode=cfg.value("diagnostic_density_projection",std::string("off"));
        if(projectionMode!="off" && projectionMode!="v1" && projectionMode!="v2")
            throw std::invalid_argument("Density projection must be off, v1 or v2");
        const bool localProjection=projectionMode!="off" && solveMode=="coupled";
        const bool naturalDamping=cfg.value("diagnostic_natural_damping",false);
        unsigned naturalSolves=0;Real naturalSolveSeconds=0.;
        const bool adaptiveJacobian=cfg.value("diagnostic_adaptive_jacobian",false);
        const bool pseudoTransient=cfg.value("diagnostic_pseudo_transient",false) && solveMode=="coupled";
        const std::string pseudoAcceptance=cfg.value("diagnostic_pseudo_acceptance",std::string("steady"));
        if(pseudoAcceptance!="steady" && pseudoAcceptance!="defect_ser" && pseudoAcceptance!="defect_model")
            throw std::invalid_argument("Unknown pseudo transient acceptance");
        if(solveMode=="coupled" && pseudoAcceptance!="steady" && !pseudoTransient)
            throw std::invalid_argument("Defect acceptance requires pseudo transient");
        const bool defectAcceptance=pseudoTransient && pseudoAcceptance!="steady";
        const Real pseudoScale=cfg.value("diagnostic_pseudo_time_scale",1.);
        if(!(pseudoScale>0. && std::isfinite(pseudoScale)))throw std::invalid_argument("Invalid pseudo time scale");
        if(pseudoTransient && (adaptiveJacobian || naturalDamping || cfg.value("diagnostic_ngmres_recovery",false)))
            throw std::invalid_argument("Pseudo transient must be isolated from other solver experiments");
        const VectorXd continuityArea=cfg.contains("recombination_area_m2")?values(cfg.at("recombination_area_m2")):siliconArea;
        Real pseudoTau=0.,initialPseudoTau=0.,massSeconds=0.;unsigned massAssemblies=0,pseudoRetries=0;
        VectorXd pseudoRows;json pseudoSteps=json::array();Real defectSeconds=0.;unsigned defectEvaluations=0;
        const bool nearSwitch=cfg.value("diagnostic_near_steady_qf_switch",false);
        const bool nearRebase=cfg.value("diagnostic_near_steady_qf_rebase",false);
        const bool contactConsistency=cfg.value("diagnostic_near_steady_contact_consistency",false);
        const bool localQfLimiter=cfg.value("diagnostic_local_qf_limiter",false);
        if(localQfLimiter && (nearSwitch || nearRebase || pseudoTransient || projectionMode!="off" || densityIterations ||
           naturalDamping || adaptiveJacobian || cfg.value("diagnostic_ngmres_recovery",false) ||
           solveMode!="coupled" || initialization!="provided_state"))
            throw std::invalid_argument("Local QF limiter requires isolated coupled provided-state QF Newton");
        if(contactConsistency && (nearSwitch || nearRebase || pseudoTransient || projectionMode!="off" || densityIterations ||
           naturalDamping || adaptiveJacobian || cfg.value("diagnostic_ngmres_recovery",false) ||
           solveMode!="coupled" || initialization!="provided_state" || !cfg.value("use_qf_references",true)))
            throw std::invalid_argument("Contact consistency requires isolated referenced coupled provided-state QF Newton");
        bool contactAttempted=false;Real contactSeconds=0.;json contactEvents=json::array();
        if(nearRebase && (nearSwitch || pseudoTransient || projectionMode!="off" || densityIterations ||
           naturalDamping || adaptiveJacobian || cfg.value("diagnostic_ngmres_recovery",false) ||
           solveMode!="coupled" || initialization!="provided_state" || !cfg.value("use_qf_references",true)))
            throw std::invalid_argument("Near-steady rebase requires isolated referenced coupled provided-state QF Newton");
        if(nearSwitch && (!pseudoTransient || projectionMode!="v1" || pseudoAcceptance!="defect_model" ||
           initialization!="provided_state" || !cfg.value("use_qf_references",true) ||
           cfg.value("diagnostic_pseudo_direction_audit",false)))
            throw std::invalid_argument("Near-steady switch requires referenced coupled provided-state PTC V1 defect_model without direction audit");
        bool nearSwitched=false;Real nearSwitchSeconds=0.;json nearSwitchEvents=json::array();
        const bool directionAudit=cfg.value("diagnostic_pseudo_direction_audit",false);
        if(directionAudit && (!pseudoTransient || maximum!=1))
            throw std::invalid_argument("Pseudo direction audit requires a single coupled pseudo iteration");
        json pseudoAudit=json::array();Real auditSeconds=0.;
        VectorXd factorRows;bool lagEligible=false,forceFreshJacobian=false;
        unsigned jacobianAge=0,laggedSolves=0,freshRetries=0;
        const bool iterationTrace=cfg.value("diagnostic_iteration_trace",false);
        Real traceSeconds=0.;
        std::unique_ptr<SiliconThermalPhysics> tracePhysics;
        if(iterationTrace)tracePhysics=std::make_unique<SiliconThermalPhysics>(siliconParameters);
        const auto traceGates=[&](const ElectrothermalAssembly& value) {
            const auto gate=rowGate(value);json violations=json::array();
            for(const auto& row:gate.violations)violations.push_back({{"node",row.nodeId},{"carrier",row.carrier},{"ratio",row.ratio}});
            return json{{"blocks",blockGates(value)},{"row",{{"satisfied",gate.satisfied},{"max_ratio",gate.maxRatio},
                {"eps_row",gate.epsRow},{"qualified_rows",gate.qualifiedRowCount},{"violations",violations}}}};
        };
        std::unique_ptr<SiliconThermalPhysics> densityPhysics;
        std::vector<SiliconThermalPhysics::DopingPreparation> densityDoping;
        Real densityEvaluationSeconds=0.;unsigned densityEvaluations=0;
        if((densityIterations || localProjection || pseudoTransient) && solveMode=="coupled") {
            densityPhysics=std::make_unique<SiliconThermalPhysics>(siliconParameters);
            densityDoping.reserve(mesh.numNodes());
            for(Index i=0;i<mesh.numNodes();++i)densityDoping.push_back(densityPhysics->prepareDoping(nd[i],na[i]));
        }
        const auto densityProperties=[&](Index i,const SiliconThermalState& state) {
            const auto start=Clock::now();++densityEvaluations;
            try {
                auto result=densityPhysics->evaluate(state,densityPhysics->prepareTemperature(state.temperature_K,densityDoping[i]));
                densityEvaluationSeconds+=seconds(start);return result;
            } catch(...) {densityEvaluationSeconds+=seconds(start);throw;}
        };
        const Real voltageUpdateLimit=cfg.value("diagnostic_voltage_update_limit_V",.2);
        if(!(voltageUpdateLimit>0. && std::isfinite(voltageUpdateLimit)))
            throw std::invalid_argument("Voltage update limit must be positive and finite");
        const bool reuseSymbolic=cfg.value("reuse_sparselu_symbolic",true);
        if(cfg.contains("diagnostic_stagnation_window") && !cfg.at("diagnostic_stagnation_window").is_number_unsigned())
            throw std::invalid_argument("Stagnation window must be a nonnegative integer");
        const auto linearBackend=cfg.value("electrothermal_linear_solver",experimental::ElectrothermalDirectSolver::defaultBackend());
        const bool retainLinear=context && reuseLinear && reuseSymbolic;
        if(context && (!retainLinear || (context->linear_ && context->linear_->backend()!=linearBackend)))
            context->clearLinearContext();
        const bool linearObjectReused=retainLinear && bool(context->linear_);
        auto linear=linearObjectReused?context->linear_:
            std::make_shared<experimental::ElectrothermalDirectSolver>(linearBackend);
        if(retainLinear)context->linear_=linear;
        auto& lu=*linear;
        const auto analysesBefore=lu.analyses();
        experimental::ElectrothermalStagnationWatch stagnation(cfg.value("diagnostic_stagnation_window",0u));
        const bool residualRecovery=cfg.value("diagnostic_ngmres_recovery",false);
        if(residualRecovery && solveMode!="coupled")throw std::invalid_argument("NGMRES recovery requires coupled equations");
        std::vector<VectorXd> recoveryStates,recoveryResiduals,recoveryEReferences,recoveryHReferences;
        json recoveryHistory=json::array();unsigned recoveryReassemblies=0,recoveryUpdates=0;
        json history=json::array();std::string stop="frozen_state";
        VectorXd columns(x.size());for(int i=0;i<x.size();++i)columns[i]=i%4==3?10.:.025;
        for(unsigned iteration=0;iteration<maximum;++iteration){
            const bool recentered=recenter();
            const bool lagged=adaptiveJacobian && lagEligible && !recentered && !forceFreshJacobian && jacobianAge<2;
            forceFreshJacobian=false;
            if(recentered || (!lagged && a.jacobian.rows()==0))a=assemble(x);
            VectorXd rows=VectorXd::Zero(x.size());
            for(int k=0;k<a.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(a.jacobian,k);it;++it)
                rows[it.row()]+=std::abs(it.value()*columns[it.col()]);
            rows=rows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
            if(lagged)rows=factorRows;
            const auto merit=[&](const ElectrothermalAssembly& r){return (r.residual.array()*rows.array()).matrix().norm();};
            Real before=merit(a);if(before<1e-9 && rowGate(a).satisfied && blocksSatisfied(a)){stop="diagnostic_scaled_residual";break;}
            if(contactConsistency && !contactAttempted && !bc.holeRecombination.empty() &&
               before<1e-9 && blocksSatisfied(a)) {
                // Terminal repair only: evaluate the same algebraic target at
                // the actual current T. Do not modify QFs, T, their references,
                // the Jacobian formula, or the ordinary Newton fallback.
                contactAttempted=true;const auto started=Clock::now();
                VectorXd candidate=x;Real maxChange=0.;unsigned changed=0;
                for(const auto& [i,contact]:bc.holeRecombination) {
                    candidate[4*i]=assembler.neutralPotential(i,bc.neutralContactBias_V.at(i),x[4*i+3]).first;
                    const Real delta=std::abs(candidate[4*i]-x[4*i]);
                    maxChange=std::max(maxChange,delta);changed+=delta!=0.;
                }
                json event={{"next_iteration",iteration+1},{"merit_before",before},
                    {"row_ratio_before",rowGate(a).maxRatio},{"changed_nodes",changed},
                    {"max_potential_change_V",maxChange},{"accepted",false}};
                // A near-steady algebraic repair must stay within the existing
                // representative-point state consistency budget (not a new gate).
                if(changed && candidate.allFinite() && maxChange<=1e-8) {
                    auto trial=assemble(candidate);VectorXd trialRows=VectorXd::Zero(x.size());
                    for(int k=0;k<trial.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(trial.jacobian,k);it;++it)
                        trialRows[it.row()]+=std::abs(it.value()*columns[it.col()]);
                    trialRows=trialRows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
                    const Real trialMerit=(trial.residual.array()*trialRows.array()).matrix().norm();
                    const auto gate=rowGate(trial);
                    const bool accepted=trialMerit<1e-9 && gate.satisfied && blocksSatisfied(trial);
                    event.update({{"merit_after",trialMerit},{"row_ratio_after",gate.maxRatio},
                        {"blocks_after",blockGates(trial)},{"accepted",accepted},{"reassembled",true}});
                    if(accepted){x=std::move(candidate);a=std::move(trial);rows=std::move(trialRows);before=trialMerit;}
                } else event.update({{"reassembled",false},{"reason",changed?"state_change_guard":"already_consistent"}});
                const Real elapsed=seconds(started);contactSeconds+=elapsed;event["seconds"]=elapsed;
                contactEvents.push_back(event);
                // Rejected candidates leave x/a/scales, references, stagnation
                // history and the remaining Newton iteration budget untouched.
                if(event.at("accepted").get<bool>()){stop="diagnostic_scaled_residual";break;}
            }
            if((nearSwitch || nearRebase) && !nearSwitched && before<1e-9 && blocksSatisfied(a)) {
                const auto started=Clock::now();const Real oldMerit=before;
                const auto oldGate=rowGate(a);
                const auto rebased=experimental::rebaseSmallElectrothermalQf(x,eReference,hReference);
                nearSwitched=true;
                // References are captured by assemble; always rebuild the full
                // residual/Jacobian and its scales before the first QF step.
                a=assemble(x);rows.setZero();
                for(int k=0;k<a.jacobian.outerSize();++k)for(SparseMatrixd::InnerIterator it(a.jacobian,k);it;++it)
                    rows[it.row()]+=std::abs(it.value()*columns[it.col()]);
                rows=rows.unaryExpr([](Real v){return 1./std::max(v,1e-100);});
                before=merit(a);
                // The standalone R7 experiment changes representation only;
                // retain its existing stagnation history and iteration budget.
                if(nearSwitch)stagnation=experimental::ElectrothermalStagnationWatch(cfg.value("diagnostic_stagnation_window",0u));
                const Real elapsed=seconds(started);nearSwitchSeconds+=elapsed;
                nearSwitchEvents.push_back({{"next_iteration",iteration+1},{"merit_before",oldMerit},{"merit_after",before},
                    {"row_ratio_before",oldGate.maxRatio},{"row_ratio_after",rowGate(a).maxRatio},
                    {"changed_references",rebased.changed},{"rounding_estimate_V",rebased.maxRoundingEstimate_V},
                    {"retired_tau_s",pseudoTau},{"seconds",elapsed}});
                if(before<1e-9 && rowGate(a).satisfied && blocksSatisfied(a)){stop="diagnostic_scaled_residual";break;}
            }
            const bool pseudoStep=pseudoTransient && !nearSwitched;
            const bool projectionStep=localProjection && !nearSwitched;
            const bool defectStep=defectAcceptance && !nearSwitched;
            if(residualRecovery) {
                if(recoveryStates.size()==4){
                    recoveryStates.erase(recoveryStates.begin());recoveryResiduals.erase(recoveryResiduals.begin());
                    recoveryEReferences.erase(recoveryEReferences.begin());recoveryHReferences.erase(recoveryHReferences.begin());
                }
                recoveryStates.push_back(x);recoveryResiduals.push_back(a.residual);
                recoveryEReferences.push_back(eReference);recoveryHReferences.push_back(hReference);
            }
            SparseMatrixd matrix=a.jacobian;
            Real fixedBefore=0.;
            VectorXd pseudoOldDensity;std::vector<bool> pseudoActive;
            if(defectStep){pseudoOldDensity=VectorXd::Zero(2*mesh.numNodes());pseudoActive.assign(2*mesh.numNodes(),false);}
            if(pseudoStep) {
                const auto started=Clock::now();
                if(pseudoRows.size()==0)pseudoRows=rows;
                fixedBefore=(a.residual.array()*pseudoRows.array()).matrix().norm();
                std::vector<Eigen::Triplet<Real>> entries;std::vector<Real> timeScales;
                for(Index i=0;i<mesh.numNodes();++i)if(siliconArea[i]>0.) {
                    const std::array<bool,2> active{
                        !bc.neutralContactBias_V.contains(i)&&!bc.electronQf_V.contains(i),
                        !bc.holeQf_V.contains(i)&&(!bc.neutralContactBias_V.contains(i)||bc.holeRecombination.contains(i))};
                    if(!active[0]&&!active[1])continue;
                    const SiliconThermalState state{x[4*i],x[4*i+1],x[4*i+2],x[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                    const auto properties=densityProperties(i,state);
                    const auto storage=experimental::electrothermalCarrierStorage(properties,continuityArea[i],active);
                    if(defectStep)for(int k=0;k<2;++k) {
                        pseudoActive[2*i+k]=active[k] && continuityArea[i]>0.;
                        pseudoOldDensity[2*i+k]=k==0?properties.electrons_m3.value:properties.holes_m3.value;
                    }
                    for(int k=1;k<=2;++k)if(active[k-1]) {
                        for(int j=0;j<4;++j)if(storage[k].derivative[j]!=0.)
                            entries.emplace_back(4*i+k,4*i+j,storage[k].derivative[j]);
                        const Real diagonal=std::abs(a.jacobian.coeff(4*i+k,4*i+k));
                        const Real time=std::abs(storage[k].derivative[k])/diagonal;
                        if(time>0. && std::isfinite(time))timeScales.push_back(time);
                    }
                }
                if(pseudoTau==0.) {
                    if(timeScales.empty())pseudoTau=1.; // No free storage rows: exact algebraic no-op.
                    else {std::sort(timeScales.begin(),timeScales.end());pseudoTau=pseudoScale*timeScales[timeScales.size()/2];}
                    if(!(pseudoTau>0. && std::isfinite(pseudoTau)))throw std::invalid_argument("Nonfinite initial pseudo time");
                    initialPseudoTau=pseudoTau;
                }
                SparseMatrixd mass(x.size(),x.size());mass.setFromTriplets(entries.begin(),entries.end());
                matrix+=mass/pseudoTau;
                ++massAssemblies;massSeconds+=seconds(started);
                pseudoSteps.push_back({{"iteration",iteration+1},{"tau_s",pseudoTau},{"mass_nonzeros",mass.nonZeros()},
                    {"fixed_residual_before",fixedBefore}});
                if(directionAudit) {
                    const auto auditStart=Clock::now();
                    for(int frozen=0;frozen<4;++frozen) {
                        const auto isolated=frozen==0?a:assembler.assemble(x,bc,eReference,hReference,true,false,frozen&1,frozen&2);
                        if(!(isolated.residual.array()==a.residual.array()).all())
                            throw std::runtime_error("Direction audit changed the steady residual");
                        for(bool regularized:{false,true}) {
                            SparseMatrixd operatorMatrix=isolated.jacobian;
                            if(regularized)operatorMatrix+=mass/pseudoTau;
                            for(int k=0;k<operatorMatrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(operatorMatrix,k);it;++it)
                                it.valueRef()*=rows[it.row()]*columns[it.col()];
                            experimental::ElectrothermalDirectSolver auditSolver(cfg.value("electrothermal_linear_solver",experimental::ElectrothermalDirectSolver::defaultBackend()));
                            auditSolver.compute(operatorMatrix,false);
                            json record={{"frozen_mobility_derivatives",bool(frozen&1)},{"frozen_recombination_derivatives",bool(frozen&2)},
                                {"regularized",regularized},{"steady_residual_exact",true}};
                            if(auditSolver.info()!=Eigen::Success){record["error"]="factorization_failed";pseudoAudit.push_back(record);continue;}
                            const VectorXd rhsAudit=-(a.residual.array()*rows.array()).matrix();
                            VectorXd d=auditSolver.solve(rhsAudit);
                            if(auditSolver.info()!=Eigen::Success || !d.allFinite()){record["error"]="solve_failed";pseudoAudit.push_back(record);continue;}
                            record["relative_linear_residual"]=(operatorMatrix*d-rhsAudit).norm()/std::max(rhsAudit.norm(),Real(1e-300));
                            d.array()*=columns.array();
                            const VectorXd r=(a.residual.array()*pseudoRows.array()).matrix();
                            const VectorXd jd=((a.jacobian*d).array()*pseudoRows.array()).matrix();
                            record["true_steady_merit_directional_slope"]=r.dot(jd)/std::max(r.squaredNorm(),Real(1e-300));
                            std::array<unsigned,2> negative{};std::array<Real,2> minimum{0.,0.};std::array<int,2> worst{-1,-1};
                            VectorXd old=VectorXd::Zero(2*mesh.numNodes()),relative=old;std::vector<bool> active(old.size(),false);
                            for(Index i=0;i<mesh.numNodes();++i)if(siliconArea[i]>0.) {
                                const SiliconThermalState state{x[4*i],x[4*i+1],x[4*i+2],x[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                                const auto p=densityProperties(i,state);
                                for(int k=0;k<2;++k) {
                                    const bool free=k==0?(!bc.neutralContactBias_V.contains(i)&&!bc.electronQf_V.contains(i)):
                                        (!bc.holeQf_V.contains(i)&&(!bc.neutralContactBias_V.contains(i)||bc.holeRecombination.contains(i)));
                                    const auto& q=k==0?p.electrons_m3:p.holes_m3;if(!free || q.value<=0.)continue;
                                    Real delta=0.;for(int j=0;j<4;++j)delta+=q.derivative[j]*d[4*i+j];
                                    const Real rel=delta/q.value;old[2*i+k]=q.value;relative[2*i+k]=rel;active[2*i+k]=true;
                                    if(rel<=-1.)++negative[k];if(rel<minimum[k]){minimum[k]=rel;worst[k]=static_cast<int>(i);}
                                }
                            }
                            record["nonpositive_density_targets"]=negative;record["minimum_relative_density_change"]=minimum;
                            record["worst_density_nodes"]=worst;
                            Real cap=1.;for(int j=0;j<x.size();++j)if((j%4==0 || j%4==3) && d[j]!=0.)
                                cap=std::min(cap,(j%4==3?30.:voltageUpdateLimit)/std::abs(d[j]));
                            record["projected_trials"]=json::array();
                            for(Real fraction:{1.,1e-4}) {
                                const Real alpha=cap*fraction;json trial={{"alpha",alpha}};
                                try {
                                    VectorXd candidate=x+alpha*d;
                                    for(Index i=0;i<mesh.numNodes();++i)if(active[2*i]||active[2*i+1]) {
                                        const SiliconThermalState state{candidate[4*i],x[4*i+1],x[4*i+2],candidate[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                                        const auto p=densityProperties(i,state);
                                        for(int k=0;k<2;++k)if(active[2*i+k])candidate[4*i+k+1]=experimental::electrothermalDensityQf(p,state,
                                            experimental::electrothermalProjectedDensity(old[2*i+k],relative[2*i+k],alpha),k==1);
                                    }
                                    const auto probe=assembler.assemble(candidate,bc,eReference,hReference,false);
                                    trial["true_fixed_merit_ratio"]=(probe.residual.array()*pseudoRows.array()).matrix().norm()/fixedBefore;
                                } catch(const std::exception& error){trial["error"]=error.what();}
                                record["projected_trials"].push_back(trial);
                            }
                            pseudoAudit.push_back(record);
                        }
                    }
                    auditSeconds+=seconds(auditStart);
                }
            }
            for(int k=0;k<matrix.outerSize();++k)for(SparseMatrixd::InnerIterator it(matrix,k);it;++it)
                it.valueRef()*=rows[it.row()]*columns[it.col()];
            if(!lagged) {
                const auto factorStart=Clock::now();lu.compute(matrix,reuseSymbolic);
                factorizationSeconds+=seconds(factorStart);++factorizations;
                if(adaptiveJacobian){factorRows=rows;jacobianAge=0;}
            } else {++jacobianAge;++laggedSolves;}
            if(lu.info()!=Eigen::Success){stop="factorization_failed";break;}
            const VectorXd rhs=-(a.residual.array()*rows.array()).matrix();
            const auto solveStart=Clock::now();
            VectorXd direction=lu.solve(rhs);direction.array()*=columns.array();
            solveSeconds+=seconds(solveStart);
            if(lu.info()!=Eigen::Success || !direction.allFinite()){stop="linear_solve_failed";break;}
            Real alpha=1.;int limiter=-1;VectorXd trialDirection;
            std::array<unsigned,2> locallyLimited{};
            if(localQfLimiter)trialDirection=direction;
            std::array<int,4> maxDirectionIndex{-1,-1,-1,-1};
            std::array<Real,4> maxDirection{};
            for(int i=0;i<x.size();++i)if(direction[i]!=0.){
                const Real cap=(i%4==3?30.:voltageUpdateLimit)/std::abs(direction[i]);
                if(localQfLimiter && (i%4==1 || i%4==2)) {
                    if(cap<1.){trialDirection[i]=std::copysign(voltageUpdateLimit,direction[i]);++locallyLimited[i%4-1];}
                } else {
                    if(cap<alpha)limiter=i;
                    alpha=std::min(alpha,cap);
                }
                if(std::abs(direction[i])>maxDirection[i%4])maxDirectionIndex[i%4]=i;
                maxDirection[i%4]=std::max(maxDirection[i%4],std::abs(direction[i]));
            }
            const Real globalAlpha=alpha;
            bool densityAttempt=false,densityFallback=false;std::string densitySkipReason;
            int densityLimiter=-1;std::string densityLimiterKind="none";
            VectorXd densities,relativeDensityChange;
            std::array<unsigned,2> nonpositiveTargets{};
            std::vector<bool> densityActive;
            if((projectionStep || (!nearSwitched && iteration<densityIterations && before>1e-6)) && solveMode=="coupled") {
                densities=VectorXd::Zero(2*mesh.numNodes());relativeDensityChange=densities;
                densityActive.assign(2*mesh.numNodes(),false);
                Real densityAlpha=1.;bool valid=true;Real maximumActiveQfVt=0.;
                for(Index i=0;i<mesh.numNodes();++i)if(siliconArea[i]>0.) {
                    const SiliconThermalState state{x[4*i],x[4*i+1],x[4*i+2],x[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                    const auto properties=densityProperties(i,state);
                    for(int k=0;k<2;++k) {
                        const bool constrained=k==0?(bc.neutralContactBias_V.contains(i)||bc.electronQf_V.contains(i)):
                            (bc.holeQf_V.contains(i)||(bc.neutralContactBias_V.contains(i)&&!bc.holeRecombination.contains(i)));
                        const auto& quantity=k==0?properties.electrons_m3:properties.holes_m3;
                        if(constrained)continue;
                        if(quantity.value<=0.) {
                            if(projectionStep){valid=false;densitySkipReason="nonpositive_old_density";}
                            continue;
                        }
                        Real change=0.;for(int j=0;j<4;++j)change+=quantity.derivative[j]*direction[4*i+j];
                        const Real relative=change/quantity.value;
                        if(relative<=-1.)++nonpositiveTargets[k];
                        if(!std::isfinite(relative)){valid=false;densitySkipReason="nonfinite_relative_density_change";continue;}
                        densities[2*i+k]=quantity.value;relativeDensityChange[2*i+k]=relative;densityActive[2*i+k]=true;
                        maximumActiveQfVt=std::max(maximumActiveQfVt,std::abs(direction[4*i+k+1])/(constants::kb/constants::q*x[4*i+3]));
                        if(relative<0. && !projectionStep) {
                            const Real bound=-.99/relative;
                            if(bound<densityAlpha){densityLimiter=4*i+k+1;densityLimiterKind="density_positivity";}
                            densityAlpha=std::min(densityAlpha,bound);
                        }
                    }
                }
                for(int i=0;i<x.size();++i)if(direction[i]!=0.) {
                    const int k=i%4;
                    if((k==1||k==2)&&(projectionStep || densityActive[2*(i/4)+k-1]))continue;
                    const Real bound=(k==3?30.:voltageUpdateLimit)/std::abs(direction[i]);
                    if(bound<densityAlpha){densityLimiter=i;densityLimiterKind="component_cap";}
                    densityAlpha=std::min(densityAlpha,bound);
                }
                if(projectionMode=="v2" && maximumActiveQfVt<.01){valid=false;densitySkipReason="near_qf_direction";}
                if(valid && densityAlpha>0. && std::isfinite(densityAlpha) && std::any_of(densityActive.begin(),densityActive.end(),[](bool v){return v;})) {alpha=densityAlpha;densityAttempt=true;}
            }
            json trace;
            if(iterationTrace) {
                const auto started=Clock::now();
                trace={{"schema","vela.newton_iteration_trace.v1"},{"before_gates",traceGates(a)},
                    {"direction_maxima",json::array()},{"carrier_samples",json::array()},
                    {"density_limiter_node",densityLimiter<0?-1:densityLimiter/4},
                    {"density_limiter_component",densityLimiter<0?-1:densityLimiter%4},
                    {"density_limiter_kind",densityLimiterKind}};
                for(int k=0;k<4;++k) {
                    const int index=maxDirectionIndex[k];const int node=index<0?-1:index/4;
                    trace["direction_maxima"].push_back({{"component",k},{"node",node},
                        {"signed_direction",index<0?0.:direction[index]},
                        {"temperature_K",node<0?json(nullptr):json(x[4*node+3])}});
                }
                for(int k=1;k<=2;++k) {
                    const int index=maxDirectionIndex[k];if(index<0)continue;const Index i=index/4;
                    json sample={{"node",i},{"component",k},{"silicon",siliconArea[i]>0.}};
                    if(siliconArea[i]>0.)try {
                        const SiliconThermalState state{x[4*i],x[4*i+1],x[4*i+2],x[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                        const auto properties=tracePhysics->evaluate(state);
                        const auto& quantity=k==1?properties.electrons_m3:properties.holes_m3;
                        std::array<Real,4> parts{};Real relative=0.;
                        for(int j=0;j<4;++j){parts[j]=quantity.derivative[j]*direction[4*i+j]/quantity.value;relative+=parts[j];}
                        const Real vt=constants::kb*state.temperature_K/constants::q;
                        sample.update({{"density_before_m3",quantity.value},{"relative_linear_change",relative},
                            {"relative_change_by_component",parts},{"local_Vt_V",vt},{"qf_direction_over_Vt",direction[index]/vt}});
                    } catch(const std::exception& error){sample["error"]=error.what();}
                    trace["carrier_samples"].push_back(sample);
                }
                traceSeconds+=seconds(started);
            }
            const Real initialAlpha=alpha;
            unsigned trials=0;json projectionTrials=json::array(),naturalTrials=json::array(),trialFailures=json::array();
            bool accepted=false;Real acceptedDefectNorm=0.,acceptedModelError=0.;
            bool usedDefect=false;json defectTrials=json::array();
            for(int trial=0;trial<(densityAttempt?32:24);++trial){
                if(densityAttempt && trial==8){alpha=globalAlpha;densityFallback=true;}
                ++trials;Real nextAlpha=.5*alpha;
                VectorXd candidate=x+alpha*(localQfLimiter?trialDirection:direction);
                const auto trialStart=iterationTrace?Clock::now():Clock::time_point{};
                bool candidateJacobian=false;
                // Such a candidate must be reassembled after reference changes
                // before the next Newton solve. Its merit needs residuals only.
                try{
                    unsigned projectedElectrons=0,projectedHoles=0;
                    std::array<Real,2> projectionCharge{};Real maximumViolation=0.;
                    Real poissonResidualL1=0.;
                    if(projectionStep)for(Index i=0;i<mesh.numNodes();++i)
                        if(!bc.neutralContactBias_V.contains(i)&&!bc.potential_V.contains(i))poissonResidualL1+=std::abs(a.residual[4*i]);
                    if(densityAttempt && !densityFallback)for(Index i=0;i<mesh.numNodes();++i) {
                        if(!densityActive[2*i]&&!densityActive[2*i+1])continue;
                        const SiliconThermalState state{candidate[4*i],x[4*i+1],x[4*i+2],candidate[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                        const auto properties=densityProperties(i,state);
                        for(int k=0;k<2;++k)if(densityActive[2*i+k]) {
                            const Real raw=densities[2*i+k]*(1.+alpha*relativeDensityChange[2*i+k]);
                            const Real target=projectionStep?experimental::electrothermalProjectedDensity(densities[2*i+k],relativeDensityChange[2*i+k],alpha):raw;
                            if(target!=raw){
                                if(k==0)++projectedElectrons;else ++projectedHoles;
                                projectionCharge[k]+=constants::q*siliconArea[i]*(target-raw);
                                maximumViolation=std::max(maximumViolation,(target-raw)/densities[2*i+k]);
                            }
                            candidate[4*i+k+1]=experimental::electrothermalDensityQf(properties,state,target,k==1);
                        }
                    }
                    if(projectionStep)projectionTrials.push_back({{"alpha",alpha},{"fallback",densityFallback},
                        {"projected_electrons",projectedElectrons},{"projected_holes",projectedHoles},
                        {"max_relative_correction",maximumViolation},{"charge_correction_abs_C_per_m",projectionCharge},
                        {"poisson_residual_l1_C_per_m",poissonResidualL1},
                        {"charge_ratio",poissonResidualL1>0.?json((projectionCharge[0]+projectionCharge[1])/poissonResidualL1):json(nullptr)}});
                    candidateJacobian=!adaptiveJacobian && !(deferJacobian && willRecenter(candidate));
                    auto next=assemble(candidate,candidateJacobian);
                    bool acceptable=merit(next)<before || (before<1e-9 && merit(next)<1e-9 && rowGate(next).maxRatio<rowGate(a).maxRatio);
                    if(defectStep && before>=1e-9) {
                        const auto started=Clock::now();++defectEvaluations;
                        VectorXd defect=next.residual;
                        try {
                            for(Index i=0;i<mesh.numNodes();++i)if(pseudoActive[2*i]||pseudoActive[2*i+1]) {
                                const SiliconThermalState state{candidate[4*i],candidate[4*i+1],candidate[4*i+2],candidate[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                                const auto properties=densityProperties(i,state);
                                for(int k=0;k<2;++k)if(pseudoActive[2*i+k])
                                    defect[4*i+k+1]+=experimental::electrothermalStorageDifference(pseudoOldDensity[2*i+k],
                                        k==0?properties.electrons_m3.value:properties.holes_m3.value,continuityArea[i],k==1)/pseudoTau;
                            }
                            const auto check=experimental::electrothermalPseudoTrial(
                                (a.residual.array()*pseudoRows.array()).matrix(),(defect.array()*pseudoRows.array()).matrix(),alpha);
                            acceptable=check.accepted;
                            defectTrials.push_back({{"alpha",alpha},{"accepted",acceptable},{"defect_norm",check.norm},
                                {"model_error",check.modelError},{"steady_fixed_norm",(next.residual.array()*pseudoRows.array()).matrix().norm()}});
                            if(acceptable){acceptedDefectNorm=check.norm;acceptedModelError=check.modelError;usedDefect=true;}
                        } catch(...) {defectSeconds+=seconds(started);throw;}
                        defectSeconds+=seconds(started);
                    }
                    // Keep the original near-floor rule and all final gates.
                    // Freeze row/column scaling and the factorization for the
                    // corrector solve: no candidate-dependent residual scaling.
                    if(naturalDamping && before>=1e-9) {
                        const auto started=Clock::now();
                        const VectorXd correction=lu.solve(-(next.residual.array()*rows.array()).matrix());
                        const Real elapsed=seconds(started);naturalSolveSeconds+=elapsed;solveSeconds+=elapsed;++naturalSolves;
                        if(lu.info()!=Eigen::Success)throw std::runtime_error("Natural corrector solve failed");
                        const auto check=experimental::electrothermalNaturalTrial(
                            (direction.array()/columns.array()).matrix(),correction,alpha,projectedElectrons+projectedHoles>0);
                        acceptable=check.decreasing;nextAlpha=check.nextAlpha;
                        naturalTrials.push_back({{"alpha",alpha},{"theta",check.theta},{"accepted",acceptable},
                            {"merit",merit(next)},{"projected",projectedElectrons+projectedHoles>0}});
                    }
                    if(iterationTrace)trace["line_search_candidates"].push_back({{"trial",trial+1},
                        {"alpha",alpha},{"merit",merit(next)},{"accepted",acceptable},
                        {"built_jacobian",candidateJacobian},{"trial_seconds",seconds(trialStart)}});
                    if(acceptable){x=candidate;a=std::move(next);accepted=true;break;}}
                catch(const std::exception& error){
                    if(iterationTrace)trace["line_search_candidates"].push_back({{"trial",trial+1},
                        {"alpha",alpha},{"accepted",false},{"exception",error.what()},
                        {"built_jacobian",candidateJacobian},{"trial_seconds",seconds(trialStart)}});
                    if(projectionStep || naturalDamping || defectStep)trialFailures.push_back({{"trial",trial+1},{"alpha",alpha},{"reason",error.what()}});
                }
                alpha=nextAlpha;
            }
            history.push_back({{"iteration",iteration+1},{"scaled_l2_before",before},{"scaled_l2_after",merit(a)},{"alpha",alpha},{"accepted",accepted}});
            if(localQfLimiter)history.back()["local_qf_limiter"]={{"clipped_electrons",locallyLimited[0]},
                {"clipped_holes",locallyLimited[1]},{"cap_V",voltageUpdateLimit}};
            if(pseudoStep) {
                const Real after=(a.residual.array()*pseudoRows.array()).matrix().norm();
                pseudoSteps.back().update({{"fixed_residual_after",after},{"accepted",accepted},{"alpha",alpha},
                    {"nonpositive_density_targets",densityAttempt?json(nonpositiveTargets):json(nullptr)}});
                // SER uses the INITIAL steady residual scaling for the entire
                // solve. Never interpret changing row equilibration as progress.
                Real factor=accepted?std::clamp(fixedBefore/std::max(after,Real(1e-300)),.5,4.):.1;
                if(usedDefect && pseudoAcceptance=="defect_model")
                    factor=acceptedModelError<=.25?2.:acceptedModelError>.75?.5:1.;
                if(defectStep)pseudoSteps.back().update({{"acceptance",pseudoAcceptance},{"defect_used",usedDefect},
                    {"accepted_defect_norm",usedDefect?json(acceptedDefectNorm):json(nullptr)},
                    {"accepted_model_error",usedDefect?json(acceptedModelError):json(nullptr)},
                    {"tau_factor",factor},{"defect_trials",defectTrials}});
                pseudoTau=std::clamp(pseudoTau*factor,initialPseudoTau*1e-6,initialPseudoTau*1e12);
                history.back()["pseudo_transient"]=pseudoSteps.back();
            }
            if(nearSwitch || nearRebase)history.back()["near_steady_qf_active"]=nearSwitched;
            if(profiling){
                history.back()["initial_alpha"]=initialAlpha;
                history.back()["line_search_trials"]=trials;
                history.back()["limiter_node"]=limiter<0?-1:limiter/4;
                history.back()["limiter_component"]=limiter<0?-1:limiter%4;
                history.back()["limiter_direction"]=limiter<0?0.:direction[limiter];
                history.back()["max_abs_direction_by_block"]=maxDirection;
                if(densityIterations || projectionStep) {
                    history.back()["density_update_attempted"]=densityAttempt;
                    history.back()["density_update_fallback"]=densityFallback;
                    history.back()["global_initial_alpha"]=globalAlpha;
                }
            }
            if(projectionStep || naturalDamping || defectStep)history.back()["trial_failures"]=trialFailures;
            if(adaptiveJacobian){history.back()["lagged_jacobian"]=lagged;history.back()["jacobian_age"]=jacobianAge;}
            if(naturalDamping)history.back()["natural_trials"]=naturalTrials;
            if(projectionStep){history.back()["projection_trials"]=projectionTrials;
                history.back()["density_skip_reason"]=densitySkipReason;}
            if(iterationTrace) {
                const auto started=Clock::now();trace["after_gates"]=traceGates(a);
                for(auto& sample:trace["carrier_samples"])if(sample.contains("density_before_m3"))try {
                    const Index i=sample.at("node");const int k=sample.at("component");
                    const SiliconThermalState state{x[4*i],x[4*i+1],x[4*i+2],x[4*i+3],nd[i],na[i],eReference[i],hReference[i]};
                    const auto properties=tracePhysics->evaluate(state);
                    const Real value=k==1?properties.electrons_m3.value:properties.holes_m3.value;
                    sample["density_after_m3"]=value;
                    sample["actual_density_ratio"]=value/sample.at("density_before_m3").get<Real>();
                } catch(const std::exception& error){sample["after_error"]=error.what();}
                trace["initial_alpha"]=initialAlpha;trace["line_search_trials"]=trials;
                trace["raw_limiter_node"]=limiter<0?-1:limiter/4;trace["raw_limiter_component"]=limiter<0?-1:limiter%4;
                history.back()["iteration_trace"]=std::move(trace);traceSeconds+=seconds(started);
            }
            progress<<history.back().dump()<<std::endl;
            const bool stalled=stagnation.update(usedDefect?fixedBefore:before,usedDefect?acceptedDefectNorm:merit(a),alpha,accepted);
            if(pseudoStep && !accepted && pseudoRetries<3) {
                ++pseudoRetries;
                stagnation=experimental::ElectrothermalStagnationWatch(cfg.value("diagnostic_stagnation_window",0u));
                stop="diagnostic_iteration_limit";continue;
            }
            if(adaptiveJacobian) {
                lagEligible=accepted && alpha==1. && merit(a)<=.1*before;
                if(lagged && (!accepted || stalled)) {
                    forceFreshJacobian=true;lagEligible=false;++freshRetries;
                    stagnation=experimental::ElectrothermalStagnationWatch(cfg.value("diagnostic_stagnation_window",0u));
                    stop="diagnostic_iteration_limit";continue;
                }
            }
            if((!accepted || stalled) && residualRecovery && recoveryHistory.empty()) {
                const auto recoveryStart=Clock::now();json record={{"iteration",iteration+1},{"accepted",false},{"candidates",json::array()}};
                try {
                    // Tail steps can recenter every iteration. Convert each old
                    // referenced state and re-evaluate its residual under the
                    // CURRENT references; never mix residual coordinate epochs.
                    auto samples=recoveryStates,residuals=recoveryResiduals;
                    for(std::size_t j=0;j<samples.size();++j) {
                        if((recoveryEReferences[j].array()==eReference.array()).all() &&
                           (recoveryHReferences[j].array()==hReference.array()).all())continue;
                        for(Index i=0;i<mesh.numNodes();++i) {
                            samples[j][4*i+1]=static_cast<Real>((static_cast<long double>(recoveryEReferences[j][i])-eReference[i])+samples[j][4*i+1]);
                            samples[j][4*i+2]=static_cast<Real>((static_cast<long double>(recoveryHReferences[j][i])-hReference[i])+samples[j][4*i+2]);
                        }
                        residuals[j]=assemble(samples[j],false).residual;++recoveryReassemblies;
                    }
                    const auto mix=experimental::electrothermalResidualMix(x,a.residual,samples,residuals,rows);
                    record["rank"]=mix.rank;record["coefficient_l1"]=mix.coefficientL1;record["predicted_norm"]=mix.predictedNorm;
                    Real step=1.;for(int i=0;i<x.size();++i)if(mix.direction[i]!=0.)step=std::min(step,(i%4==3?30.:voltageUpdateLimit)/std::abs(mix.direction[i]));
                    const auto quality=[&](const ElectrothermalAssembly& value) {
                        const auto blocks=blockGates(value);std::array<Real,5> q{merit(value),rowGate(value).maxRatio,0.,0.,0.};
                        for(int k=0;k<3;++k)q[k+2]=blocks[k].at("weighted_l2");return q;
                    };
                    const auto baseline=quality(a);
                    const std::array<Real,5> floor{1e-9,1.,gateConfig.blockAbsoluteConvergence.psiResidualCeiling,
                        gateConfig.blockAbsoluteConvergence.electronResidualCeiling,gateConfig.blockAbsoluteConvergence.holeResidualCeiling};
                    record["baseline_quality"]=baseline;
                    for(unsigned trial=0;trial<4;++trial,step*=.5) {
                        json check={{"step",step},{"accepted",false}};
                        try {
                            VectorXd candidate=x+step*mix.direction;
                            if(!candidate.allFinite())throw std::invalid_argument("Nonfinite recovery state");
                            for(Index i=0;i<mesh.numNodes();++i)if(!(candidate[4*i+3]>50. && candidate[4*i+3]<5000.))
                                throw std::invalid_argument("Invalid recovery temperature");
                            const auto next=assemble(candidate,false);const auto q=quality(next);check["quality"]=q;
                            if(q[0]<.99*baseline[0] && experimental::improvesElectrothermalPredictor(q,baseline,floor)) {
                                // Commit only after the full Jacobian is available;
                                // failed preparation must preserve the old state.
                                auto full=assemble(candidate);x=std::move(candidate);a=std::move(full);
                                check["accepted"]=true;record["accepted"]=true;++recoveryUpdates;
                            }
                        } catch(const std::exception& error){check["rejected_reason"]=error.what();}
                        record["candidates"].push_back(check);if(record.at("accepted").get<bool>())break;
                    }
                } catch(const std::exception& error){record["fallback_reason"]=error.what();}
                record["wall_seconds"]=seconds(recoveryStart);recoveryHistory.push_back(record);
                progress<<json{{"ngmres_recovery",record}}.dump()<<std::endl;
                if(record.at("accepted").get<bool>()) {
                    stagnation=experimental::ElectrothermalStagnationWatch(cfg.value("diagnostic_stagnation_window",0u));
                    recoveryStates.clear();recoveryResiduals.clear();recoveryEReferences.clear();recoveryHReferences.clear();
                    stop="diagnostic_iteration_limit";continue;
                }
            }
            if(!accepted){stop="line_search_failed";break;}
            if(stalled){stop="diagnostic_stagnation_reject";break;}
            stop="diagnostic_iteration_limit";
        }

        if(a.jacobian.rows()==0)a=assemble(x);
        const auto rows=rowGate(a);
        json holeAudit=json::array();
        if(cfg.contains("diagnostic_hole_row_audit_nodes")) {
            std::vector<ElectrothermalHoleRowAudit> selected;
            for(Index node:cfg.at("diagnostic_hole_row_audit_nodes").get<std::vector<Index>>()){
                ElectrothermalHoleRowAudit row;row.node=node;selected.push_back(row);
            }
            const auto checked=assembler.assemble(x,bc,eReference,hReference,false,false,false,false,&selected);
            if(!(checked.residual.array()==a.residual.array()).all() ||
               !(checked.holeFluxAbs_A_per_m.array()==a.holeFluxAbs_A_per_m.array()).all())
                throw std::runtime_error("Read-only row audit changed assembly");
            const auto stateJson=[](const SiliconThermalState& s){return json{{"psi",s.potential_V},{"fn",s.electronQf_V},
                {"fp",s.holeQf_V},{"T",s.temperature_K},{"rn",s.electronQfReference_V},{"rp",s.holeQfReference_V}};};
            const auto propertyJson=[](const SiliconThermalResult& p){return json{{"n",p.electrons_m3.value},{"p",p.holes_m3.value},
                {"Nc",p.Nc_m3.value},{"Nv",p.Nv_m3.value},{"en",p.electronEta.value},{"ep",p.holeEta.value},
                {"ni",p.effectiveNi_m3.value},{"tn",p.electronLifetime_s.value},{"tp",p.holeLifetime_s.value},
                {"Cn",p.augerElectron_m6_per_s.value},{"Cp",p.augerHole_m6_per_s.value},
                {"srh",p.srhRate_m3_per_s.value},{"auger",p.augerRate_m3_per_s.value}};};
            for(const auto& row:selected){
                json edges=json::array();
                for(const auto& e:row.edges){
                    const Real lf=(std::log(e.pb.holes_m3.value)-std::log(e.pb.Nv_m3.value))-
                        (std::log(e.pa.holes_m3.value)-std::log(e.pa.Nv_m3.value));
                    const Real mid=.5*(e.pa.holeEta.value+e.pb.holeEta.value);
                    const auto limit=[&](){return std::max(1.,fermiDiracHalf(mid)/fermiDiracHalfDerivative(mid));};
                    Real g=std::abs(lf)>1e-8?(e.pb.holeEta.value-e.pa.holeEta.value)/lf:limit();
                    if(!(g>0.) || !std::isfinite(g))g=limit();
                    edges.push_back({{"id",e.id},{"a",e.a},{"b",e.b},{"sa",stateJson(e.sa)},{"sb",stateJson(e.sb)},
                        {"pa",propertyJson(e.pa)},{"pb",propertyJson(e.pb)},{"mu",e.mobility},{"weight",e.weight},
                        {"current",e.current},{"logF",lf},{"g",g}});
                }
                holeAudit.push_back({{"node",row.node},{"state",stateJson(row.state)},{"properties",propertyJson(row.properties)},
                    {"q_area",row.sourceAreaCharge},{"source",row.source},{"residual",row.residual},{"flux_abs",row.fluxAbs},
                    {"edges",edges},{"current_scale",currentScale},{"q",constants::q},{"kb",constants::kb}});
                if(row.finiteContact)holeAudit.back()["contact"]={{"bias",row.contactBias},{"neutral_potential",row.neutralPotential},
                    {"neutral_dT",row.neutralTemperatureDerivative},{"vt",row.contactVt},{"delta",row.contactDelta},
                    {"coefficient",row.contactCoefficient},{"outward",row.contactOutward},{"Nv",row.equilibriumNv},
                    {"eta",row.equilibriumEta},{"p0",row.equilibriumDensity},{"df",row.fermiDerivative},{"ddf",row.fermiSecondDerivative}};
            }
        }
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
            {"scope","Point solver stages including enabled predictor/recovery work; factorization includes symbolic analysis; nested preparation timers overlap these stages"},
            {"assembly_seconds",assemblySeconds},{"assembly_calls",assemblyCalls},
            {"full_assembly_seconds",fullAssemblySeconds},{"residual_assembly_seconds",residualAssemblySeconds},
            {"residual_only_calls",residualOnlyCalls},
            {"factorization_seconds",factorizationSeconds},{"factorizations",factorizations},
            {"linear_solve_seconds",solveSeconds},
            {"fermi_half_calls",physicsCallCounters.fermiDiracHalf-countersBefore.fermiDiracHalf},
            {"fermi_half_derivative_calls",physicsCallCounters.fermiDiracHalfDerivative-countersBefore.fermiDiracHalfDerivative},
            {"inverse_fermi_half_calls",physicsCallCounters.inverseFermiDiracHalf-countersBefore.inverseFermiDiracHalf}};
        if(nearSwitch)result["near_steady_qf_switch"]={{"triggered",nearSwitched},{"events",nearSwitchEvents},
            {"seconds",nearSwitchSeconds},{"merit_trigger",1e-9},{"small_qf_limit_V",1e-3},
            {"scope","One-way representation change and retirement of mass/density updates; original steady gates"}};
        if(nearRebase)result["near_steady_qf_rebase"]={{"triggered",nearSwitched},{"events",nearSwitchEvents},
            {"seconds",nearSwitchSeconds},{"merit_trigger",1e-9},{"small_qf_limit_V",1e-3},
            {"scope","One-time QF representation only; original gates, stagnation watch and iteration budget"}};
        if(contactConsistency)result["near_steady_contact_consistency"]={{"triggered",contactAttempted},
            {"events",contactEvents},{"seconds",contactSeconds},{"merit_trigger",1e-9},{"maximum_change_V",1e-8},
            {"scope","One-time finite-contact potential repair at current T; accept only original steady gates, otherwise unchanged Newton fallback; inclusive preparation/reassembly time"}};
        if(iterationTrace)result["iteration_trace_seconds"]=traceSeconds;
        if(cfg.contains("diagnostic_hole_row_audit_nodes"))result["hole_row_audit"]=holeAudit;
        if(localProjection)result["density_projection"]={{"variant",projectionMode},{"relative_floor",.01},
            {"absolute_floor_m3",1e-250},{"floor_capped_by_old_density",true},{"scope","All free silicon carrier unknowns; constrained rows preserved"}};
        if(pseudoTransient)result["pseudo_transient"]={{"initial_tau_s",initialPseudoTau},{"next_tau_s",pseudoTau},
            {"mass_assemblies",massAssemblies},{"mass_seconds",massSeconds},{"failed_step_retries",pseudoRetries},
            {"time_scale",pseudoScale},{"steps",pseudoSteps},
            {"scope","Carrier storage only; continuity source areas; fixed initial residual scaling for SER; original steady gates"}};
        if(defectAcceptance)result["pseudo_transient"].update({{"acceptance",pseudoAcceptance},
            {"defect_evaluations",defectEvaluations},{"defect_seconds",defectSeconds},
            {"scope","Actual backward-Euler defect Armijo trial; original near-floor rule and steady terminal gates; not transient qualification"}});
        if(directionAudit)result["pseudo_direction_audit"]={{"seconds",auditSeconds},{"directions",pseudoAudit},
            {"scope","Read-only first-state directions; omitted coefficient derivatives are never used for state updates"}};
        if(adaptiveJacobian)result["adaptive_jacobian"]={{"lagged_solves",laggedSolves},{"fresh_retries",freshRetries},
            {"maximum_lag",2},{"scope","Reuse only after an accepted full step with residual ratio <=0.1; recentering forces refresh"}};
        if(naturalDamping)result["natural_damping"]={{"corrector_solves",naturalSolves},{"corrector_solve_seconds",naturalSolveSeconds},
            {"scope","NLEQ_ERR-type current-Jacobian corrector test, fixed scaling; original terminal gates retained"}};
        if(profiling){
            auto& perf=result["performance"];
            perf["point_preparation_seconds"]=pointPreparationSeconds;
            perf["preparation_key_seconds"]=preparationKeySeconds;
            perf["preparation_mesh_seconds"]=preparationHit?0.:prepared->meshSeconds;
            perf["preparation_input_seconds"]=preparationHit?0.:prepared->inputSeconds;
            perf["preparation_ialmob_geometry_seconds"]=preparationHit?0.:prepared->geometrySeconds;
            perf["preparation_assembler_seconds"]=assemblerPreparationSeconds;
            perf["static_preparation_reused"]=preparationHit;
            perf["jacobian_structure_reuse_enabled"]=reuseStructure;
            if(reuseStructure){
                perf["jacobian_structure_builds"]=prepared->structure->coupled->builds-structureBefore[0];
                perf["jacobian_structure_hits"]=prepared->structure->coupled->hits-structureBefore[1];
                perf["heat_structure_builds"]=prepared->structure->heat->builds-structureBefore[2];
                perf["heat_structure_hits"]=prepared->structure->heat->hits-structureBefore[3];
                perf["jacobian_structure_nonzeros"]=prepared->structure->coupled->pattern.nonZeros();
            }
            perf["ialmob_pass_calls"]=ialKernelProfile.passes;
            perf["ialmob_pass_seconds"]=ialKernelProfile.seconds;
            perf["ialmob_high_field_evaluations"]=ialKernelProfile.highFieldEvaluations;
            perf["ialmob_high_field_reuses"]=ialKernelProfile.highFieldReuses;
            perf["ialmob_screening_cache_requests"]=ialKernelProfile.screeningRequests;
            perf["ialmob_screening_cache_hits"]=ialKernelProfile.screeningHits;
            perf["ialmob_local_preparation_hits"]=ialKernelProfile.localPreparationHits;
            perf["ialmob_local_preparation_builds"]=ialKernelProfile.localPreparationBuilds;
            perf["ialmob_screening_method"]=ialScreeningMethodName(ial_json::electrothermalScreeningMethod(cfg));
            perf["ialmob_screening_candidate_calls"]=ialKernelProfile.screeningCandidateCalls;
            perf["ialmob_screening_function_evaluations"]=ialKernelProfile.screeningFunctionEvaluations;
            perf["ialmob_screening_fallbacks"]=ialKernelProfile.screeningFallbacks;
        }
        if(electrothermalCostProfile.mode!=ElectrothermalCostProfile::Off){
            auto& cost=result["subcost_diagnostics"];
            cost["mode"]=cfg.at("diagnostic_electrothermal_cost");
            const std::array<const char*,ElectrothermalCostProfile::Count> names{
                "residual_heat_total","residual_heat_matrix_prepare","residual_heat_matrix_fill",
                "residual_heat_matrix_finish","screening_root","solver_create","symbolic_analysis","numeric_factorization"};
            for(std::size_t i=0;i<names.size();++i)cost[names[i]]={
                {"seconds",electrothermalCostProfile.seconds[i]},{"calls",electrothermalCostProfile.calls[i]}};
        }
        if(profiling)result["performance"]["preparation_counts"]=assembler.preparationCounts();
        if(profiling)result["performance"]["neutral_root_counts"]=assembler.neutralRootCounts();
        if(profiling)result["performance"]["neutral_root_iteration_counts"]=assembler.neutralRootIterationCounts();
        if(profiling)result["performance"]["symbolic_analyses"]=lu.analyses()-analysesBefore+tangentAnalyses;
        if(profiling){
            result["performance"]["linear_solver"]=lu.backend();
            result["performance"]["linear_analysis_reuse_enabled"]=reuseLinear;
            result["performance"]["linear_object_reused"]=linearObjectReused;
        }
        if(profiling && (densityIterations || localProjection || pseudoTransient)) {
            result["performance"]["density_coordinate_evaluation_seconds"]=densityEvaluationSeconds;
            result["performance"]["density_coordinate_evaluations"]=densityEvaluations;
        }
        if(profiling)result["performance"]["ialmob_screening_minimum_solves"]=physicsCallCounters.ialScreeningMinimumSolves-countersBefore.ialScreeningMinimumSolves;
        if(!predictorSelection.is_null())result["predictor_selection"]=predictorSelection;
        if(!tangentPreparation.is_null())result["tangent_preparation"]=tangentPreparation;
        if(residualRecovery)result["ngmres_recovery"]={{"updates",recoveryUpdates},{"reference_reassemblies",recoveryReassemblies},{"attempts",recoveryHistory}};
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
        linearGuard.accepted=stop=="diagnostic_scaled_residual";
        return result;
}
