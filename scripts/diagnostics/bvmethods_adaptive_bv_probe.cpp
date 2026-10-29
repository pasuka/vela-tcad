// Uninstalled controlled experiment. The physical residual is always the
// production residual; only the trial linearization and finite update vary.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/core/PerformanceProfiler.h"
#include "BvMethodsPopulationNewton.h"
#include "BvMethodsDerivativeSwitch.h"
#include "BvMethodsAdaptiveStep.h"
#include "BvMethodsExperimentJournal.h"

namespace {
using namespace vela;
using Json = nlohmann::json;
namespace fs = std::filesystem;
void matrixFile(const fs::path& path, const SparseMatrixd& a) {
    std::ofstream f(path); f << std::setprecision(17) << "row,col,value\n";
    for (int k=0; k<a.outerSize(); ++k)
        for (SparseMatrixd::InnerIterator v(a,k); v; ++v)
            f << v.row() << ',' << v.col() << ',' << v.value() << '\n';
}
}

int main(int argc, char** argv) { try {
    if (argc != 3) throw std::invalid_argument("manifest.json output-directory");
    Json m, cfg;
    { std::ifstream f(argv[1]); f >> m; }
    const std::string config=m.at("config");
    { std::ifstream f(config); f >> cfg; }
    if(m.contains("state_file")||m.contains("fields"))throw std::invalid_argument("autonomous BV forbids checkpoint input");
    const auto innerControl=diagnostics::resolveInnerStepLimit(cfg,m);
    const auto outerControl=diagnostics::resolveOuterStepControl(cfg,m);
    const std::string predictorMode=m.value("inner_predictor_mode",std::string{"enforce"});
    if(predictorMode!="enforce"&&predictorMode!="report_only")
        throw std::invalid_argument("inner_predictor_mode must be enforce or report_only");
    const bool reportPredictorOnly=predictorMode=="report_only";
    const std::string restorationMode=m.value("derivative_restoration_mode",std::string{"legacy"});
    if(restorationMode!="legacy"&&restorationMode!="automatic")
        throw std::invalid_argument("derivative_restoration_mode must be legacy or automatic");
    const bool automatic=restorationMode=="automatic";
    const std::string protection=m.value("population_protection",std::string{"fraction_to_boundary"});
    if(protection!="fraction_to_boundary"&&protection!="local_nonpositive"&&protection!="adaptive")
        throw std::invalid_argument("population_protection must be fraction_to_boundary, local_nonpositive or adaptive");
    const bool localProtection=protection=="local_nonpositive";
    const bool adaptiveProtection=protection=="adaptive";
    const double retainedFraction=m.value("local_retained_fraction",.01);
    if(!(retainedFraction>0.&&retainedFraction<1.))
        throw std::invalid_argument("local_retained_fraction must be between zero and one");
    const std::string selectionMetric=m.value("population_selection_metric",std::string{"target_normalized"});
    if(selectionMetric!="target_normalized"&&selectionMetric!="fixed_rows_l2")
        throw std::invalid_argument("population_selection_metric must be target_normalized or fixed_rows_l2");
    if(adaptiveProtection&&!automatic)throw std::invalid_argument("adaptive population protection requires automatic derivative restoration");
    const bool protectionReplay=m.contains("protection_checkpoint");
    const bool diagnosticIterationStates=m.value("diagnostic_iteration_states",false);
    const auto diagnosticNodes=m.value("diagnostic_nodes",std::vector<int>{});
    const Json optimization=m.value("evaluation_optimization",Json::object());
    const bool reuseCandidate=optimization.value("candidate_diagnostics",false);
    const bool reusePhysical=optimization.value("physical_reuse",false);
    const bool sharePorts=optimization.value("shared_ports",false);
    const bool boundaryParameter=optimization.value("boundary_parameter",false);
    const bool reuseDirection=optimization.value("direction_reuse",false);
    const bool reuseIteration=optimization.value("iteration_linearization_reuse",false);
    std::size_t directionBuilds=0,directionHits=0,qfDefectHits=0,identicalCandidates=0,linearizationHits=0,residualHits=0;
    const Json controls=m.value("derivative_switch",Json::object());
    diagnostics::DerivativeSwitchConfig sc;
    sc.minimumTrialSteps=controls.value("minimum_trial_steps",2);
    sc.maximumTrialSteps=controls.value("maximum_trial_steps",4);
    sc.failureConfirmations=controls.value("failure_confirmations",2);
    sc.bootstrapProbeInterval=controls.value("bootstrap_probe_interval",2);
    sc.cooldownSteps=controls.value("cooldown_steps",1);
    sc.minimumModelRatio=controls.value("minimum_model_ratio",.1);
    if(automatic){const diagnostics::DerivativeSwitch validated(sc);}
    cfg["sweep"].erase("initial_state_file");cfg["solver"]["warm_start"]=false;
    cfg=canonicalizeDeck(cfg);
    const fs::path out=argv[2]; fs::create_directories(out);
    const auto costStarted=std::chrono::steady_clock::now();
    PerformanceProfiler cost({m.value("cost_profiling",false),(out/"cost_profile.json").string()});
    ActivePerformanceProfilerScope costScope(&cost);
    auto timed=[&](std::string_view name,auto&& task){ScopedPerformanceTimer timer(name);return task();};
    const double inner=0.;double outer=0.;
    const int budget=m.value("max_iterations",40);
    if (budget<1 || budget>100) throw std::invalid_argument("bounded iteration budget required");
    for (auto& c:cfg["contacts"]) if(c["name"]=="drain") c["bias"]=inner;
    auto p=loadNewtonProblem(config,cfg);
    DDStateCliScope scope(p,cfg,configDirectory(config));
    if(p.newton.coupledOperator!="bvmethods_tri3_qf") throw std::invalid_argument("BVMethods profile required");
    p.newton.diagnosticBvMethodsEvaluationReuse=reusePhysical;
    p.newton.diagnosticBvMethodsSharedPorts=sharePorts;
    auto solver=makeNewtonSolver(p);
    auto sys=solver.makeArclengthSystem("drain",1e-4,boundaryParameter);
    auto unpack=solver.makeArclengthStateUnpacker();
    auto gates=solver.makeArclengthConvergenceEvaluator("drain");
    const int n=p.mesh.numNodes(),d=3*n;

    const auto& u=p.newton.inputScaling.unitSystem();
    auto cm=[&](double a){return u.internalConcentrationToM3(a)/1e6;};
    const double vt=constants::kb*p.newton.temperature_K/constants::q;
    const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*1e-6;
    auto port=solver.makeArclengthContactCurrentFunctional("drain",line);
    std::vector<ArclengthScalarFunctional> ports;
    for(const auto& c:p.mesh.contacts()) ports.push_back(solver.makeArclengthContactCurrentFunctional(c.name,line));
    const Json stateMetadata={{"mode","dd"},{"mesh_sha256",stateMeshIdentity(p.mesh,p.newton.inputScaling)},
        {"potential_origin_V",cfg.value("potential_origin_V",0.)},{"source_config_sha256",stateSha256(cfg.dump())},
        {"input_file_sha256",stateInputProvenance(cfg,configDirectory(config))},{"contact_biases_V",p.biases}};
    auto saveState=[&](const fs::path& path,const DDSolution& state,double bias){
        ScopedPerformanceTimer timer("driver.state_write");
        auto metadata=stateMetadata;metadata["contact_biases_V"]["drain"]=bias;
        DDStateArchiveScope dynamicScope(metadata);writeDDSolutionState(path,state,p.newton.inputScaling);
    };
    const auto poisson=timed("driver.cold_poisson",[&]{return solver.solvePoissonOnly(solver.buildPoissonBlockInitialization().coldInitial);});
    if(!poisson.converged)throw std::runtime_error("cold Poisson initialization failed");
    DDSolution seed=poisson.solution;
    saveState(out/"poisson_state.h5",seed,0.);
    VectorXd initial=solver.packArclengthState(seed); seed=unpack(initial);
    double replayBias=0.;
    if(protectionReplay) {
        const auto& checkpoint=m.at("protection_checkpoint");
        const auto current=readDDSolutionState(checkpoint.at("current_state_file").get<std::string>(),n,p.newton.inputScaling);
        const auto previous=readDDSolutionState(checkpoint.at("previous_state_file").get<std::string>(),n,p.newton.inputScaling);
        const double currentOuter=checkpoint.at("current_outer_V"),previousOuter=checkpoint.at("previous_outer_V");
        outer=checkpoint.at("target_outer_V");
        if(!std::isfinite(outer)||!(currentOuter>previousOuter)||!(outer>currentOuter))
            throw std::invalid_argument("ordered finite replay voltages required");
        const double ratio=checkpoint.value("predictor",true)?(outer-currentOuter)/(currentOuter-previousOuter):0.;
        const VectorXd currentPacked=solver.packArclengthState(current),previousPacked=solver.packArclengthState(previous);
        initial=currentPacked+ratio*(currentPacked-previousPacked);
        const double currentInner=checkpoint.at("current_inner_V"),previousInner=checkpoint.at("previous_inner_V");
        replayBias=currentInner+ratio*(currentInner-previousInner);seed=unpack(initial);
    }
    const double v0=solver.evaluateResidual(seed).potentialScale;
    double continuity=1.;
    if(p.newton.inputScaling.isUnitScaling()) {
        double epsr=1.,ni=1.;
        for(const auto& r:p.mesh.regions()) {
            const auto mat=p.matdb.getMaterial(r.material,p.newton.temperature_K);
            epsr=std::max(epsr,mat.eps_r); ni=std::max(ni,mat.ni);
        }
        const auto autoInputs=UnitScalingSystem::autoInputsFrom(p.mesh,p.doping,p.matdb,ni);
        const auto scaling=UnitScalingSystem::fromInputs(p.newton.temperature_K,constants::eps0*epsr,
            autoInputs,p.newton.unitScalingRefs,u);
        continuity=scaling.C0()*scaling.D0();
        if(std::abs(scaling.V0()/v0-1.)>1e-13) throw std::runtime_error("state scaling mismatch");
    }
    std::vector<double> nc(n),nv(n);
    std::vector<bool> active(n,false),efree(n,false),hfree(n,false);
    for(const auto& t:solver.evaluatePoissonTerms(seed).rows) {
        nc[t.nodeId]=t.electronDensityOfStates;nv[t.nodeId]=t.holeDensityOfStates;
        active[t.nodeId]=t.intrinsicDensity>0.;
    }
    for(const auto& t:solver.evaluateCarrierTermDiagnostics(seed,true).rows) {
        efree[t.nodeId]=t.electronContinuityActive; hfree[t.nodeId]=t.holeContinuityActive;
    }
    detail::BvMethodsTri3Context context{p.mesh,p.matdb,p.doping,
        {p.newton.mobility,p.newton.impactIonization,p.newton.inputScaling,CarrierStatisticsModel::FermiDiracJoyceDixon}};
    detail::BvMethodsQfState zero{VectorXd::Zero(n),VectorXd::Zero(n),VectorXd::Zero(n),seed.n,seed.p};
    detail::BvMethodsTri3Operator element(context,zero,nc,nv,v0,line*1e6*continuity,vt);
    VectorXd weights=VectorXd::Zero(d);
    for(const auto& c:p.mesh.contacts()) if(c.name=="drain") for(auto i:c.node_ids) {
        weights[n+i]-=constants::q;weights[2*n+i]+=constants::q;
    }
    const auto& circuit=cfg.at("sweep").at("external_circuit");
    const double resistance=circuit.at("resistance_ohm_um");
    const double eqtol=circuit.value("coupled_equation_tolerance",1e-8);
    const double loadtol=circuit.value("residual_tolerance_V",1e-6);
    DensityUpdateConvergenceConfig updateConfig;
    auto updateError=solver.makeArclengthUpdateErrorEvaluator("drain",updateConfig);
    auto coordinates=[&](const DDSolution& s,double bias) {
        VectorXd en(n),hp(n),es(n),hs(n);
        for(int i=0;i<n;++i){en[i]=cm(s.n[i]);hp[i]=cm(s.p[i]);es[i]=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.n[i],nc[i]))/vt:0.;hs[i]=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.p[i],nv[i]))/vt:0.;}
        return diagnostics::populationCoordinates(s.psi,en,hp,es,hs,active,v0,vt,bias);
    };
    auto border=[&](const SparseMatrixd& jac,const VectorXd& fv,const ArclengthScalarLinearization& t) {
        std::vector<Eigen::Triplet<double>> e;
        for(int k=0;k<jac.outerSize();++k)for(SparseMatrixd::InnerIterator v(jac,k);v;++v)e.emplace_back(v.row(),v.col(),v.value());
        for(int i=0;i<d;++i){if(fv[i]!=0.)e.emplace_back(i,d,fv[i]);if(t.stateDerivative[i]!=0.)e.emplace_back(d,i,resistance*t.stateDerivative[i]);}
        e.emplace_back(d,d,1.+resistance*t.parameterDerivative);
        SparseMatrixd a(d+1,d+1);a.setFromTriplets(e.begin(),e.end());return a;
    };
    auto project=[&](const SparseMatrixd& a) {
        std::vector<Eigen::Triplet<double>> e;
        for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator v(a,k);v;++v)
            if((v.row()>=n&&v.row()<2*n&&efree[v.row()-n])||(v.row()>=2*n&&hfree[v.row()-2*n]))e.emplace_back(v.row(),v.col(),v.value());
        SparseMatrixd b(d,d);b.setFromTriplets(e.begin(),e.end());return b;
    };
    VectorXd fixedRows=VectorXd::Ones(d+1);
    auto fullResidual=[&](const VectorXd& x,double bias) {ScopedPerformanceTimer timer("driver.full_residual");VectorXd f(d+1);f.head(d)=timed("driver.device_residual",[&]{return sys.residual(x,bias);});f[d]=bias+resistance*timed("driver.port_value",[&]{return port.value(x,bias);})-outer;return f;};
    auto inspect=[&](const VectorXd& x,double bias,const VectorXd* knownResidual=nullptr) {
        ScopedPerformanceTimer timer("driver.inspect");
        const VectorXd f=knownResidual?*knownResidual:fullResidual(x,bias);if(!f.allFinite())throw std::runtime_error("nonfinite physical RHS");
        const auto g=timed("driver.inspect.gates",[&]{return gates(x,bias);});const auto s=unpack(x);double kcl=0.;for(const auto& t:ports)kcl+=timed("driver.port_value",[&]{return t.value(x,bias);});
        const double id=timed("driver.port_value",[&]{return port.value(x,bias);}),factor=constants::q*line*continuity;
        auto carrier=[&](const NewtonGlobalContinuityCarrierClosure& c){return Json{{"ratio",c.ratio},{"qualified",c.qualified},
            {"contact_flux_A_per_um",c.contactFlux*factor},{"source_A_per_um",c.integratedSource*factor},
            {"signed_mismatch_A_per_um",c.mismatch*factor},{"net_denominator_A_per_um",c.netScale*factor},
            {"gross_contact_scale_A_per_um",c.contactNodeFluxAbsSum*factor},{"gross_source_scale_A_per_um",c.integratedSourceAbsSum*factor}};};
        Json result={{"equation_inf",f.head(d).lpNorm<Eigen::Infinity>()},{"load_V",f[d]},
            {"inner_V",bias},{"current_A_per_um",id},{"KCL_A_per_um",kcl},
            {"fixed_scaled_residual_inf",fixedRows.cwiseProduct(f).lpNorm<Eigen::Infinity>()},
            {"electron",carrier(g.globalContinuity.electron)},{"hole",carrier(g.globalContinuity.hole)},
            {"row_max_ratio",g.carrierRows.maxRatio},{"row_worst_node",g.carrierRows.maxRatioNode},{"row_worst_carrier",g.carrierRows.maxRatioCarrier},
            {"rows_satisfied",g.carrierRows.satisfied},{"global_satisfied",g.globalContinuity.satisfied},
            {"rows_enforced",g.carrierRows.enforced},{"global_enforced",g.globalContinuity.enforced}};
        result["original_gates_pass"]=f.head(d).lpNorm<Eigen::Infinity>()<=eqtol&&std::abs(f[d])<=loadtol&&
            std::abs(kcl)<=std::max(1e-12,1e-5*std::abs(id))&&(!g.carrierRows.enforced||g.carrierRows.satisfied)&&(!g.globalContinuity.enforced||g.globalContinuity.satisfied);
        ScopedPerformanceTimer ledgerTimer("driver.inspect.mp_ledger");
        using MP=linear_refinement::MP;std::array<MP,2> cf{},src{},allFlux{},freeResidual{},gross{};
        for(const auto& t:solver.evaluateCarrierTermDiagnostics(s,true).rows)for(int b=0;b<2;++b) {
            const double flux=b?t.holeFlux:t.electronFlux,rec=b?t.holeRecombination:t.electronRecombination,av=b?t.holeImpact:t.electronImpact;
            allFlux[b]+=MP(flux);
            if(active[t.nodeId]&&!(b?hfree[t.nodeId]:efree[t.nodeId])){cf[b]+=MP(flux);gross[b]+=abs(MP(flux));}
            else {src[b]+=MP(rec)+MP(av);if(b?hfree[t.nodeId]:efree[t.nodeId])freeResidual[b]+=MP(b?t.holeResidual:t.electronResidual);}
        }
        for(int b=0;b<2;++b){const MP mismatch=cf[b]-src[b];const MP floor=p.newton.globalContinuityClosure.sourceFloor;
            MP den=abs(cf[b]);if(abs(src[b])>den)den=abs(src[b]);if(floor>den)den=floor;
            result[b?"hole_mp_ledger":"electron_mp_ledger"]={{"signed_mismatch_A_per_um",static_cast<double>(mismatch*factor)},
                {"ratio",static_cast<double>(abs(mismatch)/den)},{"all_node_transport_sum_A_per_um",static_cast<double>(allFlux[b]*factor)},
                {"sum_free_equations_A_per_um",static_cast<double>(freeResidual[b]*factor)},
                {"ledger_rounding_defect_A_per_um",static_cast<double>((mismatch+freeResidual[b]-allFlux[b])*factor)}};
        }
        return result;
    };
    const double tau=m.value("fraction_tau",.99);
    if(!(tau>0.&&tau<1.))throw std::invalid_argument("invalid fraction margin");
    const double voltageUpdateLimit=circuit.value("max_inner_voltage_step_V",.1);
    const bool dumpSystems=m.value("dump_linear_systems",false);
    const double stopVoltage=cfg.at("sweep").value("stop",100000.);
    const double stopCurrent=cfg.at("sweep").at("breakdown").value("current_abs_A_per_um",.0013);
    const int attemptBudget=m.value("max_outer_attempts",1024);
    const bool extrapolate=cfg.at("sweep").value("predictor",std::string{"linear"})=="linear";
    Json ph=Json::array();for(const auto& t:poisson.history)ph.push_back({{"iteration",t.iter},{"residual",t.residualNorm},{"alpha",t.dampingFactor}});
    Json report={{"schema","vela.adaptive_derivative_bv.v1"},{"from_cold_poisson",!protectionReplay},{"imported_checkpoint",protectionReplay},
        {"backend",LinearSolver::selectedBackend()},{"poisson",{{"converged",poisson.converged},{"iterations",poisson.iters},{"residual",poisson.finalResidualNorm},{"trace",ph}}},
        {"resistance_ohm_um",resistance},{"equation_tolerance",eqtol},{"load_tolerance_V",loadtol},{"max_newton_updates",budget},
        {"population_protection",protection},{"population_selection_metric",selectionMetric},{"local_retained_fraction",retainedFraction},
        {"fraction_tau",tau},{"inner_predictor_limit_V",innerControl.volts},{"inner_predictor_source",innerControl.source},{"inner_predictor_reason",innerControl.reason},
        {"corrector_voltage_update_limit_V",voltageUpdateLimit},{"extrapolate",extrapolate},
        {"outer_initial_step_V",outerControl.initial},{"outer_min_step_V",outerControl.minimum},{"outer_max_step_V",outerControl.maximum},
        {"outer_growth",outerControl.growth},{"outer_shrink",outerControl.shrink},{"stop_current_A_per_um",stopCurrent},{"stop_outer_V",stopVoltage},
        {"inner_predictor_mode",predictorMode},{"derivative_restoration_mode",restorationMode},
        {"derivative_switch",controls},{"predictor_steps",Json::array()},
        {"evaluation_optimization",{{"candidate_diagnostics",reuseCandidate},{"physical_reuse",reusePhysical},
            {"shared_ports",sharePorts},{"boundary_parameter",boundaryParameter},
            {"direction_reuse",reuseDirection},{"iteration_linearization_reuse",reuseIteration}}},
        {"attempts",Json::array()},{"points",Json::array()},{"complete",false}};
    auto save=[&](){
        report["candidate_work"]={{"direction_builds",directionBuilds},{"direction_reuse_hits",directionHits},
            {"qf_defect_reuse_hits",qfDefectHits},{"identical_candidates_avoided",identicalCandidates},
            {"initial_linearization_reuse_hits",linearizationHits},{"residual_reuse_hits",residualHits}};
        diagnostics::writeExperimentJson(out/"summary.json",report);
    };save();
    struct Candidate {double phiBefore=0.,phiAfter=0.,modelAfter=0.;bool valid=false;VectorXd z,dx,dn,dp,x,F;double bias=0.,theta=0.,alpha=1.,merit=0.,rms=0.;Json log,diagnostics;};
    bool carriedFull=protectionReplay&&m.at("protection_checkpoint").value("initial_full_mode",true);
    auto correct=[&](const ArclengthState& start,double target,const std::string& label) {
        outer=target;VectorXd x=start.x;double bias=start.lambda;
        const fs::path folder=out/label;fs::create_directories(folder);
        diagnostics::DerivativeSwitch controller(sc,automatic&&carriedFull);
        bool bootstrap=false;Json attempt={{"target_outer_V",target},{"label",label},{"iterations",Json::array()},{"converged",false}};
        const auto first=unpack(x);const auto co=coordinates(first,bias);
        const detail::BvMethodsQfState refState{first.psi,first.phin,first.phip,first.n,first.p};
        const auto lt=element.transport(refState,true,true,true),ft=element.transport(refState,true,true,true,true),av0=element.evaluate(refState,true);
        const auto initialTerminal=port.linearize(x,bias);
        auto pt=initialTerminal;pt.stateDerivative+=line*continuity*((ft.jacobian-lt.jacobian).transpose()*weights);
        const auto initialJac=sys.jacobian(x,bias);const auto initialParameter=sys.parameterDerivative(x,bias);
        const SparseMatrixd ref=border(initialJac+project(ft.jacobian-lt.jacobian-av0.jacobian),initialParameter,pt)*co.qfFromScaledPopulation;
        fixedRows=diagnostics::equilibratePopulationSystem(ref).rowWeights;
        VectorXd cachedResidual;
        if(reuseIteration){cachedResidual=fullResidual(x,bias);attempt["initial"]=inspect(x,bias,&cachedResidual);}
        else attempt["initial"]=inspect(x,bias);
        // Preserve the original residual-only zero-bias initialization rule.
        if(target==0.&&attempt["initial"]["original_gates_pass"].get<bool>()) {
            attempt["converged"]=true;attempt["reason"]="original residual-only zero anchor";
            attempt["final_derivative_mode"]="bootstrap";
            attempt["last"]=attempt["initial"];return std::make_tuple(ArclengthState{x,bias},attempt);
        }
        LinearSolver linear;
        for(int iteration=1;iteration<=budget;++iteration)try {
            const auto s=unpack(x);const auto c=coordinates(s,bias);
            if(diagnosticIterationStates)saveState(folder/("before_"+std::to_string(iteration)+".h5"),s,bias);
            const VectorXd F=reuseIteration?cachedResidual:fullResidual(x,bias);
            if(reuseIteration)++residualHits;
            const bool initialReuse=reuseIteration&&iteration==1;
            if(initialReuse)++linearizationHits;
            const detail::BvMethodsQfState state{s.psi,s.phin,s.phip,s.n,s.p};
            const auto live=timed("driver.live_transport",[&]{return initialReuse?lt:element.transport(state,true,true,true);}),frozen=timed("driver.frozen_transport",[&]{return initialReuse?ft:element.transport(state,true,true,true,true);});
            if(!(live.source.array()==frozen.source.array()).all())throw std::runtime_error("frozen tangent changed RHS");
            const SparseMatrixd mu=frozen.jacobian-live.jacobian,av=timed("driver.av_jacobian",[&]{return initialReuse?av0:element.evaluate(state,true);}).jacobian;
            const auto fullJac=timed("driver.device_jacobian",[&]{return initialReuse?initialJac:sys.jacobian(x,bias);});const auto fv=timed("driver.parameter_derivative",[&]{return initialReuse?initialParameter:sys.parameterDerivative(x,bias);});const auto terminal=timed("driver.port_linearize",[&]{return initialReuse?initialTerminal:port.linearize(x,bias);});
            const VectorXd portMu=line*continuity*(mu.transpose()*weights);
            const SparseMatrixd physicalJac=automatic?border(fullJac,fv,terminal):SparseMatrixd{};
            if(dumpSystems&&iteration==1) {
                matrixFile(folder/"physical_bordered.csv",border(fullJac,fv,terminal));
                std::ofstream audit(folder/"physical_rhs_parameter.csv");
                audit<<std::setprecision(17)<<"row,F,dF_dV\n";
                for(int k=0;k<=d;++k)audit<<k<<','<<F[k]<<','<<(k<d?fv[k]:1.+resistance*terminal.parameterDerivative)<<'\n';
            }
            struct DirectionEntry {Candidate candidate;SparseMatrixd jq;};
            std::map<std::pair<double,bool>,DirectionEntry> directionCache;
            std::map<std::pair<double,bool>,double> defectCache;
            auto calculatePolicy=[&](double theta,bool avDerivatives,bool projection) {
                const auto key=std::make_pair(theta,avDerivatives);
                Candidate q;q.theta=theta;q.log={{"theta",theta},{"valid",false}};
                try {
                    SparseMatrixd jq;
                    const auto found=directionCache.find(key);
                    if(reuseDirection&&found!=directionCache.end()){q=found->second.candidate;jq=found->second.jq;++directionHits;}
                    else {
                    ++directionBuilds;
                    SparseMatrixd correction=(1.-theta)*mu;if(!avDerivatives)correction-=av;
                    const SparseMatrixd jac=fullJac+project(correction);auto pt=terminal;pt.stateDerivative+=(1.-theta)*portMu;
                    const auto constraintRows=diagonalConstraintRows(jac);
                    jq=border(jac,fv,pt);const SparseMatrixd jp=jq*c.qfFromScaledPopulation;
                    const std::string tag="iter"+std::to_string(iteration)+"_theta"+std::to_string(theta)+"_av"+std::to_string(avDerivatives);
                    if(dumpSystems&&iteration==1){matrixFile(folder/(tag+"_qf.csv"),jq);matrixFile(folder/(tag+"_population.csv"),jp);matrixFile(folder/"coordinate_map.csv",c.qfFromScaledPopulation);}
                    const auto eq=timed("driver.equilibration",[&]{return diagnostics::equilibratePopulationSystem(jp);});
                    q.z=eq.columnWeights.cwiseProduct(timed("driver.linear_solve",[&]{return linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(F)).eval(),constraintRows);}));
                    auto defect=timed("driver.population_defect",[&]{return linear_refinement::highPrecisionDefect(jp,F,q.z);});
                    if(defect.backwardMax>1e-12){q.z+=eq.columnWeights.cwiseProduct(timed("driver.linear_solve",[&]{return linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(defect.values)).eval(),constraintRows);}));defect=timed("driver.population_defect",[&]{return linear_refinement::highPrecisionDefect(jp,F,q.z);});}
                    q.log["linear_reduced_dimension"]=linear.lastReducedDimension();
                    q.log["population_linear_backward_error"]=defect.backwardMax;
                    if(defect.backwardMax>1e-10){
                        matrixFile(folder/(tag+"_failed_population.csv"),jp);matrixFile(folder/(tag+"_failed_qf.csv"),jq);
                        std::ofstream rf(folder/(tag+"_failed_rhs.csv"));rf<<std::setprecision(17)<<"row,F\n";for(int k=0;k<=d;++k)rf<<k<<','<<F[k]<<'\n';
                        throw std::runtime_error("population linear solve defect too large");}
                    q.dx=c.qfFromScaledPopulation*q.z;
                    if(dumpSystems&&iteration==1){std::ofstream raw(folder/(tag+"_direction.csv"));raw<<std::setprecision(17)<<"row,scaled_population,qf\n";for(int i=0;i<=d;++i)raw<<i<<','<<q.z[i]<<','<<q.dx[i]<<'\n';}
                    q.dn=c.electronScale_cm3.cwiseProduct(q.z.segment(n,n));q.dp=c.holeScale_cm3.cwiseProduct(q.z.segment(2*n,n));
                    double tangentRoundtrip=0.;
                    for(int i=0;i<n;++i)if(active[i]) {
                        const double dn=cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.n[i],nc[i]))*(q.dx[i]-q.dx[n+i])*v0/vt;
                        const double dp=cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.p[i],nv[i]))*(q.dx[2*n+i]-q.dx[i])*v0/vt;
                        tangentRoundtrip=std::max({tangentRoundtrip,std::abs(dn-q.dn[i])/c.electronScale_cm3[i],std::abs(dp-q.dp[i])/c.holeScale_cm3[i]});
                    }
                    q.log["scaled_tangent_roundtrip_error"]=tangentRoundtrip;
                    if(reuseDirection)directionCache.emplace(key,DirectionEntry{q,jq});
                    }
                    VectorXd en(n),hp(n);for(int i=0;i<n;++i){en[i]=cm(s.n[i]);hp[i]=cm(s.p[i]);}
                    auto fraction=diagnostics::FractionToBoundary{};
                    if(!projection)fraction=diagnostics::populationFractionToBoundary(en,hp,q.dn,q.dp,active,tau);
                    q.alpha=std::min(fraction.alpha,std::abs(q.dx[d])>voltageUpdateLimit?voltageUpdateLimit/std::abs(q.dx[d]):1.);q.bias=bias+q.alpha*q.dx[d];
                    DDSolution next=s;next.packedState.resize(0);next.packedLow.resize(0);next.electronQfReference.resize(0);next.holeQfReference.resize(0);next.phinIncrement.resize(0);next.phipIncrement.resize(0);
                    int negativeN=0,negativeP=0,protectedN=0,protectedP=0;double relativeN=0.,relativeP=0.;
                    double minimumN=std::numeric_limits<double>::infinity(),minimumP=minimumN;Json protectedNodes=Json::array();
                    for(int i=0;i<n;++i){next.psi[i]+=q.alpha*q.dx[i]*v0;
                        if(!active[i]){next.phin[i]+=q.alpha*q.dx[n+i]*v0;next.phip[i]+=q.alpha*q.dx[2*n+i]*v0;continue;}
                        double nn=en[i]+q.alpha*q.dn[i],pp=hp[i]+q.alpha*q.dp[i];
                        negativeN+=en[i]+q.dn[i]<=0.;negativeP+=hp[i]+q.dp[i]<=0.;
                        if(projection){
                            const auto a=diagnostics::localNonpositivePopulation(en[i],q.dn[i],q.alpha,retainedFraction);
                            const auto b=diagnostics::localNonpositivePopulation(hp[i],q.dp[i],q.alpha,retainedFraction);
                            nn=a.value;pp=b.value;protectedN+=a.protectedNonpositive;protectedP+=b.protectedNonpositive;
                            if(a.protectedNonpositive||b.protectedNonpositive)protectedNodes.push_back({{"node",i},{"electron",a.protectedNonpositive},{"hole",b.protectedNonpositive}});
                        }
                        minimumN=std::min(minimumN,nn/en[i]);minimumP=std::min(minimumP,pp/hp[i]);
                        if(!(nn>0.)||!(pp>0.)||!std::isfinite(nn)||!std::isfinite(pp))throw std::runtime_error("invalid population candidate");
                        next.n[i]=u.m3ToInternalConcentration(nn*1e6);next.p[i]=u.m3ToInternalConcentration(pp*1e6);
                        relativeN=std::max(relativeN,std::abs(q.dn[i]/en[i]));relativeP=std::max(relativeP,std::abs(q.dp[i]/hp[i]));
                    }
                    q.log.update({{"protected_electrons",protectedN},{"protected_holes",protectedP},
                        {"protected_nodes",protectedNodes},{"minimum_electron_ratio",minimumN},{"minimum_hole_ratio",minimumP},
                        {"population_protection",projection?"local_nonpositive":"fraction_to_boundary"},{"positivity_alpha",fraction.alpha}});
                    for(const auto& t:solver.evaluatePoissonTerms(next).rows){next.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;next.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
                    q.x=solver.packArclengthState(next);const auto check=unpack(q.x);double roundtrip=0.;
                    for(int i=0;i<n;++i)if(active[i])roundtrip=std::max({roundtrip,std::abs(check.n[i]/next.n[i]-1.),std::abs(check.p[i]/next.p[i]-1.)});
                    if(roundtrip>1e-10)throw std::runtime_error("population inverse roundtrip failed");
                    // A positive finite candidate must also be usable at the
                    // next Newton entry; do not accept a poisoned restart.
                    (void)coordinates(check,q.bias);
                    q.log["next_coordinate_map_finite"]=true;
                    q.log["population_column_scaling"]="relative_population";
                    q.F=fullResidual(q.x,q.bias);if(!q.F.allFinite())throw std::runtime_error("nonfinite candidate physical RHS");
                    q.merit=fixedRows.cwiseProduct(q.F).lpNorm<Eigen::Infinity>();
                    if(automatic){
                        const VectorXd model=F+q.alpha*(physicalJac*q.dx);
                        q.phiBefore=.5*fixedRows.cwiseProduct(F).squaredNorm();
                        q.phiAfter=.5*fixedRows.cwiseProduct(q.F).squaredNorm();
                        q.modelAfter=.5*fixedRows.cwiseProduct(model).squaredNorm();
                        q.log.update({{"phi_before",q.phiBefore},{"phi_after",q.phiAfter},{"model_phi_after",q.modelAfter}});
                    }
                    const auto er=timed("driver.update_error",[&]{return updateError(x,q.dx.head(d),q.x,q.bias,q.dx[d]);});
                    if(!er.valid||!std::isfinite(er.rms))throw std::runtime_error("invalid raw update RMS");q.rms=er.rms;
                    if(diagnosticIterationStates){q.log["update_block_RMS"]=er.blockRms;q.log["update_counts"]=er.counts;}
                    double qfBackward;
                    const auto existing=defectCache.find(key);
                    if(reuseDirection&&existing!=defectCache.end()){qfBackward=existing->second;++qfDefectHits;}
                    else {qfBackward=timed("driver.qf_defect",[&]{return linear_refinement::highPrecisionDefect(jq,F,q.dx);}).backwardMax;
                        if(reuseDirection)defectCache.emplace(key,qfBackward);}

                    q.log.update({{"valid",true},{"alpha",q.alpha},{"limiting_node",fraction.limitingNode},{"limiting_carrier",fraction.limitingCarrier},
                        {"RMS",q.rms},{"fixed_scaled_residual_inf",q.merit},{"equation_inf",q.F.head(d).lpNorm<Eigen::Infinity>()},{"load_V",q.F[d]},
                        {"qf_mapped_linear_backward_error",qfBackward},{"population_roundtrip_relative",roundtrip},
                        {"raw_dpsi_inf_V",q.dx.head(n).lpNorm<Eigen::Infinity>()*v0},{"raw_dfn_inf_V",q.dx.segment(n,n).lpNorm<Eigen::Infinity>()*v0},
                        {"raw_dfp_inf_V",q.dx.segment(2*n,n).lpNorm<Eigen::Infinity>()*v0},{"raw_voltage_V",q.dx[d]},
                        {"negative_electron_raw_candidates",negativeN},{"negative_hole_raw_candidates",negativeP},
                        {"max_relative_electron_tangent",relativeN},{"max_relative_hole_tangent",relativeP}});
                    q.valid=true;
                    if(diagnosticIterationStates){
                        Json nodes=Json::array();
                        for(int i:diagnosticNodes){
                            if(i<0||i>=n||!active[i])throw std::invalid_argument("diagnostic node must be in semiconductor");
                            nodes.push_back({{"node",i},{"n_before_cm3",en[i]},{"n_after_cm3",cm(check.n[i])},
                                {"p_before_cm3",hp[i]},{"p_after_cm3",cm(check.p[i])},
                                {"raw_dn_cm3",q.dn[i]},{"raw_dp_cm3",q.dp[i]},
                                {"delta_psi_V",check.psi[i]-s.psi[i]},
                                {"delta_fn_V",check.phin[i]-s.phin[i]},{"delta_fp_V",check.phip[i]-s.phip[i]},
                                {"electron_F_before",F[n+i]},{"electron_F_after",q.F[n+i]},
                                {"hole_F_before",F[2*n+i]},{"hole_F_after",q.F[2*n+i]}});
                        }
                        if(!diagnosticNodes.empty())q.log["node_updates"]=nodes;
                        const auto tag="candidate_"+std::to_string(iteration)+"_"+std::to_string(theta)+"_"+
                            std::to_string(avDerivatives)+"_"+(projection?"local":"fraction");
                        saveState(folder/(tag+".h5"),check,q.bias);
                        q.log["candidate_inner_V"]=q.bias;
                        diagnostics::writeExperimentJson(folder/(tag+".json"),q.log);
                    }
                }catch(const std::exception& e){q.log["error"]=e.what();}
                return q;
            };
            auto calculate=[&](double theta,bool avDerivatives) {
                if(!adaptiveProtection)return calculatePolicy(theta,avDerivatives,localProtection);
                auto fraction=calculatePolicy(theta,avDerivatives,false);
                const double circuitAlpha=std::abs(fraction.dx.size()?fraction.dx[d]:0.)>voltageUpdateLimit?
                    voltageUpdateLimit/std::abs(fraction.dx[d]):1.;
                if(fraction.valid&&fraction.alpha==circuitAlpha){
                    ++identicalCandidates;
                    fraction.log["population_selection"]={{"selected","fraction_to_boundary"},{"reason","identical_positive_candidate"}};
                    return fraction;
                }
                auto local=calculatePolicy(theta,avDerivatives,true);
                auto quality=[&](Candidate& q){
                    bool joint=false;
                    if(q.valid&&q.rms<1.&&q.F.head(d).lpNorm<Eigen::Infinity>()<=eqtol&&std::abs(q.F[d])<=loadtol){
                        q.diagnostics=inspect(q.x,q.bias,reusePhysical?&q.F:nullptr);
                        joint=q.diagnostics["original_gates_pass"].get<bool>();
                    }
                    const double score=selectionMetric=="target_normalized"&&q.valid?
                        diagnostics::populationTargetMerit(q.F.head(d).lpNorm<Eigen::Infinity>(),q.F[d],eqtol,loadtol):q.phiAfter;
                    return diagnostics::PopulationCandidateQuality{q.valid,joint,score};
                };
                const auto a=quality(fraction),b=quality(local);
                const auto choice=diagnostics::selectPopulationCandidate(a,b);
                Json selection={{"fraction",{{"valid",a.valid},{"joint",a.joint},{"merit",std::isfinite(a.merit)?Json(a.merit):Json(nullptr)},
                    {"alpha",fraction.alpha}}},{"local",{{"valid",b.valid},{"joint",b.joint},{"merit",std::isfinite(b.merit)?Json(b.merit):Json(nullptr)},
                    {"alpha",local.alpha}}},{"metric",selectionMetric},{"reason",a.joint||b.joint?"original_joint_gates":"physical_scaled_residual"}};
                selection["selected"]=choice==diagnostics::PopulationChoice::Local?"local_nonpositive":choice==diagnostics::PopulationChoice::Fraction?"fraction_to_boundary":"none";
                auto selected=choice==diagnostics::PopulationChoice::Local?std::move(local):std::move(fraction);
                if(choice==diagnostics::PopulationChoice::None){selected.valid=false;selected.log["valid"]=false;selected.log["error"]="no usable population candidate metric";}
                selected.log["population_selection"]=selection;
                return selected;
            };
            Json trials=Json::array();Candidate chosen;Json switchLog;
            if(automatic) {
                const auto beforeMode=controller.mode();
                const Json& beforeDiag=attempt.contains("last")?attempt["last"]:attempt["initial"];
                const bool requested=controller.requestFull(beforeDiag["original_gates_pass"].get<bool>());
                const std::string entryReason=controller.reason;
                bool fallback=false;
                if(requested) {
                    auto candidate=calculate(1.,true);trials.push_back(candidate.log);
                    diagnostics::DerivativeObservation observation;
                    observation.valid=candidate.valid;
                    if(candidate.valid) {
                        const auto dg=reuseCandidate&&!candidate.diagnostics.is_null()?candidate.diagnostics:inspect(candidate.x,candidate.bias,reusePhysical?&candidate.F:nullptr);
                        if(reuseCandidate)candidate.diagnostics=dg;
                        observation.physicalGates=dg["original_gates_pass"].get<bool>();
                        observation.before=candidate.phiBefore;observation.after=candidate.phiAfter;observation.modelAfter=candidate.modelAfter;
                    }
                    const auto quality=diagnostics::derivativeQuality(observation,sc.minimumModelRatio);
                    const bool applied=controller.fullCandidate(observation);
                    switchLog={{"before",diagnostics::derivativeModeName(beforeMode)},{"requested_full",true},
                        {"entry_reason",entryReason},{"full_candidate_applied",applied},
                        {"predicted_decrease",quality.predicted},{"actual_decrease",quality.actual},
                        {"ratio_reliable",quality.reliable},{"rho",quality.reliable?Json(quality.ratio):Json(nullptr)},
                        {"candidate_physical_gates",observation.physicalGates},{"reason",controller.reason}};
                    if(applied)chosen=std::move(candidate);
                    else {fallback=true;switchLog["rejected_candidate_written_to_state"]=false;}
                } else switchLog={{"before",diagnostics::derivativeModeName(beforeMode)},{"requested_full",false},{"reason",entryReason}};
                if(!requested||fallback){chosen=calculate(0.,false);trials.push_back(chosen.log);if(chosen.valid)controller.bootstrapApplied();}
                bootstrap=controller.mode()==diagnostics::DerivativeMode::Bootstrap;
                switchLog["after"]=diagnostics::derivativeModeName(controller.mode());
                switchLog["trial_steps"]=controller.trialSteps();switchLog["bad_steps"]=controller.badSteps();
            }
            if(!automatic&&!bootstrap) {
                chosen=calculate(0.,true);trials.push_back(chosen.log);
                for(double theta:{1.,.5,.25}) {
                    auto q=calculate(theta,true);const bool within=q.valid&&diagnostics::restorationWithinTargets(q.F.head(d).lpNorm<Eigen::Infinity>(),q.F[d],q.rms,eqtol,loadtol);
                    q.log["within_existing_targets"]=within;trials.push_back(q.log);
                    if(q.valid&&(!chosen.valid||within||diagnostics::restorationNoWorse(q.merit,q.rms,chosen.merit,chosen.rms))){chosen=std::move(q);break;}
                }
                const bool within=chosen.valid&&diagnostics::restorationWithinTargets(chosen.F.head(d).lpNorm<Eigen::Infinity>(),chosen.F[d],chosen.rms,eqtol,loadtol);
                if(diagnostics::derivativeBootstrapRequired(chosen.valid,chosen.alpha,chosen.merit,fixedRows.cwiseProduct(F).lpNorm<Eigen::Infinity>(),within)) {
                    bootstrap=true;attempt["bootstrap_used"]=true;
                    trials.push_back({{"phase_change","combination bootstrap from current own state"}});
                }
            }
            if(!automatic&&bootstrap){chosen=calculate(0.,false);trials.push_back(chosen.log);}
            if(!chosen.valid){attempt["failure"]="no valid population candidate";attempt["failed_trials"]=trials;break;}
            auto diag=reuseCandidate&&!chosen.diagnostics.is_null()?chosen.diagnostics:
                inspect(chosen.x,chosen.bias,reusePhysical?&chosen.F:nullptr);
            const bool joint=diag["original_gates_pass"].get<bool>()&&chosen.rms<1.;
            Json row={{"iteration",iteration},{"phase",bootstrap?"combination_bootstrap":"full_AV_adaptive_mobility"},
                {"theta",chosen.theta},{"alpha",chosen.alpha},{"raw_RMS",chosen.rms},{"diagnostics",diag},{"trials",trials},{"switch",switchLog},
                {"population_update",{{"selection",chosen.log.value("population_selection",Json::object())},{"policy",chosen.log["population_protection"]},{"protected_electrons",chosen.log["protected_electrons"]},{"protected_holes",chosen.log["protected_holes"]},
                    {"protected_nodes",chosen.log["protected_nodes"]},{"minimum_electron_ratio",chosen.log["minimum_electron_ratio"]},
                    {"minimum_hole_ratio",chosen.log["minimum_hole_ratio"]},{"positivity_alpha",chosen.log["positivity_alpha"]}}}};
            if(reuseIteration)cachedResidual=chosen.F;
            x=chosen.x;bias=chosen.bias;attempt["iterations"].push_back(row);attempt["last"]=diag;
            diagnostics::writeExperimentJson(folder/"attempt.json",attempt);
            diagnostics::writeExperimentJson(out/"progress.json",{{"attempt",label},{"target_outer_V",target},
                {"iteration",iteration},{"raw_RMS",chosen.rms},{"theta",chosen.theta},{"alpha",chosen.alpha},
                {"diagnostics",diag},{"switch",switchLog}});
            const bool restored=automatic?controller.mode()==diagnostics::DerivativeMode::Full:!bootstrap;
            if(joint&&restored){attempt["converged"]=true;attempt["reason"]="original gates and reported population RMS after AV restoration";break;}
            if(joint&&bootstrap&&!automatic){bootstrap=false;attempt["AV_restoration_requested"]=true;}
        }catch(const std::exception& e){attempt["failure"]=e.what();break;}
        if(!attempt["converged"].get<bool>()&&!attempt.contains("failure"))attempt["failure"]="Newton update budget exhausted";
        if(automatic)attempt["final_derivative_mode"]=diagnostics::derivativeModeName(controller.mode());
        saveState(folder/"final_state.h5",unpack(x),bias);
        return std::make_tuple(ArclengthState{x,bias},attempt);
    };
    auto attempt=[&](const ArclengthState& s,double target,const std::string& kind) {
        const auto label="attempt_"+std::to_string(report["attempts"].size())+"_"+kind;
        ArclengthState result=s;Json rec;
        try{std::tie(result,rec)=correct(s,target,label);}catch(const std::exception& e){rec={{"converged",false},{"failure",e.what()},{"target_outer_V",target},{"label",label}};}
        diagnostics::writeExperimentJson(out/label/"attempt.json",rec);
        report["attempts"].push_back(rec);save();return std::make_tuple(result,rec);
    };
    if(protectionReplay) {
        saveState(out/"shared_initial.h5",unpack(initial),replayBias);
        auto [last,rec]=attempt({initial,replayBias},outer,"protection_checkpoint");
        report["checkpoint_converged"]=rec["converged"];
        double acceptedReplayOuter=outer;
        report["complete"]=false;report["stop_reason"]="bounded population-protection checkpoint comparison";
        bool replayPassed=rec["converged"].get<bool>();
        if(replayPassed){
            report["points"].push_back({{"outer_V",outer},{"inner_V",last.lambda},{"updates",rec["iterations"].size()},{"diagnostics",rec["last"]}});
            saveState(out/"checkpoint_state_0.h5",unpack(last.x),last.lambda);
        }
        const auto& cp=m.at("protection_checkpoint");
        auto previousSeed=readDDSolutionState(cp.at("current_state_file").get<std::string>(),n,p.newton.inputScaling);
        ArclengthState previousReplay{solver.packArclengthState(previousSeed),cp.at("current_inner_V").get<double>()};
        double previousReplayOuter=cp.at("current_outer_V");
        const auto followup=m.value("checkpoint_followup_targets_V",std::vector<double>{});
        for(double target:followup) {
            if(!replayPassed)break;
            if(!std::isfinite(target)||!(target>acceptedReplayOuter))throw std::invalid_argument("increasing finite followup targets required");
            carriedFull=rec.value("final_derivative_mode",std::string{})=="full";
            const double ratio=(target-acceptedReplayOuter)/(acceptedReplayOuter-previousReplayOuter);
            ArclengthState prediction{last.x+ratio*(last.x-previousReplay.x),last.lambda+ratio*(last.lambda-previousReplay.lambda)};
            auto [next,nextRec]=attempt(prediction,target,"followup_predictor");
            if(!nextRec["converged"].get<bool>())std::tie(next,nextRec)=attempt(last,target,"followup_accepted_retry");
            replayPassed=nextRec["converged"].get<bool>();
            if(!replayPassed){report["followup_failure"]=nextRec.value("failure",std::string{"unknown"});break;}
            previousReplay=last;previousReplayOuter=acceptedReplayOuter;last=next;rec=nextRec;acceptedReplayOuter=target;
            report["points"].push_back({{"outer_V",target},{"inner_V",last.lambda},{"updates",rec["iterations"].size()},{"diagnostics",rec["last"]}});
            saveState(out/("checkpoint_state_"+std::to_string(report["points"].size()-1)+".h5"),unpack(last.x),last.lambda);save();
        }
        report["checkpoint_followup_complete"]=replayPassed;
        report["last_outer_V"]=acceptedReplayOuter;report["last_inner_V"]=last.lambda;
        save();cost.recordStage("driver.total",std::chrono::steady_clock::now()-costStarted);cost.writeJson();
        return replayPassed?0:3;
    }
    auto [current,zeroRec]=attempt({initial,0.},0.,"zero");
    if(!zeroRec["converged"].get<bool>()){report["failure"]="zero-bias coupled initialization failed";save();return 2;}
    double acceptedOuter=0.,previousOuter=0.;ArclengthState previous=current;bool hasPrevious=false;
    report["points"].push_back({{"outer_V",0.},{"inner_V",current.lambda},{"diagnostics",zeroRec["last"]},{"updates",zeroRec["iterations"].size()}});
    saveState(out/"state_0.h5",unpack(current.x),current.lambda);
    std::ofstream iv(out/"iv.csv");iv<<std::setprecision(17)<<"outer_V,inner_V,current_A_per_um,updates,electron_closure,hole_closure\n";
    auto ivrow=[&](double v,const ArclengthState& s,const Json& rec){iv<<v<<','<<s.lambda<<','<<rec["last"]["current_A_per_um"]<<','<<rec["iterations"].size()<<','<<rec["last"]["electron"]["ratio"]<<','<<rec["last"]["hole"]["ratio"]<<'\n';iv.flush();};
    ivrow(0.,current,zeroRec);double step=outerControl.initial;int retries=0;
    while(report["attempts"].size()<static_cast<std::size_t>(attemptBudget)) {
        const auto decision=hasPrevious?diagnostics::predictorStepDecision(step,outerControl,innerControl.volts,
            current.lambda-previous.lambda,acceptedOuter-previousOuter,reportPredictorOnly):
            diagnostics::PredictorStepDecision{step,step};
        const auto allowed=decision.applied;
        Json stepRecord={{"accepted_outer_V",acceptedOuter},{"proposed_outer_step_V",step},{"mode",predictorMode},
            {"guard_below_minimum",!decision.guarded.has_value()}};
        stepRecord["guarded_outer_step_V"]=decision.guarded?Json(*decision.guarded):Json(nullptr);
        stepRecord["applied_outer_step_V"]=allowed?Json(*allowed):Json(nullptr);
        stepRecord["predicted_abs_inner_step_V"]=hasPrevious&&allowed?
            Json(std::abs(*allowed*(current.lambda-previous.lambda)/(acceptedOuter-previousOuter))):Json(nullptr);
        report["predictor_steps"].push_back(stepRecord);
        if(!allowed){report["failure"]="inner predictor cap requires sub-minimum outer step";break;}
        const double target=std::min(stopVoltage,acceptedOuter+*allowed);ArclengthState prediction=current;
        if(hasPrevious&&extrapolate) {
            const double ratio=(target-acceptedOuter)/(acceptedOuter-previousOuter);
            prediction.x=current.x+ratio*(current.x-previous.x);prediction.lambda=current.lambda+ratio*(current.lambda-previous.lambda);
        }
        auto [candidate,rec]=attempt(prediction,target,"predictor");bool fallback=false;
        if(!rec["converged"].get<bool>()&&hasPrevious&&extrapolate){std::tie(candidate,rec)=attempt(current,target,"accepted_state_retry");fallback=true;}
        if(rec["converged"].get<bool>()) {
            if(automatic)carriedFull=rec.value("final_derivative_mode",std::string{})=="full";
            const double used=target-acceptedOuter;previous=current;previousOuter=acceptedOuter;current=candidate;acceptedOuter=target;hasPrevious=true;retries=0;
            const auto& dg=rec["last"];const int updates=rec["iterations"].size();
            report["points"].push_back({{"outer_V",target},{"inner_V",current.lambda},{"updates",updates},{"diagnostics",dg},{"fallback_used",fallback}});
            saveState(out/("state_"+std::to_string(report["points"].size()-1)+".h5"),unpack(current.x),current.lambda);ivrow(target,current,rec);
            std::cout<<"accepted outer="<<std::setprecision(17)<<target<<" inner="<<current.lambda<<" I="<<dg["current_A_per_um"]<<" updates="<<updates<<std::endl;
            if(std::abs(dg["current_A_per_um"].get<double>())>=stopCurrent||target>=stopVoltage){report["complete"]=true;report["stop_reason"]=target>=stopVoltage?"outer limit":"current threshold";break;}
            const double growth=fallback?1.:(outerControl.growthMode=="fixed"?outerControl.growth:1.+(outerControl.growth-1.)*std::max(0.,1.-std::max(updates-1,0)/(.75*budget)));
            step=std::clamp(used*growth,outerControl.minimum,outerControl.maximum);save();
        } else {
            std::cout<<"rejected outer="<<std::setprecision(17)<<target<<" reason="<<rec.value("failure",std::string{"unknown"})<<std::endl;
            if(*allowed<=outerControl.minimum*(1+1e-12)||++retries>outerControl.maxRetries){report["failure"]="minimum outer step or retries exhausted";break;}
            step=std::max(outerControl.minimum,*allowed*outerControl.shrink);
        }
    }
    report["last_outer_V"]=acceptedOuter;report["last_inner_V"]=current.lambda;
    if(!report["complete"].get<bool>()&&!report.contains("failure"))report["failure"]="outer attempt budget exhausted";
    save();cost.recordStage("driver.total",std::chrono::steady_clock::now()-costStarted);cost.writeJson();return report["complete"].get<bool>()?0:3;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
