// Uninstalled controlled experiment. The physical residual is always the
// production residual; only the trial linearization and finite update vary.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "BvMethodsPopulationNewton.h"
#include "BvMethodsAdaptiveStep.h"

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
    cfg["sweep"].erase("initial_state_file");cfg["solver"]["warm_start"]=false;
    cfg=canonicalizeDeck(cfg);
    const fs::path out=argv[2]; fs::create_directories(out);
    const double inner=0.;double outer=0.;
    const int budget=m.value("max_iterations",40);
    if (budget<1 || budget>100) throw std::invalid_argument("bounded iteration budget required");
    for (auto& c:cfg["contacts"]) if(c["name"]=="drain") c["bias"]=inner;
    auto p=loadNewtonProblem(config,cfg);
    DDStateCliScope scope(p,cfg,configDirectory(config));
    if(p.newton.coupledOperator!="bvmethods_tri3_qf") throw std::invalid_argument("BVMethods profile required");
    auto solver=makeNewtonSolver(p);
    auto sys=solver.makeArclengthSystem("drain");
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
        auto metadata=stateMetadata;metadata["contact_biases_V"]["drain"]=bias;
        DDStateArchiveScope dynamicScope(metadata);writeDDSolutionState(path,state,p.newton.inputScaling);
    };
    const auto poisson=solver.solvePoissonOnly(solver.buildPoissonBlockInitialization().coldInitial);
    if(!poisson.converged)throw std::runtime_error("cold Poisson initialization failed");
    DDSolution seed=poisson.solution;
    saveState(out/"poisson_state.h5",seed,0.);
    const VectorXd initial=solver.packArclengthState(seed); seed=unpack(initial);
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
    auto fullResidual=[&](const VectorXd& x,double bias) {VectorXd f(d+1);f.head(d)=sys.residual(x,bias);f[d]=bias+resistance*port.value(x,bias)-outer;return f;};
    auto inspect=[&](const VectorXd& x,double bias) {
        const VectorXd f=fullResidual(x,bias);if(!f.allFinite())throw std::runtime_error("nonfinite physical RHS");
        const auto g=gates(x,bias);const auto s=unpack(x);double kcl=0.;for(const auto& t:ports)kcl+=t.value(x,bias);
        const double id=port.value(x,bias),factor=constants::q*line*continuity;
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
    Json report={{"schema","vela.adaptive_derivative_bv.v1"},{"from_cold_poisson",true},{"imported_checkpoint",false},
        {"backend",LinearSolver::selectedBackend()},{"poisson",{{"converged",poisson.converged},{"iterations",poisson.iters},{"residual",poisson.finalResidualNorm},{"trace",ph}}},
        {"resistance_ohm_um",resistance},{"equation_tolerance",eqtol},{"load_tolerance_V",loadtol},{"max_newton_updates",budget},
        {"fraction_tau",tau},{"inner_predictor_limit_V",innerControl.volts},{"inner_predictor_source",innerControl.source},{"inner_predictor_reason",innerControl.reason},
        {"corrector_voltage_update_limit_V",voltageUpdateLimit},{"extrapolate",extrapolate},
        {"outer_initial_step_V",outerControl.initial},{"outer_min_step_V",outerControl.minimum},{"outer_max_step_V",outerControl.maximum},
        {"outer_growth",outerControl.growth},{"outer_shrink",outerControl.shrink},{"stop_current_A_per_um",stopCurrent},{"stop_outer_V",stopVoltage},
        {"attempts",Json::array()},{"points",Json::array()},{"complete",false}};
    auto save=[&](){std::ofstream(out/"summary.json")<<report.dump(2)<<'\n';};save();
    struct Candidate {bool valid=false;VectorXd z,dx,dn,dp,x,F;double bias=0.,theta=0.,alpha=1.,merit=0.,rms=0.;Json log;};
    auto correct=[&](const ArclengthState& start,double target,const std::string& label) {
        outer=target;VectorXd x=start.x;double bias=start.lambda;
        const fs::path folder=out/label;fs::create_directories(folder);
        bool bootstrap=false;Json attempt={{"target_outer_V",target},{"label",label},{"iterations",Json::array()},{"converged",false}};
        const auto first=unpack(x);const auto co=coordinates(first,bias);
        const detail::BvMethodsQfState refState{first.psi,first.phin,first.phip,first.n,first.p};
        const auto lt=element.transport(refState,true,true,true),ft=element.transport(refState,true,true,true,true),av0=element.evaluate(refState,true);
        auto pt=port.linearize(x,bias);pt.stateDerivative+=line*continuity*((ft.jacobian-lt.jacobian).transpose()*weights);
        const SparseMatrixd ref=border(sys.jacobian(x,bias)+project(ft.jacobian-lt.jacobian-av0.jacobian),sys.parameterDerivative(x,bias),pt)*co.qfFromScaledPopulation;
        fixedRows=diagnostics::equilibratePopulationSystem(ref).rowWeights;
        attempt["initial"]=inspect(x,bias);
        // Preserve the original residual-only zero-bias initialization rule.
        if(target==0.&&attempt["initial"]["original_gates_pass"].get<bool>()) {
            attempt["converged"]=true;attempt["reason"]="original residual-only zero anchor";
            attempt["last"]=attempt["initial"];return std::make_tuple(ArclengthState{x,bias},attempt);
        }
        LinearSolver linear;
        for(int iteration=1;iteration<=budget;++iteration)try {
            const auto s=unpack(x);const auto c=coordinates(s,bias);const VectorXd F=fullResidual(x,bias);
            const detail::BvMethodsQfState state{s.psi,s.phin,s.phip,s.n,s.p};
            const auto live=element.transport(state,true,true,true),frozen=element.transport(state,true,true,true,true);
            if(!(live.source.array()==frozen.source.array()).all())throw std::runtime_error("frozen tangent changed RHS");
            const SparseMatrixd mu=frozen.jacobian-live.jacobian,av=element.evaluate(state,true).jacobian;
            const auto fullJac=sys.jacobian(x,bias);const auto fv=sys.parameterDerivative(x,bias);const auto terminal=port.linearize(x,bias);
            const VectorXd portMu=line*continuity*(mu.transpose()*weights);
            auto calculate=[&](double theta,bool avDerivatives) {
                const bool projection=false;Candidate q;q.theta=theta;q.log={{"theta",theta},{"valid",false}};
                try {
                    SparseMatrixd correction=(1.-theta)*mu;if(!avDerivatives)correction-=av;
                    const SparseMatrixd jac=fullJac+project(correction);auto pt=terminal;pt.stateDerivative+=(1.-theta)*portMu;
                    const SparseMatrixd jq=border(jac,fv,pt),jp=jq*c.qfFromScaledPopulation;
                    const std::string tag="iter"+std::to_string(iteration)+"_theta"+std::to_string(theta)+"_av"+std::to_string(avDerivatives);
                    if(dumpSystems&&iteration==1){matrixFile(folder/(tag+"_qf.csv"),jq);matrixFile(folder/(tag+"_population.csv"),jp);matrixFile(folder/"coordinate_map.csv",c.qfFromScaledPopulation);}
                    const auto eq=diagnostics::equilibratePopulationSystem(jp);
                    q.z=eq.columnWeights.cwiseProduct(linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(F)).eval()));
                    auto defect=linear_refinement::highPrecisionDefect(jp,F,q.z);
                    if(defect.backwardMax>1e-12){q.z+=eq.columnWeights.cwiseProduct(linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(defect.values)).eval()));defect=linear_refinement::highPrecisionDefect(jp,F,q.z);}
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
                    VectorXd en(n),hp(n);for(int i=0;i<n;++i){en[i]=cm(s.n[i]);hp[i]=cm(s.p[i]);}
                    auto fraction=diagnostics::FractionToBoundary{};
                    if(!projection)fraction=diagnostics::populationFractionToBoundary(en,hp,q.dn,q.dp,active,tau);
                    q.alpha=std::min(fraction.alpha,std::abs(q.dx[d])>voltageUpdateLimit?voltageUpdateLimit/std::abs(q.dx[d]):1.);q.bias=bias+q.alpha*q.dx[d];
                    DDSolution next=s;next.packedState.resize(0);next.packedLow.resize(0);next.electronQfReference.resize(0);next.holeQfReference.resize(0);next.phinIncrement.resize(0);next.phipIncrement.resize(0);
                    int negativeN=0,negativeP=0;double relativeN=0.,relativeP=0.;
                    for(int i=0;i<n;++i){next.psi[i]+=q.alpha*q.dx[i]*v0;
                        if(!active[i]){next.phin[i]+=q.alpha*q.dx[n+i]*v0;next.phip[i]+=q.alpha*q.dx[2*n+i]*v0;continue;}
                        double nn=en[i]+q.alpha*q.dn[i],pp=hp[i]+q.alpha*q.dp[i];
                        negativeN+=en[i]+q.dn[i]<=0.;negativeP+=hp[i]+q.dp[i]<=0.;
                        if(projection){if(nn<=0.)nn=.01*en[i];if(pp<=0.)pp=.01*hp[i];}
                        if(!(nn>0.)||!(pp>0.)||!std::isfinite(nn)||!std::isfinite(pp))throw std::runtime_error("invalid population candidate");
                        next.n[i]=u.m3ToInternalConcentration(nn*1e6);next.p[i]=u.m3ToInternalConcentration(pp*1e6);
                        relativeN=std::max(relativeN,std::abs(q.dn[i]/en[i]));relativeP=std::max(relativeP,std::abs(q.dp[i]/hp[i]));
                    }
                    for(const auto& t:solver.evaluatePoissonTerms(next).rows){next.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;next.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
                    q.x=solver.packArclengthState(next);const auto check=unpack(q.x);double roundtrip=0.;
                    for(int i=0;i<n;++i)if(active[i])roundtrip=std::max({roundtrip,std::abs(check.n[i]/next.n[i]-1.),std::abs(check.p[i]/next.p[i]-1.)});
                    if(roundtrip>1e-10)throw std::runtime_error("population inverse roundtrip failed");
                    q.F=fullResidual(q.x,q.bias);if(!q.F.allFinite())throw std::runtime_error("nonfinite candidate physical RHS");
                    q.merit=fixedRows.cwiseProduct(q.F).lpNorm<Eigen::Infinity>();
                    const auto er=updateError(x,q.dx.head(d),q.x,q.bias,q.dx[d]);
                    if(!er.valid||!std::isfinite(er.rms))throw std::runtime_error("invalid raw update RMS");q.rms=er.rms;
                    const auto qfDefect=linear_refinement::highPrecisionDefect(jq,F,q.dx);
                    q.log.update({{"valid",true},{"alpha",q.alpha},{"limiting_node",fraction.limitingNode},{"limiting_carrier",fraction.limitingCarrier},
                        {"RMS",q.rms},{"fixed_scaled_residual_inf",q.merit},{"equation_inf",q.F.head(d).lpNorm<Eigen::Infinity>()},{"load_V",q.F[d]},
                        {"qf_mapped_linear_backward_error",qfDefect.backwardMax},{"population_roundtrip_relative",roundtrip},
                        {"raw_dpsi_inf_V",q.dx.head(n).lpNorm<Eigen::Infinity>()*v0},{"raw_dfn_inf_V",q.dx.segment(n,n).lpNorm<Eigen::Infinity>()*v0},
                        {"raw_dfp_inf_V",q.dx.segment(2*n,n).lpNorm<Eigen::Infinity>()*v0},{"raw_voltage_V",q.dx[d]},
                        {"negative_electron_raw_candidates",negativeN},{"negative_hole_raw_candidates",negativeP},
                        {"max_relative_electron_tangent",relativeN},{"max_relative_hole_tangent",relativeP}});
                    q.valid=true;
                }catch(const std::exception& e){q.log["error"]=e.what();}
                return q;
            };
            Json trials=Json::array();Candidate chosen;
            if(!bootstrap) {
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
            if(bootstrap){chosen=calculate(0.,false);trials.push_back(chosen.log);}
            if(!chosen.valid){attempt["failure"]="no valid population candidate";attempt["failed_trials"]=trials;break;}
            auto diag=inspect(chosen.x,chosen.bias);
            const bool joint=diag["original_gates_pass"].get<bool>()&&chosen.rms<1.;
            Json row={{"iteration",iteration},{"phase",bootstrap?"combination_bootstrap":"full_AV_adaptive_mobility"},
                {"theta",chosen.theta},{"alpha",chosen.alpha},{"raw_RMS",chosen.rms},{"diagnostics",diag},{"trials",trials}};
            x=chosen.x;bias=chosen.bias;attempt["iterations"].push_back(row);attempt["last"]=diag;
            if(joint&&!bootstrap){attempt["converged"]=true;attempt["reason"]="original gates and reported population RMS after AV restoration";break;}
            if(joint&&bootstrap){bootstrap=false;attempt["AV_restoration_requested"]=true;}
        }catch(const std::exception& e){attempt["failure"]=e.what();break;}
        if(!attempt["converged"].get<bool>()&&!attempt.contains("failure"))attempt["failure"]="Newton update budget exhausted";
        saveState(folder/"final_state.h5",unpack(x),bias);
        return std::make_tuple(ArclengthState{x,bias},attempt);
    };
    auto attempt=[&](const ArclengthState& s,double target,const std::string& kind) {
        const auto label="attempt_"+std::to_string(report["attempts"].size())+"_"+kind;
        ArclengthState result=s;Json rec;
        try{std::tie(result,rec)=correct(s,target,label);}catch(const std::exception& e){rec={{"converged",false},{"failure",e.what()},{"target_outer_V",target},{"label",label}};}
        report["attempts"].push_back(rec);save();return std::make_tuple(result,rec);
    };
    auto [current,zeroRec]=attempt({initial,0.},0.,"zero");
    if(!zeroRec["converged"].get<bool>()){report["failure"]="zero-bias coupled initialization failed";save();return 2;}
    double acceptedOuter=0.,previousOuter=0.;ArclengthState previous=current;bool hasPrevious=false;
    report["points"].push_back({{"outer_V",0.},{"inner_V",current.lambda},{"diagnostics",zeroRec["last"]},{"updates",zeroRec["iterations"].size()}});
    saveState(out/"state_0.h5",unpack(current.x),current.lambda);
    std::ofstream iv(out/"iv.csv");iv<<std::setprecision(17)<<"outer_V,inner_V,current_A_per_um,updates,electron_closure,hole_closure\n";
    auto ivrow=[&](double v,const ArclengthState& s,const Json& rec){iv<<v<<','<<s.lambda<<','<<rec["last"]["current_A_per_um"]<<','<<rec["iterations"].size()<<','<<rec["last"]["electron"]["ratio"]<<','<<rec["last"]["hole"]["ratio"]<<'\n';iv.flush();};
    ivrow(0.,current,zeroRec);double step=outerControl.initial;int retries=0;
    while(report["attempts"].size()<static_cast<std::size_t>(attemptBudget)) {
        auto allowed=hasPrevious?diagnostics::guardedOuterStep(step,outerControl,innerControl.volts,current.lambda-previous.lambda,acceptedOuter-previousOuter):std::optional<double>(step);
        if(!allowed){report["failure"]="inner predictor cap requires sub-minimum outer step";break;}
        const double target=std::min(stopVoltage,acceptedOuter+*allowed);ArclengthState prediction=current;
        if(hasPrevious&&extrapolate) {
            const double ratio=(target-acceptedOuter)/(acceptedOuter-previousOuter);
            prediction.x=current.x+ratio*(current.x-previous.x);prediction.lambda=current.lambda+ratio*(current.lambda-previous.lambda);
        }
        auto [candidate,rec]=attempt(prediction,target,"predictor");bool fallback=false;
        if(!rec["converged"].get<bool>()&&hasPrevious&&extrapolate){std::tie(candidate,rec)=attempt(current,target,"accepted_state_retry");fallback=true;}
        if(rec["converged"].get<bool>()) {
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
    save();return report["complete"].get<bool>()?0:3;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
