// Uninstalled checkpoint restoration control. Derived from the cold-start diagnostic.
// The physical residual is always the
// production residual; only the trial linearization and finite update vary.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "BvMethodsPopulationNewton.h"
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

    const auto innerControl=diagnostics::resolveInnerStepLimit(cfg,m);
    const auto outerControl=diagnostics::resolveOuterStepControl(cfg,m);
    const std::string predictorMode=m.value("inner_predictor_mode",std::string{"enforce"});
    if(predictorMode!="enforce"&&predictorMode!="report_only")
        throw std::invalid_argument("inner_predictor_mode must be enforce or report_only");
    const bool reportPredictorOnly=predictorMode=="report_only";
    cfg["sweep"].erase("initial_state_file");cfg["solver"]["warm_start"]=false;
    cfg=canonicalizeDeck(cfg);
    const fs::path out=argv[2]; fs::create_directories(out);
    const double inner=0.;double outer=m.at("outer_V");
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
    DDSolution currentSeed=readDDSolutionState(m.at("current_state_file").get<std::string>(),n,p.newton.inputScaling);
    DDSolution previousSeed=readDDSolutionState(m.at("previous_state_file").get<std::string>(),n,p.newton.inputScaling);
    const VectorXd currentPacked=solver.packArclengthState(currentSeed), previousPacked=solver.packArclengthState(previousSeed);
    const double ratio=(outer-m.at("current_outer_V").get<double>())/
        (m.at("current_outer_V").get<double>()-m.at("previous_outer_V").get<double>());
    const double initialBias=m.at("current_inner_V").get<double>()+ratio*
        (m.at("current_inner_V").get<double>()-m.at("previous_inner_V").get<double>());
    const VectorXd initial=currentPacked+ratio*(currentPacked-previousPacked);
    DDSolution seed=unpack(initial);
    saveState(out/"predicted_initial.h5",seed,initialBias);
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
    Json report={{"schema","vela.restoration_checkpoint.v1"},{"backend",LinearSolver::selectedBackend()},
        {"dimension",d+1},{"outer_V",outer},{"initial_inner_V",initialBias},{"potential_scale_V",v0},
        {"equation_tolerance",eqtol},{"load_tolerance_V",loadtol},{"fraction_tau",tau},
        {"corrector_voltage_update_limit_V",voltageUpdateLimit},{"max_updates_per_arm",budget},
        {"physical_RHS_unchanged",true},{"arms",Json::array()},{"complete",false}};
    auto save=[&](){diagnostics::writeExperimentJson(out/"summary.json",report);};save();
    struct Candidate {bool valid=false;VectorXd z,dx,dn,dp,x,F;double bias=0.,theta=0.,alpha=1.,merit=0.,rms=0.;Json log;};
    auto correct=[&](const ArclengthState& start,double target,const std::string& label,
                     const std::string& mode,int maxUpdates,const VectorXd* inheritedRows) {
        outer=target;VectorXd x=start.x;double bias=start.lambda;
        const fs::path folder=out/label;fs::create_directories(folder);
        bool bootstrap=false;Json attempt={{"target_outer_V",target},{"label",label},{"iterations",Json::array()},{"converged",false}};
        const auto first=unpack(x);const auto co=coordinates(first,bias);
        const detail::BvMethodsQfState refState{first.psi,first.phin,first.phip,first.n,first.p};
        const auto lt=element.transport(refState,true,true,true),ft=element.transport(refState,true,true,true,true),av0=element.evaluate(refState,true);
        auto pt=port.linearize(x,bias);pt.stateDerivative+=line*continuity*((ft.jacobian-lt.jacobian).transpose()*weights);
        const SparseMatrixd ref=border(sys.jacobian(x,bias)+project(ft.jacobian-lt.jacobian-av0.jacobian),sys.parameterDerivative(x,bias),pt)*co.qfFromScaledPopulation;
        fixedRows=inheritedRows?*inheritedRows:diagnostics::equilibratePopulationSystem(ref).rowWeights;
        {std::ofstream f(folder/"fixed_rows.csv");f<<std::setprecision(17)<<"row,weight\n";
         for(int i=0;i<=d;++i)f<<i<<','<<fixedRows[i]<<'\n';}
        attempt["initial"]=inspect(x,bias);
        // Preserve the original residual-only zero-bias initialization rule.
        if(target==0.&&attempt["initial"]["original_gates_pass"].get<bool>()) {
            attempt["converged"]=true;attempt["reason"]="original residual-only zero anchor";
            attempt["last"]=attempt["initial"];return std::make_tuple(ArclengthState{x,bias},attempt);
        }
        LinearSolver linear;
        for(int iteration=1;iteration<=maxUpdates;++iteration)try {
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
                    Json blocks=Json::array();std::array<long double,4> sums{};std::array<int,4> counts{};
                    std::array<double,4> maximum{};std::array<Json,4> worst;
                    auto metric=[&](int b,int node,double after,double delta,double reference) {
                        const long double value=1e5L*static_cast<long double>(delta)/(std::abs(static_cast<long double>(after))+reference);
                        sums[b]+=value*value;++counts[b];
                        if(std::abs(value)>maximum[b]){maximum[b]=std::abs(value);
                            worst[b]={{"node",node},{"after",after},{"raw_delta",delta},{"reference",reference},
                                {"normalized_magnitude",maximum[b]},{"psi_V",check.psi[node]},
                                {"fn_V",check.phin[node]},{"fp_V",check.phip[node]},
                                {"electron_cm3",cm(check.n[node])},{"hole_cm3",cm(check.p[node])}};}
                    };
                    for(int i=0;i<n;++i){metric(0,i,q.x[i]*v0,q.dx[i]*v0,.025852);
                        if(active[i]&&efree[i])metric(1,i,cm(check.n[i]),q.dn[i],1e10);
                        if(active[i]&&hfree[i])metric(2,i,cm(check.p[i]),q.dp[i],1e10);}
                    for(const auto& contact:p.mesh.contacts())if(contact.name=="drain")
                        for(int i:contact.node_ids)metric(3,i,q.bias,q.dx[d],.025852);
                    long double total=0.;int metricCount=0;
                    for(int b=0;b<4;++b){total+=sums[b];metricCount+=counts[b];
                        blocks.push_back({{"block",b},{"count",er.counts[b]},{"RMS",er.blockRms[b]},
                            {"independent_reconstruction_RMS",counts[b]?std::sqrt(sums[b]/counts[b]):0.L},{"worst",worst[b]}});}
                    q.log["RMS_blocks"]=blocks;
                    q.log["independent_RMS_relative_difference"]=std::abs(static_cast<double>(std::sqrt(total/metricCount))-er.rms)/std::max(1.,er.rms);
                    if(dumpSystems&&iteration==1){
                        std::ofstream f(folder/(tag+"_rhs.csv"));f<<std::setprecision(17)<<"row,F\n";
                        for(int k=0;k<=d;++k)f<<k<<','<<F[k]<<'\n';
                    }
                    if(mode=="full"&&iteration==1){
                        Json checks=Json::array();const VectorXd analytic=jq*q.dx;
                        for(double h:{.1,.01,.001}){
                            const VectorXd plus=fullResidual((x+h*q.dx.head(d)).eval(),bias+h*q.dx[d]);
                            const VectorXd minus=fullResidual((x-h*q.dx.head(d)).eval(),bias-h*q.dx[d]);
                            const VectorXd fd=(plus-minus)/(2*h),diff=fd-analytic;
                            Eigen::Index worstRow=0;diff.cwiseAbs().maxCoeff(&worstRow);
                            Json bc=Json::array();for(int b=0;b<4;++b){const int offset=b*n,len=b==3?1:n;
                                bc.push_back({{"block",b},{"relative_l2",diff.segment(offset,len).norm()/std::max(analytic.segment(offset,len).norm(),1e-300)},
                                    {"absolute_inf",diff.segment(offset,len).lpNorm<Eigen::Infinity>()},
                                    {"analytic_inf",analytic.segment(offset,len).lpNorm<Eigen::Infinity>()}});}
                            checks.push_back({{"h",h},{"relative_l2",diff.norm()/std::max(analytic.norm(),1e-300)},
                                {"relative_inf",diff.lpNorm<Eigen::Infinity>()/std::max(analytic.lpNorm<Eigen::Infinity>(),1e-300)},
                                {"worst_row",worstRow},{"blocks",bc}});
                            std::ofstream f(folder/("jv_h"+std::to_string(h)+".csv"));f<<std::setprecision(17)<<"row,analytic,finite_difference\n";
                            for(int k=0;k<=d;++k)f<<k<<','<<analytic[k]<<','<<fd[k]<<'\n';
                        }
                        q.log["physical_Jv"]=checks;
                        // A near-null Newton direction has an unusually small
                        // Jv signal. Audit independent physical columns too.
                        Json columns=Json::array();
                        for(int node:{1129,1239})for(int block:{0,1,2}) {
                            VectorXd direction=VectorXd::Zero(d+1);
                            direction[block*n+node]=vt/v0;
                            const VectorXd expected=jq*direction;
                            for(double h:{.01,.001,.0001}) {
                                const VectorXd plus=fullResidual((x+h*direction.head(d)).eval(),bias);
                                const VectorXd minus=fullResidual((x-h*direction.head(d)).eval(),bias);
                                const VectorXd numerical=(plus-minus)/(2*h),error=numerical-expected;
                                Eigen::Index worst=0;error.cwiseAbs().maxCoeff(&worst);
                                columns.push_back({{"node",node},{"block",block},{"h",h},{"physical_column_scale_V",vt},
                                    {"relative_l2",error.norm()/std::max(expected.norm(),1e-300)},
                                    {"relative_inf",error.lpNorm<Eigen::Infinity>()/std::max(expected.lpNorm<Eigen::Infinity>(),1e-300)},
                                    {"absolute_inf",error.lpNorm<Eigen::Infinity>()},{"analytic_inf",expected.lpNorm<Eigen::Infinity>()},
                                    {"worst_row",worst}});
                            }
                        }
                        q.log["physical_columns"]=columns;
                    }
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
            if(mode=="current"&&!bootstrap) {
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
            if(mode=="current"&&bootstrap){chosen=calculate(0.,false);trials.push_back(chosen.log);}
            if(mode=="full"||mode=="av_only") {chosen=calculate(mode=="full"?1.:0.,true);trials.push_back(chosen.log);}
            if(!chosen.valid){attempt["failure"]="no valid population candidate";attempt["failed_trials"]=trials;break;}
            auto diag=inspect(chosen.x,chosen.bias);
            const bool joint=diag["original_gates_pass"].get<bool>()&&chosen.rms<1.;
            Json row={{"iteration",iteration},{"phase",bootstrap?"combination_bootstrap":"full_AV_adaptive_mobility"},
                {"theta",chosen.theta},{"alpha",chosen.alpha},{"raw_RMS",chosen.rms},{"diagnostics",diag},{"trials",trials},
                {"RMS_blocks",chosen.log["RMS_blocks"]},{"before_merit",fixedRows.cwiseProduct(F).lpNorm<Eigen::Infinity>()}};
            x=chosen.x;bias=chosen.bias;attempt["iterations"].push_back(row);attempt["last"]=diag;
            diagnostics::writeExperimentJson(folder/"attempt.json",attempt);
            saveState(folder/("step"+std::to_string(iteration)+".h5"),unpack(x),bias);
            std::cout<<label<<" iteration="<<iteration<<" RMS="<<chosen.rms<<" F="<<diag["equation_inf"]<<" bootstrap="<<bootstrap<<std::endl;
            if(joint&&!bootstrap){if(!attempt.contains("first_joint_restored_iteration"))attempt["first_joint_restored_iteration"]=iteration;
                if(iteration>=m.value("minimum_diagnostic_iterations",8)&&label!="replay"){
                    attempt["converged"]=true;attempt["reason"]="original gates and raw population RMS after AV restoration";break;}}
            if(joint&&bootstrap){bootstrap=false;attempt["AV_restoration_requested"]=true;}
        }catch(const std::exception& e){attempt["failure"]=e.what();break;}
        if(!attempt["converged"].get<bool>()&&!attempt.contains("failure"))attempt["failure"]="Newton update budget exhausted";
        saveState(folder/"final_state.h5",unpack(x),bias);
        return std::make_tuple(ArclengthState{x,bias},attempt);
    };
    const int replayUpdates=m.value("replay_updates",10);
    auto [checkpoint,replay]=correct({initial,initialBias},outer,"replay","current",replayUpdates,nullptr);
    report["replay"]=replay;const VectorXd sharedRows=fixedRows;
    const double expected=m.at("expected_replay_RMS");
    const double actual=replay.at("iterations").back().at("raw_RMS");
    report["replay_RMS_relative_difference"]=std::abs(actual-expected)/std::max(1.,expected);
    if(report["replay_RMS_relative_difference"].get<double>()>1e-8)throw std::runtime_error("checkpoint replay differs from original trajectory");
    report["checkpoint"]=inspect(checkpoint.x,checkpoint.lambda);saveState(out/"shared_checkpoint.h5",unpack(checkpoint.x),checkpoint.lambda);save();
    for(const std::string mode:{"current","full","av_only"}){
        auto [last,rec]=correct(checkpoint,outer,mode,mode,budget,&sharedRows);
        rec["mode"]=mode;report["arms"].push_back(rec);save();
    }
    report["complete"]=true;save();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
