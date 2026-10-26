// Uninstalled controlled experiment. The physical residual is always the
// production residual; only the trial linearization and finite update vary.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "BvMethodsPopulationNewton.h"

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
    cfg=canonicalizeDeck(cfg);
    const fs::path out=argv[2]; fs::create_directories(out);
    const double inner=m.at("inner_V"), outer=m.at("outer_V");
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
    if(n<=312)throw std::invalid_argument("this BVMethods diagnostic requires its recorded node IDs");
    const auto& u=p.newton.inputScaling.unitSystem();
    auto cm=[&](double a){return u.internalConcentrationToM3(a)/1e6;};
    const double vt=constants::kb*p.newton.temperature_K/constants::q;
    const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*1e-6;
    auto port=solver.makeArclengthContactCurrentFunctional("drain",line);
    std::vector<ArclengthScalarFunctional> ports;
    for(const auto& c:p.mesh.contacts()) ports.push_back(solver.makeArclengthContactCurrentFunctional(c.name,line));
    DDSolution seed;
    if(m.contains("fields")) {
        seed.psi.resize(n);seed.phin=seed.phip=seed.n=seed.p=VectorXd::Zero(n);
        const auto& f=m.at("fields");
        for(int i=0;i<n;++i) seed.psi[i]=f.at("psi_V")[i];
        for(const auto& v:f.at("electron")) seed.n[v[0].get<int>()]=u.m3ToInternalConcentration(v[1].get<double>()*1e6);
        for(const auto& v:f.at("hole")) seed.p[v[0].get<int>()]=u.m3ToInternalConcentration(v[1].get<double>()*1e6);
        for(const auto& t:solver.evaluatePoissonTerms(seed).rows) {
            seed.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;
            seed.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;
        }
    } else seed=readDDSolutionState(m.at("state_file").get<std::string>(),n,p.newton.inputScaling);
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
    const double tau=m.value("fraction_tau",.99);
    const int minIterations=m.value("min_iterations",8);
    if(!(tau>0.&&tau<1.)||minIterations<1||minIterations>budget)throw std::invalid_argument("invalid experiment controls");
    std::vector<double> thetaCandidates=m.value("theta_candidates",std::vector<double>{1.,.5,.25,0.});
    for(double t:thetaCandidates)if(!std::isfinite(t)||t<0.||t>1.)throw std::invalid_argument("theta outside [0,1]");
    std::sort(thetaCandidates.begin(),thetaCandidates.end(),std::greater<double>());
    thetaCandidates.erase(std::unique(thetaCandidates.begin(),thetaCandidates.end()),thetaCandidates.end());
    if(thetaCandidates.empty()||thetaCandidates.back()!=0.)throw std::invalid_argument("theta=0 baseline required");
    auto coordinates=[&](const DDSolution& s) {
        VectorXd en(n),hp(n),es(n),hs(n);
        for(int i=0;i<n;++i){en[i]=cm(s.n[i]);hp[i]=cm(s.p[i]);es[i]=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.n[i],nc[i]))/vt:0.;hs[i]=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.p[i],nv[i]))/vt:0.;}
        return diagnostics::populationCoordinates(s.psi,en,hp,es,hs,active,v0,vt,inner);
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
    const auto initialCoordinates=coordinates(seed);
    const detail::BvMethodsQfState initialState{seed.psi,seed.phin,seed.phip,seed.n,seed.p};
    const auto lt=element.transport(initialState,true,true,true),ft=element.transport(initialState,true,true,true,true);
    const auto av0=element.evaluate(initialState,true);auto tp=port.linearize(initial,inner);
    tp.stateDerivative+=line*continuity*((ft.jacobian-lt.jacobian).transpose()*weights);
    const SparseMatrixd refDevice=sys.jacobian(initial,inner)+project(ft.jacobian-lt.jacobian-av0.jacobian);
    const SparseMatrixd reference=border(refDevice,sys.parameterDerivative(initial,inner),tp)*initialCoordinates.qfFromScaledPopulation;
    const VectorXd fixedRows=diagnostics::equilibratePopulationSystem(reference).rowWeights;
    std::ofstream scales(out/"fixed_residual_rows.csv");scales<<std::setprecision(17)<<"row,weight\n";for(int i=0;i<=d;++i)scales<<i<<','<<fixedRows[i]<<'\n';scales.close();
    Json report={{"backend",LinearSolver::selectedBackend()},{"dimension",d+1},{"potential_scale_V",v0},
        {"continuity_scale",continuity},{"current_per_normalized_flux",constants::q*line*continuity},
        {"tau",tau},{"minimum_diagnostic_iterations",minIterations},{"maximum_iterations",budget},
        {"theta_candidates",thetaCandidates},{"equation_tolerance",eqtol},{"load_tolerance_V",loadtol},
        {"outer_V",outer},{"initial_inner_V",inner},{"arms",Json::array()}};
    auto save=[&](){std::ofstream(out/"summary.json")<<report.dump(2)<<'\n';};
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
    const VectorXd initialF=fullResidual(initial,inner);
    // Independent physical-population finite differences, including psi changes
    // at fixed n/p. Do not perturb the already transformed matrix itself.
    Json fd=Json::array();
    const SparseMatrixd initialPopulationJac=border(sys.jacobian(initial,inner),sys.parameterDerivative(initial,inner),port.linearize(initial,inner))*initialCoordinates.qfFromScaledPopulation;
    for(int block:{0,1,2})for(int node:{308,310,311}) {
        if(!active[node])continue;
        VectorXd direction=VectorXd::Zero(d+1);
        const double amount=block==0?vt:block==1?cm(seed.n[node]):cm(seed.p[node]);
        const double scale=block==0?initialCoordinates.psiScale_V[node]:block==1?initialCoordinates.electronScale_cm3[node]:initialCoordinates.holeScale_cm3[node];
        direction[block*n+node]=amount/scale;const VectorXd analytic=initialPopulationJac*direction;
        for(double eps:{1e-3,1e-5}) {
            auto value=[&](double sign) {
                DDSolution state=seed;state.packedState.resize(0);state.packedLow.resize(0);
                state.electronQfReference.resize(0);state.holeQfReference.resize(0);
                state.phinIncrement.resize(0);state.phipIncrement.resize(0);
                if(block==0)state.psi[node]+=sign*eps*amount;
                else if(block==1)state.n[node]*=1.+sign*eps;
                else state.p[node]*=1.+sign*eps;
                for(const auto& t:solver.evaluatePoissonTerms(state).rows){state.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;state.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
                return fullResidual(solver.packArclengthState(state),inner);
            };
            const VectorXd numerical=(value(1.)-value(-1.))/(2*eps);
            Json entries=Json::array();for(int row:{308,310,311,312})for(int b:{1,2}) {
                const int i=b*n+row;entries.push_back({{"row",i},{"analytic",analytic[i]},{"central",numerical[i]},
                    {"relative_difference",std::abs(analytic[i]-numerical[i])/std::max({std::abs(analytic[i]),std::abs(numerical[i]),1e-300})}});
            }
            fd.push_back({{"column_block",block},{"column_node",node},{"relative_step",eps},{"entries",entries}});
        }
    }
    std::ofstream(out/"population_fd.json")<<fd.dump(2)<<'\n';
    std::ofstream initialFile(out/"initial_rhs.csv");initialFile<<std::setprecision(17)<<"row,rhs\n";for(int i=0;i<=d;++i)initialFile<<i<<','<<-initialF[i]<<'\n';initialFile.close();
    struct Arm{const char* name;double theta;bool projection;int adaptive;bool noAv;};
    const std::array<Arm,7> arms{{{"scaled_full_projection",1.,true,false,false},{"scaled_frozen_projection",0.,true,false,false},
        {"scaled_full_fraction",1.,false,false,false},{"scaled_frozen_fraction",0.,false,false,false},
        {"scaled_adaptive_fraction",0.,false,true,false},{"scaled_combo_fraction",0.,false,false,true},
        {"scaled_target_adaptive_fraction",0.,false,2,false}}};
    struct Candidate {bool valid=false;VectorXd z,dx,dn,dp,x,F;double bias=0.,theta=0.,alpha=1.,merit=0.,rms=0.;Json log;};
    for(const auto& arm:arms) {
        if(m.contains("selected_arms")){const auto names=m.at("selected_arms").get<std::vector<std::string>>();if(std::find(names.begin(),names.end(),arm.name)==names.end())continue;}
        const fs::path folder=out/arm.name;fs::create_directories(folder);
        VectorXd x=initial;double bias=inner;
        Json a={{"name",arm.name},{"projection_1pct",arm.projection},{"adaptive",arm.adaptive},{"AV_derivatives",!arm.noAv},
            {"status","budget_exhausted"},{"initial",inspect(x,bias)},{"iterations",Json::array()}};
        LinearSolver linear;
        for(int iteration=1;iteration<=budget;++iteration)try {
            const auto s=unpack(x);const auto c=coordinates(s);const VectorXd F=fullResidual(x,bias);
            if(iteration==1&&!(F.array()==initialF.array()).all())throw std::runtime_error("different initial RHS");
            const detail::BvMethodsQfState state{s.psi,s.phin,s.phip,s.n,s.p};
            const auto live=element.transport(state,true,true,true),frozen=element.transport(state,true,true,true,true);
            if(!(live.source.array()==frozen.source.array()).all())throw std::runtime_error("frozen tangent changed transport value");
            const SparseMatrixd mu=frozen.jacobian-live.jacobian;
            const SparseMatrixd av=element.evaluate(state,true).jacobian;
            const auto fullJac=sys.jacobian(x,bias);const auto fv=sys.parameterDerivative(x,bias);const auto terminal=port.linearize(x,bias);
            const SparseMatrixd full=border(fullJac,fv,terminal);
            const VectorXd portMu=line*continuity*(mu.transpose()*weights);
            auto calculate=[&](double theta) {
                Candidate q;q.theta=theta;q.log={{"theta",theta},{"valid",false}};
                try {
                    SparseMatrixd correction=(1.-theta)*mu;if(arm.noAv)correction-=av;
                    const SparseMatrixd jac=fullJac+project(correction);auto pt=terminal;pt.stateDerivative+=(1.-theta)*portMu;
                    const SparseMatrixd jq=border(jac,fv,pt),jp=jq*c.qfFromScaledPopulation;
                    const std::string tag="theta"+std::to_string(theta);
                    if(iteration==1){matrixFile(folder/(tag+"_qf.csv"),jq);matrixFile(folder/(tag+"_population.csv"),jp);matrixFile(folder/"coordinate_map.csv",c.qfFromScaledPopulation);}
                    const auto eq=diagnostics::equilibratePopulationSystem(jp);
                    q.z=eq.columnWeights.cwiseProduct(linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(F)).eval()));
                    auto defect=linear_refinement::highPrecisionDefect(jp,F,q.z);
                    if(defect.backwardMax>1e-12){q.z+=eq.columnWeights.cwiseProduct(linear.solve(eq.matrix,(-eq.rowWeights.cwiseProduct(defect.values)).eval()));defect=linear_refinement::highPrecisionDefect(jp,F,q.z);}
                    q.log["population_linear_backward_error"]=defect.backwardMax;
                    if(defect.backwardMax>1e-10)throw std::runtime_error("population linear solve defect too large");
                    q.dx=c.qfFromScaledPopulation*q.z;
                    if(iteration==1){std::ofstream raw(folder/(tag+"_direction.csv"));raw<<std::setprecision(17)<<"row,scaled_population,qf\n";for(int i=0;i<=d;++i)raw<<i<<','<<q.z[i]<<','<<q.dx[i]<<'\n';}
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
                    if(!arm.projection)fraction=diagnostics::populationFractionToBoundary(en,hp,q.dn,q.dp,active,tau);
                    q.alpha=fraction.alpha;q.bias=bias+q.alpha*q.dx[d];
                    DDSolution next=s;next.packedState.resize(0);next.packedLow.resize(0);next.electronQfReference.resize(0);next.holeQfReference.resize(0);next.phinIncrement.resize(0);next.phipIncrement.resize(0);
                    int negativeN=0,negativeP=0;double relativeN=0.,relativeP=0.;
                    for(int i=0;i<n;++i){next.psi[i]+=q.alpha*q.dx[i]*v0;
                        if(!active[i]){next.phin[i]+=q.alpha*q.dx[n+i]*v0;next.phip[i]+=q.alpha*q.dx[2*n+i]*v0;continue;}
                        double nn=en[i]+q.alpha*q.dn[i],pp=hp[i]+q.alpha*q.dp[i];
                        negativeN+=en[i]+q.dn[i]<=0.;negativeP+=hp[i]+q.dp[i]<=0.;
                        if(arm.projection){if(nn<=0.)nn=.01*en[i];if(pp<=0.)pp=.01*hp[i];}
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
            Candidate chosen=calculate(arm.theta);Json trials=Json::array();trials.push_back(chosen.log);
            if(arm.adaptive)for(double t:thetaCandidates)if(t!=0.) {
                Candidate candidate=calculate(t);trials.push_back(candidate.log);
                const bool withinTargets=arm.adaptive==2&&candidate.valid&&diagnostics::restorationWithinTargets(candidate.F.head(d).lpNorm<Eigen::Infinity>(),candidate.F[d],candidate.rms,eqtol,loadtol);
                candidate.log["within_existing_targets"]=withinTargets;
                trials.back()=candidate.log;
                if(candidate.valid&&(!chosen.valid||withinTargets||diagnostics::restorationNoWorse(candidate.merit,candidate.rms,chosen.merit,chosen.rms))){chosen=std::move(candidate);break;}
            }
            if(!chosen.valid){a["iterations"].push_back({{"iteration",iteration},{"trials",trials},{"error","no valid population direction"}});a["status"]="invalid_or_failed_candidate";break;}
            Json row={{"iteration",iteration},{"theta",chosen.theta},{"alpha",chosen.alpha},{"raw_RMS",chosen.rms},
                {"before_fixed_merit",fixedRows.cwiseProduct(F).lpNorm<Eigen::Infinity>()},{"trials",trials},
                {"diagnostics",inspect(chosen.x,chosen.bias)}};
            writeDDSolutionState(folder/("step"+std::to_string(iteration)+".h5"),unpack(chosen.x),p.newton.inputScaling);
            x=chosen.x;bias=chosen.bias;a["iterations"].push_back(row);
            if(row["diagnostics"]["original_gates_pass"].get<bool>()) {
                if(!a.contains("first_original_gate_iteration"))a["first_original_gate_iteration"]=iteration;
                if(chosen.rms<1.){if(!a.contains("first_original_and_RMS_iteration"))a["first_original_and_RMS_iteration"]=iteration;
                    if(iteration>=minIterations){a["status"]="bounded_diagnostics_passed";break;}}
            }
            std::cout<<arm.name<<' '<<iteration<<" theta="<<chosen.theta<<" alpha="<<chosen.alpha<<" RMS="<<chosen.rms<<" F="<<row["diagnostics"]["equation_inf"]<<std::endl;
        }catch(const std::exception& e){a["status"]="execution_error";a["error"]=e.what();break;}
        report["arms"].push_back(a);save();
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
