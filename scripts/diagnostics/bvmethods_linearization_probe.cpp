// Uninstalled controlled experiment. The physical residual is always the
// production residual; only the trial linearization and finite update vary.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/simulation/CoupledLoadLine.h"

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
    Json report={{"backend",LinearSolver::selectedBackend()},{"dimension",d+1},
        {"potential_scale_V",v0},{"thermal_voltage_V",vt},{"continuity_scale",continuity},
        {"outer_V",outer},{"initial_inner_V",inner},{"equation_tolerance",eqtol},
        {"load_tolerance_V",loadtol},{"update_RMS_enforced",false},{"arms",Json::array()}};
    auto save=[&](){std::ofstream(out/"summary.json")<<report.dump(2)<<'\n';};
    struct Arm {const char* name;bool freeze,noav,population;};
    const std::array<Arm,6> arms{{{"full_qf",false,false,false},{"approx_qf",true,true,false},
        {"full_population",false,false,true},{"approx_population",true,true,true},
        {"frozen_mobility_population",true,false,true},{"noav_derivative_population",false,true,true}}};
    const VectorXd initialRhs=sys.residual(initial,inner);
    for(const auto& arm:arms) {
        Json result={{"name",arm.name},{"freeze_mobility_derivatives",arm.freeze},
            {"omit_AV_derivatives",arm.noav},{"population_candidate",arm.population},
            {"status","budget_exhausted"},{"iterations",Json::array()}};
        VectorXd x=initial;double bias=inner;
        const fs::path folder=out/arm.name;fs::create_directories(folder);
        auto inspect=[&](const VectorXd& state,double voltage) {
            const VectorXd f=sys.residual(state,voltage); const auto s=unpack(state);
            if(!f.allFinite()) throw std::runtime_error("nonfinite physical residual");
            const auto g=gates(state,voltage);double kcl=0.;
            for(const auto& terminal:ports) kcl+=terminal.value(state,voltage);
            const double id=port.value(state,voltage),load=voltage+resistance*id-outer;
            bool positive=true;double minN=1e300,minP=1e300;
            for(int i=0;i<n;++i) if(active[i]) {
                positive=positive&&s.n[i]>0.&&s.p[i]>0.&&std::isfinite(s.n[i])&&std::isfinite(s.p[i]);
                minN=std::min(minN,cm(s.n[i]));minP=std::min(minP,cm(s.p[i]));
            }
            const bool physical=positive&&std::isfinite(kcl)&&std::abs(kcl)<=std::max(1e-12,1e-5*std::abs(id))&&
                (!g.carrierRows.enforced||g.carrierRows.satisfied)&&(!g.globalContinuity.enforced||g.globalContinuity.satisfied);
            Eigen::Index worst;const double inf=f.cwiseAbs().maxCoeff(&worst);
            return Json{{"equation_inf",inf},{"equation_worst_row",worst},{"load_V",load},
                {"psi_inf",g.blockNorms.psi},{"electron_inf",g.blockNorms.phin},{"hole_inf",g.blockNorms.phip},
                {"electron_closure",g.globalContinuity.electron.ratio},{"hole_closure",g.globalContinuity.hole.ratio},
                {"carrier_rows_enforced",g.carrierRows.enforced},{"carrier_rows_satisfied",g.carrierRows.satisfied},
                {"carrier_row_max_ratio",g.carrierRows.maxRatio},
                {"closure_enforced",g.globalContinuity.enforced},{"closure_satisfied",g.globalContinuity.satisfied},
                {"current_A_per_um",id},{"KCL_A_per_um",kcl},{"positive_finite",positive},
                {"min_electron_cm3",minN},{"min_hole_cm3",minP},{"inner_V",voltage},
                {"original_gates_pass",inf<=eqtol&&std::abs(load)<=loadtol&&physical}};
        };
        result["initial"]=inspect(x,bias);
        for(int iteration=1;iteration<=budget;++iteration) {
            Json row={{"iteration",iteration}};
            try {
                const auto s=unpack(x);
                const detail::BvMethodsQfState state{s.psi,s.phin,s.phip,s.n,s.p};
                const auto live=element.transport(state,true,true,true);
                const auto frozen=element.transport(state,true,true,true,true);
                if(!(live.source.array()==frozen.source.array()).all()) throw std::runtime_error("mobility freeze changed RHS");
                const auto av=element.evaluate(state,true);
                const auto ledger=solver.evaluateCarrierTermDiagnostics(s,true).rows;
                double transportError=0.,avError=0.,transportNorm=0.,avNorm=0.;
                for(const auto& t:ledger) {
                    const int i=t.nodeId;
                    transportError=std::max({transportError,std::abs(live.source[n+i]-t.electronFlux),std::abs(live.source[2*n+i]-t.holeFlux)});
                    avError=std::max(avError,std::abs(av.source[n+i]-t.electronImpact));
                    transportNorm=std::max({transportNorm,std::abs(t.electronFlux),std::abs(t.holeFlux)});
                    avNorm=std::max(avNorm,std::abs(t.electronImpact));
                }
                if(transportError>1e-12*std::max(transportNorm,1e-300)||avError>1e-12*std::max(avNorm,1e-300))
                    throw std::runtime_error("diagnostic element differs from production ledger");
                const auto terminal=port.linearize(x,bias);
                const VectorXd independentPortGradient=line*continuity*(live.jacobian.transpose()*weights);
                const double portError=(independentPortGradient-terminal.stateDerivative).norm()/std::max(terminal.stateDerivative.norm(),1e-300);
                if(portError>1e-12) throw std::runtime_error("terminal gradient mismatch");
                SparseMatrixd correction(d,d);if(arm.freeze) correction=frozen.jacobian-live.jacobian;
                if(arm.noav) correction-=av.jacobian;
                std::vector<Eigen::Triplet<double>> entries;
                for(int k=0;k<correction.outerSize();++k) for(SparseMatrixd::InnerIterator v(correction,k);v;++v)
                    if((v.row()>=n&&v.row()<2*n&&efree[v.row()-n])||(v.row()>=2*n&&hfree[v.row()-2*n]))
                        entries.emplace_back(v.row(),v.col(),v.value());
                SparseMatrixd projected(d,d);projected.setFromTriplets(entries.begin(),entries.end());
                const SparseMatrixd fullJac=sys.jacobian(x,bias);
                const SparseMatrixd jac=fullJac+projected;
                auto controlled=sys;controlled.jacobian=[&](const VectorXd&,double){return jac;};
                auto approximatePort=port;
                auto linearPort=terminal;
                if(arm.freeze) linearPort.stateDerivative+=line*continuity*((frozen.jacobian-live.jacobian).transpose()*weights);
                approximatePort.linearize=[&](const VectorXd&,double){return linearPort;};
                const VectorXd beforeRhs=sys.residual(x,bias);
                if(iteration==1&&!(beforeRhs.array()==initialRhs.array()).all()) throw std::runtime_error("arm RHS differs at checkpoint");
                if(iteration==1) {
                    const VectorXd fv=sys.parameterDerivative(x,bias);
                    std::vector<Eigen::Triplet<double>> border;
                    for(int k=0;k<jac.outerSize();++k)for(SparseMatrixd::InnerIterator v(jac,k);v;++v)border.emplace_back(v.row(),v.col(),v.value());
                    for(int i=0;i<d;++i){if(fv[i]!=0.)border.emplace_back(i,d,fv[i]);if(linearPort.stateDerivative[i]!=0.)border.emplace_back(d,i,resistance*linearPort.stateDerivative[i]);}
                    border.emplace_back(d,d,1.+resistance*linearPort.parameterDerivative);
                    SparseMatrixd bordered(d+1,d+1);bordered.setFromTriplets(border.begin(),border.end());matrixFile(folder/"matrix.csv",bordered);
                    std::ofstream rhs(folder/"rhs.csv");rhs<<std::setprecision(17)<<"row,rhs\n";
                    for(int i=0;i<d;++i)rhs<<i<<','<<-beforeRhs[i]<<'\n';rhs<<d<<','<<-(bias+resistance*terminal.value-outer)<<'\n';
                }
                CoupledLoadLineConfig c;c.maxIterations=1;c.equationTolerance=1e-30;c.loadLineTolerance_V=1e-30;
                c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.linearizationAuditEnabled=true;
                c.applyDeviceUpdateLimit=false;c.maxVoltageUpdate_V=1e100;
                const auto solved=CoupledLoadLineNewton(controlled,port.value,c,approximatePort).solve({x,bias},outer,resistance);
                if(solved.trace.empty()||solved.trace[0].rawStateUpdate.size()!=d)throw std::runtime_error("raw linear direction missing: "+solved.failureReason);
                const auto& t=solved.trace[0];const VectorXd delta=t.rawStateUpdate;const double dv=t.rawVoltageUpdate_V;
                // Retain the solved direction even if its finite QF candidate
                // subsequently leaves the positive finite population domain.
                if(iteration==1){std::ofstream v(folder/"raw_direction.csv");v<<std::setprecision(17)<<"row,value\n";for(int i=0;i<d;++i)v<<i<<','<<delta[i]<<'\n';v<<d<<','<<dv<<'\n';}
                const VectorXd trueLinearResidual=beforeRhs+fullJac*delta+sys.parameterDerivative(x,bias)*dv;
                row["true_full_linearized_equation_inf"]=trueLinearResidual.lpNorm<Eigen::Infinity>();
                row["true_full_linearized_load_V"]=bias+resistance*terminal.value-outer+dv+
                    resistance*(terminal.stateDerivative.dot(delta)+terminal.parameterDerivative*dv);
                if(iteration==1) {
                    std::ofstream f(folder/"fields.csv");f<<std::setprecision(17)<<"node,psi_V,fn_V,fp_V,n_cm3,p_cm3,dn_dV_cm3,dp_dV_cm3\n";
                    for(int i=0;i<n;++i) {
                        const double dn=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.n[i],nc[i]))/vt:0.;
                        const double dp=active[i]?cm(experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.p[i],nv[i]))/vt:0.;
                        f<<i<<','<<s.psi[i]<<','<<s.phin[i]<<','<<s.phip[i]<<','<<cm(s.n[i])<<','<<cm(s.p[i])<<','<<dn<<','<<dp<<'\n';
                    }
                }
                VectorXd candidate=x+delta;double roundtrip=0.;int negativeN=0,negativeP=0;
                double maxDn=0.,maxDp=0.;
                if(arm.population) {
                    DDSolution next=s;next.packedState.resize(0);next.packedLow.resize(0);
                    next.electronQfReference.resize(0);next.holeQfReference.resize(0);
                    next.phinIncrement.resize(0);next.phipIncrement.resize(0);
                    for(int i=0;i<n;++i) {
                        next.psi[i]+=delta[i]*v0;
                        if(active[i]) {
                            const double dn=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.n[i],nc[i])*(delta[i]-delta[n+i])*v0/vt;
                            const double dp=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(s.p[i],nv[i])*(delta[2*n+i]-delta[i])*v0/vt;
                            maxDn=std::max(maxDn,std::abs(dn/s.n[i]));maxDp=std::max(maxDp,std::abs(dp/s.p[i]));
                            negativeN+=s.n[i]+dn<=0.;negativeP+=s.p[i]+dp<=0.;
                            next.n[i]=s.n[i]+dn>0.?s.n[i]+dn:.01*s.n[i];next.p[i]=s.p[i]+dp>0.?s.p[i]+dp:.01*s.p[i];
                        }
                    }
                    for(const auto& v:solver.evaluatePoissonTerms(next).rows){next.phin[v.nodeId]=v.electronRequiredQuasiFermiPotential;next.phip[v.nodeId]=v.holeRequiredQuasiFermiPotential;}
                    candidate=solver.packArclengthState(next);const auto back=unpack(candidate);
                    for(int i=0;i<n;++i)if(active[i])roundtrip=std::max({roundtrip,std::abs(back.n[i]/next.n[i]-1.),std::abs(back.p[i]/next.p[i]-1.)});
                }
                row["raw_dpsi_inf_V"]=delta.head(n).lpNorm<Eigen::Infinity>()*v0;
                row["raw_dfn_inf_V"]=delta.segment(n,n).lpNorm<Eigen::Infinity>()*v0;
                row["raw_dfp_inf_V"]=delta.tail(n).lpNorm<Eigen::Infinity>()*v0;
                row["raw_inner_update_V"]=dv;row["linear_backward_error"]=t.componentwiseBackwardError;
                row["negative_electron_candidates"]=negativeN;row["negative_hole_candidates"]=negativeP;
                row["max_relative_population_tangent_n"]=maxDn;row["max_relative_population_tangent_p"]=maxDp;
                row["population_roundtrip_relative"]=roundtrip;row["port_gradient_relative_check"]=portError;
                row["transport_ledger_relative_check"]=transportError/std::max(transportNorm,1e-300);
                row["AV_ledger_relative_check"]=avError/std::max(avNorm,1e-300);
                row["diagnostics"]=inspect(candidate,bias+dv);
                const auto error=updateError(x,delta,candidate,bias+dv,dv);
                row["raw_tangent_update_RMS"]=error.rms;row["update_RMS_valid"]=error.valid;
                row["accepted_delta_fp_inf_V"]=(candidate-x).tail(n).lpNorm<Eigen::Infinity>()*v0;
                writeDDSolutionState(folder/("step"+std::to_string(iteration)+".h5"),unpack(candidate),p.newton.inputScaling);
                x=candidate;bias+=dv;
                result["iterations"].push_back(row);
                std::cout<<arm.name<<' '<<iteration<<' '<<row["diagnostics"].dump()<<std::endl;
                if(row["diagnostics"]["original_gates_pass"].get<bool>()) {
                    if(!result.contains("first_original_gate_iteration"))result["first_original_gate_iteration"]=iteration;
                    // Continue the diagnostic beyond an original residual-only
                    // stop to see whether the independently reported RMS settles.
                    // This does not change the production acceptance rule.
                    if(error.valid&&error.rms<1.){result["status"]="original_gates_and_reported_RMS_passed";break;}
                }
            } catch(const std::exception& e) {
                row["error"]=e.what();result["iterations"].push_back(row);result["status"]="invalid_or_failed_candidate";break;
            }
        }
        report["arms"].push_back(result);save();
    }
    return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
