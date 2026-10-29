// Uninstalled regression driver. Build with the same Release compiler and
// link flags as vela_example_runner. Inputs contain physical QF state only;
// native SG caches, source values and Jacobians are not accepted.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"
#include <filesystem>

int main(int argc,char**argv) {try {
    using namespace vela;
    if(argc!=3)throw std::invalid_argument("usage: bvmethods_profile_probe manifest.json output-directory");
    nlohmann::json manifest;{std::ifstream f(argv[1]);f>>manifest;}
    nlohmann::json config;{std::ifstream f(manifest.at("config").get<std::string>());f>>config;}
    config=canonicalizeDeck(config);
    config["solver"]["coupled_operator"]="bvmethods_tri3_qf";
    config["solver"]["carrier_statistics"]={{"model","fermi_dirac_joyce_dixon"}};
    config["solver"]["impact_ionization"]["source_jacobian"]="local_ad";
    const std::filesystem::path root=argv[2];std::filesystem::create_directories(root);
    const auto matrixFile=[](const auto& path,const SparseMatrixd&M){std::ofstream f(path);f<<std::setprecision(17)<<"row,col,value\n";
        for(int k=0;k<M.outerSize();++k)for(SparseMatrixd::InnerIterator v(M,k);v;++v)f<<v.row()<<','<<v.col()<<','<<v.value()<<'\n';};
    nlohmann::json report={{"schema","vela.bvmethods_profile_regression.v1"},{"operator","bvmethods_tri3_qf"},
        {"coordinates","scaled physical psi, fn, fp; populations reconstructed autonomously"},
        {"SG_arithmetic","log-density/expm1 QF; no native values"},{"frames",nlohmann::json::array()},{"sweep",nlohmann::json::array()},{"attempts",nlohmann::json::array()},{"continuation",nlohmann::json::array()}};
    for(const auto&frame:manifest.at("frames")) {
        const std::string label=frame.at("label");const double inner=frame.at("inner_V"),outer=frame.at("outer_V");
        auto cfg=config;for(auto&c:cfg["contacts"])if(c["name"]=="drain")c["bias"]=inner;
        auto p=loadNewtonProblem(manifest.at("config").get<std::string>(),cfg);auto solver=makeNewtonSolver(p);
        const int n=p.mesh.numNodes(),d=3*n;DDSolution input;
        input.psi.resize(n);input.phin.resize(n);input.phip.resize(n);input.n=VectorXd::Zero(n);input.p=VectorXd::Zero(n);
        for(int i=0;i<n;++i){input.psi[i]=frame.at("psi_V")[i];input.phin[i]=frame.at("fn_V")[i];input.phip[i]=frame.at("fp_V")[i];}
        auto x=solver.packArclengthState(input);auto sys=solver.makeArclengthSystem("drain");auto unpack=solver.makeArclengthStateUnpacker();
        const auto&u=p.newton.inputScaling.unitSystem();const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*1e-6;
        auto port=solver.makeArclengthContactCurrentFunctional("drain",line);
        auto fields=unpack(x);auto ev=solver.evaluateResidual(fields);const double vscale=ev.potentialScale;
        const auto device=sys.jacobian(x,inner);const auto fv=sys.parameterDerivative(x,inner);const auto terminal=port.linearize(x,inner);
        std::vector<Eigen::Triplet<double>> entries;
        for(int k=0;k<device.outerSize();++k)for(SparseMatrixd::InnerIterator v(device,k);v;++v)entries.emplace_back(v.row(),v.col(),v.value());
        for(int i=0;i<d;++i){if(fv[i]!=0)entries.emplace_back(i,d,fv[i]);if(terminal.stateDerivative[i]!=0)entries.emplace_back(d,i,1e7*terminal.stateDerivative[i]);}
        entries.emplace_back(d,d,1+1e7*terminal.parameterDerivative);SparseMatrixd M(d+1,d+1);M.setFromTriplets(entries.begin(),entries.end());
        VectorXd F(d+1);F.head(d)=sys.residual(x,inner);F[d]=inner+1e7*terminal.value-outer;
        auto folder=root/label;std::filesystem::create_directories(folder);matrixFile(folder/"matrix.csv",M);
        {std::ofstream f(folder/"rhs.csv");f<<std::setprecision(17)<<"row,rhs,x\n";for(int i=0;i<=d;++i)f<<i<<','<<-F[i]<<','<<(i==d?inner:x[i])<<'\n';}
        const auto rows=solver.evaluateCarrierTermDiagnostics(fields,true).rows;double ledgerError=0.;
        {std::ofstream f(folder/"fields.csv");f<<std::setprecision(17)<<"node,psi,fn,fp,n_cm3,p_cm3,e_flux,h_flux,e_srh,h_srh,e_AV,h_AV\n";
         for(const auto&t:rows){int i=t.nodeId;f<<i<<','<<fields.psi[i]<<','<<fields.phin[i]<<','<<fields.phip[i]<<','<<u.internalConcentrationToM3(fields.n[i])/1e6<<','<<u.internalConcentrationToM3(fields.p[i])/1e6<<','<<t.electronFlux<<','<<t.holeFlux<<','<<t.electronRecombination<<','<<t.holeRecombination<<','<<t.electronImpact<<','<<t.holeImpact<<'\n';
          ledgerError=std::max({ledgerError,std::abs(t.electronResidual-F[n+i]),std::abs(t.holeResidual-F[2*n+i])});}}
        auto meta=nlohmann::json{{"label",label},{"dimension",d+1},{"potential_scale_V",vscale},{"thermal_voltage_V",constants::kb*p.newton.temperature_K/constants::q},
            {"current_A_per_um",terminal.value},{"ledger_residual_difference_Linf",ledgerError},{"Jv",nlohmann::json::array()}};
        for(int dir=0;dir<3;++dir){VectorXd v=VectorXd::Zero(d);double dv=0.;
            if(dir<2)v[(dir==0?n:2*n)+2587]=1e-3/vscale;
            else{for(int i=0;i<d;++i)v[i]=3e-4*std::sin(i*.37)/vscale;dv=1e-4;}
            const VectorXd an=device*v+fv*dv;
            for(double eps:{1e-3,1e-5,1e-7,1e-9}){const VectorXd fd=(sys.residual(x+eps*v,inner+eps*dv)-sys.residual(x-eps*v,inner-eps*dv))/(2*eps);
                VectorXd a=an,b=fd;for(int i=0;i<d;++i){a[i]/=device.coeff(i,i);b[i]/=device.coeff(i,i);}
                meta["Jv"].push_back({{"direction",dir},{"epsilon",eps},{"relative_L2",(a-b).norm()/std::max(a.norm(),b.norm())}});
            }
        }
        std::ofstream(folder/"identity.json")<<meta.dump(2);report["frames"].push_back(meta);
        if(label=="base"&&manifest.contains("outer_targets_V")){
            CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.maxIterations=40;c.equationTolerance=1e-8;c.loadLineTolerance_V=1e-6;
            c.applyDeviceUpdateLimit=false;c.updateConvergence.mode="enforce";
            CoupledLoadLineNewton coupled(sys,port.value,c,port,solver.makeArclengthUpdateErrorEvaluator("drain",c.updateConvergence));
            ArclengthState current;current.x=x;current.lambda=inner;
            double acceptedOuter=outer;
            const auto& controls=cfg.at("sweep").at("external_circuit");
            double nextStep=controls.value("coupled_initial_outer_step_V",.01);
            const double minimumStep=controls.value("coupled_min_outer_step_V",1e-5);
            const double maximumStep=controls.value("coupled_max_outer_step_V",2500.);
            const double growth=controls.value("coupled_outer_growth_factor",1.41);
            const double shrink=controls.value("coupled_outer_shrink_factor",.5);
            const int maxRetries=controls.value("coupled_max_step_retries",20);
            bool failed=false;
            for(double goal:manifest.at("outer_targets_V")){
                double target=std::min(goal,acceptedOuter+nextStep);int retries=0;
                while(true){
                    if(report["attempts"].size()>=128)throw std::runtime_error("validation continuation exceeded 128 attempts");
                    const auto result=coupled.solve(current,target,1e7);nlohmann::json trace=nlohmann::json::array();
                    for(const auto&t:result.trace)trace.push_back({{"iteration",t.iteration},{"equation_inf",t.equationResidualNorm},{"load_V",t.loadLineResidual_V},{"alpha",t.acceptedAlpha},{"raw_update_norm",t.rawStateUpdateNorm},{"update_rms",t.updateError.rms}});
                    nlohmann::json record={{"outer_V",target},{"requested_goal_V",goal},{"from_outer_V",acceptedOuter},{"inner_V",result.state.lambda},{"current_A_per_um",result.directedCurrent_A_per_um},
                        {"converged",result.converged},{"iterations",result.iterations},{"equation_inf",result.equationResidualNorm},{"load_V",result.loadLineResidual_V},
                        {"reason",result.convergenceReason},{"failure",result.failureReason},{"trace",trace}};
                    report["attempts"].push_back(record);
                    if(result.converged){
                        const double step=target-acceptedOuter;
                        current=result.state;acceptedOuter=target;report["continuation"].push_back(record);
                        if(step>0.)nextStep=std::min(maximumStep,step*growth);
                        if(target==goal){report["sweep"].push_back(record);break;}
                        target=std::min(goal,acceptedOuter+nextStep);retries=0;
                    }else{
                        const double step=(target-acceptedOuter)*shrink;
                        if(step<minimumStep || ++retries>maxRetries){failed=true;report["sweep"].push_back(record);break;}
                        // Retry from the last accepted state; a failed iterate
                        // never becomes a continuation seed. Same 40-iteration
                        // budget and tolerances at every smaller outer step.
                        target=acceptedOuter+step;
                    }
                    std::ofstream(root/"summary.json")<<report.dump(2);
                }
                std::ofstream(root/"summary.json")<<report.dump(2);if(failed)break;
            }
        }
        std::cout<<label<<" exported"<<std::endl;
    }
    std::ofstream(root/"summary.json")<<report.dump(2);
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
