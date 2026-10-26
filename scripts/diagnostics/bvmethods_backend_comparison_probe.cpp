// Uninstalled diagnostic: autonomous Poisson and matched-bias continuation,
// or independent same-state native replays. No native operator enters a solve.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace {
using namespace vela;
using Json=nlohmann::json;
namespace fs=std::filesystem;
void jsonFile(const fs::path& p,const Json& j){std::ofstream(p)<<j.dump(2)<<'\n';}
void matrixFile(const fs::path& p,const SparseMatrixd& a){
    std::ofstream f(p);f<<std::setprecision(17)<<"row,col,value\n";
    for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator v(a,k);v;++v)
        f<<v.row()<<','<<v.col()<<','<<v.value()<<'\n';
}
}
int main(int argc,char**argv){try{
    if(argc!=3)throw std::invalid_argument("manifest output-directory");
    Json m,cfg;{std::ifstream f(argv[1]);f>>m;}{std::ifstream f(m.at("config").get<std::string>());f>>cfg;}
    const fs::path out=argv[2];fs::create_directories(out);
    cfg=canonicalizeDeck(cfg);cfg["solver"]["max_iter"]=100;
    auto problem=loadNewtonProblem(m.at("config").get<std::string>(),cfg);
    DDStateCliScope stateScope(problem,cfg,configDirectory(m.at("config").get<std::string>()));
    auto owner=makeNewtonSolver(problem);const int n=problem.mesh.numNodes(),d=3*n;
    const auto& units=problem.newton.inputScaling.unitSystem();
    auto sys=owner.makeArclengthSystem("drain");auto unpack=owner.makeArclengthStateUnpacker();
    const double line=units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6;
    auto port=owner.makeArclengthContactCurrentFunctional("drain",line);
    PerformanceProfiler profiler({true,(out/"performance.json").string()});
    ActivePerformanceProfilerScope profiling(&profiler);
    Json report={{"backend",LinearSolver::selectedBackend()},{"points",Json::array()}};
    auto fields=[&](const fs::path& folder,const VectorXd& x,double inner){
        fs::create_directories(folder);const auto s=unpack(x);
        writeDDSolutionState(folder/"state.h5",s,problem.newton.inputScaling);
        const auto terms=owner.evaluatePoissonTerms(s).rows;
        std::ofstream f(folder/"fields.csv");f<<std::setprecision(17)<<"node,psi_V,fn_V,fp_V,n_cm3,p_cm3,Nc_cm3,Nv_cm3,ni_cm3\n";
        for(const auto&t:terms){int i=t.nodeId;auto cm=[&](double a){return units.internalConcentrationToM3(a)/1e6;};
            f<<i<<','<<s.psi[i]<<','<<s.phin[i]<<','<<s.phip[i]<<','<<cm(s.n[i])<<','<<cm(s.p[i])<<','<<cm(t.electronDensityOfStates)<<','<<cm(t.holeDensityOfStates)<<','<<cm(t.intrinsicDensity)<<'\n';}
        jsonFile(folder/"state.json",{{"inner_V",inner},{"x",std::vector<double>(x.data(),x.data()+x.size())}});
    };
    auto exportOperator=[&](const fs::path& folder,const VectorXd& x,double inner,double outer,const SparseMatrixd& device){
        fields(folder,x,inner);const auto fv=sys.parameterDerivative(x,inner);const auto p=port.linearize(x,inner);
        std::vector<Eigen::Triplet<double>> entries;
        for(int k=0;k<device.outerSize();++k)for(SparseMatrixd::InnerIterator v(device,k);v;++v)entries.emplace_back(v.row(),v.col(),v.value());
        for(int i=0;i<d;++i){if(fv[i]!=0)entries.emplace_back(i,d,fv[i]);if(p.stateDerivative[i]!=0)entries.emplace_back(d,i,1e7*p.stateDerivative[i]);}
        entries.emplace_back(d,d,1+1e7*p.parameterDerivative);SparseMatrixd a(d+1,d+1);a.setFromTriplets(entries.begin(),entries.end());
        matrixFile(folder/"matrix.csv",a);const VectorXd f=sys.residual(x,inner);
        std::ofstream r(folder/"rhs.csv");r<<std::setprecision(17)<<"row,rhs\n";for(int i=0;i<d;++i)r<<i<<','<<-f[i]<<'\n';
        r<<d<<','<<-(inner+1e7*p.value-outer)<<'\n';
        const double scale=owner.evaluateResidual(unpack(x)).potentialScale;
        jsonFile(folder/"identity.json",{{"inner_V",inner},{"outer_V",outer},{"dimension",d+1},{"potential_scale_V",scale},{"current_A_per_um",p.value},{"backend",LinearSolver::selectedBackend()}});
    };
    CoupledLoadLineConfig c;c.maxIterations=m.value("max_iterations",40);c.equationTolerance=1e-8;c.loadLineTolerance_V=1e-6;
    c.linearSolver="direct_bordered";c.lineSearchMode=m.value("line_search_mode",std::string{"merit"});c.filterEnvelopeFactor=1e12;
    c.maxVoltageUpdate_V=.1;c.updateConvergence.mode="off";
    if(m.contains("frames")){
        for(const auto& frame:m.at("frames")){
            auto s=owner.buildPoissonBlockInitialization().coldInitial;
            s.packedState.resize(0);s.packedLow.resize(0);s.electronQfReference.resize(0);s.holeQfReference.resize(0);s.phinIncrement.resize(0);s.phipIncrement.resize(0);
            for(int i=0;i<n;++i)s.psi[i]=frame.at("psi_V")[i];
            for(const auto&v:frame.at("all_electron"))s.n[v[0].get<int>()]=units.m3ToInternalConcentration(v[1].get<double>()*1e6);
            for(const auto&v:frame.at("all_hole"))s.p[v[0].get<int>()]=units.m3ToInternalConcentration(v[1].get<double>()*1e6);
            for(const auto&t:owner.evaluatePoissonTerms(s).rows){s.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;s.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
            const VectorXd x=owner.packArclengthState(s);const double v=frame.at("inner_V"),u=frame.at("outer_V");
            const auto folder=out/frame.at("label").get<std::string>();exportOperator(folder,x,v,u,sys.jacobian(x,v));
            // Force one direction for audit only, never count this as a converged trajectory.
            auto cc=c;cc.maxIterations=1;cc.equationTolerance=1e-30;cc.loadLineTolerance_V=1e-30;cc.linearizationAuditEnabled=true;
            auto result=CoupledLoadLineNewton(sys,port.value,cc,port).solve({x,v},u,1e7);
            Json direction={{"converged",result.converged},{"failure",result.failureReason}};
            if(!result.trace.empty()&&result.trace.front().rawStateUpdate.size()==d){
                const auto&t=result.trace.front();const double scale=owner.evaluateResidual(unpack(x)).potentialScale;
                std::ofstream f(folder/"direction.csv");f<<std::setprecision(17)<<"node,dpsi_V,dfn_V,dfp_V\n";
                for(int i=0;i<n;++i)f<<i<<','<<t.rawStateUpdate[i]*scale<<','<<t.rawStateUpdate[n+i]*scale<<','<<t.rawStateUpdate[2*n+i]*scale<<'\n';
                direction["raw_inner_V"]=t.rawVoltageUpdate_V;direction["backward_error"]=t.componentwiseBackwardError;
            }
            jsonFile(folder/"direction.json",direction);std::cout<<frame.at("label")<<" exported"<<std::endl;
        }
    }else{
        const auto poisson=owner.solvePoissonOnly(owner.buildPoissonBlockInitialization().coldInitial);
        auto x=owner.packArclengthState(poisson.solution);fields(out/"poisson",x,0.);
        Json ph=Json::array();for(const auto&t:poisson.history)ph.push_back({{"iteration",t.iter},{"residual",t.residualNorm},{"alpha",t.dampingFactor},{"raw_step",t.rawStepNorm}});
        report["poisson"]={{"converged",poisson.converged},{"iterations",poisson.iters},{"residual",poisson.finalResidualNorm},{"trace",ph}};
        jsonFile(out/"summary.json",report);if(!poisson.converged)throw std::runtime_error("Poisson did not converge");
        ArclengthState current{x,0.},previous=current;bool hasPrevious=false;
        int point=0;
        for(double outer:m.at("outer_targets_V")){
            const auto folder=out/("point"+std::to_string(point));int jacCalls=0;
            auto traced=sys;const auto jac=sys.jacobian;
            traced.jacobian=[&](const VectorXd& z,double v){auto a=jac(z,v);exportOperator(folder/("linearization"+std::to_string(jacCalls++)),z,v,outer,a);return a;};
            CoupledLoadLineNewton corrector(traced,port.value,c,port);
            auto seed=hasPrevious&&m.value("extrapolate",true)?corrector.secantPredict(previous,current,outer,1e7):current;
            auto result=corrector.solve(seed,outer,1e7);Json trace=Json::array();
            for(const auto&t:result.trace){Json row={{"iteration",t.iteration},{"equation_inf",t.equationResidualNorm},{"load_V",t.loadLineResidual_V},{"alpha",t.acceptedAlpha},{"reported_alpha",t.reportedAlpha},{"full_step_merit",t.fullStepMerit},{"raw_inner_V",t.rawVoltageUpdate_V}};
                if(t.rawStateUpdate.size()==d){const double scale=owner.evaluateResidual(unpack(seed.x)).potentialScale;row["raw_psi_inf_V"]=t.rawStateUpdate.head(n).lpNorm<Eigen::Infinity>()*scale;row["raw_fn_inf_V"]=t.rawStateUpdate.segment(n,n).lpNorm<Eigen::Infinity>()*scale;row["raw_fp_inf_V"]=t.rawStateUpdate.tail(n).lpNorm<Eigen::Infinity>()*scale;}
                trace.push_back(row);}
            Json rec={{"outer_V",outer},{"inner_V",result.state.lambda},{"current_A_per_um",result.directedCurrent_A_per_um},{"converged",result.converged},{"iterations",result.iterations},{"failure",result.failureReason},{"equation_inf",result.equationResidualNorm},{"load_V",result.loadLineResidual_V},{"trace",trace},{"jacobian_calls",jacCalls}};
            report["points"].push_back(rec);jsonFile(out/"summary.json",report);
            exportOperator(folder/"final",result.state.x,result.state.lambda,outer,sys.jacobian(result.state.x,result.state.lambda));
            std::cout<<"point "<<point<<" "<<rec.dump()<<std::endl;
            if(!result.converged)break;
            previous=current;current=result.state;hasPrevious=true;++point;
        }
    }
    profiler.writeJson();jsonFile(out/"summary.json",report);return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
