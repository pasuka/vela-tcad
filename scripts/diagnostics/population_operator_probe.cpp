// Validation client of the core mode, never a second nonlinear implementation.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace {
using namespace vela;
using Json=nlohmann::json;
namespace fs=std::filesystem;
void matrix(const fs::path& path,const SparseMatrixd& a) {
    std::ofstream f(path);f<<std::setprecision(17)<<"row,col,value\n";
    for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator e(a,k);e;++e)
        f<<e.row()<<','<<e.col()<<','<<e.value()<<'\n';
}
Json support(const SparseMatrixd& a) {
    int nonzero=0;for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator e(a,k);e;++e)nonzero+=e.value()!=0.;
    return {{"dimension",a.rows()},{"stored",a.nonZeros()},{"strict_nonzero",nonzero}};
}
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::invalid_argument("manifest.json output-directory");
    Json manifest,cfg;{std::ifstream f(argv[1]);f>>manifest;}
    const std::string config=manifest.at("config");{std::ifstream f(config);f>>cfg;}
    const fs::path out=argv[2];fs::create_directories(out);Json reports=Json::array();
    for(const auto& frame:manifest.at("frames")) {
        auto control=cfg;control["solver"]["state_variables"]="potential_quasi_fermi";
        control["solver"].erase("population_update");
        const double bias=frame.at("inner_V");
        for(auto& contact:control["contacts"])if(contact["name"]=="drain")contact["bias"]=bias;
        control=canonicalizeDeck(control);auto qp=loadNewtonProblem(config,control);auto q=makeNewtonSolver(qp);
        auto direct=control;direct["solver"]["state_variables"]="potential_density";
        direct["solver"]["quasi_fermi_update_limit_V"]=0.;direct["solver"]["quasi_fermi_update_limit_minority_V"]=0.;
        auto pp=loadNewtonProblem(config,direct);auto p=makeNewtonSolver(pp);DDStateCliScope stateScope(pp,direct,configDirectory(config));
        auto physical=readDDSolutionState(frame.at("state_file").get<std::string>(),pp.mesh.numNodes(),pp.newton.inputScaling);
        // Define the common comparison state by the QF operator's actual
        // forward population, then store that population exactly in the new mode.
        physical=p.unpackArclengthState(p.packArclengthState(physical));
        const auto x=q.packArclengthState(physical);physical=q.unpackArclengthState(x);
        const auto y=p.packArclengthState(physical);const auto ps=p.makeArclengthSystem("drain"),qs=q.makeArclengthSystem("drain");
        const auto F=ps.residual(y,bias),Q=qs.residual(x,bias);const auto J=ps.jacobian(y,bias),K=qs.jacobian(x,bias);
        const int n=pp.mesh.numNodes(),d=3*n;const double vt=constants::kb*pp.newton.temperature_K/constants::q;
        const double v0=p.evaluateResidual(physical).potentialScale;std::vector<double> nc(n),nv(n);std::vector<bool> active(n);
        for(const auto& t:p.evaluatePoissonTerms(physical).rows){nc[t.nodeId]=t.electronDensityOfStates;nv[t.nodeId]=t.holeDensityOfStates;active[t.nodeId]=t.intrinsicDensity>0.;}
        std::vector<Eigen::Triplet<double>> entries;VectorXd scale=VectorXd::Ones(d);
        for(int i=0;i<n;++i) {
            entries.emplace_back(i,i,1.);
            if(active[i]) {
                entries.emplace_back(n+i,i,1.);entries.emplace_back(2*n+i,i,1.);
                entries.emplace_back(n+i,n+i,-vt/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(physical.n[i],nc[i])/v0);
                entries.emplace_back(2*n+i,2*n+i,vt/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(physical.p[i],nv[i])/v0);
                scale[n+i]=physical.n[i];scale[2*n+i]=physical.p[i];
            }else {entries.emplace_back(n+i,n+i,1.);entries.emplace_back(2*n+i,2*n+i,1.);}
        }
        SparseMatrixd T(d,d);T.setFromTriplets(entries.begin(),entries.end());const SparseMatrixd transformed=K*T;
        const SparseMatrixd scaled=J*scale.asDiagonal(),expected=transformed*scale.asDiagonal(),diff=scaled-expected;
        const auto folder=out/frame.at("label").get<std::string>();fs::create_directories(folder);
        matrix(folder/"population.csv",J);matrix(folder/"qf_chain.csv",transformed);
        {std::ofstream f(folder/"rhs.csv");f<<std::setprecision(17)<<"row,population,qf\n";for(int i=0;i<d;++i)f<<i<<','<<F[i]<<','<<Q[i]<<'\n';}
        Json report={{"label",frame.at("label")},{"population",support(J)},{"chain",support(transformed)},
            {"reduced_dimension",d+1-static_cast<int>(ps.prescribedRows(y,bias).size())},
            {"rhs_absolute_inf",(F-Q).lpNorm<Eigen::Infinity>()},
            {"scaled_matrix_relative_frobenius",diff.norm()/std::max(expected.norm(),1e-300)},
            {"fd",Json::array()}};
        for(const auto& direction:manifest.value("directions",Json::array({{{"kind","global"}}}))) {
            VectorXd v=VectorXd::Zero(d);const std::string kind=direction.at("kind");
            if(kind=="global")for(int i=0;i<n;++i) {
                v[i]=.01*std::sin(i*.71);
                if(active[i]){v[n+i]=physical.n[i]*.03*std::sin(i*.37);v[2*n+i]=physical.p[i]*.03*std::cos(i*.43);}
            }else {
                const int node=direction.at("node");if(node<0||node>=n||!active[node])throw std::invalid_argument("invalid FD node");
                if(kind!="electron"&&kind!="hole")throw std::invalid_argument("invalid FD kind");
                const int k=(kind=="electron"?n:2*n)+node;v[k]=y[k];
            }
            const VectorXd analytic=J*v;
            for(double h:manifest.value("relative_steps",std::vector<double>{.01,.001,.0001})) {
                if(!(h>0.&&h<.5))throw std::invalid_argument("invalid FD step");
                const VectorXd plus=ps.residual(y+h*v,bias),minus=ps.residual(y-h*v,bias);
                const VectorXd forward=(plus-F)/h,backward=(F-minus)/h,central=(plus-minus)/(2*h);
                if(kind!="global" && h==.001) {
                    std::ofstream f(folder/(kind+"_fd.csv"));f<<std::setprecision(17)<<"row,analytic,forward,backward,base\n";
                    for(int i=0;i<d;++i)if(analytic[i]!=0.||forward[i]!=0.||backward[i]!=0.)f<<i<<','<<analytic[i]<<','<<forward[i]<<','<<backward[i]<<','<<F[i]<<'\n';
                }
                auto relative=[&](const VectorXd& value){return (value-analytic).norm()/std::max({value.norm(),analytic.norm(),1e-300});};
                Json blocks=Json::array();
                for(int b=0;b<3;++b){const auto a=analytic.segment(b*n,n),f=central.segment(b*n,n);
                    blocks.push_back({{"analytic_l2",a.norm()},{"absolute_error_l2",(f-a).norm()},
                        {"relative",(f-a).norm()/std::max({f.norm(),a.norm(),1e-300})}});}
                report["fd"].push_back({{"direction",direction},{"h",h},{"analytic_l2",analytic.norm()},
                    {"forward_relative",relative(forward)},{"backward_relative",relative(backward)},{"central_relative",relative(central)},{"blocks",blocks}});
            }
        }
        reports.push_back(report);std::ofstream f(out/"summary.json");f<<reports.dump(2)<<'\n';
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
