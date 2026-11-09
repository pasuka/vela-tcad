// Read-only frozen-state audit. No nonlinear updates or acceptance changes.
#define main velaEmbeddedStandardMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/DCSweepExperiment.h"
#include <cmath>

namespace {
using Json=nlohmann::json;
Json readJson(const std::filesystem::path& file) {
    std::ifstream stream(file);Json value;
    if(!stream)throw std::runtime_error("Cannot open frozen audit input");
    stream>>value;return value;
}
void writeJson(const std::filesystem::path& file,const Json& value) {
    std::ofstream stream(file);stream<<value.dump(2)<<'\n';
    if(!stream)throw std::runtime_error("Cannot write frozen audit output");
}
Json norm(const vela::VectorXd& v) {
    return {{"finite",v.allFinite()},{"inf",v.cwiseAbs().maxCoeff()},
            {"l2",v.stableNorm()}};
}
vela::VectorXd multiplyWide(const vela::SparseMatrixd& j,const vela::VectorXd& z) {
    std::vector<long double> sum(j.rows());
    for(int c=0;c<j.outerSize();++c)
        for(vela::SparseMatrixd::InnerIterator it(j,c);it;++it)
            sum[it.row()]+=(long double)it.value()*z[c];
    vela::VectorXd result(j.rows());
    for(int i=0;i<result.size();++i)result[i]=double(sum[i]);
    return result;
}
int nonfinite(const vela::SparseMatrixd& j) {
    int count=0;
    for(int c=0;c<j.outerSize();++c)
        for(vela::SparseMatrixd::InnerIterator it(j,c);it;++it)
            count+=!std::isfinite(it.value());
    return count;
}
}
int main(int argc,char** argv) {
    using namespace vela;namespace fs=std::filesystem;
    try {
        if(argc!=2)throw std::invalid_argument("Usage: vela_physics_frozen_audit request.json");
        const auto request=readJson(argv[1]);const fs::path out=request.at("out");
        fs::create_directories(out);const std::string config=request.at("config");
        const VectorXd x=experimental::readPackedExperiment(request.at("packed_binary"));
        const double inner=request.at("inner");
        Json report={{"diagnostic_only",true},{"direction","fixed relative population; forward positive ray"}};
        std::vector<VectorXd> residuals,directions,fluxes,avalanches,srhs;
        for(const std::string arithmetic:{"quasi_fermi","band_population"}) {
            auto cfg=canonicalizeDeck(readJson(config));
            cfg["solver"]["population_sg_current"]=arithmetic;
            // Only Full matrices are compared. This disables no live tangent;
            // the bootstrap-only control is incompatible with band arithmetic.
            if(arithmetic=="band_population")
                cfg["solver"]["population_sg_bootstrap_linearization"]="unchanged";
            auto problem=loadNewtonProblem(config,cfg);
            auto owner=makeNewtonSolver(problem);
            DDStateCliScope scope(problem,cfg,configDirectory(config));
            auto system=owner.makeArclengthSystem("drain",1e-4);
            const int n=problem.mesh.numNodes();
            if(x.size()!=3*n)throw std::invalid_argument("Frozen vector dimension mismatch");
            VectorXd columns=VectorXd::Ones(x.size()),z(x.size());
            for(int i=0;i<z.size();++i) {
                z[i]=i<n?1e-3*std::sin(i+1.):1e-4*(.1+std::abs(std::sin(i+1.)));
                if(i>=n) {
                    // Inactive (oxide) carrier slots are exactly zero. Keep a
                    // positive coordinate scale, but do not perturb those slots.
                    columns[i]=x[i]>0.?x[i]:1.;
                    if(x[i]<=0.)z[i]=0.;
                }
            }
            const VectorXd physicalDirection=columns.cwiseProduct(z);
            auto f=system.residual(x,inner);
            auto j=system.jacobianInColumns(x,inner,columns);
            const auto jv=multiplyWide(j,z);
            const auto unpack=owner.makeArclengthStateUnpacker();
            const auto terms=owner.evaluateCarrierTermDiagnostics(unpack(x),true).rows;
            VectorXd flux=VectorXd::Zero(2*n),av=flux,srh=flux;
            for(const auto& t:terms) {
                flux[t.nodeId]=t.electronFlux;flux[n+t.nodeId]=t.holeFlux;
                av[t.nodeId]=t.electronImpact;av[n+t.nodeId]=t.holeImpact;
                srh[t.nodeId]=t.electronRecombination;srh[n+t.nodeId]=t.holeRecombination;
            }
            const auto directory=out/arithmetic;fs::create_directories(directory);
            experimental::writePackedExperiment(directory/"residual.bin",f);
            experimental::writePackedExperiment(directory/"jv.bin",jv);
            experimental::writePackedExperiment(directory/"columns.bin",columns);
            experimental::writePackedExperiment(directory/"direction.bin",physicalDirection);
            experimental::writePackedExperiment(directory/"transport.bin",flux);
            experimental::writePackedExperiment(directory/"av.bin",av);
            experimental::writePackedExperiment(directory/"srh.bin",srh);
            Json measured={{"F",norm(f)},{"J_nonfinite",nonfinite(j)},
                {"Jv",norm(jv)},{"transport",norm(flux)},
                {"av",norm(av)},{"srh",norm(srh)},{"forward_fd",Json::array()}};
            if(arithmetic=="quasi_fermi") {
                const auto boot=system.bootstrapJacobianInColumns(x,inner,columns);
                measured["bootstrap_nonfinite"]=nonfinite(boot);
            }
            const auto constrained=system.prescribedRows(x,inner);
            auto free=[&](VectorXd value) {for(int row:constrained)value[row]=0.;return value;};
            measured["prescribed_rows"]=constrained;
            measured["free_F"]=norm(free(f));measured["free_Jv"]=norm(free(jv));
            measured["blocks"]=Json::array();
            const auto freeF=free(f),freeJv=free(jv);
            for(int block=0;block<3;++block)
                measured["blocks"].push_back({{"block",block},
                    {"F",norm(freeF.segment(block*n,n))},
                    {"Jv",norm(freeJv.segment(block*n,n))}});
            for(double h:{1.,.1,.01}) {
                const VectorXd trial=x+h*physicalDirection;
                const VectorXd fd=(system.residual(trial,inner)-f)/h;
                experimental::writePackedExperiment(directory/("fd_"+std::to_string(h)+".bin"),fd);
                measured["forward_fd"].push_back({{"h",h},
                    {"free_error",norm(free(fd-jv))},{"free_Jv",norm(free(jv))}});
            }
            measured["dimension"]=x.size();measured["nodes"]=n;
            report[arithmetic]=measured;
            residuals.push_back(f);directions.push_back(jv);
            fluxes.push_back(flux);avalanches.push_back(av);srhs.push_back(srh);
        }
        report["difference"]={{"F",norm(residuals[1]-residuals[0])},
            {"Jv",norm(directions[1]-directions[0])},
            {"transport",norm(fluxes[1]-fluxes[0])},
            {"av",norm(avalanches[1]-avalanches[0])},
            {"srh",norm(srhs[1]-srhs[0])}};
        writeJson(out/"audit.json",report);
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"Frozen physics audit failed: "<<e.what()<<'\n';return 1;
    }
}
