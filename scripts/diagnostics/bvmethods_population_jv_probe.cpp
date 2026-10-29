// Uninstalled fixed-population one-sided checks. No solver policy changes.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace {
using namespace vela;
using Json=nlohmann::json;
namespace fs=std::filesystem;
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::invalid_argument("manifest.json output-directory");
    Json m,original;{std::ifstream f(argv[1]);f>>m;}
    const std::string config=m.at("config");{std::ifstream f(config);f>>original;}
    const fs::path out=argv[2];fs::create_directories(out);
    Json summary=Json::array();
    const auto rows=m.at("row_nodes").get<std::vector<int>>();
    const auto hs=m.value("relative_steps",std::vector<double>{.1,.01,.001,.0001});
    for(double h:hs)if(!(h>0.&&h<.5))throw std::invalid_argument("positive relative steps below 0.5 required");
    for(const auto& frame:m.at("frames"))for(const std::string arm:{"full","noav","nosrh","transport"}){
        auto cfg=original;
        if(arm=="noav"||arm=="transport")cfg["solver"]["impact_ionization"]["model"]="none";
        if(arm=="nosrh"||arm=="transport")cfg["solver"]["recombination"]=Json::array();
        const double inner=frame.at("inner_V");
        for(auto& c:cfg["contacts"])if(c["name"]=="drain")c["bias"]=inner;
        cfg=canonicalizeDeck(cfg);auto p=loadNewtonProblem(config,cfg);
        DDStateCliScope scope(p,cfg,configDirectory(config));auto solver=makeNewtonSolver(p);
        const int n=p.mesh.numNodes();const auto& units=p.newton.inputScaling.unitSystem();
        const double vt=constants::kb*p.newton.temperature_K/constants::q;
        auto state=readDDSolutionState(frame.at("state_file").get<std::string>(),n,p.newton.inputScaling);
        auto unpack=solver.makeArclengthStateUnpacker();const VectorXd x=solver.packArclengthState(state);state=unpack(x);
        auto sys=solver.makeArclengthSystem("drain");const VectorXd F=sys.residual(x,inner);const auto J=sys.jacobian(x,inner);
        const double v0=solver.evaluateResidual(state).potentialScale;
        std::vector<double> nc(n),nv(n);for(const auto& t:solver.evaluatePoissonTerms(state).rows){nc[t.nodeId]=t.electronDensityOfStates;nv[t.nodeId]=t.holeDensityOfStates;}
        const auto folder=out/frame.at("label").get<std::string>()/arm;fs::create_directories(folder);
        std::ofstream ledger(folder/"ledger.csv");ledger<<std::setprecision(17)<<"node,e_flux,e_srh,e_av,e_rhs,h_flux,h_srh,h_av,h_rhs\n";
        for(const auto& t:solver.evaluateCarrierTermDiagnostics(state,true).rows)
            ledger<<t.nodeId<<','<<t.electronFlux<<','<<t.electronRecombination<<','<<t.electronImpact<<','<<t.electronResidual<<','<<t.holeFlux<<','<<t.holeRecombination<<','<<t.holeImpact<<','<<t.holeResidual<<'\n';
        std::ofstream fields(folder/"fields.csv");fields<<std::setprecision(17)<<"node,psi_V,fn_V,fp_V,n_cm3,p_cm3\n";
        for(int i:rows){if(i<0||i>=n)throw std::invalid_argument("row node outside mesh");fields<<i<<','<<state.psi[i]<<','<<state.phin[i]<<','<<state.phip[i]<<','<<units.internalConcentrationToM3(state.n[i])*1e-6<<','<<units.internalConcentrationToM3(state.p[i])*1e-6<<'\n';}
        std::ofstream fd(folder/"one_sided.csv");fd<<std::setprecision(17)<<"kind,column_node,row_block,row_node,h,F0,Fplus,Fminus,analytic,forward,backward,central,noise_estimate,roundtrip\n";
        double maxRoundtrip=0.,anchorChange=0.;int evaluations=0;
        for(const auto& direction:m.at("directions")){
            const int node=direction.at("node");const std::string kind=direction.at("kind");
            if(node<0||node>=n||!(state.n[node]>0.)||!(state.p[node]>0.))throw std::invalid_argument("positive semiconductor direction node required");
            double roundtrip=0.;
            auto populationPoint=[&](double offset){
                DDSolution next=state;next.packedState.resize(0);next.packedLow.resize(0);next.electronQfReference.resize(0);next.holeQfReference.resize(0);next.phinIncrement.resize(0);next.phipIncrement.resize(0);
                if(kind=="electron")next.n[node]*=1.+offset;
                else if(kind=="hole")next.p[node]*=1.+offset;
                else if(kind=="potential")next.psi[node]+=offset*vt;
                else throw std::invalid_argument("unknown population direction kind");
                for(const auto& t:solver.evaluatePoissonTerms(next).rows)if(t.nodeId==node){next.phin[node]=t.electronRequiredQuasiFermiPotential;next.phip[node]=t.holeRequiredQuasiFermiPotential;}
                const VectorXd xx=solver.packArclengthState(next);const auto back=unpack(xx);
                for(int i=0;i<n;++i)if(state.n[i]>0.)roundtrip=std::max({roundtrip,std::abs(back.n[i]/next.n[i]-1.),std::abs(back.p[i]/next.p[i]-1.)});
                return xx;
            };
            // Invert the zero perturbation too: a serialized QF roundtrip must
            // not be mistaken for a finite-difference response as h decreases.
            const VectorXd baseX=populationPoint(0.);const auto baseState=unpack(baseX);
            const VectorXd baseF=sys.residual(baseX,inner);const auto baseJ=sys.jacobian(baseX,inner);
            anchorChange=std::max(anchorChange,(baseF-F).lpNorm<Eigen::Infinity>());
            VectorXd tangent=VectorXd::Zero(3*n);
            if(kind=="electron")tangent[n+node]=-state.n[node]/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(baseState.n[node],nc[node])*vt/v0;
            else if(kind=="hole")tangent[2*n+node]=state.p[node]/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(baseState.p[node],nv[node])*vt/v0;
            else if(kind=="potential"){tangent[node]=vt/v0;tangent[n+node]=vt/v0;tangent[2*n+node]=vt/v0;}
            else throw std::invalid_argument("unknown population direction kind");
            const VectorXd analytic=baseJ*tangent;
            for(double h:hs){
                auto perturb=[&](double sign){
                    const VectorXd xx=populationPoint(sign*h);
                    ++evaluations;return VectorXd(sys.residual(xx,inner));
                };
                const VectorXd plus=perturb(1.),minus=perturb(-1.);maxRoundtrip=std::max(maxRoundtrip,roundtrip);
                if(!plus.allFinite()||!minus.allFinite()||roundtrip>1e-10)throw std::runtime_error("invalid population perturbation roundtrip");
                for(int block:{1,2})for(int i:rows){const int r=block*n+i;
                    const double noise=64.*std::numeric_limits<double>::epsilon()*std::max({std::abs(baseF[r]),std::abs(plus[r]),std::abs(minus[r]),std::numeric_limits<double>::min()});
                    fd<<kind<<','<<node<<','<<block<<','<<i<<','<<h<<','<<baseF[r]<<','<<plus[r]<<','<<minus[r]<<','<<analytic[r]<<','<<(plus[r]-baseF[r])/h<<','<<(baseF[r]-minus[r])/h<<','<<(plus[r]-minus[r])/(2*h)<<','<<noise<<','<<roundtrip<<'\n';
                }
            }
        }
        Json record={{"label",frame.at("label")},{"arm",arm},{"inner_V",inner},{"equation_inf",F.lpNorm<Eigen::Infinity>()},{"inverse_anchor_rhs_change",anchorChange},{"max_population_roundtrip",maxRoundtrip},{"residual_evaluations",evaluations},{"potential_scale_V",v0}};
        summary.push_back(record);std::ofstream(out/"summary.json")<<summary.dump(2)<<'\n';std::cout<<record.dump()<<std::endl;
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
