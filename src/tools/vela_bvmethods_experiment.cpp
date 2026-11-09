// Independent experimental entry: production JSON schemas and defaults are unchanged.
#define main velaStandardEntryMain
#include "vela_example_runner.cpp"
#undef main
#include "vela/simulation/DCSweepExperiment.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/core/PerformanceProfiler.h"

namespace {
using Json=nlohmann::json;
Json readExperimentJson(const std::string& file) {
    std::ifstream input(file);Json value;
    if(!input)throw std::runtime_error("cannot open experiment request");
    input>>value;return value;
}
int checkpointExperiment(const std::string& file) {
#ifndef VELA_ENABLE_HDF5_STATE
    throw std::runtime_error("checkpoint experiments require HDF5 state support");
#else
    using namespace vela;
    namespace fs=std::filesystem;
    const auto req=readExperimentJson(file);
    const std::string configFile=req.at("config");
    const fs::path out=req.at("out").get<std::string>();fs::create_directories(out);
    auto cfg=canonicalizeDeck(readExperimentJson(configFile));
    auto problem=loadNewtonProblem(configFile,cfg);
    auto owner=makeNewtonSolver(problem);
    DDStateCliScope stateScope(problem,cfg,configDirectory(configFile));
    auto system=owner.makeArclengthSystem("drain",1e-4);
    const auto units=problem.newton.inputScaling.unitSystem();
    auto port=owner.makeArclengthContactCurrentFunctional("drain",
        units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6);
    ArclengthState state{experimental::readPackedExperiment(req.at("packed_binary")),req.at("inner")};
    if(state.x.size()!=3*problem.mesh.numNodes())throw std::runtime_error("checkpoint dimensions differ");
    experimental::writePackedExperiment(out/"initial.bin",state.x);
    experimental::writePackedExperiment(out/"initial_rhs.bin",system.residual(state.x,state.lambda));
    CoupledLoadLineConfig c;c.maxIterations=req.value("budget",40);
    c.applyVoltageUpdateLimit=false;c.applyDeviceUpdateLimit=false;
    c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.bootstrapModelGuardEnabled=req.value("guard",false);
    c.reportModelQuality=req.value("report_model_quality",true);
    c.equationTolerance=1e-8;c.loadLineTolerance_V=1e-6;
    c.updateConvergence.mode="enforce_joint";c.enforceLoadLineConvergence=true;
    auto metric=owner.makeArclengthUpdateErrorEvaluator("drain",c.updateConvergence);
    CoupledLoadLineNewton solver(system,port.value,c,port,metric);
    PerformanceProfilingConfig pc;pc.enabled=true;PerformanceProfiler profiler(pc);
    ActivePerformanceProfilerScope profileScope(&profiler);
    Json points=Json::array();bool passed=true;const auto started=std::chrono::steady_clock::now();
    for(double target:req.at("targets").get<std::vector<double>>()) {
        const auto begin=std::chrono::steady_clock::now();
        // Each target resets history, exactly as the G2/G3 sweep controls.
        auto result=solver.solve(state,target,1e7);
        Json trace=Json::array();
        for(const auto& t:result.trace)trace.push_back({
            {"iteration",t.iteration},{"mode",t.derivativeMode},{"reason",t.derivativeReason},
            {"raw_dv",t.rawVoltageUpdate_V},{"applied_alpha",t.acceptedAlpha},
            {"protected_n",t.protectedElectrons},{"protected_p",t.protectedHoles},
            {"linear_absolute_defect",t.linearGateAbsolute},
            {"linear_absolute_tolerance",t.linearGateAbsoluteTolerance},
            {"linear_componentwise_error",t.linearGateBackward},
            {"linear_backward_tolerance",t.linearGateBackwardTolerance},
            {"phi_before",t.modelReportBefore},{"phi_model",t.modelReportPredicted},
            {"phi_actual",t.modelReportActual},{"rho_reliable",t.modelReportReliable},
            {"rho",t.modelReportRatio},{"model_report_valid",t.modelReportValid},
            {"rejected",t.candidateRejectedByRecovery},{"model_review",t.bootstrapModelRetry},
            {"linear_review",!t.bootstrapLinearFailureReason.empty()},
            {"rms",t.updateError.rms},{"candidate_rms",t.candidateUpdateError.rms}});
        long double kcl=0.;
        for(const auto& contact:problem.mesh.contacts()) {
            const auto current=owner.makeArclengthContactCurrentFunctional(contact.name,
                units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6);
            kcl+=current.value(result.state.x,result.state.lambda);
        }
        const double current=port.value(result.state.x,result.state.lambda);
        const bool kclPass=std::isfinite(double(kcl))&&std::abs(kcl)<=std::max(1e-12,1e-5*std::abs(current));
        const bool ok=result.converged&&result.equationResidualNorm<=1e-8&&
            std::abs(result.loadLineResidual_V)<=1e-6&&result.updateError.valid&&
            result.updateError.rms<1&&kclPass&&
            (!system.convergenceSatisfied||system.convergenceSatisfied(result.state.x,result.state.lambda));
        points.push_back({{"target",target},{"strict_gates",ok},{"converged",result.converged},
            {"iterations",result.iterations},{"rhs",result.equationResidualNorm},
            {"load",result.loadLineResidual_V},{"rms",result.updateError.rms},
            {"inner",result.state.lambda},{"current",current},{"kcl",double(kcl)},
            {"failure",result.failureReason},{"trace",trace},
            {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()}});
        experimental::writePackedExperiment(out/("point_"+std::to_string(points.size()-1)+"_final.bin"),result.state.x);
        if(!ok){passed=false;break;}state=result.state;
    }
    Json report={{"passed",passed},{"points",points},{"profile",profiler.toJson()},
        {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()}};
    std::ofstream(out/"report.json")<<report.dump(2);
    std::cout<<Json{{"experiment_completed",true},{"passed",passed},{"points",points.size()}}.dump()<<'\n';
    return 0; // Numeric qualification is explicit in report.json, not inferred from process success.
#endif
}
} // namespace

int main(int argc,char** argv) {
    try {
        if(argc==3&&std::string(argv[1])=="--checkpoint")return checkpointExperiment(argv[2]);
        vela::experimental::DCSweepExperiment options;
        std::vector<char*> forwarded;forwarded.push_back(argv[0]);bool enabled=false;
        for(int i=1;i<argc;++i) {
            if(std::string(argv[i])=="--experiment") {
                if(++i>=argc||enabled)throw std::invalid_argument("one --experiment request is required");
                const auto req=readExperimentJson(argv[i]);enabled=true;
                options.bootstrapModelGuardEnabled=req.value("guard",false);
                options.reportModelQuality=req.value("report_model_quality",false);
                options.directory=req.value("out",std::string{});
                options.initialPackedStateFile=req.value("packed_binary",std::string{});
                options.initialInnerVoltage_V=req.value("inner",0.);
            } else forwarded.push_back(argv[i]);
        }
        if(!enabled)throw std::invalid_argument("use --experiment options.json --config deck.json, or --checkpoint request.json");
        vela::experimental::DCSweepExperimentScope scope(options);
        return velaStandardEntryMain(static_cast<int>(forwarded.size()),forwarded.data());
    } catch(const std::exception& error) {
        std::cerr<<"BVMethods experiment failed: "<<error.what()<<'\n';return 1;
    }
}
