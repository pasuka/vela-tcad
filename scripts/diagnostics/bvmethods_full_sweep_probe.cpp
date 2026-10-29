// Case-specific BVMethods continuation, based on the 2026-10-25 full sweep.
// This uninstalled diagnostic accepts locally generated physical-potential
// checkpoints. It does not add a CMD/CLI feature or change production defaults.
// See docs/bvmethods_operator_contract.md for manifest and adaptive controls.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"
#include "BvMethodsAdaptiveStep.h"
#include <filesystem>
int main(int argc,char**argv){try{
 using namespace vela;namespace fs=std::filesystem;
 if(argc!=3)throw std::invalid_argument("manifest outputdir");
 nlohmann::json m;{std::ifstream f(argv[1]);f>>m;}
 nlohmann::json cfg;{std::ifstream f(m.at("config").get<std::string>());f>>cfg;}
 const auto innerControl=diagnostics::resolveInnerStepLimit(cfg,m);
 const auto outerControl=diagnostics::resolveOuterStepControl(cfg,m);
 const auto circuit=diagnostics::stepObject(diagnostics::stepObject(cfg,"sweep"),"external_circuit");
 const double correctorVoltageLimit=circuit.contains("max_inner_voltage_step_V")
  ? diagnostics::positiveStepNumber(circuit,"max_inner_voltage_step_V") : .1;
 cfg=canonicalizeDeck(cfg);const fs::path out=argv[2];fs::create_directories(out);
 const bool extrapolate=m.value("extrapolate",diagnostics::stepObject(cfg,"sweep").value("predictor",std::string{"linear"})=="linear");
 const double start=m["current"]["inner_V"];
 for(auto&c:cfg["contacts"])if(c["name"]=="drain")c["bias"]=start;
 auto p=loadNewtonProblem(m["config"],cfg);auto solver=makeNewtonSolver(p);const int n=p.mesh.numNodes();
 auto pack=[&](const nlohmann::json&f){DDSolution s;s.psi.resize(n);s.phin.resize(n);s.phip.resize(n);s.n=s.p=VectorXd::Zero(n);
  for(int i=0;i<n;++i){s.psi[i]=f["psi_V"][i];s.phin[i]=f["fn_V"][i];s.phip[i]=f["fp_V"][i];}
  return ArclengthState{solver.packArclengthState(s),f["inner_V"].get<double>()};};
 auto current=pack(m["current"]),previous=pack(m["previous"]);double outer=m["current"]["outer_V"];double previousOuter=m["previous"]["outer_V"];
 auto sys=solver.makeArclengthSystem("drain");auto unpack=solver.makeArclengthStateUnpacker();auto evaluate=solver.makeArclengthConvergenceEvaluator("drain");
 const auto&u=p.newton.inputScaling.unitSystem();const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*1e-6;
 auto drain=solver.makeArclengthContactCurrentFunctional("drain",line);
 std::vector<ArclengthScalarFunctional> ports;for(const auto&c:p.mesh.contacts())ports.push_back(solver.makeArclengthContactCurrentFunctional(c.name,line));
 CoupledLoadLineConfig c;c.maxIterations=40;c.equationTolerance=1e-8;c.loadLineTolerance_V=1e-6;c.maxVoltageUpdate_V=correctorVoltageLimit;
 c.linearSolver=circuit.value("coupled_linear_solver",std::string{"direct_bordered"});c.lineSearchMode=circuit.value("coupled_line_search_mode",std::string{"merit"});c.applyDeviceUpdateLimit=true;c.filterEnvelopeFactor=circuit.value("coupled_filter_envelope_factor",1e12);c.updateConvergence.mode="off";
 CoupledLoadLineNewton corrector(sys,drain.value,c,drain,solver.makeArclengthUpdateErrorEvaluator("drain",c.updateConvergence));
 nlohmann::json report={{"schema","vela.full_bv_fallback.v1"},{"source","own accepted baseline states"},{"start_outer_V",outer},
  {"attempts",nlohmann::json::array()},{"points",nlohmann::json::array()},{"complete",false},
  {"inner_step_guard_V",innerControl.volts},{"inner_step_guard_source",innerControl.source},
  {"inner_step_guard_reason",innerControl.reason},{"corrector_voltage_update_limit_V",correctorVoltageLimit},
  {"line_search_mode",c.lineSearchMode},{"filter_envelope_factor",c.filterEnvelopeFactor},{"extrapolate",extrapolate},
  {"outer_step_control",{{"initial_step_V",outerControl.initial},{"min_step_V",outerControl.minimum},
   {"max_step_V",outerControl.maximum},{"growth_factor",outerControl.growth},
   {"shrink_factor",outerControl.shrink},{"growth_mode",outerControl.growthMode},{"max_retries",outerControl.maxRetries}}}};
 auto save=[&](){std::ofstream(out/"summary.json")<<report.dump(2);};
 auto attempt=[&](const ArclengthState&seed,double target,const std::string&kind){
  auto result=corrector.solve(seed,target,1e7);double kcl=0.,id=0.,load=0.;bool physical=false;
  nlohmann::json diag=nlohmann::json::object();
  if(result.converged){const auto e=evaluate(result.state.x,result.state.lambda);id=drain.value(result.state.x,result.state.lambda);
   bool finite=true;for(auto&port:ports){double v=port.value(result.state.x,result.state.lambda);finite=finite&&std::isfinite(v);kcl+=v;}
   load=result.state.lambda+1e7*id-target;
   physical=finite&&std::isfinite(e.blockNorms.psi)&&std::isfinite(e.blockNorms.phin)&&std::isfinite(e.blockNorms.phip)&&
    std::isfinite(e.carrierRows.maxRatio)&&std::isfinite(e.globalContinuity.electron.ratio)&&std::isfinite(e.globalContinuity.hole.ratio)&&
    (!e.carrierRows.enforced||e.carrierRows.satisfied)&&(!e.globalContinuity.enforced||e.globalContinuity.satisfied)&&
    std::abs(kcl)<=std::max(1e-12,1e-5*std::abs(id))&&std::abs(load)<=1e-6;
   diag={{"psi_inf",e.blockNorms.psi},{"electron_inf",e.blockNorms.phin},{"hole_inf",e.blockNorms.phip},
    {"electron_closure",e.globalContinuity.electron.ratio},{"hole_closure",e.globalContinuity.hole.ratio}};
  }
  nlohmann::json record={{"outer_V",target},{"seed",kind},{"converged",result.converged},{"physical_acceptance",physical},
   {"iterations",result.iterations},{"equation_inf",result.equationResidualNorm},{"load_V",result.loadLineResidual_V},
   {"failure",result.failureReason},{"reason",result.convergenceReason},{"inner_V",result.state.lambda},
   {"current_A_per_um",id},{"KCL_A_per_um",kcl},{"actual_load_V",load},{"diagnostics",diag}};
  report["attempts"].push_back(record);save();
  return std::make_pair(std::move(result),physical);
 };
 auto anchor=attempt(current,outer,"restart_anchor");if(!anchor.second)throw std::runtime_error("Own anchor failed acceptance");
 current=anchor.first.state;double step=outerControl.initial,minimum=outerControl.minimum;int retries=0;bool complete=false;
 while(report["attempts"].size()<512){
  const auto guardedStep=diagnostics::guardedOuterStep(step,outerControl,innerControl.volts,
   current.lambda-previous.lambda,outer-previousOuter);
  if(!guardedStep){report["failure"]="predicted inner-step limit requires an outer step below min_step";break;}
  const double attemptedStep=*guardedStep;
  const double target=std::min(100000.,outer+attemptedStep);
  const auto prediction=extrapolate?corrector.secantPredict(previous,current,target,1e7):current;
  const bool changed=prediction.lambda!=current.lambda || !(prediction.x.array()==current.x.array()).all();
  auto candidate=attempt(prediction,target,changed?"secant":"constant");bool fallback=false;
  if(!candidate.second&&changed){candidate=attempt(current,target,"constant_fallback");fallback=true;}
  if(candidate.second){
   const double acceptedStep=target-outer;previousOuter=outer;previous=current;current=candidate.first.state;outer=target;
   auto record=report["attempts"].back();record["fallback_used"]=fallback;report["points"].push_back(record);
   const auto s=unpack(current.x);auto vec=[](const VectorXd&v){return std::vector<double>(v.data(),v.data()+v.size());};
   std::ofstream(out/("state_"+std::to_string(report["points"].size())+".json"))<<nlohmann::json{{"outer_V",outer},{"inner_V",current.lambda},
    {"psi_V",vec(s.psi)},{"fn_V",vec(s.phin)},{"fp_V",vec(s.phip)}}.dump();
   std::cout<<"accepted outer="<<std::setprecision(17)<<outer<<" inner="<<current.lambda<<" I="<<record["current_A_per_um"]<<" N="<<candidate.first.iterations<<" fallback="<<fallback<<std::endl;
   if(std::abs(record["current_A_per_um"].get<double>())>=.0013||outer>=100000.){complete=true;break;}
   const double growth=fallback?1.:(outerControl.growthMode=="fixed"?outerControl.growth:
    1.+(outerControl.growth-1.)*std::max(0.,1.-std::max(candidate.first.iterations-1,0)/(.75*c.maxIterations)));
   step=std::clamp(acceptedStep*growth,minimum,outerControl.maximum);retries=0;save();
  }else{
   if(attemptedStep<=minimum*(1+1e-12)||++retries>outerControl.maxRetries){report["failure"]="minimum outer step or retry limit";break;}
   step=std::max(minimum,attemptedStep*outerControl.shrink);
  }
 }
 report["complete"]=complete;report["last_outer_V"]=outer;report["last_inner_V"]=current.lambda;save();return complete?0:1;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
