#include "vela/simulation/ElectrothermalSimulation.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace vela {
namespace {
using json = nlohmann::json;
namespace fs = std::filesystem;
constexpr std::array<const char*,4> stateKeys{"state_interleaved",
    "referenced_state_interleaved", "electron_qf_reference_V", "hole_qf_reference_V"};
json read(const fs::path& p) {
    std::ifstream f(p); if(!f) throw std::runtime_error("Cannot read " + p.string());
    return json::parse(f);
}
void save(const fs::path& p, const json& value) {
    auto temporary=p; temporary += ".tmp";
    { std::ofstream f(temporary); if(!f) throw std::runtime_error("Cannot write " + temporary.string());
      f << value.dump(2) << '\n'; f.flush(); if(!f) throw std::runtime_error("Failed writing " + temporary.string()); }
#ifdef _WIN32
    if(!MoveFileExW(temporary.c_str(),p.c_str(),MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("Cannot replace checkpoint " + p.string());
#else
    fs::rename(temporary,p);
#endif
}
fs::path resolve(const fs::path& base,const std::string& name) {
    fs::path p(name); return fs::absolute(p.is_absolute()?p:base/p).lexically_normal();
}
json gate(const json& r,double bias) {
    json reasons=json::array();
    if(r.value("diagnostic_stop","")!="diagnostic_scaled_residual") reasons.push_back("four_equation_convergence");
    const auto blocks=r.value("electrical_block_gates",json::array());
    if(blocks.size()!=3 || !std::all_of(blocks.begin(),blocks.end(),[](const json& b){return b.value("satisfied",false);}))
        reasons.push_back("electrical_blocks");
    if(!r.value("carrier_row_gate",json::object()).value("satisfied",false)) reasons.push_back("carrier_rows");
    std::set<std::string> names; double sum=0.,largest=1e-24;
    for(const auto& c:r.at("contacts")) {
        names.insert(c.at("contact")); double i=c.at("total_outflow_A_per_m");
        if(!std::isfinite(i)) reasons.push_back("nonfinite_current");
        sum+=i; largest=std::max(largest,std::abs(i));
    }
    if(names!=std::set<std::string>{"gate","drain","source","substrate"} || r.at("contacts").size()!=4)
        reasons.push_back("terminal_set");
    double kcl=std::abs(sum)/largest;
    if(!std::isfinite(kcl)||kcl>1e-3) reasons.push_back("kcl");
    const auto t=r.at("temperature_K").get<std::vector<double>>();
    if(t.empty()||!std::all_of(t.begin(),t.end(),[](double v){return std::isfinite(v)&&v>=50.;})) reasons.push_back("temperature");
    double source=r.at("lattice_source_W_per_m"),boundary=r.at("boundary_heat_W_per_m");
    json balance=nullptr;
    if(bias==0.) {
        if(source!=0.||boundary!=0.||!std::all_of(t.begin(),t.end(),[](double v){return v==300.;}))
            reasons.push_back("zero_power_equilibrium");
    } else {
        double scale=std::max(std::abs(source),std::abs(boundary));
        double error=scale>0.?std::abs(source-boundary)/scale:1.; balance=error;
        if(!std::isfinite(error)||error>1e-3) reasons.push_back("heat_balance");
    }
    return {{"pass_gate",reasons.empty()},{"reasons",reasons},{"kcl_relative",kcl},{"heat_balance_relative",balance}};
}
void useState(json& cfg,const json& state) {
    for(const auto* key:stateKeys) cfg[key]=state.at(key);
}
bool predict(json& cfg,const json& previous,double ratio) {
    if(!(ratio>0.&&ratio<=3.)||cfg.value("potential_origin_V",0.)!=previous.value("potential_origin_V",0.)) return false;
    auto x=cfg.at("referenced_state_interleaved").get<std::vector<double>>();
    const auto old=previous.at("referenced_state_interleaved").get<std::vector<double>>();
    if(x.size()!=old.size()||x.empty()||x.size()%4) return false;
    for(std::size_t i=0;i<x.size()/4;++i) {
        for(int k:{0,3}) x[4*i+k]+=ratio*(x[4*i+k]-old[4*i+k]);
        for(int k:{1,2}) {
            const auto* key=k==1?"electron_qf_reference_V":"hole_qf_reference_V";
            const long double change=static_cast<long double>(cfg.at(key).at(i).get<double>())-
                previous.at(key).at(i).get<double>()+x[4*i+k]-old[4*i+k];
            x[4*i+k]+=ratio*static_cast<double>(change);
        }
        if(!(x[4*i+3]>50.&&x[4*i+3]<5000.)) return false;
    }
    if(!std::all_of(x.begin(),x.end(),[](double v){return std::isfinite(v);})) return false;
    cfg["referenced_state_interleaved"]=x;
    for(std::size_t i=0;i<x.size()/4;++i) for(int k:{1,2})
        x[4*i+k]=static_cast<double>(static_cast<long double>(x[4*i+k])+cfg.at(k==1?"electron_qf_reference_V":"hole_qf_reference_V").at(i).get<double>());
    cfg["state_interleaved"]=x; return true;
}
json curveSummary(const json& r) {
    double current=0.,peak=0.,integral=0.,area=0.;json terminals=json::array();
    for(const auto& c:r.at("contacts")) {
        const double value=c.at("total_outflow_A_per_m").get<double>()*1e-6;
        const std::string name=c.at("contact");if(name=="drain")current=value;
        terminals.push_back({{"contact",name},{"current_A_per_um",value}});
    }
    for(std::size_t i=0;i<r.at("temperature_K").size();++i) {
        double t=r.at("temperature_K")[i],a=r.at("nodal_area_m2")[i];
        peak=std::max(peak,t);area+=a;integral+=a*t;
    }
    return {{"current_A_per_um",current},{"peak_K",peak},{"mean_K",integral/area},{"terminals",terminals}};
}
void curves(const fs::path& root,const json& ledger) {
    std::ofstream f(root/"curve.csv"),terminal(root/"terminal_balance.csv");
    if(!f||!terminal) throw std::runtime_error("Cannot write electrothermal curves");
    f << std::setprecision(17) << "bias_V,current_total_A_per_um,peak_temperature_K,mean_temperature_K\n";
    terminal << std::setprecision(17) << "point_index,bias_V,contact,current_total_A_per_um\n";
    std::size_t index=0;
    for(const auto& point:ledger.at("exact_points")) {
        const auto s=point.contains("curve_summary")?point.at("curve_summary"):
            curveSummary(read(point.at("result").get<std::string>()));
        for(const auto& c:s.at("terminals"))
            terminal << index << ',' << point.at("bias_V").get<double>() << ',' << c.at("contact").get<std::string>() << ',' << c.at("current_A_per_um").get<double>() << '\n';
        f << point.at("bias_V").get<double>() << ',' << s.at("current_A_per_um").get<double>() << ',' << s.at("peak_K").get<double>() << ',' << s.at("mean_K").get<double>() << '\n'; ++index;
    }
}
}

json runElectrothermalSweep(const json& deck,const fs::path& configFile) {
    const auto start=std::chrono::steady_clock::now();
    const fs::path base=fs::absolute(configFile).parent_path();
    const auto inputPath=resolve(base,deck.at("input_file"));
    const auto root=resolve(base,deck.at("output_directory"));
    auto input=read(inputPath);
    input["mesh_file"]=resolve(inputPath.parent_path(),input.at("mesh_file")).string();
    const auto mesh=read(input.at("mesh_file").get<std::string>());
    const auto control=deck.at("sweep");
    const auto biases=control.at("bias_points_V").get<std::vector<double>>();
    const double initial=control.value("initial_step_V",.1),minimum=control.value("minimum_step_V",1e-4),maximum=control.value("maximum_step_V",4./3.);
    const int budget=control.value("max_newton",60),growth=control.value("growth_newton",8);
    const auto predictor=control.value("predictor",std::string("linear"));
    const double densityBiasCeiling=control.value("density_update_maximum_bias_V",std::numeric_limits<double>::infinity());
    if(control.contains("density_update_maximum_bias_V") && !(densityBiasCeiling>=0.&&std::isfinite(densityBiasCeiling)))
        throw std::invalid_argument("Density update bias ceiling must be finite and nonnegative");
    if(control.contains("density_update_requires_prediction") && !control.at("density_update_requires_prediction").is_boolean())
        throw std::invalid_argument("Density prediction guard must be boolean");
    const bool densityNeedsPrediction=control.value("density_update_requires_prediction",false);
    if(biases.empty()||biases.front()!=0.||!std::all_of(biases.begin(),biases.end(),[](double v){return std::isfinite(v);})||
       std::adjacent_find(biases.begin(),biases.end(),std::greater_equal<double>())!=biases.end())
        throw std::invalid_argument("Bias points must increase strictly from zero");
    if(!(0.<minimum&&minimum<=initial&&initial<=maximum&&std::isfinite(maximum))||budget<1||growth<1||growth>budget)
        throw std::invalid_argument("Invalid electrothermal step/Newton budget");
    if(predictor!="linear"&&predictor!="none") throw std::invalid_argument("Unknown electrothermal predictor");
    const auto initialization=deck.value("initialization",json::object());
    const auto initializeMode=initialization.value("mode",std::string("provided_state"));
    if(initializeMode!="provided_state"&&initializeMode!="neutral_300K") throw std::invalid_argument("Unknown sweep initialization mode");
    if(initializeMode=="provided_state") for(const auto* key:stateKeys)
        if(!input.contains(key)) throw std::invalid_argument(std::string("Explicit initial state required: ")+key);
    const auto& gates=input.at("electrical_gate_solver");
    if(gates.at("carrier_row_convergence").at("mode")!="enforce"||gates.at("block_absolute_convergence").at("mode")!="enforce")
        throw std::invalid_argument("Electrothermal production sweep requires enforced row and block gates");
    std::set<std::size_t> drain;
    for(const auto& c:mesh.at("contacts")) if(c.at("name")=="drain") for(const auto& i:c.at("node_ids")) drain.insert(i.get<std::size_t>());
    if(drain.empty()) throw std::invalid_argument("Missing drain contact");
    double origin=input.value("potential_origin_V",0.);
    for(const auto& b:input.at("boundaries")) if(drain.contains(b.at("node"))&&b.at("kind")=="neutral_contact"&&b.at("value").get<double>()+origin!=0.)
        throw std::invalid_argument("Initial drain bias must be zero");
    json ledger,state=input; double current=0.,step=initial;std::size_t index=0;
    if(deck.value("resume",false)) {
        ledger=read(root/"ledger.json");
        if(read(root/"input_snapshot.json")!=input || read(root/"mesh_snapshot.json")!=mesh || ledger.at("sweep")!=control ||
           ledger.value("initialization",json::object())!=initialization)
            throw std::invalid_argument("Restart input, mesh or sweep differs from checkpoint");
        if(ledger.at("status")=="complete") return ledger;
        current=ledger.at("accepted_bias_V");step=ledger.at("next_step_V");index=ledger.at("exact_points").size();
        if(!ledger.at("accepted_result").is_null()) state=read(ledger.at("accepted_result").get<std::string>());
        else if(ledger.contains("initialized_result")) state=read(ledger.at("initialized_result").get<std::string>());
    } else {
        if(fs::exists(root)) throw std::invalid_argument("Output directory already exists; use explicit resume for a checkpoint");
        fs::create_directories(root);save(root/"input_snapshot.json",input);save(root/"mesh_snapshot.json",mesh);
        ledger={{"schema","vela.electrothermal_dc_sweep.v1"},{"scope","Explicit audited silicon electrothermal DC; reference acceptance is scored separately"},
            {"input_file",inputPath.string()},{"sweep",control},{"initialization",initialization},{"status","running"},{"runs",json::array()},
            {"exact_points",json::array()},{"accepted_bias_V",0.},{"accepted_result",nullptr},{"next_step_V",step},{"wall_seconds",0.}};
    }
    const double priorWall=ledger.value("wall_seconds",0.);
    auto checkpoint=[&] {
        ledger["next_step_V"]=step;
        ledger["wall_seconds"]=priorWall+std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        save(root/"ledger.json",ledger);
    };
    ledger["status"]="running";checkpoint();
    if(initializeMode=="neutral_300K"&&!ledger.contains("initialized_result")) {
        const double finalGate=initialization.at("gate_voltage_V");
        if(!(finalGate>0.&&std::isfinite(finalGate))) throw std::invalid_argument("Positive final gate voltage required");
        std::set<std::size_t> gateNodes;
        for(const auto& c:mesh.at("contacts")) if(c.at("name")=="gate") for(const auto& i:c.at("node_ids")) gateNodes.insert(i.get<std::size_t>());
        if(gateNodes.empty()) throw std::invalid_argument("Missing gate contact for prebias");
        const auto initRoot=root/"initialization";
        if(fs::exists(initRoot)) throw std::runtime_error("Incomplete initialization exists; preserve evidence and choose a new output directory");
        fs::create_directory(initRoot);
        json initRuns=json::array();double gateVoltage=0.;std::size_t initAttempt=0;
        const auto solveInitial=[&](double voltage,bool neutral,bool poisson) {
            auto cfg=input;
            // Initialization has no accepted drain-history predictor. Apply the
            // same opt-in guard as for unpredicted sweep steps.
            if(densityNeedsPrediction) cfg["diagnostic_density_update_iterations"]=0;
            for(auto& b:cfg.at("boundaries")) if(gateNodes.contains(b.at("node"))&&b.at("kind")=="psi")
                b["value"]=b.at("value").get<double>()-finalGate+voltage;
            cfg["solve_mode"]=poisson?"poisson":"coupled";
            cfg["diagnostic_newton_max_iterations"]=initialization.value("max_newton",100);
            if(neutral) {for(const auto* key:stateKeys) cfg.erase(key);cfg["initialization"]="neutral_300K";}
            else {for(const auto* key:stateKeys) cfg[key]=state.at(key);}
            std::ostringstream name;name<<"stage_"<<std::setw(4)<<std::setfill('0')<<initAttempt++;
            auto directory=initRoot/name.str();fs::create_directory(directory);save(directory/"input.json",cfg);
            std::ofstream log(directory/"run.log");auto result=solveElectrothermalPoint(cfg,log);save(directory/"output.json",result);
            auto acceptance=gate(result,0.);
            initRuns.push_back({{"gate_V",voltage},{"solve_mode",cfg["solve_mode"]},{"result",(directory/"output.json").string()},
                {"gate",acceptance},{"newton_updates",result.at("newton_updates")}});
            ledger["initialization_runs"]=initRuns;checkpoint();
            if(acceptance.at("pass_gate").get<bool>()) {state=result;return true;}return false;
        };
        if(!solveInitial(0.,true,true)||!solveInitial(0.,false,false)) {
            ledger["status"]="initialization_failed";checkpoint();return ledger;
        }
        // Original gate prebias: Poisson-only, relative InitialStep .01,
        // Increment 1.35 and MaxStep .4, through 4 V before 8 V.
        std::vector<double> gateGoals;if(finalGate>4.)gateGoals.push_back(4.);gateGoals.push_back(finalGate);
        for(double goal:gateGoals) {
            const double span=goal-gateVoltage;double increment=.01*span;
            while(gateVoltage<goal) {
                double target=std::min(goal,gateVoltage+increment);
                if(solveInitial(target,false,true)) {gateVoltage=target;increment=std::min(.4*span,increment*1.35);}
                else {increment=(target-gateVoltage)*.5;if(increment<1e-4*span){ledger["status"]="initialization_failed";checkpoint();return ledger;}}
            }
        }
        if(!solveInitial(finalGate,false,false)) {ledger["status"]="initialization_failed";checkpoint();return ledger;}
        ledger["initialized_result"]=initRuns.back().at("result");checkpoint();
    }
    std::size_t invocationAttempts=0;
    const auto pauseAfter=deck.value("pause_after_attempts",std::size_t(0));
    while(index<biases.size()) {
        if(fs::exists(root/"STOP")||(pauseAfter&&invocationAttempts>=pauseAfter)) {ledger["status"]="stopped_at_checkpoint";checkpoint();return ledger;}
        const double target=index==0?0.:std::min(biases[index],current+step);
        json cfg=input; useState(cfg,state);
        cfg["solve_mode"]="coupled";cfg["initialization"]="provided_state";
        cfg["diagnostic_newton_max_iterations"]=budget;
        if(target>densityBiasCeiling)cfg["diagnostic_density_update_iterations"]=0;
        for(auto& b:cfg.at("boundaries")) if(drain.contains(b.at("node"))&&b.at("kind")=="neutral_contact") b["value"]=target-origin;
        json prediction={{"mode",predictor},{"used",false}};
        if(predictor=="linear"&&target>current) {
            for(auto it=ledger["runs"].rbegin();it!=ledger["runs"].rend();++it) {
                const double previousBias=it->at("bias_V"),difference=current-previousBias;
                if(it->at("gate").at("pass_gate").get<bool>()&&difference>1e-12&&(target-current)/difference<=3.) {
                    auto previous=read(fs::path(it->at("directory").get<std::string>())/"output.json");
                    bool used=predict(cfg,previous,(target-current)/difference);
                    prediction={{"mode",predictor},{"used",used},{"prior_bias_V",previousBias},{"ratio",(target-current)/difference}};break;
                }
            }
        }
        if(densityNeedsPrediction && !prediction.at("used").get<bool>())cfg["diagnostic_density_update_iterations"]=0;
        std::ostringstream name;name<<"step_"<<std::setw(4)<<std::setfill('0')<<ledger["runs"].size();
        const auto directory=root/name.str();
        if(!fs::create_directory(directory)) throw std::runtime_error("Uncheckpointed attempt exists: "+directory.string());
        save(directory/"input.json",cfg);std::ofstream log(directory/"run.log");
        const auto solveStart=std::chrono::steady_clock::now();
        json result=json::object(),acceptance;int returncode=0;
        try {result=solveElectrothermalPoint(cfg,log);save(directory/"output.json",result);acceptance=gate(result,target);}
        catch(const std::exception& e) {returncode=1;log<<e.what()<<'\n';acceptance={{"pass_gate",false},{"reasons",{"solver_error"}},{"message",e.what()}};}
        ledger["runs"].push_back({{"parent_bias_V",current},{"bias_V",target},{"directory",directory.string()},
            {"returncode",returncode},{"gate",acceptance},{"prediction",prediction},{"newton_updates",result.value("newton_updates",0)},
            {"wall_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-solveStart).count()}});
        ++invocationAttempts;
        if(acceptance.at("pass_gate").get<bool>()) {
            current=target;state=result;ledger["accepted_bias_V"]=current;ledger["accepted_result"]=(directory/"output.json").string();
            if(std::abs(current-biases[index])<1e-12) {ledger["exact_points"].push_back({{"bias_V",biases[index]},{"result",(directory/"output.json").string()},{"curve_summary",curveSummary(result)}});++index;curves(root,ledger);}
            const int updates=result.at("newton_updates");
            if(updates<=growth) step=std::min(maximum,step*1.5);else if(updates>20) step=std::max(minimum,step*.5);
        } else {
            step=(target-current)*.5;
            if(target==0.||step<minimum) {ledger["status"]="failed";checkpoint();return ledger;}
        }
        checkpoint();
    }
    ledger["status"]="complete";checkpoint();return ledger;
}
}
