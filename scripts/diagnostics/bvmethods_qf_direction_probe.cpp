// Bounded, uninstalled same-state QF-cap and bordered-direction audit.
// It reads Vela's own state; no native populations, sources or matrices enter.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"
#include <filesystem>
#include <set>

namespace {
using namespace vela;
using Json = nlohmann::json;
namespace fs = std::filesystem;
void saveJson(const fs::path& path, const Json& value) {
    std::ofstream(path) << value.dump(2) << '\n';
}
void saveMatrix(const fs::path& path, const SparseMatrixd& matrix) {
    std::ofstream f(path); f << std::setprecision(17) << "row,col,value\n";
    for (int k=0;k<matrix.outerSize();++k)
        for (SparseMatrixd::InnerIterator v(matrix,k);v;++v)
            f << v.row() << ',' << v.col() << ',' << v.value() << '\n';
}
Real inf(const VectorXd& v) { return v.lpNorm<Eigen::Infinity>(); }
}

int main(int argc,char** argv) { try {
    if (argc!=3) throw std::invalid_argument("manifest.json output-directory");
    Json m, original;
    {std::ifstream f(argv[1]);f>>m;}
    {std::ifstream f(m.at("config").get<std::string>());f>>original;}
    const fs::path out=argv[2];fs::create_directories(out);
    const Real target=m.at("outer_V"), resistance=original.at("sweep").at("external_circuit").at("resistance_ohm_um");
    Json summary={{"schema","vela.bvmethods_qf_direction_audit.v1"},{"outer_V",target},
        {"resistance_ohm_um",resistance},{"arms",Json::array()}};
    for (const auto& arm:m.at("arms")) {
        const std::string name=arm.at("name");const fs::path folder=out/name;fs::create_directories(folder);
        auto cfg=original;
        // Same-state operator ablations only; never a production physics fix.
        if(arm.value("disable_av",false))cfg["solver"]["impact_ionization"]["model"]="none";
        if(arm.value("disable_srh",false))cfg["solver"]["recombination"]=Json::array();
        cfg["solver"]["quasi_fermi_update_limit_V"]=arm.at("global_limit_V");
        cfg["solver"]["quasi_fermi_update_limit_minority_V"]=arm.at("minority_limit_V");
        Real initial=m.at("inner_V");
        for(auto& contact:cfg["contacts"])if(contact["name"]=="drain")contact["bias"]=initial;
        cfg=canonicalizeDeck(cfg);
        auto p=loadNewtonProblem(m.at("config").get<std::string>(),cfg);
        DDStateCliScope stateScope(p,cfg,configDirectory(m.at("config").get<std::string>()));
        auto solver=makeNewtonSolver(p);const int n=static_cast<int>(p.mesh.numNodes()), d=3*n;
        ArclengthState start;
        if(m.contains("packed_state_file")) {
            Json state;{std::ifstream f(m.at("packed_state_file").get<std::string>());f>>state;}
            const auto values=state.at("x").get<std::vector<double>>();
            if(static_cast<int>(values.size())!=d)throw std::invalid_argument("packed state shape");
            start.x=Eigen::Map<const VectorXd>(values.data(),d);start.lambda=state.at("inner_V");
        } else {
            const auto state=readDDSolutionState(m.at("state_file").get<std::string>(),n,p.newton.inputScaling);
            start={solver.packArclengthState(state),initial};
        }
        if(!start.x.allFinite()||!std::isfinite(start.lambda))throw std::invalid_argument("nonfinite seed");
        auto sys=solver.makeArclengthSystem("drain");
        auto unpack=solver.makeArclengthStateUnpacker();
        const Real scale=solver.evaluateResidual(unpack(start.x)).potentialScale;
        const auto& units=p.newton.inputScaling.unitSystem();
        const Real line=units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6;
        auto port=solver.makeArclengthContactCurrentFunctional("drain",line);
        if(m.contains("trial_direction_file")){
            Json trialInput;{std::ifstream f(m.at("trial_direction_file").get<std::string>());f>>trialInput;}
            const auto values=trialInput.at("delta").get<std::vector<Real>>();
            if(static_cast<int>(values.size())!=d)throw std::invalid_argument("trial direction shape");
            const VectorXd delta=Eigen::Map<const VectorXd>(values.data(),d);
            Json trials=Json::array();
            for(Real alpha:trialInput.at("alphas").get<std::vector<Real>>()){
                const VectorXd x=start.x+alpha*delta;
                const Real v=start.lambda+alpha*trialInput.at("delta_inner_V").get<Real>();
                Json record={{"alpha",alpha},{"inner_V",v}};
                try{
                    const auto s=unpack(x);record["p273"]=s.p[273];record["fp273_V"]=s.phip[273];
                    const auto f=sys.residual(x,v);record["residual_finite"]=f.allFinite();record["equation_inf"]=inf(f);
                    const Real id=port.value(x,v);record["current_A_per_um"]=id;record["current_finite"]=std::isfinite(id);
                    record["load_V"]=v+resistance*id-target;
                }catch(const std::exception& e){record["exception"]=e.what();}
                trials.push_back(record);
            }
            const auto ledger=solver.evaluateCarrierTermDiagnostics(unpack(start.x),true);
            std::ofstream f(folder/"seed_ledger.csv");f<<std::setprecision(17)<<"node,hole_flux,hole_SRH,hole_AV,hole_rhs,electron_flux,electron_SRH,electron_AV,electron_rhs\n";
            for(const auto& a:ledger.rows)f<<a.nodeId<<','<<a.holeFlux<<','<<a.holeRecombination<<','<<a.holeImpact<<','<<a.holeResidual<<','
                <<a.electronFlux<<','<<a.electronRecombination<<','<<a.electronImpact<<','<<a.electronResidual<<'\n';
            saveJson(folder/"trial_replay.json",trials);summary["arms"].push_back({{"name",name},{"trials",trials}});
            saveJson(out/"summary.json",summary);continue;
        }
        Json invalid=Json::array();
        auto saveState=[&](const fs::path& dir,const VectorXd& x,Real v){
            fs::create_directories(dir);
            saveJson(dir/"packed_state.json",{{"x",std::vector<Real>(x.data(),x.data()+x.size())},
                {"inner_V",v},{"potential_scale_V",scale}});
            const auto state=unpack(x);
            if(m.value("write_state_fields",false))
                writeDDSolutionState(dir/"state.h5",state,p.newton.inputScaling);
            std::ofstream f(dir/"populations.csv");f<<std::setprecision(17)<<"node,psi,fn,fp,n,p\n";
            for(int i=0;i<n;++i)f<<i<<','<<state.psi[i]<<','<<state.phin[i]<<','<<state.phip[i]<<','
                <<state.n[i]<<','<<state.p[i]<<'\n';
        };
        std::vector<ArclengthState> checkpoints;
        bool seedJacobianSaved=false,seedPortSaved=false;
        const auto jacobian=sys.jacobian;
        sys.jacobian=[&](const VectorXd& x, Real v) {
            if(checkpoints.empty()||checkpoints.back().lambda!=v||
               !(checkpoints.back().x.array()==x.array()).all())checkpoints.push_back({x,v});
            auto J=jacobian(x,v);Json bad=Json::array();
            for(int col=0;col<J.outerSize();++col)for(SparseMatrixd::InnerIterator a(J,col);a;++a)
                if(!std::isfinite(a.value()))bad.push_back({a.row(),a.col()});
            if(!seedJacobianSaved){
                const auto dir=folder/"seed_operator";saveState(dir,x,v);saveMatrix(dir/"device_matrix.csv",J);
                const VectorXd f=sys.residual(x,v);std::ofstream csv(dir/"rhs.csv");csv<<std::setprecision(17)<<"row,rhs\n";
                for(int i=0;i<f.size();++i)csv<<i<<','<<-f[i]<<'\n';
                saveJson(dir/"matrix_status.json",{{"rows",J.rows()},{"cols",J.cols()},{"nnz",J.nonZeros()},
                    {"nonfinite_entries",bad.size()},{"residual_finite",f.allFinite()},{"residual_inf",inf(f)}});
                seedJacobianSaved=true;
            }
            if(!bad.empty()){
                const auto dir=folder/("invalid_jacobian_"+std::to_string(invalid.size()));
                saveState(dir,x,v);saveMatrix(dir/"device_matrix.csv",J);
                invalid.push_back({{"kind","device_jacobian"},{"rows",J.rows()},{"cols",J.cols()},
                    {"bad_entries",bad},{"directory",dir.string()}});
            }
            return J;
        };
        const auto linearize=port.linearize;
        port.linearize=[&](const VectorXd& x,Real v){
            auto a=linearize(x,v);
            if(!seedPortSaved){
                const auto dir=folder/"seed_operator";fs::create_directories(dir);
                std::ofstream f(dir/"port_gradient.csv");f<<std::setprecision(17)<<"row,derivative\n";
                for(int i=0;i<a.stateDerivative.size();++i)f<<i<<','<<a.stateDerivative[i]<<'\n';
                saveJson(dir/"port_status.json",{{"size",a.stateDerivative.size()},
                    {"state_derivative_finite",a.stateDerivative.allFinite()},
                    {"parameter_derivative",a.parameterDerivative},{"value",a.value}});seedPortSaved=true;
            }
            if(a.stateDerivative.size()!=d||!a.stateDerivative.allFinite()||!std::isfinite(a.parameterDerivative)){
                const auto dir=folder/("invalid_port_"+std::to_string(invalid.size()));saveState(dir,x,v);
                Json bad=Json::array();
                for(int i=0;i<a.stateDerivative.size();++i)if(!std::isfinite(a.stateDerivative[i]))bad.push_back(i);
                invalid.push_back({{"kind","port_derivative"},{"size",a.stateDerivative.size()},
                    {"bad_entries",bad},{"parameter_derivative",a.parameterDerivative},{"directory",dir.string()}});
            }
            return a;
        };
        CoupledLoadLineConfig c;
        c.maxIterations=m.value("max_iterations",40);c.equationTolerance=1e-8;c.loadLineTolerance_V=1e-6;
        c.maxVoltageUpdate_V=.1;c.linearSolver=m.value("linear_solver",std::string{"direct_bordered"});c.lineSearchMode="merit";
        c.filterEnvelopeFactor=1e12;c.applyDeviceUpdateLimit=true;
        c.updateConvergence.mode="off";c.linearizationAuditEnabled=true;
        CoupledLoadLineNewton corrector(sys,port.value,c,port);
        const auto result=corrector.solve(start,target,resistance);
        Json report={{"name",name},{"global_limit_V",arm.at("global_limit_V")},
            {"minority_limit_V",arm.at("minority_limit_V")},{"potential_scale_V",scale},
            {"converged",result.converged},{"iterations",result.iterations},
            {"line_search_trials",result.lineSearchAttempts},{"failure",result.failureReason},
            {"equation_inf",result.equationResidualNorm},{"load_V",result.loadLineResidual_V},
            {"inner_V",result.state.lambda},{"trace",Json::array()},{"snapshots",Json::array()}};
        report["invalid_linearizations"]=invalid;
        saveState(folder/"final_state",result.state.x,result.state.lambda);
        int lastDirection=-1;
        std::ofstream directions(folder/"directions.csv");directions<<std::setprecision(17)
            <<"iteration,node,raw_delta_psi_V,raw_delta_fn_V,raw_delta_fp_V\n";
        for(std::size_t k=0;k<result.trace.size();++k) {
            const auto& t=result.trace[k];Json trials=Json::array();
            for(const auto& v:t.trials)trials.push_back({{"alpha",v.alpha},{"accepted",v.accepted},
                {"predicted_device_inf",inf(v.predictedDeviceResidual)},
                {"actual_device_inf",inf(v.actualDeviceResidual)},
                {"predicted_load_V",v.predictedLoadResidual_V},{"actual_load_V",v.actualLoadResidual_V},
                {"applied_delta_inner_V",v.appliedVoltageUpdate_V}});
            Json row={{"iteration",t.iteration},{"equation_inf",t.equationResidualNorm},
                {"load_V",t.loadLineResidual_V},{"alpha",t.acceptedAlpha},
                {"raw_delta_inner_V",t.rawVoltageUpdate_V},{"limited_delta_inner_V",t.proposedVoltageUpdate_V},
                {"linear_device_defect_inf",t.deviceLinearizationResidualNorm},
                {"linear_circuit_defect_V",t.circuitLinearizationResidual},
                {"componentwise_backward_error",t.componentwiseBackwardError},
                {"audit_current_Jv_relative_error",t.currentDirectionalRelativeError},
                {"audit_device_Jv_relative_error",t.deviceJvpRelativeError},{"trials",trials}};
            if(t.rawStateUpdate.size()==d) {
                lastDirection=static_cast<int>(k);
                row["raw_delta_psi_inf_V"]=inf(t.rawStateUpdate.head(n))*scale;
                row["raw_delta_fn_inf_V"]=inf(t.rawStateUpdate.segment(n,n))*scale;
                row["raw_delta_fp_inf_V"]=inf(t.rawStateUpdate.tail(n))*scale;
                for(int i=0;i<n;++i)directions<<k<<','<<i<<','<<t.rawStateUpdate[i]*scale<<','
                    <<t.rawStateUpdate[n+i]*scale<<','<<t.rawStateUpdate[2*n+i]*scale<<'\n';
            }
            report["trace"].push_back(row);
        }
        report["has_raw_direction"]=lastDirection>=0;
        const auto snapshotIterations=lastDirection<0?std::set<int>{}:std::set<int>{0,lastDirection};
        for(int k:snapshotIterations) {
            if(k>=static_cast<int>(checkpoints.size()))throw std::runtime_error("state/trace mismatch");
            const auto& state=checkpoints[k];const auto& t=result.trace[k];
            const auto J=jacobian(state.x,state.lambda);const auto fv=sys.parameterDerivative(state.x,state.lambda);
            const auto current=port.linearize(state.x,state.lambda);
            std::vector<Eigen::Triplet<Real>> entries;
            for(int col=0;col<J.outerSize();++col)for(SparseMatrixd::InnerIterator v(J,col);v;++v)
                entries.emplace_back(v.row(),v.col(),v.value());
            for(int i=0;i<d;++i) {
                if(fv[i]!=0)entries.emplace_back(i,d,fv[i]);
                if(current.stateDerivative[i]!=0)entries.emplace_back(d,i,resistance*current.stateDerivative[i]);
            }
            entries.emplace_back(d,d,1+resistance*current.parameterDerivative);
            SparseMatrixd M(d+1,d+1);M.setFromTriplets(entries.begin(),entries.end());
            VectorXd F(d+1),delta(d+1);F.head(d)=sys.residual(state.x,state.lambda);
            F[d]=state.lambda+resistance*current.value-target;
            delta.head(d)=t.rawStateUpdate;delta[d]=t.rawVoltageUpdate_V;
            const fs::path snapshot=folder/("iteration_"+std::to_string(k));fs::create_directories(snapshot);
            saveMatrix(snapshot/"matrix.csv",M);
            std::ofstream f(snapshot/"system.csv");f<<std::setprecision(17)<<"row,residual,rhs,state,raw_delta\n";
            for(int i=0;i<=d;++i)f<<i<<','<<F[i]<<','<<-F[i]<<','<<(i<d?state.x[i]:state.lambda)<<','<<delta[i]<<'\n';
            saveJson(snapshot/"packed_state.json",{{"x",std::vector<Real>(state.x.data(),state.x.data()+d)},
                {"inner_V",state.lambda},{"potential_scale_V",scale}});
            Json meta={{"iteration",k},{"dimension",d+1},{"nnz",M.nonZeros()},
                {"linear_defect_inf",inf(M*delta+F)},{"FD",Json::array()}};
            VectorXd diagonal=J.diagonal().cwiseAbs().cwiseMax(1e-30);
            for(int mode=0;mode<2;++mode) {
                VectorXd v=VectorXd::Zero(d);Real dv=0;
                if(mode==0) {
                    Real magnitude=std::max(inf(delta.head(d))*scale,std::abs(delta[d]));
                    v=delta.head(d)*(.001/magnitude);dv=delta[d]*(.001/magnitude);
                } else v[2*n+m.value("focus_node",342)]=.001/scale;
                const VectorXd an=J*v+fv*dv;
                const Real ia=current.stateDerivative.dot(v)+current.parameterDerivative*dv;
                for(Real h:{.1,.01,.001,.0001}) {
                    const VectorXd plus=sys.residual(state.x+h*v,state.lambda+h*dv);
                    const VectorXd minus=sys.residual(state.x-h*v,state.lambda-h*dv);
                    const VectorXd central=(plus-minus)/(2*h), forward=(plus-F.head(d))/h;
                    const Real ic=(port.value(state.x+h*v,state.lambda+h*dv)-port.value(state.x-h*v,state.lambda-h*dv))/(2*h);
                    const auto normalized=[&](const VectorXd& value){return value.cwiseQuotient(diagonal).eval();};
                    meta["FD"].push_back({{"mode",mode==0?std::string{"raw_newton"}:"fp_node"+std::to_string(m.value("focus_node",342))},{"h",h},
                        {"maximum_physical_perturbation_V",h*.001},
                        {"central_relative_L2",normalized(central-an).norm()/std::max(normalized(an).norm(),Real{1e-30})},
                        {"forward_relative_L2",normalized(forward-an).norm()/std::max(normalized(an).norm(),Real{1e-30})},
                        {"current_analytic",ia},{"current_central",ic},
                        {"current_relative_error",std::abs(ia-ic)/std::max({std::abs(ia),std::abs(ic),Real{1e-30}})}});
                }
            }
            saveJson(snapshot/"audit.json",meta);report["snapshots"].push_back(meta);
        }
        const auto evaluation=solver.makeArclengthConvergenceEvaluator("drain")(result.state.x,result.state.lambda);
        Real kcl=0;bool finite=true;
        for(const auto& contact:p.mesh.contacts()) {
            const Real value=solver.makeArclengthContactCurrentFunctional(contact.name,line).value(result.state.x,result.state.lambda);
            kcl+=value;finite=finite&&std::isfinite(value);
        }
        const Real id=port.value(result.state.x,result.state.lambda);
        report["current_A_per_um"]=id;report["KCL_A_per_um"]=kcl;
        const bool finiteDiagnostics=std::isfinite(evaluation.blockNorms.psi)&&std::isfinite(evaluation.blockNorms.phin)&&
            std::isfinite(evaluation.blockNorms.phip)&&std::isfinite(evaluation.carrierRows.maxRatio)&&
            std::isfinite(evaluation.globalContinuity.electron.ratio)&&std::isfinite(evaluation.globalContinuity.hole.ratio);
        report["physical_acceptance"]=result.converged&&finite&&finiteDiagnostics&&std::abs(kcl)<=std::max(Real{1e-12},1e-5*std::abs(id))&&
            std::abs(result.state.lambda+resistance*id-target)<=1e-6&&
            (!evaluation.carrierRows.enforced||evaluation.carrierRows.satisfied)&&
            (!evaluation.globalContinuity.enforced||evaluation.globalContinuity.satisfied);
        saveJson(folder/"summary.json",report);summary["arms"].push_back(report);
        saveJson(out/"summary.json",summary);
        std::cout<<name<<" converged="<<result.converged<<" iterations="<<result.iterations
            <<" equation="<<result.equationResidualNorm<<" load="<<result.loadLineResidual_V<<std::endl;
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
