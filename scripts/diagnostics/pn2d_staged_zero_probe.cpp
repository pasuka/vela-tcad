// Diagnostic only: independent Poisson solve, then zero-bias coupled circuit.
#define main embeddedRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"

int main(int argc, char** argv) {
    try {
        using namespace vela;
        if (argc != 3) throw std::runtime_error("config output_directory required");
        nlohmann::json cfg; {std::ifstream f(argv[1]); f >> cfg;}
        cfg = canonicalizeDeck(cfg);
        auto problem = loadNewtonProblem(argv[1], cfg);
        DDStateCliScope scope(problem, cfg, configDirectory(argv[1]));
        auto owner = makeNewtonSolver(problem);
        const std::string contact = cfg.value("staged_contact", std::string("Anode"));
        const auto units = problem.newton.inputScaling.unitSystem();
        const int n = problem.mesh.numNodes();
        std::filesystem::path out = argv[2];
        std::filesystem::create_directories(out);
        const auto fields = [&](const std::string& name, const DDSolution& state) {
            std::ofstream f(out / (name + ".csv"));
            f << std::setprecision(17) << "node,psi_V,fn_V,fp_V,n_cm3,p_cm3\n";
            for (int i=0;i<n;++i) f << i << ',' << state.psi[i] << ',' << state.phin[i]
                << ',' << state.phip[i] << ',' << units.internalConcentrationToM3(state.n[i])/1e6
                << ',' << units.internalConcentrationToM3(state.p[i])/1e6 << '\n';
        };
        // No external/native solution is used to initialize either stage.
        auto initialization = owner.buildPoissonBlockInitialization();
        fields("cold", initialization.coldInitial);
        auto poisson = owner.solvePoissonOnly(initialization.coldInitial);
        fields("poisson", poisson.solution);
        nlohmann::json result;
        result["poisson"] = {{"converged",poisson.converged},{"iterations",poisson.iters},
            {"initial_residual",poisson.initialResidualNorm},{"final_residual",poisson.finalResidualNorm},
            {"reason",poisson.convergenceReason},{"failure",poisson.failureDiagnostics.failureReason}};
        result["poisson"]["history"] = nlohmann::json::array();
        for (const auto& h:poisson.history) result["poisson"]["history"].push_back(
            {{"iteration",h.iter},{"residual",h.residualNorm},{"damping",h.dampingFactor}});
        auto sys = owner.makeArclengthSystem(contact);
        if (cfg.contains("poisson_reference_psi")) {
            for (const auto& reference:cfg["poisson_reference_psi"]) {
                DDSolution state = initialization.coldInitial;
                state.packedState.resize(0);state.packedLow.resize(0);
                state.electronQfReference.resize(0);state.holeQfReference.resize(0);
                state.phinIncrement.resize(0);state.phipIncrement.resize(0);
                state.phin.setZero();state.phip.setZero();
                for(int i=0;i<n;++i)state.psi[i]=reference["psi_V"][i];
                auto x=owner.packArclengthState(state);
                auto J=sys.jacobian(x,0.);auto F=sys.residual(x,0.);
                const std::string label=reference["label"];
                std::ofstream m(out/(label+"_matrix.csv")),r(out/(label+"_rhs.csv"));
                m<<std::setprecision(17)<<"row,col,value\n";
                r<<std::setprecision(17)<<"row,rhs\n";
                for(int k=0;k<J.outerSize();++k)for(SparseMatrixd::InnerIterator it(J,k);it;++it)
                    if(it.row()<n && it.col()<n)m<<it.row()<<','<<it.col()<<','<<it.value()<<'\n';
                for(int i=0;i<n;++i)r<<i<<','<<-F[i]<<'\n';
            }
        }
        if (poisson.converged) {
            auto x=owner.packArclengthState(poisson.solution);
            const double scale=units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6;
            auto current=owner.makeArclengthContactCurrentFunctional(contact,scale);
            if (cfg.value("staged_export_coupled_operator", false)) {
                auto J=sys.jacobian(x,0.);auto F=sys.residual(x,0.);
                auto border=sys.parameterDerivative(x,0.);
                auto terminal=current.linearize(x,0.);
                std::ofstream m(out/"D0_matrix.csv"),r(out/"D0_rhs.csv");
                m<<std::setprecision(17)<<"row,col,value\n";
                r<<std::setprecision(17)<<"row,rhs\n";
                const int size=x.size();
                for(int k=0;k<J.outerSize();++k)for(SparseMatrixd::InnerIterator it(J,k);it;++it)
                    m<<it.row()<<','<<it.col()<<','<<it.value()<<'\n';
                for(int i=0;i<size;++i) {
                    if(border[i]!=0.)m<<i<<','<<size<<','<<border[i]<<'\n';
                    if(terminal.stateDerivative[i]!=0.)m<<size<<','<<i<<','<<1e7*terminal.stateDerivative[i]<<'\n';
                    r<<i<<','<<-F[i]<<'\n';
                }
                m<<size<<','<<size<<','<<1.+1e7*terminal.parameterDerivative<<'\n';
                r<<size<<','<<-1e7*terminal.value<<'\n';
                result["coupled_operator_dimension"]=size+1;
            }
            CoupledLoadLineConfig cc;
            const auto circuit=cfg["sweep"]["external_circuit"];
            cc.linearSolver=circuit.value("coupled_linear_solver",std::string("direct_bordered"));
            cc.equationTolerance=circuit.value("coupled_equation_tolerance",1e-8);
            cc.loadLineTolerance_V=circuit.value("residual_tolerance_V",1e-6);
            CoupledLoadLineNewton coupled(sys,current.value,cc,current);
            auto solved=coupled.solve({x,0.},0.,1e7);
            auto state=owner.makeArclengthStateUnpacker()(solved.state.x);
            fields("coupled",state);
            auto gate=owner.makeArclengthConvergenceEvaluator(contact)(solved.state.x,solved.state.lambda);
            result["coupled"]={{"scalar_converged",solved.converged},{"iterations",solved.iterations},
                {"failure",solved.failureReason},{"device_residual",solved.equationResidualNorm},
                {"load_residual_V",solved.loadLineResidual_V},{"inner_V",solved.state.lambda},
                {"current_A_per_um",solved.directedCurrent_A_per_um},
                {"carrier_rows",carrierRowConvergenceJson(gate.carrierRows)},
                {"global_continuity",globalContinuityClosureJson(gate.globalContinuity)},
                {"physical_accepted",solved.converged && gate.carrierRows.satisfied && gate.globalContinuity.satisfied}};
            // Diagnostic device-only reclosure distinguishes circuit effects.
            auto device=owner.solve(poisson.solution);
            fields("device_coupled",device.solution);
            result["device_coupled"]={{"converged",device.converged},{"iterations",device.iters},
                {"reason",device.convergenceReason},{"failure",device.failureDiagnostics.failureReason}};
        }
        result["n_nodes"]=n;result["R_ohm_um"]=1e7;result["outer_V"]=0.;
        std::ofstream(out/"summary.json")<<result.dump(2);
        std::cout<<result.dump(2)<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
