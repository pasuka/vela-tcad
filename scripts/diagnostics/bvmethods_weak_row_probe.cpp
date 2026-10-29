// Uninstalled same-population source decomposition. No production update or gate changes.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace {
using namespace vela;
using Json = nlohmann::json;
namespace fs = std::filesystem;
void matrixFile(const fs::path& path, const SparseMatrixd& a) {
    std::ofstream f(path);
    f << std::setprecision(17) << "row,col,value\n";
    for (int k = 0; k < a.outerSize(); ++k)
        for (SparseMatrixd::InnerIterator v(a, k); v; ++v)
            f << v.row() << ',' << v.col() << ',' << v.value() << '\n';
}
}

int main(int argc, char** argv) { try {
    if (argc != 3) throw std::invalid_argument("manifest.json output-directory");
    Json m, original;
    { std::ifstream f(argv[1]); f >> m; }
    const std::string config = m.at("config");
    { std::ifstream f(config); f >> original; }
    const fs::path out = argv[2]; fs::create_directories(out);
    Json summary = Json::array();
    for (const auto& frame : m.at("frames")) {
        for (const std::string name : {"full", "noav", "nosrh", "transport"}) {
            auto cfg = original;
            if (name == "noav" || name == "transport") cfg["solver"]["impact_ionization"]["model"] = "none";
            if (name == "nosrh" || name == "transport") cfg["solver"]["recombination"] = Json::array();
            const double inner = frame.at("inner_V"), outer = frame.at("outer_V");
            for (auto& c : cfg["contacts"]) if (c["name"] == "drain") c["bias"] = inner;
            cfg = canonicalizeDeck(cfg);
            auto problem = loadNewtonProblem(config, cfg);
            DDStateCliScope scope(problem, cfg, configDirectory(config));
            auto solver = makeNewtonSolver(problem);
            const int n = problem.mesh.numNodes(), d = 3 * n;
            const auto& units = problem.newton.inputScaling.unitSystem();
            const auto cm = [&](double x) { return units.internalConcentrationToM3(x) / 1e6; };
            DDSolution s;
            s.psi.resize(n); s.phin = s.phip = VectorXd::Zero(n); s.n = s.p = VectorXd::Zero(n);
            const auto& fields = frame.at("fields");
            for (int i = 0; i < n; ++i) s.psi[i] = fields.at("psi_V")[i];
            for (const auto& v : fields.at("electron")) s.n[v[0].get<int>()] = units.m3ToInternalConcentration(v[1].get<double>() * 1e6);
            for (const auto& v : fields.at("hole")) s.p[v[0].get<int>()] = units.m3ToInternalConcentration(v[1].get<double>() * 1e6);
            for (const auto& t : solver.evaluatePoissonTerms(s).rows) {
                s.phin[t.nodeId] = t.electronRequiredQuasiFermiPotential;
                s.phip[t.nodeId] = t.holeRequiredQuasiFermiPotential;
            }
            const VectorXd x = solver.packArclengthState(s);
            auto unpack = solver.makeArclengthStateUnpacker(); s = unpack(x);
            auto sys = solver.makeArclengthSystem("drain");
            const double scale = solver.evaluateResidual(s).potentialScale;
            const double vt = constants::kb * problem.newton.temperature_K / constants::q;
            const auto folder = out / frame.at("label").get<std::string>() / name;
            fs::create_directories(folder);
            const VectorXd residual = sys.residual(x, inner);
            const SparseMatrixd jacobian = sys.jacobian(x, inner);
            if (!residual.allFinite()) throw std::runtime_error("nonfinite baseline residual");
            const auto fv = sys.parameterDerivative(x, inner);
            auto port = solver.makeArclengthContactCurrentFunctional("drain", units.currentDensityAM2PerInternal() * units.lengthMPerInternal() * 1e-6);
            const auto terminal = port.linearize(x, inner);
            std::vector<Eigen::Triplet<double>> entries;
            for (int k = 0; k < jacobian.outerSize(); ++k)
                for (SparseMatrixd::InnerIterator v(jacobian, k); v; ++v) entries.emplace_back(v.row(), v.col(), v.value());
            for (int i = 0; i < d; ++i) {
                if (fv[i] != 0.) entries.emplace_back(i, d, fv[i]);
                if (terminal.stateDerivative[i] != 0.) entries.emplace_back(d, i, 1e7 * terminal.stateDerivative[i]);
            }
            entries.emplace_back(d, d, 1 + 1e7 * terminal.parameterDerivative);
            SparseMatrixd bordered(d + 1, d + 1); bordered.setFromTriplets(entries.begin(), entries.end());
            matrixFile(folder / "matrix.csv", bordered);
            std::ofstream rf(folder / "rhs.csv"); rf << std::setprecision(17) << "row,rhs\n";
            for (int i = 0; i < d; ++i) rf << i << ',' << -residual[i] << '\n';
            rf << d << ',' << -(inner + 1e7 * terminal.value - outer) << '\n';
            std::ofstream state(folder / "fields.csv");
            state << std::setprecision(17) << "node,psi_V,fn_V,fp_V,n_cm3,p_cm3,Nc_cm3,Nv_cm3,dn_dV_cm3,dp_dV_cm3\n";
            for (const auto& t : solver.evaluatePoissonTerms(s).rows) {
                const int i = t.nodeId; double dn = 0., dp = 0.;
                if (t.intrinsicDensity > 0.) {
                    dn = experimental::JoyceDixonStatistics::densitySlopeFromPopulation(cm(s.n[i]), cm(t.electronDensityOfStates)) / vt;
                    dp = experimental::JoyceDixonStatistics::densitySlopeFromPopulation(cm(s.p[i]), cm(t.holeDensityOfStates)) / vt;
                }
                state << i << ',' << s.psi[i] << ',' << s.phin[i] << ',' << s.phip[i] << ',' << cm(s.n[i]) << ',' << cm(s.p[i]) << ',' << cm(t.electronDensityOfStates) << ',' << cm(t.holeDensityOfStates) << ',' << dn << ',' << dp << '\n';
            }
            std::ofstream ledger(folder / "ledger.csv");
            ledger << std::setprecision(17) << "node,e_flux,h_flux,e_srh,h_srh,e_AV,h_AV,e_residual,h_residual\n";
            for (const auto& t : solver.evaluateCarrierTermDiagnostics(s, true).rows)
                ledger << t.nodeId << ',' << t.electronFlux << ',' << t.holeFlux << ',' << t.electronRecombination << ',' << t.holeRecombination << ',' << t.electronImpact << ',' << t.holeImpact << ',' << t.electronResidual << ',' << t.holeResidual << '\n';
            Json jv = Json::array();
            for (int block : {0, 1, 2}) for (int node : {308, 311, 312}) {
                const int col = block * n + node;
                for (double h : {1e-4, 1e-6, 1e-8}) {
                    VectorXd trial = x; trial[col] += h / scale;
                    const VectorXd diff = (sys.residual(trial, inner) - residual) / h;
                    for (int rowNode : {307, 308, 309, 311, 312, 318}) {
                        const int row = 2 * n + rowNode;
                        jv.push_back({{"column_block", block}, {"column_node", node}, {"row_node", rowNode}, {"h_V", h}, {"analytic_per_V", jacobian.coeff(row, col) / scale}, {"one_sided_per_V", diff[row]}});
                    }
                }
            }
            Json identity = {{"label", frame.at("label")}, {"arm", name}, {"potential_scale_V", scale}, {"equation_inf", residual.lpNorm<Eigen::Infinity>()}, {"inner_V", inner}, {"outer_V", outer}, {"current_A_per_um", terminal.value}};
            std::ofstream(folder / "identity.json") << identity.dump(2) << '\n';
            std::ofstream(folder / "one_sided.json") << jv.dump(2) << '\n';
            summary.push_back(identity); std::cout << identity.dump() << std::endl;
        }
    }
    std::ofstream(out / "summary.json") << summary.dump(2) << '\n';
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
