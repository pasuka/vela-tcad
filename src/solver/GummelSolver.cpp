#include "vela/solver/GummelSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/DDAssembler.h"
#include "vela/solver/LinearSolver.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/io/VTKWriter.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace vela {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

/// Thermal voltage at temperature T [K].
inline double thermalVoltage(double T = 300.0)
{
    return constants::kb * T / constants::q;
}

/// Charge-neutral equilibrium electron concentration from net doping and ni.
/// Solves: n - p = Ndop,  n*p = ni²  →  n = Ndop/2 + sqrt((Ndop/2)² + ni²)
inline double nEq(double Ndop, double ni)
{
    const double half = 0.5 * Ndop;
    return half + std::sqrt(half * half + ni * ni);
}

/// Compute per-node ni from the material database.
std::vector<double> buildNiVector(const DeviceMesh&       mesh,
                                  const MaterialDatabase& matdb)
{
    const Index N = mesh.numNodes();
    std::vector<double> ni_v(N, 0.0);

    // For each node, search cells that contain it and pick the material ni.
    // For interface nodes, take the first match (sufficient for prototype).
    std::vector<bool> found(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        double      ni_mat = 0.0;
        if (matdb.hasMaterial(region.material))
            ni_mat = matdb.getMaterial(region.material).ni;
        for (Index nid : cell.node_ids) {
            if (!found[nid]) {
                ni_v[nid]   = ni_mat;
                found[nid]  = true;
            }
        }
    }
    return ni_v;
}

} // anonymous namespace


GummelConfig gummelConfigFromJson(const nlohmann::json& json)
{
    GummelConfig cfg;
    cfg.maxIter = json.value("max_iter", cfg.maxIter);
    cfg.reltol = json.value("reltol", cfg.reltol);
    cfg.dampingPsi = json.value("damping_psi", cfg.dampingPsi);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    cfg.mobility = json.value("mobility", cfg.mobility);

    if (json.contains("recombination")) {
        const auto& value = json.at("recombination");
        if (value.is_array())
            cfg.recombination = value.get<std::vector<std::string>>();
        else if (value.is_string())
            cfg.recombination = {value.get<std::string>()};
        else
            throw std::invalid_argument(
                "gummelConfigFromJson: recombination must be a string or string array.");
    }

    return cfg;
}

// ---------------------------------------------------------------------------
// runGummel
// ---------------------------------------------------------------------------

namespace {

DDSolution runGummelImpl(const DeviceMesh&                          mesh,
                         const MaterialDatabase&                     matdb,
                         const DopingModel&                          doping,
                         const std::unordered_map<std::string, Real>& contactBiases,
                         const GummelConfig&                          cfg,
                         const DDSolution*                           initialGuess)
{
    const Index  N   = mesh.numNodes();
    const double Vt  = thermalVoltage();

    // Per-node ni
    const std::vector<double> ni_v = buildNiVector(mesh, matdb);

    // ------------------------------------------------------------------
    // Build Ohmic contact Dirichlet BCs
    // ------------------------------------------------------------------
    //   psi_contact  = V_bias + Vt * ln(n_eq / ni)  (built-in potential)
    //   n_contact    = n_eq
    //   p_contact    = ni² / n_eq
    //   phin_contact = V_bias
    //   phip_contact = V_bias
    // ------------------------------------------------------------------
    std::unordered_map<Index, Real> psiBC;
    std::unordered_map<Index, Real> nBC;
    std::unordered_map<Index, Real> pBC;
    std::unordered_map<Index, Real> phinBC;
    std::unordered_map<Index, Real> phipBC;

    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& contact = mesh.getContact(c);
        auto it = contactBiases.find(contact.name);
        if (it == contactBiases.end()) continue;
        const double Vbias = it->second;

        for (Index nid : contact.node_ids) {
            const double ni_node  = ni_v[nid];
            const double Ndop     = doping.netDoping(nid);
            const double n_eq_val = nEq(Ndop, ni_node);
            const double p_eq_val = (ni_node > 0.0)
                                        ? ni_node * ni_node / n_eq_val
                                        : 0.0;
            // built-in potential (0 if ni == 0 or n_eq == 0)
            double psi_bi = 0.0;
            if (ni_node > 0.0 && n_eq_val > 0.0)
                psi_bi = Vt * std::log(n_eq_val / ni_node);

            psiBC [nid] = Vbias + psi_bi;
            nBC   [nid] = n_eq_val;
            pBC   [nid] = p_eq_val;
            phinBC[nid] = Vbias;
            phipBC[nid] = Vbias;
        }
    }

    // ------------------------------------------------------------------
    // Initial guess: solve linear Poisson (no carriers) for ψ
    // ------------------------------------------------------------------
    MobilityModelConfig mobilityConfig = mobilityModelConfig(cfg.mobility);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg.recombination, cfg.taun, cfg.taup);
    DDAssembler assembler(mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    VectorXd psi_init  = VectorXd::Zero(static_cast<int>(N));
    VectorXd n_init    = VectorXd::Zero(static_cast<int>(N));
    VectorXd p_init    = VectorXd::Zero(static_cast<int>(N));
    VectorXd phin      = VectorXd::Zero(static_cast<int>(N));
    VectorXd phip      = VectorXd::Zero(static_cast<int>(N));

    const bool useInitialGuess = initialGuess != nullptr &&
        initialGuess->psi.size() == static_cast<int>(N) &&
        initialGuess->phin.size() == static_cast<int>(N) &&
        initialGuess->phip.size() == static_cast<int>(N) &&
        initialGuess->n.size() == static_cast<int>(N) &&
        initialGuess->p.size() == static_cast<int>(N);

    if (useInitialGuess) {
        psi_init = initialGuess->psi;
        phin = initialGuess->phin;
        phip = initialGuess->phip;
        n_init = initialGuess->n;
        p_init = initialGuess->p;
    } else {
        // Override contact nodes before solving an initial Poisson problem.
        for (const auto& [nid, val] : psiBC)
            psi_init(static_cast<int>(nid)) = val;

        // Solve initial Poisson (no free carriers)
        {
            VectorXd n_zero = VectorXd::Zero(static_cast<int>(N));
            VectorXd p_zero = VectorXd::Zero(static_cast<int>(N));
            assembler.assemblePoissonWithCarriers(n_zero, p_zero, psi_init);
            assembler.applyDirichlet(psiBC);
            LinearSolver ls;
            psi_init = ls.solve(assembler.matrix(), assembler.rhs());
        }

        for (const auto& [nid, val] : phinBC) phin(static_cast<int>(nid)) = val;
        for (const auto& [nid, val] : phipBC) phip(static_cast<int>(nid)) = val;

        for (Index i = 0; i < N; ++i) {
            const int    ii     = static_cast<int>(i);
            const double ni_i   = ni_v[i];
            n_init(ii) = electronDensity(ni_i, psi_init(ii), phin(ii), Vt);
            p_init(ii) = holeDensity    (ni_i, psi_init(ii), phip(ii), Vt);
        }
    }

    // Enforce the new bias point's Ohmic contact values even when the
    // previous bias point is used as the initial guess.
    for (const auto& [nid, val] : psiBC)  psi_init(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : phinBC) phin(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : phipBC) phip(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : nBC)    n_init(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : pBC)    p_init(static_cast<int>(nid)) = val;

    // ------------------------------------------------------------------
    // Working copies
    // ------------------------------------------------------------------
    VectorXd psi = psi_init;
    VectorXd n   = n_init;
    VectorXd p   = p_init;

    LinearSolver ls;
    int iters = 0;
    bool converged = false;

    // ------------------------------------------------------------------
    // Gummel iteration
    // ------------------------------------------------------------------
    for (int iter = 0; iter < cfg.maxIter; ++iter) {
        ++iters;

        // ---- a. Solve Poisson with current carriers ----
        assembler.assemblePoissonWithCarriers(n, p, psi);
        assembler.applyDirichlet(psiBC);
        VectorXd psi_new = ls.solve(assembler.matrix(), assembler.rhs());

        // Damped update: compute full correction then apply with damping factor
        const VectorXd dpsi         = psi_new - psi;
        const VectorXd dpsi_applied = cfg.dampingPsi * dpsi;
        psi += dpsi_applied;

        // ---- b. Solve electron continuity for n ----
        assembler.assembleElectronContinuity(psi, n, p);
        assembler.applyDirichlet(nBC);
        VectorXd n_new = ls.solve(assembler.matrix(), assembler.rhs());

        // Enforce positivity (guard against small negative artefacts)
        for (int ii = 0; ii < static_cast<int>(N); ++ii)
            if (n_new(ii) < 0.0) { n_new(ii) = 0.0; }

        // ---- c. Solve hole continuity for p ----
        assembler.assembleHoleContinuity(psi, n_new, p);
        assembler.applyDirichlet(pBC);
        VectorXd p_new = ls.solve(assembler.matrix(), assembler.rhs());

        for (int ii = 0; ii < static_cast<int>(N); ++ii)
            if (p_new(ii) < 0.0) { p_new(ii) = 0.0; }

        // ---- d. Update quasi-Fermi potentials ----
        for (Index i = 0; i < N; ++i) {
            const int    ii   = static_cast<int>(i);
            const double ni_i = ni_v[i];
            if (ni_i > 0.0 && n_new(ii) > 0.0)
                phin(ii) = psi(ii) - Vt * std::log(n_new(ii) / ni_i);
            if (ni_i > 0.0 && p_new(ii) > 0.0)
                phip(ii) = psi(ii) + Vt * std::log(p_new(ii) / ni_i);
        }
        // Restore Dirichlet quasi-Fermi values at contacts
        for (const auto& [nid, val] : phinBC) phin(static_cast<int>(nid)) = val;
        for (const auto& [nid, val] : phipBC) phip(static_cast<int>(nid)) = val;

        // ---- e. Convergence check ----
        // Base residuals on the actual applied updates so that damping < 1
        // does not prevent detection of convergence.
        const double dpsi_norm = dpsi_applied.norm();
        const double psi_norm  = psi.norm();
        const double rel_err   = (psi_norm > 1.0e-30)
                                     ? dpsi_norm / psi_norm
                                     : dpsi_norm;

        const VectorXd dn = n_new - n;
        const double dn_norm = dn.norm();
        const double n_norm  = n.norm();
        const double rel_n   = (n_norm > 1.0e-30) ? dn_norm / n_norm : dn_norm;

        const VectorXd dp = p_new - p;
        const double dp_norm = dp.norm();
        const double p_norm  = p.norm();
        const double rel_p   = (p_norm > 1.0e-30) ? dp_norm / p_norm : dp_norm;

        n = n_new;
        p = p_new;

        if (rel_err < cfg.reltol && rel_n < cfg.reltol && rel_p < cfg.reltol) {
            converged = true;
            break;
        }
    }

    DDSolution sol;
    sol.psi   = psi;
    sol.phin  = phin;
    sol.phip  = phip;
    sol.n     = n;
    sol.p     = p;
    sol.iters = iters;
    sol.converged = converged;
    return sol;
}

} // anonymous namespace

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, cfg, nullptr);
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, cfg, &initialGuess);
}

// ---------------------------------------------------------------------------
// writeDDSolutionVTK
// ---------------------------------------------------------------------------

void writeDDSolutionVTK(const std::string& filename,
                        const DeviceMesh&  mesh,
                        const DopingModel& doping,
                        const DDSolution&  sol)
{
    const Index N = mesh.numNodes();

    VTKWriter writer(filename, mesh);
    writer.write();

    auto toVec = [&](const VectorXd& v) {
        std::vector<Real> out(N);
        for (Index i = 0; i < N; ++i)
            out[i] = v(static_cast<int>(i));
        return out;
    };

    writer.addNodeScalar("Potential",            toVec(sol.psi));
    writer.addNodeScalar("ElectronQuasiFermi",   toVec(sol.phin));
    writer.addNodeScalar("HoleQuasiFermi",       toVec(sol.phip));
    writer.addNodeScalar("Electrons",            toVec(sol.n));
    writer.addNodeScalar("Holes",                toVec(sol.p));

    std::vector<Real> netDop(N);
    for (Index i = 0; i < N; ++i)
        netDop[i] = doping.netDoping(i);
    writer.addNodeScalar("NetDoping", netDop);
}

} // namespace vela
