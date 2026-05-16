#include "vela/solver/GummelSolver.h"
#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/DDAssembler.h"
#include "vela/solver/LinearSolver.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/io/VTKWriter.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <optional>

namespace vela {



// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

/// Thermal voltage at temperature T [K].
inline double thermalVoltage(double T)
{
    if (T <= 0.0)
        throw std::invalid_argument("thermalVoltage: temperature_K must be positive.");
    return constants::kb * T / constants::q;
}

/// Charge-neutral equilibrium electron concentration from net doping and ni.
/// Solves: n - p = Ndop,  n*p = ni^2  ->  n = Ndop/2 + hypot(Ndop/2, ni)
inline double nEq(double Ndop, double ni)
{
    const double half = 0.5 * Ndop;
    const double root = std::hypot(half, ni);

    if (Ndop >= 0.0)
        return half + root;

    // std::hypot avoids overflow in the radicand. For p-type contacts with
    // |Ndop| >> ni, half + root still suffers severe cancellation, so use
    // the algebraically equivalent minority-carrier form
    // n = ni^2 / p, with p = root - half, to preserve high-doping accuracy.
    const double p_eq = root - half;
    return (p_eq > 0.0) ? (ni * ni / p_eq) : 0.0;
}

/// Compute per-node ni from the material database.
std::vector<double> buildNiVector(const DeviceMesh&       mesh,
                                  const MaterialDatabase& matdb,
                                  Real                    cfgTemperature_K)
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
            ni_mat = matdb.getMaterial(region.material, cfgTemperature_K).ni;
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
    cfg.abstol = json.value("abstol", cfg.abstol);
    cfg.temperature_K = json.value("temperature_K", cfg.temperature_K);
    cfg.dampingPsi = json.value("damping_psi", cfg.dampingPsi);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    cfg.mobility = json.value("mobility", cfg.mobility);
    if (json.contains("bandgap_narrowing")) {
        const auto& value = json.at("bandgap_narrowing");
        if (value.is_string()) {
            cfg.bandgapNarrowing.model = value.get<std::string>();
        } else if (value.is_object()) {
            cfg.bandgapNarrowing.model = value.value("model", cfg.bandgapNarrowing.model);
            cfg.bandgapNarrowing.referenceDoping = value.value(
                "reference_doping_m3", cfg.bandgapNarrowing.referenceDoping);
            cfg.bandgapNarrowing.coefficient = value.value(
                "coefficient_eV", cfg.bandgapNarrowing.coefficient);
            cfg.bandgapNarrowing.smoothing = value.value(
                "smoothing", cfg.bandgapNarrowing.smoothing);
        } else {
            throw std::invalid_argument(
                "gummelConfigFromJson: bandgap_narrowing must be a string or object.");
        }
    }

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


    if (json.contains("impact_ionization")) {
        const auto& value = json.at("impact_ionization");
        if (value.is_string()) {
            cfg.impactIonization.model = value.get<std::string>();
        } else if (value.is_object()) {
            cfg.impactIonization.model = value.value("model", cfg.impactIonization.model);
            cfg.impactIonization.electronA = value.value(
                "electron_A_m_inv", cfg.impactIonization.electronA);
            cfg.impactIonization.electronB = value.value(
                "electron_B_V_m", cfg.impactIonization.electronB);
            cfg.impactIonization.holeA = value.value(
                "hole_A_m_inv", cfg.impactIonization.holeA);
            cfg.impactIonization.holeB = value.value(
                "hole_B_V_m", cfg.impactIonization.holeB);
            cfg.impactIonization.carrierVelocity = value.value(
                "carrier_velocity_m_s", cfg.impactIonization.carrierVelocity);
        } else {
            throw std::invalid_argument(
                "gummelConfigFromJson: impact_ionization must be a string or object.");
        }
    }

    if (cfg.temperature_K <= 0.0)
        throw std::invalid_argument("gummelConfigFromJson: temperature_K must be positive.");

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
                         const ContactSpecsMap&                       contactSpecs,
                         const GummelConfig&                          cfg,
                         const DDSolution*                           initialGuess)
{
    const Index  N   = mesh.numNodes();
    const double Vt  = thermalVoltage(cfg.temperature_K);

    // Per-node effective ni, including optional bandgap narrowing.
    std::vector<double> ni_v = buildNiVector(mesh, matdb, cfg.temperature_K);
    const auto bgn = makeBandgapNarrowingModel(cfg.bandgapNarrowing);
    for (Index i = 0; i < N; ++i) {
        const double deltaEg = bgn->deltaEg(doping.netDoping(i), 0.0, 0.0);
        ni_v[i] = effectiveIntrinsicDensity(ni_v[i], Vt, deltaEg);
    }

    // Look up a contact-region material for Schottky barrier helpers.
    auto contactMaterial = [&](const Contact& contact) -> std::optional<Material> {
        if (contact.region_id < mesh.numRegions()) {
            const Region& region = mesh.getRegion(contact.region_id);
            if (matdb.hasMaterial(region.material))
                return matdb.getMaterial(region.material, cfg.temperature_K);
        }
        // Fallback: pick any cell that contains a contact node.
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            for (Index nid : cell.node_ids) {
                for (Index cnid : contact.node_ids) {
                    if (nid == cnid) {
                        const Region& region = mesh.getRegion(cell.region_id);
                        if (matdb.hasMaterial(region.material))
                            return matdb.getMaterial(region.material, cfg.temperature_K);
                    }
                }
            }
        }
        return std::nullopt;
    };

    // ------------------------------------------------------------------
    // Build contact Dirichlet BCs
    // ------------------------------------------------------------------
    //   Ohmic (default):
    //     psi_contact  = V_bias + Vt * ln(n_eq / ni)
    //     n_contact    = n_eq
    //     p_contact    = ni^2 / n_eq
    //     phin = phip  = V_bias
    //   Schottky (prototype Dirichlet barrier):
    //     ContactState from computeSchottkyContactState()
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

        auto specIt = contactSpecs.find(contact.name);
        const ContactBoundarySpec* spec =
            (specIt != contactSpecs.end()) ? &specIt->second : nullptr;
        const bool isSchottky =
            spec != nullptr && spec->type == ContactType::Schottky;

        if (isSchottky) {
            ContactBoundarySpec effSpec = *spec;
            effSpec.bias = Vbias; // honour the swept bias point

            const auto materialOpt = contactMaterial(contact);
            const Real bandgap_eV = materialOpt && materialOpt->bandgap_eV
                ? *materialOpt->bandgap_eV
                : std::numeric_limits<Real>::quiet_NaN();
            const Real affinity_eV = materialOpt && materialOpt->electron_affinity_eV
                ? *materialOpt->electron_affinity_eV
                : std::numeric_limits<Real>::quiet_NaN();
            const Real NcVal = materialOpt && materialOpt->Nc_m3
                ? *materialOpt->Nc_m3 : 0.0;
            const Real NvVal = materialOpt && materialOpt->Nv_m3
                ? *materialOpt->Nv_m3 : 0.0;

            for (Index nid : contact.node_ids) {
                const double ni_node = ni_v[nid];
                const double Ndop = doping.netDoping(nid);
                const ContactState state = computeSchottkyContactState(
                    effSpec, ni_node, NcVal, NvVal,
                    bandgap_eV, affinity_eV, Ndop, cfg.temperature_K);
                psiBC [nid] = state.psi;
                nBC   [nid] = state.n;
                pBC   [nid] = state.p;
                phinBC[nid] = state.phin;
                phipBC[nid] = state.phip;
            }
            continue;
        }

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
    // Initial guess: solve linear Poisson (no carriers) for psi
    // ------------------------------------------------------------------
    MobilityModelConfig mobilityConfig = mobilityModelConfig(cfg.mobility);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg.recombination, cfg.taun, cfg.taup);
    DDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg.bandgapNarrowing,
        cfg.impactIonization);

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

        const bool relativeConverged =
            rel_err < cfg.reltol && rel_n < cfg.reltol && rel_p < cfg.reltol;
        const bool absoluteConverged = cfg.abstol > 0.0 &&
            dpsi_norm <= cfg.abstol && dn_norm <= cfg.abstol && dp_norm <= cfg.abstol;
        if (relativeConverged || absoluteConverged) {
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
    ContactSpecsMap emptySpecs;
    return runGummelImpl(mesh, matdb, doping, contactBiases, emptySpecs, cfg, nullptr);
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess)
{
    ContactSpecsMap emptySpecs;
    return runGummelImpl(mesh, matdb, doping, contactBiases, emptySpecs, cfg, &initialGuess);
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, contactSpecs, cfg, nullptr);
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, contactSpecs, cfg, &initialGuess);
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
