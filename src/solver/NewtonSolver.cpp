#include "vela/solver/NewtonSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace vela {
namespace {

inline double thermalVoltage(double T)
{
    if (T <= 0.0)
        throw std::invalid_argument("thermalVoltage: temperature_K must be positive.");
    return constants::kb * T / constants::q;
}

inline double nEq(double Ndop, double ni)
{
    const double half = 0.5 * Ndop;
    return half + std::sqrt(half * half + ni * ni);
}

void validateResidualWeights(
    Real psi,
    Real phin,
    Real phip,
    const char* context)
{
    const bool allFinite = std::isfinite(psi) && std::isfinite(phin) && std::isfinite(phip);
    const bool allNonnegative = psi >= 0.0 && phin >= 0.0 && phip >= 0.0;
    const bool anyEnabled = psi > 0.0 || phin > 0.0 || phip > 0.0;
    if (!allFinite || !allNonnegative || !anyEnabled) {
        throw std::invalid_argument(
            std::string(context)
            + ": residual_weights values must be finite, nonnegative, "
              "and leave at least one block enabled.");
    }
}

ResidualBlockWeights residualWeightsFromConfig(const NewtonConfig& cfg)
{
    return {cfg.residualWeightPsi, cfg.residualWeightPhin, cfg.residualWeightPhip};
}

ResidualBlockNormValue residualScalesFromConfig(
    const NewtonConfig& cfg,
    const ResidualBlockNormValue& initialBlocks)
{
    ResidualBlockNormValue scales;
    scales.psi = std::max(initialBlocks.psi, 1.0);
    scales.phin = std::max(initialBlocks.phin, 1.0);
    scales.phip = std::max(initialBlocks.phip, 1.0);
    if (cfg.residualScalePsi > 0.0)
        scales.psi = cfg.residualScalePsi;
    if (cfg.residualScalePhin > 0.0)
        scales.phin = cfg.residualScalePhin;
    if (cfg.residualScalePhip > 0.0)
        scales.phip = cfg.residualScalePhip;
    scales.combined = std::sqrt(scales.psi * scales.psi
        + scales.phin * scales.phin
        + scales.phip * scales.phip);
    return scales;
}

Real maxRelativePermittivityAcrossRegions(const DeviceMesh& mesh,
                                          const MaterialDatabase& matdb,
                                          Real temperature_K)
{
    Real maxEpsr = 0.0;
    for (const Region& region : mesh.regions()) {
        const Material& material = matdb.getMaterial(region.material, temperature_K);
        maxEpsr = std::max(maxEpsr, material.eps_r);
    }
    return std::max(maxEpsr, 1.0);
}

Real maxIntrinsicDensityAcrossRegions(const DeviceMesh& mesh,
                                      const MaterialDatabase& matdb,
                                      Real temperature_K)
{
    Real maxNi = 1.0;
    for (const Region& region : mesh.regions()) {
        const Material& material = matdb.getMaterial(region.material, temperature_K);
        maxNi = std::max(maxNi, material.ni);
    }
    return maxNi;
}

} // namespace


NewtonConfig newtonConfigFromJson(const nlohmann::json& json, UnitScalingConfig scaling)
{
    NewtonConfig cfg;
    cfg.inputScaling = scaling;
    cfg.maxIter = json.value("max_iter", cfg.maxIter);
    cfg.reltol = json.value("reltol", cfg.reltol);
    cfg.abstol = json.value("abstol", cfg.abstol);
    cfg.temperature_K = json.value("temperature_K", cfg.temperature_K);
    cfg.dampingFactor = json.value("damping_factor", cfg.dampingFactor);
    cfg.lineSearch = json.value("line_search", cfg.lineSearch);
    cfg.verbose = json.value("verbose", cfg.verbose);
    cfg.warmStart = json.value("warm_start", cfg.warmStart);
    cfg.diagnostics = json.value("diagnostics", cfg.diagnostics);
    cfg.diagnostics = json.value("diagnostic_history", cfg.diagnostics);
    cfg.finiteDifferenceStep = json.value("finite_difference_step", cfg.finiteDifferenceStep);
    cfg.jacobian = json.value("jacobian", cfg.jacobian);
    cfg.residualNorm = json.value("residual_norm", cfg.residualNorm);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    if (json.contains("mobility"))
        cfg.mobility = mobilityModelConfigFromJson(json.at("mobility"), scaling);
    if (json.contains("bandgap_narrowing")) {
        const auto& value = json.at("bandgap_narrowing");
        if (value.is_string()) {
            cfg.bandgapNarrowing.model = value.get<std::string>();
        } else if (value.is_object()) {
            cfg.bandgapNarrowing.model = value.value("model", cfg.bandgapNarrowing.model);
            if (value.contains("reference_doping_m3")) {
                cfg.bandgapNarrowing.referenceDoping = scaling.concentrationToSI(
                    value.at("reference_doping_m3").get<Real>());
            }
            cfg.bandgapNarrowing.coefficient = value.value(
                "coefficient_eV", cfg.bandgapNarrowing.coefficient);
            cfg.bandgapNarrowing.smoothing = value.value(
                "smoothing", cfg.bandgapNarrowing.smoothing);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: bandgap_narrowing must be a string or object.");
        }
    }
    if (json.contains("residual_weights")) {
        const auto& weights = json.at("residual_weights");
        cfg.residualWeightPsi = weights.value("psi", cfg.residualWeightPsi);
        cfg.residualWeightPhin = weights.value("phin", cfg.residualWeightPhin);
        cfg.residualWeightPhip = weights.value("phip", cfg.residualWeightPhip);
    }
    if (json.contains("residual_scales")) {
        const auto& scales = json.at("residual_scales");
        cfg.residualScalePsi = scales.value("psi", cfg.residualScalePsi);
        cfg.residualScalePhin = scales.value("phin", cfg.residualScalePhin);
        cfg.residualScalePhip = scales.value("phip", cfg.residualScalePhip);
    }
    if (json.contains("recombination")) {
        const auto& value = json.at("recombination");
        if (value.is_array())
            cfg.recombination = value.get<std::vector<std::string>>();
        else if (value.is_string())
            cfg.recombination = {value.get<std::string>()};
        else
            throw std::invalid_argument(
                "newtonConfigFromJson: recombination must be a string or string array.");
    }

    if (json.contains("impact_ionization")) {
        const auto& value = json.at("impact_ionization");
        if (value.is_string()) {
            cfg.impactIonization.model = value.get<std::string>();
        } else if (value.is_object()) {
            cfg.impactIonization.model = value.value("model", cfg.impactIonization.model);
            if (value.contains("electron_A_m_inv")) {
                cfg.impactIonization.electronA = scaling.inverseLengthToSI(
                    value.at("electron_A_m_inv").get<Real>());
            }
            if (value.contains("electron_B_V_m")) {
                cfg.impactIonization.electronB = scaling.electricFieldToSI(
                    value.at("electron_B_V_m").get<Real>());
            }
            if (value.contains("hole_A_m_inv")) {
                cfg.impactIonization.holeA = scaling.inverseLengthToSI(
                    value.at("hole_A_m_inv").get<Real>());
            }
            if (value.contains("hole_B_V_m")) {
                cfg.impactIonization.holeB = scaling.electricFieldToSI(
                    value.at("hole_B_V_m").get<Real>());
            }
            cfg.impactIonization.carrierVelocity = value.value(
                "carrier_velocity_m_s", cfg.impactIonization.carrierVelocity);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: impact_ionization must be a string or object.");
        }
    }

    if (cfg.jacobian != "analytic" && cfg.jacobian != "finite_difference")
        throw std::invalid_argument(
            "newtonConfigFromJson: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg.inputScaling.isUnitScaling())
        cfg.jacobian = "finite_difference";
    if (cfg.residualNorm != "block" && cfg.residualNorm != "l2")
        throw std::invalid_argument(
            "newtonConfigFromJson: residual_norm must be 'block' or 'l2'.");
    validateResidualWeights(
        cfg.residualWeightPsi,
        cfg.residualWeightPhin,
        cfg.residualWeightPhip,
        "newtonConfigFromJson");
    if (cfg.temperature_K <= 0.0)
        throw std::invalid_argument("newtonConfigFromJson: temperature_K must be positive.");

    return cfg;
}

NewtonSolver::NewtonSolver(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    NewtonConfig cfg,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , contactBiases_(contactBiases)
    , cfg_(cfg)
    , fixedCharges_(std::move(fixedCharges))
    , sheetCharges_(std::move(sheetCharges))
{
    if (cfg_.jacobian != "analytic" && cfg_.jacobian != "finite_difference")
        throw std::invalid_argument(
            "NewtonSolver: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg_.residualNorm != "block" && cfg_.residualNorm != "l2")
        throw std::invalid_argument(
            "NewtonSolver: residual_norm must be 'block' or 'l2'.");
    validateResidualWeights(
        cfg_.residualWeightPsi,
        cfg_.residualWeightPhin,
        cfg_.residualWeightPhip,
        "NewtonSolver");
}

DDScalingSpec NewtonSolver::buildScalingSpec() const
{
    DDScalingSpec scaling;
    if (!cfg_.inputScaling.isUnitScaling())
        return scaling;

    const Real epsRef = constants::eps0 *
        maxRelativePermittivityAcrossRegions(mesh_, matdb_, cfg_.temperature_K);
    const Real niFloor =
        maxIntrinsicDensityAcrossRegions(mesh_, matdb_, cfg_.temperature_K);
    const UnitScalingSystem::AutoInputs autoInputs =
        UnitScalingSystem::autoInputsFrom(mesh_, doping_, matdb_, niFloor);
    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        cfg_.temperature_K, epsRef, autoInputs, cfg_.unitScalingRefs);

    scaling.enabled = true;
    scaling.V0 = sc.V0();
    scaling.C0 = sc.C0();
    scaling.mu0 = sc.mu0();
    scaling.D0 = sc.D0();
    scaling.L0 = sc.L0();
    scaling.permittivityReference_F_per_m = epsRef;
    return scaling;
}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler) const
{
    CoupledDDBoundaryConditions bcs;
    const auto& ni = assembler.intrinsicDensity();
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;

    for (Index c = 0; c < mesh_.numContacts(); ++c) {
        const Contact& contact = mesh_.getContact(c);
        auto it = contactBiases_.find(contact.name);
        if (it == contactBiases_.end()) continue;

        const double Vbias = it->second;
        for (Index nid : contact.node_ids) {
            const double niNode = ni[nid];
            const double neq = nEq(doping_.netDoping(nid), niNode);
            double psiBuiltIn = 0.0;
            if (niNode > 0.0 && neq > 0.0)
                psiBuiltIn = Vt * std::log(neq / niNode);

            bcs.psi[nid] = (Vbias + psiBuiltIn) / potentialScale;
            bcs.phin[nid] = Vbias / potentialScale;
            bcs.phip[nid] = Vbias / potentialScale;
        }
    }
    return bcs;
}

DDSolution NewtonSolver::buildInitialGuess(
    const CoupledDDAssembler&, const CoupledDDBoundaryConditions&) const
{
    GummelConfig gcfg;
    gcfg.maxIter = 1;
    gcfg.reltol = 0.0;
    gcfg.temperature_K = cfg_.temperature_K;
    gcfg.dampingPsi = 0.5;
    gcfg.taun = cfg_.taun;
    gcfg.taup = cfg_.taup;
    gcfg.inputScaling = cfg_.inputScaling;
    gcfg.unitScalingRefs = cfg_.unitScalingRefs;
    gcfg.mobility = cfg_.mobility;
    gcfg.recombination = cfg_.recombination;
    gcfg.bandgapNarrowing = cfg_.bandgapNarrowing;
    gcfg.impactIonization = cfg_.impactIonization;
    DDSolution sol = runGummel(mesh_, matdb_, doping_, contactBiases_, ContactSpecsMap{}, gcfg, fixedCharges_, sheetCharges_);

    // The default cold-start path removes tiny quasi-Fermi numerical noise left
    // by the one-step Gummel initializer.  A caller can opt into warm_start when
    // the supplied/constructed quasi-Fermi potentials should be used as-is.
    if (!cfg_.warmStart) {
        const int N = static_cast<int>(mesh_.numNodes());
        sol.phin = VectorXd::Zero(N);
        sol.phip = VectorXd::Zero(N);
    }
    return sol;
}

DDSolution NewtonSolver::makeSolution(const CoupledDDAssembler& assembler,
                                      const VectorXd& x,
                                      int iters) const
{
    CoupledDDState state = assembler.unpack(x);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    DDSolution sol;
    sol.psi = state.psi * potentialScale;
    sol.phin = state.phin * potentialScale;
    sol.phip = state.phip * potentialScale;
    sol.n = assembler.electronDensity(x);
    sol.p = assembler.holeDensity(x);
    sol.iters = iters;
    return sol;
}

NewtonResult NewtonSolver::solve() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    const DDScalingSpec scaling = buildScalingSpec();
    CoupledDDAssembler assembler(
        mesh_,
        matdb_,
        doping_,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg_.bandgapNarrowing,
        cfg_.impactIonization,
        fixedCharges_,
        sheetCharges_,
        scaling);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    return solve(buildInitialGuess(assembler, bcs));
}

NewtonResult NewtonSolver::solve(const DDSolution& initial) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    const DDScalingSpec scaling = buildScalingSpec();
    CoupledDDAssembler assembler(
        mesh_,
        matdb_,
        doping_,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg_.bandgapNarrowing,
        cfg_.impactIonization,
        fixedCharges_,
        sheetCharges_,
        scaling);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);

    // By default Newton uses a conservative cold start for quasi-Fermi
    // potentials: interior phin/phip are reset to equilibrium values because
    // tiny external-initializer noise can be strongly amplified by the balanced
    // Scharfetter-Gummel flux.  Set warm_start=true to preserve the supplied
    // quasi-Fermi potentials, which is useful for continuation runs where the
    // previous bias point is already a high-quality initial guess.
    VectorXd phinInit = initial.phin;
    VectorXd phipInit = initial.phip;
    const int N = static_cast<int>(mesh_.numNodes());
    if (!cfg_.warmStart) {
        for (int i = 0; i < N; ++i) {
            const Index nid = static_cast<Index>(i);
            if (bcs.phin.find(nid) == bcs.phin.end()) {
                phinInit(i) = 0.0;
                phipInit(i) = 0.0;
            }
        }
    }

    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    VectorXd x = assembler.pack({
        initial.psi / potentialScale,
        phinInit / potentialScale,
        phipInit / potentialScale});
    VectorXd r = assembler.residual(x, bcs);
    const ResidualBlockNormValue initialBlocks =
        ResidualNorm::computeBlocks(r, mesh_.numNodes());
    const ResidualBlockNormValue residualScales =
        residualScalesFromConfig(cfg_, initialBlocks);
    const ResidualBlockWeights residualWeights = residualWeightsFromConfig(cfg_);
    const auto residualNormFn = [&](const VectorXd& residual) {
        if (cfg_.residualNorm == "l2")
            return residual.norm();
        return ResidualNorm::normalizedBlockL2(
            ResidualNorm::computeBlocks(residual, mesh_.numNodes()),
            residualScales,
            residualWeights);
    };
    const Real initialNorm = residualNormFn(r);

    NewtonResult result;
    result.solution = initial;
    result.initialResidualNorm = initialNorm;
    result.finalResidualNorm = initialNorm;

    if (cfg_.verbose) {
        std::cout << "Newton iter 0 residual=" << initialNorm
                  << " step=0 damping=0\n";
    }

    if (initialNorm <= cfg_.abstol) {
        result.converged = true;
        result.solution = makeSolution(assembler, x, 0);
        return result;
    }

    LinearSolver linearSolver;
    LineSearchConfig lscfg;
    lscfg.enabled = cfg_.lineSearch;
    lscfg.initialDamping = cfg_.dampingFactor;
    lscfg.recordHistory = cfg_.diagnostics;
    BacktrackingLineSearch lineSearch(lscfg);

    VectorXd acceptedX = x;
    VectorXd acceptedR = r;
    int acceptedIters = 0;

    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
            ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
            : assembler.assembleJacobian(x, bcs);
        VectorXd step;
        try {
            step = linearSolver.solve(J, -r);
        } catch (const std::runtime_error&) {
            result.finalResidualNorm = residualNormFn(acceptedR);
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            if (cfg_.verbose) {
                std::cerr << "Newton failed at iter " << iter
                          << ": residual=" << residualNormFn(r)
                          << " damping=0 step=0 (linear solve failed)\n";
            }
            return result;
        }
        const Real stepNorm = step.norm();

        auto ls = lineSearch.search(
            x, step, r,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            },
            residualNormFn);

        if (!ls.accepted) {
            result.finalResidualNorm = residualNormFn(acceptedR);
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            if (cfg_.verbose) {
                std::cerr << "Newton failed at iter " << iter
                          << ": residual=" << ls.residualNorm
                          << " damping=" << ls.damping
                          << " step=" << stepNorm
                          << " (line search rejected step)\n";
            }
            return result;
        }

        x = ls.x;
        r = ls.residual;
        acceptedX = x;
        acceptedR = r;
        acceptedIters = iter;

        // Record the norm of the actually applied update (damped step) so that
        // per-iteration metrics are consistent with the accepted solution.
        const Real appliedStepNorm = ls.damping * stepNorm;
        const Real residualNorm = ls.residualNorm;
        NewtonIterationInfo info;
        info.iter = iter;
        info.residualNorm = residualNorm;
        info.stepNorm = appliedStepNorm;
        info.dampingFactor = ls.damping;
        info.relativeResidualNorm = ResidualNorm::relative(residualNorm, initialNorm);
        info.rawStepNorm = stepNorm;
        info.lineSearchAttempts = ls.attempts;
        info.lineSearchAccepted = ls.accepted;
        if (cfg_.diagnostics)
            info.lineSearchHistory = std::move(ls.history);
        result.history.push_back(std::move(info));
        if (cfg_.verbose) {
            std::cout << "Newton iter " << iter
                      << " residual=" << residualNorm
                      << " step=" << appliedStepNorm
                      << " damping=" << ls.damping << '\n';
        }

        const Real rel = result.history.back().relativeResidualNorm;
        if (residualNorm <= cfg_.abstol || rel <= cfg_.reltol) {
            result.converged = true;
            result.iters = iter;
            result.finalResidualNorm = residualNorm;
            result.solution = makeSolution(assembler, x, iter);
            return result;
        }
    }

    result.converged = false;
    result.iters = acceptedIters;
    result.finalResidualNorm = residualNormFn(acceptedR);
    result.solution = makeSolution(assembler, acceptedX, acceptedIters);
    if (cfg_.verbose) {
        std::cerr << "Newton failed after " << cfg_.maxIter
                  << " iterations: residual=" << result.finalResidualNorm
                  << " damping="
                  << (result.history.empty() ? 0.0 : result.history.back().dampingFactor)
                  << " step="
                  << (result.history.empty() ? 0.0 : result.history.back().stepNorm)
                  << '\n';
    }
    return result;
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg)
{
    return NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve();
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg)
{
    return NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve(initial);
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges)
{
    return NewtonSolver(
        mesh,
        matdb,
        doping,
        contactBiases,
        cfg,
        std::move(fixedCharges),
        std::move(sheetCharges)).solve();
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges)
{
    return NewtonSolver(
        mesh,
        matdb,
        doping,
        contactBiases,
        cfg,
        std::move(fixedCharges),
        std::move(sheetCharges)).solve(initial);
}

} // namespace vela
