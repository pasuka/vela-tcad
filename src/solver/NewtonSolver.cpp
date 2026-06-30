#include "vela/solver/NewtonSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace vela {
namespace {

void parseImpactIonizationDrivingForceInterpolation(
    const nlohmann::json& value,
    const UnitScalingConfig& scaling,
    ImpactIonizationModelConfig& config,
    const char* context)
{
    if (!value.contains("driving_force_interpolation"))
        return;

    const auto& interpolation = value.at("driving_force_interpolation");
    if (interpolation.is_string()) {
        config.drivingForceInterpolation = interpolation.get<std::string>();
        return;
    }
    if (!interpolation.is_object()) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation must be a string or object.");
    }

    config.drivingForceInterpolation = interpolation.value(
        "mode", config.drivingForceInterpolation);
    if (interpolation.contains("electron_ref_density_m3")) {
        config.electronDrivingForceRefDensity = scaling.concentrationToSI(
            interpolation.at("electron_ref_density_m3").get<Real>());
    }
    if (interpolation.contains("hole_ref_density_m3")) {
        config.holeDrivingForceRefDensity = scaling.concentrationToSI(
            interpolation.at("hole_ref_density_m3").get<Real>());
    }
}

std::string normalizeToken(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        if (ch == '-' || std::isspace(ch))
            return static_cast<char>('_');
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string canonicalContactBoundaryReconstruction(const std::string& raw)
{
    const std::string mode = normalizeToken(raw);
    if (mode == "dominant_signed_contact_mean" || mode == "dominant_signed_region")
        return "dominant_signed_contact_mean";
    if (mode == "legacy_node_local" || mode == "node_local")
        return "legacy_node_local";
    throw std::invalid_argument(
        "newtonConfigFromJson: contact_boundary_reconstruction must be "
        "'dominant_signed_contact_mean' or 'legacy_node_local'.");
}

std::string canonicalMinorityRelaxationContactSide(const std::string& raw)
{
    const std::string mode = normalizeToken(raw);
    if (mode == "p_contact_only" || mode == "p_contact" || mode == "p")
        return "p_contact_only";
    if (mode == "n_contact_only" || mode == "n_contact" || mode == "n")
        return "n_contact_only";
    if (mode == "both_contacts" || mode == "both" || mode == "all_contacts")
        return "both_contacts";
    throw std::invalid_argument(
        "newtonConfigFromJson: "
        "contact_boundary_minority_electron_relaxation_contact_side must be "
        "'p_contact_only', 'n_contact_only', or 'both_contacts'.");
}

Real clampMinorityRelaxationStrength(Real value)
{
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw std::invalid_argument(
            "newtonConfigFromJson: contact_boundary_minority_electron_relaxation_strength "
            "must be finite and lie in [0, 1].");
    }
    return value;
}

inline double thermalVoltage(double T)
{
    if (T <= 0.0)
        throw std::invalid_argument("thermalVoltage: temperature_K must be positive.");
    return constants::kb * T / constants::q;
}

inline double nEq(double Ndop, double ni)
{
    const double half = 0.5 * Ndop;
    const double root = std::hypot(half, ni);
    if (Ndop >= 0.0)
        return half + root;

    const double pEq = root - half;
    return (pEq > 0.0) ? (ni * ni / pEq) : 0.0;
}

Real ohmicContactNetDoping(const DeviceMesh& mesh,
                           const DopingModel& doping,
                           const Contact& contact,
                           Index nodeId,
                           bool dominantSignedContactMean)
{
    const Real local = doping.netDoping(nodeId);
    if (!dominantSignedContactMean)
        return local;
    if (contact.node_ids.empty())
        return local;

    Real sum = 0.0;
    for (Index nid : contact.node_ids)
        sum += doping.netDoping(nid);
    const Real mean = sum / static_cast<Real>(contact.node_ids.size());

    // Contact nodes in imported meshes can include compensated/tie points whose
    // node-owned signed doping flips relative to the contact side. For Ohmic
    // BC reconstruction, align those outliers with the contact-local average
    // polarity to avoid injecting an opposite-type built-in potential on one
    // endpoint of the same terminal.
    if (mean == 0.0 || local == 0.0)
        return mean != 0.0 ? mean : local;
    if ((local > 0.0 && mean > 0.0) || (local < 0.0 && mean < 0.0))
        return local;
    return mean;
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

NewtonBlockResidualInfo blockResidualInfo(const VectorXd& residual, Index nodeCount)
{
    const ResidualBlockNormValue blocks = ResidualNorm::computeBlocks(residual, nodeCount);
    return {blocks.psi, blocks.phin, blocks.phip, blocks.combined};
}

std::vector<int> jacobianAuditRows(const std::string& block,
                                   int nodeCount,
                                   const CoupledDDBoundaryConditions& bcs = {})
{
    std::vector<int> rows;
    if (block == "poisson") {
        for (int i = 0; i < nodeCount; ++i)
            rows.push_back(i);
    } else if (block == "transport" ||
               block == "srh_auger" ||
               block == "sg_avalanche") {
        for (int i = nodeCount; i < 3 * nodeCount; ++i)
            rows.push_back(i);
    } else if (block == "dirichlet_or_gauge") {
        for (const auto& [node, value] : bcs.psi) {
            (void)value;
            if (node < static_cast<Index>(nodeCount))
                rows.push_back(static_cast<int>(node));
        }
        for (const auto& [node, value] : bcs.phin) {
            (void)value;
            if (node < static_cast<Index>(nodeCount))
                rows.push_back(nodeCount + static_cast<int>(node));
        }
        for (const auto& [node, value] : bcs.phip) {
            (void)value;
            if (node < static_cast<Index>(nodeCount))
                rows.push_back(2 * nodeCount + static_cast<int>(node));
        }
    }
    return rows;
}

Real restrictedSparseNorm(const SparseMatrixd& matrix,
                          const std::vector<int>& rows)
{
    std::vector<char> rowMask(static_cast<std::size_t>(matrix.rows()), 0);
    for (int row : rows) {
        if (row >= 0 && row < matrix.rows())
            rowMask[static_cast<std::size_t>(row)] = 1;
    }

    Real sum = 0.0;
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator it(matrix, outer); it; ++it) {
            if (rowMask[static_cast<std::size_t>(it.row())]) {
                const Real value = it.value();
                sum += value * value;
            }
        }
    }
    return std::sqrt(sum);
}

NewtonJacobianBlockAuditRow jacobianAuditRow(
    const std::string& block,
    const SparseMatrixd& analytic,
    const SparseMatrixd& fd,
    const std::vector<int>& rows)
{
    const SparseMatrixd diff = analytic - fd;
    NewtonJacobianBlockAuditRow row;
    row.block = block;
    row.analyticNorm = restrictedSparseNorm(analytic, rows);
    row.fdNorm = restrictedSparseNorm(fd, rows);
    row.diffNorm = restrictedSparseNorm(diff, rows);
    row.relDiff = row.diffNorm / std::max<Real>(1.0, row.fdNorm);
    return row;
}

SparseMatrixd sparseBlock(const SparseMatrixd& matrix,
                          int rowStart,
                          int colStart,
                          int rows,
                          int cols)
{
    SparseMatrixd block = matrix.block(rowStart, colStart, rows, cols);
    block.makeCompressed();
    return block;
}

Real addCarrierRowRegularization(SparseMatrixd& matrix,
                                 int nodeCount,
                                 Real regularizationScale)
{
    if (regularizationScale <= 0.0)
        return 0.0;

    std::vector<Real> rowAbsSums(static_cast<std::size_t>(2 * nodeCount), 0.0);
    for (int col = 0; col < 3 * nodeCount; ++col) {
        for (SparseMatrixd::InnerIterator it(matrix, col); it; ++it) {
            const int row = static_cast<int>(it.row());
            if (row >= nodeCount && row < 3 * nodeCount) {
                rowAbsSums[static_cast<std::size_t>(row - nodeCount)] +=
                    std::abs(it.value());
            }
        }
    }

    Real diagonalNormSq = 0.0;
    for (int localRow = 0; localRow < 2 * nodeCount; ++localRow) {
        const int row = nodeCount + localRow;
        const Real diagonal = matrix.coeff(row, row);
        const Real sign = diagonal < 0.0 ? -1.0 : 1.0;
        const Real addition =
            sign * regularizationScale * rowAbsSums[static_cast<std::size_t>(localRow)];
        matrix.coeffRef(row, row) += addition;
        diagonalNormSq += addition * addition;
    }
    matrix.makeCompressed();
    return std::sqrt(diagonalNormSq);
}

bool clampQuasiFermiStep(VectorXd& step,
                         const DopingModel& doping,
                         Real globalLimit,
                         Real minorityLimit,
                         int nodeCount)
{
    if (globalLimit <= 0.0 && minorityLimit <= 0.0)
        return false;

    const auto resolveLimit = [&](bool isMinority) -> Real {
        if (isMinority && minorityLimit > 0.0)
            return globalLimit > 0.0 ? std::min(globalLimit, minorityLimit)
                                     : minorityLimit;
        return globalLimit;
    };
    const auto clampEntry = [&](int index, Real limit) -> bool {
        if (limit <= 0.0)
            return false;
        if (step(index) > limit) {
            step(index) = limit;
            return true;
        }
        if (step(index) < -limit) {
            step(index) = -limit;
            return true;
        }
        return false;
    };

    bool clippedQuasiFermi = false;
    for (int i = 0; i < nodeCount; ++i) {
        const Real net = doping.netDoping(i);
        const bool electronMinority = net < 0.0; // p-type node
        const bool holeMinority = net > 0.0;     // n-type node
        clippedQuasiFermi |= clampEntry(nodeCount + i, resolveLimit(electronMinority));
        clippedQuasiFermi |= clampEntry(2 * nodeCount + i, resolveLimit(holeMinority));
    }
    return clippedQuasiFermi;
}

bool applyConfiguredQuasiFermiStepCaps(VectorXd& step,
                                       const NewtonConfig& cfg,
                                       int nodeCount,
                                       Real potentialScale,
                                       const DopingModel& doping)
{
    const Real globalLimit = cfg.quasiFermiUpdateLimit_V > 0.0
        ? cfg.quasiFermiUpdateLimit_V / potentialScale
        : 0.0;
    const Real minorityLimit = cfg.quasiFermiUpdateLimitMinority_V > 0.0
        ? cfg.quasiFermiUpdateLimitMinority_V / potentialScale
        : 0.0;
    return clampQuasiFermiStep(step, doping, globalLimit, minorityLimit, nodeCount);
}

bool applyConfiguredStepCaps(VectorXd& step,
                             const NewtonConfig& cfg,
                             int nodeCount,
                             Real potentialScale,
                             const DopingModel& doping)
{
    if (cfg.maxUpdate > 0.0) {
        const Real maxAbsStep = step.cwiseAbs().maxCoeff();
        if (maxAbsStep > cfg.maxUpdate)
            step *= cfg.maxUpdate / maxAbsStep;
    }

    return applyConfiguredQuasiFermiStepCaps(
        step, cfg, nodeCount, potentialScale, doping);
}

void recorrectPoissonStepForClippedQuasiFermi(VectorXd& step,
                                              const SparseMatrixd& J,
                                              const VectorXd& residual,
                                              int nodeCount)
{
    if (nodeCount <= 0)
        return;

    const SparseMatrixd poissonBlock = sparseBlock(J, 0, 0, nodeCount, nodeCount);
    const SparseMatrixd poissonCarrierCoupling = sparseBlock(
        J, 0, nodeCount, nodeCount, 2 * nodeCount);
    const VectorXd qfStep = step.segment(nodeCount, 2 * nodeCount);
    const VectorXd rhs = -residual.segment(0, nodeCount) - poissonCarrierCoupling * qfStep;
    LinearSolver linearSolver;
    step.segment(0, nodeCount) = linearSolver.solve(poissonBlock, rhs);
}

void applyConfiguredStepCapsAndPoissonRecorrection(VectorXd& step,
                                                   const SparseMatrixd& J,
                                                   const VectorXd& residual,
                                                   const NewtonConfig& cfg,
                                                   int nodeCount,
                                                   Real potentialScale,
                                                   const DopingModel& doping)
{
    const bool clippedQuasiFermi = applyConfiguredStepCaps(
        step, cfg, nodeCount, potentialScale, doping);
    if (clippedQuasiFermi)
        recorrectPoissonStepForClippedQuasiFermi(step, J, residual, nodeCount);
}

NewtonCarrierDiagnostics carrierDiagnostics(const CoupledDDAssembler& assembler,
                                            const VectorXd& x)
{
    NewtonCarrierDiagnostics diagnostics;
    const VectorXd n = assembler.electronDensity(x);
    const VectorXd p = assembler.holeDensity(x);
    diagnostics.minElectronDensity = n.size() > 0
        ? std::numeric_limits<Real>::infinity()
        : 0.0;
    diagnostics.minHoleDensity = p.size() > 0
        ? std::numeric_limits<Real>::infinity()
        : 0.0;

    for (int i = 0; i < n.size(); ++i) {
        if (!std::isfinite(n(i)))
            ++diagnostics.nonfiniteElectronCount;
        else
            diagnostics.minElectronDensity = std::min(diagnostics.minElectronDensity, n(i));
        if (!std::isfinite(p(i)))
            ++diagnostics.nonfiniteHoleCount;
        else
            diagnostics.minHoleDensity = std::min(diagnostics.minHoleDensity, p(i));
        if (!(n(i) > 0.0))
            ++diagnostics.nonpositiveElectronCount;
        if (!(p(i) > 0.0))
            ++diagnostics.nonpositiveHoleCount;
    }

    if (!std::isfinite(diagnostics.minElectronDensity))
        diagnostics.minElectronDensity = 0.0;
    if (!std::isfinite(diagnostics.minHoleDensity))
        diagnostics.minHoleDensity = 0.0;
    diagnostics.positiveFinite = diagnostics.nonfiniteElectronCount == 0 &&
        diagnostics.nonfiniteHoleCount == 0 &&
        diagnostics.nonpositiveElectronCount == 0 &&
        diagnostics.nonpositiveHoleCount == 0;
    return diagnostics;
}

std::vector<NewtonTopResidualNode> topPoissonResidualNodes(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const CoupledDDAssembler& assembler,
    const VectorXd& residual,
    std::size_t limit = 10)
{
    struct RankedNode {
        Index nodeId = 0;
        Real absResidual = 0.0;
    };

    const Index nodeCount = mesh.numNodes();
    std::vector<RankedNode> ranked;
    ranked.reserve(static_cast<std::size_t>(nodeCount));
    for (Index nodeId = 0; nodeId < nodeCount; ++nodeId) {
        const Real value = residual(static_cast<int>(nodeId));
        ranked.push_back({nodeId, std::abs(value)});
    }
    std::sort(ranked.begin(), ranked.end(), [](const RankedNode& a, const RankedNode& b) {
        if (a.absResidual == b.absResidual)
            return a.nodeId < b.nodeId;
        return a.absResidual > b.absResidual;
    });

    const std::vector<Real>& ni = assembler.intrinsicDensity();
    const std::size_t count = std::min(limit, ranked.size());
    std::vector<NewtonTopResidualNode> out;
    out.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const Index nodeId = ranked[index].nodeId;
        const Node& node = mesh.getNode(nodeId);
        const Real poissonResidual = residual(static_cast<int>(nodeId));
        out.push_back({
            nodeId,
            node.x,
            node.y,
            poissonResidual,
            std::abs(poissonResidual),
            doping.donors(nodeId),
            doping.acceptors(nodeId),
            doping.netDoping(nodeId),
            nodeId < static_cast<Index>(ni.size()) ? ni[static_cast<std::size_t>(nodeId)] : 0.0});
    }
    return out;
}

NewtonFailureDiagnostics buildFailureDiagnostics(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const CoupledDDAssembler& assembler,
    const VectorXd& x,
    const VectorXd& residual,
    const std::string& failureReason,
    int failedIteration,
    Real residualNorm,
    Real stepNorm,
    Real dampingFactor,
    int lineSearchAttempts,
    const std::string& lineSearchFailureReason,
    std::vector<LineSearchIterationInfo> lineSearchHistory = {})
{
    NewtonFailureDiagnostics diagnostics;
    diagnostics.failureReason = failureReason;
    diagnostics.failedIteration = failedIteration;
    diagnostics.residualNorm = residualNorm;
    diagnostics.stepNorm = stepNorm;
    diagnostics.dampingFactor = dampingFactor;
    diagnostics.lineSearchAttempts = lineSearchAttempts;
    diagnostics.lineSearchFailureReason = lineSearchFailureReason;
    diagnostics.blockResiduals = blockResidualInfo(residual, mesh.numNodes());
    diagnostics.carrierDiagnostics = carrierDiagnostics(assembler, x);
    diagnostics.lineSearchHistory = std::move(lineSearchHistory);
    diagnostics.topPoissonResidualNodes = topPoissonResidualNodes(mesh, doping, assembler, residual);
    return diagnostics;
}

void printFailureDiagnostics(const NewtonFailureDiagnostics& diagnostics)
{
    std::cerr << "  failure_reason=" << diagnostics.failureReason
              << " line_search_reason=" << diagnostics.lineSearchFailureReason
              << " blocks=(" << diagnostics.blockResiduals.psi << ','
              << diagnostics.blockResiduals.phin << ','
              << diagnostics.blockResiduals.phip << ")"
              << " carriers_positive_finite="
              << (diagnostics.carrierDiagnostics.positiveFinite ? "1" : "0")
              << '\n';
    if (!diagnostics.topPoissonResidualNodes.empty()) {
        const NewtonTopResidualNode& node = diagnostics.topPoissonResidualNodes.front();
        std::cerr << "  top_poisson_residual_node=" << node.nodeId
                  << " x=" << node.x
                  << " y=" << node.y
                  << " residual=" << node.poissonResidual
                  << " donors=" << node.donors
                  << " acceptors=" << node.acceptors
                  << " net=" << node.netDoping
                  << " ni_eff=" << node.effectiveIntrinsicDensity
                  << '\n';
    }
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
    cfg.maxUpdate = json.value("max_update", cfg.maxUpdate);
    cfg.quasiFermiUpdateLimit_V = json.value(
        "quasi_fermi_update_limit_V",
        cfg.quasiFermiUpdateLimit_V);
    cfg.quasiFermiUpdateLimitMinority_V = json.value(
        "quasi_fermi_update_limit_minority_V",
        cfg.quasiFermiUpdateLimitMinority_V);
    cfg.stallResidualFloor = json.value("stall_residual_floor", cfg.stallResidualFloor);
    cfg.carrierRegularizationScale = json.value(
        "carrier_regularization_scale",
        cfg.carrierRegularizationScale);
    cfg.finiteDifferenceStep = json.value("finite_difference_step", cfg.finiteDifferenceStep);
    cfg.jacobian = json.value("jacobian", cfg.jacobian);
    cfg.residualNorm = json.value("residual_norm", cfg.residualNorm);
    cfg.contactBoundaryReconstruction =
        json.value("contact_boundary_reconstruction", cfg.contactBoundaryReconstruction);
    cfg.contactBoundaryMinorityElectronRelaxation = json.value(
        "contact_boundary_minority_electron_relaxation",
        cfg.contactBoundaryMinorityElectronRelaxation);
    cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V = json.value(
        "contact_boundary_minority_electron_relaxation_bias_threshold_V",
        cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V);
    cfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly = json.value(
        "contact_boundary_minority_electron_relaxation_two_terminal_only",
        cfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
    cfg.contactBoundaryMinorityElectronRelaxationContactSide = json.value(
        "contact_boundary_minority_electron_relaxation_contact_side",
        cfg.contactBoundaryMinorityElectronRelaxationContactSide);
    cfg.contactBoundaryMinorityElectronRelaxationStrength = json.value(
        "contact_boundary_minority_electron_relaxation_strength",
        cfg.contactBoundaryMinorityElectronRelaxationStrength);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    cfg.augerCn = json.value("auger_cn_m6_per_s", cfg.augerCn);
    cfg.augerCp = json.value("auger_cp_m6_per_s", cfg.augerCp);
    if (json.contains("mobility"))
        cfg.mobility = mobilityModelConfigFromJson(json.at("mobility"), scaling);
    if (json.contains("bandgap_narrowing")) {
        const auto& value = json.at("bandgap_narrowing");
        if (value.is_string()) {
            cfg.bandgapNarrowing = bandgapNarrowingConfig(value.get<std::string>());
        } else if (value.is_object()) {
            cfg.bandgapNarrowing = bandgapNarrowingConfig(
                value.value("model", cfg.bandgapNarrowing.model));
            if (value.contains("reference_doping_m3")) {
                cfg.bandgapNarrowing.referenceDoping = scaling.concentrationToSI(
                    value.at("reference_doping_m3").get<Real>());
            }
            cfg.bandgapNarrowing.coefficient = value.value(
                "coefficient_eV", cfg.bandgapNarrowing.coefficient);
            cfg.bandgapNarrowing.smoothing = value.value(
                "smoothing", cfg.bandgapNarrowing.smoothing);
            cfg.bandgapNarrowing.offset = value.value(
                "offset_eV", cfg.bandgapNarrowing.offset);
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
            cfg.impactIonization.parameterSet = value.value(
                "parameter_set", cfg.impactIonization.parameterSet);
            cfg.impactIonization.drivingForce = value.value(
                "driving_force", cfg.impactIonization.drivingForce);
            cfg.impactIonization.generation = value.value(
                "generation", cfg.impactIonization.generation);
            cfg.impactIonization.currentApproximation = value.value(
                "current_approximation", cfg.impactIonization.currentApproximation);
            cfg.impactIonization.currentMagnitudeMode = value.value(
                "current_magnitude_mode", cfg.impactIonization.currentMagnitudeMode);
            cfg.impactIonization.quasiFermiGradientDiscretization = value.value(
                "quasi_fermi_gradient_discretization",
                cfg.impactIonization.quasiFermiGradientDiscretization);
            parseImpactIonizationDrivingForceInterpolation(
                value, scaling, cfg.impactIonization, "newtonConfigFromJson");
            cfg.impactIonization.sourceGeometryScale = value.value(
                "source_geometry_scale", cfg.impactIonization.sourceGeometryScale);
            cfg.impactIonization.sourceVolumePolicy = value.value(
                "source_volume_policy", cfg.impactIonization.sourceVolumePolicy);
            cfg.impactIonization.sourceVolumeFactor = value.value(
                "source_volume_factor", cfg.impactIonization.sourceVolumeFactor);
            cfg.impactIonization.sourceMappingMode = value.value(
                "source_mapping_mode", cfg.impactIonization.sourceMappingMode);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_truncation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_trucation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.minimumField = scaling.electricFieldToSI(value.value(
                "minimum_field_V_m", cfg.impactIonization.minimumField));
            cfg.impactIonization.debugRawVanOverstraeten = value.value(
                "debug_raw_vanoverstraeten",
                cfg.impactIonization.debugRawVanOverstraeten);
            cfg.impactIonization.aScale = value.value(
                "A_scale", cfg.impactIonization.aScale);
            cfg.impactIonization.bScale = value.value(
                "B_scale", cfg.impactIonization.bScale);
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
            if (value.contains("electron_a_low_m_inv")) {
                cfg.impactIonization.electronALow = scaling.inverseLengthToSI(
                    value.at("electron_a_low_m_inv").get<Real>());
            }
            if (value.contains("electron_a_high_m_inv")) {
                cfg.impactIonization.electronAHigh = scaling.inverseLengthToSI(
                    value.at("electron_a_high_m_inv").get<Real>());
            }
            if (value.contains("electron_b_low_V_m")) {
                cfg.impactIonization.electronBLow = scaling.electricFieldToSI(
                    value.at("electron_b_low_V_m").get<Real>());
            }
            if (value.contains("electron_b_high_V_m")) {
                cfg.impactIonization.electronBHigh = scaling.electricFieldToSI(
                    value.at("electron_b_high_V_m").get<Real>());
            }
            if (value.contains("hole_a_low_m_inv")) {
                cfg.impactIonization.holeALow = scaling.inverseLengthToSI(
                    value.at("hole_a_low_m_inv").get<Real>());
            }
            if (value.contains("hole_a_high_m_inv")) {
                cfg.impactIonization.holeAHigh = scaling.inverseLengthToSI(
                    value.at("hole_a_high_m_inv").get<Real>());
            }
            if (value.contains("hole_b_low_V_m")) {
                cfg.impactIonization.holeBLow = scaling.electricFieldToSI(
                    value.at("hole_b_low_V_m").get<Real>());
            }
            if (value.contains("hole_b_high_V_m")) {
                cfg.impactIonization.holeBHigh = scaling.electricFieldToSI(
                    value.at("hole_b_high_V_m").get<Real>());
            }
            if (value.contains("switch_field_V_m")) {
                cfg.impactIonization.switchField = scaling.electricFieldToSI(
                    value.at("switch_field_V_m").get<Real>());
            }
            cfg.impactIonization.phononEnergy = value.value(
                "phonon_energy_eV", cfg.impactIonization.phononEnergy);
            cfg.impactIonization.referenceTemperature_K = value.value(
                "reference_temperature_K", cfg.impactIonization.referenceTemperature_K);
            cfg.impactIonization.temperature_K = value.value(
                "temperature_K", cfg.impactIonization.temperature_K);
            cfg.impactIonization.carrierVelocity = value.value(
                "carrier_velocity_m_s", cfg.impactIonization.carrierVelocity);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: impact_ionization must be a string or object.");
        }
    }
    detail::validateImpactIonizationDrivingForce(cfg.impactIonization, "newtonConfigFromJson");

    if (cfg.jacobian != "analytic" && cfg.jacobian != "finite_difference")
        throw std::invalid_argument(
            "newtonConfigFromJson: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg.maxUpdate < 0.0 || !std::isfinite(cfg.maxUpdate))
        throw std::invalid_argument(
            "newtonConfigFromJson: max_update must be non-negative and finite.");
    if (cfg.quasiFermiUpdateLimit_V < 0.0 || !std::isfinite(cfg.quasiFermiUpdateLimit_V))
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_update_limit_V must be non-negative and finite.");
    if (cfg.quasiFermiUpdateLimitMinority_V < 0.0 ||
        !std::isfinite(cfg.quasiFermiUpdateLimitMinority_V))
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_update_limit_minority_V must be non-negative and finite.");
    if (cfg.stallResidualFloor < 0.0 || !std::isfinite(cfg.stallResidualFloor))
        throw std::invalid_argument(
            "newtonConfigFromJson: stall_residual_floor must be non-negative and finite.");
    if (cfg.carrierRegularizationScale < 0.0 || !std::isfinite(cfg.carrierRegularizationScale))
        throw std::invalid_argument(
            "newtonConfigFromJson: carrier_regularization_scale must be non-negative and finite.");
    if (cfg.finiteDifferenceStep <= 0.0 || !std::isfinite(cfg.finiteDifferenceStep))
        throw std::invalid_argument(
            "newtonConfigFromJson: finite_difference_step must be positive and finite.");
    if (cfg.residualNorm != "block" && cfg.residualNorm != "l2")
        throw std::invalid_argument(
            "newtonConfigFromJson: residual_norm must be 'block' or 'l2'.");
    cfg.contactBoundaryReconstruction =
        canonicalContactBoundaryReconstruction(cfg.contactBoundaryReconstruction);
    cfg.contactBoundaryMinorityElectronRelaxationContactSide =
        canonicalMinorityRelaxationContactSide(
            cfg.contactBoundaryMinorityElectronRelaxationContactSide);
    cfg.contactBoundaryMinorityElectronRelaxationStrength =
        clampMinorityRelaxationStrength(
            cfg.contactBoundaryMinorityElectronRelaxationStrength);
    if (!std::isfinite(cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V)
        || cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V < 0.0) {
        throw std::invalid_argument(
            "newtonConfigFromJson: "
            "contact_boundary_minority_electron_relaxation_bias_threshold_V "
            "must be finite and non-negative.");
    }
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
    if (cfg_.stallResidualFloor < 0.0 || !std::isfinite(cfg_.stallResidualFloor)) {
        throw std::invalid_argument(
            "NewtonSolver: stall_residual_floor must be non-negative and finite.");
    }
    if (cfg_.carrierRegularizationScale < 0.0 ||
        !std::isfinite(cfg_.carrierRegularizationScale)) {
        throw std::invalid_argument(
            "NewtonSolver: carrier_regularization_scale must be non-negative and finite.");
    }
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
    return buildBoundaryConditions(assembler, contactBiases_);
}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler,
    const std::unordered_map<std::string, Real>& contactBiases) const
{
    CoupledDDBoundaryConditions bcs;
    const auto& ni = assembler.intrinsicDensity();
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const bool dominantSignedContactMean =
        cfg_.contactBoundaryReconstruction == "dominant_signed_contact_mean";
    const bool enableMinorityElectronRelaxation =
        cfg_.contactBoundaryMinorityElectronRelaxation;
    const Real relaxationBiasThreshold =
        cfg_.contactBoundaryMinorityElectronRelaxationBiasThreshold_V;
    const bool twoTerminalOnly =
        cfg_.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly;
    const std::string contactSidePolicy =
        cfg_.contactBoundaryMinorityElectronRelaxationContactSide;
    const bool allowsMinorityRelaxation =
        !twoTerminalOnly || mesh_.numContacts() == 2;

    for (Index c = 0; c < mesh_.numContacts(); ++c) {
        const Contact& contact = mesh_.getContact(c);
        auto it = contactBiases.find(contact.name);
        if (it == contactBiases.end()) continue;

        const double Vbias = it->second;
        const bool relaxByBiasAndTopology =
            enableMinorityElectronRelaxation
            && allowsMinorityRelaxation
            && (std::abs(Vbias) >= relaxationBiasThreshold);
        const bool relaxMinorityOnPContact =
            relaxByBiasAndTopology
            && (contactSidePolicy == "p_contact_only"
                || contactSidePolicy == "both_contacts");
        const bool relaxMinorityOnNContact =
            relaxByBiasAndTopology
            && (contactSidePolicy == "n_contact_only"
                || contactSidePolicy == "both_contacts");
            const Real relaxedMinorityBias =
                (1.0 - cfg_.contactBoundaryMinorityElectronRelaxationStrength) * Vbias;
        for (Index nid : contact.node_ids) {
            const double niNode = ni[nid];
            const double Ndop = ohmicContactNetDoping(
                mesh_, doping_, contact, nid, dominantSignedContactMean);
            const double neq = nEq(Ndop, niNode);
            double psiBuiltIn = 0.0;
            if (niNode > 0.0 && neq > 0.0)
                psiBuiltIn = Vt * std::log(neq / niNode);

            bcs.psi[nid] = (Vbias + psiBuiltIn) / potentialScale;
            if (Ndop >= 0.0) {
                bcs.phin[nid] = Vbias / potentialScale;
                if (relaxMinorityOnNContact)
                    bcs.phip[nid] = relaxedMinorityBias / potentialScale;
                else
                    bcs.phip[nid] = Vbias / potentialScale;
            } else {
                bcs.phip[nid] = Vbias / potentialScale;
                if (relaxMinorityOnPContact)
                    bcs.phin[nid] = relaxedMinorityBias / potentialScale;
                else
                    bcs.phin[nid] = Vbias / potentialScale;
            }
        }
    }
    return bcs;
}

DDSolution NewtonSolver::buildInitialGuess(
    const CoupledDDAssembler& assembler, const CoupledDDBoundaryConditions& bcs) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    const Real Vt = thermalVoltage(cfg_.temperature_K);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const auto& ni = assembler.intrinsicDensity();

    DDSolution sol;
    sol.psi = VectorXd::Zero(N);
    sol.phin = VectorXd::Zero(N);
    sol.phip = VectorXd::Zero(N);
    sol.n = VectorXd::Zero(N);
    sol.p = VectorXd::Zero(N);
    sol.iters = 0;
    sol.converged = false;

    for (int i = 0; i < N; ++i) {
        const Index nid = static_cast<Index>(i);
        const Real niNode = ni[nid];
        const Real neq = nEq(doping_.netDoping(nid), niNode);
        Real psiBuiltIn = 0.0;
        if (niNode > 0.0 && neq > 0.0)
            psiBuiltIn = Vt * std::log(neq / niNode);
        sol.psi(i) = psiBuiltIn;
    }

    for (const auto& [nid, value] : bcs.psi)
        sol.psi(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phin)
        sol.phin(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phip)
        sol.phip(static_cast<int>(nid)) = value * potentialScale;

    for (int i = 0; i < N; ++i) {
        const Index nid = static_cast<Index>(i);
        sol.n(i) = electronDensity(ni[nid], sol.psi(i), sol.phin(i), Vt);
        sol.p(i) = holeDensity(ni[nid], sol.psi(i), sol.phip(i), Vt);
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

std::shared_ptr<CoupledDDAssembler> NewtonSolver::makeArclengthAssembler() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    const DDScalingSpec scaling = buildScalingSpec();
    return std::make_shared<CoupledDDAssembler>(
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
}

ArclengthSystem NewtonSolver::makeArclengthSystem(const std::string& activeContact,
                                                  Real biasFiniteDifferenceStep_V) const
{
    if (!(biasFiniteDifferenceStep_V > 0.0) ||
        !std::isfinite(biasFiniteDifferenceStep_V)) {
        throw std::invalid_argument(
            "NewtonSolver::makeArclengthSystem: biasFiniteDifferenceStep_V must be "
            "finite and positive.");
    }
    if (contactBiases_.find(activeContact) == contactBiases_.end()) {
        throw std::invalid_argument(
            "NewtonSolver::makeArclengthSystem: active contact '" + activeContact +
            "' is not present in the contact bias map.");
    }

    auto assembler = makeArclengthAssembler();
    const NewtonSolver* self = this;
    const std::unordered_map<std::string, Real> baseBiases = contactBiases_;
    const Real h = biasFiniteDifferenceStep_V;
    const int nodeCount = static_cast<int>(mesh_.numNodes());
    const Real potentialScale = assembler->usesScaledState()
        ? assembler->potentialScale()
        : 1.0;

    auto biasesAt = [baseBiases, activeContact](Real lambda) {
        std::unordered_map<std::string, Real> biases = baseBiases;
        biases[activeContact] = lambda;
        return biases;
    };

    ArclengthSystem system;
    system.residual = [self, assembler, biasesAt](const VectorXd& x, Real lambda) {
        const CoupledDDBoundaryConditions bcs =
            self->buildBoundaryConditions(*assembler, biasesAt(lambda));
        return assembler->residual(x, bcs);
    };
    system.parameterDerivative =
        [self, assembler, biasesAt, h](const VectorXd& x, Real lambda) {
            const CoupledDDBoundaryConditions bcsPlus =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda + h));
            const CoupledDDBoundaryConditions bcsMinus =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda - h));
            const VectorXd fPlus = assembler->residual(x, bcsPlus);
            const VectorXd fMinus = assembler->residual(x, bcsMinus);
            return VectorXd((fPlus - fMinus) / (2.0 * h));
        };
    system.solveJacobian =
        [self, assembler, biasesAt](const VectorXd& x, Real lambda,
                                    const VectorXd& b, VectorXd& y) {
            const CoupledDDBoundaryConditions bcs =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda));
            const SparseMatrixd jacobian = assembler->assembleJacobian(x, bcs);
            Eigen::SparseLU<SparseMatrixd> lu;
            lu.compute(jacobian);
            if (lu.info() != Eigen::Success)
                return false;
            y = lu.solve(b);
            if (lu.info() != Eigen::Success)
                return false;
            return y.allFinite();
        };
    if (cfg_.quasiFermiUpdateLimit_V > 0.0 ||
        cfg_.quasiFermiUpdateLimitMinority_V > 0.0) {
        system.limitUpdate = [self, nodeCount, potentialScale](
            const VectorXd&, VectorXd& deltaX, Real&) {
            applyConfiguredQuasiFermiStepCaps(
                deltaX, self->cfg_, nodeCount, potentialScale, self->doping_);
        };
    }
    return system;
}

VectorXd NewtonSolver::packArclengthState(const DDSolution& state) const
{
    auto assembler = makeArclengthAssembler();
    const Real potentialScale =
        assembler->usesScaledState() ? assembler->potentialScale() : 1.0;
    return assembler->pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
}

DDSolution NewtonSolver::unpackArclengthState(const VectorXd& x) const
{
    auto assembler = makeArclengthAssembler();
    return makeSolution(*assembler, x, 0);
}

NewtonResidualEvaluation NewtonSolver::evaluateResidual(const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);
    const NewtonBlockResidualInfo blocks = blockResidualInfo(raw, mesh_.numNodes());

    NewtonResidualEvaluation evaluation;
    evaluation.raw = raw;
    evaluation.blockNorms = blocks;
    evaluation.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.scaledState = assembler.usesScaledState();
    evaluation.potentialScale = potentialScale;
    return evaluation;
}

NewtonStepEvaluation NewtonSolver::evaluateStep(const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    LinearSolver linearSolver;
    VectorXd step = linearSolver.solve(J, -raw);
    const Real rawStepNorm = step.norm();

    const int N = static_cast<int>(mesh_.numNodes());
    applyConfiguredStepCapsAndPoissonRecorrection(
        step, J, raw, cfg_, N, potentialScale, doping_);

    const VectorXd trialX = x + step;
    const VectorXd trialRaw = assembler.residual(trialX, bcs);

    NewtonStepEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.trialResidual.raw = trialRaw;
    evaluation.trialResidual.blockNorms = blockResidualInfo(trialRaw, mesh_.numNodes());
    evaluation.trialResidual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.trialResidual.scaledState = assembler.usesScaledState();
    evaluation.trialResidual.potentialScale = potentialScale;
    evaluation.trialSolution = makeSolution(assembler, trialX, 1);
    evaluation.deltaPsi = step.segment(0, N) * potentialScale;
    evaluation.deltaPhin = step.segment(N, N) * potentialScale;
    evaluation.deltaPhip = step.segment(2 * N, N) * potentialScale;
    evaluation.rawStepNorm = rawStepNorm;
    evaluation.stepNorm = step.norm();
    return evaluation;
}

NewtonDirectionalDerivativeEvaluation NewtonSolver::evaluateDirectionalDerivative(
    const DDSolution& state,
    const DDSolution& physicalPerturbation) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N ||
        physicalPerturbation.psi.size() != N ||
        physicalPerturbation.phin.size() != N ||
        physicalPerturbation.phip.size() != N) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateDirectionalDerivative: state and perturbation sizes "
            "must match the mesh node count.");
    }

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd dx = assembler.pack({
        physicalPerturbation.psi / potentialScale,
        physicalPerturbation.phin / potentialScale,
        physicalPerturbation.phip / potentialScale});
    if (dx.norm() == 0.0)
        throw std::invalid_argument(
            "NewtonSolver::evaluateDirectionalDerivative: perturbation must be non-zero.");

    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = assembler.assembleJacobian(x, bcs);
    const VectorXd analytic = J * dx;
    const VectorXd forward = assembler.residual(x + dx, bcs);
    const VectorXd backward = assembler.residual(x - dx, bcs);
    const VectorXd finiteDifference = 0.5 * (forward - backward);
    const VectorXd error = analytic - finiteDifference;

    NewtonDirectionalDerivativeEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.perturbationPsi = physicalPerturbation.psi;
    evaluation.perturbationPhin = physicalPerturbation.phin;
    evaluation.perturbationPhip = physicalPerturbation.phip;
    evaluation.analyticJv = analytic;
    evaluation.finiteDifferenceJv = finiteDifference;
    evaluation.forwardResidual = forward;
    evaluation.backwardResidual = backward;
    evaluation.perturbationNorm = dx.norm();
    evaluation.analyticNorm = analytic.norm();
    evaluation.finiteDifferenceNorm = finiteDifference.norm();
    evaluation.absoluteError = error.norm();
    evaluation.relativeError = evaluation.absoluteError /
        std::max<Real>(1.0, evaluation.finiteDifferenceNorm);
    return evaluation;
}

NewtonBlockStepEvaluation NewtonSolver::evaluateBlockStep(
    const DDSolution& state,
    const std::string& mode) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd step = VectorXd::Zero(3 * N);
    LinearSolver linearSolver;
    if (mode == "poisson_only") {
        const SparseMatrixd block = sparseBlock(J, 0, 0, N, N);
        step.segment(0, N) = linearSolver.solve(block, -raw.segment(0, N));
    } else if (mode == "carrier_only") {
        const SparseMatrixd block = sparseBlock(J, N, N, 2 * N, 2 * N);
        step.segment(N, 2 * N) = linearSolver.solve(block, -raw.segment(N, 2 * N));
    } else {
        throw std::invalid_argument(
            "NewtonSolver::evaluateBlockStep: mode must be 'poisson_only' "
            "or 'carrier_only'.");
    }

    const Real rawStepNorm = step.norm();
    applyConfiguredStepCaps(step, cfg_, N, potentialScale, doping_);
    const VectorXd trialX = x + step;
    const VectorXd trialRaw = assembler.residual(trialX, bcs);

    NewtonBlockStepEvaluation evaluation;
    evaluation.mode = mode;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.trialResidual.raw = trialRaw;
    evaluation.trialResidual.blockNorms = blockResidualInfo(trialRaw, mesh_.numNodes());
    evaluation.trialResidual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.trialResidual.scaledState = assembler.usesScaledState();
    evaluation.trialResidual.potentialScale = potentialScale;
    evaluation.trialSolution = makeSolution(assembler, trialX, 1);
    evaluation.deltaPsi = step.segment(0, N) * potentialScale;
    evaluation.deltaPhin = step.segment(N, N) * potentialScale;
    evaluation.deltaPhip = step.segment(2 * N, N) * potentialScale;
    evaluation.rawStepNorm = rawStepNorm;
    evaluation.stepNorm = step.norm();
    return evaluation;
}

NewtonRegularizedCarrierStepEvaluation NewtonSolver::evaluateRegularizedCarrierStep(
    const DDSolution& state,
    Real regularizationScale) const
{
    if (!std::isfinite(regularizationScale) || regularizationScale < 0.0) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateRegularizedCarrierStep: "
            "regularization scale must be finite and non-negative.");
    }

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    const int N = static_cast<int>(mesh_.numNodes());
    const SparseMatrixd carrierBlock = sparseBlock(J, N, N, 2 * N, 2 * N);
    std::vector<Real> rowAbsSums(static_cast<std::size_t>(2 * N), 0.0);
    for (int col = 0; col < carrierBlock.outerSize(); ++col) {
        for (SparseMatrixd::InnerIterator it(carrierBlock, col); it; ++it) {
            rowAbsSums[static_cast<std::size_t>(it.row())] += std::abs(it.value());
        }
    }

    SparseMatrixd regularizedBlock = carrierBlock;
    Real regularizationDiagonalNormSq = 0.0;
    for (int row = 0; row < 2 * N; ++row) {
        const Real diagonal = carrierBlock.coeff(row, row);
        const Real sign = diagonal < 0.0 ? -1.0 : 1.0;
        const Real addition =
            sign * regularizationScale * rowAbsSums[static_cast<std::size_t>(row)];
        regularizedBlock.coeffRef(row, row) += addition;
        regularizationDiagonalNormSq += addition * addition;
    }
    regularizedBlock.makeCompressed();

    VectorXd step = VectorXd::Zero(3 * N);
    LinearSolver linearSolver;
    step.segment(N, 2 * N) =
        linearSolver.solve(regularizedBlock, -raw.segment(N, 2 * N));

    const Real rawStepNorm = step.norm();
    applyConfiguredStepCaps(step, cfg_, N, potentialScale, doping_);
    const VectorXd trialX = x + step;
    const VectorXd trialRaw = assembler.residual(trialX, bcs);

    NewtonRegularizedCarrierStepEvaluation evaluation;
    evaluation.regularizationScale = regularizationScale;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.trialResidual.raw = trialRaw;
    evaluation.trialResidual.blockNorms = blockResidualInfo(trialRaw, mesh_.numNodes());
    evaluation.trialResidual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.trialResidual.scaledState = assembler.usesScaledState();
    evaluation.trialResidual.potentialScale = potentialScale;
    evaluation.trialSolution = makeSolution(assembler, trialX, 1);
    evaluation.deltaPsi = step.segment(0, N) * potentialScale;
    evaluation.deltaPhin = step.segment(N, N) * potentialScale;
    evaluation.deltaPhip = step.segment(2 * N, N) * potentialScale;
    evaluation.rawStepNorm = rawStepNorm;
    evaluation.stepNorm = step.norm();
    evaluation.regularizationDiagonalNorm = std::sqrt(regularizationDiagonalNormSq);
    return evaluation;
}

NewtonCarrierRowDiagnosticsEvaluation NewtonSolver::evaluateCarrierRowDiagnostics(
    const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    const int N = static_cast<int>(mesh_.numNodes());
    const SparseMatrixd carrierBlock = sparseBlock(J, N, N, 2 * N, 2 * N);
    LinearSolver linearSolver;
    VectorXd rawStep = VectorXd::Zero(3 * N);
    rawStep.segment(N, 2 * N) =
        linearSolver.solve(carrierBlock, -raw.segment(N, 2 * N));
    VectorXd cappedStep = rawStep;
    applyConfiguredStepCaps(cappedStep, cfg_, N, potentialScale, doping_);

    std::vector<Real> electronRowAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> holeRowAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronRowL2Sq(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> holeRowL2Sq(static_cast<std::size_t>(N), 0.0);
    for (int col = 0; col < J.outerSize(); ++col) {
        for (SparseMatrixd::InnerIterator it(J, col); it; ++it) {
            const int row = static_cast<int>(it.row());
            const Real value = it.value();
            if (row >= N && row < 2 * N) {
                const std::size_t node = static_cast<std::size_t>(row - N);
                electronRowAbs[node] += std::abs(value);
                electronRowL2Sq[node] += value * value;
            } else if (row >= 2 * N && row < 3 * N) {
                const std::size_t node = static_cast<std::size_t>(row - 2 * N);
                holeRowAbs[node] += std::abs(value);
                holeRowL2Sq[node] += value * value;
            }
        }
    }

    NewtonCarrierRowDiagnosticsEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.potentialScale = potentialScale;
    evaluation.rawCarrierStepNorm = rawStep.norm();
    evaluation.cappedCarrierStepNorm = cappedStep.norm();
    evaluation.rows.reserve(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        const int eRow = N + i;
        const int hRow = 2 * N + i;
        const Real eDiag = J.coeff(eRow, eRow);
        const Real hDiag = J.coeff(hRow, hRow);
        NewtonCarrierRowDiagnostic row;
        row.nodeId = static_cast<Index>(i);
        row.electronResidual = raw(eRow);
        row.holeResidual = raw(hRow);
        row.electronDiagonal = eDiag;
        row.holeDiagonal = hDiag;
        row.electronRowAbsSum = electronRowAbs[static_cast<std::size_t>(i)];
        row.holeRowAbsSum = holeRowAbs[static_cast<std::size_t>(i)];
        row.electronOffdiagAbsSum = row.electronRowAbsSum - std::abs(eDiag);
        row.holeOffdiagAbsSum = row.holeRowAbsSum - std::abs(hDiag);
        row.electronRowL2Norm =
            std::sqrt(electronRowL2Sq[static_cast<std::size_t>(i)]);
        row.holeRowL2Norm =
            std::sqrt(holeRowL2Sq[static_cast<std::size_t>(i)]);
        row.rawDeltaPhin_V = rawStep(N + i) * potentialScale;
        row.rawDeltaPhip_V = rawStep(2 * N + i) * potentialScale;
        row.cappedDeltaPhin_V = cappedStep(N + i) * potentialScale;
        row.cappedDeltaPhip_V = cappedStep(2 * N + i) * potentialScale;
        evaluation.rows.push_back(row);
    }
    return evaluation;
}

NewtonCarrierTermDiagnosticsEvaluation NewtonSolver::evaluateCarrierTermDiagnostics(
    const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = assembler.residual(x, bcs);

    NewtonCarrierTermDiagnosticsEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.rows = assembler.carrierContinuityTermDiagnostics(x, bcs);
    return evaluation;
}

std::vector<NewtonJacobianBlockAuditRow> NewtonSolver::evaluateJacobianBlockAudit(
    const DDSolution& state,
    Real finiteDifferenceStep,
    std::vector<std::string> blocks) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateJacobianBlockAudit: state size must match the mesh node count.");
    }
    if (finiteDifferenceStep <= 0.0 || !std::isfinite(finiteDifferenceStep)) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateJacobianBlockAudit: finite difference step must be positive.");
    }

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    const DDScalingSpec scaling = buildScalingSpec();

    const auto makeRecombinationConfig =
        [&](const std::vector<std::string>& models) {
            RecombinationModelConfig config =
                recombinationModelConfig(models, cfg_.taun, cfg_.taup);
            config.augerCn = cfg_.augerCn;
            config.augerCp = cfg_.augerCp;
            return config;
        };
    const RecombinationModelConfig noRecombinationConfig =
        makeRecombinationConfig({"none"});
    const RecombinationModelConfig recombinationConfig =
        makeRecombinationConfig(cfg_.recombination);
    const ImpactIonizationModelConfig noImpactConfig{};

    const auto makeAssembler =
        [&](const RecombinationModelConfig& recombination,
            const ImpactIonizationModelConfig& impact) {
            return CoupledDDAssembler(
                mesh_,
                matdb_,
                doping_,
                Vt,
                mobilityConfig,
                recombination,
                cfg_.bandgapNarrowing,
                impact,
                fixedCharges_,
                sheetCharges_,
                scaling);
        };

    CoupledDDAssembler baseAssembler =
        makeAssembler(noRecombinationConfig, noImpactConfig);
    const CoupledDDBoundaryConditions bcs =
        buildBoundaryConditions(baseAssembler);
    const Real potentialScale =
        baseAssembler.usesScaledState() ? baseAssembler.potentialScale() : 1.0;
    const VectorXd x = baseAssembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});

    const auto matrixPair =
        [&](CoupledDDAssembler& assembler) {
            return std::pair<SparseMatrixd, SparseMatrixd>{
                assembler.assembleJacobian(x, bcs),
                assembler.finiteDifferenceJacobian(x, bcs, finiteDifferenceStep),
            };
        };

    const std::vector<std::string> defaultBlocks = {
        "poisson",
        "transport",
        "srh_auger",
        "sg_avalanche",
        "dirichlet_or_gauge",
    };
    if (blocks.empty())
        blocks = defaultBlocks;

    const auto wants = [&](const std::string& block) {
        return std::find(blocks.begin(), blocks.end(), block) != blocks.end();
    };
    const bool needsBase =
        wants("poisson") || wants("transport") || wants("srh_auger") ||
        wants("sg_avalanche") || wants("dirichlet_or_gauge");
    const bool needsRecombination = wants("srh_auger");
    const bool needsImpact = wants("sg_avalanche");

    std::optional<std::pair<SparseMatrixd, SparseMatrixd>> base;
    std::optional<std::pair<SparseMatrixd, SparseMatrixd>> withRecombination;
    std::optional<std::pair<SparseMatrixd, SparseMatrixd>> withImpact;
    if (needsBase)
        base = matrixPair(baseAssembler);
    if (needsRecombination) {
        CoupledDDAssembler recombinationAssembler =
            makeAssembler(recombinationConfig, noImpactConfig);
        withRecombination = matrixPair(recombinationAssembler);
    }
    if (needsImpact) {
        CoupledDDAssembler impactAssembler =
            makeAssembler(noRecombinationConfig, cfg_.impactIonization);
        withImpact = matrixPair(impactAssembler);
    }

    std::vector<NewtonJacobianBlockAuditRow> rows;
    rows.reserve(blocks.size());
    for (const std::string& block : blocks) {
        if (block == "poisson") {
            rows.push_back(jacobianAuditRow(
                block, base->first, base->second, jacobianAuditRows(block, N)));
        } else if (block == "transport") {
            rows.push_back(jacobianAuditRow(
                block, base->first, base->second, jacobianAuditRows(block, N)));
        } else if (block == "srh_auger") {
            rows.push_back(jacobianAuditRow(
                block,
                withRecombination->first - base->first,
                withRecombination->second - base->second,
                jacobianAuditRows(block, N)));
        } else if (block == "sg_avalanche") {
            rows.push_back(jacobianAuditRow(
                block,
                withImpact->first - base->first,
                withImpact->second - base->second,
                jacobianAuditRows(block, N)));
        } else if (block == "dirichlet_or_gauge") {
            rows.push_back(jacobianAuditRow(
                block, base->first, base->second, jacobianAuditRows(block, N, bcs)));
        } else {
            throw std::invalid_argument(
                "NewtonSolver::evaluateJacobianBlockAudit: unknown block '" + block + "'.");
        }
    }
    return rows;
}

std::vector<CoupledDDEdgeFluxDiagnostic> NewtonSolver::evaluateSgEdgeFluxDiagnostics(
    const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    return assembler.sgEdgeFluxDiagnostics(x, bcs);
}

NewtonResult NewtonSolver::solve() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
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
    VectorXd psiInit = initial.psi;
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
    for (const auto& [nid, value] : bcs.psi)
        psiInit(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phin)
        phinInit(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phip)
        phipInit(static_cast<int>(nid)) = value * potentialScale;
    VectorXd x = assembler.pack({
        psiInit / potentialScale,
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
        const ResidualBlockNormValue blocks = ResidualNorm::computeBlocks(r, mesh_.numNodes());
        std::cout << "Newton iter 0 residual=" << initialNorm
                  << " step=0 damping=0"
                  << " blocks=(" << blocks.psi << ',' << blocks.phin << ','
                  << blocks.phip << ")\n";
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

    // Stall-recovery threshold: when the damped Newton step cannot reduce the
    // residual but the residual already sits at or below this normalized floor,
    // the state is effectively solved (the line search is fighting numerical
    // noise) and convergence is reported instead of failing.
    const Real stallResidualFloor = cfg_.stallResidualFloor;


    auto capConfiguredStep = [&](VectorXd& candidateStep,
                                 const SparseMatrixd& jacobian,
                                 const VectorXd& residual) {
        applyConfiguredStepCapsAndPoissonRecorrection(
            candidateStep, jacobian, residual, cfg_, N, potentialScale, doping_);
    };

    VectorXd acceptedX = x;
    VectorXd acceptedR = r;
    int acceptedIters = 0;

    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        SparseMatrixd J = (cfg_.jacobian == "finite_difference")
            ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
            : assembler.assembleJacobian(x, bcs);
        addCarrierRowRegularization(J, N, cfg_.carrierRegularizationScale);
        VectorXd step;
        try {
            step = linearSolver.solve(J, -r);
        } catch (const std::runtime_error&) {
            result.finalResidualNorm = residualNormFn(acceptedR);
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            result.failureDiagnostics = buildFailureDiagnostics(
                mesh_,
                doping_,
                assembler,
                acceptedX,
                acceptedR,
                "linear_solve_failed",
                iter,
                result.finalResidualNorm,
                0.0,
                0.0,
                0,
                {});
            if (cfg_.verbose) {
                std::cerr << "Newton failed at iter " << iter
                          << ": residual=" << residualNormFn(r)
                          << " damping=0 step=0 (linear solve failed)\n";
                printFailureDiagnostics(result.failureDiagnostics);
            }
            return result;
        }
        capConfiguredStep(step, J, r);
        Real stepNorm = step.norm();

        auto ls = lineSearch.search(
            x, step, r,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            },
            residualNormFn);

        if (!ls.accepted) {
            const Real stalledNorm = residualNormFn(acceptedR);
            // Effectively-solved state: the residual already sits at the
            // numerical floor, so the rejected step is only fighting noise.
            // Declaring convergence here avoids spurious failures when the
            // Newton iterate has already reached the achievable precision.
            if (stalledNorm <= stallResidualFloor) {
                result.converged = true;
                result.iters = acceptedIters;
                result.finalResidualNorm = stalledNorm;
                result.solution = makeSolution(assembler, acceptedX, acceptedIters);
                return result;
            }
            result.finalResidualNorm = stalledNorm;
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            result.failureDiagnostics = buildFailureDiagnostics(
                mesh_,
                doping_,
                assembler,
                acceptedX,
                acceptedR,
                ls.failureReason.empty() ? std::string("line_search_rejected") : ls.failureReason,
                iter,
                result.finalResidualNorm,
                stepNorm,
                ls.damping,
                ls.attempts,
                ls.failureReason,
                std::move(ls.history));
            if (cfg_.verbose) {
                std::cerr << "Newton failed at iter " << iter
                          << ": residual=" << ls.residualNorm
                          << " damping=" << ls.damping
                          << " step=" << stepNorm
                          << " (line search rejected step; reason="
                          << result.failureDiagnostics.failureReason << ")\n";
                printFailureDiagnostics(result.failureDiagnostics);
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
        info.blockResiduals = blockResidualInfo(r, mesh_.numNodes());
        if (cfg_.diagnostics)
            info.lineSearchHistory = std::move(ls.history);
        result.history.push_back(std::move(info));
        if (cfg_.verbose) {
            const NewtonBlockResidualInfo& blocks = result.history.back().blockResiduals;
            std::cout << "Newton iter " << iter
                      << " residual=" << residualNorm
                      << " step=" << appliedStepNorm
                      << " damping=" << ls.damping
                      << " blocks=(" << blocks.psi << ',' << blocks.phin << ','
                      << blocks.phip << ")\n";
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

    result.iters = acceptedIters;
    result.finalResidualNorm = residualNormFn(acceptedR);
    result.solution = makeSolution(assembler, acceptedX, acceptedIters);
    if (result.finalResidualNorm <= stallResidualFloor) {
        result.converged = true;
        return result;
    }

    result.converged = false;
    result.failureDiagnostics = buildFailureDiagnostics(
        mesh_,
        doping_,
        assembler,
        acceptedX,
        acceptedR,
        "max_iterations",
        cfg_.maxIter,
        result.finalResidualNorm,
        result.history.empty() ? 0.0 : result.history.back().stepNorm,
        result.history.empty() ? 0.0 : result.history.back().dampingFactor,
        result.history.empty() ? 0 : result.history.back().lineSearchAttempts,
        {});
    if (cfg_.verbose) {
        std::cerr << "Newton failed after " << cfg_.maxIter
                  << " iterations: residual=" << result.finalResidualNorm
                  << " damping="
                  << (result.history.empty() ? 0.0 : result.history.back().dampingFactor)
                  << " step="
                  << (result.history.empty() ? 0.0 : result.history.back().stepNorm)
                  << '\n';
        printFailureDiagnostics(result.failureDiagnostics);
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
