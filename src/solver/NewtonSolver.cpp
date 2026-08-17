#include "vela/solver/NewtonSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/RuntimeLog.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace vela {
namespace {

void emitVerboseLine(const std::string& line)
{
    if (runtimeLogEnabled())
        runtimeLogInfo(line);
    else
        std::cout << line << '\n';
}

void emitVerboseErrorLine(const std::string& line)
{
    if (runtimeLogEnabled())
        runtimeLogError(line);
    else
        std::cerr << line << '\n';
}

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
        config.electronDrivingForceRefDensity = scaling.concentrationToInternal(
            interpolation.at("electron_ref_density_m3").get<Real>());
    }
    if (interpolation.contains("hole_ref_density_m3")) {
        config.holeDrivingForceRefDensity = scaling.concentrationToInternal(
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

NewtonResidualPeak residualPeak(
    const VectorXd& residual,
    Eigen::Index offset,
    Index nodeCount)
{
    NewtonResidualPeak peak;
    peak.nodeId = -1;
    for (Index node = 0; node < nodeCount; ++node) {
        const Real value =
            residual(offset + static_cast<Eigen::Index>(node));
        const Real magnitude = std::abs(value);
        if (peak.nodeId < 0 || magnitude > peak.absoluteResidual) {
            peak.nodeId = node;
            peak.signedResidual = value;
            peak.absoluteResidual = magnitude;
        }
    }
    return peak;
}

void fillNewtonIterationResidualTrace(
    NewtonIterationInfo& info,
    const VectorXd& state,
    const VectorXd& residual,
    const VectorXd& rowWeights,
    Index nodeCount,
    const CoupledDDAssembler& assembler,
    bool captureActiveBranches)
{
    info.blockResiduals = blockResidualInfo(residual, nodeCount);
    info.rowScaledBlockResiduals =
        blockResidualInfo(residual.cwiseProduct(rowWeights), nodeCount);
    info.topPoissonResidual = residualPeak(residual, 0, nodeCount);
    info.topElectronResidual = residualPeak(
        residual, static_cast<Eigen::Index>(nodeCount), nodeCount);
    info.topHoleResidual = residualPeak(
        residual, 2 * static_cast<Eigen::Index>(nodeCount), nodeCount);
    info.sourceJacobianActiveBranchFingerprint = captureActiveBranches
        ? assembler.impactIonizationActiveBranchFingerprint(state)
        : assembler.impactIonizationConfigurationFingerprint();
}

bool isPoissonLineSearchStall(const LineSearchResult& lineSearch,
                              const NewtonBlockResidualInfo& blocks,
                              Real stalledNorm,
                              Real maxContactMajorityQfDrop,
                              const NewtonConfig& cfg)
{
    if (lineSearch.failureReason != "line_search_non_decrease" ||
        !lineSearch.bestRejectedCandidate ||
        !std::isfinite(stalledNorm) ||
        !std::isfinite(lineSearch.bestRejectedResidualNorm) ||
        stalledNorm > cfg.poissonLineSearchStallResidualFloor) {
        return false;
    }

    const Real allowedResidual = cfg.poissonLineSearchStallResidualFloor *
        (1.0 + cfg.poissonLineSearchStallRelativeIncrease);
    if (lineSearch.bestRejectedResidualNorm > allowedResidual)
        return false;

    if (!std::isfinite(blocks.psi) || !std::isfinite(blocks.phin) ||
        !std::isfinite(blocks.phip) || !std::isfinite(blocks.combined)) {
        return false;
    }
    if (blocks.psi > cfg.poissonLineSearchStallResidualFloor)
        return false;
    if (blocks.phin > cfg.poissonLineSearchStallCarrierResidualFloor ||
        blocks.phip > cfg.poissonLineSearchStallCarrierResidualFloor) {
        return false;
    }
    if (cfg.poissonLineSearchStallContactMajorityQfDropLimit_V > 0.0 &&
        (!std::isfinite(maxContactMajorityQfDrop) ||
         maxContactMajorityQfDrop > cfg.poissonLineSearchStallContactMajorityQfDropLimit_V)) {
        return false;
    }

    return blocks.combined <= 0.0 || blocks.psi >= 0.5 * blocks.combined;
}

bool isCarrierRowQualifiedLineSearchStall(
    const LineSearchResult& lineSearch,
    const NewtonBlockResidualInfo& blocks,
    Real stalledNorm,
    Real maxContactMajorityQfDrop,
    const NewtonCarrierRowConvergenceEvaluation& rowEvaluation,
    const NewtonConfig& cfg)
{
    if (!cfg.carrierRowQualifiedStallAcceptance ||
        !rowEvaluation.enforced || !rowEvaluation.satisfied ||
        lineSearch.failureReason != "line_search_non_decrease" ||
        !lineSearch.bestRejectedCandidate ||
        !std::isfinite(stalledNorm) ||
        !(stalledNorm > 0.0) ||
        !std::isfinite(lineSearch.bestRejectedResidualNorm)) {
        return false;
    }
    const Real allowedResidual = stalledNorm *
        (1.0 + cfg.poissonLineSearchStallRelativeIncrease);
    if (lineSearch.bestRejectedResidualNorm > allowedResidual)
        return false;
    // Carrier rows are governed by the enforced source-relative row metric.
    // Retaining an absolute carrier-block ceiling here would reintroduce the
    // scale dependence this opt-in acceptance path is intended to remove.
    if (!std::isfinite(blocks.psi) || !std::isfinite(blocks.phin) ||
        !std::isfinite(blocks.phip) ||
        blocks.psi > cfg.poissonLineSearchStallResidualFloor) {
        return false;
    }
    if (cfg.poissonLineSearchStallContactMajorityQfDropLimit_V > 0.0 &&
        (!std::isfinite(maxContactMajorityQfDrop) ||
         maxContactMajorityQfDrop >
             cfg.poissonLineSearchStallContactMajorityQfDropLimit_V)) {
        return false;
    }
    return true;
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
                          const std::vector<int>& rows,
                          int columnStart = 0,
                          int columnCount = -1)
{
    std::vector<char> rowMask(static_cast<std::size_t>(matrix.rows()), 0);
    for (int row : rows) {
        if (row >= 0 && row < matrix.rows())
            rowMask[static_cast<std::size_t>(row)] = 1;
    }
    const int matrixColumns = static_cast<int>(matrix.cols());
    const int columnEnd = columnCount < 0
        ? matrixColumns
        : std::min(matrixColumns, columnStart + columnCount);

    Real sum = 0.0;
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator it(matrix, outer); it; ++it) {
            if (it.col() >= columnStart && it.col() < columnEnd &&
                rowMask[static_cast<std::size_t>(it.row())]) {
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
    const Real referenceNorm = std::max(row.analyticNorm, row.fdNorm);
    row.relDiff = referenceNorm > 0.0 ? row.diffNorm / referenceNorm : 0.0;

    const int nodeCount = analytic.cols() / 3;
    auto setColumnBlock = [&](int blockIndex,
                              Real& analyticNorm,
                              Real& fdNorm,
                              Real& diffNorm,
                              Real& relativeDiff) {
        const int start = blockIndex * nodeCount;
        analyticNorm = restrictedSparseNorm(analytic, rows, start, nodeCount);
        fdNorm = restrictedSparseNorm(fd, rows, start, nodeCount);
        diffNorm = restrictedSparseNorm(diff, rows, start, nodeCount);
        const Real blockReference = std::max(analyticNorm, fdNorm);
        relativeDiff = blockReference > 0.0 ? diffNorm / blockReference : 0.0;
    };
    setColumnBlock(0, row.analyticPsiColumnNorm, row.fdPsiColumnNorm,
                   row.diffPsiColumnNorm, row.relPsiColumnDiff);
    setColumnBlock(1, row.analyticPhinColumnNorm, row.fdPhinColumnNorm,
                   row.diffPhinColumnNorm, row.relPhinColumnDiff);
    setColumnBlock(2, row.analyticPhipColumnNorm, row.fdPhipColumnNorm,
                   row.diffPhipColumnNorm, row.relPhipColumnDiff);

    const auto carrierRows = [&](int blockIndex) {
        const int start = blockIndex * nodeCount;
        const int end = start + nodeCount;
        std::vector<int> selected;
        selected.reserve(rows.size());
        for (const int candidate : rows) {
            if (candidate >= start && candidate < end)
                selected.push_back(candidate);
        }
        return selected;
    };
    const auto setSubBlock = [&](const std::vector<int>& selectedRows,
                                 int columnBlock,
                                 Real& analyticNorm,
                                 Real& fdNorm,
                                 Real& diffNorm,
                                 Real& relativeDiff) {
        const int columnStart = columnBlock * nodeCount;
        analyticNorm = restrictedSparseNorm(
            analytic, selectedRows, columnStart, nodeCount);
        fdNorm = restrictedSparseNorm(fd, selectedRows, columnStart, nodeCount);
        diffNorm = restrictedSparseNorm(diff, selectedRows, columnStart, nodeCount);
        const Real blockReference = std::max(analyticNorm, fdNorm);
        relativeDiff = blockReference > 0.0 ? diffNorm / blockReference : 0.0;
    };
    const std::vector<int> electronRows = carrierRows(1);
    const std::vector<int> holeRows = carrierRows(2);
    setSubBlock(electronRows, 1,
                row.analyticElectronPhinNorm,
                row.fdElectronPhinNorm,
                row.diffElectronPhinNorm,
                row.relElectronPhinDiff);
    setSubBlock(electronRows, 2,
                row.analyticElectronPhipNorm,
                row.fdElectronPhipNorm,
                row.diffElectronPhipNorm,
                row.relElectronPhipDiff);
    setSubBlock(holeRows, 1,
                row.analyticHolePhinNorm,
                row.fdHolePhinNorm,
                row.diffHolePhinNorm,
                row.relHolePhinDiff);
    setSubBlock(holeRows, 2,
                row.analyticHolePhipNorm,
                row.fdHolePhipNorm,
                row.diffHolePhipNorm,
                row.relHolePhipDiff);
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

NewtonMatrixConditionEstimate matrixConditionEstimate(
    const Eigen::MatrixXd& matrix)
{
    NewtonMatrixConditionEstimate estimate;
    estimate.rows = static_cast<Index>(matrix.rows());
    estimate.columns = static_cast<Index>(matrix.cols());
    if (matrix.rows() == 0 || matrix.cols() == 0)
        return estimate;

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VectorXd singular = svd.singularValues();
    if (singular.size() == 0)
        return estimate;
    estimate.largestSingularValue = singular(0);
    const Real threshold =
        std::numeric_limits<Real>::epsilon() *
        static_cast<Real>(std::max(matrix.rows(), matrix.cols())) *
        estimate.largestSingularValue;
    for (int i = 0; i < singular.size(); ++i) {
        if (singular(i) > threshold) {
            ++estimate.numericalRank;
            estimate.smallestResolvedSingularValue = singular(i);
        }
    }
    if (estimate.smallestResolvedSingularValue > 0.0) {
        estimate.resolvedConditionNumber =
            estimate.largestSingularValue /
            estimate.smallestResolvedSingularValue;
    }
    return estimate;
}

Eigen::MatrixXd l2Equilibrated(const Eigen::MatrixXd& matrix)
{
    Eigen::MatrixXd equilibrated = matrix;
    for (int row = 0; row < equilibrated.rows(); ++row) {
        const Real norm = equilibrated.row(row).norm();
        if (norm > 0.0)
            equilibrated.row(row) /= norm;
    }
    for (int column = 0; column < equilibrated.cols(); ++column) {
        const Real norm = equilibrated.col(column).norm();
        if (norm > 0.0)
            equilibrated.col(column) /= norm;
    }
    return equilibrated;
}

Real relativeVectorDifference(
    const VectorXd& left,
    const VectorXd& right)
{
    const Real reference = std::max<Real>(
        1.0, std::max(left.norm(), right.norm()));
    return (left - right).norm() / reference;
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
    const std::vector<Real>& ni = assembler.intrinsicDensity();
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
        const bool transportCarrierNode = ni[static_cast<std::size_t>(i)] > 0.0;
        if (transportCarrierNode && !(n(i) > 0.0))
            ++diagnostics.nonpositiveElectronCount;
        if (transportCarrierNode && !(p(i) > 0.0))
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

std::vector<NewtonTopCarrierResidualNode> topCarrierResidualNodes(
    const DeviceMesh& mesh,
    const VectorXd& residual,
    Index blockOffset,
    std::size_t limit = 10)
{
    std::vector<NewtonTopCarrierResidualNode> ranked;
    ranked.reserve(static_cast<std::size_t>(mesh.numNodes()));
    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        const Node& node = mesh.getNode(nodeId);
        const Real value = residual(static_cast<int>(blockOffset + nodeId));
        ranked.push_back({nodeId, node.x, node.y, value, std::abs(value)});
    }
    std::sort(
        ranked.begin(), ranked.end(),
        [](const NewtonTopCarrierResidualNode& a,
           const NewtonTopCarrierResidualNode& b) {
            if (a.absResidual == b.absResidual)
                return a.nodeId < b.nodeId;
            return a.absResidual > b.absResidual;
        });
    if (ranked.size() > limit)
        ranked.resize(limit);
    return ranked;
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
    Real maxContactMajorityQfDrop = 0.0,
    Real bestRejectedContactMajorityQfDrop = 0.0,
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
    diagnostics.maxContactMajorityQfDrop = maxContactMajorityQfDrop;
    diagnostics.bestRejectedContactMajorityQfDrop = bestRejectedContactMajorityQfDrop;
    diagnostics.lineSearchHistory = std::move(lineSearchHistory);
    diagnostics.topPoissonResidualNodes = topPoissonResidualNodes(mesh, doping, assembler, residual);
    diagnostics.topElectronResidualNodes = topCarrierResidualNodes(
        mesh, residual, mesh.numNodes());
    diagnostics.topHoleResidualNodes = topCarrierResidualNodes(
        mesh, residual, 2 * mesh.numNodes());
    return diagnostics;
}

void printFailureDiagnostics(const NewtonFailureDiagnostics& diagnostics)
{
    emitVerboseErrorLine(
        "  failure_reason=" + diagnostics.failureReason +
        " line_search_reason=" + diagnostics.lineSearchFailureReason +
        " blocks=(" + std::to_string(diagnostics.blockResiduals.psi) + "," +
        std::to_string(diagnostics.blockResiduals.phin) + "," +
        std::to_string(diagnostics.blockResiduals.phip) + ")" +
        " max_contact_majority_qf_drop_V=" +
        std::to_string(diagnostics.maxContactMajorityQfDrop) +
        " best_rejected_contact_majority_qf_drop_V=" +
        std::to_string(diagnostics.bestRejectedContactMajorityQfDrop) +
        " carriers_positive_finite=" +
        std::string(diagnostics.carrierDiagnostics.positiveFinite ? "1" : "0"));
    if (!diagnostics.topPoissonResidualNodes.empty()) {
        const NewtonTopResidualNode& node = diagnostics.topPoissonResidualNodes.front();
        emitVerboseErrorLine(
            "  top_poisson_residual_node=" + std::to_string(node.nodeId) +
            " x=" + std::to_string(node.x) +
            " y=" + std::to_string(node.y) +
            " residual=" + std::to_string(node.poissonResidual) +
            " donors=" + std::to_string(node.donors) +
            " acceptors=" + std::to_string(node.acceptors) +
            " net=" + std::to_string(node.netDoping) +
            " ni_eff=" + std::to_string(node.effectiveIntrinsicDensity));
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

NewtonCarrierRowConvergenceEvaluation evaluateCarrierRowConvergence(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const NewtonCarrierRowConvergenceConfig& cfg)
{
    NewtonCarrierRowConvergenceEvaluation evaluation;
    evaluation.enabled = cfg.mode != "off";
    evaluation.enforced = cfg.mode == "enforce";
    evaluation.epsRow = cfg.epsRow;
    if (!evaluation.enabled)
        return evaluation;

    const Real eps = cfg.epsRow;
    const Real floor = std::max<Real>(cfg.scaleFloor, 0.0);
    auto checkCarrier = [&](const CoupledDDCarrierTermDiagnostic& row,
                            const std::string& carrier,
                            Real residual,
                            Real flux,
                            Real fluxAbsSum,
                            Real recombination,
                            Real impact) {
        const Real sourceScale = std::max(std::abs(recombination), std::abs(impact));
        const Real scale = std::max({std::abs(fluxAbsSum), sourceScale, floor});
        const bool sourceQualified = scale > 0.0 &&
            sourceScale >= cfg.minSourceScale &&
            sourceScale >= cfg.minSourceScaleFraction * scale;
        const Real ratio = scale > 0.0 ? std::abs(residual) / scale : 0.0;
        if (sourceQualified && ratio > evaluation.maxRatio) {
            evaluation.maxRatio = ratio;
            evaluation.maxRatioNode = row.nodeId;
            evaluation.maxRatioCarrier = carrier;
        }
        if (sourceQualified && ratio > eps) {
            NewtonCarrierRowConvergenceViolation violation;
            violation.nodeId = row.nodeId;
            violation.carrier = carrier;
            violation.residual = residual;
            violation.scale = scale;
            violation.ratio = ratio;
            violation.flux = flux;
            violation.recombination = recombination;
            violation.impact = impact;
            evaluation.violations.push_back(std::move(violation));
        }
    };

    for (const CoupledDDCarrierTermDiagnostic& row : rows) {
        checkCarrier(row, "electron", row.electronResidual, row.electronFlux,
                     row.electronFluxAbsSum, row.electronRecombination, row.electronImpact);
        checkCarrier(row, "hole", row.holeResidual, row.holeFlux,
                     row.holeFluxAbsSum, row.holeRecombination, row.holeImpact);
    }
    evaluation.satisfied = evaluation.violations.empty();
    return evaluation;
}

NewtonGlobalContinuityClosureEvaluation evaluateGlobalContinuityClosure(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes,
    const NewtonGlobalContinuityClosureConfig& cfg)
{
    NewtonGlobalContinuityClosureEvaluation evaluation;
    evaluation.enabled = cfg.mode != "off";
    evaluation.enforced = cfg.mode == "enforce";
    if (!evaluation.enabled)
        return evaluation;

    std::vector<bool> electronContact(rows.size(), false);
    std::vector<bool> holeContact(rows.size(), false);
    for (Index node : electronContactNodes) {
        if (node >= 0 && static_cast<std::size_t>(node) < rows.size())
            electronContact[static_cast<std::size_t>(node)] = true;
    }
    for (Index node : holeContactNodes) {
        if (node >= 0 && static_cast<std::size_t>(node) < rows.size())
            holeContact[static_cast<std::size_t>(node)] = true;
    }

    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        if (electronContact[i])
            evaluation.electron.contactFlux += row.electronFlux;
        else
            evaluation.electron.integratedSource +=
                row.electronRecombination + row.electronImpact;
        if (holeContact[i])
            evaluation.hole.contactFlux += row.holeFlux;
        else
            evaluation.hole.integratedSource +=
                row.holeRecombination + row.holeImpact;
    }

    auto finishCarrier = [&](NewtonGlobalContinuityCarrierClosure& carrier) {
        carrier.mismatch = carrier.contactFlux - carrier.integratedSource;
        const Real sourceMagnitude = std::abs(carrier.integratedSource);
        carrier.qualified = sourceMagnitude >= cfg.sourceFloor;
        const Real scale = std::max({
            std::abs(carrier.contactFlux), sourceMagnitude, cfg.sourceFloor});
        carrier.ratio = scale > 0.0 ? std::abs(carrier.mismatch) / scale : 0.0;
        if (carrier.qualified && carrier.ratio > cfg.tolerance)
            evaluation.satisfied = false;
    };
    finishCarrier(evaluation.electron);
    finishCarrier(evaluation.hole);
    return evaluation;
}

namespace {

VectorXd continuityRowWeights(
    const CoupledDDAssembler& assembler,
    const VectorXd& state,
    const CoupledDDBoundaryConditions& bcs,
    const NewtonContinuityRowScalingConfig& cfg)
{
    const int N = static_cast<int>(assembler.numNodes());
    VectorXd weights = VectorXd::Ones(3 * N);
    if (!cfg.enabled)
        return weights;

    const auto terms =
        assembler.carrierContinuityEquationTermDiagnostics(state, bcs);
    auto rowWeight = [&](Real fluxAbsSum, Real recombination, Real impact) {
        const Real sourceScale = std::max(std::abs(recombination), std::abs(impact));
        if (sourceScale < cfg.minSourceScale ||
            sourceScale < cfg.fluxFraction * std::abs(fluxAbsSum)) {
            return 1.0;
        }
        const Real scale = std::max(sourceScale, cfg.scaleFloor);
        return std::clamp(1.0 / scale, cfg.minWeight, cfg.maxWeight);
    };

    for (int i = 0; i < N; ++i) {
        const Index node = static_cast<Index>(i);
        if (bcs.phin.find(node) == bcs.phin.end()) {
            const auto& term = terms[node];
            weights(N + i) = rowWeight(
                term.electronFluxAbsSum,
                term.electronRecombination,
                term.electronImpact);
        }
        if (bcs.phip.find(node) == bcs.phip.end()) {
            const auto& term = terms[node];
            weights(2 * N + i) = rowWeight(
                term.holeFluxAbsSum,
                term.holeRecombination,
                term.holeImpact);
        }
    }
    return weights;
}

SparseMatrixd leftScaleRows(const SparseMatrixd& matrix, const VectorXd& weights)
{
    SparseMatrixd scaled = matrix;
    for (int outer = 0; outer < scaled.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(scaled, outer); entry; ++entry)
            entry.valueRef() *= weights(entry.row());
    }
    return scaled;
}

DDScalingSpec buildRecoveryScalingSpec(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       const NewtonConfig& cfg)
{
    DDScalingSpec scaling;
    if (!cfg.inputScaling.isUnitScaling())
        return scaling;

    const Real epsRef = constants::eps0 *
        maxRelativePermittivityAcrossRegions(mesh, matdb, cfg.temperature_K);
    const Real niFloor =
        maxIntrinsicDensityAcrossRegions(mesh, matdb, cfg.temperature_K);
    const UnitScalingSystem::AutoInputs autoInputs =
        UnitScalingSystem::autoInputsFrom(mesh, doping, matdb, niFloor);
    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        cfg.temperature_K, epsRef, autoInputs, cfg.unitScalingRefs, cfg.inputScaling.unitSystem());

    scaling.enabled = true;
    scaling.V0 = sc.V0();
    scaling.C0 = sc.C0();
    scaling.mu0 = sc.mu0();
    scaling.D0 = sc.D0();
    scaling.L0 = sc.L0();
    scaling.permittivityReference_F_per_m = epsRef;
    scaling.unitSystem = cfg.inputScaling.unitSystem();
    scaling.chargeAreaFactor = cfg.inputScaling.unitSystem().chargeAreaFactor();
    scaling.chargeLineFactor = cfg.inputScaling.unitSystem().chargeLineFactor();
    scaling.fieldFromCoordinateDeltaFactor = cfg.inputScaling.unitSystem().fieldFromCoordinateDeltaFactor();
    scaling.currentDensityLineIntegralFactor =
        cfg.inputScaling.unitSystem().currentDensityAM2PerInternal() *
        cfg.inputScaling.unitSystem().lengthMPerInternal();
    return scaling;
}

bool hasContactBiasForRecovery(const std::unordered_map<std::string, Real>& biases,
                               const std::string& name)
{
    if (biases.find(name) != biases.end())
        return true;
    auto lower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return s;
    };
    const std::string target = lower(name);
    for (const auto& [key, value] : biases) {
        (void)value;
        if (lower(key) == target)
            return true;
    }
    return false;
}

} // namespace

NewtonCarrierRowRecoveryResult recoverCarrierRowsWithGummelDensity(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    const NewtonConfig& cfg,
    const DDSolution& state,
    const std::vector<NewtonCarrierRowConvergenceViolation>& violations,
    const NewtonCarrierRowRecoveryConfig& recovery,
    const ContactSpecsMap& contactSpecs)
{
    NewtonCarrierRowRecoveryResult result;
    result.solution = state;
    result.mode = recovery.mode;
    if (recovery.mode != "gummel_density" || violations.empty())
        return result;
    result.attempted = true;
    result.cyclesAttempted = 1;

    const double Vt = thermalVoltage(cfg.temperature_K);
    MobilityModelConfig mobilityConfig = cfg.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg.recombination, cfg.taun, cfg.taup, cfg.srhDopingDependence);
    recombinationConfig.augerCn = cfg.augerCn;
    recombinationConfig.augerCp = cfg.augerCp;
    recombinationConfig.bandToBand = cfg.bandToBand;
    const DDScalingSpec scaling = buildRecoveryScalingSpec(mesh, matdb, doping, cfg);

    DDAssembler densityAssembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg.bandgapNarrowing,
        cfg.impactIonization,
        {},
        {},
        scaling,
        cfg.carrierStatistics);
    CoupledDDAssembler qfAssembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg.bandgapNarrowing,
        cfg.impactIonization,
        {},
        {},
        scaling,
        cfg.carrierDiagonalFloor,
        cfg.carrierStatistics);

    const int nNodes = static_cast<int>(mesh.numNodes());
    if (state.psi.size() != nNodes || state.phin.size() != nNodes ||
        state.phip.size() != nNodes || state.n.size() != nNodes ||
        state.p.size() != nNodes) {
        throw std::invalid_argument(
            "recoverCarrierRowsWithGummelDensity: state size does not match mesh.");
    }

    VectorXd psi = scaling.enabled ? (state.psi / scaling.V0) : state.psi;
    VectorXd nOld = scaling.enabled ? (state.n / scaling.C0) : state.n;
    VectorXd pOld = scaling.enabled ? (state.p / scaling.C0) : state.p;

    std::unordered_map<Index, Real> nBC;
    std::unordered_map<Index, Real> pBC;
    const auto& ni = qfAssembler.intrinsicDensity();
    const auto& Nc = qfAssembler.electronDensityOfStates();
    const auto& Nv = qfAssembler.holeDensityOfStates();
    const bool dominantSignedContactMean =
        cfg.contactBoundaryReconstruction == "dominant_signed_contact_mean";
    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& contact = mesh.getContact(c);
        if (!hasContactBiasForRecovery(contactBiases, contact.name))
            continue;
        const auto specIt = contactSpecs.find(contact.name);
        if (specIt != contactSpecs.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        for (Index node : contact.node_ids) {
            const Real niNode = ni[static_cast<std::size_t>(node)];
            const Real ndop = ohmicContactNetDoping(
                mesh, doping, contact, node, dominantSignedContactMean);
            const EquilibriumCarrierState equilibrium = equilibriumCarrierState(
                ndop, niNode, Nc[static_cast<std::size_t>(node)],
                Nv[static_cast<std::size_t>(node)], Vt, cfg.carrierStatistics);
            nBC[node] = scaling.enabled ? (equilibrium.n / scaling.C0) : equilibrium.n;
            pBC[node] = scaling.enabled ? (equilibrium.p / scaling.C0) : equilibrium.p;
        }
    }

    LinearSolver linearSolver;
    VectorXd nNew = nOld;
    VectorXd pNew = pOld;
    const int densityPasses = std::max(1, recovery.maxAttempts);
    const Real densityChangeReltol = std::max<Real>(0.0, recovery.densityChangeReltol);
    for (int pass = 0; pass < densityPasses; ++pass) {
        const VectorXd nBefore = nNew;
        const VectorXd pBefore = pNew;
        densityAssembler.assembleElectronContinuity(psi, nNew, pNew);
        densityAssembler.applyDirichlet(nBC);
        VectorXd nCandidate = linearSolver.solve(densityAssembler.matrix(), densityAssembler.rhs());
        for (int i = 0; i < nNodes; ++i) {
            if (nBC.find(static_cast<Index>(i)) == nBC.end() && nCandidate(i) <= 0.0)
                nCandidate(i) = std::numeric_limits<Real>::min();
        }

        densityAssembler.assembleHoleContinuity(psi, nCandidate, pNew);
        densityAssembler.applyDirichlet(pBC);
        VectorXd pCandidate = linearSolver.solve(densityAssembler.matrix(), densityAssembler.rhs());
        for (int i = 0; i < nNodes; ++i) {
            if (pBC.find(static_cast<Index>(i)) == pBC.end() && pCandidate(i) <= 0.0)
                pCandidate(i) = std::numeric_limits<Real>::min();
        }

        Real maxRelChange = 0.0;
        for (int i = 0; i < nNodes; ++i) {
            const Real nDenom = std::max({std::abs(nBefore(i)), std::abs(nCandidate(i)), std::numeric_limits<Real>::min()});
            const Real pDenom = std::max({std::abs(pBefore(i)), std::abs(pCandidate(i)), std::numeric_limits<Real>::min()});
            maxRelChange = std::max(maxRelChange, std::abs(nCandidate(i) - nBefore(i)) / nDenom);
            maxRelChange = std::max(maxRelChange, std::abs(pCandidate(i) - pBefore(i)) / pDenom);
        }
        nNew = std::move(nCandidate);
        pNew = std::move(pCandidate);
        result.maxDensityRelativeChange = maxRelChange;
        result.densityPasses = pass + 1;
        if (maxRelChange <= densityChangeReltol) {
            result.densityConverged = true;
            break;
        }
    }

    bool recoverElectrons = false;
    bool recoverHoles = false;
    for (const auto& violation : violations) {
        if (violation.nodeId >= mesh.numNodes())
            continue;
        if (violation.carrier == "electron")
            recoverElectrons = true;
        else if (violation.carrier == "hole")
            recoverHoles = true;
    }

    for (int i = 0; i < nNodes; ++i) {
        const Index node = static_cast<Index>(i);
        const Real psiSi = state.psi(i);
        const Real niNode = ni[static_cast<std::size_t>(i)];
        if (recoverElectrons && nBC.find(node) == nBC.end() && niNode > 0.0) {
            const Real recoveredN = scaling.enabled ? (nNew(i) * scaling.C0) : nNew(i);
            if (recoveredN > 0.0 && std::isfinite(recoveredN)) {
                const Real before = std::max(std::abs(state.n(i)), std::numeric_limits<Real>::min());
                result.solution.n(i) = recoveredN;
                result.solution.phin(i) = electronQuasiFermiPotential(
                    niNode, Nc[static_cast<std::size_t>(i)], psiSi,
                    recoveredN, Vt, cfg.carrierStatistics);
                result.maxCarrierDensityRatio = std::max(
                    result.maxCarrierDensityRatio, recoveredN / before);
                ++result.electronRowsUpdated;
            }
        }
        if (recoverHoles && pBC.find(node) == pBC.end() && niNode > 0.0) {
            const Real recoveredP = scaling.enabled ? (pNew(i) * scaling.C0) : pNew(i);
            if (recoveredP > 0.0 && std::isfinite(recoveredP)) {
                const Real before = std::max(std::abs(state.p(i)), std::numeric_limits<Real>::min());
                result.solution.p(i) = recoveredP;
                result.solution.phip(i) = holeQuasiFermiPotential(
                    niNode, Nv[static_cast<std::size_t>(i)], psiSi,
                    recoveredP, Vt, cfg.carrierStatistics);
                result.maxCarrierDensityRatio = std::max(
                    result.maxCarrierDensityRatio, recoveredP / before);
                ++result.holeRowsUpdated;
            }
        }
    }
    result.maxPsiDelta_V = 0.0;
    return result;
}
NewtonConfig newtonConfigFromJson(const nlohmann::json& json, UnitScalingConfig scaling)
{
    NewtonConfig cfg;
    cfg.inputScaling = scaling;
    cfg.maxIter = json.value("max_iter", cfg.maxIter);
    cfg.reltol = json.value("reltol", cfg.reltol);
    cfg.abstol = json.value("abstol", cfg.abstol);
    cfg.temperature_K = json.value("temperature_K", cfg.temperature_K);
    if (json.contains("electron_quantum_potential")) {
        const auto& value = json.at("electron_quantum_potential");
        if (value.is_boolean()) {
            cfg.electronQuantumPotential.enabled = value.get<bool>();
        } else if (value.is_object()) {
            cfg.electronQuantumPotential.enabled = value.value(
                "enabled", cfg.electronQuantumPotential.enabled);
            cfg.electronQuantumPotential.couplingMode = value.value(
                "coupling_mode", cfg.electronQuantumPotential.couplingMode);
            cfg.electronQuantumPotential.formulation = value.value(
                "formulation", cfg.electronQuantumPotential.formulation);
            cfg.electronQuantumPotential.gamma = value.value(
                "gamma", cfg.electronQuantumPotential.gamma);
            cfg.electronQuantumPotential.effectiveMassRatio = value.value(
                "effective_mass_ratio", cfg.electronQuantumPotential.effectiveMassRatio);
            if (value.contains("coefficient_mass_ratio")) {
                cfg.electronQuantumPotential.coefficientMassRatio =
                    value.at("coefficient_mass_ratio").get<Real>();
            }
            cfg.electronQuantumPotential.includeInsulators = value.value(
                "include_insulators", cfg.electronQuantumPotential.includeInsulators);
            cfg.electronQuantumPotential.insulatorGamma = value.value(
                "insulator_gamma", cfg.electronQuantumPotential.insulatorGamma);
            cfg.electronQuantumPotential.insulatorEffectiveMassRatio = value.value(
                "insulator_effective_mass_ratio",
                cfg.electronQuantumPotential.insulatorEffectiveMassRatio);
            if (value.contains("insulator_coefficient_mass_ratio")) {
                cfg.electronQuantumPotential.insulatorCoefficientMassRatio =
                    value.at("insulator_coefficient_mass_ratio").get<Real>();
            }
            cfg.electronQuantumPotential.interfaceBoundary = value.value(
                "interface_boundary", cfg.electronQuantumPotential.interfaceBoundary);
            cfg.electronQuantumPotential.theta = value.value(
                "theta", cfg.electronQuantumPotential.theta);
            cfg.electronQuantumPotential.conductionBandNarrowingFraction = value.value(
                "conduction_band_narrowing_fraction",
                cfg.electronQuantumPotential.conductionBandNarrowingFraction);
            cfg.electronQuantumPotential.maxIterations = value.value(
                "max_iterations", cfg.electronQuantumPotential.maxIterations);
            cfg.electronQuantumPotential.outerMaxIterations = value.value(
                "outer_max_iterations", cfg.electronQuantumPotential.outerMaxIterations);
            cfg.electronQuantumPotential.relativeTolerance = value.value(
                "relative_tolerance", cfg.electronQuantumPotential.relativeTolerance);
            cfg.electronQuantumPotential.absoluteTolerance_V = value.value(
                "absolute_tolerance_V", cfg.electronQuantumPotential.absoluteTolerance_V);
            if (value.contains("outer_absolute_tolerance_V")) {
                cfg.electronQuantumPotential.outerAbsoluteTolerance_V =
                    value.at("outer_absolute_tolerance_V").get<Real>();
            }
            cfg.electronQuantumPotential.damping = value.value(
                "damping", cfg.electronQuantumPotential.damping);
            cfg.electronQuantumPotential.maxUpdate_V = value.value(
                "max_update_V", cfg.electronQuantumPotential.maxUpdate_V);
            cfg.electronQuantumPotential.outerAcceleration = value.value(
                "outer_acceleration", cfg.electronQuantumPotential.outerAcceleration);
            cfg.electronQuantumPotential.outerRelaxation = value.value(
                "outer_relaxation", cfg.electronQuantumPotential.outerRelaxation);
            cfg.electronQuantumPotential.outerRelaxationMin = value.value(
                "outer_relaxation_min", cfg.electronQuantumPotential.outerRelaxationMin);
            cfg.electronQuantumPotential.outerRelaxationMax = value.value(
                "outer_relaxation_max", cfg.electronQuantumPotential.outerRelaxationMax);
            cfg.electronQuantumPotential.residualDiagnosticPrefix = value.value(
                "residual_diagnostic_prefix",
                cfg.electronQuantumPotential.residualDiagnosticPrefix);
            cfg.electronQuantumPotential.residualDiagnosticUseInitialState = value.value(
                "residual_diagnostic_use_initial_state",
                cfg.electronQuantumPotential.residualDiagnosticUseInitialState);
            cfg.electronQuantumPotential.globalDiscretization = value.value(
                "global_discretization",
                cfg.electronQuantumPotential.globalDiscretization);
            cfg.electronQuantumPotential.sentaurusInterfaceInsulatorHalfJumpOffset =
                value.value(
                    "sentaurus_interface_insulator_half_jump_offset",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceInsulatorHalfJumpOffset);
            cfg.electronQuantumPotential.sentaurusInterfaceSiliconHalfJumpOffset =
                value.value("sentaurus_interface_silicon_half_jump_offset",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceSiliconHalfJumpOffset);
            cfg.electronQuantumPotential.sentaurusInterfacePolysiliconHalfJumpOffset =
                value.value("sentaurus_interface_polysilicon_half_jump_offset",
                    cfg.electronQuantumPotential.
                        sentaurusInterfacePolysiliconHalfJumpOffset);
            cfg.electronQuantumPotential.sentaurusInterfaceSiliconReactionWeight =
                value.value("sentaurus_interface_silicon_reaction_weight",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceSiliconReactionWeight);
            cfg.electronQuantumPotential.sentaurusInterfacePolysiliconReactionWeight =
                value.value("sentaurus_interface_polysilicon_reaction_weight",
                    cfg.electronQuantumPotential.
                        sentaurusInterfacePolysiliconReactionWeight);
            cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtSiliconReactionWeight =
                value.value(
                    "sentaurus_interface_insulator_at_silicon_reaction_weight",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceInsulatorAtSiliconReactionWeight);
            cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtPolysiliconReactionWeight =
                value.value(
                    "sentaurus_interface_insulator_at_polysilicon_reaction_weight",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceInsulatorAtPolysiliconReactionWeight);
            cfg.electronQuantumPotential.sentaurusInterfaceSiliconReactionOffset_V =
                value.value("sentaurus_interface_silicon_reaction_offset_V",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceSiliconReactionOffset_V);
            cfg.electronQuantumPotential.
                sentaurusInterfacePolysiliconReactionOffset_V =
                value.value("sentaurus_interface_polysilicon_reaction_offset_V",
                    cfg.electronQuantumPotential.
                        sentaurusInterfacePolysiliconReactionOffset_V);
            cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtSiliconReactionOffset_V =
                value.value(
                    "sentaurus_interface_insulator_at_silicon_reaction_offset_V",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceInsulatorAtSiliconReactionOffset_V);
            cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V =
                value.value(
                    "sentaurus_interface_insulator_at_polysilicon_reaction_offset_V",
                    cfg.electronQuantumPotential.
                        sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V);
            cfg.electronQuantumPotential.
                sentaurusInsulatorReentrantCornerReactionWeight =
                value.value(
                    "sentaurus_insulator_reentrant_corner_reaction_weight",
                    cfg.electronQuantumPotential.
                        sentaurusInsulatorReentrantCornerReactionWeight);
            cfg.electronQuantumPotential.oxideBoundary = value.value(
                "oxide_boundary", cfg.electronQuantumPotential.oxideBoundary);
            cfg.electronQuantumPotential.oxideQuantumMassRatio = value.value(
                "oxide_quantum_mass_ratio",
                cfg.electronQuantumPotential.oxideQuantumMassRatio);
            cfg.electronQuantumPotential.oxideBarrierMassRatio = value.value(
                "oxide_barrier_mass_ratio",
                cfg.electronQuantumPotential.oxideBarrierMassRatio);
            cfg.electronQuantumPotential.oxideBarrierHeight_V = value.value(
                "oxide_barrier_height_V",
                cfg.electronQuantumPotential.oxideBarrierHeight_V);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: electron_quantum_potential must be a boolean or object.");
        }
    }
    if (json.contains("carrier_statistics")) {
        const auto& statistics = json.at("carrier_statistics");
        if (statistics.is_string())
            cfg.carrierStatistics.model = statistics.get<std::string>();
        else if (statistics.is_object())
            cfg.carrierStatistics.model = statistics.value("model", cfg.carrierStatistics.model);
        else
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_statistics must be a string or object.");
    }
    if (json.contains("statistics"))
        cfg.carrierStatistics.model = json.at("statistics").get<std::string>();
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
    cfg.poissonLineSearchStallResidualFloor = json.value(
        "poisson_line_search_stall_residual_floor",
        cfg.poissonLineSearchStallResidualFloor);
    cfg.poissonLineSearchStallRelativeIncrease = json.value(
        "poisson_line_search_stall_relative_increase",
        cfg.poissonLineSearchStallRelativeIncrease);
    cfg.poissonLineSearchStallCarrierResidualFloor = json.value(
        "poisson_line_search_stall_carrier_residual_floor",
        cfg.poissonLineSearchStallCarrierResidualFloor);
    cfg.poissonLineSearchStallContactMajorityQfDropLimit_V = json.value(
        "poisson_line_search_stall_contact_majority_qf_drop_limit_V",
        cfg.poissonLineSearchStallContactMajorityQfDropLimit_V);
    cfg.carrierRowQualifiedStallAcceptance = json.value(
        "carrier_row_qualified_stall_acceptance",
        cfg.carrierRowQualifiedStallAcceptance);
    cfg.carrierRegularizationScale = json.value(
        "carrier_regularization_scale",
        cfg.carrierRegularizationScale);
    auto parseCarrierDiagonalFloor = [&](const nlohmann::json& value) {
        if (value.is_boolean()) {
            cfg.carrierDiagonalFloor.enabled = value.get<bool>();
            return;
        }
        if (!value.is_object()) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_diagonal_floor must be a boolean or object.");
        }
        cfg.carrierDiagonalFloor.enabled = value.value(
            "enabled", cfg.carrierDiagonalFloor.enabled);
        cfg.carrierDiagonalFloor.scale = value.value(
            "scale", cfg.carrierDiagonalFloor.scale);
        cfg.carrierDiagonalFloor.minorityDensityRatio = value.value(
            "minority_density_ratio", cfg.carrierDiagonalFloor.minorityDensityRatio);
    };
    if (json.contains("carrier_diagonal_floor_regularization"))
        parseCarrierDiagonalFloor(json.at("carrier_diagonal_floor_regularization"));
    if (json.contains("carrier_diagonal_floor"))
        parseCarrierDiagonalFloor(json.at("carrier_diagonal_floor"));
    if (json.contains("carrier_row_convergence")) {
        const auto& value = json.at("carrier_row_convergence");
        if (value.is_boolean()) {
            cfg.carrierRowConvergence.mode = value.get<bool>() ? "report" : "off";
        } else if (value.is_object()) {
            cfg.carrierRowConvergence.mode = value.value(
                "mode", cfg.carrierRowConvergence.mode);
            if (value.contains("enabled") && !value.at("enabled").get<bool>())
                cfg.carrierRowConvergence.mode = "off";
            cfg.carrierRowConvergence.epsRow = value.value(
                "eps_row", cfg.carrierRowConvergence.epsRow);
            cfg.carrierRowConvergence.scaleFloor = value.value(
                "scale_floor", cfg.carrierRowConvergence.scaleFloor);
            cfg.carrierRowConvergence.minSourceScaleFraction = value.value(
                "min_source_scale_fraction", cfg.carrierRowConvergence.minSourceScaleFraction);
            cfg.carrierRowConvergence.minSourceScale = value.value(
                "min_source_scale", cfg.carrierRowConvergence.minSourceScale);
            cfg.carrierRowConvergence.minEnforceMaxIter = value.value(
                "min_newton_max_iter", cfg.carrierRowConvergence.minEnforceMaxIter);
            cfg.carrierRowConvergence.diagnosticCsvFile = value.value(
                "diagnostic_csv", cfg.carrierRowConvergence.diagnosticCsvFile);
            cfg.carrierRowConvergence.traceCsvFile = value.value(
                "trace_csv", cfg.carrierRowConvergence.traceCsvFile);
            cfg.carrierRowConvergence.traceNodes = value.value(
                "trace_nodes", cfg.carrierRowConvergence.traceNodes);
            cfg.carrierRowConvergence.traceFirstIterations = value.value(
                "trace_first_iterations", cfg.carrierRowConvergence.traceFirstIterations);
            cfg.carrierRowConvergence.traceEveryIterations = value.value(
                "trace_every_iterations", cfg.carrierRowConvergence.traceEveryIterations);
            if (value.contains("recovery")) {
                const auto& recovery = value.at("recovery");
                if (recovery.is_string()) {
                    cfg.carrierRowRecovery.mode = recovery.get<std::string>();
                } else if (recovery.is_object()) {
                    cfg.carrierRowRecovery.mode = recovery.value(
                        "mode", cfg.carrierRowRecovery.mode);
                    cfg.carrierRowRecovery.maxAttempts = recovery.value(
                        "max_attempts", cfg.carrierRowRecovery.maxAttempts);
                    cfg.carrierRowRecovery.maxCycles = recovery.value(
                        "max_cycles", cfg.carrierRowRecovery.maxCycles);
                    cfg.carrierRowRecovery.densityChangeReltol = recovery.value(
                        "density_change_reltol", cfg.carrierRowRecovery.densityChangeReltol);
                } else if (recovery.is_boolean()) {
                    cfg.carrierRowRecovery.mode = recovery.get<bool>() ? "gummel_density" : "off";
                } else {
                    throw std::invalid_argument(
                        "newtonConfigFromJson: carrier_row_convergence.recovery must be a string, boolean, or object.");
                }
            }
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence must be a boolean or object.");
        }
        if (cfg.carrierRowConvergence.mode != "off" &&
            cfg.carrierRowConvergence.mode != "report" &&
            cfg.carrierRowConvergence.mode != "enforce") {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.mode must be 'off', 'report', or 'enforce'.");
        }
        if (!(cfg.carrierRowConvergence.epsRow > 0.0) ||
            !std::isfinite(cfg.carrierRowConvergence.epsRow)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.eps_row must be positive and finite.");
        }
        if (cfg.carrierRowConvergence.scaleFloor < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.scaleFloor)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.scale_floor must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.minSourceScaleFraction < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minSourceScaleFraction)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_source_scale_fraction must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.minSourceScale < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minSourceScale)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_source_scale must be finite and nonnegative.");
        }
        if (cfg.carrierRowRecovery.mode != "off" &&
            cfg.carrierRowRecovery.mode != "gummel_density") {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.recovery.mode must be 'off' or 'gummel_density'.");
        }
        if (cfg.carrierRowRecovery.maxAttempts < 0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.recovery.max_attempts must be nonnegative.");
        }
        if (cfg.carrierRowRecovery.maxCycles < 0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.recovery.max_cycles must be nonnegative.");
        }
        if (cfg.carrierRowRecovery.densityChangeReltol < 0.0 ||
            !std::isfinite(cfg.carrierRowRecovery.densityChangeReltol)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.recovery.density_change_reltol must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.mode == "enforce" &&
            cfg.maxIter < cfg.carrierRowConvergence.minEnforceMaxIter) {
            cfg.maxIter = cfg.carrierRowConvergence.minEnforceMaxIter;
        }
    }
    if (json.contains("continuity_row_scaling")) {
        const auto& value = json.at("continuity_row_scaling");
        if (value.is_boolean()) {
            cfg.continuityRowScaling.enabled = value.get<bool>();
        } else if (value.is_object()) {
            cfg.continuityRowScaling.enabled =
                value.value("enabled", cfg.continuityRowScaling.enabled);
            cfg.continuityRowScaling.fluxFraction =
                value.value("flux_fraction", cfg.continuityRowScaling.fluxFraction);
            cfg.continuityRowScaling.scaleFloor =
                value.value("scale_floor", cfg.continuityRowScaling.scaleFloor);
            cfg.continuityRowScaling.minSourceScale =
                value.value("min_source_scale", cfg.continuityRowScaling.minSourceScale);
            cfg.continuityRowScaling.minWeight =
                value.value("min_weight", cfg.continuityRowScaling.minWeight);
            cfg.continuityRowScaling.maxWeight =
                value.value("max_weight", cfg.continuityRowScaling.maxWeight);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: continuity_row_scaling must be a boolean or object.");
        }
        if (cfg.continuityRowScaling.fluxFraction < 0.0 ||
            !std::isfinite(cfg.continuityRowScaling.fluxFraction) ||
            cfg.continuityRowScaling.scaleFloor <= 0.0 ||
            !std::isfinite(cfg.continuityRowScaling.scaleFloor) ||
            cfg.continuityRowScaling.minSourceScale < 0.0 ||
            !std::isfinite(cfg.continuityRowScaling.minSourceScale) ||
            cfg.continuityRowScaling.minWeight <= 0.0 ||
            !std::isfinite(cfg.continuityRowScaling.minWeight) ||
            cfg.continuityRowScaling.maxWeight < cfg.continuityRowScaling.minWeight ||
            !std::isfinite(cfg.continuityRowScaling.maxWeight)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: invalid continuity_row_scaling bounds.");
        }
    }
    if (json.contains("global_continuity_closure")) {
        const auto& value = json.at("global_continuity_closure");
        if (value.is_boolean()) {
            cfg.globalContinuityClosure.mode = value.get<bool>() ? "report" : "off";
        } else if (value.is_object()) {
            cfg.globalContinuityClosure.mode =
                value.value("mode", cfg.globalContinuityClosure.mode);
            if (value.contains("enabled") && !value.at("enabled").get<bool>())
                cfg.globalContinuityClosure.mode = "off";
            cfg.globalContinuityClosure.tolerance =
                value.value("tolerance", cfg.globalContinuityClosure.tolerance);
            cfg.globalContinuityClosure.sourceFloor =
                value.value("source_floor", cfg.globalContinuityClosure.sourceFloor);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: global_continuity_closure must be a boolean or object.");
        }
        if (cfg.globalContinuityClosure.mode != "off" &&
            cfg.globalContinuityClosure.mode != "report" &&
            cfg.globalContinuityClosure.mode != "enforce") {
            throw std::invalid_argument(
                "newtonConfigFromJson: global_continuity_closure.mode must be 'off', 'report', or 'enforce'.");
        }
        if (!(cfg.globalContinuityClosure.tolerance > 0.0) ||
            !std::isfinite(cfg.globalContinuityClosure.tolerance) ||
            cfg.globalContinuityClosure.sourceFloor < 0.0 ||
            !std::isfinite(cfg.globalContinuityClosure.sourceFloor)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: invalid global_continuity_closure bounds.");
        }
    }
    cfg.finiteDifferenceStep = json.value("finite_difference_step", cfg.finiteDifferenceStep);
    cfg.jacobian = json.value("jacobian", cfg.jacobian);
    cfg.quasiFermiReference =
        json.value("quasi_fermi_reference", cfg.quasiFermiReference);
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
    if (json.contains("srh_doping_dependence")) {
        cfg.srhDopingDependence = srhDopingDependenceConfigFromJson(
            json.at("srh_doping_dependence"), scaling);
    }
    cfg.srhDopingDependence.temperature_K = cfg.temperature_K;
    // The compiled Auger defaults are SI literals; express them in the active
    // internal unit before applying any deck override, which is already
    // written in internal units.
    cfg.augerCn = scaling.unitSystem().m6PerSToInternalAugerCoefficient(cfg.augerCn);
    cfg.augerCp = scaling.unitSystem().m6PerSToInternalAugerCoefficient(cfg.augerCp);
    cfg.augerCn = scaling.augerCoefficientToInternal(
        json.value("auger_cn_m6_per_s", cfg.augerCn));
    cfg.augerCp = scaling.augerCoefficientToInternal(
        json.value("auger_cp_m6_per_s", cfg.augerCp));
    if (json.contains("mobility"))
        cfg.mobility = mobilityModelConfigFromJson(json.at("mobility"), scaling);
    if (json.contains("bandgap_narrowing")) {
        const auto& value = json.at("bandgap_narrowing");
        if (value.is_string()) {
            cfg.bandgapNarrowing =
                bandgapNarrowingConfig(value.get<std::string>(), scaling);
        } else if (value.is_object()) {
            cfg.bandgapNarrowing = bandgapNarrowingConfig(
                value.value("model", cfg.bandgapNarrowing.model), scaling);
            if (value.contains("reference_doping_m3")) {
                cfg.bandgapNarrowing.referenceDoping = scaling.concentrationToInternal(
                    value.at("reference_doping_m3").get<Real>());
            }
            cfg.bandgapNarrowing.coefficient = value.value(
                "coefficient_eV", cfg.bandgapNarrowing.coefficient);
            cfg.bandgapNarrowing.smoothing = value.value(
                "smoothing", cfg.bandgapNarrowing.smoothing);
            cfg.bandgapNarrowing.offset = value.value(
                "offset_eV", cfg.bandgapNarrowing.offset);
            cfg.bandgapNarrowing.fermiStatisticsCorrection = value.value(
                "fermi_statistics_correction",
                cfg.bandgapNarrowing.fermiStatisticsCorrection);
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
    if (json.contains("band_to_band")) {
        cfg.bandToBand = bandToBandTunnelingConfigFromJson(
            json.at("band_to_band"), "newtonConfigFromJson");
    }

    if (json.contains("impact_ionization")) {
        const auto& value = json.at("impact_ionization");
        if (value.is_string()) {
            cfg.impactIonization = impactIonizationModelConfig(
                value.get<std::string>(), scaling);
        } else if (value.is_object()) {
            cfg.impactIonization = impactIonizationModelConfig(
                value.value("model", cfg.impactIonization.model), scaling);
            cfg.impactIonization.couplingMode = value.value(
                "coupling_mode", cfg.impactIonization.couplingMode);
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
            cfg.impactIonization.eparallelFieldRecovery = value.value(
                "eparallel_field_recovery", cfg.impactIonization.eparallelFieldRecovery);
            cfg.impactIonization.cellReconstructedMidpointDensity = value.value(
                "cell_reconstructed_midpoint_density",
                cfg.impactIonization.cellReconstructedMidpointDensity);
            cfg.impactIonization.quasiFermiGradientDiscretization = value.value(
                "quasi_fermi_gradient_discretization",
                cfg.impactIonization.quasiFermiGradientDiscretization);
            cfg.impactIonization.contactElectricFieldFallback = value.value(
                "contact_electric_field_fallback",
                cfg.impactIonization.contactElectricFieldFallback);
            cfg.impactIonization.contactElectricFieldFallbackScope = value.value(
                "contact_electric_field_fallback_scope",
                cfg.impactIonization.contactElectricFieldFallbackScope);
            cfg.impactIonization.contactElectricFieldFallbackMode = value.value(
                "contact_electric_field_fallback_mode",
                cfg.impactIonization.contactElectricFieldFallbackMode);
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
            cfg.impactIonization.edgeSourcePartition = value.value(
                "edge_source_partition", cfg.impactIonization.edgeSourcePartition);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_truncation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_trucation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.minimumField = scaling.electricFieldToInternal(value.value(
                "minimum_field_V_m", cfg.impactIonization.minimumField));
            cfg.impactIonization.debugRawVanOverstraeten = value.value(
                "debug_raw_vanoverstraeten",
                cfg.impactIonization.debugRawVanOverstraeten);
            cfg.impactIonization.aScale = value.value(
                "A_scale", cfg.impactIonization.aScale);
            cfg.impactIonization.bScale = value.value(
                "B_scale", cfg.impactIonization.bScale);
            if (value.contains("electron_A_m_inv")) {
                cfg.impactIonization.electronA = scaling.inverseLengthToInternal(
                    value.at("electron_A_m_inv").get<Real>());
            }
            if (value.contains("electron_B_V_m")) {
                cfg.impactIonization.electronB = scaling.electricFieldToInternal(
                    value.at("electron_B_V_m").get<Real>());
            }
            if (value.contains("hole_A_m_inv")) {
                cfg.impactIonization.holeA = scaling.inverseLengthToInternal(
                    value.at("hole_A_m_inv").get<Real>());
            }
            if (value.contains("hole_B_V_m")) {
                cfg.impactIonization.holeB = scaling.electricFieldToInternal(
                    value.at("hole_B_V_m").get<Real>());
            }
            if (value.contains("electron_a_low_m_inv")) {
                cfg.impactIonization.electronALow = scaling.inverseLengthToInternal(
                    value.at("electron_a_low_m_inv").get<Real>());
            }
            if (value.contains("electron_a_high_m_inv")) {
                cfg.impactIonization.electronAHigh = scaling.inverseLengthToInternal(
                    value.at("electron_a_high_m_inv").get<Real>());
            }
            if (value.contains("electron_b_low_V_m")) {
                cfg.impactIonization.electronBLow = scaling.electricFieldToInternal(
                    value.at("electron_b_low_V_m").get<Real>());
            }
            if (value.contains("electron_b_high_V_m")) {
                cfg.impactIonization.electronBHigh = scaling.electricFieldToInternal(
                    value.at("electron_b_high_V_m").get<Real>());
            }
            if (value.contains("hole_a_low_m_inv")) {
                cfg.impactIonization.holeALow = scaling.inverseLengthToInternal(
                    value.at("hole_a_low_m_inv").get<Real>());
            }
            if (value.contains("hole_a_high_m_inv")) {
                cfg.impactIonization.holeAHigh = scaling.inverseLengthToInternal(
                    value.at("hole_a_high_m_inv").get<Real>());
            }
            if (value.contains("hole_b_low_V_m")) {
                cfg.impactIonization.holeBLow = scaling.electricFieldToInternal(
                    value.at("hole_b_low_V_m").get<Real>());
            }
            if (value.contains("hole_b_high_V_m")) {
                cfg.impactIonization.holeBHigh = scaling.electricFieldToInternal(
                    value.at("hole_b_high_V_m").get<Real>());
            }
            if (value.contains("switch_field_V_m")) {
                cfg.impactIonization.switchField = scaling.electricFieldToInternal(
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
    (void)BandToBandTunnelingModel(cfg.bandToBand);

    if (cfg.jacobian != "analytic" && cfg.jacobian != "finite_difference")
        throw std::invalid_argument(
            "newtonConfigFromJson: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg.quasiFermiReference != "none" &&
        cfg.quasiFermiReference != "contact_majority") {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_reference must be 'none' or "
            "'contact_majority'.");
    }
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
    if (cfg.poissonLineSearchStallResidualFloor < 0.0 ||
        !std::isfinite(cfg.poissonLineSearchStallResidualFloor))
        throw std::invalid_argument(
            "newtonConfigFromJson: poisson_line_search_stall_residual_floor must be non-negative and finite.");
    if (cfg.poissonLineSearchStallRelativeIncrease < 0.0 ||
        !std::isfinite(cfg.poissonLineSearchStallRelativeIncrease))
        throw std::invalid_argument(
            "newtonConfigFromJson: poisson_line_search_stall_relative_increase must be non-negative and finite.");
    if (cfg.poissonLineSearchStallCarrierResidualFloor < 0.0 ||
        !std::isfinite(cfg.poissonLineSearchStallCarrierResidualFloor))
        throw std::invalid_argument(
            "newtonConfigFromJson: poisson_line_search_stall_carrier_residual_floor must be non-negative and finite.");
    if (cfg.poissonLineSearchStallContactMajorityQfDropLimit_V < 0.0 ||
        !std::isfinite(cfg.poissonLineSearchStallContactMajorityQfDropLimit_V))
        throw std::invalid_argument(
            "newtonConfigFromJson: poisson_line_search_stall_contact_majority_qf_drop_limit_V must be non-negative and finite.");
    if (cfg.carrierRegularizationScale < 0.0 || !std::isfinite(cfg.carrierRegularizationScale))
        throw std::invalid_argument(
            "newtonConfigFromJson: carrier_regularization_scale must be non-negative and finite.");
    if (cfg.carrierDiagonalFloor.scale < 0.0 ||
        !std::isfinite(cfg.carrierDiagonalFloor.scale))
        throw std::invalid_argument(
            "newtonConfigFromJson: carrier_diagonal_floor scale must be non-negative and finite.");
    if (cfg.carrierDiagonalFloor.minorityDensityRatio < 0.0 ||
        !std::isfinite(cfg.carrierDiagonalFloor.minorityDensityRatio))
        throw std::invalid_argument(
            "newtonConfigFromJson: carrier_diagonal_floor minority_density_ratio must be non-negative and finite.");
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
    if (cfg.electronQuantumPotential.couplingMode != "outer" &&
        cfg.electronQuantumPotential.couplingMode != "frozen") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.coupling_mode "
            "must be 'outer' or 'frozen'.");
    }
    if (cfg.electronQuantumPotential.formulation != "potential_based" &&
        cfg.electronQuantumPotential.formulation != "density_based") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.formulation "
            "must be 'potential_based' or 'density_based'.");
    }
    if (!(cfg.electronQuantumPotential.maxUpdate_V >= 0.0) ||
        !std::isfinite(cfg.electronQuantumPotential.maxUpdate_V)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.max_update_V "
            "must be finite and non-negative.");
    }
    (void)densityGradientCoefficientVm2(
        cfg.electronQuantumPotential.gamma,
        cfg.electronQuantumPotential.coefficientMassRatio.value_or(
            cfg.electronQuantumPotential.effectiveMassRatio));
    (void)densityGradientCoefficientVm2(
        cfg.electronQuantumPotential.insulatorGamma,
        cfg.electronQuantumPotential.insulatorCoefficientMassRatio.value_or(
            cfg.electronQuantumPotential.insulatorEffectiveMassRatio));
    if (cfg.electronQuantumPotential.interfaceBoundary != "homogeneous_neumann" &&
        cfg.electronQuantumPotential.interfaceBoundary != "sentaurus_step") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.interface_boundary "
            "must be 'homogeneous_neumann' or 'sentaurus_step'.");
    }
    if (cfg.electronQuantumPotential.globalDiscretization !=
            "exponential_fitted" &&
        cfg.electronQuantumPotential.globalDiscretization != "p1_direct" &&
        cfg.electronQuantumPotential.globalDiscretization != "cvfem_full" &&
        cfg.electronQuantumPotential.globalDiscretization !=
            "p1_lambda_direct" &&
        cfg.electronQuantumPotential.globalDiscretization !=
            "gss_potentiallike_fitted" &&
        cfg.electronQuantumPotential.globalDiscretization !=
            "sentaurus_box" &&
        cfg.electronQuantumPotential.globalDiscretization !=
            "conservative_sqrt_fitted" &&
        cfg.electronQuantumPotential.globalDiscretization !=
            "gss_density_fitted") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "global_discretization must be 'exponential_fitted', "
            "'p1_direct', 'cvfem_full', 'p1_lambda_direct', "
            "'gss_potentiallike_fitted', 'sentaurus_box', "
            "'conservative_sqrt_fitted', or "
            "'gss_density_fitted'.");
    }
    for (const Real offset : {
             cfg.electronQuantumPotential.
                 sentaurusInterfaceInsulatorHalfJumpOffset,
             cfg.electronQuantumPotential.sentaurusInterfaceSiliconHalfJumpOffset,
             cfg.electronQuantumPotential.
                 sentaurusInterfacePolysiliconHalfJumpOffset}) {
        if (!std::isfinite(offset) || std::abs(offset) > 1.0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: sentaurus_box interface half-jump "
                "offsets must be finite and have absolute value at most one.");
        }
    }
    for (const Real weight : {
             cfg.electronQuantumPotential.sentaurusInterfaceSiliconReactionWeight,
             cfg.electronQuantumPotential.sentaurusInterfacePolysiliconReactionWeight,
             cfg.electronQuantumPotential.
                 sentaurusInterfaceInsulatorAtSiliconReactionWeight,
             cfg.electronQuantumPotential.
                 sentaurusInterfaceInsulatorAtPolysiliconReactionWeight}) {
        if (!(weight > 0.0) || !std::isfinite(weight)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: sentaurus_box interface reaction "
                "weights must be finite and positive.");
        }
    }
    for (const Real offset : {
             cfg.electronQuantumPotential.
                 sentaurusInterfaceSiliconReactionOffset_V,
             cfg.electronQuantumPotential.
                 sentaurusInterfacePolysiliconReactionOffset_V,
             cfg.electronQuantumPotential.
                 sentaurusInterfaceInsulatorAtSiliconReactionOffset_V,
             cfg.electronQuantumPotential.
                 sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V}) {
        if (!std::isfinite(offset) || std::abs(offset) > 1.0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: sentaurus_box interface reaction "
                "offsets must be finite and have absolute value at most one V.");
        }
    }
    const Real insulatorCornerWeight = cfg.electronQuantumPotential.
        sentaurusInsulatorReentrantCornerReactionWeight;
    if (!(insulatorCornerWeight > 0.0) ||
        !std::isfinite(insulatorCornerWeight)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: sentaurus_box insulator re-entrant "
            "corner reaction weight must be finite and positive.");
    }
    if (cfg.electronQuantumPotential.oxideBoundary != "none" &&
        cfg.electronQuantumPotential.oxideBoundary != "devsim_wkb") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "oxide_boundary must be 'none' or 'devsim_wkb'.");
    }
    if (cfg.electronQuantumPotential.includeInsulators &&
        cfg.electronQuantumPotential.oxideBoundary != "none") {
        throw std::invalid_argument(
            "newtonConfigFromJson: explicit insulator density-gradient solve "
            "cannot be combined with an oxide WKB truncation boundary.");
    }
    if (!(cfg.electronQuantumPotential.oxideQuantumMassRatio > 0.0) ||
        !(cfg.electronQuantumPotential.oxideBarrierMassRatio > 0.0) ||
        !(cfg.electronQuantumPotential.oxideBarrierHeight_V > 0.0)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential oxide WKB "
            "masses and barrier height must be positive.");
    }
    if (!(cfg.electronQuantumPotential.theta > 0.0) ||
        !std::isfinite(cfg.electronQuantumPotential.theta)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.theta must be "
            "finite and positive.");
    }
    if (cfg.electronQuantumPotential.conductionBandNarrowingFraction < 0.0 ||
        cfg.electronQuantumPotential.conductionBandNarrowingFraction > 1.0 ||
        !std::isfinite(
            cfg.electronQuantumPotential.conductionBandNarrowingFraction)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "conduction_band_narrowing_fraction must be finite in [0,1].");
    }
    if (cfg.electronQuantumPotential.outerAbsoluteTolerance_V.has_value() &&
        (!(cfg.electronQuantumPotential.outerAbsoluteTolerance_V.value() >= 0.0) ||
         !std::isfinite(
             cfg.electronQuantumPotential.outerAbsoluteTolerance_V.value()))) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "outer_absolute_tolerance_V must be finite and non-negative.");
    }
    if (cfg.electronQuantumPotential.outerAcceleration != "none" &&
        cfg.electronQuantumPotential.outerAcceleration != "aitken") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.outer_acceleration "
            "must be 'none' or 'aitken'.");
    }
    if (!(cfg.electronQuantumPotential.outerRelaxation > 0.0) ||
        !std::isfinite(cfg.electronQuantumPotential.outerRelaxation) ||
        !(cfg.electronQuantumPotential.outerRelaxationMin > 0.0) ||
        !std::isfinite(cfg.electronQuantumPotential.outerRelaxationMin) ||
        !(cfg.electronQuantumPotential.outerRelaxationMax >=
          cfg.electronQuantumPotential.outerRelaxationMin) ||
        !std::isfinite(cfg.electronQuantumPotential.outerRelaxationMax) ||
        cfg.electronQuantumPotential.outerRelaxation <
            cfg.electronQuantumPotential.outerRelaxationMin ||
        cfg.electronQuantumPotential.outerRelaxation >
            cfg.electronQuantumPotential.outerRelaxationMax) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential outer relaxation "
            "must be finite, positive, and within [outer_relaxation_min, "
            "outer_relaxation_max].");
    }
    (void)usesFermiDirac(cfg.carrierStatistics);

    return cfg;
}

NewtonSolver::NewtonSolver(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    NewtonConfig cfg,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    ContactSpecsMap contactSpecs)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , contactBiases_(contactBiases)
    , contactSpecs_(std::move(contactSpecs))
    , cfg_(cfg)
    , fixedCharges_(std::move(fixedCharges))
    , sheetCharges_(std::move(sheetCharges))
{
    if (cfg_.jacobian != "analytic" && cfg_.jacobian != "finite_difference")
        throw std::invalid_argument(
            "NewtonSolver: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg_.quasiFermiReference != "none" &&
        cfg_.quasiFermiReference != "contact_majority") {
        throw std::invalid_argument(
            "NewtonSolver: quasi_fermi_reference must be 'none' or "
            "'contact_majority'.");
    }
    if (cfg_.residualNorm != "block" && cfg_.residualNorm != "l2")
        throw std::invalid_argument(
            "NewtonSolver: residual_norm must be 'block' or 'l2'.");
    if (cfg_.stallResidualFloor < 0.0 || !std::isfinite(cfg_.stallResidualFloor)) {
        throw std::invalid_argument(
            "NewtonSolver: stall_residual_floor must be non-negative and finite.");
    }
    if (cfg_.poissonLineSearchStallResidualFloor < 0.0 ||
        !std::isfinite(cfg_.poissonLineSearchStallResidualFloor)) {
        throw std::invalid_argument(
            "NewtonSolver: poisson_line_search_stall_residual_floor must be non-negative and finite.");
    }
    if (cfg_.poissonLineSearchStallRelativeIncrease < 0.0 ||
        !std::isfinite(cfg_.poissonLineSearchStallRelativeIncrease)) {
        throw std::invalid_argument(
            "NewtonSolver: poisson_line_search_stall_relative_increase must be non-negative and finite.");
    }
    if (cfg_.poissonLineSearchStallCarrierResidualFloor < 0.0 ||
        !std::isfinite(cfg_.poissonLineSearchStallCarrierResidualFloor)) {
        throw std::invalid_argument(
            "NewtonSolver: poisson_line_search_stall_carrier_residual_floor must be non-negative and finite.");
    }
    if (cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V < 0.0 ||
        !std::isfinite(cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V)) {
        throw std::invalid_argument(
            "NewtonSolver: poisson_line_search_stall_contact_majority_qf_drop_limit_V must be non-negative and finite.");
    }
    if (cfg_.carrierRegularizationScale < 0.0 ||
        !std::isfinite(cfg_.carrierRegularizationScale)) {
        throw std::invalid_argument(
            "NewtonSolver: carrier_regularization_scale must be non-negative and finite.");
    }
    if (cfg_.carrierDiagonalFloor.scale < 0.0 ||
        !std::isfinite(cfg_.carrierDiagonalFloor.scale)) {
        throw std::invalid_argument(
            "NewtonSolver: carrier_diagonal_floor scale must be non-negative and finite.");
    }
    if (cfg_.carrierDiagonalFloor.minorityDensityRatio < 0.0 ||
        !std::isfinite(cfg_.carrierDiagonalFloor.minorityDensityRatio)) {
        throw std::invalid_argument(
            "NewtonSolver: carrier_diagonal_floor minority_density_ratio must be non-negative and finite.");
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
        cfg_.temperature_K, epsRef, autoInputs, cfg_.unitScalingRefs, cfg_.inputScaling.unitSystem());

    scaling.enabled = true;
    scaling.V0 = sc.V0();
    scaling.C0 = sc.C0();
    scaling.mu0 = sc.mu0();
    scaling.D0 = sc.D0();
    scaling.L0 = sc.L0();
    scaling.permittivityReference_F_per_m = epsRef;
    scaling.unitSystem = cfg_.inputScaling.unitSystem();
    scaling.chargeAreaFactor = cfg_.inputScaling.unitSystem().chargeAreaFactor();
    scaling.chargeLineFactor = cfg_.inputScaling.unitSystem().chargeLineFactor();
    scaling.fieldFromCoordinateDeltaFactor = cfg_.inputScaling.unitSystem().fieldFromCoordinateDeltaFactor();
    scaling.currentDensityLineIntegralFactor =
        cfg_.inputScaling.unitSystem().currentDensityAM2PerInternal() *
        cfg_.inputScaling.unitSystem().lengthMPerInternal();
    return scaling;
}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler) const
{
    return buildBoundaryConditions(assembler, contactBiases_);
}

void NewtonSolver::configureQuasiFermiReferences(
    CoupledDDAssembler& assembler) const
{
    if (cfg_.quasiFermiReference == "none")
        return;

    Real electronReference = 0.0;
    Real holeReference = 0.0;
    Real strongestElectronDoping = 0.0;
    Real strongestHoleDoping = 0.0;
    for (const Contact& contact : mesh_.contacts()) {
        const auto biasIt = contactBiases_.find(contact.name);
        if (biasIt == contactBiases_.end() || contact.node_ids.empty())
            continue;
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }

        Real meanNetDoping = 0.0;
        for (const Index node : contact.node_ids)
            meanNetDoping += doping_.netDoping(node);
        meanNetDoping /= static_cast<Real>(contact.node_ids.size());
        if (meanNetDoping > strongestElectronDoping) {
            strongestElectronDoping = meanNetDoping;
            electronReference = biasIt->second;
        }
        if (meanNetDoping < strongestHoleDoping) {
            strongestHoleDoping = meanNetDoping;
            holeReference = biasIt->second;
        }
    }
    assembler.setQuasiFermiReferences(electronReference, holeReference);
}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler,
    const std::unordered_map<std::string, Real>& contactBiases) const
{
    CoupledDDBoundaryConditions bcs;
    const auto& ni = assembler.intrinsicDensity();
    const auto& Nc = assembler.electronDensityOfStates();
    const auto& Nv = assembler.holeDensityOfStates();
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
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            for (Index nid : contact.node_ids)
                bcs.psi[nid] = Vbias / potentialScale;
            continue;
        }
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
            const EquilibriumCarrierState equilibrium = equilibriumCarrierState(
                Ndop, niNode, Nc[nid], Nv[nid], Vt, cfg_.carrierStatistics);
            const double psiBuiltIn = equilibrium.potential;

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
    const auto& Nc = assembler.electronDensityOfStates();
    const auto& Nv = assembler.holeDensityOfStates();

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
        const EquilibriumCarrierState equilibrium = equilibriumCarrierState(
            doping_.netDoping(nid), niNode, Nc[nid], Nv[nid], Vt,
            cfg_.carrierStatistics);
        sol.psi(i) = equilibrium.potential;
    }

    for (const auto& [nid, value] : bcs.psi)
        sol.psi(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phin)
        sol.phin(static_cast<int>(nid)) = value * potentialScale;
    for (const auto& [nid, value] : bcs.phip)
        sol.phip(static_cast<int>(nid)) = value * potentialScale;

    for (int i = 0; i < N; ++i) {
        const Index nid = static_cast<Index>(i);
        sol.n(i) = electronDensity(
            ni[nid], Nc[nid], sol.psi(i), sol.phin(i), Vt,
            cfg_.carrierStatistics);
        sol.p(i) = holeDensity(
            ni[nid], Nv[nid], sol.psi(i), sol.phip(i), Vt,
            cfg_.carrierStatistics);
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
    sol.electronQuantumPotential = assembler.electronQuantumPotential();
    const int N = static_cast<int>(assembler.numNodes());
    sol.phinIncrement = x.segment(N, N) * potentialScale;
    sol.phipIncrement = x.segment(2 * N, N) * potentialScale;
    sol.electronQfReference_V = assembler.electronQuasiFermiReference();
    sol.holeQfReference_V = assembler.holeQuasiFermiReference();
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    const DDScalingSpec scaling = buildScalingSpec();
    auto assembler = std::make_shared<CoupledDDAssembler>(
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(*assembler);
    return assembler;
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

Real NewtonSolver::maxContactMajorityQuasiFermiDrop(const DDSolution& state) const
{
    if (state.phin.size() < static_cast<int>(mesh_.numNodes()) ||
        state.phip.size() < static_cast<int>(mesh_.numNodes()) ||
        state.n.size() < static_cast<int>(mesh_.numNodes()) ||
        state.p.size() < static_cast<int>(mesh_.numNodes())) {
        return std::numeric_limits<Real>::infinity();
    }

    Real maxDrop = 0.0;
    for (const Contact& contact : mesh_.contacts()) {
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        std::vector<bool> isContactNode(mesh_.numNodes(), false);
        Real netDopingSum = 0.0;
        int netDopingCount = 0;
        for (Index node : contact.node_ids) {
            if (node >= mesh_.numNodes())
                continue;
            isContactNode[node] = true;
            netDopingSum += doping_.netDoping(node);
            ++netDopingCount;
        }
        if (netDopingCount == 0)
            continue;

        const Real meanNetDoping = netDopingSum / static_cast<Real>(netDopingCount);
        const bool electronMajority = meanNetDoping > 0.0;
        const bool holeMajority = meanNetDoping < 0.0;
        for (Index e = 0; e < mesh_.numEdges(); ++e) {
            const Edge& edge = mesh_.getEdge(e);
            const bool n0Contact = isContactNode[edge.n0];
            const bool n1Contact = isContactNode[edge.n1];
            if (n0Contact == n1Contact)
                continue;
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            const bool electronTransportEdge =
                state.n(i) > 0.0 && state.n(j) > 0.0 &&
                std::isfinite(state.n(i)) && std::isfinite(state.n(j));
            const bool holeTransportEdge =
                state.p(i) > 0.0 && state.p(j) > 0.0 &&
                std::isfinite(state.p(i)) && std::isfinite(state.p(j));
            if ((electronMajority || !holeMajority) && electronTransportEdge)
                maxDrop = std::max(maxDrop, std::abs(state.phin(i) - state.phin(j)));
            if ((holeMajority || !electronMajority) && holeTransportEdge)
                maxDrop = std::max(maxDrop, std::abs(state.phip(i) - state.phip(j)));
        }
    }
    return maxDrop;
}
NewtonResidualEvaluation NewtonSolver::evaluateResidual(const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    if (state.electronQuantumPotential.size() ==
        static_cast<int>(mesh_.numNodes())) {
        assembler.setElectronQuantumPotential(state.electronQuantumPotential);
    }
    configureQuasiFermiReferences(assembler);
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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

std::vector<NewtonFeedbackSubstitutionEvaluation>
NewtonSolver::evaluateFeedbackSubstitutions(
    const DDSolution& state,
    const DDSolution& replacementState) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    const auto requireSize = [N](const VectorXd& values, const char* name) {
        if (values.size() != N) {
            throw std::invalid_argument(
                std::string("NewtonSolver::evaluateFeedbackSubstitutions: ") +
                name + " size mismatch.");
        }
    };
    requireSize(state.psi, "state psi");
    requireSize(state.phin, "state phin");
    requireSize(state.phip, "state phip");
    requireSize(replacementState.phin, "replacement phin");
    requireSize(replacementState.phip, "replacement phip");
    requireSize(replacementState.n, "replacement electron density");
    requireSize(replacementState.p, "replacement hole density");

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});

    const VectorXd baselineElectronDensity = assembler.electronDensity(x);
    const VectorXd baselineHoleDensity = assembler.holeDensity(x);
    CoupledDDFeedbackStateSubstitution replacement;
    replacement.electronDensity = replacementState.n;
    replacement.holeDensity = replacementState.p;
    replacement.electronQuasiFermi_V = replacementState.phin;
    replacement.holeQuasiFermi_V = replacementState.phip;

    // Keep every contact value on the baseline state.  The intervention is
    // confined to interior operator inputs, so paired boundary rows are
    // exactly unchanged.
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        const int i = static_cast<int>(node);
        replacement.electronDensity(i) = baselineElectronDensity(i);
        replacement.electronQuasiFermi_V(i) = state.phin(i);
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        const int i = static_cast<int>(node);
        replacement.holeDensity(i) = baselineHoleDensity(i);
        replacement.holeQuasiFermi_V(i) = state.phip(i);
    }

    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);
    const SparseMatrixd carrierBlock =
        sparseBlock(J, N, N, 2 * N, 2 * N);
    VectorXd targetStep = VectorXd::Zero(3 * N);
    targetStep.segment(N, N) =
        (replacement.electronQuasiFermi_V - state.phin) / potentialScale;
    targetStep.segment(2 * N, N) =
        (replacement.holeQuasiFermi_V - state.phip) / potentialScale;
    const VectorXd desiredResidual = -(J * targetStep);

    struct Variant {
        const char* name;
        bool electronDensity;
        bool holeDensity;
        bool electronQuasiFermi;
        bool holeQuasiFermi;
    };
    const std::vector<Variant> variants = {
        {"baseline", false, false, false, false},
        {"electron_density_only", true, false, false, false},
        {"hole_density_only", false, true, false, false},
        {"density_only", true, true, false, false},
        {"electron_qfp_only", false, false, true, false},
        {"hole_qfp_only", false, false, false, true},
        {"qfp_only", false, false, true, true},
        {"density_qfp", true, true, true, true},
    };

    std::vector<NewtonFeedbackSubstitutionEvaluation> evaluations;
    evaluations.reserve(variants.size());
    LinearSolver linearSolver;
    for (const Variant& variant : variants) {
        CoupledDDFeedbackStateSubstitution active = replacement;
        active.replaceElectronDensity = variant.electronDensity;
        active.replaceHoleDensity = variant.holeDensity;
        active.replaceElectronQuasiFermi = variant.electronQuasiFermi;
        active.replaceHoleQuasiFermi = variant.holeQuasiFermi;
        const bool replacesDensity =
            variant.electronDensity || variant.holeDensity;
        const bool replacesQuasiFermi =
            variant.electronQuasiFermi || variant.holeQuasiFermi;

        const VectorXd raw = (replacesDensity || replacesQuasiFermi)
            ? assembler.feedbackSubstitutionResidual(x, bcs, active)
            : assembler.residual(x, bcs);
        const auto terms = (replacesDensity || replacesQuasiFermi)
            ? assembler.feedbackSubstitutionCarrierContinuityTermDiagnostics(
                x, bcs, active)
            : assembler.carrierContinuityTermDiagnostics(x, bcs);
        VectorXd step = linearSolver.solve(J, -raw);
        const Real rawStepNorm = step.norm();
        applyConfiguredStepCapsAndPoissonRecorrection(
            step, J, raw, cfg_, N, potentialScale, doping_);
        VectorXd carrierOnlyStep = VectorXd::Zero(3 * N);
        carrierOnlyStep.segment(N, 2 * N) = linearSolver.solve(
            carrierBlock, -raw.segment(N, 2 * N));
        const Real carrierOnlyRawStepNorm = carrierOnlyStep.norm();
        applyConfiguredStepCaps(
            carrierOnlyStep, cfg_, N, potentialScale, doping_);
        const VectorXd productionTrialRaw = assembler.residual(x + step, bcs);

        NewtonFeedbackSubstitutionEvaluation evaluation;
        evaluation.variant = variant.name;
        evaluation.replacesDensity = replacesDensity;
        evaluation.replacesQuasiFermi = replacesQuasiFermi;
        evaluation.residual.raw = raw;
        evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
        evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
        evaluation.residual.scaledState = assembler.usesScaledState();
        evaluation.residual.potentialScale = potentialScale;
        evaluation.productionTrialResidual.raw = productionTrialRaw;
        evaluation.productionTrialResidual.blockNorms =
            blockResidualInfo(productionTrialRaw, mesh_.numNodes());
        evaluation.productionTrialResidual.intrinsicDensity =
            assembler.intrinsicDensity();
        evaluation.productionTrialResidual.scaledState =
            assembler.usesScaledState();
        evaluation.productionTrialResidual.potentialScale = potentialScale;
        evaluation.carrierTerms = terms;
        evaluation.desiredResidual = desiredResidual;
        evaluation.deltaPsi = step.segment(0, N) * potentialScale;
        evaluation.deltaPhin = step.segment(N, N) * potentialScale;
        evaluation.deltaPhip = step.segment(2 * N, N) * potentialScale;
        evaluation.carrierOnlyDeltaPhin =
            carrierOnlyStep.segment(N, N) * potentialScale;
        evaluation.carrierOnlyDeltaPhip =
            carrierOnlyStep.segment(2 * N, N) * potentialScale;
        evaluation.rawStepNorm = rawStepNorm;
        evaluation.stepNorm = step.norm();
        evaluation.carrierOnlyRawStepNorm = carrierOnlyRawStepNorm;
        evaluation.carrierOnlyStepNorm = carrierOnlyStep.norm();
        evaluations.push_back(std::move(evaluation));
    }
    return evaluations;
}

NewtonPoissonQfpCrossBlockEvaluation
NewtonSolver::evaluatePoissonQfpCrossBlockDecomposition(
    const DDSolution& state,
    const DDSolution& replacementState) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    const auto requireSize = [N](const VectorXd& values, const char* name) {
        if (values.size() != N) {
            throw std::invalid_argument(
                std::string(
                    "NewtonSolver::evaluatePoissonQfpCrossBlockDecomposition: ") +
                name + " size mismatch.");
        }
    };
    requireSize(state.psi, "state psi");
    requireSize(state.phin, "state phin");
    requireSize(state.phip, "state phip");
    requireSize(replacementState.phin, "replacement phin");
    requireSize(replacementState.phip, "replacement phip");

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});

    CoupledDDFeedbackStateSubstitution replacement;
    replacement.electronQuasiFermi_V = replacementState.phin;
    replacement.holeQuasiFermi_V = replacementState.phip;
    replacement.replaceElectronQuasiFermi = true;
    replacement.replaceHoleQuasiFermi = true;
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        replacement.electronQuasiFermi_V(static_cast<int>(node)) =
            state.phin(static_cast<int>(node));
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        replacement.holeQuasiFermi_V(static_cast<int>(node)) =
            state.phip(static_cast<int>(node));
    }

    const VectorXd raw =
        assembler.feedbackSubstitutionResidual(x, bcs, replacement);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);
    const SparseMatrixd A = sparseBlock(J, 0, 0, N, N);
    const SparseMatrixd B = sparseBlock(J, 0, N, N, 2 * N);
    const SparseMatrixd C = sparseBlock(J, N, 0, 2 * N, N);
    const SparseMatrixd D = sparseBlock(J, N, N, 2 * N, 2 * N);
    const VectorXd rPsi = raw.segment(0, N);
    const VectorXd rQfp = raw.segment(N, 2 * N);

    LinearSolver linearSolver;
    const VectorXd independentPsi = linearSolver.solve(A, -rPsi);
    const VectorXd independentQfp = linearSolver.solve(D, -rQfp);
    const VectorXd noPsiQfpPsi = independentPsi;
    const VectorXd noPsiQfpQfp =
        linearSolver.solve(D, -rQfp - C * noPsiQfpPsi);
    const VectorXd noQfpPsiQfp = independentQfp;
    const VectorXd noQfpPsiPsi =
        linearSolver.solve(A, -rPsi - B * noQfpPsiQfp);

    const Eigen::MatrixXd denseB = Eigen::MatrixXd(B);
    Eigen::MatrixXd aInverseB(N, 2 * N);
    for (int column = 0; column < 2 * N; ++column)
        aInverseB.col(column) = linearSolver.solve(A, denseB.col(column));
    const Eigen::MatrixXd schur =
        Eigen::MatrixXd(D) - Eigen::MatrixXd(C) * aInverseB;
    const Eigen::MatrixXd effectiveLoop =
        Eigen::MatrixXd(C) * aInverseB;
    SparseMatrixd sparseSchur = schur.sparseView();
    sparseSchur.makeCompressed();
    const VectorXd schurQfp = linearSolver.solve(
        sparseSchur, -rQfp - C * independentPsi);
    const VectorXd schurPsi =
        linearSolver.solve(A, -rPsi - B * schurQfp);

    const VectorXd fullRaw = linearSolver.solve(J, -raw);
    VectorXd fullCapped = fullRaw;
    applyConfiguredStepCapsAndPoissonRecorrection(
        fullCapped, J, raw, cfg_, N, potentialScale, doping_);
    VectorXd schurStep(3 * N);
    schurStep.segment(0, N) = schurPsi;
    schurStep.segment(N, 2 * N) = schurQfp;

    const auto makeDiagnosticAssembler =
        [&](const RecombinationModelConfig& recombination,
            const ImpactIonizationModelConfig& impact) {
            CoupledDDAssembler diagnostic(
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
                scaling,
                cfg_.carrierDiagonalFloor,
                cfg_.carrierStatistics);
            configureQuasiFermiReferences(diagnostic);
            return diagnostic;
        };
    RecombinationModelConfig noRecombination =
        recombinationModelConfig({"none"}, cfg_.taun, cfg_.taup);
    const ImpactIonizationModelConfig noImpact{};
    CoupledDDAssembler transportAssembler =
        makeDiagnosticAssembler(noRecombination, noImpact);
    CoupledDDAssembler recombinationAssembler =
        makeDiagnosticAssembler(recombinationConfig, noImpact);
    CoupledDDAssembler impactAssembler =
        makeDiagnosticAssembler(noRecombination, cfg_.impactIonization);
    const SparseMatrixd transportJ =
        transportAssembler.assembleJacobian(x, bcs);
    const SparseMatrixd recombinationJ =
        recombinationAssembler.assembleJacobian(x, bcs);
    const SparseMatrixd impactJ =
        impactAssembler.assembleJacobian(x, bcs);
    const SparseMatrixd transportC =
        sparseBlock(transportJ, N, 0, 2 * N, N);
    SparseMatrixd recombinationC = sparseBlock(
        recombinationJ - transportJ, N, 0, 2 * N, N);
    SparseMatrixd impactC = sparseBlock(
        impactJ - transportJ, N, 0, 2 * N, N);
    recombinationC.makeCompressed();
    impactC.makeCompressed();

    std::vector<NewtonSchurLoopComponentEvaluation> loopComponents;
    loopComponents.reserve(3);
    SparseMatrixd componentSum(2 * N, N);
    componentSum.setZero();
    const auto addLoopComponent =
        [&](const std::string& name, const SparseMatrixd& componentC) {
            const Eigen::MatrixXd componentLoop =
                Eigen::MatrixXd(componentC) * aInverseB;
            const Eigen::MatrixXd leaveOutC =
                Eigen::MatrixXd(C) - Eigen::MatrixXd(componentC);
            SparseMatrixd leaveOutSchur =
                (Eigen::MatrixXd(D) - leaveOutC * aInverseB).sparseView();
            leaveOutSchur.makeCompressed();
            const VectorXd leaveOutQfp = linearSolver.solve(
                leaveOutSchur,
                -rQfp - leaveOutC * independentPsi);
            SparseMatrixd onlySchur =
                (Eigen::MatrixXd(D)
                 - Eigen::MatrixXd(componentC) * aInverseB).sparseView();
            onlySchur.makeCompressed();
            const VectorXd onlyQfp = linearSolver.solve(
                onlySchur,
                -rQfp - Eigen::MatrixXd(componentC) * independentPsi);

            NewtonSchurLoopComponentEvaluation component;
            component.name = name;
            component.jacobianQfpPsi = componentC;
            component.effectiveLoop = componentLoop.sparseView();
            component.effectiveLoop.makeCompressed();
            component.leaveOutDeltaPhin =
                leaveOutQfp.segment(0, N) * potentialScale;
            component.leaveOutDeltaPhip =
                leaveOutQfp.segment(N, N) * potentialScale;
            component.onlyDeltaPhin =
                onlyQfp.segment(0, N) * potentialScale;
            component.onlyDeltaPhip =
                onlyQfp.segment(N, N) * potentialScale;
            loopComponents.push_back(std::move(component));
            componentSum += componentC;
        };
    addLoopComponent("transport_boundary", transportC);
    addLoopComponent("srh_auger", recombinationC);
    addLoopComponent("sg_avalanche", impactC);
    const Eigen::MatrixXd componentLoopClosure =
        effectiveLoop - Eigen::MatrixXd(componentSum) * aInverseB;

    VectorXd qfpDirection = independentQfp;
    if (qfpDirection.norm() == 0.0) {
        qfpDirection.segment(0, N) =
            (replacement.electronQuasiFermi_V - state.phin) / potentialScale;
        qfpDirection.segment(N, N) =
            (replacement.holeQuasiFermi_V - state.phip) / potentialScale;
    }
    if (qfpDirection.norm() == 0.0) {
        throw std::runtime_error(
            "NewtonSolver::evaluatePoissonQfpCrossBlockDecomposition: "
            "cannot form a non-zero QFP finite-difference direction.");
    }
    qfpDirection.normalize();
    VectorXd psiDirection =
        -linearSolver.solve(A, B * qfpDirection);
    if (psiDirection.norm() == 0.0) {
        throw std::runtime_error(
            "NewtonSolver::evaluatePoissonQfpCrossBlockDecomposition: "
            "cannot form a non-zero Poisson-response finite-difference direction.");
    }
    psiDirection.normalize();
    constexpr Real directionalStep = 1.0e-7;
    VectorXd qfpPerturbation = VectorXd::Zero(3 * N);
    qfpPerturbation.segment(N, 2 * N) = qfpDirection;
    const VectorXd psiQfpFiniteDifference =
        (assembler.residual(x + directionalStep * qfpPerturbation, bcs)
         - assembler.residual(x - directionalStep * qfpPerturbation, bcs))
        / (2.0 * directionalStep);
    VectorXd psiPerturbation = VectorXd::Zero(3 * N);
    psiPerturbation.segment(0, N) = psiDirection;
    const VectorXd qfpPsiFiniteDifference =
        (assembler.residual(x + directionalStep * psiPerturbation, bcs)
         - assembler.residual(x - directionalStep * psiPerturbation, bcs))
        / (2.0 * directionalStep);
    const VectorXd analyticPsiQfp = B * qfpDirection;
    const VectorXd analyticQfpPsi = C * psiDirection;
    const VectorXd finiteDifferencePsiQfp =
        psiQfpFiniteDifference.segment(0, N);
    const VectorXd finiteDifferenceQfpPsi =
        qfpPsiFiniteDifference.segment(N, 2 * N);

    NewtonPoissonQfpCrossBlockEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.jacobianPsiPsi = A;
    evaluation.jacobianPsiQfp = B;
    evaluation.jacobianQfpPsi = C;
    evaluation.jacobianQfpQfp = D;
    evaluation.effectiveSchurLoop = effectiveLoop.sparseView();
    evaluation.effectiveSchurLoop.makeCompressed();
    evaluation.loopComponents = std::move(loopComponents);
    evaluation.targetDeltaPhin =
        replacement.electronQuasiFermi_V - state.phin;
    evaluation.targetDeltaPhip =
        replacement.holeQuasiFermi_V - state.phip;
    evaluation.independentDeltaPsi = independentPsi * potentialScale;
    evaluation.independentDeltaPhin =
        independentQfp.segment(0, N) * potentialScale;
    evaluation.independentDeltaPhip =
        independentQfp.segment(N, N) * potentialScale;
    evaluation.noPsiQfpDeltaPsi = noPsiQfpPsi * potentialScale;
    evaluation.noPsiQfpDeltaPhin =
        noPsiQfpQfp.segment(0, N) * potentialScale;
    evaluation.noPsiQfpDeltaPhip =
        noPsiQfpQfp.segment(N, N) * potentialScale;
    evaluation.noQfpPsiDeltaPsi = noQfpPsiPsi * potentialScale;
    evaluation.noQfpPsiDeltaPhin =
        noQfpPsiQfp.segment(0, N) * potentialScale;
    evaluation.noQfpPsiDeltaPhip =
        noQfpPsiQfp.segment(N, N) * potentialScale;
    evaluation.schurDeltaPsi = schurPsi * potentialScale;
    evaluation.schurDeltaPhin =
        schurQfp.segment(0, N) * potentialScale;
    evaluation.schurDeltaPhip =
        schurQfp.segment(N, N) * potentialScale;
    evaluation.fullRawDeltaPsi =
        fullRaw.segment(0, N) * potentialScale;
    evaluation.fullRawDeltaPhin =
        fullRaw.segment(N, N) * potentialScale;
    evaluation.fullRawDeltaPhip =
        fullRaw.segment(2 * N, N) * potentialScale;
    evaluation.fullCappedDeltaPsi =
        fullCapped.segment(0, N) * potentialScale;
    evaluation.fullCappedDeltaPhin =
        fullCapped.segment(N, N) * potentialScale;
    evaluation.fullCappedDeltaPhip =
        fullCapped.segment(2 * N, N) * potentialScale;
    evaluation.psiQfpProduct = B * schurQfp;
    evaluation.qfpPsiProduct = C * schurPsi;
    evaluation.qfpFiniteDifferenceDirectionPhin =
        qfpDirection.segment(0, N) * potentialScale;
    evaluation.qfpFiniteDifferenceDirectionPhip =
        qfpDirection.segment(N, N) * potentialScale;
    evaluation.psiFiniteDifferenceDirection =
        psiDirection * potentialScale;
    evaluation.analyticPsiQfpDirectionalDerivative = analyticPsiQfp;
    evaluation.finiteDifferencePsiQfpDirectionalDerivative =
        finiteDifferencePsiQfp;
    evaluation.analyticQfpPsiDirectionalDerivative = analyticQfpPsi;
    evaluation.finiteDifferenceQfpPsiDirectionalDerivative =
        finiteDifferenceQfpPsi;
    evaluation.jacobianPsiPsiCondition =
        matrixConditionEstimate(Eigen::MatrixXd(A));
    evaluation.jacobianPsiPsiEquilibratedCondition =
        matrixConditionEstimate(l2Equilibrated(Eigen::MatrixXd(A)));
    evaluation.jacobianQfpQfpCondition =
        matrixConditionEstimate(Eigen::MatrixXd(D));
    evaluation.jacobianQfpQfpEquilibratedCondition =
        matrixConditionEstimate(l2Equilibrated(Eigen::MatrixXd(D)));
    evaluation.schurCondition = matrixConditionEstimate(schur);
    evaluation.schurEquilibratedCondition =
        matrixConditionEstimate(l2Equilibrated(schur));
    evaluation.effectiveSchurLoopCondition =
        matrixConditionEstimate(effectiveLoop);
    evaluation.jacobianPsiPsiNorm = A.norm();
    evaluation.jacobianPsiQfpNorm = B.norm();
    evaluation.jacobianQfpPsiNorm = C.norm();
    evaluation.jacobianQfpQfpNorm = D.norm();
    evaluation.fullLinearClosureNorm = (J * fullRaw + raw).norm();
    evaluation.schurClosureNorm = (J * schurStep + raw).norm();
    evaluation.schurRelativeClosure =
        evaluation.schurClosureNorm / std::max<Real>(raw.norm(), 1.0);
    evaluation.loopComponentClosureNorm = componentLoopClosure.norm();
    evaluation.finiteDifferenceRelativeStep = directionalStep;
    evaluation.psiQfpDirectionalDerivativeRelativeError =
        relativeVectorDifference(analyticPsiQfp, finiteDifferencePsiQfp);
    evaluation.qfpPsiDirectionalDerivativeRelativeError =
        relativeVectorDifference(analyticQfpPsi, finiteDifferenceQfpPsi);
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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

NewtonCarrierBlockDecompositionEvaluation
NewtonSolver::evaluateCarrierBlockDecomposition(const DDSolution& state) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateCarrierBlockDecomposition: state size must "
            "match the mesh node count.");
    }

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    const DDScalingSpec scaling = buildScalingSpec();
    const auto makeRecombinationConfig =
        [&](const std::vector<std::string>& models) {
            RecombinationModelConfig config =
                recombinationModelConfig(
                    models, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
            config.augerCn = cfg_.augerCn;
            config.augerCp = cfg_.augerCp;
            if (models.size() != 1 || models.front() != "none")
                config.bandToBand = cfg_.bandToBand;
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
            CoupledDDAssembler assembler(
                mesh_, matdb_, doping_, Vt, mobilityConfig, recombination,
                cfg_.bandgapNarrowing, impact, fixedCharges_, sheetCharges_,
                scaling, cfg_.carrierDiagonalFloor, cfg_.carrierStatistics);
            configureQuasiFermiReferences(assembler);
            return assembler;
        };

    CoupledDDAssembler fullAssembler =
        makeAssembler(recombinationConfig, cfg_.impactIonization);
    CoupledDDAssembler noRecombinationAssembler =
        makeAssembler(noRecombinationConfig, cfg_.impactIonization);
    CoupledDDAssembler noImpactAssembler =
        makeAssembler(recombinationConfig, noImpactConfig);
    CoupledDDAssembler transportAssembler =
        makeAssembler(noRecombinationConfig, noImpactConfig);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(fullAssembler);
    const Real potentialScale =
        fullAssembler.usesScaledState() ? fullAssembler.potentialScale() : 1.0;
    const VectorXd x = fullAssembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const VectorXd raw = fullAssembler.residual(x, bcs);
    const VectorXd noRecombinationRaw =
        noRecombinationAssembler.residual(x, bcs);
    const VectorXd noImpactRaw = noImpactAssembler.residual(x, bcs);
    const VectorXd transportRaw = transportAssembler.residual(x, bcs);

    const auto carrierBlock = [&](CoupledDDAssembler& assembler) {
        return sparseBlock(
            assembler.assembleJacobian(x, bcs), N, N, 2 * N, 2 * N);
    };
    const SparseMatrixd fullBlock = carrierBlock(fullAssembler);
    const SparseMatrixd noRecombinationBlock =
        carrierBlock(noRecombinationAssembler);
    const SparseMatrixd noImpactBlock = carrierBlock(noImpactAssembler);
    const SparseMatrixd transportBlock = carrierBlock(transportAssembler);
    const VectorXd carrierResidual = raw.segment(N, 2 * N);
    const VectorXd allRowWeights = continuityRowWeights(
        fullAssembler, x, bcs, cfg_.continuityRowScaling);
    const VectorXd carrierRowWeights = allRowWeights.segment(N, 2 * N);

    std::vector<int> freeUnknowns;
    freeUnknowns.reserve(static_cast<std::size_t>(2 * N));
    for (int node = 0; node < N; ++node) {
        if (bcs.phin.find(static_cast<Index>(node)) == bcs.phin.end())
            freeUnknowns.push_back(node);
    }
    const int freeElectronUnknowns = static_cast<int>(freeUnknowns.size());
    for (int node = 0; node < N; ++node) {
        if (bcs.phip.find(static_cast<Index>(node)) == bcs.phip.end())
            freeUnknowns.push_back(N + node);
    }
    const int freeHoleUnknowns =
        static_cast<int>(freeUnknowns.size()) - freeElectronUnknowns;

    const auto reduceMatrix = [&](const SparseMatrixd& matrix) {
        Eigen::MatrixXd reduced = Eigen::MatrixXd::Zero(
            static_cast<int>(freeUnknowns.size()),
            static_cast<int>(freeUnknowns.size()));
        for (int row = 0; row < reduced.rows(); ++row) {
            for (int col = 0; col < reduced.cols(); ++col) {
                reduced(row, col) = matrix.coeff(
                    freeUnknowns[static_cast<std::size_t>(row)],
                    freeUnknowns[static_cast<std::size_t>(col)]);
            }
        }
        return reduced;
    };
    const auto reduceVector = [&](const VectorXd& vector) {
        VectorXd reduced(static_cast<int>(freeUnknowns.size()));
        for (int row = 0; row < reduced.size(); ++row)
            reduced(row) = vector(freeUnknowns[static_cast<std::size_t>(row)]);
        return reduced;
    };
    const Eigen::MatrixXd fullReduced = reduceMatrix(fullBlock);
    const Eigen::MatrixXd noRecombinationReduced =
        reduceMatrix(noRecombinationBlock);
    const Eigen::MatrixXd noImpactReduced = reduceMatrix(noImpactBlock);
    const Eigen::MatrixXd transportReduced = reduceMatrix(transportBlock);
    const Eigen::MatrixXd recombinationReduced =
        noImpactReduced - transportReduced;
    const Eigen::MatrixXd avalancheReduced = fullReduced - noImpactReduced;
    Eigen::MatrixXd avalancheDiagonalReduced = avalancheReduced;
    avalancheDiagonalReduced.block(
        0, freeElectronUnknowns,
        freeElectronUnknowns, freeHoleUnknowns).setZero();
    avalancheDiagonalReduced.block(
        freeElectronUnknowns, 0,
        freeHoleUnknowns, freeElectronUnknowns).setZero();
    const Eigen::MatrixXd avalancheCrossReduced =
        avalancheReduced - avalancheDiagonalReduced;
    const VectorXd transportReducedRhs = -reduceVector(transportRaw.segment(N, 2 * N));
    const VectorXd recombinationReducedRhs = -reduceVector(
        noImpactRaw.segment(N, 2 * N) - transportRaw.segment(N, 2 * N));
    const VectorXd avalancheReducedRhs = -reduceVector(
        raw.segment(N, 2 * N) - noImpactRaw.segment(N, 2 * N));
    const VectorXd reducedRowWeights = reduceVector(carrierRowWeights);
    Eigen::MatrixXd rowScaledReduced = fullReduced;
    for (int row = 0; row < rowScaledReduced.rows(); ++row)
        rowScaledReduced.row(row) *= reducedRowWeights(row);

    NewtonCarrierBlockDecompositionEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = fullAssembler.intrinsicDensity();
    evaluation.residual.scaledState = fullAssembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.freeElectronUnknowns = static_cast<Index>(freeElectronUnknowns);
    evaluation.freeHoleUnknowns = static_cast<Index>(freeHoleUnknowns);
    evaluation.rawCondition = matrixConditionEstimate(fullReduced);
    evaluation.rowScaledCondition = matrixConditionEstimate(rowScaledReduced);
    evaluation.l2EquilibratedCondition =
        matrixConditionEstimate(l2Equilibrated(fullReduced));

    const auto subBlockNorm = [&](const Eigen::MatrixXd& matrix,
                                  int rowStart,
                                  int rowCount,
                                  int colStart,
                                  int colCount) {
        if (rowCount == 0 || colCount == 0)
            return Real{0.0};
        return matrix.block(rowStart, colStart, rowCount, colCount).norm();
    };
    evaluation.electronElectronNorm = subBlockNorm(
        fullReduced, 0, freeElectronUnknowns, 0, freeElectronUnknowns);
    evaluation.electronHoleNorm = subBlockNorm(
        fullReduced, 0, freeElectronUnknowns,
        freeElectronUnknowns, freeHoleUnknowns);
    evaluation.holeElectronNorm = subBlockNorm(
        fullReduced, freeElectronUnknowns, freeHoleUnknowns,
        0, freeElectronUnknowns);
    evaluation.holeHoleNorm = subBlockNorm(
        fullReduced, freeElectronUnknowns, freeHoleUnknowns,
        freeElectronUnknowns, freeHoleUnknowns);
    const auto crossNorm = [&](const Eigen::MatrixXd& matrix) {
        return std::hypot(
            subBlockNorm(matrix, 0, freeElectronUnknowns,
                         freeElectronUnknowns, freeHoleUnknowns),
            subBlockNorm(matrix, freeElectronUnknowns, freeHoleUnknowns,
                         0, freeElectronUnknowns));
    };
    evaluation.crossCarrierNormFraction = fullReduced.norm() > 0.0
        ? crossNorm(fullReduced) / fullReduced.norm()
        : 0.0;
    evaluation.recombinationCrossNorm =
        crossNorm(fullReduced - noRecombinationReduced);
    evaluation.avalancheCrossNorm = crossNorm(fullReduced - noImpactReduced);
    evaluation.transportCrossNorm = crossNorm(transportReduced);

    const auto positiveSpread = [](const VectorXd& values) {
        Real smallest = std::numeric_limits<Real>::infinity();
        Real largest = 0.0;
        for (int i = 0; i < values.size(); ++i) {
            if (values(i) > 0.0 && std::isfinite(values(i))) {
                smallest = std::min(smallest, values(i));
                largest = std::max(largest, values(i));
            }
        }
        return std::isfinite(smallest) && smallest > 0.0
            ? largest / smallest
            : Real{0.0};
    };
    VectorXd columnNorms(fullReduced.cols());
    VectorXd rowNorms(fullReduced.rows());
    for (int i = 0; i < fullReduced.cols(); ++i)
        columnNorms(i) = fullReduced.col(i).norm();
    for (int i = 0; i < fullReduced.rows(); ++i)
        rowNorms(i) = fullReduced.row(i).norm();
    evaluation.freeColumnNormSpread = positiveSpread(columnNorms);
    evaluation.freeRowNormSpread = positiveSpread(rowNorms);
    evaluation.rowWeightSpread = positiveSpread(reducedRowWeights);

    LinearSolver linearSolver;
    const auto solve = [&](const std::string& name,
                           const SparseMatrixd& matrix,
                           const VectorXd& equationResidual) {
        NewtonCarrierBlockSolveVariantEvaluation result;
        result.name = name;
        const VectorXd delta = linearSolver.solve(matrix, -equationResidual);
        result.deltaPhin = delta.head(N) * potentialScale;
        result.deltaPhip = delta.tail(N) * potentialScale;
        result.scaledStepNorm = delta.norm();
        result.physicalStepNorm_V = delta.norm() * potentialScale;
        result.relativeLinearClosure = equationResidual.norm() > 0.0
            ? (matrix * delta + equationResidual).norm() /
                equationResidual.norm()
            : (matrix * delta).norm();
        return std::pair<NewtonCarrierBlockSolveVariantEvaluation, VectorXd>{
            std::move(result), delta};
    };
    auto [fullSolve, fullDelta] = solve("full", fullBlock, carrierResidual);
    evaluation.solveVariants.push_back(std::move(fullSolve));
    auto appendVariant = [&](NewtonCarrierBlockSolveVariantEvaluation variant,
                             const VectorXd& delta) {
        const Real fullNorm = fullDelta.norm();
        variant.relativeDifferenceFromFull = fullNorm > 0.0
            ? (delta - fullDelta).norm() / fullNorm
            : delta.norm();
        const Real denominator = delta.norm() * fullNorm;
        variant.cosineWithFull = denominator > 0.0
            ? delta.dot(fullDelta) / denominator
            : 0.0;
        evaluation.solveVariants.push_back(std::move(variant));
    };
    evaluation.solveVariants.front().relativeDifferenceFromFull = 0.0;
    evaluation.solveVariants.front().cosineWithFull = 1.0;

    const SparseMatrixd rowScaledBlock =
        leftScaleRows(fullBlock, carrierRowWeights);
    const VectorXd rowScaledResidual =
        carrierResidual.cwiseProduct(carrierRowWeights);
    auto [rowScaledSolve, rowScaledDelta] =
        solve("row_scaled", rowScaledBlock, rowScaledResidual);
    appendVariant(std::move(rowScaledSolve), rowScaledDelta);

    Eigen::MatrixXd noCrossDense(fullBlock);
    noCrossDense.block(0, N, N, N).setZero();
    noCrossDense.block(N, 0, N, N).setZero();
    SparseMatrixd noCrossBlock = noCrossDense.sparseView();
    auto [noCrossSolve, noCrossDelta] =
        solve("no_cross_carrier", noCrossBlock, carrierResidual);
    appendVariant(std::move(noCrossSolve), noCrossDelta);
    auto [noRecombinationSolve, noRecombinationDelta] =
        solve("no_recombination", noRecombinationBlock, carrierResidual);
    appendVariant(std::move(noRecombinationSolve), noRecombinationDelta);
    auto [noAvalancheSolve, noAvalancheDelta] =
        solve("no_avalanche", noImpactBlock, carrierResidual);
    appendVariant(std::move(noAvalancheSolve), noAvalancheDelta);
    auto [transportSolve, transportDelta] =
        solve("transport_only", transportBlock, carrierResidual);
    appendVariant(std::move(transportSolve), transportDelta);

    evaluation.columns.reserve(freeUnknowns.size());
    for (int reducedColumn = 0;
         reducedColumn < fullReduced.cols(); ++reducedColumn) {
        const int carrierDof =
            freeUnknowns[static_cast<std::size_t>(reducedColumn)];
        const bool electron = carrierDof < N;
        const Real columnNorm = fullReduced.col(reducedColumn).norm();
        const Real electronRows = freeElectronUnknowns > 0
            ? fullReduced.col(reducedColumn).head(freeElectronUnknowns).norm()
            : 0.0;
        const Real holeRows = freeHoleUnknowns > 0
            ? fullReduced.col(reducedColumn).tail(freeHoleUnknowns).norm()
            : 0.0;
        NewtonCarrierBlockColumnDiagnostic column;
        column.carrier = electron ? "electron" : "hole";
        column.nodeId = static_cast<Index>(electron ? carrierDof : carrierDof - N);
        column.reducedColumn = static_cast<Index>(reducedColumn);
        column.diagonal = fullReduced(reducedColumn, reducedColumn);
        column.columnL2Norm = columnNorm;
        column.electronRowL2Norm = electronRows;
        column.holeRowL2Norm = holeRows;
        column.diagonalFraction = columnNorm > 0.0
            ? std::abs(column.diagonal) / columnNorm
            : 0.0;
        column.crossCarrierRowFraction = columnNorm > 0.0
            ? (electron ? holeRows : electronRows) / columnNorm
            : 0.0;
        column.continuityRowWeight = carrierRowWeights(carrierDof);
        column.residual = carrierResidual(carrierDof);
        column.fullDeltaQfp_V = fullDelta(carrierDof) * potentialScale;
        evaluation.columns.push_back(std::move(column));
    }

    if (!freeUnknowns.empty()) {
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            fullReduced, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const VectorXd singular = svd.singularValues();
        const VectorXd reducedRhs = -reduceVector(carrierResidual);
        const VectorXd reducedDelta = reduceVector(fullDelta);
        const VectorXd reducedNoCrossDelta = reduceVector(noCrossDelta);
        const VectorXd reducedNoRecombinationDelta =
            reduceVector(noRecombinationDelta);
        const VectorXd reducedNoAvalancheDelta = reduceVector(noAvalancheDelta);
        const VectorXd reducedTransportDelta = reduceVector(transportDelta);
        const Real rhsEnergy = reducedRhs.squaredNorm();
        const Real stepEnergy = reducedDelta.squaredNorm();
        const Real largest = singular.size() > 0 ? singular(0) : 0.0;
        evaluation.singularModes.reserve(
            static_cast<std::size_t>(singular.size()));
        for (int mode = 0; mode < singular.size(); ++mode) {
            const VectorXd right = svd.matrixV().col(mode);
            const VectorXd left = svd.matrixU().col(mode);
            const Real rhsProjection = left.dot(reducedRhs);
            // Project the step produced by the production sparse solver rather
            // than reconstructing it as U^T b / sigma.  The raw carrier block
            // spans many scale-separated directions below the dense-SVD rank
            // threshold; dropping those modes would hide where the actual
            // Newton update lives.
            const Real stepAmplitude = right.dot(reducedDelta);
            Eigen::Index topRight = 0;
            Eigen::Index topLeft = 0;
            right.cwiseAbs().maxCoeff(&topRight);
            left.cwiseAbs().maxCoeff(&topLeft);
            const int topRightDof =
                freeUnknowns[static_cast<std::size_t>(topRight)];
            const int topLeftDof =
                freeUnknowns[static_cast<std::size_t>(topLeft)];
            NewtonCarrierBlockSingularModeDiagnostic record;
            record.modeIndex = static_cast<Index>(mode);
            record.singularValue = singular(mode);
            record.relativeSingularValue = largest > 0.0
                ? singular(mode) / largest
                : 0.0;
            record.rhsProjection = rhsProjection;
            record.rhsEnergyFraction = rhsEnergy > 0.0
                ? rhsProjection * rhsProjection / rhsEnergy
                : 0.0;
            record.stepAmplitude = stepAmplitude;
            record.stepEnergyFraction = stepEnergy > 0.0
                ? stepAmplitude * stepAmplitude / stepEnergy
                : 0.0;
            record.transportJacobianProjection =
                left.dot(transportReduced * right);
            record.recombinationJacobianProjection =
                left.dot(recombinationReduced * right);
            record.avalancheDiagonalJacobianProjection =
                left.dot(avalancheDiagonalReduced * right);
            record.avalancheCrossJacobianProjection =
                left.dot(avalancheCrossReduced * right);
            record.jacobianProjectionClosure =
                record.transportJacobianProjection +
                record.recombinationJacobianProjection +
                record.avalancheDiagonalJacobianProjection +
                record.avalancheCrossJacobianProjection - singular(mode);
            record.transportRhsProjection = left.dot(transportReducedRhs);
            record.recombinationRhsProjection =
                left.dot(recombinationReducedRhs);
            record.avalancheRhsProjection = left.dot(avalancheReducedRhs);
            record.rhsProjectionClosure =
                record.transportRhsProjection +
                record.recombinationRhsProjection +
                record.avalancheRhsProjection - rhsProjection;
            record.noCrossCarrierStepAmplitude =
                right.dot(reducedNoCrossDelta);
            record.noRecombinationStepAmplitude =
                right.dot(reducedNoRecombinationDelta);
            record.noAvalancheStepAmplitude =
                right.dot(reducedNoAvalancheDelta);
            record.transportOnlyStepAmplitude =
                right.dot(reducedTransportDelta);
            record.rightElectronFraction = freeElectronUnknowns > 0
                ? right.head(freeElectronUnknowns).squaredNorm()
                : 0.0;
            record.leftElectronFraction = freeElectronUnknowns > 0
                ? left.head(freeElectronUnknowns).squaredNorm()
                : 0.0;
            record.topRightCarrier = topRightDof < N ? "electron" : "hole";
            record.topRightNode = static_cast<Index>(
                topRightDof < N ? topRightDof : topRightDof - N);
            record.topRightValue = right(topRight);
            record.topLeftCarrier = topLeftDof < N ? "electron" : "hole";
            record.topLeftNode = static_cast<Index>(
                topLeftDof < N ? topLeftDof : topLeftDof - N);
            record.topLeftValue = left(topLeft);
            evaluation.singularModes.push_back(std::move(record));
        }
    }
    return evaluation;
}

NewtonCarrierTermDiagnosticsEvaluation NewtonSolver::evaluateCarrierTermDiagnostics(
    const DDSolution& state) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
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
    std::vector<std::string> blocks,
    const std::string& finiteDifferenceMode) const
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
    if (finiteDifferenceMode != "double_symmetric" &&
        finiteDifferenceMode != "multiprecision_branch_resolved") {
        throw std::invalid_argument(
            "NewtonSolver::evaluateJacobianBlockAudit: unknown finite difference "
            "mode '" + finiteDifferenceMode + "'.");
    }

    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    const DDScalingSpec scaling = buildScalingSpec();

    const auto makeRecombinationConfig =
        [&](const std::vector<std::string>& models) {
            RecombinationModelConfig config =
                recombinationModelConfig(
                    models, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
            config.augerCn = cfg_.augerCn;
            config.augerCp = cfg_.augerCp;
            if (models.size() != 1 || models.front() != "none")
                config.bandToBand = cfg_.bandToBand;
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
            CoupledDDAssembler assembler(
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
                scaling,
                cfg_.carrierDiagonalFloor,
                cfg_.carrierStatistics);
            configureQuasiFermiReferences(assembler);
            return assembler;
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
    const auto matrixDifferencePair =
        [&](CoupledDDAssembler& withTerm,
            CoupledDDAssembler& withoutTerm,
            const std::string& termName) {
            const bool directElementEdgeImpact =
                termName == "sg_avalanche" &&
                detail::usesElementEdgeGssLauxAvalancheSource(
                    cfg_.impactIonization);
            if (directElementEdgeImpact) {
                const SparseMatrixd reference =
                    finiteDifferenceMode == "multiprecision_branch_resolved"
                    ? withTerm
                        .impactIonizationSourceBranchResolvedFiniteDifferenceJacobian(
                            x, bcs, finiteDifferenceStep)
                    : withTerm.impactIonizationSourceFiniteDifferenceJacobian(
                        x, bcs, finiteDifferenceStep);
                return std::pair<SparseMatrixd, SparseMatrixd>{
                    withTerm.impactIonizationSourceJacobian(x, bcs),
                    reference};
            }
            SparseMatrixd analytic =
                withTerm.assembleJacobian(x, bcs) -
                withoutTerm.assembleJacobian(x, bcs);
            const int unknowns = static_cast<int>(x.size());
            auto diagnosticResidual =
                [&](CoupledDDAssembler& assembler,
                    const VectorXd& values,
                    const std::string& selectedTerm) {
                    if (selectedTerm == "sg_avalanche")
                        return assembler.impactIonizationSourceResidual(values, bcs);
                    const auto terms =
                        assembler.carrierContinuityTermDiagnostics(values, bcs);
                    VectorXd residual = VectorXd::Zero(unknowns);
                    for (int node = 0; node < N; ++node) {
                        const auto& term = terms[static_cast<std::size_t>(node)];
                        if (selectedTerm == "sg_avalanche") {
                            residual(N + node) = term.electronImpact;
                            residual(2 * N + node) = term.holeImpact;
                        } else if (selectedTerm == "srh_auger") {
                            residual(N + node) = term.electronRecombination;
                            residual(2 * N + node) = term.holeRecombination;
                        } else {
                            residual(N + node) = term.electronGauge;
                            residual(2 * N + node) = term.holeGauge;
                        }
                    }
                    return residual;
                };
            auto centralDifference =
                [&](auto&& residualEvaluator, Real relativeStep) {
                    std::vector<Eigen::Triplet<double>> triplets;
                    triplets.reserve(static_cast<std::size_t>(unknowns) * 7);
                    for (int col = 0; col < unknowns; ++col) {
                        const Real step = relativeStep *
                            std::max<Real>(1.0, std::abs(x(col)));
                        VectorXd plus = x;
                        VectorXd minus = x;
                        plus(col) += step;
                        minus(col) -= step;
                        // Materialize both residuals before subtracting them.
                        // Some evaluators combine temporary Eigen vectors; keeping
                        // their results alive here prevents a returned lazy
                        // expression from retaining dangling references.
                        const VectorXd plusResidual = residualEvaluator(plus);
                        const VectorXd minusResidual = residualEvaluator(minus);
                        const VectorXd derivative =
                            (plusResidual - minusResidual) / (2.0 * step);
                        for (int row = 0; row < unknowns; ++row) {
                            if (derivative(row) != 0.0)
                                triplets.emplace_back(row, col, derivative(row));
                        }
                    }
                    SparseMatrixd matrix(unknowns, unknowns);
                    matrix.setFromTriplets(triplets.begin(), triplets.end());
                    return matrix;
                };
            auto termResidual = [&](const VectorXd& values) {
                return diagnosticResidual(withTerm, values, termName);
            };
            SparseMatrixd finiteDifference =
                centralDifference(termResidual, finiteDifferenceStep);
            auto gaugeResidual = [&](const VectorXd& values) -> VectorXd {
                VectorXd residual =
                    diagnosticResidual(withTerm, values, "gauge");
                residual -= diagnosticResidual(
                    withoutTerm, values, "gauge");
                return residual;
            };
            analytic -= centralDifference(gaugeResidual, finiteDifferenceStep);
            return std::pair<SparseMatrixd, SparseMatrixd>{
                analytic, finiteDifference};
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
    std::string impactConfigurationFingerprint;
    std::string impactActiveBranchFingerprint;
    if (needsBase)
        base = matrixPair(baseAssembler);
    if (needsRecombination) {
        CoupledDDAssembler recombinationAssembler =
            makeAssembler(recombinationConfig, noImpactConfig);
        withRecombination = matrixDifferencePair(
            recombinationAssembler, baseAssembler, "srh_auger");
    }
    if (needsImpact) {
        CoupledDDAssembler impactAssembler =
            makeAssembler(noRecombinationConfig, cfg_.impactIonization);
        impactConfigurationFingerprint =
            impactAssembler.impactIonizationConfigurationFingerprint();
        impactActiveBranchFingerprint =
            impactAssembler.impactIonizationActiveBranchFingerprint(x);
        withImpact = matrixDifferencePair(
            impactAssembler, baseAssembler, "sg_avalanche");
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
                withRecombination->first,
                withRecombination->second,
                jacobianAuditRows(block, N)));
        } else if (block == "sg_avalanche") {
            NewtonJacobianBlockAuditRow row = jacobianAuditRow(
                block,
                withImpact->first,
                withImpact->second,
                jacobianAuditRows(block, N));
            row.configurationFingerprint = impactConfigurationFingerprint;
            row.activeBranchFingerprint = impactActiveBranchFingerprint;
            rows.push_back(std::move(row));
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
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    return assembler.sgEdgeFluxDiagnostics(x, bcs);
}

std::vector<CoupledDDTransportEdgeJacobianDiagnostic>
NewtonSolver::evaluateTransportEdgeJacobianDiagnostics(
    const DDSolution& state,
    Real physicalFiniteDifferenceStep_V) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, cfg_.mobility, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, buildScalingSpec(), cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    std::vector<CoupledDDTransportEdgeJacobianDiagnostic> records =
        assembler.transportEdgeJacobianDiagnostics(
            x, bcs, physicalFiniteDifferenceStep_V);
    const VectorXd rowWeights = continuityRowWeights(
        assembler, x, bcs, cfg_.continuityRowScaling);
    const int N = static_cast<int>(mesh_.numNodes());
    for (auto& record : records) {
        const int offset = record.carrier == "electron" ? N : 2 * N;
        record.continuityRowWeight =
            rowWeights(offset + static_cast<int>(record.rowNode));
        record.solverProductionEdgeDerivative =
            record.continuityRowWeight *
            record.contactEliminatedProductionEdgeDerivative;
        record.solverContactIdentityEntry =
            record.continuityRowWeight * record.contactIdentityEntry;
    }
    return records;
}

NewtonResult NewtonSolver::solve() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    return solve(buildInitialGuess(assembler, bcs));
}

NewtonPoissonBlockInitialization NewtonSolver::buildPoissonBlockInitialization() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);

    NewtonPoissonBlockInitialization out;
    out.coldInitial = buildInitialGuess(assembler, bcs);
    const NewtonBlockStepEvaluation step =
        evaluateBlockStep(out.coldInitial, "poisson_only");
    out.poissonBlockInitial = step.trialSolution;
    out.coldBlockResiduals = step.residual.blockNorms;
    out.poissonBlockResiduals = step.trialResidual.blockNorms;
    out.rawStepNorm = step.rawStepNorm;
    out.stepNorm = step.stepNorm;
    return out;
}

NewtonResult NewtonSolver::solvePoissonOnly(const DDSolution& initial) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, cfg_.mobility, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, buildScalingSpec(), cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureQuasiFermiReferences(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const int N = static_cast<int>(mesh_.numNodes());
    if (initial.psi.size() != N || initial.phin.size() != N ||
        initial.phip.size() != N || initial.n.size() != N ||
        initial.p.size() != N) {
        throw std::invalid_argument(
            "NewtonSolver::solvePoissonOnly: initial state size does not match mesh.");
    }

    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    struct ContactAnchor {
        Real x = 0.0;
        Real y = 0.0;
        Real bias = 0.0;
        Real meanNetDoping = 0.0;
    };
    std::vector<ContactAnchor> electronAnchors;
    std::vector<ContactAnchor> holeAnchors;
    for (const Contact& contact : mesh_.contacts()) {
        const auto biasIt = contactBiases_.find(contact.name);
        if (biasIt == contactBiases_.end() || contact.node_ids.empty())
            continue;
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        ContactAnchor anchor;
        for (Index node : contact.node_ids) {
            const Node& geometry = mesh_.getNode(node);
            anchor.x += geometry.x;
            anchor.y += geometry.y;
            anchor.meanNetDoping += doping_.netDoping(node);
        }
        const Real inverseCount = 1.0 / static_cast<Real>(contact.node_ids.size());
        anchor.x *= inverseCount;
        anchor.y *= inverseCount;
        anchor.meanNetDoping *= inverseCount;
        anchor.bias = biasIt->second;
        if (anchor.meanNetDoping > 0.0)
            electronAnchors.push_back(anchor);
        if (anchor.meanNetDoping < 0.0)
            holeAnchors.push_back(anchor);
    }
    const auto nearestBias = [&](Index node,
                                 const std::vector<ContactAnchor>& anchors,
                                 Real fallback) {
        if (anchors.empty())
            return fallback;
        const Node& geometry = mesh_.getNode(node);
        Real bestDistance = std::numeric_limits<Real>::infinity();
        Real bestBias = fallback;
        for (const ContactAnchor& anchor : anchors) {
            const Real dx = geometry.x - anchor.x;
            const Real dy = geometry.y - anchor.y;
            const Real distance = dx * dx + dy * dy;
            if (distance < bestDistance) {
                bestDistance = distance;
                bestBias = anchor.bias;
            }
        }
        return bestBias;
    };

    VectorXd psi = initial.psi;
    VectorXd phin = initial.phin;
    VectorXd phip = initial.phip;
    for (int i = 0; i < N; ++i) {
        const Index node = static_cast<Index>(i);
        phin(i) = nearestBias(node, electronAnchors, initial.phin(i));
        phip(i) = nearestBias(node, holeAnchors, initial.phip(i));
    }
    for (const auto& [node, value] : bcs.psi)
        psi(static_cast<int>(node)) = value * potentialScale;
    for (const auto& [node, value] : bcs.phin)
        phin(static_cast<int>(node)) = value * potentialScale;
    for (const auto& [node, value] : bcs.phip)
        phip(static_cast<int>(node)) = value * potentialScale;

    VectorXd x = assembler.pack({
        psi / potentialScale,
        phin / potentialScale,
        phip / potentialScale});
    VectorXd residual = assembler.residual(x, bcs);
    const auto poissonNorm = [N](const VectorXd& value) {
        return value.head(N).norm();
    };
    const Real initialNorm = poissonNorm(residual);

    NewtonResult result;
    result.initialResidualNorm = initialNorm;
    auto finish = [&](bool converged,
                      const std::string& reason,
                      int iterations,
                      const VectorXd& state,
                      const VectorXd& finalResidual) {
        result.converged = converged;
        result.iters = iterations;
        result.finalResidualNorm = poissonNorm(finalResidual);
        result.finalBlockNorms = blockResidualInfo(finalResidual, mesh_.numNodes());
        result.solution = makeSolution(assembler, state, iterations);
        if (converged) {
            result.convergenceReason = reason;
        } else {
            result.failureDiagnostics.failureReason = reason;
            result.failureDiagnostics.failedIteration = iterations;
            result.failureDiagnostics.residualNorm = result.finalResidualNorm;
            result.failureDiagnostics.blockResiduals = result.finalBlockNorms;
        }
    };

    if (initialNorm <= cfg_.abstol) {
        finish(true, "poisson_only_contact_basin_initial_abstol", 0, x, residual);
        return result;
    }

    LineSearchConfig lineSearchConfig;
    lineSearchConfig.enabled = cfg_.lineSearch;
    lineSearchConfig.initialDamping = cfg_.dampingFactor;
    lineSearchConfig.recordHistory = cfg_.diagnostics;
    BacktrackingLineSearch lineSearch(lineSearchConfig);
    LinearSolver linearSolver;
    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        const SparseMatrixd jacobian = cfg_.jacobian == "finite_difference"
            ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
            : assembler.assembleJacobian(x, bcs);
        const SparseMatrixd poissonBlock = sparseBlock(jacobian, 0, 0, N, N);
        VectorXd step = VectorXd::Zero(3 * N);
        try {
            step.head(N) = linearSolver.solve(poissonBlock, -residual.head(N));
        } catch (const std::runtime_error&) {
            finish(false, "poisson_only_linear_solve_failed", iter - 1, x, residual);
            return result;
        }
        applyConfiguredStepCaps(step, cfg_, N, potentialScale, doping_);
        const Real rawStepNorm = step.norm();
        auto searched = lineSearch.search(
            x, step, residual,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            },
            poissonNorm);
        if (!searched.accepted) {
            finish(false, "poisson_only_line_search_rejected", iter - 1, x, residual);
            result.failureDiagnostics.stepNorm = rawStepNorm;
            return result;
        }
        x = std::move(searched.x);
        residual = std::move(searched.residual);
        const Real norm = poissonNorm(residual);
        NewtonIterationInfo info;
        info.iter = iter;
        info.residualNorm = norm;
        info.relativeResidualNorm = ResidualNorm::relative(norm, initialNorm);
        info.rawStepNorm = rawStepNorm;
        info.stepNorm = searched.damping * rawStepNorm;
        info.dampingFactor = searched.damping;
        info.lineSearchAttempts = searched.attempts;
        info.lineSearchAccepted = true;
        info.event = "poisson_only_contact_basin_iteration";
        info.blockResiduals = blockResidualInfo(residual, mesh_.numNodes());
        result.history.push_back(info);
        result.trace.push_back(std::move(info));
        if (norm <= cfg_.abstol ||
            ResidualNorm::relative(norm, initialNorm) <= cfg_.reltol) {
            finish(true, "poisson_only_contact_basin_converged", iter, x, residual);
            return result;
        }
    }
    finish(false, "poisson_only_max_iterations", cfg_.maxIter, x, residual);
    return result;
}

NewtonResult NewtonSolver::solve(const DDSolution& initial) const
{
    if (cfg_.electronQuantumPotential.enabled) {
        const int nodeCount = static_cast<int>(mesh_.numNodes());
        VectorXd quantumPotential =
            initial.electronQuantumPotential.size() == nodeCount
            ? initial.electronQuantumPotential
            : VectorXd::Zero(nodeCount);
        VectorXd quantumPotentialLike =
            initial.electronQuantumPotentialLike.size() == nodeCount
            ? initial.electronQuantumPotentialLike
            : VectorXd{};
        const bool continuousLambdaState =
            cfg_.electronQuantumPotential.globalDiscretization ==
                "gss_density_fitted" ||
            cfg_.electronQuantumPotential.globalDiscretization ==
                "p1_lambda_direct";
        if (continuousLambdaState)
            quantumPotentialLike.resize(0);
        if (cfg_.electronQuantumPotential.couplingMode == "frozen") {
            NewtonResult result = solveClassicalWithFrozenElectronQuantumPotential(
                initial, quantumPotential);
            result.solution.electronQuantumPotential = quantumPotential;
            result.solution.electronQuantumPotentialLike = quantumPotentialLike;
            if (result.converged)
                result.convergenceReason += "+frozen_electron_density_gradient";
            return result;
        }
        DDSolution outerState = initial;
        NewtonResult last;
        std::vector<DensityGradientOuterIterationInfo> outerHistory;
        VectorXd previousRawUpdate;
        Real previousRelaxation = cfg_.electronQuantumPotential.outerRelaxation;
        const bool potentialBased =
            cfg_.electronQuantumPotential.formulation == "potential_based";
        const Real Vt = thermalVoltage(cfg_.temperature_K);
        std::vector<bool> activeNodes(static_cast<std::size_t>(nodeCount), false);
        std::vector<bool> transportNodes(static_cast<std::size_t>(nodeCount), false);
        std::vector<bool> nodeOwned(static_cast<std::size_t>(nodeCount), false);
        std::vector<bool> nodeTransportOwned(
            static_cast<std::size_t>(nodeCount), false);
        std::vector<bool> cellTransport(
            static_cast<std::size_t>(mesh_.numCells()), false);
        std::vector<Real> cellElectronAffinity_eV(
            static_cast<std::size_t>(mesh_.numCells()), 0.0);
        std::vector<Real> cellQuantumGamma(
            static_cast<std::size_t>(mesh_.numCells()), 0.0);
        std::vector<Real> cellQuantumCoefficientMassRatio(
            static_cast<std::size_t>(mesh_.numCells()), 0.0);
        VectorXd materialCoefficient = VectorXd::Zero(nodeCount);
        VectorXd materialBandDrive = VectorXd::Zero(nodeCount);
        std::vector<Real> cellDosMassDrive_V(
            static_cast<std::size_t>(mesh_.numCells()), 0.0);
        std::unordered_map<Index, Real> quantumDirichlet;
        for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
            const Cell& cell = mesh_.getCell(cellId);
            const Region& region = mesh_.getRegion(cell.region_id);
            const Material material = matdb_.getMaterial(region.material);
            const bool transport = material.ni > 0.0 && material.mun > 0.0;
            cellTransport[cellId] = transport;
            cellElectronAffinity_eV[cellId] =
                material.electron_affinity_eV.value_or(0.0);
            const bool active = transport ||
                cfg_.electronQuantumPotential.includeInsulators;
            if (!active)
                continue;
            const Real gamma = material.electron_quantum_gamma.value_or(
                transport ? cfg_.electronQuantumPotential.gamma
                          : cfg_.electronQuantumPotential.insulatorGamma);
            const Real dosMassRatio =
                material.electron_quantum_dos_mass_ratio.value_or(
                    transport
                        ? cfg_.electronQuantumPotential.effectiveMassRatio
                        : cfg_.electronQuantumPotential.
                            insulatorEffectiveMassRatio);
            const Real coefficientMassRatio =
                material.electron_quantum_coefficient_mass_ratio.value_or(
                    transport
                        ? cfg_.electronQuantumPotential.coefficientMassRatio.
                            value_or(dosMassRatio)
                        : cfg_.electronQuantumPotential.
                            insulatorCoefficientMassRatio.value_or(dosMassRatio));
            cellQuantumGamma[cellId] = gamma;
            cellQuantumCoefficientMassRatio[cellId] = coefficientMassRatio;
            const Real coefficient = densityGradientCoefficientVm2(
                gamma, coefficientMassRatio);
            // Ec/q = -psi-chi.  For xi=eta=1 (semiconductor), the
            // Eq. 231 drive is psi-phin+chi.  In insulators xi=eta=0;
            // the explicit (eta-1)q*grad(phi) term cancels Ec's
            // electrostatic contribution, leaving chi plus the DOS term.
            const Real dosMassDrive =
                1.5 * Vt * std::log(dosMassRatio);
            cellDosMassDrive_V[cellId] = dosMassDrive;
            for (Index node : cell.node_ids) {
                activeNodes[node] = true;
                transportNodes[node] = transportNodes[node] || transport;
                if (!nodeOwned[node] || (transport && !nodeTransportOwned[node])) {
                    materialCoefficient(static_cast<int>(node)) = coefficient;
                    materialBandDrive(static_cast<int>(node)) =
                        material.electron_affinity_eV.value_or(0.0) +
                        dosMassDrive;
                    nodeOwned[node] = true;
                    nodeTransportOwned[node] = transport;
                }
            }
        }
        for (const Contact& contact : mesh_.contacts()) {
            const auto spec = contactSpecs_.find(contact.name);
            if (spec != contactSpecs_.end() &&
                spec->second.type == ContactType::MetalGate &&
                !cfg_.electronQuantumPotential.includeInsulators)
                continue;
            for (Index node : contact.node_ids) {
                if (activeNodes[node])
                    quantumDirichlet[node] = 0.0;
            }
        }
        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)])
                quantumPotential(i) = 0.0;
        }

        struct QuantumStepInterfaceDescriptor {
            Index edgeId = 0;
            Index solvedCellId = 0;
            Real solvedAffinity_eV = 0.0;
            Real unsolvedAffinity_eV = 0.0;
        };
        std::vector<QuantumStepInterfaceDescriptor> stepInterfaceDescriptors;
        if (cfg_.electronQuantumPotential.interfaceBoundary == "sentaurus_step") {
            if (cfg_.electronQuantumPotential.includeInsulators) {
                throw std::invalid_argument(
                    "electron_quantum_potential sentaurus_step requires "
                    "include_insulators=false.");
            }
            const auto edgeCells = detail::buildEdgeCellMap(mesh_);
            for (Index edgeId = 0; edgeId < mesh_.numEdges(); ++edgeId) {
                if (edgeCells[edgeId].size() != 2)
                    continue;
                const Index cell0 = edgeCells[edgeId][0];
                const Index cell1 = edgeCells[edgeId][1];
                if (cellTransport[cell0] == cellTransport[cell1])
                    continue;
                const Index solvedCell = cellTransport[cell0] ? cell0 : cell1;
                const Index unsolvedCell = cellTransport[cell0] ? cell1 : cell0;
                stepInterfaceDescriptors.push_back({
                    edgeId,
                    solvedCell,
                    cellElectronAffinity_eV[solvedCell],
                    cellElectronAffinity_eV[unsolvedCell],
                });
            }
        }
        const auto quantumBgn = makeBandgapNarrowingModel(cfg_.bandgapNarrowing);
        // Sentaurus OldSlotboom evaluates the Eq. 231 band drive from total
        // ionized impurity.  Carrier-density dominance belongs to Vela's
        // generic BGN model but must not leak into this quantum equation.
        const auto quantumBandgapNarrowing = [&](Index node, bool transport) {
            return transport
                ? quantumBgn->deltaEg(
                    doping_.totalImpurity(node), 0.0, 0.0)
                : Real{0.0};
        };

        for (int outer = 1;
             outer <= cfg_.electronQuantumPotential.outerMaxIterations; ++outer) {
            last = solveClassicalWithFrozenElectronQuantumPotential(
                outerState, quantumPotential);
            if (!last.converged) {
                last.electronQuantumOuterHistory = std::move(outerHistory);
                return last;
            }

            VectorXd classicalDensity(nodeCount);
            VectorXd drivingPotential = VectorXd::Zero(nodeCount);
            if (potentialBased) {
                // For default xi=eta=1 and Boltzmann statistics, Sentaurus'
                // potential formula has auxiliary variable proportional to
                // exp((psi-phin-Lambda)/(2*Vt)).  The normalization constant
                // drops out of the homogeneous DG equation.
                Real referencePotential = -std::numeric_limits<Real>::infinity();
                for (int i = 0; i < nodeCount; ++i) {
                    if (!activeNodes[static_cast<std::size_t>(i)]) {
                        drivingPotential(i) = 0.0;
                        continue;
                    }
                    drivingPotential(i) = materialBandDrive(i);
                    if (transportNodes[static_cast<std::size_t>(i)]) {
                        drivingPotential(i) +=
                            last.solution.psi(i) - last.solution.phin(i);
                    }
                    referencePotential = std::max(referencePotential, drivingPotential(i));
                }
                for (int i = 0; i < nodeCount; ++i) {
                    classicalDensity(i) = activeNodes[static_cast<std::size_t>(i)]
                        ? std::exp(std::clamp(
                            (drivingPotential(i) - referencePotential) / Vt,
                            Real{-700.0}, Real{0.0}))
                        : 0.0;
                }
            } else {
                // Density-based audit form: recover n_classical from the
                // converged density relation before the nonlinear DG solve.
                for (int i = 0; i < nodeCount; ++i) {
                    classicalDensity(i) = last.solution.n(i) * std::exp(std::clamp(
                        quantumPotential(i) / Vt, Real{-700.0}, Real{700.0}));
                }
            }
            std::vector<DensityGradientStepBoundary> stepBoundaries;
            stepBoundaries.reserve(stepInterfaceDescriptors.size());
            for (const auto& descriptor : stepInterfaceDescriptors) {
                const Edge& edge = mesh_.getEdge(descriptor.edgeId);
                const Cell& solvedCell = mesh_.getCell(descriptor.solvedCellId);
                bool gradientValid = false;
                Real solvedCellArea = 0.0;
                const Point2 drivingGradient = detail::cellScalarGradient(
                    mesh_, solvedCell,
                    [&](Index node) {
                        return drivingPotential(static_cast<int>(node));
                    },
                    gradientValid, solvedCellArea);
                const Node& p0 = mesh_.getNode(edge.n0);
                const Node& p1 = mesh_.getNode(edge.n1);
                Point2 centroid = Point2::Zero();
                for (Index node : solvedCell.node_ids) {
                    const Node& point = mesh_.getNode(node);
                    centroid += Point2{point.x, point.y};
                }
                centroid /= static_cast<Real>(solvedCell.node_ids.size());
                const Point2 midpoint{
                    0.5 * (p0.x + p1.x), 0.5 * (p0.y + p1.y)};
                const Point2 outwardDelta = midpoint - centroid;
                Point2 outwardNormal = Point2::Zero();
                if (outwardDelta.norm() > 0.0)
                    outwardNormal = outwardDelta / outwardDelta.norm();
                const Real normalDriveGradient = gradientValid
                    ? drivingGradient.dot(outwardNormal) /
                        cfg_.inputScaling.unitSystem().lengthMPerInternal()
                    : 0.0;
                const auto barrierAt = [&](Index node) {
                    const Real narrowing = quantumBandgapNarrowing(node, true);
                    return descriptor.solvedAffinity_eV -
                        descriptor.unsolvedAffinity_eV +
                        cfg_.electronQuantumPotential.
                            conductionBandNarrowingFraction * narrowing;
                };
                stepBoundaries.push_back({
                    descriptor.edgeId,
                    barrierAt(edge.n0),
                    barrierAt(edge.n1),
                    cfg_.electronQuantumPotential.insulatorEffectiveMassRatio,
                    cfg_.electronQuantumPotential.gamma,
                    cfg_.electronQuantumPotential.theta,
                    1.0,
                    normalDriveGradient,
                });
            }
            DensityGradientQuantumPotentialResult quantum;
            if (cfg_.electronQuantumPotential.includeInsulators) {
                const DDSolution& densityGradientState =
                    cfg_.electronQuantumPotential.residualDiagnosticUseInitialState
                    ? initial : last.solution;
                std::vector<DensityGradientCellMaterial> cellMaterials;
                cellMaterials.reserve(mesh_.numCells());
                VectorXd nodeOutputShift = VectorXd::Zero(nodeCount);
                for (int i = 0; i < nodeCount; ++i) {
                    const Real narrowing = quantumBandgapNarrowing(
                        static_cast<Index>(i),
                        transportNodes[static_cast<std::size_t>(i)]);
                    nodeOutputShift(i) =
                        materialBandDrive(i) + densityGradientState.psi(i) +
                        cfg_.electronQuantumPotential.
                            conductionBandNarrowingFraction * narrowing;
                }
                if (!continuousLambdaState &&
                    quantumPotentialLike.size() != nodeCount)
                    quantumPotentialLike = quantumPotential - nodeOutputShift;
                for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
                    const Cell& cell = mesh_.getCell(cellId);
                    if (cell.node_ids.size() != 3)
                        continue;
                    const bool transport = cellTransport[cellId];
                    DensityGradientCellMaterial data;
                    data.cellId = cellId;
                    data.isTransport = transport;
                    data.coefficientVm2 = densityGradientCoefficientVm2(
                        cellQuantumGamma[cellId],
                        cellQuantumCoefficientMassRatio[cellId]);
                    for (int local = 0; local < 3; ++local) {
                        const int node = static_cast<int>(cell.node_ids[local]);
                        const Real narrowing = quantumBandgapNarrowing(
                            cell.node_ids[local], transport);
                        // Phi/q is the continuous Eq. 231 primary unknown,
                        // whereas Lambda = Phi/q - Ec/q - Phi_m/q is a
                        // region-side quantity at a material interface.  The
                        // merged nodeOutputShift is only the trace convention
                        // used for restart/output Lambda; the reaction term
                        // must reconstruct Lambda with this cell's affinity
                        // and DOS mass.
                        data.materialBandDrive_V[local] =
                            densityGradientState.psi(node) +
                            cellElectronAffinity_eV[cellId] +
                            cellDosMassDrive_V[cellId] +
                            cfg_.electronQuantumPotential.
                                conductionBandNarrowingFraction * narrowing;
                        data.dynamicDrivingPotential_V[local] = transport
                            ? -densityGradientState.phin(node)
                            : -densityGradientState.psi(node);
                        data.initialLambda_V[local] = quantumPotential(node);
                    }
                    cellMaterials.push_back(data);
                }
                quantum = solveElectronDensityGradientPotentialLikeGlobal(
                    mesh_, cellMaterials, nodeOutputShift, activeNodes,
                    quantumDirichlet, Vt, cfg_.inputScaling.unitSystem(),
                    cfg_.electronQuantumPotential, quantumPotential,
                    quantumPotentialLike);
            } else {
                quantum = solveElectronDensityGradientPotential(
                    mesh_, classicalDensity, materialCoefficient,
                    activeNodes, quantumDirichlet, stepBoundaries,
                    Vt,
                    cfg_.inputScaling.unitSystem(), cfg_.electronQuantumPotential,
                    quantumPotential);
            }
            VectorXd rawUpdate = quantum.potential_V - quantumPotential;
            for (int i = 0; i < nodeCount; ++i) {
                if (!activeNodes[static_cast<std::size_t>(i)])
                    rawUpdate(i) = 0.0;
            }
            const Real change = rawUpdate.lpNorm<Eigen::Infinity>();
            Real relaxation = cfg_.electronQuantumPotential.outerRelaxation;
            if (cfg_.electronQuantumPotential.outerAcceleration == "aitken" &&
                previousRawUpdate.size() == rawUpdate.size()) {
                const VectorXd residualDelta = rawUpdate - previousRawUpdate;
                const Real denominator = residualDelta.squaredNorm();
                if (denominator > std::numeric_limits<Real>::min()) {
                    relaxation = std::clamp(
                        -previousRelaxation *
                            previousRawUpdate.dot(residualDelta) / denominator,
                        cfg_.electronQuantumPotential.outerRelaxationMin,
                        cfg_.electronQuantumPotential.outerRelaxationMax);
                }
            }
            const VectorXd appliedUpdate = relaxation * rawUpdate;
            DensityGradientOuterIterationInfo outerInfo;
            outerInfo.iteration = outer;
            outerInfo.innerIterations = quantum.iterations;
            outerInfo.innerConverged = quantum.converged;
            outerInfo.innerResidualInfinityNorm = quantum.residualInfinityNorm;
            outerInfo.innerLastUpdateInfinityNorm_V =
                quantum.lastUpdateInfinityNorm_V;
            outerInfo.innerMaxUpdateNode = quantum.maxUpdateNode;
            outerInfo.innerMaxUpdateNodeValue_V = quantum.maxUpdateNodeValue_V;
            outerInfo.rawChangeInfinityNorm_V = change;
            outerInfo.appliedChangeInfinityNorm_V =
                appliedUpdate.lpNorm<Eigen::Infinity>();
            outerInfo.potentialInfinityNorm_V =
                quantum.potential_V.lpNorm<Eigen::Infinity>();
            outerInfo.relaxation = relaxation;
            outerHistory.push_back(outerInfo);
            if (cfg_.verbose || runtimeLogEnabled()) {
                emitVerboseLine(
                    "Density-gradient outer iter " + std::to_string(outer) +
                    " inner_iters=" + std::to_string(quantum.iterations) +
                    " inner_converged=" + std::to_string(quantum.converged ? 1 : 0) +
                    " inner_residual=" + std::to_string(quantum.residualInfinityNorm) +
                    " inner_last_update_V=" +
                        std::to_string(quantum.lastUpdateInfinityNorm_V) +
                    " inner_max_update_node=" +
                        std::to_string(quantum.maxUpdateNode) +
                    " inner_max_update_signed_V=" +
                        std::to_string(quantum.maxUpdateNodeValue_V) +
                    " raw_change_V=" + std::to_string(change) +
                    " applied_change_V=" +
                        std::to_string(outerInfo.appliedChangeInfinityNorm_V) +
                    " relaxation=" + std::to_string(relaxation) +
                    " potential_max_V=" +
                        std::to_string(outerInfo.potentialInfinityNorm_V));
            }
            const Real scale = std::max<Real>(
                quantum.potential_V.lpNorm<Eigen::Infinity>(), 1.0);
            const Real outerAbsoluteTolerance_V =
                cfg_.electronQuantumPotential.outerAbsoluteTolerance_V.
                    value_or(cfg_.electronQuantumPotential.absoluteTolerance_V);
            const bool outerConverged = quantum.converged &&
                change <= outerAbsoluteTolerance_V +
                    cfg_.electronQuantumPotential.relativeTolerance * scale;
            quantumPotential = outerConverged
                ? std::move(quantum.potential_V)
                : quantumPotential + appliedUpdate;
            if (quantum.potentialLike_V.size() == nodeCount) {
                // The potential-like state is the primary global Eq. 231
                // unknown. Apply the same outer relaxation used for Lambda so
                // the restart remains consistent with the accepted iterate.
                if (quantumPotentialLike.size() != nodeCount)
                    quantumPotentialLike = quantum.potentialLike_V;
                else
                    quantumPotentialLike += relaxation *
                        (quantum.potentialLike_V - quantumPotentialLike);
            } else if (continuousLambdaState) {
                quantumPotentialLike.resize(0);
            }
            previousRawUpdate = rawUpdate;
            previousRelaxation = relaxation;
            outerState = last.solution;
            outerState.electronQuantumPotential = quantumPotential;
            outerState.electronQuantumPotentialLike = quantumPotentialLike;
            if (outerConverged) {
                // One final frozen solve makes the returned carrier state
                // consistent with the converged quantum potential.
                last = solveClassicalWithFrozenElectronQuantumPotential(
                    outerState, quantumPotential);
                last.electronQuantumOuterHistory = std::move(outerHistory);
                last.solution.electronQuantumPotential = quantumPotential;
                last.solution.electronQuantumPotentialLike = quantumPotentialLike;
                if (last.converged)
                    last.convergenceReason += "+electron_density_gradient";
                return last;
            }
        }
        last.solution.electronQuantumPotential = quantumPotential;
        last.solution.electronQuantumPotentialLike = quantumPotentialLike;
        last.electronQuantumOuterHistory = std::move(outerHistory);
        last.converged = false;
        last.failureDiagnostics.failureReason =
            "electron_density_gradient_max_iterations";
        return last;
    }
    return solveClassicalWithFrozenElectronQuantumPotential(
        initial, VectorXd::Zero(static_cast<int>(mesh_.numNodes())));
}

NewtonResult NewtonSolver::solveClassicalWithFrozenElectronQuantumPotential(
    const DDSolution& initial,
    const VectorXd& electronQuantumPotential_V) const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.bandToBand = cfg_.bandToBand;
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
        scaling,
        cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics,
        cfg_.electronQuantumPotential);
    assembler.setElectronQuantumPotential(electronQuantumPotential_V);
    configureQuasiFermiReferences(assembler);
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
    VectorXd r;
    {
        ScopedPerformanceTimer timer("newton.initial_residual");
        r = assembler.residual(x, bcs);
    }
    VectorXd activeRowWeights;
    {
        ScopedPerformanceTimer timer("newton.row_weights");
        activeRowWeights = continuityRowWeights(
            assembler, x, bcs, cfg_.continuityRowScaling);
    }
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
    result.finalBlockNorms = blockResidualInfo(r, mesh_.numNodes());

    auto carrierRowEval = [&](const VectorXd& state) {
        ScopedPerformanceTimer timer("newton.carrier_row_evaluation");
        if (cfg_.carrierRowConvergence.mode == "off")
            return NewtonCarrierRowConvergenceEvaluation{};
        return evaluateCarrierRowConvergence(
            assembler.carrierContinuityEquationTermDiagnostics(state, bcs),
            cfg_.carrierRowConvergence);
    };
    auto carrierRowsAcceptConvergence = [](const NewtonCarrierRowConvergenceEvaluation& evaluation) {
        return !evaluation.enforced || evaluation.satisfied;
    };
    std::vector<Index> electronContactNodes;
    std::vector<Index> holeContactNodes;
    electronContactNodes.reserve(bcs.phin.size());
    holeContactNodes.reserve(bcs.phip.size());
    for (const auto& [node, _] : bcs.phin)
        electronContactNodes.push_back(node);
    for (const auto& [node, _] : bcs.phip)
        holeContactNodes.push_back(node);
    auto globalClosureEval = [&](const VectorXd& state) {
        ScopedPerformanceTimer timer("newton.global_continuity_evaluation");
        if (cfg_.globalContinuityClosure.mode == "off")
            return NewtonGlobalContinuityClosureEvaluation{};
        return evaluateGlobalContinuityClosure(
            assembler.carrierContinuityEquationTermDiagnostics(
                state, CoupledDDBoundaryConditions{}),
            electronContactNodes,
            holeContactNodes,
            cfg_.globalContinuityClosure);
    };
    auto globalClosureAcceptsConvergence =
        [](const NewtonGlobalContinuityClosureEvaluation& evaluation) {
            return !evaluation.enforced || evaluation.satisfied;
        };
    Real activeGlobalElectronScale = cfg_.globalContinuityClosure.sourceFloor;
    Real activeGlobalHoleScale = cfg_.globalContinuityClosure.sourceFloor;
    auto globalClosureLineSearchNorm = [&](const VectorXd& residual) {
        const Real baseNorm = residualNormFn(residual);
        if (cfg_.globalContinuityClosure.mode == "off")
            return baseNorm;
        Real electronSum = 0.0;
        Real holeSum = 0.0;
        for (int i = 0; i < N; ++i) {
            const Index node = static_cast<Index>(i);
            if (bcs.phin.find(node) == bcs.phin.end())
                electronSum += residual(N + i);
            if (bcs.phip.find(node) == bcs.phip.end())
                holeSum += residual(2 * N + i);
        }
        const Real electronRatio =
            electronSum / std::max(activeGlobalElectronScale, Real{1.0e-300});
        const Real holeRatio =
            holeSum / std::max(activeGlobalHoleScale, Real{1.0e-300});
        return std::sqrt(
            baseNorm * baseNorm
            + electronRatio * electronRatio
            + holeRatio * holeRatio);
    };
    auto writeCarrierRowDiagnosticCsv = [&](const NewtonCarrierRowConvergenceEvaluation& evaluation,
                                            int iteration,
                                            const std::string& event) {
        if (!evaluation.enabled || cfg_.carrierRowConvergence.diagnosticCsvFile.empty())
            return;
        const std::filesystem::path path(cfg_.carrierRowConvergence.diagnosticCsvFile);
        if (!path.parent_path().empty())
            std::filesystem::create_directories(path.parent_path());
        const bool writeHeader = !std::filesystem::exists(path);
        std::ofstream out(path, std::ios::app);
        if (writeHeader) {
            out << "event,iteration,node_id,carrier,residual,scale,ratio,flux,srh,impact\n";
        }
        out << std::setprecision(17);
        for (const auto& row : evaluation.violations) {
            out << event << ',' << iteration << ',' << row.nodeId << ',' << row.carrier << ','
                << row.residual << ',' << row.scale << ',' << row.ratio << ','
                << row.flux << ',' << row.recombination << ',' << row.impact << '\n';
        }
    };
    auto shouldWriteCarrierRowTrace = [&](int iteration) {
        if (cfg_.carrierRowConvergence.traceCsvFile.empty() ||
            cfg_.carrierRowConvergence.traceNodes.empty())
            return false;
        if (iteration <= cfg_.carrierRowConvergence.traceFirstIterations)
            return true;
        const int every = std::max(1, cfg_.carrierRowConvergence.traceEveryIterations);
        return iteration % every == 0;
    };
    auto writeCarrierRowTraceCsv = [&](const VectorXd& state,
                                       const NewtonCarrierRowConvergenceEvaluation& evaluation,
                                       int iteration,
                                       Real residualNorm,
                                       const std::string& event) {
        if (!evaluation.enabled || !shouldWriteCarrierRowTrace(iteration))
            return;
        const std::filesystem::path path(cfg_.carrierRowConvergence.traceCsvFile);
        if (!path.parent_path().empty())
            std::filesystem::create_directories(path.parent_path());
        const bool writeHeader = !std::filesystem::exists(path);
        std::ofstream out(path, std::ios::app);
        if (writeHeader) {
            out << "event,iteration,residual_norm,node_id,psi_V,phin_V,phip_V,n_m3,p_m3,"
                << "electron_residual,electron_scale,electron_ratio,electron_flux,electron_srh,electron_impact,"
                << "hole_residual,hole_scale,hole_ratio,hole_flux,hole_srh,hole_impact\n";
        }
        const auto terms =
            assembler.carrierContinuityEquationTermDiagnostics(state, bcs);
        const VectorXd n = assembler.electronDensity(state);
        const VectorXd p = assembler.holeDensity(state);
        const CoupledDDState unpacked = assembler.unpack(state);
        auto scaleFor = [&](Real flux, Real recombination, Real impact) {
            return std::max({std::abs(flux), std::abs(recombination), std::abs(impact),
                             std::max<Real>(cfg_.carrierRowConvergence.scaleFloor, 0.0)});
        };
        out << std::setprecision(17);
        for (Index node : cfg_.carrierRowConvergence.traceNodes) {
            if (node < 0 || node >= mesh_.numNodes())
                continue;
            const auto& term = terms[static_cast<std::size_t>(node)];
            const Real eScale = scaleFor(term.electronFluxAbsSum, term.electronRecombination, term.electronImpact);
            const Real hScale = scaleFor(term.holeFluxAbsSum, term.holeRecombination, term.holeImpact);
            const Real eRatio = eScale > 0.0 ? std::abs(term.electronResidual) / eScale : 0.0;
            const Real hRatio = hScale > 0.0 ? std::abs(term.holeResidual) / hScale : 0.0;
            out << event << ',' << iteration << ',' << residualNorm << ',' << node << ','
                << unpacked.psi(static_cast<int>(node)) * potentialScale << ','
                << unpacked.phin(static_cast<int>(node)) * potentialScale << ','
                << unpacked.phip(static_cast<int>(node)) * potentialScale << ','
                << n(static_cast<int>(node)) << ',' << p(static_cast<int>(node)) << ','
                << term.electronResidual << ',' << eScale << ',' << eRatio << ','
                << term.electronFlux << ',' << term.electronRecombination << ',' << term.electronImpact << ','
                << term.holeResidual << ',' << hScale << ',' << hRatio << ','
                << term.holeFlux << ',' << term.holeRecombination << ',' << term.holeImpact << '\n';
        }
    };
    auto finishConverged = [&](const std::string& reason,
                               const VectorXd& state,
                               const VectorXd& residual,
                               int iterations,
                               Real norm,
                               const NewtonCarrierRowConvergenceEvaluation& rowEval,
                               const NewtonGlobalContinuityClosureEvaluation& globalEval) {
        result.converged = true;
        result.iters = iterations;
        result.finalResidualNorm = norm;
        result.convergenceReason = reason;
        result.finalBlockNorms = blockResidualInfo(residual, mesh_.numNodes());
        result.finalCarrierRowConvergence = rowEval;
        result.finalGlobalContinuityClosure = globalEval;
        result.solution = makeSolution(assembler, state, iterations);
        writeCarrierRowDiagnosticCsv(rowEval, iterations, reason);
        writeCarrierRowTraceCsv(state, rowEval, iterations, norm, reason);
    };

    auto retryAfterCarrierRowRecovery = [&](const VectorXd& state,
                                            int iterations,
                                            const NewtonCarrierRowConvergenceEvaluation& rowEval)
        -> std::optional<NewtonResult> {
        if (cfg_.carrierRowRecovery.mode == "off" || cfg_.carrierRowRecovery.maxAttempts <= 0 ||
            cfg_.carrierRowRecovery.maxCycles <= 0 || !rowEval.enabled || rowEval.satisfied ||
            rowEval.violations.empty()) {
            return std::nullopt;
        }
        DDSolution base = makeSolution(assembler, state, iterations);
        NewtonCarrierRowRecoveryResult recovery = recoverCarrierRowsWithGummelDensity(
            mesh_, matdb_, doping_, contactBiases_, cfg_, base, rowEval.violations,
            cfg_.carrierRowRecovery, contactSpecs_);
        if (!recovery.attempted ||
            (recovery.electronRowsUpdated == 0 && recovery.holeRowsUpdated == 0)) {
            return std::nullopt;
        }
        NewtonConfig retryCfg = cfg_;
        retryCfg.carrierRowRecovery.maxCycles = cfg_.carrierRowRecovery.maxCycles - 1;
        retryCfg.warmStart = true;
        NewtonSolver retrySolver(
            mesh_, matdb_, doping_, contactBiases_, retryCfg, fixedCharges_, sheetCharges_,
            contactSpecs_);
        NewtonResult retried = retrySolver.solve(recovery.solution);
        std::vector<NewtonIterationInfo> combinedTrace = result.trace;
        combinedTrace.insert(
            combinedTrace.end(),
            retried.trace.begin(),
            retried.trace.end());
        retried.trace = std::move(combinedTrace);
        if (retried.carrierRowRecovery.attempted) {
            recovery.electronRowsUpdated += retried.carrierRowRecovery.electronRowsUpdated;
            recovery.holeRowsUpdated += retried.carrierRowRecovery.holeRowsUpdated;
            recovery.densityPasses += retried.carrierRowRecovery.densityPasses;
            recovery.cyclesAttempted += retried.carrierRowRecovery.cyclesAttempted;
            recovery.densityConverged = recovery.densityConverged &&
                retried.carrierRowRecovery.densityConverged;
            recovery.maxDensityRelativeChange = std::max(
                recovery.maxDensityRelativeChange,
                retried.carrierRowRecovery.maxDensityRelativeChange);
            recovery.maxCarrierDensityRatio = std::max(
                recovery.maxCarrierDensityRatio,
                retried.carrierRowRecovery.maxCarrierDensityRatio);
            recovery.maxPsiDelta_V = std::max(
                recovery.maxPsiDelta_V,
                retried.carrierRowRecovery.maxPsiDelta_V);
            recovery.solution = retried.carrierRowRecovery.solution;
        }
        retried.carrierRowRecovery = recovery;
        return retried;
    };


    if (cfg_.verbose) {
        const ResidualBlockNormValue blocks = ResidualNorm::computeBlocks(r, mesh_.numNodes());
        emitVerboseLine(
            "Newton iter 0 residual=" + std::to_string(initialNorm) +
            " step=0 damping=0 blocks=(" +
            std::to_string(blocks.psi) + "," +
            std::to_string(blocks.phin) + "," +
            std::to_string(blocks.phip) + ")");
    }

    NewtonCarrierRowConvergenceEvaluation initialRowEval = carrierRowEval(x);
    NewtonGlobalContinuityClosureEvaluation initialGlobalEval = globalClosureEval(x);
    NewtonIterationInfo initialTrace;
    initialTrace.iter = 0;
    initialTrace.residualNorm = initialNorm;
    initialTrace.relativeResidualNorm = 1.0;
    initialTrace.lineSearchAccepted = true;
    initialTrace.event = "initial";
    initialTrace.carrierRowConvergence = initialRowEval;
    fillNewtonIterationResidualTrace(
        initialTrace, x, r, activeRowWeights, mesh_.numNodes(), assembler,
        cfg_.diagnostics);
    result.trace.push_back(std::move(initialTrace));
    writeCarrierRowTraceCsv(x, initialRowEval, 0, initialNorm, "initial");
    if (initialNorm <= cfg_.abstol &&
        carrierRowsAcceptConvergence(initialRowEval) &&
        globalClosureAcceptsConvergence(initialGlobalEval)) {
        finishConverged(
            "initial_abstol", x, r, 0, initialNorm,
            initialRowEval, initialGlobalEval);
        return result;
    }
    if (initialNorm <= cfg_.abstol && initialRowEval.enforced &&
        !initialRowEval.satisfied) {
        writeCarrierRowDiagnosticCsv(
            initialRowEval, 0, "carrier_row_convergence_initial_abstol_rejected");
        writeCarrierRowTraceCsv(
            x, initialRowEval, 0, initialNorm,
            "carrier_row_convergence_initial_abstol_rejected");
        if (auto recovered = retryAfterCarrierRowRecovery(x, 0, initialRowEval))
            return *recovered;
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
    NewtonGlobalContinuityClosureEvaluation activeGlobalEval =
        initialGlobalEval;

    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        ScopedPerformanceTimer iterationTimer("newton.iteration");
        activeGlobalElectronScale = std::max(
            std::abs(activeGlobalEval.electron.integratedSource),
            cfg_.globalContinuityClosure.sourceFloor);
        activeGlobalHoleScale = std::max(
            std::abs(activeGlobalEval.hole.integratedSource),
            cfg_.globalContinuityClosure.sourceFloor);
        {
            ScopedPerformanceTimer timer("newton.row_weights");
            activeRowWeights = continuityRowWeights(
                assembler, x, bcs, cfg_.continuityRowScaling);
        }
        SparseMatrixd J;
        {
            ScopedPerformanceTimer timer("newton.jacobian");
            J = (cfg_.jacobian == "finite_difference")
                ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
                : assembler.assembleJacobian(x, bcs);
            addCarrierRowRegularization(J, N, cfg_.carrierRegularizationScale);
        }
        VectorXd step;
        try {
            if (cfg_.continuityRowScaling.enabled) {
                ScopedPerformanceTimer timer("newton.linear_row_scaling");
                const SparseMatrixd scaledJ = leftScaleRows(J, activeRowWeights);
                step = linearSolver.solve(
                    scaledJ,
                    -r.cwiseProduct(activeRowWeights));
            } else {
                step = linearSolver.solve(J, -r);
            }
        } catch (const std::runtime_error&) {
            NewtonIterationInfo failedTrace;
            failedTrace.iter = iter;
            failedTrace.residualNorm = residualNormFn(r);
            failedTrace.relativeResidualNorm =
                ResidualNorm::relative(failedTrace.residualNorm, initialNorm);
            failedTrace.lineSearchAccepted = false;
            failedTrace.event = "linear_solve_failed";
            failedTrace.carrierRowConvergence = carrierRowEval(x);
            fillNewtonIterationResidualTrace(
                failedTrace, x, r, activeRowWeights, mesh_.numNodes(), assembler,
                cfg_.diagnostics);
            result.trace.push_back(std::move(failedTrace));
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
                emitVerboseErrorLine(
                    "Newton failed at iter " + std::to_string(iter) +
                    ": residual=" + std::to_string(residualNormFn(r)) +
                    " damping=0 step=0 (linear solve failed)");
                printFailureDiagnostics(result.failureDiagnostics);
            }
            return result;
        }
        {
            ScopedPerformanceTimer timer("newton.update_constraints");
            capConfiguredStep(step, J, r);
        }
        Real stepNorm = step.norm();

        auto ls = lineSearch.search(
            x, step, r,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            },
            globalClosureLineSearchNorm);

        if (!ls.accepted) {
            const Real stalledNorm = residualNormFn(acceptedR);
            // Effectively-solved state: the residual already sits at the
            // numerical floor, so the rejected step is only fighting noise.
            // Declaring convergence here avoids spurious failures when the
            // Newton iterate has already reached the achievable precision.
            const NewtonCarrierRowConvergenceEvaluation stalledRowEval = carrierRowEval(acceptedX);
            const NewtonGlobalContinuityClosureEvaluation stalledGlobalEval =
                globalClosureEval(acceptedX);
            const NewtonBlockResidualInfo stalledBlocks = blockResidualInfo(acceptedR, mesh_.numNodes());
            NewtonIterationInfo rejectedTrace;
            rejectedTrace.iter = iter;
            rejectedTrace.residualNorm = stalledNorm;
            rejectedTrace.relativeResidualNorm =
                ResidualNorm::relative(stalledNorm, initialNorm);
            rejectedTrace.rawStepNorm = stepNorm;
            rejectedTrace.stepNorm = 0.0;
            rejectedTrace.dampingFactor = ls.damping;
            rejectedTrace.lineSearchAttempts = ls.attempts;
            rejectedTrace.lineSearchAccepted = false;
            rejectedTrace.event = ls.failureReason.empty()
                ? "line_search_rejected"
                : ls.failureReason;
            rejectedTrace.carrierRowConvergence = stalledRowEval;
            if (cfg_.diagnostics)
                rejectedTrace.lineSearchHistory = ls.history;
            fillNewtonIterationResidualTrace(
                rejectedTrace,
                acceptedX,
                acceptedR,
                activeRowWeights,
                mesh_.numNodes(),
                assembler,
                cfg_.diagnostics);
            result.trace.push_back(std::move(rejectedTrace));
            if (stalledNorm <= stallResidualFloor &&
                carrierRowsAcceptConvergence(stalledRowEval) &&
                globalClosureAcceptsConvergence(stalledGlobalEval)) {
                finishConverged("stall_residual_floor", acceptedX, acceptedR,
                                acceptedIters, stalledNorm, stalledRowEval,
                                stalledGlobalEval);
                return result;
            }
            const DDSolution stalledSolution = makeSolution(assembler, acceptedX, acceptedIters);
            const Real stalledContactMajorityQfDrop = maxContactMajorityQuasiFermiDrop(stalledSolution);
            Real bestRejectedContactMajorityQfDrop = std::numeric_limits<Real>::infinity();
            if (ls.bestRejectedCandidate) {
                const DDSolution bestRejectedSolution = makeSolution(assembler, ls.bestRejectedX, acceptedIters);
                bestRejectedContactMajorityQfDrop = maxContactMajorityQuasiFermiDrop(bestRejectedSolution);
            }

            if (isPoissonLineSearchStall(ls, stalledBlocks, stalledNorm, stalledContactMajorityQfDrop, cfg_) &&
                carrierRowsAcceptConvergence(stalledRowEval) &&
                globalClosureAcceptsConvergence(stalledGlobalEval)) {
                finishConverged("poisson_line_search_stall_floor", acceptedX, acceptedR,
                                acceptedIters, stalledNorm, stalledRowEval,
                                stalledGlobalEval);
                return result;
            }
            if (isCarrierRowQualifiedLineSearchStall(
                    ls, stalledBlocks, stalledNorm,
                    stalledContactMajorityQfDrop, stalledRowEval, cfg_) &&
                globalClosureAcceptsConvergence(stalledGlobalEval)) {
                finishConverged(
                    "carrier_row_qualified_stall_floor",
                    acceptedX, acceptedR, acceptedIters, stalledNorm,
                    stalledRowEval, stalledGlobalEval);
                return result;
            }
            if (stalledRowEval.enforced && !stalledRowEval.satisfied) {
                writeCarrierRowDiagnosticCsv(
                    stalledRowEval, acceptedIters,
                    "carrier_row_convergence_line_search_rejected");
                writeCarrierRowTraceCsv(
                    acceptedX, stalledRowEval, acceptedIters, stalledNorm,
                    "carrier_row_convergence_line_search_rejected");
                if (auto recovered =
                        retryAfterCarrierRowRecovery(acceptedX, acceptedIters, stalledRowEval))
                    return *recovered;
            }
            result.finalResidualNorm = stalledNorm;
            result.iters = acceptedIters;
            result.finalBlockNorms = stalledBlocks;
            result.finalCarrierRowConvergence = stalledRowEval;
            result.finalGlobalContinuityClosure = stalledGlobalEval;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            const std::string rejectedFailureReason =
                (stalledGlobalEval.enforced && !stalledGlobalEval.satisfied)
                    ? std::string("global_continuity_closure_line_search_rejected")
                    : (stalledRowEval.enforced && !stalledRowEval.satisfied)
                    ? std::string("carrier_row_convergence_line_search_rejected")
                    : (ls.failureReason.empty()
                        ? std::string("line_search_rejected")
                        : ls.failureReason);
            result.failureDiagnostics = buildFailureDiagnostics(
                mesh_,
                doping_,
                assembler,
                acceptedX,
                acceptedR,
                rejectedFailureReason,
                iter,
                result.finalResidualNorm,
                stepNorm,
                ls.damping,
                ls.attempts,
                ls.failureReason,
                stalledContactMajorityQfDrop,
                bestRejectedContactMajorityQfDrop,
                std::move(ls.history));
            if (cfg_.verbose) {
                emitVerboseErrorLine(
                    "Newton failed at iter " + std::to_string(iter) +
                    ": residual=" + std::to_string(ls.residualNorm) +
                    " damping=" + std::to_string(ls.damping) +
                    " step=" + std::to_string(stepNorm) +
                    " (line search rejected step; reason=" +
                    result.failureDiagnostics.failureReason + ")");
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
        const Real residualNorm = residualNormFn(r);
        NewtonIterationInfo info;
        info.iter = iter;
        info.residualNorm = residualNorm;
        info.stepNorm = appliedStepNorm;
        info.dampingFactor = ls.damping;
        info.relativeResidualNorm = ResidualNorm::relative(residualNorm, initialNorm);
        info.rawStepNorm = stepNorm;
        info.lineSearchAttempts = ls.attempts;
        info.lineSearchAccepted = ls.accepted;
        info.event = "accepted_iteration";
        fillNewtonIterationResidualTrace(
            info, x, r, activeRowWeights, mesh_.numNodes(), assembler,
            cfg_.diagnostics);
        info.carrierRowConvergence = carrierRowEval(x);
        writeCarrierRowTraceCsv(x, info.carrierRowConvergence, iter, residualNorm, "iteration");
        if (cfg_.diagnostics)
            info.lineSearchHistory = std::move(ls.history);
        result.trace.push_back(info);
        result.history.push_back(std::move(info));
        if (cfg_.verbose) {
            const NewtonBlockResidualInfo& blocks = result.history.back().blockResiduals;
            emitVerboseLine(
                "Newton iter " + std::to_string(iter) +
                " residual=" + std::to_string(residualNorm) +
                " step=" + std::to_string(appliedStepNorm) +
                " damping=" + std::to_string(ls.damping) +
                " blocks=(" + std::to_string(blocks.psi) + "," +
                std::to_string(blocks.phin) + "," +
                std::to_string(blocks.phip) + ")");
        }

        const Real rel = result.history.back().relativeResidualNorm;
        const bool absoluteConverged = residualNorm <= cfg_.abstol;
        const bool relativeConverged = rel <= cfg_.reltol;
        const NewtonCarrierRowConvergenceEvaluation& rowEval =
            result.history.back().carrierRowConvergence;
        const NewtonGlobalContinuityClosureEvaluation globalEval =
            globalClosureEval(x);
        activeGlobalEval = globalEval;
        if ((absoluteConverged || relativeConverged) &&
            carrierRowsAcceptConvergence(rowEval) &&
            globalClosureAcceptsConvergence(globalEval)) {
            finishConverged(
                absoluteConverged ? "abstol" : "reltol",
                x,
                r,
                iter,
                residualNorm,
                rowEval,
                globalEval);
            return result;
        }
        if (cfg_.carrierRowQualifiedStallAcceptance &&
            rowEval.enforced && rowEval.satisfied &&
            residualNorm <= stallResidualFloor &&
            globalClosureAcceptsConvergence(globalEval)) {
            const DDSolution floorSolution =
                makeSolution(assembler, x, iter);
            const Real contactMajorityQfDrop =
                maxContactMajorityQuasiFermiDrop(floorSolution);
            if (cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V <= 0.0 ||
                (std::isfinite(contactMajorityQfDrop) &&
                 contactMajorityQfDrop <=
                     cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V)) {
                finishConverged(
                    "carrier_row_qualified_residual_floor", x, r, iter,
                    residualNorm, rowEval, globalEval);
                return result;
            }
        }
    }

    result.iters = acceptedIters;
    result.finalResidualNorm = residualNormFn(acceptedR);
    result.solution = makeSolution(assembler, acceptedX, acceptedIters);
    const NewtonCarrierRowConvergenceEvaluation finalRowEval = carrierRowEval(acceptedX);
    const NewtonGlobalContinuityClosureEvaluation& finalGlobalEval =
        activeGlobalEval;
    result.finalBlockNorms = blockResidualInfo(acceptedR, mesh_.numNodes());
    result.finalCarrierRowConvergence = finalRowEval;
    result.finalGlobalContinuityClosure = finalGlobalEval;
    if (result.finalResidualNorm <= stallResidualFloor &&
        carrierRowsAcceptConvergence(finalRowEval) &&
        globalClosureAcceptsConvergence(finalGlobalEval)) {
        finishConverged("max_iter_stall_residual_floor", acceptedX, acceptedR,
                        acceptedIters, result.finalResidualNorm, finalRowEval,
                        finalGlobalEval);
        return result;
    }

    result.converged = false;
    const std::string finalFailureReason =
        (finalGlobalEval.enforced && !finalGlobalEval.satisfied)
            ? std::string("global_continuity_closure")
            : (finalRowEval.enforced && !finalRowEval.satisfied)
            ? std::string("carrier_row_convergence")
            : std::string("max_iterations");
    writeCarrierRowDiagnosticCsv(finalRowEval, acceptedIters, finalFailureReason);
    writeCarrierRowTraceCsv(acceptedX, finalRowEval, acceptedIters,
                            result.finalResidualNorm, finalFailureReason);
    if (finalRowEval.enforced && !finalRowEval.satisfied) {
        if (auto recovered = retryAfterCarrierRowRecovery(acceptedX, acceptedIters, finalRowEval))
            return *recovered;
    }
    result.failureDiagnostics = buildFailureDiagnostics(
        mesh_,
        doping_,
        assembler,
        acceptedX,
        acceptedR,
        finalFailureReason,
        cfg_.maxIter,
        result.finalResidualNorm,
        result.history.empty() ? 0.0 : result.history.back().stepNorm,
        result.history.empty() ? 0.0 : result.history.back().dampingFactor,
        result.history.empty() ? 0 : result.history.back().lineSearchAttempts,
        {});
    if (cfg_.verbose) {
        emitVerboseErrorLine(
            "Newton failed after " + std::to_string(cfg_.maxIter) +
            " iterations: residual=" + std::to_string(result.finalResidualNorm) +
            " damping=" + std::to_string(
                result.history.empty() ? 0.0 : result.history.back().dampingFactor) +
            " step=" + std::to_string(
                result.history.empty() ? 0.0 : result.history.back().stepNorm));
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
    const auto started = std::chrono::steady_clock::now();
    NewtonResult result = NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve();
    if (PerformanceProfiler* profiler = activePerformanceProfiler()) {
        profiler->recordNewtonSolve(
            result.converged, result.iters, result.initialResidualNorm,
            result.finalResidualNorm,
            result.converged ? result.convergenceReason
                             : result.failureDiagnostics.failureReason,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started));
    }
    return result;
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg)
{
    const auto started = std::chrono::steady_clock::now();
    NewtonResult result =
        NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve(initial);
    if (PerformanceProfiler* profiler = activePerformanceProfiler()) {
        profiler->recordNewtonSolve(
            result.converged, result.iters, result.initialResidualNorm,
            result.finalResidualNorm,
            result.converged ? result.convergenceReason
                             : result.failureDiagnostics.failureReason,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started));
    }
    return result;
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges,
                       ContactSpecsMap contactSpecs)
{
    const auto started = std::chrono::steady_clock::now();
    NewtonResult result = NewtonSolver(
        mesh,
        matdb,
        doping,
        contactBiases,
        cfg,
        std::move(fixedCharges),
        std::move(sheetCharges),
        std::move(contactSpecs)).solve();
    if (PerformanceProfiler* profiler = activePerformanceProfiler()) {
        profiler->recordNewtonSolve(
            result.converged, result.iters, result.initialResidualNorm,
            result.finalResidualNorm,
            result.converged ? result.convergenceReason
                             : result.failureDiagnostics.failureReason,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started));
    }
    return result;
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges,
                       ContactSpecsMap contactSpecs)
{
    const auto started = std::chrono::steady_clock::now();
    NewtonResult result = NewtonSolver(
        mesh,
        matdb,
        doping,
        contactBiases,
        cfg,
        std::move(fixedCharges),
        std::move(sheetCharges),
        std::move(contactSpecs)).solve(initial);
    if (PerformanceProfiler* profiler = activePerformanceProfiler()) {
        profiler->recordNewtonSolve(
            result.converged, result.iters, result.initialResidualNorm,
            result.finalResidualNorm,
            result.converged ? result.convergenceReason
                             : result.failureDiagnostics.failureReason,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started));
    }
    return result;
}

NewtonResult runNewtonPoissonOnly(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    const DDSolution& initial,
    const NewtonConfig& cfg,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    ContactSpecsMap contactSpecs)
{
    const auto started = std::chrono::steady_clock::now();
    NewtonResult result = NewtonSolver(
        mesh,
        matdb,
        doping,
        contactBiases,
        cfg,
        std::move(fixedCharges),
        std::move(sheetCharges),
        std::move(contactSpecs)).solvePoissonOnly(initial);
    if (PerformanceProfiler* profiler = activePerformanceProfiler()) {
        profiler->recordNewtonSolve(
            result.converged,
            result.iters,
            result.initialResidualNorm,
            result.finalResidualNorm,
            result.converged ? result.convergenceReason
                             : result.failureDiagnostics.failureReason,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started));
    }
    return result;
}

} // namespace vela
