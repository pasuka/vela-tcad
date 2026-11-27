#include "vela/core/ArithmeticPrecision.h"
#include "vela/numerics/WeightedSparseRows.h"
#include "vela/numerics/StableMeritComparison.h"
#include "vela/equation/SplitDDRuntime.h"
#include "vela/equation/ExtendedPoissonResidual.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/numerics/SplitContactVoltage.h"
#include "vela/simulation/CommonRelativeLinearSolve.h"
#include "vela/simulation/PrimaryVariableCoordinates.h"
#include "vela/solver/experimental/NewtonCapture.h"
#include "vela/solver/experimental/NewtonAudit.h"
#include "vela/core/experimental/CacheAudit.h"
#include "vela/solver/LinearRefinement.h"
#include "vela/solver/ShadowConvergenceDiagnostics.h"
#include "vela/solver/SlotPopulationNumerics.h"
#include "vela/solver/PopulationWriteback.h"
#include "PopulationStepControl.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/RuntimeLog.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/LinearSolver.h"
#include "vela/solver/detail/ContinuityTermCache.h"
#include "vela/solver/detail/NewtonTrendStop.h"
#include "LinearCapture.h"
#include "vela/io/DDSolutionCsv.h"
#include "vela/io/DDSolutionState.h"
#include <nlohmann/json.hpp>
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace vela {
namespace {

VectorXd packReferencedSolution(
    CoupledDDAssembler& assembler,
    const DDSolution& state,
    const CoupledDDBoundaryConditions& bcs,
    bool checkConsistency = false)
{
    (void)bcs;
    if(assembler.usesPopulationState()) {
        if(state.packedLow.size())throw std::invalid_argument("population state cannot import split packed coordinates");
        const double scale=assembler.usesScaledState()?assembler.potentialScale():1.;
        CoupledDDState normalized;
        normalized.psi=state.normalizedPopulationPotential(scale);
        normalized.n=state.n;normalized.p=state.p;
        return assembler.pack(normalized);
    }
    if(state.packedLow.size()!=0) {
        if(!state.hasConsistentSplitPackedState())
            throw std::invalid_argument("Split checkpoint field mismatch");
        // Read-only packing must use the checkpoint's coordinate frame. A
        // terminal-nearest frame selected from its final physical fields can
        // differ from the one fixed at the beginning of the Newton solve.
        // Continuation rebasing remains explicit in restoreForContinuation().
        bool matchingFrame=true;
        for(Index i=0;i<assembler.numNodes();++i)
            matchingFrame=matchingFrame &&
                assembler.electronQuasiFermiReferenceAt(i)==state.electronQuasiFermiReferenceAt(i) &&
                assembler.holeQuasiFermiReferenceAt(i)==state.holeQuasiFermiReferenceAt(i);
        if(!matchingFrame)
            assembler.setQuasiFermiReferenceFields(state.electronQfReference,state.holeQfReference);
        auto runtime=assembler.splitRuntime();
        if(!runtime)throw std::invalid_argument("Split checkpoint requires split_dd_state");
        return runtime->restore(state);
    }
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    VectorXd packed = assembler.pack({
        state.psi / potentialScale,
        state.phin / potentialScale,
        state.phip / potentialScale});
    const int N = static_cast<int>(assembler.numNodes());
    // Subtract physical references before nondimensionalizing even for
    // six-column external states; scaled subtraction loses tiny QF drops.
    for (int i = 0; i < N; ++i) {
        packed(N + i) = static_cast<Real>((static_cast<vela::ArithmeticReal>(state.phin(i)) -
            assembler.electronQuasiFermiReferenceAt(static_cast<Index>(i))) / potentialScale);
        packed(2 * N + i) = static_cast<Real>((static_cast<vela::ArithmeticReal>(state.phip(i)) -
            assembler.holeQuasiFermiReferenceAt(static_cast<Index>(i))) / potentialScale);
    }
    const auto consistentCoordinate = [](Real physical, Real reconstructed) {
        const Real tolerance =
            32.0 * std::numeric_limits<Real>::epsilon() *
            std::max({Real{1.0}, std::abs(physical), std::abs(reconstructed)});
        return std::abs(physical - reconstructed) <= tolerance;
    };
    if (state.hasReferencedElectronQuasiFermi()) {
        for (int i = 0; i < N; ++i) {
            // Match ordinary Newton restart semantics: an intentional edit
            // of the absolute field takes precedence over stale coordinates.
            if (checkConsistency && !consistentCoordinate(
                    state.phin(i), state.electronQuasiFermiAt(i)))
                continue;
            // Cancel matching references before adding the tiny coordinate.
            // Even vela::ArithmeticReal can lose it if added to the bias first.
            const vela::ArithmeticReal increment =
                (static_cast<vela::ArithmeticReal>(
                     state.electronQuasiFermiReferenceAt(i)) -
                 static_cast<vela::ArithmeticReal>(
                     assembler.electronQuasiFermiReferenceAt(
                         static_cast<Index>(i)))) +
                static_cast<vela::ArithmeticReal>(state.phinIncrement(i));
            packed(N + i) = static_cast<Real>(
                increment / static_cast<vela::ArithmeticReal>(potentialScale));
        }
    }
    if (state.hasReferencedHoleQuasiFermi()) {
        for (int i = 0; i < N; ++i) {
            if (checkConsistency && !consistentCoordinate(
                    state.phip(i), state.holeQuasiFermiAt(i)))
                continue;
            const vela::ArithmeticReal increment =
                (static_cast<vela::ArithmeticReal>(
                     state.holeQuasiFermiReferenceAt(i)) -
                 static_cast<vela::ArithmeticReal>(
                     assembler.holeQuasiFermiReferenceAt(
                         static_cast<Index>(i)))) +
                static_cast<vela::ArithmeticReal>(state.phipIncrement(i));
            packed(2 * N + i) = static_cast<Real>(
                increment / static_cast<vela::ArithmeticReal>(potentialScale));
        }
    }
    if(state.hasConsistentPackedState() && state.packedPotentialScale_V==potentialScale) {
        packed.head(N)=state.packedState.head(N);
        for(int i=0;i<N;++i) {
            if(state.electronQuasiFermiReferenceAt(i)==assembler.electronQuasiFermiReferenceAt(i))
                packed(N+i)=state.packedState(N+i);
            if(state.holeQuasiFermiReferenceAt(i)==assembler.holeQuasiFermiReferenceAt(i))
                packed(2*N+i)=state.packedState(2*N+i);
        }
    }
    return packed;
}

void restoreElectronQuantumPotential(
    CoupledDDAssembler& assembler,
    const DDSolution& state)
{
    if (state.electronQuantumPotential.size() ==
        static_cast<int>(assembler.numNodes())) {
        assembler.setElectronQuantumPotential(state.electronQuantumPotential);
    }
}

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

std::string scientificReal(Real value)
{
    std::ostringstream out;
    out << std::scientific << std::setprecision(9) << value;
    return out.str();
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
                                 Real regularizationScale,
                                 const CoupledDDBoundaryConditions& bcs)
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
        // Boundary identity equations must stay exact. Stabilize only the
        // physical free-carrier rows, including when exact updates are used.
        if ((localRow < nodeCount && bcs.phin.contains(localRow)) ||
            (localRow >= nodeCount && bcs.phip.contains(localRow - nodeCount)))
            continue;
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

Real configuredQuasiFermiUniformStepScale(const VectorXd& step,
                                          const NewtonConfig& cfg,
                                          int nodeCount,
                                          Real potentialScale,
                                          const DopingModel& doping)
{
    VectorXd componentLimited = step;
    if (!applyConfiguredQuasiFermiStepCaps(
            componentLimited, cfg, nodeCount, potentialScale, doping)) {
        return 1.0;
    }

    Real scale = 1.0;
    for (int i = nodeCount; i < 3 * nodeCount; ++i) {
        const Real magnitude = std::abs(step(i));
        if (magnitude > 0.0) {
            scale = std::min(scale, std::abs(componentLimited(i)) / magnitude);
        }
    }
    return scale;
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

    if (cfg.quasiFermiUpdateLimitMode == "uniform_trust_region") {
        const Real scale = configuredQuasiFermiUniformStepScale(
            step, cfg, nodeCount, potentialScale, doping);
        if (scale < 1.0)
            step *= scale;
        // The complete Newton vector was scaled by one scalar, so it remains
        // a descent direction for the linearized residual.  No Poisson
        // sub-block recorrection is needed (and doing one would destroy that
        // direction-preserving property).
        return false;
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
    if (clippedQuasiFermi) {
        recorrectPoissonStepForClippedQuasiFermi(step, J, residual, nodeCount);
        // Poisson recorrection replaces the psi block after max_update was
        // applied above and can therefore reintroduce an arbitrarily large
        // full-vector update.  Reapply the public whole-step contract to the
        // recorrected direction; uniform rescaling preserves its block
        // balance and keeps globalization effective.
        if (cfg.maxUpdate > 0.0) {
            const Real maxAbsStep = step.cwiseAbs().maxCoeff();
            if (maxAbsStep > cfg.maxUpdate)
                step *= cfg.maxUpdate / maxAbsStep;
        }
    }
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
    const NewtonCarrierRowConvergenceConfig& cfg,
    std::vector<NewtonCarrierRowConvergenceViolation>* activeRowAudit)
{
    if (activeRowAudit) activeRowAudit->clear();
    NewtonCarrierRowConvergenceEvaluation evaluation;
    if (cfg.floorLimitedAction != "enforce" && cfg.floorLimitedAction != "report")
        throw std::invalid_argument("carrier_row_convergence.floor_limited_action must be enforce or report");
    evaluation.enabled = cfg.mode != "off";
    evaluation.enforced = cfg.mode == "enforce";
    evaluation.epsRow = cfg.epsRow;
    evaluation.scaleFloor = std::max<Real>(cfg.scaleFloor, 0.0);
    evaluation.floorLimitedAction = cfg.floorLimitedAction;
    if (!evaluation.enabled)
        return evaluation;

    const Real eps = cfg.epsRow;
    const Real floor = std::max<Real>(cfg.scaleFloor, 0.0);
    Real maxElectronFluxScale = 0.0;
    Real maxHoleFluxScale = 0.0;
    Real maxElectronSourceScale = 0.0;
    Real maxHoleSourceScale = 0.0;
    for (const CoupledDDCarrierTermDiagnostic& row : rows) {
        maxElectronFluxScale = std::max(maxElectronFluxScale, std::abs(row.electronFluxAbsSum));
        maxHoleFluxScale = std::max(maxHoleFluxScale, std::abs(row.holeFluxAbsSum));
        maxElectronSourceScale = std::max(maxElectronSourceScale,
            std::max(std::abs(row.electronRecombination), std::abs(row.electronImpact)));
        maxHoleSourceScale = std::max(maxHoleSourceScale,
            std::max(std::abs(row.holeRecombination), std::abs(row.holeImpact)));
    }
    auto checkCarrier = [&](const CoupledDDCarrierTermDiagnostic& row,
                            const std::string& carrier,
                            bool continuityActive,
                            Real carrierDensity,
                            Real residual,
                            Real flux,
                            Real fluxAbsSum,
                            Real recombination,
                            Real impact,
                            Real maxFluxScale,
                            Real maxSourceScale) {
        if (!continuityActive) {
            ++evaluation.ignoredRowCount;
            return;
        }
        const Real sourceScale = std::max(std::abs(recombination), std::abs(impact));
        const Real physicalScale = std::max(std::abs(fluxAbsSum), sourceScale);
        const bool floorLimited = floor > 0.0 && physicalScale <= floor;
        const Real scale = std::max({std::abs(fluxAbsSum), sourceScale, floor});
        const bool sourceQualified = sourceScale > 0.0 && scale > 0.0 &&
            sourceScale >= cfg.minSourceScale &&
            sourceScale >= cfg.minSourceScaleFraction * scale &&
            sourceScale >= cfg.minSourceGlobalFraction * maxSourceScale;
        const bool densityQualified = cfg.minCarrierDensity_m3 > 0.0 &&
            carrierDensity >= cfg.minCarrierDensity_m3;
        const bool fluxQualificationEnabled =
            cfg.minFluxScale > 0.0 || cfg.minFluxScaleFraction > 0.0;
        const Real fluxThreshold = std::max(
            cfg.minFluxScale,
            cfg.minFluxScaleFraction * maxFluxScale);
        const bool fluxQualified = fluxQualificationEnabled &&
            fluxAbsSum > 0.0 && fluxAbsSum >= fluxThreshold;
        const bool qualified = sourceQualified || densityQualified || fluxQualified;
        const Real ratio = scale > 0.0 ? std::abs(residual) / scale : 0.0;
        if (qualified) {
            ++evaluation.qualifiedRowCount;
            auto& group = floorLimited ? evaluation.floorLimited : evaluation.nonFloorLimited;
            ++group.rowCount;
            group.violationCount += ratio > eps ? 1 : 0;
            group.maxRatio = std::max(group.maxRatio, ratio);
            group.maxAbsoluteResidual = std::max(group.maxAbsoluteResidual, std::abs(residual));
        } else
            ++evaluation.ignoredRowCount;
        if (qualified && ratio > evaluation.maxRatio) {
            evaluation.maxRatio = ratio;
            evaluation.maxRatioNode = row.nodeId;
            evaluation.maxRatioCarrier = carrier;
        }
        if ((qualified && ratio > eps) || activeRowAudit) {
            NewtonCarrierRowConvergenceViolation violation;
            violation.qualified = qualified;
            violation.nodeId = row.nodeId;
            violation.carrier = carrier;
            violation.residual = residual;
            violation.scale = scale;
            violation.ratio = ratio;
            violation.physicalScale = physicalScale;
            violation.floorLimited = floorLimited;
            violation.blocksAcceptance = !(floorLimited && cfg.floorLimitedAction == "report");
            violation.carrierDensity_m3 = carrierDensity;
            violation.flux = flux;
            violation.recombination = recombination;
            violation.impact = impact;
            violation.densityQualified = densityQualified;
            violation.fluxQualified = fluxQualified;
            violation.sourceQualified = sourceQualified;
            if (activeRowAudit) {
                auto audited = violation;
                audited.blocksAcceptance = qualified && ratio > eps && violation.blocksAcceptance;
                activeRowAudit->push_back(std::move(audited));
            }
            if (qualified && ratio > eps) {
                evaluation.blockingViolationCount += violation.blocksAcceptance ? 1 : 0;
                evaluation.violations.push_back(std::move(violation));
            }
        }
    };

    for (const CoupledDDCarrierTermDiagnostic& row : rows) {
        checkCarrier(row, "electron", row.electronContinuityActive,
                     row.electronDensity_m3,
                     row.electronResidual, row.electronFlux,
                     row.electronFluxAbsSum, row.electronRecombination,
                     row.electronImpact, maxElectronFluxScale,
                     maxElectronSourceScale);
        checkCarrier(row, "hole", row.holeContinuityActive,
                     row.holeDensity_m3,
                     row.holeResidual, row.holeFlux,
                     row.holeFluxAbsSum, row.holeRecombination,
                     row.holeImpact, maxHoleFluxScale,
                     maxHoleSourceScale);
    }
    evaluation.satisfied = evaluation.blockingViolationCount == 0;
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

    long double electronContactFlux = 0.0L;
    long double electronIntegratedSource = 0.0L;
    long double holeContactFlux = 0.0L;
    long double holeIntegratedSource = 0.0L;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        if (electronContact[i])
            electronContactFlux += static_cast<long double>(row.electronFlux);
        else
            electronIntegratedSource +=
                static_cast<long double>(row.electronRecombination) +
                static_cast<long double>(row.electronImpact);
        if (holeContact[i])
            holeContactFlux += static_cast<long double>(row.holeFlux);
        else
            holeIntegratedSource +=
                static_cast<long double>(row.holeRecombination) +
                static_cast<long double>(row.holeImpact);
        if(electronContact[i])evaluation.electron.contactNodeFluxAbsSum+=std::abs(row.electronFlux);
        else evaluation.electron.integratedSourceAbsSum+=std::abs(row.electronRecombination+row.electronImpact);
        if(holeContact[i])evaluation.hole.contactNodeFluxAbsSum+=std::abs(row.holeFlux);
        else evaluation.hole.integratedSourceAbsSum+=std::abs(row.holeRecombination+row.holeImpact);
    }
    evaluation.electron.contactFlux = static_cast<Real>(electronContactFlux);
    evaluation.electron.integratedSource =
        static_cast<Real>(electronIntegratedSource);
    evaluation.hole.contactFlux = static_cast<Real>(holeContactFlux);
    evaluation.hole.integratedSource =
        static_cast<Real>(holeIntegratedSource);

    auto finishCarrier = [&](NewtonGlobalContinuityCarrierClosure& carrier) {
        carrier.mismatch = carrier.contactFlux - carrier.integratedSource;
        const Real sourceMagnitude = std::abs(carrier.integratedSource);
        const Real physicalScale =
            std::max(std::abs(carrier.contactFlux), sourceMagnitude);
        carrier.qualified = physicalScale >= cfg.sourceFloor;
        const Real scale = std::max(physicalScale, cfg.sourceFloor);
        carrier.netScale=scale;
        carrier.ratio = scale > 0.0 ? std::abs(carrier.mismatch) / scale : 0.0;
        carrier.satisfied = carrier.qualified
            ? carrier.ratio <= cfg.tolerance
            : std::abs(carrier.mismatch) <= cfg.absoluteTolerance;
        if (!carrier.satisfied)
            evaluation.satisfied = false;
    };
    finishCarrier(evaluation.electron);
    finishCarrier(evaluation.hole);
    return evaluation;
}

namespace {

std::array<int, 2> selectGlobalContinuityReplacementRows(
    const SparseMatrixd&,
    int nodeCount,
    const std::vector<bool>& electronContact,
    const std::vector<bool>& holeContact)
{
    auto firstFreeRow = [&](int rowOffset,
                            const std::vector<bool>& contact) {
        for (int node = 0; node < nodeCount; ++node) {
            if (!contact[static_cast<std::size_t>(node)])
                return rowOffset + node;
        }
        return -1;
    };

    const int electronRow = firstFreeRow(nodeCount, electronContact);
    const int holeRow = firstFreeRow(2 * nodeCount, holeContact);
    return {electronRow, holeRow};
}

} // namespace

bool stabilizeGlobalContinuityLinearSystem(
    SparseMatrixd& matrix,
    VectorXd& rhs,
    int nodeCount,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes)
{
    if (nodeCount <= 0 || matrix.rows() != 3 * nodeCount ||
        matrix.cols() != 3 * nodeCount || rhs.size() != 3 * nodeCount) {
        throw std::invalid_argument(
            "stabilizeGlobalContinuityLinearSystem: expected a 3N square system.");
    }

    std::vector<bool> electronContact(static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> holeContact(static_cast<std::size_t>(nodeCount), false);
    for (Index node : electronContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            electronContact[static_cast<std::size_t>(node)] = true;
    }
    for (Index node : holeContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            holeContact[static_cast<std::size_t>(node)] = true;
    }

    const std::array<int, 2> replacementRows =
        selectGlobalContinuityReplacementRows(
            matrix, nodeCount, electronContact, holeContact);
    const int electronReplacementRow = replacementRows[0];
    const int holeReplacementRow = replacementRows[1];
    if (electronReplacementRow < 0 && holeReplacementRow < 0)
        return false;

    const int systemSize = 3 * nodeCount;
    std::vector<long double> electronRowSum(
        static_cast<std::size_t>(systemSize), 0.0L);
    std::vector<long double> holeRowSum(
        static_cast<std::size_t>(systemSize), 0.0L);
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry) {
            const int row = entry.row();
            if (row >= nodeCount && row < 2 * nodeCount &&
                !electronContact[static_cast<std::size_t>(row - nodeCount)]) {
                electronRowSum[static_cast<std::size_t>(entry.col())] +=
                    static_cast<long double>(entry.value());
            } else if (row >= 2 * nodeCount && row < 3 * nodeCount &&
                       !holeContact[static_cast<std::size_t>(row - 2 * nodeCount)]) {
                holeRowSum[static_cast<std::size_t>(entry.col())] +=
                    static_cast<long double>(entry.value());
            }
        }
    }

    long double electronRhsSum = 0.0L;
    long double holeRhsSum = 0.0L;
    for (int node = 0; node < nodeCount; ++node) {
        if (!electronContact[static_cast<std::size_t>(node)])
            electronRhsSum += static_cast<long double>(rhs(nodeCount + node));
        if (!holeContact[static_cast<std::size_t>(node)])
            holeRhsSum += static_cast<long double>(rhs(2 * nodeCount + node));
    }

    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(static_cast<std::size_t>(matrix.nonZeros()) +
                     2U * static_cast<std::size_t>(systemSize));
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry) {
            if (entry.row() != electronReplacementRow &&
                entry.row() != holeReplacementRow) {
                triplets.emplace_back(entry.row(), entry.col(), entry.value());
            }
        }
    }
    auto appendReplacementRow = [&](int row,
                                    const std::vector<long double>& values) {
        if (row < 0)
            return;
        for (int col = 0; col < systemSize; ++col) {
            const Real value =
                static_cast<Real>(values[static_cast<std::size_t>(col)]);
            if (value != 0.0)
                triplets.emplace_back(row, col, value);
        }
    };
    appendReplacementRow(electronReplacementRow, electronRowSum);
    appendReplacementRow(holeReplacementRow, holeRowSum);

    SparseMatrixd stabilized(systemSize, systemSize);
    stabilized.setFromTriplets(triplets.begin(), triplets.end());
    stabilized.makeCompressed();
    matrix = std::move(stabilized);
    if (electronReplacementRow >= 0)
        rhs(electronReplacementRow) = static_cast<Real>(electronRhsSum);
    if (holeReplacementRow >= 0)
        rhs(holeReplacementRow) = static_cast<Real>(holeRhsSum);
    return true;
}

VectorXd solveGlobalContinuityStabilizedSystem(
    LinearSolver& linearSolver,
    const SparseMatrixd& matrix,
    const VectorXd& rhs,
    int nodeCount,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes)
{
    SparseMatrixd stabilized = matrix;
    VectorXd stabilizedRhs = rhs;
    if (!stabilizeGlobalContinuityLinearSystem(
            stabilized,
            stabilizedRhs,
            nodeCount,
            electronContactNodes,
            holeContactNodes)) {
        return linearSolver.solve(matrix, rhs);
    }

    std::vector<bool> electronContact(static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> holeContact(static_cast<std::size_t>(nodeCount), false);
    for (Index node : electronContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            electronContact[static_cast<std::size_t>(node)] = true;
    }
    for (Index node : holeContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            holeContact[static_cast<std::size_t>(node)] = true;
    }

    const std::array<int, 2> selectedRows =
        selectGlobalContinuityReplacementRows(
            matrix, nodeCount, electronContact, holeContact);
    std::vector<int> replacementRows;
    for (int row : selectedRows) {
        if (row >= 0)
            replacementRows.push_back(row);
    }

    const int systemSize = 3 * nodeCount;
    const int correctionRank = static_cast<int>(replacementRows.size());
    Eigen::MatrixXd rowUpdates = Eigen::MatrixXd::Zero(
        systemSize, correctionRank);
    Eigen::VectorXd rhsUpdates(correctionRank);
    for (int update = 0; update < correctionRank; ++update) {
        const int row = replacementRows[static_cast<std::size_t>(update)];
        rhsUpdates(update) = stabilizedRhs(row) - rhs(row);
        for (int col = 0; col < systemSize; ++col) {
            rowUpdates(col, update) =
                stabilized.coeff(row, col) - matrix.coeff(row, col);
        }
    }

    const VectorXd baseSolution = linearSolver.solve(matrix, rhs);
    Eigen::MatrixXd inverseReplacementColumns(systemSize, correctionRank);
    for (int update = 0; update < correctionRank; ++update) {
        VectorXd unit = VectorXd::Zero(systemSize);
        unit(replacementRows[static_cast<std::size_t>(update)]) = 1.0;
        inverseReplacementColumns.col(update) =
            linearSolver.solve(matrix, unit);
    }

    const Eigen::MatrixXd border =
        Eigen::MatrixXd::Identity(correctionRank, correctionRank) +
        rowUpdates.transpose() * inverseReplacementColumns;
    Eigen::FullPivLU<Eigen::MatrixXd> borderSolver(border);
    if (!borderSolver.isInvertible()) {
        throw std::runtime_error(
            "solveGlobalContinuityStabilizedSystem: singular Woodbury border.");
    }
    const Eigen::VectorXd correction = borderSolver.solve(
        rhsUpdates - rowUpdates.transpose() * baseSolution);
    return baseSolution + inverseReplacementColumns * correction;
}

bool restoreGlobalContinuityLinearizedStep(
    VectorXd& step,
    const SparseMatrixd& matrix,
    const VectorXd& residual,
    int nodeCount,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes,
    nlohmann::json* diagnostics,
    const std::vector<CoupledDDCarrierTermDiagnostic>* physicalRows)
{
    if (diagnostics)
        *diagnostics = {{"status", "invalid_system"}};
    const int systemSize = 3 * nodeCount;
    if (nodeCount <= 0 || matrix.rows() != systemSize ||
        matrix.cols() != systemSize || residual.size() != systemSize ||
        step.size() != systemSize) {
        return false;
    }

    std::vector<bool> electronContact(static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> holeContact(static_cast<std::size_t>(nodeCount), false);
    for (Index node : electronContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            electronContact[static_cast<std::size_t>(node)] = true;
    }
    for (Index node : holeContactNodes) {
        if (node < static_cast<Index>(nodeCount))
            holeContact[static_cast<std::size_t>(node)] = true;
    }

    if (physicalRows) {
        if (physicalRows->size() != static_cast<std::size_t>(nodeCount)) return false;
        for (int node = 0; node < nodeCount; ++node) {
            electronContact[node] = electronContact[node] || !(*physicalRows)[node].electronContinuityActive;
            holeContact[node] = holeContact[node] || !(*physicalRows)[node].holeContinuityActive;
        }
    }

    std::array<std::vector<long double>, 2> closureRows{
        std::vector<long double>(static_cast<std::size_t>(systemSize), 0.0L),
        std::vector<long double>(static_cast<std::size_t>(systemSize), 0.0L)};
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry) {
            const int row = entry.row();
            const int col = entry.col();
            if (row >= nodeCount && row < 2 * nodeCount &&
                !electronContact[static_cast<std::size_t>(row - nodeCount)]) {
                closureRows[0][static_cast<std::size_t>(col)] +=
                    static_cast<long double>(entry.value());
            } else if (row >= 2 * nodeCount && row < 3 * nodeCount &&
                       !holeContact[static_cast<std::size_t>(row - 2 * nodeCount)]) {
                closureRows[1][static_cast<std::size_t>(col)] +=
                    static_cast<long double>(entry.value());
            }
        }
    }

    // Do not move electrostatic or Dirichlet unknowns. Both free carrier
    // blocks remain active because recombination and avalanche couple the two
    // global continuity rows.
    for (int col = 0; col < systemSize; ++col) {
        const bool freeElectronQf = col >= nodeCount && col < 2 * nodeCount &&
            !electronContact[static_cast<std::size_t>(col - nodeCount)];
        const bool freeHoleQf = col >= 2 * nodeCount && col < 3 * nodeCount &&
            !holeContact[static_cast<std::size_t>(col - 2 * nodeCount)];
        if (!freeElectronQf && !freeHoleQf) {
            closureRows[0][static_cast<std::size_t>(col)] = 0.0L;
            closureRows[1][static_cast<std::size_t>(col)] = 0.0L;
        }
    }

    std::array<long double, 2> rowNorm{0.0L, 0.0L};
    std::array<long double, 2> predictedSum{0.0L, 0.0L};
    for (int carrier = 0; carrier < 2; ++carrier) {
        for (int node = 0; node < nodeCount; ++node) {
            const bool contact = carrier == 0
                ? electronContact[static_cast<std::size_t>(node)]
                : holeContact[static_cast<std::size_t>(node)];
            if (!contact) {
                predictedSum[carrier] += static_cast<long double>(
                    residual((carrier + 1) * nodeCount + node));
            }
        }
        for (int col = 0; col < systemSize; ++col) {
            const long double value =
                closureRows[carrier][static_cast<std::size_t>(col)];
            rowNorm[carrier] = std::hypotl(rowNorm[carrier], value);
            predictedSum[carrier] +=
                value * static_cast<long double>(step(col));
        }
        if (!(rowNorm[carrier] > 0.0L) ||
            !std::isfinite(rowNorm[carrier]) ||
            !std::isfinite(predictedSum[carrier])) {
            return false;
        }
    }

    long double gram00 = 0.0L;
    long double gram01 = 0.0L;
    long double gram11 = 0.0L;
    for (int col = 0; col < systemSize; ++col) {
        const long double a0 =
            closureRows[0][static_cast<std::size_t>(col)] / rowNorm[0];
        const long double a1 =
            closureRows[1][static_cast<std::size_t>(col)] / rowNorm[1];
        gram00 += a0 * a0;
        gram01 += a0 * a1;
        gram11 += a1 * a1;
    }
    const long double determinant = gram00 * gram11 - gram01 * gram01;
    if (diagnostics) {
        // Diagnostic SVD only: preserve the existing projection arithmetic.
        Eigen::MatrixXd normalizedRows(2, systemSize);
        for (int col = 0; col < systemSize; ++col)
            for (int carrier = 0; carrier < 2; ++carrier)
                normalizedRows(carrier, col) = static_cast<Real>(
                    closureRows[carrier][static_cast<std::size_t>(col)] / rowNorm[carrier]);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalizedRows);
        const auto singular = svd.singularValues();
        *diagnostics = {
            {"status", "degenerate_gram"},
            {"row_norms", {static_cast<Real>(rowNorm[0]), static_cast<Real>(rowNorm[1])}},
            {"normalized_gram", {static_cast<Real>(gram00), static_cast<Real>(gram01), static_cast<Real>(gram11)}},
            {"sin_theta", static_cast<Real>(std::sqrt(std::max(0.0L, determinant / (gram00 * gram11))))},
            {"normalized_singular_values", {singular(0), singular(1)}},
            {"normalized_condition_number", singular(1) > 0 ? singular(0) / singular(1) : std::numeric_limits<Real>::infinity()},
            {"predicted_free_sums_before", {static_cast<Real>(predictedSum[0]), static_cast<Real>(predictedSum[1])}}};
    }
    const long double determinantScale =
        std::max(std::abs(gram00 * gram11), 1.0L);
    if (!std::isfinite(determinant) ||
        std::abs(determinant) <=
            64.0L * std::numeric_limits<long double>::epsilon() *
                determinantScale) {
        return false;
    }

    const long double rhs0 = -predictedSum[0] / rowNorm[0];
    const long double rhs1 = -predictedSum[1] / rowNorm[1];
    const long double multiplier0 =
        (rhs0 * gram11 - rhs1 * gram01) / determinant;
    const long double multiplier1 =
        (gram00 * rhs1 - gram01 * rhs0) / determinant;
    if (diagnostics)
        (*diagnostics)["multipliers"] = {static_cast<Real>(multiplier0), static_cast<Real>(multiplier1)};
    VectorXd projected = step;
    for (int col = 0; col < systemSize; ++col) {
        const long double correction =
            closureRows[0][static_cast<std::size_t>(col)] / rowNorm[0] *
                multiplier0 +
            closureRows[1][static_cast<std::size_t>(col)] / rowNorm[1] *
                multiplier1;
        const long double value =
            static_cast<long double>(projected(col)) + correction;
        if (!std::isfinite(value))
            return false;
        projected(col) = static_cast<Real>(value);
    }
    step = std::move(projected);
    if (diagnostics) {
        (*diagnostics)["status"] = "projected";
        std::vector<Index> excludedElectron, excludedHole;
        for (int node = 0; node < nodeCount; ++node) {
            if (electronContact[node]) excludedElectron.push_back(node);
            if (holeContact[node]) excludedHole.push_back(node);
        }
        (*diagnostics)["physical_rows_only"] = physicalRows != nullptr;
        (*diagnostics)["excluded_electron_rows"] = excludedElectron.size();
        (*diagnostics)["excluded_hole_rows"] = excludedHole.size();
        (*diagnostics)["predicted_free_sums_after"] = globalContinuityLinearizedSums(
            matrix, residual, step, nodeCount, excludedElectron, excludedHole);
    }
    return true;
}

std::array<Real, 2> globalContinuityLinearizedSums(
    const SparseMatrixd& matrix, const VectorXd& residual, const VectorXd& step,
    int nodeCount, const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes)
{
    if (nodeCount <= 0 || matrix.rows() != 3 * nodeCount ||
        matrix.cols() != 3 * nodeCount || residual.size() != 3 * nodeCount ||
        step.size() != 3 * nodeCount)
        throw std::invalid_argument("globalContinuityLinearizedSums: expected a 3N system.");
    std::vector<int> carrierForRow(3 * nodeCount, -1);
    for (int node = 0; node < nodeCount; ++node) {
        carrierForRow[nodeCount + node] = 0;
        carrierForRow[2 * nodeCount + node] = 1;
    }
    for (Index node : electronContactNodes)
        if (node >= 0 && node < nodeCount) carrierForRow[nodeCount + node] = -1;
    for (Index node : holeContactNodes)
        if (node >= 0 && node < nodeCount) carrierForRow[2 * nodeCount + node] = -1;
    std::array<long double, 2> sums{};
    for (int row = nodeCount; row < 3 * nodeCount; ++row)
        if (carrierForRow[row] >= 0) sums[carrierForRow[row]] += residual(row);
    // Sum columns before applying the direction to retain cancellation precision.
    std::array<std::vector<long double>, 2> rows{
        std::vector<long double>(3 * nodeCount), std::vector<long double>(3 * nodeCount)};
    for (int outer = 0; outer < matrix.outerSize(); ++outer)
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry)
            if (carrierForRow[entry.row()] >= 0)
                rows[carrierForRow[entry.row()]][entry.col()] += entry.value();
    for (int carrier = 0; carrier < 2; ++carrier)
        for (int col = 0; col < 3 * nodeCount; ++col)
            sums[carrier] += rows[carrier][col] * static_cast<long double>(step(col));
    return {static_cast<Real>(sums[0]), static_cast<Real>(sums[1])};
}

Real globalContinuityIdealAlphaLowerBound(
    const NewtonGlobalContinuityClosureEvaluation& evaluation,
    const NewtonGlobalContinuityClosureConfig& config)
{
    Real lower = 0.0;
    for (const auto* carrier : {&evaluation.electron, &evaluation.hole}) {
        const Real limit = carrier->qualified
            ? config.tolerance * std::max({std::abs(carrier->contactFlux),
                std::abs(carrier->integratedSource), config.sourceFloor})
            : config.absoluteTolerance;
        if (std::abs(carrier->mismatch) > limit)
            lower = std::max(lower, 1.0 - limit / std::abs(carrier->mismatch));
    }
    return lower;
}

namespace {

VectorXd continuityRowWeights(
    const CoupledDDAssembler& assembler,
    const VectorXd& state,
    const CoupledDDBoundaryConditions& bcs,
    const NewtonContinuityRowScalingConfig& cfg,
    const std::vector<CoupledDDCarrierTermDiagnostic>* cachedTerms = nullptr)
{
    const int N = static_cast<int>(assembler.numNodes());
    VectorXd weights = VectorXd::Ones(3 * N);
    if (!cfg.enabled)
        return weights;

    const auto computedTerms = cachedTerms == nullptr
        ? assembler.carrierContinuityEquationTermDiagnostics(state, bcs)
        : std::vector<CoupledDDCarrierTermDiagnostic>{};
    const auto& terms = cachedTerms == nullptr ? computedTerms : *cachedTerms;
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

struct TwoSidedL2EquilibratedSystem {
    SparseMatrixd matrix;
    VectorXd rhs;
    VectorXd rowScale;
    VectorXd columnScale;
};

VectorXd sparseRowL2Norms(const SparseMatrixd& matrix)
{
    VectorXd squared = VectorXd::Zero(matrix.rows());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry)
            squared(entry.row()) += entry.value() * entry.value();
    }
    return squared.cwiseSqrt();
}

VectorXd sparseColumnL2Norms(const SparseMatrixd& matrix)
{
    VectorXd squared = VectorXd::Zero(matrix.cols());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(matrix, outer); entry; ++entry)
            squared(entry.col()) += entry.value() * entry.value();
    }
    return squared.cwiseSqrt();
}

Real positiveNormSpread(const VectorXd& norms)
{
    Real minimum = std::numeric_limits<Real>::infinity();
    Real maximum = 0.0;
    for (int i = 0; i < norms.size(); ++i) {
        const Real value = norms(i);
        if (!std::isfinite(value) || value <= 0.0)
            continue;
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    return std::isfinite(minimum) && minimum > 0.0 ? maximum / minimum : 0.0;
}

TwoSidedL2EquilibratedSystem twoSidedL2Equilibrate(
    const SparseMatrixd& matrix,
    const VectorXd& rhs)
{
    TwoSidedL2EquilibratedSystem result;
    result.rowScale = VectorXd::Ones(matrix.rows());
    const VectorXd rowNorms = sparseRowL2Norms(matrix);
    for (int row = 0; row < rowNorms.size(); ++row) {
        if (std::isfinite(rowNorms(row)) && rowNorms(row) > 0.0)
            result.rowScale(row) = 1.0 / rowNorms(row);
    }

    result.matrix = leftScaleRows(matrix, result.rowScale);
    result.rhs = rhs.cwiseProduct(result.rowScale);
    result.columnScale = VectorXd::Ones(matrix.cols());
    const VectorXd columnNorms = sparseColumnL2Norms(result.matrix);
    for (int column = 0; column < columnNorms.size(); ++column) {
        if (std::isfinite(columnNorms(column)) && columnNorms(column) > 0.0)
            result.columnScale(column) = 1.0 / columnNorms(column);
    }
    for (int outer = 0; outer < result.matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator entry(result.matrix, outer); entry; ++entry)
            entry.valueRef() *= result.columnScale(entry.col());
    }
    result.matrix.makeCompressed();
    return result;
}

VectorXd solveNewtonLinearSystem(const SparseMatrixd& J, const VectorXd& residual,
                                const VectorXd& rowWeights, const NewtonConfig& cfg,
                                LinearSolver& solver)
{
    const auto solve = [&solver, &cfg](const SparseMatrixd& matrix,
                                     const VectorXd& rhs) -> VectorXd {
        if (cfg.linearEquilibration.mode != "l2_row_column")
            return solver.solve(matrix, rhs);
        ScopedPerformanceTimer timer("newton.linear_l2_row_column");
        const auto equilibrated = twoSidedL2Equilibrate(matrix, rhs);
        return solver.solve(equilibrated.matrix, equilibrated.rhs)
            .cwiseProduct(equilibrated.columnScale).eval();
    };
    VectorXd step;
    if (cfg.continuityRowScaling.enabled) {
        ScopedPerformanceTimer timer("newton.linear_row_scaling");
        step = solve(leftScaleRows(J, rowWeights),
                     (-residual.cwiseProduct(rowWeights)).eval());
    } else {
        step = solve(J, -residual);
    }
    if (cfg.linearRefinementIterations > 0) {
        ScopedPerformanceTimer timer("newton.linear_refinement");
        const VectorXd weights = cfg.continuityRowScaling.enabled
            ? rowWeights : VectorXd::Ones(residual.size());
        // Corrections must use the same equilibration as the initial solve.
        step = linear_refinement::refineWithSolve(J, residual, weights,
            std::move(step), solve, cfg.linearRefinementIterations,
            [](int, const VectorXd&, const linear_refinement::Defect&) {});
    }
    return step;
}

DDScalingSpec buildRecoveryScalingSpec(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       const NewtonConfig& cfg)
{
    DDScalingSpec scaling;
    scaling.poissonChargeVolumePolicy = cfg.poissonChargeVolumePolicy;
    scaling.poissonPerCellPermittivity = cfg.poissonPerCellPermittivity;
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

    const double Vt = thermalVoltageForInputs(cfg.temperature_K, cfg.unitScalingRefs);
    MobilityModelConfig mobilityConfig = cfg.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg.recombination, cfg.taun, cfg.taup, cfg.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg.augerCn;
    recombinationConfig.augerCp = cfg.augerCp;
    recombinationConfig.augerDensityDependence = cfg.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg.augerExcessProduct;
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
                    qfAssembler.electronReferencePrefactor()[i], Nc[static_cast<std::size_t>(i)], psiSi,
                    recoveredN, Vt, cfg.carrierStatistics);
                if (result.solution.hasReferencedElectronQuasiFermi()) {
                    result.solution.phinIncrement(i) = result.solution.phin(i)
                        - result.solution.electronQuasiFermiReferenceAt(i);
                }
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
                    qfAssembler.holeReferencePrefactor()[i], Nv[static_cast<std::size_t>(i)], psiSi,
                    recoveredP, Vt, cfg.carrierStatistics);
                if (result.solution.hasReferencedHoleQuasiFermi()) {
                    result.solution.phipIncrement(i) = result.solution.phip(i)
                        - result.solution.holeQuasiFermiReferenceAt(i);
                }
                result.maxCarrierDensityRatio = std::max(
                    result.maxCarrierDensityRatio, recoveredP / before);
                ++result.holeRowsUpdated;
            }
        }
    }
    result.maxPsiDelta_V = 0.0;
    return result;
}
namespace {
void resolveFermiSgDefault(NewtonConfig& config)
{
    auto& statistics = config.carrierStatistics;
    if (statistics.sgFormulation != FermiSgFormulation::Automatic)
        return;
    const bool classicalFermi = statistics.model == "fermi_dirac" &&
        config.coupledOperator == "legacy" && config.impactIonization.model == "none" &&
        !config.electronQuantumPotential.enabled;
    statistics.sgFormulation = classicalFermi
        ? FermiSgFormulation::EffectivePotential
        : FermiSgFormulation::GeneralizedEinstein;
}

// One coordinate policy for JSON, direct C++ construction and CMD execution.
// Positive populations and equation gates are not QF trust-region controls.
void resolvePopulationDefaults(NewtonConfig& c) {
    if(c.stateVariables.empty())c.stateVariables=c.coupledOperator=="bvmethods_tri3_qf"?"potential_density":"potential_quasi_fermi";
    const bool population=c.stateVariables=="potential_density";
    const bool bvPopulation=population && c.coupledOperator=="bvmethods_tri3_qf";
    if(c.populationMinimumDensity_cm3==-1.)c.populationMinimumDensity_cm3=bvPopulation?1e-100:0.;
    if(c.populationDerivativeRecovery.empty())c.populationDerivativeRecovery=bvPopulation?"after_protection":"full";
    if(c.primaryVariableCoordinates.empty())c.primaryVariableCoordinates=bvPopulation?"unit":"physical";
    if(c.populationStoppingRule.empty())c.populationStoppingRule=bvPopulation?"native_approx":"legacy";
    if(c.primaryVariableCoordinates!="physical"&&c.primaryVariableCoordinates!="unit")
        throw std::invalid_argument("primary_variable_coordinates must be physical or unit");
    if(c.populationStoppingRule!="legacy"&&c.populationStoppingRule!="native_approx")
        throw std::invalid_argument("population_stopping_rule must be legacy or native_approx");
    if(!population&&(c.primaryVariableCoordinates!="physical"||c.populationStoppingRule!="legacy"))
        throw std::invalid_argument("Unit primary coordinates and native_approx stopping require potential_density");
}

void validatePopulationNumerics(const NewtonConfig& c)
{
    if (c.coupledOperator != "bvmethods_tri3_qf" && c.bvMethodsNumerics != BvMethodsNumericalOptions{})
        throw std::invalid_argument("bvmethods_numerics requires bvmethods_tri3_qf");
    if (c.coupledOperator == "bvmethods_tri3_qf" && c.populationSgCurrent == "band_population" &&
        (c.carrierStatistics.model != "fermi_dirac_joyce_dixon" ||
         c.bvMethodsNumerics.avalancheCurrent != BvMethodsAvalancheCurrent::RepresentativeVertex))
        throw std::invalid_argument("band_population requires Joyce-Dixon and representative-vertex avalanche current");
    if (c.populationLinearDefect != "mp100" && c.populationLinearDefect != "eft")
        throw std::invalid_argument("population_linear_defect must be mp100 or eft");
    if (c.populationLinearDefect == "eft" &&
        (c.stateVariables != "potential_density" || c.coupledOperator != "legacy" ||
         c.primaryVariableCoordinates != "unit"))
        throw std::invalid_argument("population_linear_defect=eft requires legacy Unit population coordinates");
    if(c.populationWriteback!="independent" && c.populationWriteback!="boltzmann_psi_compensated")
        throw std::invalid_argument("unknown population_writeback");
    if(c.populationWriteback!="independent" &&
       (c.stateVariables!="potential_density" || c.coupledOperator!="legacy" ||
        c.carrierStatistics.model!="boltzmann"))
        throw std::invalid_argument("population_writeback compensation requires legacy Boltzmann potential_density");
    if (c.compensatedPopulationState &&
        (c.stateVariables != "potential_density" || c.coupledOperator != "legacy" ||
         (c.carrierStatistics.model != "boltzmann" &&
          c.electronQuantumPotential.couplingMode != "monolithic_box" &&
          !(c.carrierStatistics.model == "fermi_dirac" &&
            c.carrierStatistics.sgFormulation == FermiSgFormulation::EffectivePotential &&
            !c.electronQuantumPotential.enabled && c.impactIonization.model == "none")) ||
         !c.inputScaling.isUnitScaling() || c.primaryVariableCoordinates != "unit" ||
         c.contactBoundaryMinorityElectronRelaxation))
        throw std::invalid_argument(
            "Compensated population requires legacy Unit density coordinates and no minority-contact relaxation; classical transport requires Boltzmann statistics or effective-potential integral Fermi-Dirac without avalanche");
    if (c.electronQuantumPotential.couplingMode == "monolithic_box") {
        const auto& dg = c.electronQuantumPotential;
        if (c.coupledOperator != "legacy" || c.splitDDState ||
            (c.stateVariables == "potential_density" &&
             (c.primaryVariableCoordinates != "unit" || !c.inputScaling.isUnitScaling())) ||
            !dg.enabled || !dg.includeInsulators || dg.globalDiscretization != "sentaurus_box" ||
            dg.formulation != "potential_based" || dg.transportCoupling != "direct_band_edge" ||
            c.blockAbsoluteConvergence.mode == "enforce" ||
            c.carrierRowConvergence.mode == "enforce" || c.globalContinuityClosure.mode == "enforce" ||
            c.impactIonization.model != "none" || c.bandToBand.model != "none")
            throw std::invalid_argument("monolithic_box requires legacy DD with the global fitted DG box, no avalanche/BTBT or enforced DD block diagnostics, and Unit scaling for population coordinates");
    }
    if (c.stateVariables != "potential_density")
        return;
    for(double v:{c.populationMinimumDensity_cm3,c.populationPredictorLowDensity_cm3})
        if(!std::isfinite(v)||v<0.)throw std::invalid_argument("population density limits must be finite and nonnegative");
    if((c.populationMinimumDensity_cm3!=0.||c.populationPredictorLowDensity_cm3!=0.||c.populationUpdate.minimumDensity!=0.)&&c.stateVariables!="potential_density")
        throw std::invalid_argument("population density limits require potential_density");
    validatePopulationUpdate(c.populationUpdate);
    const auto reject=[](bool active,const char* option) {
        if(active)throw std::invalid_argument(std::string("potential_density does not support ")+option);
    };
    reject(c.splitDDState ||
        (c.electronQuantumPotential.enabled && c.electronQuantumPotential.couplingMode != "monolithic_box"),
        "split or quantum state formulations");
    reject(c.coupledOperator=="legacy" && (c.populationSgCurrent!="quasi_fermi" ||
        c.populationDerivativeRecovery!="full"), "BVMethods-only SG or derivative recovery");
    reject(c.poissonResidualPrecision!="double" || c.bankRose.enabled || c.populationStepFraction!=0.,
        "QF Poisson arithmetic, Bank-Rose or QF population-step control");
    reject(c.globalContinuityClosure.linearSystemStabilization!="off" ||
        c.globalContinuityClosure.feasibilityRestoration!="off", "QF closure replacement/restoration");
    reject(c.quasiFermiUpdateLimit_V!=0. || c.quasiFermiUpdateLimitMinority_V!=0.,"QF voltage update caps");
    reject(c.maxUpdate!=0.,"max_update on mixed potential/concentration units");
    reject(c.quasiFermiReference!="none" || c.quasiFermiRecenterOnInitialState || c.quasiFermiRecenterOnStall,
        "QF reference/recentering");
    reject(c.quasiFermiUpdateLimitMode!="componentwise_poisson_recorrect" ||
        c.quasiFermiTrustRegionGrowthFactor!=1. || c.quasiFermiTrustRegionMaxMultiplier!=1. ||
        c.quasiFermiTrustRegionExpansionThreshold!=.75 || c.quasiFermiTrustRegionShrinkFactor!=1. ||
        c.quasiFermiTrustRegionMinMultiplier!=1.,"QF trust-region controls");
    reject(c.contactMajorityQfBranchDropLimit_V!=0.,"contact-majority QF branch stopping");
    reject(c.carrierRowRecovery.mode!="off" || c.continuityRowScaling.enabled ||
        c.carrierRegularizationScale!=0. || c.carrierDiagonalFloor.enabled,
        "legacy carrier recovery/scaling/regularization");
    reject(c.linearEquilibration.mode!="off", "legacy linear_equilibration (population columns have their own scaling)");
    reject(c.jacobian!="analytic", "legacy finite-difference Jacobian selection");
    reject(c.stableMeritComparison || c.carrierRowQualifiedStallAcceptance,
        "legacy merit/stall acceptance");
    reject(c.lineSearch && c.lineSearchMode!="merit", "legacy block_filter line search");
}

void validateNumericalPrecisionOptions(const NewtonConfig& cfg)
{
    validateBankRoseConfig(cfg.bankRose);
    if (cfg.bankRose.enabled && (cfg.lineSearch || cfg.lineSearchMode != "merit" ||
        cfg.dampingFactor != 1.0 || cfg.stableMeritComparison))
        throw std::invalid_argument("bank_rose requires line_search=false, line_search_mode=merit, "
            "damping_factor=1 and stable_merit_comparison=false; its decrease test must not be replaced.");
    if (cfg.poissonResidualPrecision != "double" && cfg.poissonResidualPrecision != "binary128")
        throw std::invalid_argument("poisson_residual_precision must be double or binary128.");
    if (cfg.poissonResidualPrecision == "binary128" &&
        (cfg.carrierStatistics.model != "boltzmann" || cfg.electronQuantumPotential.enabled))
        throw std::invalid_argument("binary128 Poisson residual requires classical Boltzmann statistics without a quantum potential.");
    if (cfg.stableMeritComparison &&
        (cfg.lineSearchMode != "merit" || cfg.globalContinuityClosure.mode != "off"))
        throw std::invalid_argument("stable_merit_comparison requires merit line search and global continuity merit off.");
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
    if (json.contains("electron_quantum_potential")) {
        const auto& value = json.at("electron_quantum_potential");
        if (value.is_boolean()) {
            cfg.electronQuantumPotential.enabled = value.get<bool>();
        } else if (value.is_object()) {
            cfg.electronQuantumPotential.enabled = value.value(
                "enabled", cfg.electronQuantumPotential.enabled);
            cfg.electronQuantumPotential.couplingMode = value.value(
                "coupling_mode", cfg.electronQuantumPotential.couplingMode);
            cfg.electronQuantumPotential.transportCoupling = value.value(
                "transport_coupling",
                cfg.electronQuantumPotential.transportCoupling);
            cfg.electronQuantumPotential.transportCouplingWeight = value.value(
                "transport_coupling_weight",
                cfg.electronQuantumPotential.transportCouplingWeight);
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
        else if (statistics.is_object()) {
            cfg.carrierStatistics.model = statistics.value("model", cfg.carrierStatistics.model);
            const auto formulation = statistics.value("sg_formulation", "auto");
            if (formulation == "effective_potential")
                cfg.carrierStatistics.sgFormulation = FermiSgFormulation::EffectivePotential;
            else if (formulation == "generalized_einstein")
                cfg.carrierStatistics.sgFormulation = FermiSgFormulation::GeneralizedEinstein;
            else if (formulation != "auto")
                throw std::invalid_argument("Unknown carrier_statistics.sg_formulation: " + formulation);
            cfg.carrierStatistics.transport = statistics.value("transport", cfg.carrierStatistics.transport);
        }
        else
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_statistics must be a string or object.");
    }
    if (json.contains("statistics"))
        cfg.carrierStatistics.model = json.at("statistics").get<std::string>();
    cfg.dampingFactor = json.value("damping_factor", cfg.dampingFactor);
    if(json.contains("population_step_control")) {
        const auto& p=json.at("population_step_control");
        const auto mode=p.value("mode",std::string("legacy"));
        if(mode!="legacy"&&mode!="fraction_to_boundary")
            throw std::invalid_argument("population_step_control mode must be legacy or fraction_to_boundary");
        cfg.populationStepFraction=mode=="legacy"?0.:p.value("fraction",.99);
        if(mode!="legacy"&&(!(cfg.populationStepFraction>0.&&cfg.populationStepFraction<1.)||!std::isfinite(cfg.populationStepFraction)))
            throw std::invalid_argument("population step fraction must be finite and strictly between zero and one");
    }
    cfg.lineSearch = json.value("line_search", cfg.lineSearch);
    cfg.lineSearchMode = json.value("line_search_mode", cfg.lineSearchMode);
    if (json.contains("bank_rose")) {
        const auto& br = json.at("bank_rose");
        if (!br.is_object()) throw std::invalid_argument("bank_rose must be an object.");
        const std::vector<std::string> keys{"enabled", "not_damped", "omega", "initial_K",
            "seed_K", "growth", "decay", "min_damping", "max_trials"};
        for (auto it = br.begin(); it != br.end(); ++it)
            if (std::find(keys.begin(), keys.end(), it.key()) == keys.end())
                throw std::invalid_argument("Unknown bank_rose key: " + it.key());
        for (const auto* key : {"not_damped", "max_trials"})
            if (br.contains(key) && !br.at(key).is_number_integer())
                throw std::invalid_argument(std::string("bank_rose.") + key + " must be an integer.");
        cfg.bankRose.enabled = br.value("enabled", cfg.bankRose.enabled);
        cfg.bankRose.notDamped = br.value("not_damped", cfg.bankRose.notDamped);
        cfg.bankRose.sufficientDecrease = br.value("omega", cfg.bankRose.sufficientDecrease);
        cfg.bankRose.initialK = br.value("initial_K", cfg.bankRose.initialK);
        cfg.bankRose.seedK = br.value("seed_K", cfg.bankRose.seedK);
        cfg.bankRose.growth = br.value("growth", cfg.bankRose.growth);
        cfg.bankRose.decay = br.value("decay", cfg.bankRose.decay);
        cfg.bankRose.minDamping = br.value("min_damping", cfg.bankRose.minDamping);
        cfg.bankRose.maxTrials = br.value("max_trials", cfg.bankRose.maxTrials);
    }
    cfg.returnBestFeasibleOnFailure = json.value(
        "return_best_feasible_on_failure", cfg.returnBestFeasibleOnFailure);
    cfg.residualFilterGamma = json.value(
        "residual_filter_gamma", cfg.residualFilterGamma);
    cfg.residualFilterEnvelopeFactor = json.value(
        "residual_filter_envelope_factor", cfg.residualFilterEnvelopeFactor);
    cfg.verbose = json.value("verbose", cfg.verbose);
    cfg.warmStart = json.value("warm_start", cfg.warmStart);
    cfg.contactPreprojection = json.value("contact_preprojection", cfg.contactPreprojection);
    cfg.diagnostics = json.value("diagnostics", cfg.diagnostics);
    cfg.diagnostics = json.value("diagnostic_history", cfg.diagnostics);
    cfg.progressDirectory = json.value("progress_directory", std::string{});
    cfg.captureLinearSystem = json.value("capture_linear_system", false);
    cfg.softTimeLimitSeconds = json.value("soft_time_limit_seconds", 0.0);
    cfg.diagnosticStopPolicy = json.value("diagnostic_stop_policy", std::string("off"));
    if (cfg.diagnosticStopPolicy != "off" && cfg.diagnosticStopPolicy != "row_trend_v1")
        throw std::invalid_argument("diagnostic_stop_policy must be off or row_trend_v1");
    if (!std::isfinite(cfg.softTimeLimitSeconds) || cfg.softTimeLimitSeconds < 0)
        throw std::invalid_argument("soft_time_limit_seconds must be finite and nonnegative");
    cfg.maxUpdate = json.value("max_update", cfg.maxUpdate);
    cfg.quasiFermiUpdateLimit_V = json.value(
        "quasi_fermi_update_limit_V",
        cfg.quasiFermiUpdateLimit_V);
    cfg.quasiFermiUpdateLimitMinority_V = json.value(
        "quasi_fermi_update_limit_minority_V",
        cfg.quasiFermiUpdateLimitMinority_V);
    cfg.quasiFermiUpdateLimitMode = json.value(
        "quasi_fermi_update_limit_mode",
        cfg.quasiFermiUpdateLimitMode);
    cfg.quasiFermiTrustRegionGrowthFactor = json.value(
        "quasi_fermi_trust_region_growth_factor",
        cfg.quasiFermiTrustRegionGrowthFactor);
    cfg.quasiFermiTrustRegionMaxMultiplier = json.value(
        "quasi_fermi_trust_region_max_multiplier",
        cfg.quasiFermiTrustRegionMaxMultiplier);
    cfg.quasiFermiTrustRegionExpansionThreshold = json.value(
        "quasi_fermi_trust_region_expansion_threshold",
        cfg.quasiFermiTrustRegionExpansionThreshold);
    cfg.quasiFermiTrustRegionShrinkFactor = json.value(
        "quasi_fermi_trust_region_shrink_factor",
        cfg.quasiFermiTrustRegionShrinkFactor);
    cfg.quasiFermiTrustRegionMinMultiplier = json.value(
        "quasi_fermi_trust_region_min_multiplier",
        cfg.quasiFermiTrustRegionMinMultiplier);
    cfg.linearRefinementIterations = json.value("linear_refinement_iterations", 0);
    cfg.populationLinearDefect = json.value("population_linear_defect", cfg.populationLinearDefect);
    cfg.reuseFixedVoltageBorderedAnalysis =
        json.value("reuse_fixed_voltage_bordered_analysis", cfg.reuseFixedVoltageBorderedAnalysis);
    if (json.contains("legacy_edge_jacobian_steps")) {
        const auto& steps = json.at("legacy_edge_jacobian_steps");
        if (!steps.is_object()) throw std::invalid_argument("legacy_edge_jacobian_steps must be an object");
        cfg.legacyTransportJacobianStep = steps.value("transport", cfg.legacyTransportJacobianStep);
        cfg.legacyAvalancheJacobianStep = steps.value("avalanche", cfg.legacyAvalancheJacobianStep);
    }
    if (!std::isfinite(cfg.legacyTransportJacobianStep) || cfg.legacyTransportJacobianStep <= 0.0 ||
        !std::isfinite(cfg.legacyAvalancheJacobianStep) || cfg.legacyAvalancheJacobianStep <= 0.0)
        throw std::invalid_argument("legacy_edge_jacobian_steps must be finite and positive");
    cfg.poissonResidualPrecision = json.value("poisson_residual_precision", cfg.poissonResidualPrecision);
    cfg.stableMeritComparison = json.value("stable_merit_comparison", cfg.stableMeritComparison);
    cfg.splitDDState = json.value("split_dd_state",cfg.splitDDState);
    cfg.coupledOperator = json.value("coupled_operator",cfg.coupledOperator);
    if (json.contains("bvmethods_numerics")) {
        const auto& options = json.at("bvmethods_numerics");
        if (cfg.coupledOperator != "bvmethods_tri3_qf" || !options.is_object())
            throw std::invalid_argument("bvmethods_numerics requires a BVMethods options object");
        for (const auto& [name, value] : options.items()) {
            if (name == "avalanche_weights" && value == "coefficient_box")
                cfg.bvMethodsNumerics.avalancheWeights = BvMethodsAvalancheWeights::CoefficientBox;
            else if (name == "avalanche_current" && value == "gss_laux")
                cfg.bvMethodsNumerics.avalancheCurrent = BvMethodsAvalancheCurrent::GssLaux;
            else if (name == "avalanche_switch" && value == "sharp")
                cfg.bvMethodsNumerics.avalancheSwitch = BvMethodsAvalancheSwitch::Sharp;
            else if (!(name == "avalanche_weights" && value == "truncated_box") &&
                     !(name == "avalanche_current" && value == "representative_vertex") &&
                     !(name == "avalanche_switch" && value == "smooth"))
                throw std::invalid_argument("Invalid BVMethods numerical option: " + name);
        }
    }
    cfg.stateVariables = json.value("state_variables",cfg.stateVariables);
    cfg.primaryVariableCoordinates=json.value("primary_variable_coordinates",cfg.primaryVariableCoordinates);
    cfg.populationStoppingRule=json.value("population_stopping_rule",cfg.populationStoppingRule);
    cfg.compensatedPopulationState = json.value("compensated_population_state", false);
    resolvePopulationDefaults(cfg);
    if (json.contains("population_update_convergence")) {
        const auto& metric = json.at("population_update_convergence");
        if (cfg.stateVariables != "potential_density" || cfg.populationStoppingRule != "native_approx")
            throw std::invalid_argument("population_update_convergence requires potential_density and native_approx");
        if (!metric.is_object())
            throw std::invalid_argument("population_update_convergence must be an object");
        for (auto field = metric.begin(); field != metric.end(); ++field) {
            if (field.key() != "digits" && field.key() != "potential_reference_V" &&
                field.key() != "electron_reference_cm3" && field.key() != "hole_reference_cm3")
                throw std::invalid_argument("unknown population_update_convergence field: " + field.key());
        }
        if (metric.contains("digits") && !metric.at("digits").is_number_integer())
            throw std::invalid_argument("population_update_convergence digits must be an integer");
        auto& update = cfg.populationUpdateConvergence;
        update.digits = metric.value("digits", update.digits);
        update.potentialReference_V = metric.value("potential_reference_V", update.potentialReference_V);
        update.electronReference_cm3 = metric.value("electron_reference_cm3", update.electronReference_cm3);
        update.holeReference_cm3 = metric.value("hole_reference_cm3", update.holeReference_cm3);
        validateDensityUpdateConvergence(update);
    }
    cfg.populationSgCurrent=json.value("population_sg_current",cfg.populationSgCurrent);
    cfg.populationSgBootstrapLinearization=json.value("population_sg_bootstrap_linearization",cfg.populationSgBootstrapLinearization);
    cfg.populationWriteback=json.value("population_writeback",cfg.populationWriteback);
    if(cfg.stateVariables=="potential_density" && cfg.coupledOperator=="legacy" &&
       (!json.contains("coupled_operator") || !json.contains("carrier_statistics")))
        throw std::invalid_argument("legacy potential_density requires explicit operator and carrier statistics");
    if((cfg.populationSgCurrent!="quasi_fermi" && cfg.populationSgCurrent!="band_population") ||
       (cfg.populationSgCurrent!="quasi_fermi" && cfg.stateVariables!="potential_density"))
        throw std::invalid_argument("population_sg_current requires band_population/potential_density or quasi_fermi");
    if(json.contains("population_update")) {
        const auto& p=json.at("population_update");
        if(!p.is_object())throw std::invalid_argument("population_update must be an object");
        for(auto it=p.begin();it!=p.end();++it)
            if(it.key()!="protection" && it.key()!="local_retained_fraction" && it.key()!="fraction_tau" && it.key()!="minimum_density_cm3" && it.key()!="minimum_density_m3")
                throw std::invalid_argument("unknown population_update field: "+it.key());
        cfg.populationUpdate.protection=p.value("protection",cfg.populationUpdate.protection);
        cfg.populationUpdate.localRetainedFraction=p.value("local_retained_fraction",cfg.populationUpdate.localRetainedFraction);
        cfg.populationUpdate.fractionTau=p.value("fraction_tau",cfg.populationUpdate.fractionTau);
        if(p.contains("minimum_density_cm3")&&!p["minimum_density_cm3"].is_number())throw std::invalid_argument("minimum_density_cm3 must be numeric");
        if(p.contains("minimum_density_cm3")&&p.contains("minimum_density_m3"))throw std::invalid_argument("duplicate population minimum density units");
        if(p.contains("minimum_density_m3")&&!p["minimum_density_m3"].is_number())throw std::invalid_argument("minimum_density_m3 must be numeric");
        cfg.populationMinimumDensity_cm3=p.contains("minimum_density_m3")?scaling.concentrationToSI(p["minimum_density_m3"].get<double>())/1e6:p.value("minimum_density_cm3",cfg.populationMinimumDensity_cm3);
    }
    if(json.contains("population_predictor_low_density_limit_cm3")&&!json["population_predictor_low_density_limit_cm3"].is_number())throw std::invalid_argument("population predictor density limit must be numeric");
    if(json.contains("population_predictor_low_density_limit_cm3")&&json.contains("population_predictor_low_density_limit_m3"))throw std::invalid_argument("duplicate predictor density limit units");
    if(json.contains("population_predictor_low_density_limit_m3")&&!json["population_predictor_low_density_limit_m3"].is_number())throw std::invalid_argument("population predictor density limit must be numeric");
    cfg.populationPredictorLowDensity_cm3=json.contains("population_predictor_low_density_limit_m3")?scaling.concentrationToSI(json["population_predictor_low_density_limit_m3"].get<double>())/1e6:json.value("population_predictor_low_density_limit_cm3",0.);
    for(double v:{cfg.populationMinimumDensity_cm3,cfg.populationPredictorLowDensity_cm3})
        if(!std::isfinite(v)||v<0.)throw std::invalid_argument("population density limits must be finite and nonnegative");
    if(cfg.populationPredictorLowDensity_cm3!=0. && cfg.stateVariables!="potential_density")
        throw std::invalid_argument("population predictor density limit requires potential_density");
    if(cfg.populationSgBootstrapLinearization!="unchanged" &&
       (cfg.populationSgBootstrapLinearization!="floor_adjacent_edges" ||
        cfg.stateVariables!="potential_density" || cfg.coupledOperator!="bvmethods_tri3_qf" ||
        cfg.populationSgCurrent!="quasi_fermi" || !(cfg.populationMinimumDensity_cm3>0.)))
        throw std::invalid_argument("population_sg_bootstrap_linearization floor_adjacent_edges requires BVMethods potential_density, quasi_fermi SG and positive minimum_density_cm3");
    validatePopulationUpdate(cfg.populationUpdate);
    cfg.populationDerivativeRecovery=json.value("population_derivative_recovery",cfg.populationDerivativeRecovery);
    if(cfg.populationDerivativeRecovery!="full" && cfg.populationDerivativeRecovery!="automatic" && cfg.populationDerivativeRecovery!="after_protection")
        throw std::invalid_argument("population_derivative_recovery must be full, automatic or after_protection");
    if(cfg.populationDerivativeRecovery!="full" && cfg.stateVariables!="potential_density")
        throw std::invalid_argument("automatic population derivative recovery requires potential_density");
    if(cfg.stateVariables!="potential_quasi_fermi" && cfg.stateVariables!="potential_density")
        throw std::invalid_argument("unknown state_variables");
    if(cfg.stateVariables=="potential_density") {
        const auto method=json.value("method",std::string{"newton"});
        const bool frozenLegacy = cfg.coupledOperator == "legacy" &&
            (method == "frozen_state" || method == "diagnostic_state_replay");
        if (method != "newton" && method != "poisson_only" && !frozenLegacy)
            throw std::invalid_argument("potential_density requires Newton, Poisson initialization, or legacy frozen replay");
    }
    if(json.contains("population_update") && cfg.stateVariables!="potential_density")
        throw std::invalid_argument("population_update requires potential_density state_variables");
    if(cfg.coupledOperator!="legacy" && cfg.coupledOperator!="bvmethods_tri3_qf")
        throw std::invalid_argument("Unknown coupled_operator; expected legacy or bvmethods_tri3_qf");
    cfg.linearFactorRefinementSteps = json.value("linear_factor_refinement_steps", 0);
    if (cfg.linearFactorRefinementSteps < 0 || cfg.linearFactorRefinementSteps > 3 ||
        (cfg.linearFactorRefinementSteps && cfg.linearRefinementIterations))
        throw std::invalid_argument("linear_factor_refinement_steps requires 0..3 and a separate MP100 refinement arm.");
    if(json.contains("split_contact_bias_low_V")) {
        cfg.splitContactBiasLow_V=json.at("split_contact_bias_low_V").get<std::unordered_map<std::string,Real>>();
        if(!cfg.splitDDState && !cfg.splitContactBiasLow_V.empty())throw std::invalid_argument("Contact low parts require split_dd_state");
        for(const auto& [name,low]:cfg.splitContactBiasLow_V)
            if(name.empty() || !std::isfinite(low))throw std::invalid_argument("Invalid contact low part");
    }
    cfg.exactDirichletUpdates = json.value("exact_dirichlet_updates", cfg.exactDirichletUpdates);
    if (cfg.linearRefinementIterations < 0 || cfg.linearRefinementIterations > 10)
        throw std::invalid_argument("linear_refinement_iterations must be between 0 and 10.");
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
    if (json.contains("block_absolute_convergence")) {
        const auto& value = json.at("block_absolute_convergence");
        if (!value.is_object()) {
            throw std::invalid_argument(
                "newtonConfigFromJson: block_absolute_convergence must be an object.");
        }
        cfg.blockAbsoluteConvergence.mode = value.value(
            "mode", cfg.blockAbsoluteConvergence.mode);
        cfg.blockAbsoluteConvergence.psiResidualCeiling = value.value(
            "psi_residual_ceiling", cfg.blockAbsoluteConvergence.psiResidualCeiling);
        cfg.blockAbsoluteConvergence.electronResidualCeiling = value.value(
            "electron_residual_ceiling",
            cfg.blockAbsoluteConvergence.electronResidualCeiling);
        cfg.blockAbsoluteConvergence.holeResidualCeiling = value.value(
            "hole_residual_ceiling", cfg.blockAbsoluteConvergence.holeResidualCeiling);
    }
    cfg.contactMajorityQfBranchDropLimit_V = json.value(
        "contact_majority_qf_branch_drop_limit_V",
        cfg.contactMajorityQfBranchDropLimit_V);
    cfg.contactMajorityQfBranchGuardContacts = json.value(
        "contact_majority_qf_branch_guard_contacts",
        cfg.contactMajorityQfBranchGuardContacts);
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
    if (json.contains("region_resolved_interface_assembly")) {
        const auto& value = json.at("region_resolved_interface_assembly");
        if (value.is_boolean()) {
            const bool enabled = value.get<bool>();
            cfg.regionResolvedInterfaceAssembly.poissonEdgeCoupling = enabled;
            cfg.regionResolvedInterfaceAssembly.transportEdgeCoupling = enabled;
            cfg.regionResolvedInterfaceAssembly.transportNodeVolume = enabled;
            cfg.regionResolvedInterfaceAssembly
                .transportSignedAverageBoxNodeVolume = false;
            cfg.regionResolvedInterfaceAssembly
                .poissonElectronTransportNodeVolume = false;
            cfg.regionResolvedInterfaceAssembly
                .poissonHoleTransportNodeVolume = false;
            cfg.regionResolvedInterfaceAssembly
                .poissonDopantTransportNodeVolume = false;
        } else if (value.is_object()) {
            auto& assembly = cfg.regionResolvedInterfaceAssembly;
            const bool enabled = value.value("enabled", false);
            assembly.srhSignedTransportVolumeFraction = value.value("srh_signed_transport_volume_fraction", 0.0);
            if (!std::isfinite(assembly.srhSignedTransportVolumeFraction) ||
                std::abs(assembly.srhSignedTransportVolumeFraction) > 1.0)
                throw std::invalid_argument("srh_signed_transport_volume_fraction must be finite and in [-1,1].");
            assembly.poissonChargeNodeVolume = value.value("poisson_charge_node_volume", "inherit");
            assembly.transportEdgeGeometry = value.value("transport_edge_geometry", "inherit");
            if (assembly.transportEdgeGeometry != "inherit" && assembly.transportEdgeGeometry != "element_box")
                throw std::invalid_argument("transport_edge_geometry must be inherit or element_box.");
            assembly.transportEdgeCouplingRatios = value.value(
                "transport_edge_coupling_ratios", std::vector<Real>{});
            assembly.poissonEdgeCoupling = value.value(
                "poisson_edge_coupling", enabled);
            assembly.transportEdgeCoupling = value.value(
                "transport_edge_coupling", enabled);
            assembly.transportNodeVolume = value.value(
                "transport_node_volume", enabled);
            assembly.transportSignedAverageBoxNodeVolume = value.value(
                "transport_signed_average_box_node_volume", false);
            assembly.transportSignedAverageBoxNodeVolumeScope = value.value(
                "transport_signed_average_box_node_volume_scope",
                assembly.transportSignedAverageBoxNodeVolumeScope);
            assembly.poissonElectronTransportNodeVolume = value.value(
                "poisson_electron_transport_node_volume", false);
            assembly.poissonHoleTransportNodeVolume = value.value(
                "poisson_hole_transport_node_volume", false);
            assembly.poissonDopantTransportNodeVolume = value.value(
                "poisson_dopant_transport_node_volume", false);
            if (assembly.poissonChargeNodeVolume != "inherit" &&
                assembly.poissonChargeNodeVolume != "signed_transport")
                throw std::invalid_argument("poisson_charge_node_volume must be inherit or signed_transport.");
            if (assembly.poissonChargeNodeVolume != "inherit" &&
                (assembly.poissonElectronTransportNodeVolume || assembly.poissonHoleTransportNodeVolume ||
                 assembly.poissonDopantTransportNodeVolume))
                throw std::invalid_argument("Charge volume policy cannot combine with per-term volume flags.");
            if (!assembly.transportEdgeCouplingRatios.empty() && assembly.transportEdgeCoupling)
                throw std::invalid_argument("Explicit transport ratios cannot combine with transport_edge_coupling.");
            for (Real ratio : assembly.transportEdgeCouplingRatios)
                if (!std::isfinite(ratio) || ratio < 0)
                    throw std::invalid_argument("Transport ratios must be finite and nonnegative.");
            if (assembly.transportSignedAverageBoxNodeVolumeScope !=
                    "all_transport_nodes" &&
                assembly.transportSignedAverageBoxNodeVolumeScope !=
                    "external_boundary_contact_support") {
                throw std::invalid_argument(
                    "newtonConfigFromJson: "
                    "transport_signed_average_box_node_volume_scope must be "
                    "'all_transport_nodes' or "
                    "'external_boundary_contact_support'.");
            }
            if (assembly.transportNodeVolume &&
                assembly.transportSignedAverageBoxNodeVolume) {
                throw std::invalid_argument(
                    "newtonConfigFromJson: transport_node_volume and "
                    "transport_signed_average_box_node_volume are mutually "
                    "exclusive.");
            }
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: region_resolved_interface_assembly "
                "must be a boolean or object.");
        }
    }
    if (json.contains("local_update_diagnostics")) {
        const auto& value = json.at("local_update_diagnostics");
        if (value.is_boolean()) {
            cfg.localUpdateDiagnostics.enabled = value.get<bool>();
        } else if (value.is_object()) {
            cfg.localUpdateDiagnostics.enabled = value.value(
                "enabled", cfg.localUpdateDiagnostics.enabled);
            cfg.localUpdateDiagnostics.csvFile = value.value(
                "csv_file", cfg.localUpdateDiagnostics.csvFile);
            cfg.localUpdateDiagnostics.nodes = value.value(
                "nodes", cfg.localUpdateDiagnostics.nodes);
            cfg.localUpdateDiagnostics.firstIterations = value.value(
                "first_iterations", cfg.localUpdateDiagnostics.firstIterations);
            cfg.localUpdateDiagnostics.everyIterations = value.value(
                "every_iterations", cfg.localUpdateDiagnostics.everyIterations);
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: local_update_diagnostics must be a boolean or object.");
        }
        if (cfg.localUpdateDiagnostics.firstIterations < 0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: local_update_diagnostics.first_iterations must be non-negative.");
        }
        if (cfg.localUpdateDiagnostics.everyIterations <= 0) {
            throw std::invalid_argument(
                "newtonConfigFromJson: local_update_diagnostics.every_iterations must be positive.");
        }
        if (cfg.localUpdateDiagnostics.enabled &&
            cfg.localUpdateDiagnostics.csvFile.empty()) {
            throw std::invalid_argument(
                "newtonConfigFromJson: enabled local_update_diagnostics requires csv_file.");
        }
    }
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
            cfg.carrierRowConvergence.floorLimitedAction = value.value(
                "floor_limited_action", cfg.carrierRowConvergence.floorLimitedAction);
            cfg.carrierRowConvergence.minSourceScaleFraction = value.value(
                "min_source_scale_fraction", cfg.carrierRowConvergence.minSourceScaleFraction);
            cfg.carrierRowConvergence.minSourceScale = value.value(
                "min_source_scale", cfg.carrierRowConvergence.minSourceScale);
            cfg.carrierRowConvergence.minSourceGlobalFraction = value.value(
                "min_source_global_fraction",
                cfg.carrierRowConvergence.minSourceGlobalFraction);
            cfg.carrierRowConvergence.minCarrierDensity_m3 = value.value(
                "min_carrier_density_m3", cfg.carrierRowConvergence.minCarrierDensity_m3);
            cfg.carrierRowConvergence.minFluxScaleFraction = value.value(
                "min_flux_scale_fraction", cfg.carrierRowConvergence.minFluxScaleFraction);
            cfg.carrierRowConvergence.minFluxScale = value.value(
                "min_flux_scale", cfg.carrierRowConvergence.minFluxScale);
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
        if (cfg.carrierRowConvergence.floorLimitedAction != "enforce" &&
            cfg.carrierRowConvergence.floorLimitedAction != "report")
            throw std::invalid_argument("newtonConfigFromJson: invalid floor_limited_action");
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
        if (cfg.carrierRowConvergence.minSourceGlobalFraction < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minSourceGlobalFraction)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_source_global_fraction must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.minCarrierDensity_m3 < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minCarrierDensity_m3)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_carrier_density_m3 must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.minFluxScaleFraction < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minFluxScaleFraction)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_flux_scale_fraction must be finite and nonnegative.");
        }
        if (cfg.carrierRowConvergence.minFluxScale < 0.0 ||
            !std::isfinite(cfg.carrierRowConvergence.minFluxScale)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: carrier_row_convergence.min_flux_scale must be finite and nonnegative.");
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
    if (json.contains("linear_equilibration")) {
        const auto& value = json.at("linear_equilibration");
        if (value.is_boolean()) {
            cfg.linearEquilibration.mode =
                value.get<bool>() ? "l2_row_column" : "off";
        } else if (value.is_string()) {
            cfg.linearEquilibration.mode = value.get<std::string>();
        } else if (value.is_object()) {
            cfg.linearEquilibration.mode = value.value(
                "mode", cfg.linearEquilibration.mode);
            if (value.contains("enabled") && !value.at("enabled").get<bool>())
                cfg.linearEquilibration.mode = "off";
        } else {
            throw std::invalid_argument(
                "newtonConfigFromJson: linear_equilibration must be a boolean, string, or object.");
        }
        if (cfg.linearEquilibration.mode != "off" &&
            cfg.linearEquilibration.mode != "l2_row_column") {
            throw std::invalid_argument(
                "newtonConfigFromJson: linear_equilibration.mode must be 'off' or 'l2_row_column'.");
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
            cfg.globalContinuityClosure.absoluteTolerance = value.value(
                "absolute_tolerance",
                cfg.globalContinuityClosure.absoluteTolerance);
            cfg.globalContinuityClosure.linearSystemStabilization = value.value(
                "linear_system_stabilization",
                cfg.globalContinuityClosure.linearSystemStabilization);
            cfg.globalContinuityClosure.feasibilityRestoration = value.value(
                "feasibility_restoration",
                cfg.globalContinuityClosure.feasibilityRestoration);
            cfg.globalContinuityClosure.restorationPhysicalRowsOnly = value.value(
                "restoration_physical_rows_only", cfg.globalContinuityClosure.restorationPhysicalRowsOnly);
            cfg.globalContinuityClosure.roundoffAudit = value.value(
                "roundoff_audit", cfg.globalContinuityClosure.roundoffAudit);
            cfg.globalContinuityClosure.auditJsonlFile = value.value(
                "audit_jsonl_file", cfg.globalContinuityClosure.auditJsonlFile);
            cfg.globalContinuityClosure.restorationJacobian = value.value(
                "restoration_jacobian", cfg.globalContinuityClosure.restorationJacobian);
            cfg.globalContinuityClosure.feasibilityRestorationFactor = value.value(
                "feasibility_restoration_factor",
                cfg.globalContinuityClosure.feasibilityRestorationFactor);
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
            !std::isfinite(cfg.globalContinuityClosure.sourceFloor) ||
            cfg.globalContinuityClosure.absoluteTolerance < 0.0 ||
            !std::isfinite(cfg.globalContinuityClosure.absoluteTolerance)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: invalid global_continuity_closure bounds.");
        }
        if (cfg.globalContinuityClosure.linearSystemStabilization != "off" &&
            cfg.globalContinuityClosure.linearSystemStabilization !=
                "replace_free_row_sum") {
            throw std::invalid_argument(
                "newtonConfigFromJson: global_continuity_closure.linear_system_stabilization must be 'off' or 'replace_free_row_sum'.");
        }
        if (cfg.globalContinuityClosure.linearSystemStabilization != "off" &&
            cfg.continuityRowScaling.enabled) {
            throw std::invalid_argument(
                "newtonConfigFromJson: global continuity row stabilization and continuity_row_scaling cannot be enabled together.");
        }
        if (cfg.globalContinuityClosure.feasibilityRestoration != "off" &&
            cfg.globalContinuityClosure.feasibilityRestoration !=
                "linear_qf_projection" &&
            cfg.globalContinuityClosure.feasibilityRestoration !=
                "final_linear_qf_projection") {
            throw std::invalid_argument(
                "newtonConfigFromJson: global_continuity_closure.feasibility_restoration must be 'off', 'linear_qf_projection', or 'final_linear_qf_projection'.");
        }
        if (cfg.globalContinuityClosure.restorationJacobian != "stabilized" &&
            cfg.globalContinuityClosure.restorationJacobian != "physical")
            throw std::invalid_argument("restoration_jacobian must be stabilized or physical.");
        if (!(cfg.globalContinuityClosure.feasibilityRestorationFactor > 0.0) ||
            cfg.globalContinuityClosure.feasibilityRestorationFactor > 1.0 ||
            !std::isfinite(
                cfg.globalContinuityClosure.feasibilityRestorationFactor)) {
            throw std::invalid_argument(
                "newtonConfigFromJson: global_continuity_closure.feasibility_restoration_factor must be in (0, 1].");
        }
    }
    cfg.finiteDifferenceStep = json.value("finite_difference_step", cfg.finiteDifferenceStep);
    cfg.jacobian = json.value("jacobian", cfg.jacobian);
    cfg.quasiFermiReference =
        json.value("quasi_fermi_reference", cfg.quasiFermiReference);
    cfg.poissonOnlyQuasiFermiInitialization = json.value(
        "poisson_only_quasi_fermi_initialization",
        cfg.poissonOnlyQuasiFermiInitialization);
    cfg.quasiFermiRecenterOnInitialState = json.value(
        "quasi_fermi_recenter_on_initial_state",
        cfg.quasiFermiRecenterOnInitialState);
    cfg.quasiFermiRecenterOnStall = json.value(
        "quasi_fermi_recenter_on_stall", cfg.quasiFermiRecenterOnStall);
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
    cfg.srhDopingDependence.densityCoupling = json.value(
        "srh_density_coupling", cfg.srhDopingDependence.densityCoupling);
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
    cfg.augerWithGeneration = json.value("auger_with_generation", cfg.augerWithGeneration);
    cfg.augerExcessProduct = json.value(
        "auger_excess_product", cfg.augerExcessProduct);
    cfg.recombinationVolumePolicy = json.value("recombination_volume_policy", cfg.recombinationVolumePolicy);
    if (cfg.recombinationVolumePolicy != "legacy_node" &&
        cfg.recombinationVolumePolicy != "material_local_mixed_voronoi" &&
        cfg.recombinationVolumePolicy != "material_local_circumcentric" &&
        cfg.recombinationVolumePolicy != "material_local_box")
        throw std::invalid_argument("Unsupported recombination_volume_policy");
    if (json.contains("recombination_materials"))
        cfg.recombinationMaterials=json.at("recombination_materials").get<std::vector<std::string>>();
    if (json.contains("auger_density_dependence"))
        cfg.augerDensityDependence = augerDensityDependenceConfigFromJson(
            json.at("auger_density_dependence"), scaling);
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
            cfg.bandgapNarrowing.equalNiFluxEvaluation = value.value(
                "equal_ni_flux_evaluation",
                cfg.bandgapNarrowing.equalNiFluxEvaluation);
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
            cfg.impactIonization.ialmobBoundaryFieldProjection = value.value(
                "ialmob_boundary_field_projection",
                cfg.impactIonization.ialmobBoundaryFieldProjection);
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
            cfg.impactIonization.sourceJacobianMode = value.value(
                "source_jacobian", cfg.impactIonization.sourceJacobianMode);
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
    const auto contactField = gradQfContactFieldFromString(
        json.value("grad_qf_contact_field", std::string{"electric_field"}));
    cfg.mobility.gradQfContactField = contactField;
    cfg.impactIonization.gradQfContactField = contactField;
    detail::validateImpactIonizationDrivingForce(cfg.impactIonization, "newtonConfigFromJson");
    (void)BandToBandTunnelingModel(cfg.bandToBand);

    if (cfg.jacobian != "analytic" && cfg.jacobian != "finite_difference")
        throw std::invalid_argument(
            "newtonConfigFromJson: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg.srhDopingDependence.densityCoupling != "quantum" &&
        cfg.srhDopingDependence.densityCoupling != "sentaurus_default" &&
        !cfg.srhDopingDependence.usesQuantumBandEdge()) {
        throw std::invalid_argument(
            "newtonConfigFromJson: srh_density_coupling must be 'quantum', "
            "'sentaurus_default', or 'quantum_band_edge'.");
    }
    if (cfg.srhDopingDependence.usesQuantumBandEdge() &&
        (!cfg.electronQuantumPotential.enabled ||
         cfg.electronQuantumPotential.transportCoupling != "direct_band_edge")) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quantum_band_edge SRH requires enabled "
            "electron quantum potential with direct_band_edge transport coupling.");
    }
    if (cfg.augerExcessProduct != "generalized_fermi" &&
        cfg.augerExcessProduct != "classical_np") {
        throw std::invalid_argument(
            "newtonConfigFromJson: auger_excess_product must be "
            "'generalized_fermi' or 'classical_np'.");
    }
    if (cfg.quasiFermiReference != "none" &&
        cfg.quasiFermiReference != "contact_majority" &&
        cfg.quasiFermiReference != "contact_basin" &&
        cfg.quasiFermiReference != "contact_nearest_state") {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_reference must be 'none', "
            "'contact_majority', 'contact_basin', or 'contact_nearest_state'.");
    }
    if (cfg.poissonOnlyQuasiFermiInitialization != "contact_basin" &&
        cfg.poissonOnlyQuasiFermiInitialization != "preserve_initial") {
        throw std::invalid_argument(
            "newtonConfigFromJson: poisson_only_quasi_fermi_initialization "
            "must be 'contact_basin' or 'preserve_initial'.");
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
    if (cfg.quasiFermiUpdateLimitMode != "componentwise_poisson_recorrect" &&
        cfg.quasiFermiUpdateLimitMode != "uniform_trust_region") {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_update_limit_mode must be "
            "'componentwise_poisson_recorrect' or 'uniform_trust_region'.");
    }
    if (cfg.quasiFermiTrustRegionGrowthFactor < 1.0 ||
        !std::isfinite(cfg.quasiFermiTrustRegionGrowthFactor)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_trust_region_growth_factor "
            "must be finite and at least one.");
    }
    if (cfg.quasiFermiTrustRegionMaxMultiplier < 1.0 ||
        !std::isfinite(cfg.quasiFermiTrustRegionMaxMultiplier)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_trust_region_max_multiplier "
            "must be finite and at least one.");
    }
    if (!(cfg.quasiFermiTrustRegionExpansionThreshold > 0.0) ||
        cfg.quasiFermiTrustRegionExpansionThreshold > 1.0 ||
        !std::isfinite(cfg.quasiFermiTrustRegionExpansionThreshold)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_trust_region_expansion_threshold "
            "must be finite and in (0,1].");
    }
    if (!(cfg.quasiFermiTrustRegionShrinkFactor > 0.0) ||
        cfg.quasiFermiTrustRegionShrinkFactor > 1.0 ||
        !std::isfinite(cfg.quasiFermiTrustRegionShrinkFactor)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_trust_region_shrink_factor "
            "must be finite and in (0,1].");
    }
    if (!(cfg.quasiFermiTrustRegionMinMultiplier > 0.0) ||
        cfg.quasiFermiTrustRegionMinMultiplier >
            cfg.quasiFermiTrustRegionMaxMultiplier ||
        !std::isfinite(cfg.quasiFermiTrustRegionMinMultiplier)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: quasi_fermi_trust_region_min_multiplier "
            "must be finite, positive, and no greater than the maximum multiplier.");
    }
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
    if (cfg.blockAbsoluteConvergence.mode != "off" &&
        cfg.blockAbsoluteConvergence.mode != "enforce") {
        throw std::invalid_argument(
            "newtonConfigFromJson: block_absolute_convergence.mode must be 'off' or 'enforce'.");
    }
    const auto validateBlockCeiling = [](Real value, const char* name) {
        if (value < 0.0 || !std::isfinite(value)) {
            throw std::invalid_argument(
                std::string("newtonConfigFromJson: block_absolute_convergence.") +
                name + " must be non-negative and finite.");
        }
    };
    validateBlockCeiling(
        cfg.blockAbsoluteConvergence.psiResidualCeiling,
        "psi_residual_ceiling");
    validateBlockCeiling(
        cfg.blockAbsoluteConvergence.electronResidualCeiling,
        "electron_residual_ceiling");
    validateBlockCeiling(
        cfg.blockAbsoluteConvergence.holeResidualCeiling,
        "hole_residual_ceiling");
    if (cfg.blockAbsoluteConvergence.mode == "enforce" &&
        (!(cfg.blockAbsoluteConvergence.psiResidualCeiling > 0.0) ||
         !(cfg.blockAbsoluteConvergence.electronResidualCeiling > 0.0) ||
         !(cfg.blockAbsoluteConvergence.holeResidualCeiling > 0.0))) {
        throw std::invalid_argument(
            "newtonConfigFromJson: enforced block_absolute_convergence requires positive psi, electron, and hole ceilings.");
    }
    if (cfg.contactMajorityQfBranchDropLimit_V < 0.0 ||
        !std::isfinite(cfg.contactMajorityQfBranchDropLimit_V)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: contact_majority_qf_branch_drop_limit_V must be non-negative and finite.");
    }
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
    if (cfg.localUpdateDiagnostics.firstIterations < 0)
        throw std::invalid_argument(
            "newtonConfigFromJson: local_update_diagnostics.first_iterations must be non-negative.");
    if (cfg.localUpdateDiagnostics.everyIterations <= 0)
        throw std::invalid_argument(
            "newtonConfigFromJson: local_update_diagnostics.every_iterations must be positive.");
    if (cfg.localUpdateDiagnostics.enabled && cfg.localUpdateDiagnostics.csvFile.empty())
        throw std::invalid_argument(
            "newtonConfigFromJson: enabled local_update_diagnostics requires csv_file.");
    if (cfg.finiteDifferenceStep <= 0.0 || !std::isfinite(cfg.finiteDifferenceStep))
        throw std::invalid_argument(
            "newtonConfigFromJson: finite_difference_step must be positive and finite.");
    if (cfg.residualNorm != "block" && cfg.residualNorm != "l2")
        throw std::invalid_argument(
            "newtonConfigFromJson: residual_norm must be 'block' or 'l2'.");
    if (cfg.lineSearchMode != "merit" &&
        cfg.lineSearchMode != "block_filter" &&
        cfg.lineSearchMode != "block_closure_filter")
        throw std::invalid_argument(
            "newtonConfigFromJson: line_search_mode must be 'merit', "
            "'block_filter', or 'block_closure_filter'.");
    if (!(cfg.residualFilterGamma > 0.0) ||
        !(cfg.residualFilterGamma < 1.0) ||
        !std::isfinite(cfg.residualFilterGamma)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: residual_filter_gamma must be finite and in (0,1).");
    }
    if (cfg.residualFilterEnvelopeFactor < 1.0 ||
        !std::isfinite(cfg.residualFilterEnvelopeFactor)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: residual_filter_envelope_factor must be finite and at least one.");
    }
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
        cfg.electronQuantumPotential.couplingMode != "frozen" &&
        cfg.electronQuantumPotential.couplingMode != "monolithic_box") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential.coupling_mode "
            "must be 'outer', 'frozen', or 'monolithic_box'.");
    }
    if (cfg.electronQuantumPotential.transportCoupling != "direct_band_edge" &&
        cfg.electronQuantumPotential.transportCoupling !=
            "sentaurus_exponential") {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "transport_coupling must be 'direct_band_edge' or "
            "'sentaurus_exponential'.");
    }
    if (!(cfg.electronQuantumPotential.transportCouplingWeight >= 0.0) ||
        !(cfg.electronQuantumPotential.transportCouplingWeight <= 1.0) ||
        !std::isfinite(
            cfg.electronQuantumPotential.transportCouplingWeight)) {
        throw std::invalid_argument(
            "newtonConfigFromJson: electron_quantum_potential."
            "transport_coupling_weight must be finite and in [0,1].");
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
    if(cfg.coupledOperator=="bvmethods_tri3_qf") {
        if(cfg.splitDDState || cfg.electronQuantumPotential.enabled)
            throw std::invalid_argument("BVMethods profile requires ordinary classical QF coordinates");
        if (!usesFermiDirac(cfg.carrierStatistics))
            throw std::invalid_argument("BVMethods profile requires explicit Fermi statistics");
    } else if(cfg.carrierStatistics.model=="fermi_dirac_joyce_dixon" &&
              (cfg.coupledOperator!="legacy" || cfg.stateVariables!="potential_density" ||
               cfg.carrierStatistics.transport!="effective_band_sg")) {
        throw std::invalid_argument(
            "Legacy Joyce-Dixon Newton requires potential_density with effective_band_sg");
    }

    resolveFermiSgDefault(cfg);
    validatePopulationNumerics(cfg);
    validateNumericalPrecisionOptions(cfg);
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
    resolvePopulationDefaults(cfg_);
    resolveFermiSgDefault(cfg_);
    validateNumericalPrecisionOptions(cfg_);
    if (cfg_.linearFactorRefinementSteps < 0 || cfg_.linearFactorRefinementSteps > 3 ||
        (cfg_.linearFactorRefinementSteps && (cfg_.linearRefinementIterations || cfg_.splitDDState ||
         cfg_.continuityRowScaling.enabled || cfg_.linearEquilibration.mode != "off")))
        throw std::invalid_argument("guarded factor refinement requires a separate unscaled three-block Newton arm");
    if((cfg_.populationSgCurrent!="quasi_fermi" && cfg_.populationSgCurrent!="band_population") ||
       (cfg_.populationSgCurrent!="quasi_fermi" && cfg_.stateVariables!="potential_density"))
        throw std::invalid_argument("invalid population SG arithmetic mode");
    if(cfg_.stateVariables!="potential_quasi_fermi" && cfg_.stateVariables!="potential_density")
        throw std::invalid_argument("unknown state_variables");
    if(cfg_.populationSgBootstrapLinearization!="unchanged" &&
       (cfg_.populationSgBootstrapLinearization!="floor_adjacent_edges" ||
        cfg_.stateVariables!="potential_density" || cfg_.coupledOperator!="bvmethods_tri3_qf" ||
        cfg_.populationSgCurrent!="quasi_fermi" || !(cfg_.populationMinimumDensity_cm3>0.)))
        throw std::invalid_argument("population_sg_bootstrap_linearization floor_adjacent_edges requires BVMethods potential_density, quasi_fermi SG and positive minimum_density_cm3");
    if((cfg_.populationDerivativeRecovery!="full"&&cfg_.populationDerivativeRecovery!="automatic"&&cfg_.populationDerivativeRecovery!="after_protection") ||
       (cfg_.populationDerivativeRecovery!="full"&&cfg_.stateVariables!="potential_density"))
        throw std::invalid_argument("invalid population derivative recovery mode");
    validatePopulationNumerics(cfg_);
    if (!std::isfinite(cfg_.legacyTransportJacobianStep) || cfg_.legacyTransportJacobianStep <= 0.0 ||
        !std::isfinite(cfg_.legacyAvalancheJacobianStep) || cfg_.legacyAvalancheJacobianStep <= 0.0)
        throw std::invalid_argument("NewtonSolver: legacy edge Jacobian steps must be finite and positive");
    if (cfg_.linearRefinementIterations < 0 || cfg_.linearRefinementIterations > 10)
        throw std::invalid_argument("NewtonSolver: linear_refinement_iterations must be between 0 and 10.");
    if (cfg_.jacobian != "analytic" && cfg_.jacobian != "finite_difference")
        throw std::invalid_argument(
            "NewtonSolver: jacobian must be 'analytic' or 'finite_difference'.");
    if (cfg_.quasiFermiReference != "none" &&
        cfg_.quasiFermiReference != "contact_majority" &&
        cfg_.quasiFermiReference != "contact_basin" &&
        cfg_.quasiFermiReference != "contact_nearest_state") {
        throw std::invalid_argument(
            "NewtonSolver: quasi_fermi_reference must be 'none', "
            "'contact_majority', 'contact_basin', or 'contact_nearest_state'.");
    }
    if (cfg_.poissonOnlyQuasiFermiInitialization != "contact_basin" &&
        cfg_.poissonOnlyQuasiFermiInitialization != "preserve_initial") {
        throw std::invalid_argument(
            "NewtonSolver: poisson_only_quasi_fermi_initialization must be "
            "'contact_basin' or 'preserve_initial'.");
    }
    if (cfg_.quasiFermiRecenterOnInitialState &&
        (!cfg_.warmStart || cfg_.quasiFermiReference == "none")) {
        throw std::invalid_argument(
            "NewtonSolver: quasi_fermi_recenter_on_initial_state requires "
            "warm_start=true and a non-'none' quasi_fermi_reference.");
    }
    if (cfg_.residualNorm != "block" && cfg_.residualNorm != "l2")
        throw std::invalid_argument(
            "NewtonSolver: residual_norm must be 'block' or 'l2'.");
    if (cfg_.linearEquilibration.mode != "off" &&
        cfg_.linearEquilibration.mode != "l2_row_column") {
        throw std::invalid_argument(
            "NewtonSolver: linear_equilibration.mode must be 'off' or 'l2_row_column'.");
    }
    if (cfg_.lineSearchMode != "merit" &&
        cfg_.lineSearchMode != "block_filter" &&
        cfg_.lineSearchMode != "block_closure_filter")
        throw std::invalid_argument(
            "NewtonSolver: line_search_mode must be 'merit', 'block_filter', "
            "or 'block_closure_filter'.");
    if (cfg_.quasiFermiUpdateLimitMode != "componentwise_poisson_recorrect" &&
        cfg_.quasiFermiUpdateLimitMode != "uniform_trust_region") {
        throw std::invalid_argument(
            "NewtonSolver: quasi_fermi_update_limit_mode must be "
            "'componentwise_poisson_recorrect' or 'uniform_trust_region'.");
    }
    if (cfg_.quasiFermiTrustRegionGrowthFactor < 1.0 ||
        !std::isfinite(cfg_.quasiFermiTrustRegionGrowthFactor) ||
        cfg_.quasiFermiTrustRegionMaxMultiplier < 1.0 ||
        !std::isfinite(cfg_.quasiFermiTrustRegionMaxMultiplier) ||
        !(cfg_.quasiFermiTrustRegionExpansionThreshold > 0.0) ||
        cfg_.quasiFermiTrustRegionExpansionThreshold > 1.0 ||
        !std::isfinite(cfg_.quasiFermiTrustRegionExpansionThreshold) ||
        !(cfg_.quasiFermiTrustRegionShrinkFactor > 0.0) ||
        cfg_.quasiFermiTrustRegionShrinkFactor > 1.0 ||
        !std::isfinite(cfg_.quasiFermiTrustRegionShrinkFactor) ||
        !(cfg_.quasiFermiTrustRegionMinMultiplier > 0.0) ||
        cfg_.quasiFermiTrustRegionMinMultiplier >
            cfg_.quasiFermiTrustRegionMaxMultiplier ||
        !std::isfinite(cfg_.quasiFermiTrustRegionMinMultiplier)) {
        throw std::invalid_argument(
            "NewtonSolver: invalid quasi-Fermi trust-region expansion controls.");
    }
    if (!(cfg_.residualFilterGamma > 0.0) ||
        !(cfg_.residualFilterGamma < 1.0) ||
        !std::isfinite(cfg_.residualFilterGamma) ||
        cfg_.residualFilterEnvelopeFactor < 1.0 ||
        !std::isfinite(cfg_.residualFilterEnvelopeFactor)) {
        throw std::invalid_argument(
            "NewtonSolver: invalid block residual-filter controls.");
    }
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
    if (cfg_.blockAbsoluteConvergence.mode != "off" &&
        cfg_.blockAbsoluteConvergence.mode != "enforce") {
        throw std::invalid_argument(
            "NewtonSolver: block_absolute_convergence.mode must be 'off' or 'enforce'.");
    }
    const auto validBlockCeiling = [](Real value) {
        return value >= 0.0 && std::isfinite(value);
    };
    if (!validBlockCeiling(cfg_.blockAbsoluteConvergence.psiResidualCeiling) ||
        !validBlockCeiling(cfg_.blockAbsoluteConvergence.electronResidualCeiling) ||
        !validBlockCeiling(cfg_.blockAbsoluteConvergence.holeResidualCeiling) ||
        (cfg_.blockAbsoluteConvergence.mode == "enforce" &&
         (!(cfg_.blockAbsoluteConvergence.psiResidualCeiling > 0.0) ||
          !(cfg_.blockAbsoluteConvergence.electronResidualCeiling > 0.0) ||
          !(cfg_.blockAbsoluteConvergence.holeResidualCeiling > 0.0)))) {
        throw std::invalid_argument(
            "NewtonSolver: invalid block_absolute_convergence ceilings.");
    }
    if (cfg_.contactMajorityQfBranchDropLimit_V < 0.0 ||
        !std::isfinite(cfg_.contactMajorityQfBranchDropLimit_V)) {
        throw std::invalid_argument(
            "NewtonSolver: contact_majority_qf_branch_drop_limit_V must be non-negative and finite.");
    }
    for (const std::string& requestedContact :
         cfg_.contactMajorityQfBranchGuardContacts) {
        const bool exists = std::any_of(
            mesh_.contacts().begin(), mesh_.contacts().end(),
            [&](const Contact& contact) {
                return contact.name == requestedContact;
            });
        if (!exists) {
            throw std::invalid_argument(
                "NewtonSolver: contact_majority_qf_branch_guard_contacts contains unknown contact '" +
                requestedContact + "'.");
        }
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
    if (cfg_.localUpdateDiagnostics.firstIterations < 0 ||
        cfg_.localUpdateDiagnostics.everyIterations <= 0 ||
        (cfg_.localUpdateDiagnostics.enabled &&
         cfg_.localUpdateDiagnostics.csvFile.empty())) {
        throw std::invalid_argument(
            "NewtonSolver: invalid local_update_diagnostics configuration.");
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
    scaling.coupledOperator = cfg_.coupledOperator;
    scaling.bvMethodsNumerics = cfg_.bvMethodsNumerics;
    scaling.stateVariables = cfg_.stateVariables;
    scaling.populationSgCurrent = cfg_.populationSgCurrent;
    scaling.populationSgBootstrapFloor_cm3 = cfg_.populationSgBootstrapLinearization=="floor_adjacent_edges"?cfg_.populationMinimumDensity_cm3:0.;
    scaling.regionResolvedInterfaceAssembly =
        cfg_.regionResolvedInterfaceAssembly;
    scaling.poissonChargeVolumePolicy = cfg_.poissonChargeVolumePolicy;
    scaling.poissonPerCellPermittivity = cfg_.poissonPerCellPermittivity;
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

void NewtonSolver::configureAssemblerNumerics(
    CoupledDDAssembler& assembler, const DDSolution* state) const
{
    assembler.setBvMethodsEvaluationReuse(cfg_.diagnosticBvMethodsEvaluationReuse);
    assembler.enableSplitDDState(cfg_.splitDDState);
    assembler.setExtendedPoissonResidual(cfg_.poissonResidualPrecision == "binary128");
    assembler.setDiagnosticFreezeTransportQfMobility(cfg_.diagnosticFreezeTransportQfMobility);
    assembler.setLegacyEdgeJacobianSteps(
        cfg_.legacyTransportJacobianStep, cfg_.legacyAvalancheJacobianStep);
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
    if (cfg_.quasiFermiReference == "contact_majority") {
        assembler.setQuasiFermiReferences(electronReference, holeReference);
        return;
    }

    // Geometry alone can assign a source reference to a drain-connected
    // accumulation layer. Its tiny current-carrying QF differences then live
    // on top of the entire terminal voltage and cannot be resolved in Real.
    // This opt-in coordinate choice selects the closest terminal voltage to
    // each supplied QF. References remain fixed for the ensuing Newton solve.
    // Without a supplied state, use the geometric basins below for cold start.
    if (cfg_.quasiFermiReference == "contact_nearest_state" && state) {
        const int nodeCount = static_cast<int>(mesh_.numNodes());
        if (state->phin.size() != nodeCount || state->phip.size() != nodeCount)
            throw std::invalid_argument("contact_nearest_state: state size does not match mesh");
        std::vector<Real> anchors;
        for (const Contact& contact : mesh_.contacts()) {
            const auto bias = contactBiases_.find(contact.name);
            const auto spec = contactSpecs_.find(contact.name);
            if (bias == contactBiases_.end() || contact.node_ids.empty() ||
                (spec != contactSpecs_.end() && spec->second.type == ContactType::MetalGate))
                continue;
            anchors.push_back(bias->second);
        }
        if (anchors.empty()) anchors.push_back(0.0);
        const auto referenceField = [&](bool electron) {
            VectorXd references(nodeCount);
            const bool referenced = electron
                ? state->hasReferencedElectronQuasiFermi()
                : state->hasReferencedHoleQuasiFermi();
            for (int i = 0; i < nodeCount; ++i) {
                const Real oldReference = !referenced ? 0.0 : (electron
                    ? state->electronQuasiFermiReferenceAt(i)
                    : state->holeQuasiFermiReferenceAt(i));
                const Real increment = referenced ? (electron
                    ? state->phinIncrement(i) : state->phipIncrement(i))
                    : (electron ? state->phin(i) : state->phip(i));
                if (!std::isfinite(oldReference) || !std::isfinite(increment))
                    throw std::invalid_argument("contact_nearest_state: non-finite quasi-Fermi state");
                Real selected = anchors.front();
                long double distance = std::numeric_limits<long double>::infinity();
                for (Real anchor : anchors) {
                    const long double candidate = std::abs(
                        (static_cast<long double>(oldReference) - anchor) + increment);
                    if (candidate < distance) {
                        distance = candidate;
                        selected = anchor;
                    }
                }
                references(i) = selected;
            }
            return references;
        };
        assembler.setQuasiFermiReferenceFields(referenceField(true), referenceField(false));
        return;
    }

    struct ReferenceSeed {
        Index node = 0;
        Real reference_V = 0.0;
        std::size_t rank = 0;
    };
    std::vector<ReferenceSeed> electronSeeds;
    std::vector<ReferenceSeed> holeSeeds;
    std::size_t seedRank = 0;
    for (const Contact& contact : mesh_.contacts()) {
        const auto biasIt = contactBiases_.find(contact.name);
        if (biasIt == contactBiases_.end() || contact.node_ids.empty())
            continue;
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        for (Index node : contact.node_ids) {
            // A shorted source/body contact may contain both n+ and p+
            // segments. Its mean doping must not erase either carrier's
            // reference seeds. Follow the actual Ohmic boundary doping rule.
            const Real nodeNetDoping = ohmicContactNetDoping(
                mesh_, doping_, contact, node,
                cfg_.contactBoundaryReconstruction == "dominant_signed_contact_mean");
            if (nodeNetDoping > 0.0)
                electronSeeds.push_back({node, biasIt->second, seedRank});
            if (nodeNetDoping < 0.0)
                holeSeeds.push_back({node, biasIt->second, seedRank});
        }
        ++seedRank;
    }

    const auto& edgeCells = detail::buildEdgeCellMap(mesh_);
    std::vector<std::vector<std::pair<Index, Real>>> adjacency(mesh_.numNodes());
    for (Index edgeId = 0; edgeId < mesh_.numEdges(); ++edgeId) {
        bool activeTransport = false;
        for (Index cellId : edgeCells[edgeId]) {
            const Cell& cell = mesh_.getCell(cellId);
            const Region& region = mesh_.getRegion(cell.region_id);
            if (detail::isTransportMaterial(
                    matdb_.getMaterial(region.material, cfg_.temperature_K))) {
                activeTransport = true;
                break;
            }
        }
        if (!activeTransport)
            continue;
        const Edge& edge = mesh_.getEdge(edgeId);
        adjacency[edge.n0].push_back({edge.n1, edge.length});
        adjacency[edge.n1].push_back({edge.n0, edge.length});
    }

    const auto basinField = [&](const std::vector<ReferenceSeed>& seeds,
                                Real fallbackReference) {
        const int nodeCount = static_cast<int>(mesh_.numNodes());
        VectorXd references = VectorXd::Constant(nodeCount, fallbackReference);
        if (seeds.empty())
            return references;
        const Real infinity = std::numeric_limits<Real>::infinity();
        std::vector<Real> distance(static_cast<std::size_t>(nodeCount), infinity);
        std::vector<std::size_t> owner(
            static_cast<std::size_t>(nodeCount),
            std::numeric_limits<std::size_t>::max());
        using QueueEntry = std::tuple<Real, std::size_t, Index, Real>;
        std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                            std::greater<QueueEntry>> queue;
        for (const ReferenceSeed& seed : seeds) {
            const std::size_t node = static_cast<std::size_t>(seed.node);
            if (0.0 < distance[node] ||
                (distance[node] == 0.0 && seed.rank < owner[node])) {
                distance[node] = 0.0;
                owner[node] = seed.rank;
                references(static_cast<int>(seed.node)) = seed.reference_V;
                queue.emplace(0.0, seed.rank, seed.node, seed.reference_V);
            }
        }
        while (!queue.empty()) {
            const auto [currentDistance, currentOwner, node, reference] =
                queue.top();
            queue.pop();
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            if (currentDistance != distance[nodeIndex] ||
                currentOwner != owner[nodeIndex]) {
                continue;
            }
            for (const auto& [neighbor, length] : adjacency[node]) {
                const Real candidate = currentDistance + length;
                const std::size_t neighborIndex =
                    static_cast<std::size_t>(neighbor);
                if (candidate < distance[neighborIndex] ||
                    (candidate == distance[neighborIndex] &&
                     currentOwner < owner[neighborIndex])) {
                    distance[neighborIndex] = candidate;
                    owner[neighborIndex] = currentOwner;
                    references(static_cast<int>(neighbor)) = reference;
                    queue.emplace(candidate, currentOwner, neighbor, reference);
                }
            }
        }
        for (int i = 0; i < nodeCount; ++i) {
            if (std::isfinite(distance[static_cast<std::size_t>(i)]))
                continue;
            const Node& node = mesh_.getNode(static_cast<Index>(i));
            const ReferenceSeed* nearest = &seeds.front();
            Real nearestSquared = std::numeric_limits<Real>::infinity();
            for (const ReferenceSeed& seed : seeds) {
                const Node& seedNode = mesh_.getNode(seed.node);
                const Real dx = node.x - seedNode.x;
                const Real dy = node.y - seedNode.y;
                const Real squared = dx * dx + dy * dy;
                if (squared < nearestSquared ||
                    (squared == nearestSquared && seed.rank < nearest->rank)) {
                    nearestSquared = squared;
                    nearest = &seed;
                }
            }
            references(i) = nearest->reference_V;
        }
        return references;
    };

    assembler.setQuasiFermiReferenceFields(
        basinField(electronSeeds, electronReference),
        basinField(holeSeeds, holeReference));
}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler,
    const std::unordered_map<std::string, Real>& contactBiases,
    const std::unordered_map<std::string, Real>* contactLowOverride) const
{
    CoupledDDBoundaryConditions bcs;
    const auto& ni = assembler.intrinsicDensity();
    const auto& Nc = assembler.electronDensityOfStates();
    const auto& Nv = assembler.holeDensityOfStates();
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
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

    const auto& contactLows=contactLowOverride?*contactLowOverride:cfg_.splitContactBiasLow_V;
    for(const auto& [name,low]:contactLows) {
        if(!cfg_.splitDDState || !contactBiases.contains(name) || !std::isfinite(low) ||
           contactBiases.at(name)+low!=contactBiases.at(name) || enableMinorityElectronRelaxation)
            throw std::invalid_argument("Contact low part needs a normalized split Dirichlet voltage");
        const auto spec=contactSpecs_.find(name);
        if(spec!=contactSpecs_.end() && spec->second.type==ContactType::Schottky)
            throw std::invalid_argument("Split contact voltage does not support Schottky boundaries");
    }

    for (Index c = 0; c < mesh_.numContacts(); ++c) {
        const Contact& contact = mesh_.getContact(c);
        auto it = contactBiases.find(contact.name);
        if (it == contactBiases.end()) continue;

        const double Vbias = it->second;
        const double VbiasLow=contactLows.contains(contact.name)?contactLows.at(contact.name):0.;
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            for (Index nid : contact.node_ids) {
                bcs.psi[nid] = Vbias / potentialScale;
                if(cfg_.splitDDState)setSplitBoundaryTarget(bcs,0,nid,0.,Vbias,potentialScale,VbiasLow);
            }
            continue;
        }
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::Schottky) {
            for (Index nid : contact.node_ids)
                if (assembler.intrinsicReferencePotential()[nid] != 0.0)
                    throw std::invalid_argument(
                        "Nonzero intrinsic reference with Schottky contacts is unsupported");
            ContactBoundarySpec spec = specIt->second;
            spec.bias = Vbias;
            if (spec.emissionModel != "thermionic_robin") {
                throw std::runtime_error(
                    "NewtonSolver: Schottky contact '" + contact.name +
                    "' requires emission_model='thermionic_robin'.");
            }

            const Real electronVelocity =
                spec.electronSurfaceRecombinationVelocity.value_or(
                    spec.surfaceRecombinationVelocity.value_or(0.0));
            const Real holeVelocity =
                spec.holeSurfaceRecombinationVelocity.value_or(
                    spec.surfaceRecombinationVelocity.value_or(0.0));
            if (!(electronVelocity > 0.0) || !(holeVelocity > 0.0)) {
                throw std::runtime_error(
                    "NewtonSolver: thermionic Schottky contact '" + contact.name +
                    "' requires positive electron and hole surface velocities.");
            }

            std::optional<Material> material;
            if (contact.region_id < mesh_.numRegions()) {
                const Region& region = mesh_.getRegion(contact.region_id);
                if (matdb_.hasMaterial(region.material))
                    material = matdb_.getMaterial(region.material, cfg_.temperature_K);
            }
            if (!material) {
                for (Index cellId = 0; cellId < mesh_.numCells() && !material; ++cellId) {
                    const Cell& cell = mesh_.getCell(cellId);
                    const bool touches = std::any_of(
                        cell.node_ids.begin(), cell.node_ids.end(), [&](Index node) {
                            return std::find(contact.node_ids.begin(), contact.node_ids.end(), node)
                                != contact.node_ids.end();
                        });
                    if (!touches)
                        continue;
                    const Region& region = mesh_.getRegion(cell.region_id);
                    if (matdb_.hasMaterial(region.material))
                        material = matdb_.getMaterial(region.material, cfg_.temperature_K);
                }
            }
            const Real bandgap_eV = material && material->bandgap_eV
                ? *material->bandgap_eV
                : std::numeric_limits<Real>::quiet_NaN();
            const Real affinity_eV = material && material->electron_affinity_eV
                ? *material->electron_affinity_eV
                : std::numeric_limits<Real>::quiet_NaN();

            std::unordered_set<Index> contactNodes(
                contact.node_ids.begin(), contact.node_ids.end());
            std::unordered_map<Index, Real> nodalMeasure;
            for (Index node : contact.node_ids)
                nodalMeasure[node] = 0.0;
            for (Index edgeId = 0; edgeId < mesh_.numEdges(); ++edgeId) {
                const Edge& edge = mesh_.getEdge(edgeId);
                if (!contactNodes.contains(edge.n0) || !contactNodes.contains(edge.n1))
                    continue;
                nodalMeasure[edge.n0] += 0.5 * edge.length;
                nodalMeasure[edge.n1] += 0.5 * edge.length;
            }

            for (Index nid : contact.node_ids) {
                if (!(nodalMeasure[nid] > 0.0)) {
                    throw std::runtime_error(
                        "NewtonSolver: thermionic Schottky contact '" + contact.name +
                        "' has a node with zero boundary measure.");
                }
                const ContactState state = computeSchottkyContactState(
                    spec, ni[nid], Nc[nid], Nv[nid], bandgap_eV, affinity_eV,
                    doping_.netDoping(nid), cfg_.temperature_K);
                bcs.psi[nid] = state.psi / potentialScale;
                bcs.thermionic[nid] = CoupledDDThermionicBoundary{
                    electronVelocity,
                    holeVelocity,
                    state.n,
                    state.p,
                    nodalMeasure[nid],
                };
            }
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
            const double psiBuiltIn = equilibrium.potential + assembler.intrinsicReferencePotential()[nid];

            bcs.psi[nid] = (Vbias + psiBuiltIn) / potentialScale;
            if(cfg_.splitDDState) {
                setSplitBoundaryTarget(bcs,0,nid,psiBuiltIn,Vbias,potentialScale,VbiasLow);
                const Real enBias=(Ndop<0 && relaxMinorityOnPContact)?relaxedMinorityBias:Vbias;
                const Real hpBias=(Ndop>=0 && relaxMinorityOnNContact)?relaxedMinorityBias:Vbias;
                setSplitBoundaryTarget(bcs,1,nid,0.,enBias,potentialScale,VbiasLow);
                setSplitBoundaryTarget(bcs,2,nid,0.,hpBias,potentialScale,VbiasLow);
                continue;
            }
            if (cfg_.compensatedPopulationState ||
                cfg_.electronQuantumPotential.couplingMode == "monolithic_box") {
                bcs.ohmicPopulation.emplace(nid, std::array<Real,3>{psiBuiltIn, equilibrium.n, equilibrium.p});
                bcs.ohmicBias.emplace(nid, Vbias);
            }
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
    const Real Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
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
        sol.psi(i) = equilibrium.potential + assembler.intrinsicReferencePotential()[nid];
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
            assembler.electronReferencePrefactor()[nid], Nc[nid], sol.psi(i), sol.phin(i), Vt,
            cfg_.carrierStatistics);
        sol.p(i) = holeDensity(
            assembler.holeReferencePrefactor()[nid], Nv[nid], sol.psi(i), sol.phip(i), Vt,
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
    sol.electronQfReference = assembler.electronQuasiFermiReferenceField();
    sol.holeQfReference = assembler.holeQuasiFermiReferenceField();
    sol.n = assembler.usesPopulationState()?state.n:assembler.electronDensity(x);
    sol.p = assembler.usesPopulationState()?state.p:assembler.holeDensity(x);
    if(assembler.usesPopulationState()) {
        sol.stateVariables="potential_density";
        sol.populationPsi=x.head(N);sol.populationPotentialScale_V=potentialScale;
        const auto physical=assembler.populationPhysicalFields(x);
        sol.phin=physical.phin_V;sol.phip=physical.phip_V;
        sol.phinIncrement.resize(0);sol.phipIncrement.resize(0);
        sol.electronQfReference.resize(0);sol.holeQfReference.resize(0);
    }
    sol.iters = iters;
    if(assembler.usesExtendedPoissonResidual()) {
        sol.packedState=x;
        sol.packedPotentialScale_V=potentialScale;
    }
    if(auto runtime=assembler.splitRuntime())runtime->save(x,sol);
    return sol;
}

std::shared_ptr<CoupledDDAssembler> NewtonSolver::makeArclengthAssembler(
    const DDSolution* referenceState) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
        cfg_.carrierStatistics,
        cfg_.electronQuantumPotential.couplingMode == "monolithic_box"
            ? cfg_.electronQuantumPotential : DensityGradientQuantumPotentialConfig{});
    configureAssemblerNumerics(*assembler, referenceState);
    assembler->setFermiNodeCache(cfg_.diagnosticFermiNodeCache);
    assembler->setStructureCache(cfg_.sequentialJacobianStructure);
    return assembler;
}

ArclengthSystem NewtonSolver::makeArclengthSystem(const std::string& activeContact,
                                                  Real biasFiniteDifferenceStep_V,
                                                  bool positiveDensityTrials,
                                                  const DDSolution* referenceState,
                                                  bool boundaryOnlyParameterDerivative) const
{
    if (positiveDensityTrials &&
        (usesFermiDirac(cfg_.carrierStatistics) || cfg_.electronQuantumPotential.enabled)) {
        throw std::invalid_argument(
            "NewtonSolver::makeArclengthSystem: positive density trials require "
            "classical Boltzmann carrier statistics without quantum potential.");
    }
    if (cfg_.compensatedPopulationState)
        throw std::invalid_argument("Compensated population state requires voltage-controlled Newton");
    if(cfg_.splitDDState)throw std::invalid_argument("Split DD arclength updates are not qualified");
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

    auto assembler = makeArclengthAssembler(referenceState);
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
        [self, assembler, biasesAt, h, boundaryOnlyParameterDerivative, nodeCount](const VectorXd& x, Real lambda) {
            const CoupledDDBoundaryConditions bcsPlus =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda + h));
            const CoupledDDBoundaryConditions bcsMinus =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda - h));
            if(boundaryOnlyParameterDerivative && assembler->usesBvMethodsOperator() && !assembler->usesPopulationState() &&
               bcsPlus.thermionic.empty() && bcsMinus.thermionic.empty() &&
               assembler->electronQuasiFermiReferenceField().isZero(0.) &&
               assembler->holeQuasiFermiReferenceField().isZero(0.)) {
                const auto sameKeys=[](const auto& a,const auto& b){
                    if(a.size()!=b.size())return false;
                    for(const auto& [node,value]:a)if(!b.contains(node))return false;
                    return true;
                };
                if(sameKeys(bcsPlus.psi,bcsMinus.psi)&&sameKeys(bcsPlus.phin,bcsMinus.phin)&&sameKeys(bcsPlus.phip,bcsMinus.phip)) {
                    // Keep one full finite-state check. Only the Dirichlet
                    // replacements depend on bias in this explicit profile.
                    // Subtract the original x-value expressions, rather than
                    // changing rounding by differentiating boundary values.
                    const auto base=assembler->residual(x,self->buildBoundaryConditions(*assembler,biasesAt(lambda)));
                    if(base.allFinite()) {
                        VectorXd plus=base,minus=base;
                        const auto replace=[&](const auto& a,const auto& b,int offset){
                            for(const auto& [node,value]:a)plus[offset+node]=x[offset+node]-value;
                            for(const auto& [node,value]:b)minus[offset+node]=x[offset+node]-value;
                        };
                        replace(bcsPlus.psi,bcsMinus.psi,0);replace(bcsPlus.phin,bcsMinus.phin,nodeCount);
                        replace(bcsPlus.phip,bcsMinus.phip,2*nodeCount);
                        return VectorXd((plus-minus)/(2.0*h));
                    }
                }
            }
            const VectorXd fPlus = assembler->residual(x, bcsPlus);
            const VectorXd fMinus = assembler->residual(x, bcsMinus);
            return VectorXd((fPlus - fMinus) / (2.0 * h));
        };
    system.jacobian =
        [self, assembler, biasesAt](const VectorXd& x, Real lambda) {
            const CoupledDDBoundaryConditions bcs =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda));
            return assembler->assembleJacobian(x, bcs);
        };
    system.solveJacobian =
        [self, assembler, biasesAt](const VectorXd& x, Real lambda,
                                    const VectorXd& b, VectorXd& y) {
            const CoupledDDBoundaryConditions bcs =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda));
            const SparseMatrixd jacobian = assembler->assembleJacobian(x, bcs);
            try {
                LinearSolver linear("sparselu");
                y = linear.solve(jacobian,b);
                return y.allFinite();
            } catch (const std::exception&) { return false; }
        };
    system.solveJacobianPair =
        [self, assembler, biasesAt](const VectorXd& x,
                                    Real lambda,
                                    const VectorXd& b1,
                                    const VectorXd& b2,
                                    VectorXd& y1,
                                    VectorXd& y2) {
            const CoupledDDBoundaryConditions bcs =
                self->buildBoundaryConditions(*assembler, biasesAt(lambda));
            const SparseMatrixd jacobian = assembler->assembleJacobian(x, bcs);
            try {
                LinearSolver linear("sparselu");
                y1 = linear.solve(jacobian,b1);
                y2 = linear.solve(jacobian,b2);
                return y1.allFinite() && y2.allFinite();
            } catch (const std::exception&) { return false; }
        };
    if(assembler->usesPopulationState()) {
        system.predictState=[self,assembler](const VectorXd& previous,const VectorXd& current,Real ratio) {
            (void)assembler->unpack(previous);(void)assembler->unpack(current);
            if(ratio==0.)return current;
            const double cutoff=self->cfg_.inputScaling.unitSystem().m3ToInternalConcentration(self->cfg_.populationPredictorLowDensity_cm3*1e6);
            auto next=predictPopulation(previous,current,ratio,cutoff);
            if(cutoff>0.) {
                int skipped=0;for(int k=current.size()/3;k<current.size();++k)if(current[k]>0.&&current[k]<=cutoff)++skipped;
                if(skipped)runtimeLogInfo("population_predictor_cutoff: skipped="+std::to_string(skipped));
            }
            // Constant populations are exact stationary predictor coordinates;
            // do not round them through exp(log(c)), especially at contacts.
            for(int k=current.size()/3;k<current.size();++k)
                if(previous[k]==current[k])next[k]=current[k];
            (void)assembler->unpack(next);return next;
        };
        // Only BVMethods can assemble its population tangents directly in
        // relative carrier columns. Legacy keeps its physical Jacobian and
        // the ordinary row/column equilibration for dilute states.
        if (assembler->usesBvMethodsOperator()) {
            system.jacobianInColumns=[self,assembler,biasesAt](const VectorXd& x,Real v,const VectorXd& columns){
                return assembler->populationJacobianInColumns(x,self->buildBoundaryConditions(*assembler,biasesAt(v)),columns);
            };
            system.bootstrapJacobianInColumns=[self,assembler,biasesAt](const VectorXd& x,Real v,const VectorXd& columns){
                return assembler->populationJacobianInColumns(x,self->buildBoundaryConditions(*assembler,biasesAt(v)),columns,true);
            };
            system.automaticDerivativeRecovery=cfg_.populationDerivativeRecovery=="automatic";
            system.postProtectionDerivativeRecovery=cfg_.populationDerivativeRecovery=="after_protection";
            system.protectedJacobian=[self,assembler,biasesAt](const VectorXd& x,Real v) {
                return assembler->populationProtectedJacobian(x,self->buildBoundaryConditions(*assembler,biasesAt(v)));
            };
            system.protectedJacobianInColumns=[self,assembler,biasesAt](const VectorXd& x,Real v,const VectorXd& columns) {
                return assembler->populationProtectedJacobian(x,self->buildBoundaryConditions(*assembler,biasesAt(v)),&columns);
            };
            system.bootstrapJacobian=[self,assembler,biasesAt](const VectorXd& x,Real v) {
                return assembler->populationBootstrapJacobian(x,self->buildBoundaryConditions(*assembler,biasesAt(v)));
            };
        }
        const double vt=thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
        system.columnScales=[assembler,nodeCount,potentialScale,vt](const ArclengthState& state) {
            (void)assembler->unpack(state.x);VectorXd scale=VectorXd::Ones(3*nodeCount+1);
            scale.head(nodeCount)=state.x.head(nodeCount).cwiseAbs().array()+vt/potentialScale;
            for(int i=0;i<nodeCount;++i)if(assembler->intrinsicDensity()[i]>0.) {
                scale[nodeCount+i]=state.x[nodeCount+i];scale[2*nodeCount+i]=state.x[2*nodeCount+i];
            }
            scale[3*nodeCount]=std::abs(state.lambda)+vt;return scale;
        };
        system.prescribedRows=[self,assembler,biasesAt,nodeCount](const VectorXd&,Real bias) {
            const auto b=self->buildBoundaryConditions(*assembler,biasesAt(bias));std::vector<int> rows;
            for(const auto& [i,v]:b.psi)rows.push_back(i);
            for(int i=0;i<nodeCount;++i) {
                if(assembler->intrinsicDensity()[i]<=0. || b.phin.contains(i))rows.push_back(nodeCount+i);
                if(assembler->intrinsicDensity()[i]<=0. || b.phip.contains(i))rows.push_back(2*nodeCount+i);
            }
            return rows;
        };
        system.makeTrial=[self,assembler,biasesAt,nodeCount,potentialScale,vt](const ArclengthState& before,const VectorXd& raw,Real dv,Real alpha) {
            const auto b=self->buildBoundaryConditions(*assembler,biasesAt(before.lambda));
            std::vector<bool> en(nodeCount),hp(nodeCount);
            for(int i=0;i<nodeCount;++i) {
                en[i]=assembler->intrinsicDensity()[i]>0. && !b.phin.contains(i);
                hp[i]=assembler->intrinsicDensity()[i]>0. && !b.phip.contains(i);
            }
            auto control=self->cfg_.populationUpdate;
            if(self->cfg_.populationMinimumDensity_cm3>0.)control.minimumDensity=self->cfg_.inputScaling.unitSystem().m3ToInternalConcentration(self->cfg_.populationMinimumDensity_cm3*1e6);
            auto c=makePopulationTrial(before.x,raw,alpha,en,hp,control);
            if(!c.flooredElectrons.empty()||!c.flooredHoles.empty())runtimeLogInfo("population_trial_floor: electrons="+
                std::to_string(c.flooredElectrons.size())+" holes="+std::to_string(c.flooredHoles.size()));
            if(self->cfg_.populationWriteback=="boltzmann_psi_compensated") {
                std::vector<bool> active(nodeCount);
                for(int i=0;i<nodeCount;++i)active[i]=assembler->intrinsicDensity()[i]>0.;
                // Contact populations also consume rounded psi. Their original
                // QF identities are checked by the unchanged residual/guards.
                c=compensateBoltzmannPsiWriteback(before.x,raw,std::move(c),potentialScale,vt,active);
            }
            (void)assembler->unpack(c.state);
            return ArclengthTrial{{c.state,before.lambda+c.alpha*dv},c.alpha,
                static_cast<int>(c.protectedElectrons.size()),static_cast<int>(c.protectedHoles.size())};
        };
        const auto matrix=system.jacobian;
        const auto rows=system.prescribedRows;
        system.solveJacobian=[matrix,rows](const VectorXd& x,Real v,const VectorXd& b,VectorXd& dx) {
            try{LinearSolver linear;dx=linear.solve(matrix(x,v),b,rows(x,v));return dx.allFinite();}
            catch(const std::exception&){return false;}
        };
        const auto evaluate=makeArclengthConvergenceEvaluator(activeContact);
        const auto block=cfg_.blockAbsoluteConvergence;
        system.convergenceSatisfied=[evaluate,block](const VectorXd& x,Real v) {
            const auto g=evaluate(x,v);
            const bool ceilings=block.mode!="enforce" || (g.blockNorms.psi<=block.psiResidualCeiling &&
                g.blockNorms.phin<=block.electronResidualCeiling && g.blockNorms.phip<=block.holeResidualCeiling);
            return ceilings && (!g.carrierRows.enforced || g.carrierRows.satisfied) &&
                (!g.globalContinuity.enforced || g.globalContinuity.satisfied);
        };
    }
    if (!assembler->usesPopulationState() && positiveDensityTrials) {
        const Real scaledThermalVoltage = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs) / potentialScale;
        system.trialState = [self, nodeCount, potentialScale, scaledThermalVoltage](
            const VectorXd& x, const VectorXd& deltaX, Real alpha) {
            if (x.size() != 3 * nodeCount || deltaX.size() != x.size() ||
                !x.allFinite() || !deltaX.allFinite() ||
                !std::isfinite(alpha) || alpha < 0.0) {
                throw std::invalid_argument("Invalid positive density trial state or step.");
            }
            VectorXd mapped = alpha * deltaX;
            for (int i = 0; i < nodeCount; ++i) {
                // Express the same Newton tangent in n,p coordinates. The
                // finite trial retains at least 10% of each old density;
                // log1p maps it back to quasi-Fermi potentials without the
                // enormous voltage steps caused by tiny minority densities.
                const Real relativeN = alpha * (deltaX(i) - deltaX(nodeCount + i)) /
                    scaledThermalVoltage;
                const Real relativeP = alpha * (deltaX(2 * nodeCount + i) - deltaX(i)) /
                    scaledThermalVoltage;
                mapped(nodeCount + i) = mapped(i) - scaledThermalVoltage *
                    std::log1p(std::max(-0.9, relativeN));
                mapped(2 * nodeCount + i) = mapped(i) + scaledThermalVoltage *
                    std::log1p(std::max(-0.9, relativeP));
            }
            applyConfiguredQuasiFermiStepCaps(
                mapped, self->cfg_, nodeCount, potentialScale, self->doping_);
            return VectorXd(x + mapped);
        };
    }
    if (!assembler->usesPopulationState() && !positiveDensityTrials && (cfg_.quasiFermiUpdateLimit_V > 0.0 ||
        cfg_.quasiFermiUpdateLimitMinority_V > 0.0)) {
        system.limitUpdate = [self, nodeCount, potentialScale](
            const VectorXd&, VectorXd& deltaX, Real& deltaLambda) {
            // A component-wise quasi-Fermi clip changes the bordered Newton
            // direction and can turn it into a non-descent direction for the
            // arclength constraint.  Uniformly scale the complete bordered
            // update so J dx + F_lambda d_lambda and the constraint row retain
            // the same first-order decrease.
            const Real scale = configuredQuasiFermiUniformStepScale(
                deltaX, self->cfg_, nodeCount, potentialScale, self->doping_);
            deltaX *= scale;
            deltaLambda *= scale;
        };
    }
    if(cfg_.populationStepFraction>0.) {
        const auto previousLimit=system.limitUpdate;
        const auto constraints=buildBoundaryConditions(*assembler,biasesAt(0.));
        const double fraction=cfg_.populationStepFraction;
        system.limitUpdate=[assembler,constraints,previousLimit,fraction](
            const VectorXd& x,VectorXd& dx,Real& dv) {
            if(previousLimit)previousLimit(x,dx,dv);
            const auto f=population_step::frame(*assembler,x,dx,constraints);
            const auto bound=population_step::relativeFractionToBoundary(f.relativeN,f.relativeP,f.active,fraction);
            dx*=bound.alpha;dv*=bound.alpha;
        };
        system.trialState=[assembler,constraints](const VectorXd& x,const VectorXd& dx,Real alpha) {
            const auto f=population_step::frame(*assembler,x,dx,constraints);
            return VectorXd(x+population_step::finiteStep(*assembler,f,dx,alpha));
        };
    }
    return system;
}

PrimaryVariableScales NewtonSolver::primaryVariableScales() const {
    const Real epsRef=constants::eps0*maxRelativePermittivityAcrossRegions(mesh_,matdb_,cfg_.temperature_K);
    const Real ni=maxIntrinsicDensityAcrossRegions(mesh_,matdb_,cfg_.temperature_K);
    const auto inputs=UnitScalingSystem::autoInputsFrom(mesh_,doping_,matdb_,ni);
    const auto& units=cfg_.inputScaling.unitSystem();
    const auto scale=UnitScalingSystem::fromInputs(cfg_.temperature_K,epsRef,inputs,cfg_.unitScalingRefs,units);
    const int n=mesh_.numNodes();PrimaryVariableScales out;
    out.device=VectorXd::Constant(3*n,scale.C0());
    out.device.head(n).setConstant(cfg_.inputScaling.isUnitScaling()?1.:scale.V0());
    out.voltage=scale.V0();
    // Match the qualified equation-normalization current reference. This is
    // a numerical port scale, not a separately inferred material current.
    out.current=units.currentDensityAM2PerInternal()*units.lengthMPerInternal()*1e-6*
        scale.C0()*scale.D0()*constants::q;
    out.validate(3*n);return out;
}

ApproximateNativeRhsScales NewtonSolver::approximateNativeRhsScales() const {
    // Frozen 300 K calibration, not a general model of native residual units.
    // Unsupported normalizations retain the valid physical E branch only.
    if(cfg_.coupledOperator!="bvmethods_tri3_qf" || cfg_.stateVariables!="potential_density" ||
       !cfg_.inputScaling.isUnitScaling() || cfg_.temperature_K!=300. ||
       constants::profileName!="sentaurus_u2022_12" || mesh_.numNodes()!=2600 || mesh_.numCells()!=5020)
        return {};
    const auto s=primaryVariableScales();
    const auto near=[](Real x,Real y){return std::abs(x/y-1.)<1e-10;};
    if(!near(s.device[mesh_.numNodes()],5.089999999999982e20) ||
       !near(s.current,2.98741982845799e-5) ||
       !near(maxRelativePermittivityAcrossRegions(mesh_,matdb_,cfg_.temperature_K),11.7))return {};
    return {11.7/11.9,1.2577838410557051e8,4.210268101838678e5};
}

VectorXd NewtonSolver::packArclengthState(
    const DDSolution& state, const DDSolution* referenceState) const
{
    auto assembler = makeArclengthAssembler(referenceState);
    return packReferencedSolution(
        *assembler, state, CoupledDDBoundaryConditions{}, true);
}

DDSolution NewtonSolver::unpackArclengthState(
    const VectorXd& x, const DDSolution* referenceState) const
{
    auto assembler = makeArclengthAssembler(referenceState);
    return makeSolution(*assembler, x, 0);
}

std::function<DDSolution(const VectorXd&)>
NewtonSolver::makeArclengthStateUnpacker(const DDSolution* referenceState) const
{
    auto assembler = makeArclengthAssembler(referenceState);
    const NewtonSolver* self = this;
    return [self, assembler](const VectorXd& x) {
        return self->makeSolution(*assembler, x, 0);
    };
}

std::function<NewtonStateConvergenceEvaluation(const VectorXd&, Real)>
NewtonSolver::makeArclengthConvergenceEvaluator(const std::string& activeContact) const
{
    // Retain exactly the packed coordinates, including split QF tails. Repacking
    // via ordinary voltage vectors would erase small carrier/current signals.
    auto assembler = makeArclengthAssembler();
    const NewtonSolver* self = this;
    return [self, assembler, activeContact](const VectorXd& x, Real bias) {
        auto biases = self->contactBiases_;
        biases[activeContact] = bias;
        const auto bcs = self->buildBoundaryConditions(*assembler, biases);
        NewtonStateConvergenceEvaluation evaluation;
        std::vector<CoupledDDCarrierTermDiagnostic> constrained,unconstrained;
        if(assembler->usesBvMethodsOperator()) {
            auto physical=assembler->bvMethodsPhysicalDiagnostics(x,bcs);
            evaluation.blockNorms=blockResidualInfo(physical.residual,self->mesh_.numNodes());
            constrained=std::move(physical.constrainedTerms);unconstrained=std::move(physical.unconstrainedTerms);
        } else {
            evaluation.blockNorms=blockResidualInfo(assembler->residual(x,bcs),self->mesh_.numNodes());
            constrained=assembler->carrierContinuityEquationTermDiagnostics(x,bcs);
            unconstrained=assembler->carrierContinuityEquationTermDiagnostics(x,{});
        }
        evaluation.carrierRows = evaluateCarrierRowConvergence(constrained,self->cfg_.carrierRowConvergence);
        std::vector<Index> electronContacts, holeContacts;
        for (const auto& [node, value] : bcs.phin) electronContacts.push_back(node);
        for (const auto& [node, value] : bcs.phip) holeContacts.push_back(node);
        evaluation.globalContinuity = evaluateGlobalContinuityClosure(
            unconstrained,
            electronContacts, holeContacts, self->cfg_.globalContinuityClosure);
        return evaluation;
    };
}

DensityUpdateErrorEvaluator NewtonSolver::makeArclengthUpdateErrorEvaluator(
    const std::string& activeContact,
    const DensityUpdateConvergenceConfig& config,
    bool includeContactVoltage) const
{
    validateDensityUpdateConvergence(config);
    if (cfg_.splitDDState || cfg_.electronQuantumPotential.enabled)
        throw std::invalid_argument("density update stopping requires classical isothermal DD");
    auto assembler = makeArclengthAssembler();
    std::size_t aliases = 0;
    bool foundContact = false;
    for (const auto& contact : mesh_.contacts()) {
        if (contact.name == activeContact) {
            aliases = contact.node_ids.size();
            foundContact = true;
        }
    }
    if (!foundContact || aliases == 0)
        throw std::invalid_argument("density update stopping requires a named nonempty contact");
    if (!includeContactVoltage)
        aliases = 0;
    const Real potentialScale = assembler->usesScaledState() ? assembler->potentialScale() : 1.;
    const Real densityToCm3 = cfg_.inputScaling.unitSystem().concentrationM3PerInternal()*1e-6;
    return [this, assembler, activeContact, config, aliases, potentialScale, densityToCm3](
        const VectorXd& before, const VectorXd& direction, const VectorXd& after,
        Real voltageAfter, Real rawVoltageDirection) {
        DensityUpdateError result;
        const int N = assembler->numNodes();
        if (before.size() != 3*N || direction.size() != 3*N || after.size() != 3*N ||
            !before.allFinite() || !direction.allFinite() || !after.allFinite()) return result;
        // Some contact laws change the free-carrier mask with applied bias.
        // Do not freeze a zero-bias mask for an entire continuation sweep.
        auto biases = contactBiases_;
        biases[activeContact] = voltageAfter;
        const auto bcs = buildBoundaryConditions(*assembler, biases);
        std::vector<bool> electronFree(N, true), holeFree(N, true);
        for (const auto& [node, value] : bcs.phin) electronFree[node] = false;
        for (const auto& [node, value] : bcs.phip) holeFree[node] = false;
        const auto [dn, dp] = assembler->carrierDensityTangent(before, direction);
        const auto n = assembler->electronDensity(after), p = assembler->holeDensity(after);
        DensityUpdateAccumulator<ArithmeticReal, Real> metric(config.digits);
        bool valid = true;
        for (int i=0; i<N; ++i) {
            metric.add(0, after(i)*potentialScale, direction(i)*potentialScale,
                       config.potentialReference_V);
            if (assembler->intrinsicDensity()[i] > 0. && electronFree[i]) {
                valid = valid && n(i)>0.;
                metric.add(1, n(i)*densityToCm3, dn(i)*densityToCm3, config.electronReference_cm3);
            }
            if (assembler->intrinsicDensity()[i] > 0. && holeFree[i]) {
                valid = valid && p(i)>0.;
                metric.add(2, p(i)*densityToCm3, dp(i)*densityToCm3, config.holeReference_cm3);
            }
        }
        for (std::size_t i=0; i<aliases; ++i)
            metric.add(3, voltageAfter, rawVoltageDirection, config.contactReference_V);
        return metric.finish(valid);
    };
}

ArclengthScalarFunctional NewtonSolver::makeArclengthContactCurrentFunctional(
    const std::string& contactName,
    Real currentScale,
    const DDSolution* referenceState) const
{
    if(cfg_.splitDDState)throw std::invalid_argument("Split terminal functional requires an explicit state");
    if (!std::isfinite(currentScale)) {
        throw std::invalid_argument(
            "NewtonSolver::makeArclengthContactCurrentFunctional: currentScale "
            "must be finite.");
    }
    const Contact* contact = nullptr;
    for (const Contact& candidate : mesh_.contacts()) {
        if (candidate.name == contactName) {
            contact = &candidate;
            break;
        }
    }
    if (contact == nullptr) {
        throw std::invalid_argument(
            "NewtonSolver::makeArclengthContactCurrentFunctional: unknown contact '" +
            contactName + "'.");
    }

    auto assembler = makeArclengthAssembler(referenceState);
    const int nodeCount = static_cast<int>(assembler->numNodes());
    VectorXd residualWeights = VectorXd::Zero(3 * nodeCount);
    for (Index node : contact->node_ids) {
        if (node >= static_cast<Index>(nodeCount))
            continue;
        const int row = static_cast<int>(node);
        // ContactCurrent::computeFromResidual uses
        // I=I_n-I_p=-q*r_n*scale+q*r_p*scale.
        residualWeights(nodeCount + row) -= constants::q;
        residualWeights(2 * nodeCount + row) += constants::q;
    }
    const Real continuityScale = assembler->continuityResidualScale();
    const Real scale = currentScale * continuityScale;

    ArclengthScalarFunctional functional;
    functional.value = [assembler, residualWeights, scale](
                           const VectorXd& x, Real) {
        const VectorXd residual =
            assembler->usesBvMethodsOperator()?assembler->bvMethodsTransportResidual(x,BvMethodsTransportPurpose::terminalValue):assembler->residual(x, CoupledDDBoundaryConditions{});
        return scale * residualWeights.dot(residual);
    };
    functional.linearize = [assembler, residualWeights, scale](
                               const VectorXd& x, Real) {
        const CoupledDDBoundaryConditions noBoundaryConditions;
        const VectorXd residual = assembler->usesBvMethodsOperator()?assembler->bvMethodsTransportResidual(x,BvMethodsTransportPurpose::terminalLinearization):assembler->residual(x, noBoundaryConditions);
        const SparseMatrixd jacobian =
            assembler->usesBvMethodsOperator()?assembler->bvMethodsTransportJacobian(x):assembler->assembleJacobian(x, noBoundaryConditions);
        ArclengthScalarLinearization result;
        result.value = scale * residualWeights.dot(residual);
        Eigen::Index unrelatedNonfinite=0;
        result.stateDerivative = scale * weightedSparseRows(jacobian,residualWeights,&unrelatedNonfinite);
        if(unrelatedNonfinite)runtimeLogInfo("terminal_current_derivative: excluded_unrelated_nonfinite="+
            std::to_string(unrelatedNonfinite)+" selected_derivative_finite="+
            std::to_string(result.stateDerivative.allFinite()));
        result.parameterDerivative = 0.0;
        return result;
    };
    if(assembler->usesPopulationState()&&assembler->usesBvMethodsOperator()) {
        functional.linearizeInColumns=[assembler,residualWeights,scale](const VectorXd& x,Real,const VectorXd& columns){
            ArclengthScalarLinearization result;
            result.value=scale*residualWeights.dot(assembler->bvMethodsTransportResidual(x,BvMethodsTransportPurpose::terminalLinearization));
            result.stateDerivative=scale*weightedSparseRows(assembler->bvMethodsTransportJacobianInColumns(x,columns),residualWeights);
            result.parameterDerivative=0.;return result;
        };
    }
    return functional;
}

std::unordered_map<std::string, Real> NewtonSolver::evaluatePackedContactCurrents(
    const VectorXd& packed,Real currentScale) const {
    if(!std::isfinite(currentScale))throw std::invalid_argument("invalid current scale");
    auto assembler=makeArclengthAssembler();
    const int N=mesh_.numNodes();
    if(packed.size()!=3*N || !packed.allFinite())throw std::invalid_argument("invalid packed current state");
    const auto residual=assembler->residual(packed,{});
    const Real scale=currentScale*assembler->continuityResidualScale();
    std::unordered_map<std::string,Real> result;
    for(const auto& contact:mesh_.contacts()){
        VectorXd weights=VectorXd::Zero(3*N);
        for(auto node:contact.node_ids)if(node<static_cast<Index>(N)){
            weights(N+node)-=constants::q;weights(2*N+node)+=constants::q;
        }
        result[contact.name]=scale*weights.dot(residual);
    }
    return result;
}

Real NewtonSolver::evaluateStateTerminalCurrent(const DDSolution& state,
    const std::string& contactName,Real currentScale) const {
    if(!std::isfinite(currentScale))throw std::invalid_argument("Invalid terminal current scale");
    const Contact* contact=nullptr;
    for(const auto& c:mesh_.contacts())if(c.name==contactName)contact=&c;
    if(!contact)throw std::invalid_argument("Unknown terminal current contact");
    auto assembler=makeArclengthAssembler(cfg_.splitDDState?&state:nullptr);
    // This is read-only evaluation, not continuation. Restore the checkpoint's
    // exact QF coordinate frame instead of choosing a new geometric origin.
    if(state.packedLow.size()!=0)
        assembler->setQuasiFermiReferenceFields(state.electronQfReference,state.holeQfReference);
    if (state.populationLow.size()) {
        const auto terms = assembler->compensatedPopulationDiagnostics(
            state.psi,state.n,state.p,state.populationLow,{},false);
        vela::ArithmeticReal sum = 0.;
        for (Index node : contact->node_ids) {
            const auto& row = terms.unconstrainedTerms[node];
            sum += static_cast<vela::ArithmeticReal>(row.holeFlux) - row.electronFlux;
        }
        return static_cast<Real>(sum*constants::q*assembler->continuityResidualScale()*currentScale);
    }
    restoreElectronQuantumPotential(*assembler, state);
    const auto x=packReferencedSolution(*assembler,state,{});
    const auto residual=assembler->usesBvMethodsOperator()?assembler->bvMethodsTransportResidual(x,BvMethodsTransportPurpose::terminalValidation):assembler->residual(x,{});const int N=mesh_.numNodes();
    vela::ArithmeticReal sum=0;
    for(auto i:contact->node_ids)sum+=static_cast<vela::ArithmeticReal>(residual(2*N+i))-static_cast<vela::ArithmeticReal>(residual(N+i));
    return static_cast<Real>(sum*static_cast<vela::ArithmeticReal>(constants::q)*assembler->continuityResidualScale()*currentScale);
}

#include "ZeroLoadLineProbe.inl"

NewtonTerminalCurrentAdjointEvaluation
NewtonSolver::evaluateTerminalCurrentAdjoint(
    const DDSolution& state,
    const std::string& contactName,
    Real currentScale) const
{
    const auto biasIt = contactBiases_.find(contactName);
    if (biasIt == contactBiases_.end()) {
        throw std::invalid_argument(
            "NewtonSolver::evaluateTerminalCurrentAdjoint: contact '" +
            contactName + "' is not present in the contact bias map.");
    }

    const VectorXd x = packArclengthState(state);
    ArclengthSystem system = makeArclengthSystem(contactName);
    ArclengthScalarFunctional functional =
        makeArclengthContactCurrentFunctional(contactName, currentScale);
    const ArclengthScalarLinearization linearization =
        functional.linearize(x, biasIt->second);
    const SparseMatrixd jacobian = system.jacobian(x, biasIt->second);
    SparseMatrixd transpose = jacobian.transpose();
    transpose.makeCompressed();

    LinearSolver linearSolver;
    const VectorXd adjoint =
        linearSolver.solve(transpose, linearization.stateDerivative);
    const VectorXd residual =
        transpose * adjoint - linearization.stateDerivative;
    const Real derivativeNorm = linearization.stateDerivative.norm();

    NewtonTerminalCurrentAdjointEvaluation result;
    result.current = linearization.value;
    result.contactBias_V = biasIt->second;
    result.nodeCount = static_cast<int>(mesh_.numNodes());
    result.potentialScale_V = cfg_.inputScaling.isUnitScaling()
        ? buildScalingSpec().V0
        : 1.0;
    result.stateDerivativeNorm = derivativeNorm;
    result.adjointNorm = adjoint.norm();
    result.adjointResidualNorm = residual.norm();
    result.adjointRelativeResidual = result.adjointResidualNorm /
        std::max(derivativeNorm, std::numeric_limits<Real>::min());
    result.stateDerivative = linearization.stateDerivative;
    result.adjoint = adjoint;
    return result;
}

Real NewtonSolver::maxContactMajorityQuasiFermiDrop(const DDSolution& state) const
{
    return maxContactMajorityQuasiFermiDrop(state, {});
}

Real NewtonSolver::maxContactMajorityQuasiFermiDrop(
    const DDSolution& state,
    const std::vector<std::string>& contacts) const
{
    if (state.phin.size() < static_cast<int>(mesh_.numNodes()) ||
        state.phip.size() < static_cast<int>(mesh_.numNodes()) ||
        state.n.size() < static_cast<int>(mesh_.numNodes()) ||
        state.p.size() < static_cast<int>(mesh_.numNodes())) {
        return std::numeric_limits<Real>::infinity();
    }

    std::vector<bool> isAnyTransportContactNode(mesh_.numNodes(), false);
    for (const Contact& contact : mesh_.contacts()) {
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        for (Index node : contact.node_ids) {
            if (node < mesh_.numNodes())
                isAnyTransportContactNode[node] = true;
        }
    }

    Real maxDrop = 0.0;
    for (const Contact& contact : mesh_.contacts()) {
        if (!contacts.empty() &&
            std::find(contacts.begin(), contacts.end(), contact.name) ==
                contacts.end()) {
            continue;
        }
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        std::vector<bool> isContactNode(mesh_.numNodes(), false);
        for (Index node : contact.node_ids) {
            if (node >= mesh_.numNodes())
                continue;
            isContactNode[node] = true;
        }
        for (Index e = 0; e < mesh_.numEdges(); ++e) {
            const Edge& edge = mesh_.getEdge(e);
            const bool n0Contact = isContactNode[edge.n0];
            const bool n1Contact = isContactNode[edge.n1];
            if (n0Contact == n1Contact)
                continue;
            // A boundary edge connecting two different contacts has no
            // contact-to-interior transport drop and must not contaminate a
            // contact-scoped branch guard.
            if (isAnyTransportContactNode[edge.n0] &&
                isAnyTransportContactNode[edge.n1]) {
                continue;
            }
            // A single electrode can span both n+ and p+ semiconductor.
            // Select the equilibrium majority species at this edge's contact
            // endpoint, not from an average over unrelated contact segments.
            // Compensated/undoped endpoints retain the guard on both species.
            const Index contactNode = n0Contact ? edge.n0 : edge.n1;
            const Real contactNetDoping = doping_.netDoping(contactNode);
            const bool electronMajority = contactNetDoping > 0.0;
            const bool holeMajority = contactNetDoping < 0.0;
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            const bool electronTransportEdge =
                state.n(i) > 0.0 && state.n(j) > 0.0 &&
                std::isfinite(state.n(i)) && std::isfinite(state.n(j));
            const bool holeTransportEdge =
                state.p(i) > 0.0 && state.p(j) > 0.0 &&
                std::isfinite(state.p(i)) && std::isfinite(state.p(j));
            if ((electronMajority || !holeMajority) && electronTransportEdge) {
                const Real drop = state.hasReferencedElectronQuasiFermi()
                    ? static_cast<Real>(
                        (static_cast<vela::ArithmeticReal>(
                            state.electronQuasiFermiReferenceAt(j)) -
                         static_cast<vela::ArithmeticReal>(
                            state.electronQuasiFermiReferenceAt(i))) +
                        (static_cast<vela::ArithmeticReal>(state.phinIncrement(j)) -
                         static_cast<vela::ArithmeticReal>(state.phinIncrement(i))))
                    : state.phin(j) - state.phin(i);
                maxDrop = std::max(maxDrop, std::abs(drop));
            }
            if ((holeMajority || !electronMajority) && holeTransportEdge) {
                const Real drop = state.hasReferencedHoleQuasiFermi()
                    ? static_cast<Real>(
                        (static_cast<vela::ArithmeticReal>(
                            state.holeQuasiFermiReferenceAt(j)) -
                         static_cast<vela::ArithmeticReal>(
                            state.holeQuasiFermiReferenceAt(i))) +
                        (static_cast<vela::ArithmeticReal>(state.phipIncrement(j)) -
                         static_cast<vela::ArithmeticReal>(state.phipIncrement(i))))
                    : state.phip(j) - state.phip(i);
                maxDrop = std::max(maxDrop, std::abs(drop));
            }
        }
    }
    return maxDrop;
}
NewtonResidualEvaluation NewtonSolver::evaluateResidual(const DDSolution& state) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
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

NewtonPoissonTermEvaluation NewtonSolver::evaluatePoissonTerms(
    const DDSolution& state) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup,
            cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, mobilityConfig, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, buildScalingSpec(), cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics, cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);

    NewtonPoissonTermEvaluation evaluation;
    evaluation.rows = assembler.poissonTermDiagnostics(
        x, bcs, state.n, state.p);
    evaluation.scaledState = assembler.usesScaledState();
    evaluation.potentialScale = potentialScale;
    return evaluation;
}

nlohmann::json NewtonSolver::exportPoissonArithmeticModel(const DDSolution& state) const
{
    RecombinationModelConfig recombinationConfig=recombinationModelConfig(
        cfg_.recombination,cfg_.taun,cfg_.taup,cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials=cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy=cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn=cfg_.augerCn;
    recombinationConfig.augerCp=cfg_.augerCp;
    recombinationConfig.augerDensityDependence=cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration=cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct=cfg_.augerExcessProduct;
    recombinationConfig.bandToBand=cfg_.bandToBand;
    CoupledDDAssembler assembler(mesh_,matdb_,doping_,thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs),
        cfg_.mobility,recombinationConfig,cfg_.bandgapNarrowing,cfg_.impactIonization,
        fixedCharges_,sheetCharges_,buildScalingSpec(),cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics,cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler,state);
    configureAssemblerNumerics(assembler,&state);
    const auto bcs=buildBoundaryConditions(assembler);
    const auto x=packReferencedSolution(assembler,state,bcs);
    const auto model=assembler.poissonArithmeticModel(x,bcs);
    nlohmann::json result={{"schema","vela.poisson-arithmetic-model.v1"},
        {"read_only",true},{"qualified",false},{"carrier_quasi_fermi_fixed",true},
        {"density_internal_to_m3",cfg_.inputScaling.unitSystem().internalConcentrationToM3(1.)},
        {"nodes",nlohmann::json::array()},{"edges",nlohmann::json::array()},
        {"dirichlet",nlohmann::json::array()},{"zero_low_residual",nlohmann::json::array()},
        {"zero_low_residual_double",nlohmann::json::array()}};
    for(const auto& a:model.nodes)result["nodes"].push_back({a.psi,a.n,a.p,a.xpsi,a.xn,a.xp,a.potentialScale,
        a.eref,a.href,a.niElectron,a.niHole,a.vt,a.doping,a.voln,a.volp,a.vold,a.q,a.area,a.scale,a.interfaceRhs});
    for(const auto& e:model.edges)result["edges"].push_back({e.i,e.j,e.g});
    for(const auto& [i,value]:model.dirichlet)result["dirichlet"].push_back({i,value});
    const auto values=poisson_precision::evaluate(model);
    for(std::size_t i=0;i<values.size();++i) {
        result["zero_low_residual"].push_back(values[i].str(40,std::ios_base::scientific));
        // Preserve the actual binary64 projection separately: a rounded decimal
        // string can move an exact halfway value to the other side of a tie.
        const auto boundary=bcs.psi.find(static_cast<Index>(i));
        result["zero_low_residual_double"].push_back(boundary==bcs.psi.end()
            ? static_cast<Real>(values[i]) : x(static_cast<int>(i))-boundary->second);
    }
    return result;
}

NewtonLinearSystemEvaluation NewtonSolver::evaluateLinearSystem(const DDSolution& state,
    NewtonLinearSystemEvaluation* mobilityAlternative,int factorRefinementSteps) const
{
    return evaluateSourceDiagnosticPair(state,mobilityAlternative,factorRefinementSteps,0);
}

NewtonLinearSystemEvaluation NewtonSolver::evaluateSourceDiagnosticPair(const DDSolution& state,
    NewtonLinearSystemEvaluation* mobilityAlternative,int factorRefinementSteps,int sourceMode) const
{
    if(sourceMode<0 || sourceMode>2 || (sourceMode && !mobilityAlternative))
        throw std::invalid_argument("source linearization diagnostic requires an explicit paired output");
    if(factorRefinementSteps<0 || factorRefinementSteps>3)throw std::invalid_argument("invalid factor refinement bound");
    const auto t0=std::chrono::steady_clock::now();
    auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();};
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    const double prepared=elapsed();
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    // Match the production pre-cap solve, but never form or accept a trial.
    // Reject unsupported variants rather than silently audit another operator.
    if (cfg_.linearRefinementIterations != 0 || cfg_.continuityRowScaling.enabled ||
        cfg_.linearEquilibration.mode == "l2_row_column")
        throw std::invalid_argument("linear-system probe does not support refinement or row/equilibration scaling");
    const int N = static_cast<int>(mesh_.numNodes());
    NewtonLinearSystemEvaluation result;
    result.physicalJacobian = J;
    result.regularizedJacobian = J;
    addCarrierRowRegularization(result.regularizedJacobian, N, cfg_.carrierRegularizationScale, bcs);
    result.residual = raw;
    result.solveMatrix = result.regularizedJacobian;
    result.solveRhs = -raw;
    result.potentialScale = potentialScale;
    LinearSolver linearSolver;
    if (cfg_.globalContinuityClosure.linearSystemStabilization == "replace_free_row_sum") {
        std::vector<Index> electrons, holes;
        std::vector<bool> ec(N, false), hc(N, false);
        for (const auto& [node, value] : bcs.phin) { electrons.push_back(node); ec[node] = true; }
        for (const auto& [node, value] : bcs.phip) { holes.push_back(node); hc[node] = true; }
        for (int row : selectGlobalContinuityReplacementRows(result.regularizedJacobian, N, ec, hc))
            if (row >= 0) result.replacementRows.push_back(row);
        stabilizeGlobalContinuityLinearSystem(result.solveMatrix, result.solveRhs, N, electrons, holes);
        result.direction = solveGlobalContinuityStabilizedSystem(linearSolver,
            result.regularizedJacobian, -raw, N, electrons, holes);
    } else {
        result.direction = linearSolver.solve(result.regularizedJacobian, -raw);
    }
    if(factorRefinementSteps){
        if(cfg_.globalContinuityClosure.linearSystemStabilization!="off")throw std::invalid_argument("factor audit requires unstabilized system");
        const auto r=slot_population::refineWithFactors(result.solveMatrix,result.solveRhs,result.direction,potentialScale,N,factorRefinementSteps,
            [&](const VectorXd& rhs){return linearSolver.solve(result.solveMatrix,rhs);});
        result.timings["factor_refinement_steps"]=r.sweeps;
        result.timings["factor_refinement_before_mixed"]=r.before.mixedRatio;
        result.timings["factor_refinement_after_mixed"]=r.after.mixedRatio;
    }
    result.timings["shared_prepare_s"]=prepared;
    result.timings["original_assembly_solve_s"]=elapsed()-prepared;
    if(mobilityAlternative){
        const double alternateStart=elapsed();
        if(cfg_.diagnosticFreezeTransportQfMobility || cfg_.jacobian=="finite_difference" ||
           cfg_.carrierRegularizationScale!=0 || cfg_.globalContinuityClosure.linearSystemStabilization!="off")
            throw std::invalid_argument("shared mobility pair requires original analytic unstabilized system");
        assembler.setDiagnosticFreezeTransportQfMobility(true);
        NewtonLinearSystemEvaluation other;
        other.residual=assembler.residual(x,bcs);
        if(!(other.residual.array()==raw.array()).all())
            throw std::runtime_error("shared mobility pair changed RHS");
        other.physicalJacobian=assembler.assembleSourceDiagnosticJacobian(x,bcs,sourceMode);
        other.regularizedJacobian=other.solveMatrix=other.physicalJacobian;
        other.solveRhs=-other.residual;other.potentialScale=potentialScale;
        LinearSolver diagnosticSolver;
        other.direction=diagnosticSolver.solve(other.solveMatrix,other.solveRhs);
        if(factorRefinementSteps){
            const auto r=slot_population::refineWithFactors(other.solveMatrix,other.solveRhs,other.direction,potentialScale,N,factorRefinementSteps,
                [&](const VectorXd& rhs){return diagnosticSolver.solve(other.solveMatrix,rhs);});
            other.timings["factor_refinement_steps"]=r.sweeps;
            other.timings["factor_refinement_before_mixed"]=r.before.mixedRatio;
            other.timings["factor_refinement_after_mixed"]=r.after.mixedRatio;
        }
        *mobilityAlternative=std::move(other);
        result.timings["diagnostic_assembly_solve_s"]=elapsed()-alternateStart;
    }
    return result;
}

NewtonStepEvaluation NewtonSolver::evaluateStep(const DDSolution& state, Real directionalAuditScale) const
{
    if(cfg_.stateVariables=="potential_density")
        throw std::invalid_argument("evaluateStep QF direction fields are not population-coordinate diagnostics; use the population update report");
    if (!std::isfinite(directionalAuditScale) || directionalAuditScale < 0.)
        throw std::invalid_argument("directional audit scale must be finite and nonnegative");
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    const VectorXd raw = assembler.residual(x, bcs);
    SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    const int N = static_cast<int>(mesh_.numNodes());
    addCarrierRowRegularization(J, N, cfg_.carrierRegularizationScale, bcs);
    const VectorXd weights = continuityRowWeights(
        assembler, x, bcs, cfg_.continuityRowScaling);
    LinearSolver linearSolver;
    VectorXd step = solveNewtonLinearSystem(J, raw, weights, cfg_, linearSolver);
    if (cfg_.exactDirichletUpdates) {
        for (const auto& [node, value] : bcs.psi) step(node) = -raw(node);
        for (const auto& [node, value] : bcs.phin) step(N + node) = -raw(N + node);
        for (const auto& [node, value] : bcs.phip) step(2 * N + node) = -raw(2 * N + node);
    }
    const VectorXd rawStep = step;
    const Real rawStepNorm = step.norm();

    applyConfiguredStepCapsAndPoissonRecorrection(
        step, J, raw, cfg_, N, potentialScale, doping_);

    const VectorXd trialX = assembler.splitRuntime() ? assembler.splitRuntime()->candidate(x,step,1.,assembler.splitRuntime()->low()) : VectorXd(x + step);
    const VectorXd trialRaw = assembler.residual(trialX, bcs);

    NewtonStepEvaluation evaluation;
    const auto rawDefect = linear_refinement::highPrecisionDefect(J, raw, rawStep);
    const auto cappedDefect = linear_refinement::highPrecisionDefect(J, raw, step);
    evaluation.rawDeltaPsi = rawStep.head(N) * potentialScale;
    evaluation.rawDeltaPhin = rawStep.segment(N, N) * potentialScale;
    evaluation.rawDeltaPhip = rawStep.tail(N) * potentialScale;
    evaluation.rawLinearDefect = rawDefect.values;
    evaluation.cappedLinearDefect = cappedDefect.values;
    evaluation.rawLinearBackwardErrors = rawDefect.componentwiseBackwardErrors;
    evaluation.cappedLinearBackwardErrors = cappedDefect.componentwiseBackwardErrors;
    evaluation.directionalAuditScale = directionalAuditScale;
    evaluation.freeElectronRows = VectorXd::Zero(N);
    evaluation.freeHoleRows = VectorXd::Zero(N);
    for (const auto& terms : assembler.carrierContinuityEquationTermDiagnostics(x, bcs)) {
        const int node = static_cast<int>(terms.nodeId);
        evaluation.freeElectronRows(node) = terms.electronContinuityActive && !bcs.phin.contains(terms.nodeId);
        evaluation.freeHoleRows(node) = terms.holeContinuityActive && !bcs.phip.contains(terms.nodeId);
    }
    if (directionalAuditScale > 0.0) {
        const VectorXd direction = directionalAuditScale * step;
        evaluation.directionalAnalyticChange = J * direction;
        evaluation.directionalFiniteDifferenceChange =
            (0.5 * (assembler.residual((x + direction).eval(), bcs) -
                    assembler.residual((x - direction).eval(), bcs))).eval();
    } else {
        evaluation.directionalAnalyticChange = VectorXd::Zero(3 * N);
        evaluation.directionalFiniteDifferenceChange = VectorXd::Zero(3 * N);
    }
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);

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
    const DDSolution& replacementState,
    bool computeConditionEstimates) const
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);

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
            configureAssemblerNumerics(diagnostic, &state);
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
    if (computeConditionEstimates) {
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
    }
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
    const DDSolution& physicalPerturbation,
    bool freezeTransportMobility) const
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    // A perturbation is a tangent vector, not a physical state.  In
    // reference coordinates its quasi-Fermi components only scale; they must
    // not be translated by the state's quasi-Fermi reference.
    VectorXd dx(3 * N);
    dx.segment(0, N) = physicalPerturbation.psi / potentialScale;
    dx.segment(N, N) = physicalPerturbation.phin / potentialScale;
    dx.segment(2 * N, N) = physicalPerturbation.phip / potentialScale;
    if (dx.norm() == 0.0)
        throw std::invalid_argument(
            "NewtonSolver::evaluateDirectionalDerivative: perturbation must be non-zero.");

    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = assembler.assembleJacobian(x, bcs);
    const VectorXd analytic = J * dx;
    CoupledDDFeedbackStateSubstitution frozen;
    if (freezeTransportMobility) {
        if (cfg_.mobility.jacobianFieldDerivatives) {
            throw std::invalid_argument(
                "Frozen mobility JVP requires jacobian_field_derivatives=false.");
        }
        frozen.replaceTransportMobility = true;
        frozen.electronEdgeMobility = VectorXd::Zero(mesh_.numEdges());
        frozen.holeEdgeMobility = VectorXd::Zero(mesh_.numEdges());
        for (const auto& edge : assembler.sgEdgeFluxDiagnostics(x, bcs)) {
            frozen.electronEdgeMobility(static_cast<int>(edge.edgeId)) =
                edge.electronMobility_m2_V_s;
            frozen.holeEdgeMobility(static_cast<int>(edge.edgeId)) =
                edge.holeMobility_m2_V_s;
        }
    }
    VectorXd forward, backward, finiteDifference;
    if (auto runtime = assembler.splitRuntime()) {
        if (freezeTransportMobility)
            throw std::invalid_argument("Split DD frozen mobility JVP is unqualified");
        finiteDifference = runtime->symmetricDifference(x, dx, bcs, forward, backward);
    } else {
    forward = freezeTransportMobility
        ? assembler.feedbackSubstitutionResidual(x + dx, bcs, frozen)
        : assembler.residual(x + dx, bcs);
    backward = freezeTransportMobility
        ? assembler.feedbackSubstitutionResidual(x - dx, bcs, frozen)
        : assembler.residual(x - dx, bcs);
    finiteDifference = 0.5 * (forward - backward);
    }
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
    const std::string& mode,
    Real damping,
    const std::vector<int>& localCarrierRows) const
{
    if (!std::isfinite(damping) || damping <= 0. || damping > 1.)
        throw std::invalid_argument("Block diagnostic damping must be in (0,1]");
    if (mode != "carrier_local" && !localCarrierRows.empty())
        throw std::invalid_argument("Local rows require carrier_local mode");
    if (mode == "carrier_density" && damping != 1.)
        throw std::invalid_argument("Density diagnostic has its own backtracking");
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    const VectorXd raw = assembler.residual(x, bcs);
    const SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);

    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd step = VectorXd::Zero(3 * N);
    LinearSolver linearSolver;
    if (mode == "carrier_local") {
        // Restricted QF Newton correction. Densities are reconstructed by the
        // same assembler (including Fermi statistics), never by an exp update.
        if (localCarrierRows.empty() || localCarrierRows.size() > 64)
            throw std::invalid_argument("Local carrier diagnostic requires 1..64 DOFs");
        auto sorted = localCarrierRows;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument("Duplicate local carrier row");
        const auto terms = assembler.carrierContinuityTermDiagnostics(x,bcs);
        for (const int k : sorted)
            if (k<0 || k>=2*N || (k<N ? bcs.phin.count(k) : bcs.phip.count(k-N)) ||
                !(k<N ? terms[k].electronContinuityActive : terms[k-N].holeContinuityActive))
                throw std::invalid_argument("Local row must be an active noncontact carrier equation");
        const int M = static_cast<int>(sorted.size());
        Eigen::MatrixXd block(M,M);
        VectorXd rhs(M);
        for (int i=0;i<M;++i) {
            const int k=sorted[i];
            if (k<0 || k>=2*N || !(k<N ? terms[k].electronContinuityActive
                                                       : terms[k-N].holeContinuityActive))
                throw std::invalid_argument("Local row must be an active carrier equation");
            rhs[i]=-raw[N+k];
            for (int j=0;j<M;++j) block(i,j)=J.coeff(N+k,N+sorted[j]);
            const Real scale=block.row(i).cwiseAbs().maxCoeff();
            if (!(scale>0.) || !std::isfinite(scale))
                throw std::runtime_error("Invalid local Jacobian row scale");
            block.row(i)/=scale;rhs[i]/=scale;
        }
        const VectorXd delta=linearSolver.solve(block.sparseView(),rhs);
        if (!delta.allFinite() || (block*delta-rhs).norm()>1e-8*std::max(rhs.norm(),Real{1e-30}))
            throw std::runtime_error("Local carrier linear residual check failed");
        for(int i=0;i<M;++i) step[N+sorted[i]]=delta[i];
    } else if (mode == "poisson_only") {
        const SparseMatrixd block = sparseBlock(J, 0, 0, N, N);
        step.segment(0, N) = linearSolver.solve(block, -raw.segment(0, N));
    } else if (mode == "carrier_only" || mode == "carrier_density") {
        if (mode == "carrier_density" &&
            (usesFermiDirac(cfg_.carrierStatistics) || cfg_.electronQuantumPotential.enabled))
            throw std::invalid_argument("carrier_density candidate requires classical Boltzmann densities");
        if (mode == "carrier_density") {
            for (const auto& [node,value] : bcs.phin)
                if (std::abs(raw(N+node)) > 1e-12)
                    throw std::invalid_argument("carrier_density requires projected carrier contact values");
            for (const auto& [node,value] : bcs.phip)
                if (std::abs(raw(2*N+node)) > 1e-12)
                    throw std::invalid_argument("carrier_density requires projected carrier contact values");
        }
        SparseMatrixd block = sparseBlock(J, N, N, 2 * N, 2 * N);
        VectorXd coordinateFactor=VectorXd::Ones(2*N);
        if (mode == "carrier_density") {
            const VectorXd n=assembler.electronDensity(x),p=assembler.holeDensity(x);
            for(int node=0;node<N;++node) if(assembler.intrinsicDensity()[node]>0.) {
                if (!(n[node]>0.) || !(p[node]>0.))
                    throw std::runtime_error("Nonpositive base density in carrier_density");
                coordinateFactor[node]=-Vt/(potentialScale*n[node]);
                coordinateFactor[N+node]=Vt/(potentialScale*p[node]);
            }
            if (!coordinateFactor.allFinite())
                throw std::runtime_error("Nonfinite density coordinate transformation");
            for(int col=0;col<block.outerSize();++col)
                for(SparseMatrixd::InnerIterator entry(block,col);entry;++entry)
                    entry.valueRef()*=coordinateFactor[col];
        }
        const VectorXd solved=linearSolver.solve(block,-raw.segment(N,2*N));
        step.segment(N,2*N)=solved.cwiseProduct(coordinateFactor);
    } else {
        throw std::invalid_argument(
            "NewtonSolver::evaluateBlockStep: mode must be 'poisson_only' "
            "or 'carrier_only', 'carrier_local' or 'carrier_density'.");
    }

    const Real rawStepNorm = step.norm();
    if (mode == "carrier_density") {
        // The linear system above is solved in density coordinates. Freeze
        // physical row scales during backtracking so a larger trial flux cannot
        // artificially reduce the merit. Poisson remains frozen in this probe.
        const auto terms=assembler.carrierContinuityTermDiagnostics(x,bcs);
        VectorXd weights=VectorXd::Zero(2*N);
        for(int node=0;node<N;++node) {
            const auto& t=terms[node];
            if(t.electronContinuityActive) weights[node]=1./std::max({t.electronFluxAbsSum,
                std::abs(t.electronRecombination),std::abs(t.electronImpact),Real{1e-30}});
            if(t.holeContinuityActive) weights[N+node]=1./std::max({t.holeFluxAbsSum,
                std::abs(t.holeRecombination),std::abs(t.holeImpact),Real{1e-30}});
        }
        const auto merit=[&](const VectorXd& residual) {
            return residual.tail(2*N).cwiseProduct(weights).stableNorm();
        };
        const Real initialMerit=merit(raw);
        const VectorXd direction=step;
        bool accepted=false;
        nlohmann::json backtracks=nlohmann::json::array();
        for(int backtrack=0;backtrack<20;++backtrack) {
            const Real alpha=std::ldexp(Real{1.},-backtrack);
            step=alpha*direction;
            int depletionClips=0;
            for(int node=0;node<N;++node) if(assembler.intrinsicDensity()[node]>0.)
                for(int carrier=0;carrier<2;++carrier) {
                    const int row=(carrier+1)*N+node;
                    const Real sign=carrier==0 ? -1. : 1.;
                    const Real relativeChange=sign*step[row]*potentialScale/Vt;
                    if(!std::isfinite(relativeChange))
                        throw std::runtime_error("Nonfinite carrier density Newton candidate");
                    if(relativeChange < -.9) ++depletionClips;
                    // Limit depletion to 90 percent, then map back to physical QF.
                    step[row]=sign*Vt*std::log1p(std::max(relativeChange,Real{-0.9}))/potentialScale;
                }
            const VectorXd beforeCaps=step;
            applyConfiguredStepCaps(step,cfg_,N,potentialScale,doping_);
            int qfClips=0;
            for(int row=N;row<3*N;++row) if(step[row]!=beforeCaps[row]) ++qfClips;
            const Real trialMerit=merit(assembler.residual(x+step,bcs));
            backtracks.push_back({{"alpha",alpha},{"finite",std::isfinite(trialMerit)},
                {"trial_merit",std::isfinite(trialMerit) ? nlohmann::json(trialMerit) : nlohmann::json(nullptr)},
                {"depletion_clips",depletionClips},{"qf_clips",qfClips}});
            if(std::isfinite(trialMerit) && trialMerit<=initialMerit*(1.-1e-4*alpha)) {
                accepted=true;break;
            }
        }
        if(!accepted) {
            if(!std::isfinite(initialMerit))
                throw std::runtime_error("Nonfinite carrier density base merit");
            throw NewtonBlockNoCandidate(nlohmann::json{
                {"schema","vela.block_candidate_failure.v1"},
                {"failure_class","numerical_no_candidate"},{"mode",mode},
                {"initial_merit",initialMerit},{"backtracks",backtracks}}.dump());
        }
    }
    applyConfiguredStepCaps(step, cfg_, N, potentialScale, doping_);
    step *= damping;
    const VectorXd trialX = assembler.splitRuntime() ? assembler.splitRuntime()->candidate(x,step,1.,assembler.splitRuntime()->low()) : VectorXd(x + step);
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
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
        if ((row < N && bcs.phin.contains(row)) ||
            (row >= N && bcs.phip.contains(row - N)))
            continue;
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
    const VectorXd trialX = assembler.splitRuntime() ? assembler.splitRuntime()->candidate(x,step,1.,assembler.splitRuntime()->low()) : VectorXd(x + step);
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
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
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
    std::vector<Real> electronPsiColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhinColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhipColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronContactColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronFreeColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPsiContactColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPsiFreeColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhinContactColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhinFreeColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhipContactColumnAbs(static_cast<std::size_t>(N), 0.0);
    std::vector<Real> electronPhipFreeColumnAbs(static_cast<std::size_t>(N), 0.0);
    const std::vector<bool> contactNodes = detail::contactNodeMask(mesh_);
    for (int col = 0; col < J.outerSize(); ++col) {
        for (SparseMatrixd::InnerIterator it(J, col); it; ++it) {
            const int row = static_cast<int>(it.row());
            const Real value = it.value();
            if (row >= N && row < 2 * N) {
                const std::size_t node = static_cast<std::size_t>(row - N);
                const Real absValue = std::abs(value);
                electronRowAbs[node] += absValue;
                electronRowL2Sq[node] += value * value;
                const std::size_t columnNode =
                    static_cast<std::size_t>(col % N);
                const bool contactColumn = contactNodes[columnNode];
                if (col < N) {
                    electronPsiColumnAbs[node] += absValue;
                    (contactColumn ? electronPsiContactColumnAbs
                                   : electronPsiFreeColumnAbs)[node] += absValue;
                } else if (col < 2 * N) {
                    electronPhinColumnAbs[node] += absValue;
                    (contactColumn ? electronPhinContactColumnAbs
                                   : electronPhinFreeColumnAbs)[node] += absValue;
                } else {
                    electronPhipColumnAbs[node] += absValue;
                    (contactColumn ? electronPhipContactColumnAbs
                                   : electronPhipFreeColumnAbs)[node] += absValue;
                }
                if (contactColumn)
                    electronContactColumnAbs[node] += absValue;
                else
                    electronFreeColumnAbs[node] += absValue;
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
        row.electronPsiColumnAbsSum =
            electronPsiColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhinColumnAbsSum =
            electronPhinColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhipColumnAbsSum =
            electronPhipColumnAbs[static_cast<std::size_t>(i)];
        row.electronContactColumnAbsSum =
            electronContactColumnAbs[static_cast<std::size_t>(i)];
        row.electronFreeColumnAbsSum =
            electronFreeColumnAbs[static_cast<std::size_t>(i)];
        row.electronPsiContactColumnAbsSum =
            electronPsiContactColumnAbs[static_cast<std::size_t>(i)];
        row.electronPsiFreeColumnAbsSum =
            electronPsiFreeColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhinContactColumnAbsSum =
            electronPhinContactColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhinFreeColumnAbsSum =
            electronPhinFreeColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhipContactColumnAbsSum =
            electronPhipContactColumnAbs[static_cast<std::size_t>(i)];
        row.electronPhipFreeColumnAbsSum =
            electronPhipFreeColumnAbs[static_cast<std::size_t>(i)];
        row.rawDeltaPhin_V = rawStep(N + i) * potentialScale;
        row.rawDeltaPhip_V = rawStep(2 * N + i) * potentialScale;
        row.cappedDeltaPhin_V = cappedStep(N + i) * potentialScale;
        row.cappedDeltaPhip_V = cappedStep(2 * N + i) * potentialScale;
        evaluation.rows.push_back(row);
    }
    return evaluation;
}

NewtonPoissonLinearDiagnosticsEvaluation
NewtonSolver::evaluatePoissonLinearDiagnostics(
    const DDSolution& state,
    Index focusNode) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (focusNode < 0 || focusNode >= N)
        throw std::invalid_argument(
            "NewtonSolver::evaluatePoissonLinearDiagnostics: focus node is out of range.");

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    const DDScalingSpec scaling = buildScalingSpec();
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, mobilityConfig, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, scaling, cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics, cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    const VectorXd raw = assembler.residual(x, bcs);
    SparseMatrixd J = (cfg_.jacobian == "finite_difference")
        ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
        : assembler.assembleJacobian(x, bcs);
    addCarrierRowRegularization(J, N, cfg_.carrierRegularizationScale, bcs);

    const VectorXd rowWeights = continuityRowWeights(
        assembler, x, bcs, cfg_.continuityRowScaling);
    const SparseMatrixd productionMatrix = cfg_.continuityRowScaling.enabled
        ? leftScaleRows(J, rowWeights)
        : J;
    VectorXd productionRhs = -raw;
    if (cfg_.continuityRowScaling.enabled)
        productionRhs = -raw.cwiseProduct(rowWeights);

    LinearSolver rawSolver;
    const VectorXd rawStep = rawSolver.solve(productionMatrix, productionRhs);
    const TwoSidedL2EquilibratedSystem equilibrated =
        twoSidedL2Equilibrate(productionMatrix, productionRhs);
    LinearSolver equilibratedSolver;
    const VectorXd equilibratedCoordinates =
        equilibratedSolver.solve(equilibrated.matrix, equilibrated.rhs);
    const VectorXd equilibratedStep =
        equilibratedCoordinates.cwiseProduct(equilibrated.columnScale);

    const VectorXd rawClosure = J * rawStep + raw;
    const VectorXd equilibratedClosure = J * equilibratedStep + raw;
    const Real closureScale = std::max<Real>(1.0, raw.norm());
    const SparseMatrixd poissonBlock = sparseBlock(J, 0, 0, N, N);
    const VectorXd poissonRows = sparseRowL2Norms(poissonBlock);
    const VectorXd poissonColumns = sparseColumnL2Norms(poissonBlock);
    const VectorXd fullRows = sparseRowL2Norms(J);
    const VectorXd fullColumns = sparseColumnL2Norms(J);

    std::vector<Index> patchNodes;
    patchNodes.push_back(focusNode);
    for (SparseMatrixd::InnerIterator entry(poissonBlock, focusNode); entry; ++entry) {
        if (entry.row() != focusNode)
            patchNodes.push_back(static_cast<Index>(entry.row()));
    }
    for (int column = 0; column < poissonBlock.outerSize(); ++column) {
        if (column != focusNode && poissonBlock.coeff(focusNode, column) != 0.0)
            patchNodes.push_back(static_cast<Index>(column));
    }
    std::sort(patchNodes.begin(), patchNodes.end());
    patchNodes.erase(std::unique(patchNodes.begin(), patchNodes.end()), patchNodes.end());

    Eigen::MatrixXd patch = Eigen::MatrixXd::Zero(
        static_cast<int>(patchNodes.size()), static_cast<int>(patchNodes.size()));
    for (int row = 0; row < patch.rows(); ++row) {
        for (int column = 0; column < patch.cols(); ++column) {
            patch(row, column) = poissonBlock.coeff(
                patchNodes[static_cast<std::size_t>(row)],
                patchNodes[static_cast<std::size_t>(column)]);
        }
    }

    NewtonPoissonLinearDiagnosticsEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.focusNode = focusNode;
    evaluation.focusPatchNodes = patchNodes;
    evaluation.focusPatchRawCondition = matrixConditionEstimate(patch);
    evaluation.focusPatchEquilibratedCondition =
        matrixConditionEstimate(l2Equilibrated(patch));
    evaluation.potentialScale = potentialScale;
    evaluation.poissonRowNormSpread = positiveNormSpread(poissonRows);
    evaluation.poissonColumnNormSpread = positiveNormSpread(poissonColumns);
    evaluation.fullRowNormSpread = positiveNormSpread(fullRows);
    evaluation.fullColumnNormSpread = positiveNormSpread(fullColumns);
    evaluation.rawStepNorm = rawStep.norm();
    evaluation.equilibratedStepNorm = equilibratedStep.norm();
    evaluation.rawLinearClosureNorm = rawClosure.norm();
    evaluation.equilibratedLinearClosureNorm = equilibratedClosure.norm();
    evaluation.rawRelativeLinearClosure = rawClosure.norm() / closureScale;
    evaluation.equilibratedRelativeLinearClosure =
        equilibratedClosure.norm() / closureScale;
    evaluation.relativeStepDifference =
        (rawStep - equilibratedStep).norm() /
        std::max<Real>(1.0, std::max(rawStep.norm(), equilibratedStep.norm()));
    evaluation.rows.reserve(static_cast<std::size_t>(N));
    for (int node = 0; node < N; ++node) {
        NewtonPoissonLinearRowDiagnostic row;
        row.nodeId = static_cast<Index>(node);
        row.residual = raw(node);
        row.diagonal = J.coeff(node, node);
        row.poissonRowL2Norm = poissonRows(node);
        row.poissonColumnL2Norm = poissonColumns(node);
        row.fullRowL2Norm = fullRows(node);
        row.fullColumnL2Norm = fullColumns(node);
        row.rawDeltaPsi_V = rawStep(node) * potentialScale;
        row.equilibratedDeltaPsi_V = equilibratedStep(node) * potentialScale;
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    const DDScalingSpec scaling = buildScalingSpec();
    const auto makeRecombinationConfig =
        [&](const std::vector<std::string>& models) {
            RecombinationModelConfig config =
                recombinationModelConfig(
                    models, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
            config.activeMaterials = cfg_.recombinationMaterials;
            config.volumePolicy = cfg_.recombinationVolumePolicy;
            config.augerCn = cfg_.augerCn;
            config.augerCp = cfg_.augerCp;
            config.augerDensityDependence = cfg_.augerDensityDependence;
            config.augerWithGeneration = cfg_.augerWithGeneration;
            config.augerExcessProduct = cfg_.augerExcessProduct;
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
                scaling, cfg_.carrierDiagonalFloor, cfg_.carrierStatistics,
                cfg_.electronQuantumPotential);
            restoreElectronQuantumPotential(assembler, state);
            configureAssemblerNumerics(assembler, &state);
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
    const VectorXd x = packReferencedSolution(fullAssembler, state, bcs);
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

NewtonGlobalContinuityClosureEvaluation NewtonSolver::evaluateFrozenGlobalContinuityClosure(
    const DDSolution& state) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    // Freeze the saved coordinate representation as well as the physical state.
    // Choosing new contact basins can round away a sub-ULP saved increment.
    VectorXd electronReference(mesh_.numNodes()), holeReference(mesh_.numNodes());
    for (Index node = 0; node < mesh_.numNodes(); ++node) {
        electronReference(node) = state.hasReferencedElectronQuasiFermi()
            ? state.electronQuasiFermiReferenceAt(node)
            : assembler.electronQuasiFermiReferenceAt(node);
        holeReference(node) = state.hasReferencedHoleQuasiFermi()
            ? state.holeQuasiFermiReferenceAt(node)
            : assembler.holeQuasiFermiReferenceAt(node);
    }
    assembler.setQuasiFermiReferenceFields(electronReference, holeReference);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    // Empty BCs here preserve the saved contact reference/increment pairs.
    // Real BCs still define which physical rows contribute terminal flux.
    const VectorXd x = packReferencedSolution(assembler, state, {});
    std::vector<Index> ec, hc;
    for (const auto& [node, value] : bcs.phin) ec.push_back(node);
    for (const auto& [node, value] : bcs.phip) hc.push_back(node);
    return evaluateGlobalContinuityClosure(
        assembler.carrierContinuityEquationTermDiagnostics(x, {}),
        ec, hc, cfg_.globalContinuityClosure);
}

NewtonCarrierTermDiagnosticsEvaluation NewtonSolver::evaluateCarrierTermDiagnostics(
    const DDSolution& state,
    bool solvedEquationTerms) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    if (state.populationLow.size()) {
        auto physical = assembler.compensatedPopulationDiagnostics(
            state.psi,state.n,state.p,state.populationLow,bcs);
        NewtonCarrierTermDiagnosticsEvaluation evaluation;
        evaluation.residual.raw = std::move(physical.residual);
        evaluation.residual.blockNorms = blockResidualInfo(evaluation.residual.raw,mesh_.numNodes());
        evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
        evaluation.residual.scaledState = assembler.usesScaledState();
        evaluation.residual.potentialScale = potentialScale;
        evaluation.rows = solvedEquationTerms ? std::move(physical.constrainedTerms) :
            std::move(physical.unconstrainedTerms);
        return evaluation;
    }
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    const VectorXd raw = assembler.residual(x, bcs);

    NewtonCarrierTermDiagnosticsEvaluation evaluation;
    evaluation.residual.raw = raw;
    evaluation.residual.blockNorms = blockResidualInfo(raw, mesh_.numNodes());
    evaluation.residual.intrinsicDensity = assembler.intrinsicDensity();
    evaluation.residual.scaledState = assembler.usesScaledState();
    evaluation.residual.potentialScale = potentialScale;
    evaluation.rows = solvedEquationTerms
        ? assembler.carrierContinuityEquationTermDiagnostics(x, bcs)
        : assembler.carrierContinuityTermDiagnostics(x, bcs);
    evaluation.continuityResidualToAPerUm = constants::q
        * assembler.continuityResidualScale()
        * (scaling.enabled ? scaling.currentDensityLineIntegralFactor : 1.0)
        * 1.0e-6;
    return evaluation;
}

std::vector<NewtonJacobianBlockAuditRow> NewtonSolver::evaluateJacobianBlockAudit(
    const DDSolution& state,
    Real finiteDifferenceStep,
    std::vector<std::string> blocks,
    const std::string& finiteDifferenceMode,
    std::vector<Index> columnNodes) const
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

    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    const DDScalingSpec scaling = buildScalingSpec();

    const auto makeRecombinationConfig =
        [&](const std::vector<std::string>& models) {
            RecombinationModelConfig config =
                recombinationModelConfig(
                    models, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
            config.activeMaterials = cfg_.recombinationMaterials;
            config.volumePolicy = cfg_.recombinationVolumePolicy;
            config.augerCn = cfg_.augerCn;
            config.augerCp = cfg_.augerCp;
            config.augerDensityDependence = cfg_.augerDensityDependence;
            config.augerWithGeneration = cfg_.augerWithGeneration;
            config.augerExcessProduct = cfg_.augerExcessProduct;
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
                cfg_.carrierStatistics,
                cfg_.electronQuantumPotential);
            restoreElectronQuantumPotential(assembler, state);
            configureAssemblerNumerics(assembler, &state);
            return assembler;
        };

    CoupledDDAssembler baseAssembler =
        makeAssembler(noRecombinationConfig, noImpactConfig);
    const CoupledDDBoundaryConditions bcs =
        buildBoundaryConditions(baseAssembler);
    const Real potentialScale =
        baseAssembler.usesScaledState() ? baseAssembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(baseAssembler, state, bcs);
    std::vector<int> selectedColumns;
    if (!columnNodes.empty()) {
        std::sort(columnNodes.begin(), columnNodes.end());
        columnNodes.erase(std::unique(columnNodes.begin(), columnNodes.end()), columnNodes.end());
        selectedColumns.reserve(3 * columnNodes.size());
        for (Index node : columnNodes) {
            if (node >= mesh_.numNodes()) {
                throw std::out_of_range(
                    "NewtonSolver::evaluateJacobianBlockAudit: column node is outside the mesh.");
            }
            for (int block = 0; block < 3; ++block)
                selectedColumns.push_back(block * N + static_cast<int>(node));
        }
    }
    const auto filterSelectedColumns = [&](const SparseMatrixd& matrix) {
        if (selectedColumns.empty())
            return matrix;
        std::vector<Eigen::Triplet<double>> triplets;
        for (int col : selectedColumns) {
            for (SparseMatrixd::InnerIterator entry(matrix, col); entry; ++entry)
                triplets.emplace_back(entry.row(), entry.col(), entry.value());
        }
        SparseMatrixd filtered(matrix.rows(), matrix.cols());
        filtered.setFromTriplets(triplets.begin(), triplets.end());
        return filtered;
    };

    const auto matrixPair =
        [&](CoupledDDAssembler& assembler) {
            if (!selectedColumns.empty()) {
                std::vector<Eigen::Triplet<double>> triplets;
                for (int col : selectedColumns) {
                    const Real step = finiteDifferenceStep * std::max<Real>(1., std::abs(x(col)));
                    VectorXd plus=x, minus=x;
                    plus(col)+=step;minus(col)-=step;
                    const VectorXd rp=assembler.residual(plus,bcs);
                    const VectorXd rm=assembler.residual(minus,bcs);
                    const VectorXd derivative=(rp-rm)/(2.*step);
                    for (int row=0;row<x.size();++row)
                        if (derivative(row)!=0.) triplets.emplace_back(row,col,derivative(row));
                }
                SparseMatrixd reference(x.size(),x.size());
                reference.setFromTriplets(triplets.begin(),triplets.end());
                return std::pair<SparseMatrixd,SparseMatrixd>{
                    filterSelectedColumns(assembler.assembleJacobian(x,bcs)),reference};
            }
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
                detail::usesElementLocalAvalancheSource(
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
                    filterSelectedColumns(
                        withTerm.impactIonizationSourceJacobian(x, bcs)),
                    filterSelectedColumns(reference)};
            }
            SparseMatrixd analytic = termName == "srh_auger"
                ? withTerm.recombinationSourceJacobian(x,bcs)
                : SparseMatrixd(withTerm.assembleJacobian(x, bcs) -
                                withoutTerm.assembleJacobian(x, bcs));
            analytic = filterSelectedColumns(analytic);
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
                    std::vector<int> allColumns;
                    const std::vector<int>* columns = &selectedColumns;
                    if (selectedColumns.empty()) {
                        allColumns.resize(static_cast<std::size_t>(unknowns));
                        std::iota(allColumns.begin(), allColumns.end(), 0);
                        columns = &allColumns;
                    }
                    for (int col : *columns) {
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
            if (termName != "srh_auger")
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
        wants("poisson") || wants("transport") || wants("dirichlet_or_gauge");
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
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
    return assembler.sgEdgeFluxDiagnostics(x, bcs);
}

std::vector<CoupledDDElectronTransportFactorDiagnostic>
NewtonSolver::evaluateElectronTransportFactorDiagnostics(
    const DDSolution& state,
    const DDSolution& replacementState) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
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
        cfg_.carrierStatistics, cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const VectorXd baseline = packReferencedSolution(assembler, state, bcs);
    const VectorXd replacement =
        packReferencedSolution(assembler, replacementState, bcs);
    std::vector<CoupledDDElectronTransportFactorDiagnostic> rows =
        assembler.electronTransportFactorDiagnostics(
            baseline, replacement, bcs);
    // Freeze the Shapley endpoints to independently evaluated production
    // equation terms. The fourteen mixed interventions are diagnostic states;
    // evaluating both endpoints on fresh production paths makes their role as
    // authoritative vertices explicit and independent of intervention order.
    const auto baselineTerms = evaluateCarrierTermDiagnostics(state, true).rows;
    const auto replacementTerms =
        evaluateCarrierTermDiagnostics(replacementState, true).rows;
    for (auto& row : rows) {
        if (row.mask == 0)
            row.electronFlux = baselineTerms[row.nodeId].electronFlux;
        else if (row.mask == 15)
            row.electronFlux = replacementTerms[row.nodeId].electronFlux;
    }
    return rows;
}

CoupledDDElectronTransportSecantFactorEvaluation
NewtonSolver::evaluateElectronTransportSecantFactorDiagnostics(
    const DDSolution& state,
    const DDSolution& replacementState) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
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
        cfg_.carrierStatistics, cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const VectorXd baseline = packReferencedSolution(assembler, state, bcs);
    const VectorXd replacement =
        packReferencedSolution(assembler, replacementState, bcs);
    CoupledDDElectronTransportSecantFactorEvaluation evaluation =
        assembler.electronTransportSecantFactorDiagnostics(
            baseline, replacement, bcs);
    // Keep both Shapley endpoints authoritative and independent of the mixed
    // secant-factor arithmetic.
    const auto baselineTerms = evaluateCarrierTermDiagnostics(state, true).rows;
    const auto replacementTerms =
        evaluateCarrierTermDiagnostics(replacementState, true).rows;
    for (auto& row : evaluation.nodes) {
        if (row.mask == 0)
            row.electronFlux = baselineTerms[row.nodeId].electronFlux;
        else if (row.mask == 7)
            row.electronFlux = replacementTerms[row.nodeId].electronFlux;
    }
    return evaluation;
}

std::vector<CoupledDDTransportEdgeJacobianDiagnostic>
NewtonSolver::evaluateTransportEdgeJacobianDiagnostics(
    const DDSolution& state,
    Real physicalFiniteDifferenceStep_V) const
{
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, cfg_.mobility, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, buildScalingSpec(), cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics, cfg_.electronQuantumPotential);
    restoreElectronQuantumPotential(assembler, state);
    configureAssemblerNumerics(assembler, &state);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    const Real potentialScale =
        assembler.usesScaledState() ? assembler.potentialScale() : 1.0;
    const VectorXd x = packReferencedSolution(assembler, state, bcs);
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
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    configureAssemblerNumerics(assembler);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    return solve(buildInitialGuess(assembler, bcs));
}

NewtonPoissonBlockInitialization NewtonSolver::buildPoissonBlockInitialization() const
{
    if(cfg_.stateVariables=="potential_density") {
        auto c=cfg_;c.compensatedPopulationState=false;c.stateVariables="potential_quasi_fermi";c.populationDerivativeRecovery="full";c.primaryVariableCoordinates="physical";c.populationStoppingRule="legacy";c.populationSgCurrent="quasi_fermi";c.populationSgBootstrapLinearization="unchanged";
        c.electronQuantumPotential.couplingMode = "outer";
        c.populationLinearDefect = "mp100"; // This temporary solver only initializes Poisson.
        return NewtonSolver(mesh_,matdb_,doping_,contactBiases_,c,fixedCharges_,sheetCharges_,contactSpecs_).buildPoissonBlockInitialization();
    }
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    configureAssemblerNumerics(assembler);
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
    if(cfg_.stateVariables=="potential_density") {
        auto c=cfg_;c.compensatedPopulationState=false;c.stateVariables="potential_quasi_fermi";c.populationDerivativeRecovery="full";c.primaryVariableCoordinates="physical";c.populationStoppingRule="legacy";c.populationSgCurrent="quasi_fermi";c.populationSgBootstrapLinearization="unchanged";
        c.electronQuantumPotential.couplingMode = "outer";
        c.populationLinearDefect = "mp100";
        auto result=NewtonSolver(mesh_,matdb_,doping_,contactBiases_,c,fixedCharges_,sheetCharges_,contactSpecs_).solvePoissonOnly(initial);
        // Poisson may solve in QF coordinates, but its checkpoint belongs to
        // the caller's population contract. Keep the computed n/p untouched.
        auto& state=result.solution;state.stateVariables="potential_density";
        state.packedState.resize(0);state.packedLow.resize(0);
        state.phinIncrement.resize(0);state.phipIncrement.resize(0);
        state.electronQfReference.resize(0);state.holeQfReference.resize(0);
        return result;
    }
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
    recombinationConfig.bandToBand = cfg_.bandToBand;
    CoupledDDAssembler assembler(
        mesh_, matdb_, doping_, Vt, cfg_.mobility, recombinationConfig,
        cfg_.bandgapNarrowing, cfg_.impactIonization, fixedCharges_,
        sheetCharges_, buildScalingSpec(), cfg_.carrierDiagonalFloor,
        cfg_.carrierStatistics);
    configureAssemblerNumerics(assembler, &initial);
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

    // A nearby ohmic contact across an insulator does not set the carrier
    // reference for this semiconductor. In particular, the polysilicon
    // gate must not seed a 1 V electron quasi-Fermi potential in the channel.
    std::vector<std::vector<Index>> adjacency(N);
    const auto& edgeCells = detail::buildEdgeCellMap(mesh_);
    for (Index edgeId = 0; edgeId < mesh_.numEdges(); ++edgeId) {
        const auto& cells = edgeCells[edgeId];
        const bool conducts = std::any_of(cells.begin(), cells.end(), [&](Index cellId) {
            const auto& region = mesh_.getRegion(mesh_.getCell(cellId).region_id);
            return detail::isTransportMaterial(matdb_.getMaterial(region.material, cfg_.temperature_K));
        });
        if (conducts) {
            const auto& edge = mesh_.getEdge(edgeId);
            adjacency[edge.n0].push_back(edge.n1);
            adjacency[edge.n1].push_back(edge.n0);
        }
    }
    const auto noComponent = static_cast<Index>(N);
    std::vector<Index> component(N, noComponent);
    for (Index root = 0; root < noComponent; ++root) {
        if (component[root] != noComponent || adjacency[root].empty())
            continue;
        std::vector<Index> pending{root};
        component[root] = root;
        for (std::size_t i = 0; i < pending.size(); ++i) {
            for (Index neighbor : adjacency[pending[i]]) {
                if (component[neighbor] == noComponent) {
                    component[neighbor] = root;
                    pending.push_back(neighbor);
                }
            }
        }
    }
    struct ContactAnchor {
        Real x = 0.0;
        Real y = 0.0;
        Real bias = 0.0;
        Real meanNetDoping = 0.0;
    };
    std::vector<std::vector<ContactAnchor>> componentAnchors(N);
    for (const Contact& contact : mesh_.contacts()) {
        const auto biasIt = contactBiases_.find(contact.name);
        if (biasIt == contactBiases_.end() || contact.node_ids.empty())
            continue;
        const auto specIt = contactSpecs_.find(contact.name);
        if (specIt != contactSpecs_.end() &&
            specIt->second.type == ContactType::MetalGate) {
            continue;
        }
        std::map<Index, std::pair<ContactAnchor, std::size_t>> anchors;
        for (Index node : contact.node_ids) {
            if (component[node] == noComponent)
                continue;
            auto& [anchor, count] = anchors[component[node]];
            const Node& geometry = mesh_.getNode(node);
            anchor.x += geometry.x;
            anchor.y += geometry.y;
            anchor.meanNetDoping += doping_.netDoping(node);
            ++count;
        }
        for (auto& [owner, entry] : anchors) {
            auto& [anchor, count] = entry;
            const Real inverseCount = 1.0 / static_cast<Real>(count);
            anchor.x *= inverseCount;
            anchor.y *= inverseCount;
            anchor.meanNetDoping *= inverseCount;
            anchor.bias = biasIt->second;
            componentAnchors[owner].push_back(anchor);
        }
    }
    const auto nearestBias = [&](Index node,
                                 bool electron,
                                 Real fallback) {
        if (component[node] == noComponent)
            return fallback;
        const auto& anchors = componentAnchors[component[node]];
        if (anchors.empty())
            return fallback;
        if (std::all_of(anchors.begin(), anchors.end(), [&](const ContactAnchor& anchor) {
                return anchor.bias == anchors.front().bias;
            }))
            return anchors.front().bias;
        const Node& geometry = mesh_.getNode(node);
        Real bestDistance = std::numeric_limits<Real>::infinity();
        Real bestBias = fallback;
        for (const ContactAnchor& anchor : anchors) {
            if (electron ? anchor.meanNetDoping <= 0.0 : anchor.meanNetDoping >= 0.0)
                continue;
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
    if (cfg_.poissonOnlyQuasiFermiInitialization == "contact_basin") {
        for (int i = 0; i < N; ++i) {
            const Index node = static_cast<Index>(i);
            phin(i) = nearestBias(node, true, initial.phin(i));
            phip(i) = nearestBias(node, false, initial.phip(i));
        }
    }
    for (const auto& [node, value] : bcs.psi)
        psi(static_cast<int>(node)) = value * potentialScale;
    for (const auto& [node, value] : bcs.phin)
        phin(static_cast<int>(node)) = value * potentialScale;
    for (const auto& [node, value] : bcs.phip)
        phip(static_cast<int>(node)) = value * potentialScale;

    VectorXd x;
    if (cfg_.poissonOnlyQuasiFermiInitialization == "preserve_initial") {
        DDSolution projected = initial;
        projected.psi = psi;
        projected.phin = phin;
        projected.phip = phip;
        x = packReferencedSolution(assembler, projected, bcs, true);
    } else {
        // Contact-basin initialization intentionally replaces the imported QFs.
        x = assembler.pack({
            psi / potentialScale,
            phin / potentialScale,
            phip / potentialScale});
    }
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
    BankRoseDamping bankRose(cfg_.bankRose);
    LinearSolver linearSolver;
    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        const auto poissonBlock = [&] {
            if (cfg_.jacobian == "finite_difference")
                return sparseBlock(assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep),
                                   0, 0, N, N);
            if (assembler.usesBvMethodsOperator())
                return assembler.assemblePoissonJacobian(x, bcs);
            return sparseBlock(assembler.assembleJacobian(x, bcs), 0, 0, N, N);
        }();
        VectorXd step = VectorXd::Zero(3 * N);
        try {
            step.head(N) = linearSolver.solve(poissonBlock, -residual.head(N));
        } catch (const std::runtime_error&) {
            finish(false, "poisson_only_linear_solve_failed", iter - 1, x, residual);
            return result;
        }
        applyConfiguredStepCaps(step, cfg_, N, potentialScale, doping_);
        const Real rawStepNorm = step.norm();
        auto searched = cfg_.bankRose.enabled ? bankRose.search(
            iter, x, step, residual,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            }, poissonNorm, {}, cfg_.diagnostics) : lineSearch.search(
            x, step, residual,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            },
            poissonNorm);
        if (!searched.accepted) {
            // Poisson-only initialization keeps both quasi-Fermi potentials
            // fixed at their contact-basin values.  Once the Poisson block and
            // the induced carrier residuals are below the qualified stall
            // floors, another Newton correction can be dominated by floating-
            // point noise and fail the strict-decrease line search.  Reuse the
            // same bounded stall contract as the coupled solver instead of
            // discarding this physically safe equilibrium state.
            const NewtonBlockResidualInfo stalledBlocks =
                blockResidualInfo(residual, mesh_.numNodes());
            if (isPoissonLineSearchStall(
                    searched, stalledBlocks, poissonNorm(residual), 0.0, cfg_)) {
                finish(
                    true, "poisson_only_line_search_stall_floor", iter - 1,
                    x, residual);
                return result;
            }
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
        info.lineSearchHistory = std::move(searched.history);
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
    if (cfg_.stateVariables == "potential_density" &&
        cfg_.electronQuantumPotential.couplingMode != "monolithic_box")
        return solvePopulation(initial);
    if (cfg_.captureLinearSystem && (cfg_.progressDirectory.empty() || cfg_.splitDDState))
        throw std::invalid_argument("Newton linear capture requires classical progress snapshots");
#ifndef VELA_ENABLE_HDF5_STATE
    if (!cfg_.progressDirectory.empty())
        throw std::invalid_argument("Newton progress snapshots require HDF5 state archives");
#endif
    if (!std::isfinite(cfg_.softTimeLimitSeconds) || cfg_.softTimeLimitSeconds<0.)
        throw std::invalid_argument("Newton soft time limit must be finite and nonnegative");
    if ((!cfg_.progressDirectory.empty() || cfg_.softTimeLimitSeconds > 0) &&
        (cfg_.electronQuantumPotential.enabled || cfg_.carrierRowRecovery.mode != "off"))
        throw std::invalid_argument("Newton progress/time limit requires classical solve without recursive recovery");
    if (cfg_.diagnosticStopPolicy != "off" &&
        (cfg_.diagnosticStopPolicy != "row_trend_v1" || cfg_.progressDirectory.empty() ||
         cfg_.carrierRowConvergence.mode != "enforce" || cfg_.splitDDState ||
         cfg_.blockAbsoluteConvergence.mode != "off"))
        throw std::invalid_argument("row_trend_v1 requires classical packed progress, enforced rows and ordinary residual gates");
    if(initial.packedLow.size()!=0 && !cfg_.splitDDState)throw std::invalid_argument("Split checkpoint requires split_dd_state");
    if(cfg_.splitDDState && (!cfg_.warmStart || cfg_.electronQuantumPotential.enabled || cfg_.carrierRowRecovery.mode!="off" || cfg_.temperature_K!=300.))throw std::invalid_argument("Split DD requires a 300 K classical warm start without Gummel recovery");
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
        const Real Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
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
        const bool useSentaurusBox =
            cfg_.electronQuantumPotential.includeInsulators &&
            cfg_.electronQuantumPotential.globalDiscretization == "sentaurus_box";
        for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
            const Cell& cell = mesh_.getCell(cellId);
            const Region& region = mesh_.getRegion(cell.region_id);
            const Material material = matdb_.getMaterial(region.material);
            const bool transport = material.ni > 0.0 && material.mun > 0.0;
            cellTransport[cellId] = transport;
            const Real affinity = material.electron_affinity_eV.value_or(0.0);
            cellElectronAffinity_eV[cellId] = useSentaurusBox
                ? material.electron_quantum_affinity_eV.value_or(affinity)
                : affinity;
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
                        cellElectronAffinity_eV[cellId] + dosMassDrive;
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
            const auto& edgeCells = detail::buildEdgeCellMap(mesh_);
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

        const auto buildBoxInputs = [&](const VectorXd& psi, const VectorXd& phin,
                                        const VectorXd& oldLambda) {
            VectorXd nodeOutputShift(nodeCount);
            for (int i = 0; i < nodeCount; ++i) {
                const Real narrowing = quantumBandgapNarrowing(
                    static_cast<Index>(i), transportNodes[static_cast<std::size_t>(i)]);
                nodeOutputShift(i) = materialBandDrive(i) + psi(i) +
                    cfg_.electronQuantumPotential.conductionBandNarrowingFraction * narrowing;
            }
            std::vector<DensityGradientCellMaterial> cells;
            cells.reserve(mesh_.numCells());
            for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
                const Cell& cell = mesh_.getCell(cellId);
                if (cell.node_ids.size() != 3)
                    continue;
                const bool transport = cellTransport[cellId];
                DensityGradientCellMaterial data;
                data.cellId = cellId;
                data.isTransport = transport;
                data.coefficientVm2 = densityGradientCoefficientVm2(
                    cellQuantumGamma[cellId], cellQuantumCoefficientMassRatio[cellId]);
                for (int local = 0; local < 3; ++local) {
                    const int node = static_cast<int>(cell.node_ids[local]);
                    const Real narrowing = quantumBandgapNarrowing(cell.node_ids[local], transport);
                    // Lambda is reconstructed with the band and DOS data
                    // of this material, including at shared interface nodes.
                    data.materialBandDrive_V[local] = psi(node) +
                        cellElectronAffinity_eV[cellId] + cellDosMassDrive_V[cellId] +
                        cfg_.electronQuantumPotential.conductionBandNarrowingFraction * narrowing;
                    data.dynamicDrivingPotential_V[local] = transport ? -phin(node) : -psi(node);
                    data.initialLambda_V[local] = oldLambda(node);
                }
                cells.push_back(data);
            }
            return std::pair{std::move(cells), std::move(nodeOutputShift)};
        };

        if (cfg_.electronQuantumPotential.couplingMode == "monolithic_box") {
            const VectorXd zero = VectorXd::Zero(nodeCount);
            auto [cells, outputDrive] = buildBoxInputs(zero, zero, zero);
            return solveMonolithicDensityGradientBox(
                initial, std::move(cells), outputDrive, activeNodes, quantumDirichlet);
        }

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
                auto [cellMaterials, nodeOutputShift] = buildBoxInputs(
                    densityGradientState.psi, densityGradientState.phin, quantumPotential);
                if (!continuousLambdaState &&
                    quantumPotentialLike.size() != nodeCount)
                    quantumPotentialLike = quantumPotential - nodeOutputShift;
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
                    " raw_change_V=" + scientificReal(change) +
                    " applied_change_V=" +
                        scientificReal(outerInfo.appliedChangeInfinityNorm_V) +
                    " relaxation=" + scientificReal(relaxation) +
                    " potential_max_V=" +
                        scientificReal(outerInfo.potentialInfinityNorm_V));
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
    const auto progressStarted = std::chrono::steady_clock::now();
    const double Vt = thermalVoltageForInputs(cfg_.temperature_K, cfg_.unitScalingRefs);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg_.recombination, cfg_.taun, cfg_.taup, cfg_.srhDopingDependence);
    recombinationConfig.activeMaterials = cfg_.recombinationMaterials;
    recombinationConfig.volumePolicy = cfg_.recombinationVolumePolicy;
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    recombinationConfig.augerDensityDependence = cfg_.augerDensityDependence;
    recombinationConfig.augerWithGeneration = cfg_.augerWithGeneration;
    recombinationConfig.augerExcessProduct = cfg_.augerExcessProduct;
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
    assembler.setStructureCache(cfg_.sequentialJacobianStructure);
    assembler.setFermiNodeCache(cfg_.diagnosticFermiNodeCache);
    assembler.setElectronQuantumPotential(electronQuantumPotential_V);
    configureAssemblerNumerics(assembler, cfg_.warmStart ? &initial : nullptr);
    if (cfg_.quasiFermiRecenterOnInitialState) {
        assembler.setQuasiFermiReferenceFields(initial.phin, initial.phip);
    }
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
    // A restart written by a referenced solve carries the cancellation-free
    // coordinates explicitly.  Repacking only the physical quasi-Fermi
    // potentials would quantize increments smaller than the ULP of a biased
    // contact.  Restore the saved coordinates for unconstrained nodes when
    // their saved references match this assembler's coordinate system. Subtract
    // references first so an unchanged reference preserves the saved increment.
    if (cfg_.warmStart && initial.hasReferencedElectronQuasiFermi()) {
        for (int i = 0; i < N; ++i) {
            const Index node = static_cast<Index>(i);
            const Real reconstructed =
                initial.electronQuasiFermiAt(i);
            const Real consistencyTolerance =
                32.0 * std::numeric_limits<Real>::epsilon() *
                std::max({Real{1.0}, std::abs(initial.phin(i)),
                          std::abs(reconstructed)});
            if (bcs.phin.find(node) == bcs.phin.end() &&
                std::abs(initial.phin(i) - reconstructed) <=
                    consistencyTolerance) {
                const vela::ArithmeticReal repartitioned =
                    (static_cast<vela::ArithmeticReal>(initial.electronQuasiFermiReferenceAt(i)) -
                     static_cast<vela::ArithmeticReal>(assembler.electronQuasiFermiReferenceAt(node))) +
                    static_cast<vela::ArithmeticReal>(initial.phinIncrement(i));
                x(N + i) = static_cast<Real>(
                    repartitioned / static_cast<vela::ArithmeticReal>(potentialScale));
            }
        }
    }
    if (cfg_.warmStart && initial.hasReferencedHoleQuasiFermi()) {
        for (int i = 0; i < N; ++i) {
            const Index node = static_cast<Index>(i);
            const Real reconstructed = initial.holeQuasiFermiAt(i);
            const Real consistencyTolerance =
                32.0 * std::numeric_limits<Real>::epsilon() *
                std::max({Real{1.0}, std::abs(initial.phip(i)),
                          std::abs(reconstructed)});
            if (bcs.phip.find(node) == bcs.phip.end() &&
                std::abs(initial.phip(i) - reconstructed) <=
                    consistencyTolerance) {
                const vela::ArithmeticReal repartitioned =
                    (static_cast<vela::ArithmeticReal>(initial.holeQuasiFermiReferenceAt(i)) -
                     static_cast<vela::ArithmeticReal>(assembler.holeQuasiFermiReferenceAt(node))) +
                    static_cast<vela::ArithmeticReal>(initial.phipIncrement(i));
                x(2 * N + i) = static_cast<Real>(
                    repartitioned / static_cast<vela::ArithmeticReal>(potentialScale));
            }
        }
    }
    if(cfg_.warmStart && initial.packedLow.size()==0 && initial.hasConsistentPackedState() &&
       initial.packedPotentialScale_V==potentialScale) {
        const VectorXd saved=packReferencedSolution(assembler,initial,bcs);
        for(int i=0;i<N;++i) {
            if(!bcs.psi.contains(i))x(i)=saved(i);
            if(!bcs.phin.contains(i))x(N+i)=saved(N+i);
            if(!bcs.phip.contains(i))x(2*N+i)=saved(2*N+i);
        }
    }
    auto splitRuntime=assembler.splitRuntime();
    experimental::NewtonAuditScope::begin(!splitRuntime && !cfg_.electronQuantumPotential.enabled &&
        !cfg_.quasiFermiRecenterOnStall && cfg_.carrierRowRecovery.mode == "off");
    experimental::NewtonCaptureScope::record([&](auto& capture) {
        // Other formulations need their own state/history identity contracts.
        if (capture.initialized || splitRuntime || cfg_.electronQuantumPotential.enabled ||
            cfg_.quasiFermiRecenterOnStall || cfg_.carrierRowRecovery.mode != "off") {
            capture.unsupported = true;
            return;
        }
        capture.potentialScale = potentialScale;
        capture.electronReference = assembler.electronQuasiFermiReferenceField();
        capture.holeReference = assembler.holeQuasiFermiReferenceField();
        capture.initialized = true;
    });
    if(splitRuntime) {
        if(initial.packedLow.size()!=0)x=splitRuntime->restoreForContinuation(initial);
        splitRuntime->applyBoundary(x,bcs);
    }
    // Clear coordinate-keyed terms whenever the optional QF origin changes.
    // Share exact-state equation terms between row scaling and local gates;
    // global closure uses different boundary semantics and remains independent.
    detail::ContinuityTermCache continuityTerms(assembler, bcs);
    const auto continuityAt = [&](const VectorXd& state) -> const auto& {
        // The high coordinate alone is not a complete key for a split state.
        if (splitRuntime) continuityTerms.clear();
        return continuityTerms.evaluate(state);
    };
    VectorXd r;
    {
        ScopedPerformanceTimer timer("newton.initial_residual");
        r = assembler.residual(x, bcs);
    }
    VectorXd activeRowWeights;
    {
        ScopedPerformanceTimer timer("newton.row_weights");
        activeRowWeights = continuityRowWeights(
            assembler, x, bcs, cfg_.continuityRowScaling,
            cfg_.continuityRowScaling.enabled ? &continuityAt(x) : nullptr);
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

    const auto blockAbsoluteConvergenceSatisfied = [&](const VectorXd& residual) {
        if (cfg_.blockAbsoluteConvergence.mode != "enforce")
            return true;
        const NewtonBlockResidualInfo blocks =
            blockResidualInfo(residual, mesh_.numNodes());
        return blocks.psi <=
                   cfg_.blockAbsoluteConvergence.psiResidualCeiling &&
            blocks.phin <=
                   cfg_.blockAbsoluteConvergence.electronResidualCeiling &&
            blocks.phip <=
                   cfg_.blockAbsoluteConvergence.holeResidualCeiling;
    };
    const auto contactMajorityQfBranchDrop = [&](const VectorXd& state,
                                                  int iterations) {
        return maxContactMajorityQuasiFermiDrop(
            makeSolution(assembler, state, iterations),
            cfg_.contactMajorityQfBranchGuardContacts);
    };
    const auto contactMajorityQfBranchAccepts = [&](const VectorXd& state,
                                                    int iterations) {
        if (cfg_.contactMajorityQfBranchDropLimit_V <= 0.0)
            return true;
        const Real drop = contactMajorityQfBranchDrop(state, iterations);
        return std::isfinite(drop) &&
            drop <= cfg_.contactMajorityQfBranchDropLimit_V;
    };
    const auto updateBestAcceptedIterate = [&](const VectorXd& state,
                                               const VectorXd& residual,
                                               int iterations,
                                               Real norm) {
        const Real drop = contactMajorityQfBranchDrop(state, iterations);
        const bool branchEligible =
            cfg_.contactMajorityQfBranchDropLimit_V <= 0.0 ||
            (std::isfinite(drop) &&
             drop <= cfg_.contactMajorityQfBranchDropLimit_V);
        if (!branchEligible || !std::isfinite(norm) ||
            (result.hasBestSolution && !(norm < result.bestResidualNorm))) {
            return;
        }
        result.hasBestSolution = true;
        result.bestIteration = iterations;
        result.bestResidualNorm = norm;
        result.bestBlockNorms = blockResidualInfo(residual, mesh_.numNodes());
        result.bestContactMajorityQfDrop = drop;
        result.bestSolution = makeSolution(assembler, state, iterations);
    };
    updateBestAcceptedIterate(x, r, 0, initialNorm);

    auto carrierRowEval = [&](const VectorXd& state) {
        ScopedPerformanceTimer timer("newton.carrier_row_evaluation");
        if (cfg_.carrierRowConvergence.mode == "off")
            return NewtonCarrierRowConvergenceEvaluation{};
        return evaluateCarrierRowConvergence(
            continuityAt(state),
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
    auto freeCarrierResidualSums = [&](const VectorXd& residual) {
        long double electronSum = 0.0L;
        long double holeSum = 0.0L;
        for (int i = 0; i < N; ++i) {
            const Index node = static_cast<Index>(i);
            if (bcs.phin.find(node) == bcs.phin.end())
                electronSum += static_cast<long double>(residual(N + i));
            if (bcs.phip.find(node) == bcs.phip.end())
                holeSum += static_cast<long double>(residual(2 * N + i));
        }
        return std::array<long double, 2>{electronSum, holeSum};
    };
    Real activeGlobalElectronScale = cfg_.globalContinuityClosure.sourceFloor;
    Real activeGlobalHoleScale = cfg_.globalContinuityClosure.sourceFloor;
    auto globalClosureLineSearchNorm = [&](const VectorXd& residual) {
        const Real baseNorm = residualNormFn(residual);
        // Global closure is a convergence invariant, not an implicit path
        // constraint.  Keep ordinary merit and device-block line searches on
        // the device residual alone.  A caller that deliberately wants
        // closure to shape the nonlinear path must opt in through the
        // block_closure_filter mode.
        if (cfg_.globalContinuityClosure.mode == "off" ||
            cfg_.lineSearchMode != "block_closure_filter")
            return baseNorm;
        const auto sums = freeCarrierResidualSums(residual);
        const Real electronRatio =
            static_cast<Real>(sums[0]) /
            std::max(activeGlobalElectronScale, Real{1.0e-300});
        const Real holeRatio =
            static_cast<Real>(sums[1]) /
            std::max(activeGlobalHoleScale, Real{1.0e-300});
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
            out << "event,iteration,node_id,carrier,residual,scale,ratio,carrier_density_m3,"
                << "flux,srh,impact,density_qualified,flux_qualified,source_qualified\n";
        }
        out << std::setprecision(17);
        for (const auto& row : evaluation.violations) {
            out << event << ',' << iteration << ',' << row.nodeId << ',' << row.carrier << ','
                << row.residual << ',' << row.scale << ',' << row.ratio << ','
                << row.carrierDensity_m3 << ',' << row.flux << ','
                << row.recombination << ',' << row.impact << ','
                << (row.densityQualified ? 1 : 0) << ','
                << (row.fluxQualified ? 1 : 0) << ','
                << (row.sourceQualified ? 1 : 0) << '\n';
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
    const bool closureAuditEnabled = !cfg_.globalContinuityClosure.auditJsonlFile.empty();
    auto writeClosureAudit = [&](nlohmann::json record) {
        if (!closureAuditEnabled) return;
        const std::filesystem::path path(cfg_.globalContinuityClosure.auditJsonlFile);
        if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
        std::ofstream out(path, std::ios::app);
        if (!out) throw std::runtime_error("Cannot open global continuity audit: " + path.string());
        record["schema_version"] = 1;
        record["restoration_jacobian"] = cfg_.globalContinuityClosure.restorationJacobian;
        record["restoration_physical_rows_only"] = cfg_.globalContinuityClosure.restorationPhysicalRowsOnly;
        out << record.dump() << '\n';
    };
    auto auditReplacementRows = [&](int iteration, const VectorXd& beforeX,
                                    const VectorXd& beforeR, const VectorXd& afterX,
                                    const VectorXd& afterR) {
        if (!closureAuditEnabled || cfg_.globalContinuityClosure.linearSystemStabilization !=
            "replace_free_row_sum") return;
        std::vector<bool> ec(N, false), hc(N, false);
        for (Index node : electronContactNodes) ec[node] = true;
        for (Index node : holeContactNodes) hc[node] = true;
        const auto selected = selectGlobalContinuityReplacementRows(SparseMatrixd{}, N, ec, hc);
        const auto beforeN = assembler.electronDensity(beforeX);
        const auto beforeP = assembler.holeDensity(beforeX);
        const auto afterN = assembler.electronDensity(afterX);
        const auto afterP = assembler.holeDensity(afterX);
        for (int carrier = 0; carrier < 2; ++carrier) {
            const int row = selected[carrier];
            if (row < 0) continue;
            const int node = row - (carrier + 1) * N;
            auto rank = [&](const VectorXd& residual) {
                int value = 1;
                for (int i = 0; i < N; ++i)
                    if (std::abs(residual((carrier + 1) * N + i)) > std::abs(residual(row))) ++value;
                return value;
            };
            const Real ni = assembler.intrinsicDensity()[node];
            const Real density = carrier == 0 ? afterN(node) : afterP(node);
            const auto& position = mesh_.nodes()[node];
            writeClosureAudit({{"event", "replacement_row"}, {"iteration", iteration},
                {"carrier", carrier == 0 ? "electron" : "hole"}, {"row_index", row},
                {"node_id", node}, {"x_m", scaling.unitSystem.internalLengthToMeters(position.x)},
                {"y_m", scaling.unitSystem.internalLengthToMeters(position.y)},
                {"density_before_internal", carrier == 0 ? beforeN(node) : beforeP(node)},
                {"density_after_internal", density}, {"intrinsic_density_internal", ni},
                {"density_below_intrinsic", density < ni},
                {"residual_before", beforeR(row)}, {"residual_after", afterR(row)},
                {"rank_before", rank(beforeR)}, {"rank_after", rank(afterR)},
                {"top10_before", rank(beforeR) <= 10}, {"top10_after", rank(afterR) <= 10}});
        }
    };
    auto shouldWriteLocalUpdateTrace = [&](int iteration) {
        if (!cfg_.localUpdateDiagnostics.enabled ||
            cfg_.localUpdateDiagnostics.csvFile.empty() ||
            cfg_.localUpdateDiagnostics.nodes.empty()) {
            return false;
        }
        if (iteration <= cfg_.localUpdateDiagnostics.firstIterations)
            return true;
        return iteration % std::max(1, cfg_.localUpdateDiagnostics.everyIterations) == 0;
    };
    auto writeLocalUpdateTraceCsv = [&](int iteration,
                                        const VectorXd& state,
                                        const VectorXd& currentResidual,
                                        const SparseMatrixd& jacobian,
                                        const VectorXd& rawStep,
                                        const VectorXd& cappedStep,
                                        const VectorXd& rawLinearResidual,
                                        const VectorXd& cappedLinearResidual,
                                        const LineSearchResult& lineSearchResult,
                                        const VectorXd& auditBaseLow) {
        if (!shouldWriteLocalUpdateTrace(iteration))
            return;
        ScopedPerformanceTimer diagnosticTimer("newton.local_update_diagnostics");
        const std::filesystem::path path(cfg_.localUpdateDiagnostics.csvFile);
        if (!path.parent_path().empty())
            std::filesystem::create_directories(path.parent_path());
        const bool writeHeader = !std::filesystem::exists(path);
        std::ofstream out(path, std::ios::app);
        if (writeHeader) {
            out << "iteration,node_id,carrier,row_index,line_search_accepted,line_search_attempts,"
                << "reported_damping,selected_damping,best_rejected_candidate,"
                << "qf_reference_V,state_increment_V,state_physical_qf_V,residual,rhs,"
                << "jacobian_diagonal,row_weight,raw_linear_step_V,capped_step_V,applied_step_V,"
                << "raw_linear_residual,capped_linear_residual,selected_trial_residual,"
                << "raw_linear_residual_l2,raw_linear_residual_inf,"
                << "capped_linear_residual_l2,capped_linear_residual_inf,"
                << "current_psi_block,current_electron_block,current_hole_block,"
                << "raw_predicted_psi_block,raw_predicted_electron_block,raw_predicted_hole_block,"
                << "capped_predicted_psi_block,capped_predicted_electron_block,capped_predicted_hole_block,"
                << "selected_psi_block,selected_electron_block,selected_hole_block,"
                << "raw_step_l2,capped_step_l2,applied_step_l2,"
                << "qf_capped_dofs,qf_max_cap_delta_V,"
                << "current_electron_free_sum,current_hole_free_sum,"
                << "selected_electron_free_sum,selected_hole_free_sum,"
                << "qf_update_limit_V,qf_update_limit_minority_V,max_update,"
                << "carrier_regularization_scale,carrier_diagonal_floor_enabled,"
                << "carrier_diagonal_floor_scale,jacobian_inf_norm,raw_step_inf_norm,"
                << "rhs_inf_norm,raw_linear_relative_residual_inf,"
                << "raw_linear_backward_error_inf,raw_linear_componentwise_backward_error,"
                << "componentwise_max_row,componentwise_max_linear_residual,"
                << "componentwise_max_denominator,componentwise_max_rhs,"
                << "selected_coordinate_change_V,selected_high_change_V,low_before_internal,low_after_internal,"
                << "selected_writeback_error_V,complete_coordinate_change_available\n";
        }
        VectorXd appliedStep = VectorXd::Zero(state.size());
        VectorXd selectedResidual = currentResidual;
        Real selectedDamping = 0.0;
        if (lineSearchResult.accepted) {
            appliedStep = lineSearchResult.x - state;
            selectedResidual = lineSearchResult.residual;
            selectedDamping = lineSearchResult.damping;
        } else if (lineSearchResult.bestRejectedCandidate) {
            appliedStep = lineSearchResult.bestRejectedX - state;
            selectedResidual = lineSearchResult.bestRejectedResidual;
            selectedDamping = lineSearchResult.bestRejectedDamping;
        }
        // Rejected split candidates do not retain their low state. Never label
        // their high-only difference a complete applied coordinate change.
        const bool fullChangeAvailable=(lineSearchResult.accepted || lineSearchResult.bestRejectedCandidate) &&
            (!splitRuntime || lineSearchResult.accepted);
        const VectorXd fullChange=(splitRuntime && lineSearchResult.accepted)
            ? splitRuntime->coordinateDelta(lineSearchResult.x,state,auditBaseLow) : appliedStep;
        const VectorXd auditAfterLow=(splitRuntime && lineSearchResult.accepted)?splitRuntime->low():VectorXd{};
        const Real rawLinearResidualInf = rawLinearResidual.size() > 0
            ? rawLinearResidual.lpNorm<Eigen::Infinity>() : 0.0;
        const Real cappedLinearResidualInf = cappedLinearResidual.size() > 0
            ? cappedLinearResidual.lpNorm<Eigen::Infinity>() : 0.0;
        const ResidualBlockNormValue currentBlocks =
            ResidualNorm::computeBlocks(currentResidual, mesh_.numNodes());
        const ResidualBlockNormValue rawPredictedBlocks =
            ResidualNorm::computeBlocks(rawLinearResidual, mesh_.numNodes());
        const ResidualBlockNormValue cappedPredictedBlocks =
            ResidualNorm::computeBlocks(cappedLinearResidual, mesh_.numNodes());
        const ResidualBlockNormValue selectedBlocks =
            ResidualNorm::computeBlocks(selectedResidual, mesh_.numNodes());
        int qfCappedDofs = 0;
        Real qfMaxCapDelta_V = 0.0;
        for (int row = N; row < 3 * N; ++row) {
            const Real difference = std::abs(rawStep(row) - cappedStep(row));
            const Real scale = std::max(
                {Real{1.0}, std::abs(rawStep(row)), std::abs(cappedStep(row))});
            if (difference > 64.0 * std::numeric_limits<Real>::epsilon() * scale)
                ++qfCappedDofs;
            qfMaxCapDelta_V = std::max(
                qfMaxCapDelta_V, difference * potentialScale);
        }
        const auto currentFreeSums = freeCarrierResidualSums(currentResidual);
        const auto selectedFreeSums = freeCarrierResidualSums(selectedResidual);
        VectorXd rowAbsSums = VectorXd::Zero(jacobian.rows());
        VectorXd componentDenominator = currentResidual.cwiseAbs();
        for (int column = 0; column < jacobian.outerSize(); ++column) {
            for (SparseMatrixd::InnerIterator entry(jacobian, column); entry; ++entry) {
                const Real magnitude = std::abs(entry.value());
                rowAbsSums(entry.row()) += magnitude;
                componentDenominator(entry.row()) += magnitude * std::abs(rawStep(entry.col()));
            }
        }
        const Real jacobianInf = rowAbsSums.maxCoeff();
        const Real rawStepInf = rawStep.lpNorm<Eigen::Infinity>();
        const Real rhsInf = currentResidual.lpNorm<Eigen::Infinity>();
        const auto safeRatio = [](Real numerator, Real denominator) {
            return denominator > 0.0 ? numerator / denominator
                : (numerator == 0.0 ? 0.0 : std::numeric_limits<Real>::infinity());
        };
        Real componentBackwardError = 0.0;
        int componentMaxRow = 0;
        for (int row = 0; row < rawLinearResidual.size(); ++row) {
            const Real candidate =
                safeRatio(std::abs(rawLinearResidual(row)), componentDenominator(row));
            if (candidate > componentBackwardError) {
                componentBackwardError = candidate;
                componentMaxRow = row;
            }
        }
        out << std::setprecision(17);
        for (Index node : cfg_.localUpdateDiagnostics.nodes) {
            if (node < 0 || node >= mesh_.numNodes())
                continue;
            for (int carrierBlock = 1; carrierBlock <= 2; ++carrierBlock) {
                const int row = carrierBlock * N + static_cast<int>(node);
                const bool electron = carrierBlock == 1;
                const Real reference = electron
                    ? assembler.electronQuasiFermiReferenceAt(node)
                    : assembler.holeQuasiFermiReferenceAt(node);
                const Real stateIncrement = state(row) * potentialScale;
                const Real rowWeight = row < activeRowWeights.size()
                    ? activeRowWeights(row) : 1.0;
                out << iteration << ',' << node << ','
                    << (electron ? "electron" : "hole") << ',' << row << ','
                    << (lineSearchResult.accepted ? 1 : 0) << ','
                    << lineSearchResult.attempts << ',' << lineSearchResult.damping << ','
                    << selectedDamping << ','
                    << (lineSearchResult.bestRejectedCandidate ? 1 : 0) << ','
                    << reference << ',' << stateIncrement << ','
                    << reference + stateIncrement << ','
                    << currentResidual(row) << ',' << -currentResidual(row) << ','
                    << jacobian.coeff(row, row) << ',' << rowWeight << ','
                    << rawStep(row) * potentialScale << ','
                    << cappedStep(row) * potentialScale << ','
                    << appliedStep(row) * potentialScale << ','
                    << rawLinearResidual(row) << ','
                    << cappedLinearResidual(row) << ','
                    << selectedResidual(row) << ','
                    << rawLinearResidual.norm() << ',' << rawLinearResidualInf << ','
                    << cappedLinearResidual.norm() << ',' << cappedLinearResidualInf << ','
                    << currentBlocks.psi << ',' << currentBlocks.phin << ','
                    << currentBlocks.phip << ',' << rawPredictedBlocks.psi << ','
                    << rawPredictedBlocks.phin << ',' << rawPredictedBlocks.phip << ','
                    << cappedPredictedBlocks.psi << ',' << cappedPredictedBlocks.phin << ','
                    << cappedPredictedBlocks.phip << ',' << selectedBlocks.psi << ','
                    << selectedBlocks.phin << ',' << selectedBlocks.phip << ','
                    << rawStep.norm() << ',' << cappedStep.norm() << ','
                    << appliedStep.norm() << ',' << qfCappedDofs << ','
                    << qfMaxCapDelta_V << ','
                    << static_cast<Real>(currentFreeSums[0]) << ','
                    << static_cast<Real>(currentFreeSums[1]) << ','
                    << static_cast<Real>(selectedFreeSums[0]) << ','
                    << static_cast<Real>(selectedFreeSums[1]) << ','
                    << cfg_.quasiFermiUpdateLimit_V << ','
                    << cfg_.quasiFermiUpdateLimitMinority_V << ','
                    << cfg_.maxUpdate << ',' << cfg_.carrierRegularizationScale << ','
                    << (cfg_.carrierDiagonalFloor.enabled ? 1 : 0) << ','
                    << cfg_.carrierDiagonalFloor.scale << ',' << jacobianInf << ','
                    << rawStepInf << ',' << rhsInf << ','
                    << safeRatio(rawLinearResidualInf, rhsInf) << ','
                    << safeRatio(rawLinearResidualInf, jacobianInf * rawStepInf + rhsInf) << ','
                    << componentBackwardError << ',' << componentMaxRow << ','
                    << rawLinearResidual(componentMaxRow) << ','
                    << componentDenominator(componentMaxRow) << ','
                    << -currentResidual(componentMaxRow) << ',';
                if(fullChangeAvailable)out<<fullChange(row)*potentialScale;
                out<<','<<appliedStep(row)*potentialScale<<',';
                if(auditBaseLow.size()==state.size())out<<auditBaseLow(row);
                out<<',';
                if(auditAfterLow.size()==state.size())out<<auditAfterLow(row);
                out<<',';
                if(fullChangeAvailable)out<<(fullChange(row)-selectedDamping*cappedStep(row))*potentialScale;
                out<<','<<fullChangeAvailable<<'\n';
            }
        }
    };
    auto finishConverged = [&](const std::string& reason,
                               const VectorXd& state,
                               const VectorXd& residual,
                               int iterations,
                               Real norm,
                               const NewtonCarrierRowConvergenceEvaluation& rowEval,
                               const NewtonGlobalContinuityClosureEvaluation& globalEval) {
        updateBestAcceptedIterate(state, residual, iterations, norm);
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
    const auto contactMajorityQfAcceptsFloor = [&](const VectorXd& state,
                                                    int iterations) {
        if (cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V <= 0.0)
            return true;
        const DDSolution floorSolution =
            makeSolution(assembler, state, iterations);
        const Real drop = maxContactMajorityQuasiFermiDrop(floorSolution);
        return std::isfinite(drop) &&
            drop <= cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V;
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
        std::vector<NewtonCarrierRowConvergenceViolation> blockingRows;
        for (const auto& row : rowEval.violations)
            if (row.blocksAcceptance) blockingRows.push_back(row);
        NewtonCarrierRowRecoveryResult recovery = recoverCarrierRowsWithGummelDensity(
            mesh_, matdb_, doping_, contactBiases_, cfg_, base, blockingRows,
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
        // A recovery retry is part of this solve. Keep all accepted Newton
        // updates in the public work count, including those before recovery.
        retried.iters += iterations;
        retried.solution.iters = retried.iters;
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
    // A closure-filter run compares all trial iterates against one physical
    // current scale for this bias point. Recomputing the denominator after
    // every accepted step creates a moving target while contact current is
    // being rebalanced, so freeze it from the supplied restart state.
    activeGlobalElectronScale = std::max(
        {std::abs(initialGlobalEval.electron.contactFlux),
         std::abs(initialGlobalEval.electron.integratedSource),
         cfg_.globalContinuityClosure.sourceFloor});
    activeGlobalHoleScale = std::max(
        {std::abs(initialGlobalEval.hole.contactFlux),
         std::abs(initialGlobalEval.hole.integratedSource),
         cfg_.globalContinuityClosure.sourceFloor});
    bool hasBestFeasibleIterate =
        carrierRowsAcceptConvergence(initialRowEval) &&
        globalClosureAcceptsConvergence(initialGlobalEval);
    VectorXd bestFeasibleX = x;
    VectorXd bestFeasibleR = r;
    int bestFeasibleIters = 0;
    Real bestFeasibleNorm = initialNorm;
    NewtonCarrierRowConvergenceEvaluation bestFeasibleRowEval = initialRowEval;
    NewtonGlobalContinuityClosureEvaluation bestFeasibleGlobalEval = initialGlobalEval;
    NewtonIterationInfo initialTrace;
    initialTrace.iter = 0;
    initialTrace.residualNorm = initialNorm;
    initialTrace.relativeResidualNorm = 1.0;
    initialTrace.contactMajorityQfDrop =
        contactMajorityQfBranchDrop(x, 0);
    initialTrace.lineSearchAccepted = true;
    initialTrace.event = "initial";
    initialTrace.carrierRowConvergence = initialRowEval;
    fillNewtonIterationResidualTrace(
        initialTrace, x, r, activeRowWeights, mesh_.numNodes(), assembler,
        cfg_.diagnostics);
    detail::NewtonTrendStop trendStop;
    unsigned linearCaptureSequence = 0;
    std::string pendingTrendStop;
#include "NewtonShadowCapture.inl"
    auto writeProgress = [&](const NewtonIterationInfo& info,
                             const NewtonGlobalContinuityClosureEvaluation& closure,
                             const VectorXd& actualDelta,
                             const VectorXd& rawDirection,
                             const VectorXd& cappedDirection) {
        if (cfg_.progressDirectory.empty()) return;
        ScopedPerformanceTimer snapshotTimer("newton.progress_snapshot");
        const std::filesystem::path folder(cfg_.progressDirectory);
        std::filesystem::create_directories(folder);
        const auto stem = folder / ("iteration_" + std::to_string(info.iter));
        const auto statePath = stem.string() + ".h5";
        const auto jsonPath = stem.string() + ".json";
        if (std::filesystem::exists(statePath) || std::filesystem::exists(jsonPath))
            throw std::runtime_error("Refusing to overwrite Newton progress evidence");
        nlohmann::json doc = {{"schema", "vela.newton_progress.v2"},
            {"qualified", false}, {"iteration", info.iter}, {"event", info.event},
            {"elapsed_s", std::chrono::duration<double>(std::chrono::steady_clock::now()-progressStarted).count()},
            {"residual_norm", info.residualNorm}, {"raw_direction_norm", rawDirection.norm()},
            {"capped_direction_norm", cappedDirection.norm()}, {"potential_scale_V",potentialScale},
            {"applied_update_norm", info.stepNorm}, {"alpha", info.dampingFactor},
            {"line_search_attempts", info.lineSearchAttempts},
            {"blocks", {{"psi",info.blockResiduals.psi},{"phin",info.blockResiduals.phin},{"phip",info.blockResiduals.phip}}},
            {"row_max_ratio",info.carrierRowConvergence.maxRatio},
            {"electron_closure_ratio",closure.electron.ratio}, {"hole_closure_ratio",closure.hole.ratio},
            {"state_file",statePath}, {"violations",nlohmann::json::array()}};
        if (cfg_.bankRose.enabled) {
            doc["bank_rose"] = {{"not_damped", cfg_.bankRose.notDamped},
                {"active", info.iter > cfg_.bankRose.notDamped},
                {"trials", nlohmann::json::array()}};
            for (const auto& trial : info.lineSearchHistory)
                doc["bank_rose"]["trials"].push_back({{"K", trial.bankRoseK},
                    {"alpha", trial.damping}, {"residual_norm", trial.residualNorm},
                    {"target_norm", trial.targetResidualNorm}, {"accepted", trial.accepted},
                    {"reason", trial.rejectionReason}});
        }
        for (const auto& row : info.carrierRowConvergence.violations)
            doc["violations"].push_back({{"node",row.nodeId},{"carrier",row.carrier},
                {"residual",row.residual},{"scale",row.scale},{"ratio",row.ratio},
                {"density_m3",row.carrierDensity_m3}});
        for (int block=0;block<3;++block) {
            Eigen::Index index=0;
            const Real peak=actualDelta.segment(block*N,N).cwiseAbs().maxCoeff(&index);
            doc["actual_update_blocks"][std::to_string(block)]={{"node",index},{"max_abs_internal",peak}};
            const Real rawPeak=rawDirection.segment(block*N,N).cwiseAbs().maxCoeff(&index);
            doc["raw_direction_blocks"][std::to_string(block)]={{"node",index},{"max_abs_internal",rawPeak}};
        }
        int clipped=0,roundedAway=0;
        for(int i=0;i<x.size();++i) {
            if(rawDirection[i]!=cappedDirection[i]) ++clipped;
            if(info.dampingFactor*cappedDirection[i]!=0. && actualDelta[i]==0.) ++roundedAway;
        }
        doc["cap_changed_dofs"]=clipped;doc["requested_nonzero_but_actual_zero_dofs"]=roundedAway;
        if (cfg_.diagnosticStopPolicy == "row_trend_v1") {
            detail::NewtonTrendSample sample;
            sample.iteration=info.iter; sample.alpha=info.dampingFactor;
            sample.rawDirection=rawDirection.norm();
            for (const auto& row:info.carrierRowConvergence.violations) {
                const int offset=row.carrier=="electron" ? 0 : 3;
                sample.rows[offset]+=1;
                sample.rows[offset+1]+=std::max(0.0,std::log10(row.ratio/cfg_.carrierRowConvergence.epsRow));
                sample.rows[offset+2]=std::max(sample.rows[offset+2],row.ratio);
            }
            sample.gates[0]=std::min(info.residualNorm/cfg_.abstol,
                ResidualNorm::relative(info.residualNorm,initialNorm)/cfg_.reltol);
            auto closureRatio=[&](const NewtonGlobalContinuityCarrierClosure& c) {
                if (!closure.enforced || c.satisfied) return 0.0;
                return c.qualified ? c.ratio/cfg_.globalContinuityClosure.tolerance :
                    std::abs(c.mismatch)/cfg_.globalContinuityClosure.absoluteTolerance;
            };
            sample.gates[1]=closureRatio(closure.electron);
            sample.gates[2]=closureRatio(closure.hole);
            if (cfg_.contactMajorityQfBranchDropLimit_V>0)
                sample.gates[3]=info.contactMajorityQfDrop/cfg_.contactMajorityQfBranchDropLimit_V;
            pendingTrendStop=trendStop.observe(sample);
            doc["trend_stop"]={{"policy",cfg_.diagnosticStopPolicy},{"rows_N_E_max",sample.rows},
                {"gates_primary_eclosure_hclosure_qf",sample.gates},
                {"consecutive_stalled_windows",trendStop.consecutiveWindows()},
                {"pending_reason",pendingTrendStop},{"success_checks_have_priority",true}};
        }
        auto metadata = activeDDStateArchiveMetadata() ? *activeDDStateArchiveMetadata() : nlohmann::json::object();
        metadata["state_role"]="newton_iteration_diagnostic";
        metadata["qualified"]=false;
        DDStateArchiveScope scope(std::move(metadata));
        appendShadow(doc, info.iter, rawDirection, cappedDirection, info.dampingFactor, info, closure);
        writeDDSolutionState(statePath, makeSolution(assembler,x,info.iter),cfg_.inputScaling);
        // Unique files: a hard kill can leave an incomplete final file, never corrupt a prior snapshot.
        std::ofstream stream(jsonPath); stream << doc.dump(2) << '\n'; stream.flush();
        if (!stream) throw std::runtime_error("Newton progress write failed");
    };
    writeProgress(initialTrace, initialGlobalEval, VectorXd::Zero(x.size()),
                  VectorXd::Zero(x.size()),VectorXd::Zero(x.size()));
    result.trace.push_back(std::move(initialTrace));
    writeCarrierRowTraceCsv(x, initialRowEval, 0, initialNorm, "initial");
    const bool initialBlockAbsoluteConverged =
        cfg_.blockAbsoluteConvergence.mode == "enforce" &&
        blockAbsoluteConvergenceSatisfied(r);
    const bool initialPrimaryToleranceConverged =
        initialBlockAbsoluteConverged ||
        (cfg_.blockAbsoluteConvergence.mode == "off" &&
         initialNorm <= cfg_.abstol);
    if (initialPrimaryToleranceConverged &&
        carrierRowsAcceptConvergence(initialRowEval) &&
        globalClosureAcceptsConvergence(initialGlobalEval) &&
        contactMajorityQfBranchAccepts(x, 0)) {
        finishConverged(
            initialBlockAbsoluteConverged
                ? "initial_block_abstol"
                : "initial_abstol",
            x, r, 0, initialNorm,
            initialRowEval, initialGlobalEval);
        return result;
    }
    if (initialNorm <= cfg_.stallResidualFloor &&
        blockAbsoluteConvergenceSatisfied(r) &&
        carrierRowsAcceptConvergence(initialRowEval) &&
        globalClosureAcceptsConvergence(initialGlobalEval) &&
        contactMajorityQfBranchAccepts(x, 0) &&
        contactMajorityQfAcceptsFloor(x, 0)) {
        finishConverged(
            "initial_stall_residual_floor", x, r, 0, initialNorm,
            initialRowEval, initialGlobalEval);
        return result;
    }
    if (initialPrimaryToleranceConverged && initialRowEval.enforced &&
        !initialRowEval.satisfied) {
        writeCarrierRowDiagnosticCsv(
            initialRowEval, 0, "carrier_row_convergence_initial_abstol_rejected");
        writeCarrierRowTraceCsv(
            x, initialRowEval, 0, initialNorm,
            "carrier_row_convergence_initial_abstol_rejected");
        if (auto recovered = retryAfterCarrierRowRecovery(x, 0, initialRowEval))
            return *recovered;
    }
    std::unique_ptr<LinearSolver> localLinearSolver;
    if (!cfg_.sequentialLinearSolver) localLinearSolver = std::make_unique<LinearSolver>();
    LinearSolver& linearSolver = cfg_.sequentialLinearSolver
        ? *cfg_.sequentialLinearSolver : *localLinearSolver;
    LineSearchConfig lscfg;
    lscfg.enabled = cfg_.lineSearch;
    lscfg.initialDamping = cfg_.dampingFactor;
    lscfg.recordHistory = cfg_.diagnostics;
    BacktrackingLineSearch lineSearch(lscfg);
    BankRoseDamping bankRose(cfg_.bankRose);

    // Stall-recovery threshold: when the damped Newton step cannot reduce the
    // residual but the residual already sits at or below this normalized floor,
    // the state is effectively solved (the line search is fighting numerical
    // noise) and convergence is reported instead of failing.
    const Real stallResidualFloor = cfg_.stallResidualFloor;


    Real quasiFermiTrustRegionMultiplier = 1.0;
    auto capConfiguredStep = [&](VectorXd& candidateStep,
                                 const SparseMatrixd& jacobian,
                                 const VectorXd& residual) {
        NewtonConfig activeLimits = cfg_;
        activeLimits.quasiFermiUpdateLimit_V *= quasiFermiTrustRegionMultiplier;
        activeLimits.quasiFermiUpdateLimitMinority_V *=
            quasiFermiTrustRegionMultiplier;
        applyConfiguredStepCapsAndPoissonRecorrection(
            candidateStep, jacobian, residual, activeLimits,
            N, potentialScale, doping_);
    };
    auto capConfiguredQfOnly = [&](VectorXd& candidateStep) {
        NewtonConfig activeLimits = cfg_;
        activeLimits.quasiFermiUpdateLimit_V *= quasiFermiTrustRegionMultiplier;
        activeLimits.quasiFermiUpdateLimitMinority_V *=
            quasiFermiTrustRegionMultiplier;
        applyConfiguredQuasiFermiStepCaps(
            candidateStep, activeLimits, N, potentialScale, doping_);
    };
    auto limitPopulationStep = [&](VectorXd& direction) {
        if(cfg_.populationStepFraction==0.)return;
        const auto f=population_step::frame(assembler,x,direction,bcs);
        const auto bound=population_step::relativeFractionToBoundary(
            f.relativeN,f.relativeP,f.active,cfg_.populationStepFraction);
        direction*=bound.alpha;
        if(cfg_.diagnostics)std::clog<<"global_population_step factor="<<std::setprecision(17)
            <<bound.alpha<<" limiting_node="<<bound.limitingNode
            <<" limiting_carrier="<<bound.limitingCarrier<<'\n';
    };

    if (closureAuditEnabled && cfg_.globalContinuityClosure.roundoffAudit) {
        auto audit = auditGlobalContinuityRoundoff(assembler, x, bcs,
            cfg_.globalContinuityClosure, mesh_.numEdges());
        audit["stage"] = "initial";
        writeClosureAudit(std::move(audit));
    }
    auditReplacementRows(0, x, r, x, r);
    VectorXd acceptedX = x;
    VectorXd acceptedR = r;
    int acceptedIters = 0;
    NewtonGlobalContinuityClosureEvaluation activeGlobalEval =
        initialGlobalEval;

    int quasiFermiRecenters = 0;

    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        experimental::CacheAuditScope::iteration(iter);
        ScopedPerformanceTimer iterationTimer("newton.iteration");
        if (cfg_.softTimeLimitSeconds > 0 &&
            std::chrono::duration<double>(std::chrono::steady_clock::now()-progressStarted).count() >= cfg_.softTimeLimitSeconds) {
            result.converged=false;
            result.solution=makeSolution(assembler,x,acceptedIters);
            result.iters=acceptedIters;
            result.finalResidualNorm=residualNormFn(r);
            result.finalBlockNorms=result.trace.back().blockResiduals;
            result.finalCarrierRowConvergence=result.trace.back().carrierRowConvergence;
            result.finalGlobalContinuityClosure=activeGlobalEval;
            result.failureDiagnostics.failureReason="soft_time_limit";
            return result;
        }
        {
            ScopedPerformanceTimer timer("newton.row_weights");
            activeRowWeights = continuityRowWeights(
                assembler, x, bcs, cfg_.continuityRowScaling,
                cfg_.continuityRowScaling.enabled ? &continuityAt(x) : nullptr);
        }
        std::filesystem::path linearCaptureFolder;
        VectorXd captureReferences;
        if (cfg_.captureLinearSystem) {
            linearCaptureFolder = std::filesystem::path(cfg_.progressDirectory) /
                ("linear_" + std::to_string(++linearCaptureSequence));
            if (!std::filesystem::create_directory(linearCaptureFolder))
                throw std::runtime_error("Refusing to overwrite Newton linear capture directory");
            captureReferences = VectorXd::Zero(x.size());
            for (int i=0;i<N;++i) {
                captureReferences(N+i)=assembler.electronQuasiFermiReferenceAt(i);
                captureReferences(2*N+i)=assembler.holeQuasiFermiReferenceAt(i);
            }
            std::ofstream metadata(linearCaptureFolder/"identity.json");
            metadata.exceptions(std::ios::badbit|std::ios::failbit);
            metadata << nlohmann::json{{"schema","vela.newton_linear_capture.v1"},
                {"iteration",iter},{"sequence",linearCaptureSequence},{"potential_scale_V",potentialScale},
                {"qualified",false},{"regularization",cfg_.carrierRegularizationScale},
                {"stabilization",cfg_.globalContinuityClosure.linearSystemStabilization}}.dump(2);
        }
        const auto captureVectors = [&](const char* name, const VectorXd& state,
                                       const VectorXd& residual, const VectorXd& direction) {
            if (cfg_.captureLinearSystem)
                detail::captureNewtonVectors((linearCaptureFolder/name).string(),state,residual,direction,captureReferences);
        };
        SparseMatrixd J;
        experimental::AssemblyPopulationCapture capturedPopulation;
        {
            ScopedPerformanceTimer timer("newton.jacobian");
            experimental::AssemblyPopulationCaptureScope populationScope(
                experimental::NewtonCaptureScope::hasCapacity() && cfg_.jacobian != "finite_difference"
                    ? &capturedPopulation : nullptr);
            J = (cfg_.jacobian == "finite_difference")
                ? assembler.finiteDifferenceJacobian(x, bcs, cfg_.finiteDifferenceStep)
                : assembler.assembleJacobian(x, bcs);
            if (cfg_.captureLinearSystem)
                detail::captureNewtonMatrix((linearCaptureFolder/"physical.bin").string(),J);
            addCarrierRowRegularization(J, N, cfg_.carrierRegularizationScale, bcs);
        }
        if (cfg_.captureLinearSystem) {
            detail::captureNewtonMatrix((linearCaptureFolder/"regularized.bin").string(),J);
            SparseMatrixd diagnosticMatrix=J; VectorXd diagnosticRhs=-r;
            if (cfg_.globalContinuityClosure.linearSystemStabilization=="replace_free_row_sum")
                stabilizeGlobalContinuityLinearSystem(diagnosticMatrix,diagnosticRhs,N,electronContactNodes,holeContactNodes);
            detail::captureNewtonMatrix((linearCaptureFolder/"closure_system.bin").string(),diagnosticMatrix);
            captureVectors("assembled.csv",x,r,diagnosticRhs);
        }
        VectorXd step;
        experimental::NewtonAuditScope::system(iter,x,r,J,activeRowWeights);
        experimental::NewtonCaptureScope::record([&](auto& capture) {
            if (capture.iterations.size() >= capture.limit) {
                capture.truncated = true;
                return;
            }
            experimental::NewtonCapturedIteration item;
            item.iteration = iter;
            item.x = x;
            item.residual = r;
            item.rowWeights = activeRowWeights;
            item.jacobian = J;
            item.population = std::move(capturedPopulation);
            capture.iterations.push_back(std::move(item));
        });
        try {
            const auto solveProductionSystem = [&linearSolver, this](
                const SparseMatrixd& matrix,
                const VectorXd& rhs) -> VectorXd {
                if (cfg_.linearEquilibration.mode != "l2_row_column")
                    return linearSolver.solve(matrix, rhs);
                ScopedPerformanceTimer timer("newton.linear_l2_row_column");
                const TwoSidedL2EquilibratedSystem equilibrated =
                    twoSidedL2Equilibrate(matrix, rhs);
                const VectorXd equilibratedCoordinates =
                    linearSolver.solve(equilibrated.matrix, equilibrated.rhs);
                return equilibratedCoordinates.cwiseProduct(
                    equilibrated.columnScale).eval();
            };
            if (cfg_.globalContinuityClosure.linearSystemStabilization ==
                "replace_free_row_sum") {
                ScopedPerformanceTimer timer(
                    "newton.global_continuity_linear_stabilization");
                step = solveGlobalContinuityStabilizedSystem(
                    linearSolver,
                    J,
                    -r,
                    N,
                    electronContactNodes,
                    holeContactNodes);
            } else if (cfg_.continuityRowScaling.enabled) {
                ScopedPerformanceTimer timer("newton.linear_row_scaling");
                const SparseMatrixd scaledJ = leftScaleRows(J, activeRowWeights);
                step = solveProductionSystem(
                    scaledJ, -r.cwiseProduct(activeRowWeights));
            } else {
                step = solveProductionSystem(J, -r);
            }
            if (cfg_.linearFactorRefinementSteps > 0) {
                ScopedPerformanceTimer timer("newton.guarded_factor_refinement");
                // Audit the actual equations, including the closure replacement.
                // The unmodified initial direction is still the usual Woodbury/direct solve.
                SparseMatrixd solvedMatrix=J;
                VectorXd solvedRhs=-r;
                if (cfg_.globalContinuityClosure.linearSystemStabilization=="replace_free_row_sum")
                    stabilizeGlobalContinuityLinearSystem(solvedMatrix,solvedRhs,N,electronContactNodes,holeContactNodes);
                const auto audit=slot_population::refineWithFactors(
                    solvedMatrix,solvedRhs,step,potentialScale,N,cfg_.linearFactorRefinementSteps,
                    [&](const VectorXd& rhs){return linearSolver.solve(solvedMatrix,rhs);});
                std::cout << "VELA_LINEAR_FACTOR_AUDIT " << nlohmann::json{
                    {"iteration",iter},{"dimension",solvedMatrix.rows()},
                    {"stabilization",cfg_.globalContinuityClosure.linearSystemStabilization},
                    {"corrections",audit.sweeps},{"before_mixed",audit.before.mixedRatio},
                    {"after_mixed",audit.after.mixedRatio},{"block_normwise",audit.after.blockNormwise},
                    {"constraint_error_V",audit.after.constraintV},{"passed",audit.after.accepted}
                }.dump() << std::endl;
                if (!audit.after.accepted)
                    throw std::runtime_error("guarded linear factor refinement exhausted");
            }
            if (cfg_.linearRefinementIterations > 0) {
                ScopedPerformanceTimer timer("newton.linear_refinement");
                const VectorXd weights = cfg_.continuityRowScaling.enabled
                    ? activeRowWeights : VectorXd::Ones(r.size());
                step = linear_refinement::refine(J, r, weights, std::move(step),
                    linearSolver, cfg_.linearRefinementIterations,
                    [](int, const VectorXd&, const linear_refinement::Defect&) {});
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
            if (cfg_.returnBestFeasibleOnFailure && hasBestFeasibleIterate) {
                result.iters = bestFeasibleIters;
                result.finalResidualNorm = bestFeasibleNorm;
                result.finalBlockNorms = blockResidualInfo(
                    bestFeasibleR, mesh_.numNodes());
                result.finalCarrierRowConvergence = bestFeasibleRowEval;
                result.finalGlobalContinuityClosure = bestFeasibleGlobalEval;
                result.solution = makeSolution(
                    assembler, bestFeasibleX, bestFeasibleIters);
            }
            return result;
        }

        captureVectors("linear.csv",x,r,step);
        if (shadowCapture) {
            SparseMatrixd auditedMatrix = J;
            VectorXd auditedRhs = -r;
            const bool closureReplaced =
                cfg_.globalContinuityClosure.linearSystemStabilization == "replace_free_row_sum";
            if (closureReplaced)
                stabilizeGlobalContinuityLinearSystem(auditedMatrix, auditedRhs, N,
                    electronContactNodes, holeContactNodes);
            else if (cfg_.continuityRowScaling.enabled) {
                auditedMatrix = leftScaleRows(auditedMatrix, activeRowWeights);
                auditedRhs = auditedRhs.cwiseProduct(activeRowWeights).eval();
            }
            shadowSolvedLinear = shadow_convergence::linearQuality(auditedMatrix, -auditedRhs, step);
            shadowSolvedLinear["matrix_identity"] =
                "returned-coordinate solve equations before exact Dirichlet update; not eliminated PARDISO CSR";
            shadowSolvedLinear["closure_replacement"] = closureReplaced;
            shadowSolvedLinear["continuity_row_scaling"] = !closureReplaced && cfg_.continuityRowScaling.enabled;
            shadowSolvedLinear["linear_equilibration"] = cfg_.linearEquilibration.mode;
            shadowSolvedLinear["configured_regularization"] = cfg_.carrierRegularizationScale;
        }
        experimental::NewtonAuditScope::direction(iter,step);
        experimental::NewtonCaptureScope::record([&](auto& capture) {
            if (capture.truncated || capture.iterations.empty()) return;
            auto& item = capture.iterations.back();
            item.linearDirection = step;
            item.hasLinearDirection = true;
        });

        // These rows are exact identity equations. A sparse direct solve can
        // pollute their zero update with tiny roundoff, creating artificial
        // merit decreases after all physical rows have reached their floor.
        if (cfg_.exactDirichletUpdates) {
            for (const auto& [node, value] : bcs.psi)
                step(static_cast<int>(node)) = -r(static_cast<int>(node));
            for (const auto& [node, value] : bcs.phin)
                step(N + static_cast<int>(node)) = -r(N + static_cast<int>(node));
            for (const auto& [node, value] : bcs.phip)
                step(2 * N + static_cast<int>(node)) = -r(2 * N + static_cast<int>(node));
        }
        const VectorXd rawStep = step;
        captureVectors("raw.csv",x,r,rawStep);
        if (shadowCapture) {
            shadowLinear = shadow_convergence::linearQuality(J, r, rawStep);
            shadowLinear["matrix_identity"] = "pre_closure_replacement J, post_constraint; configured regularization";
            shadowLinear["factorized_system_quality"] =
                shadow_convergence::unavailable("internal eliminated/equilibrated factorization coordinates not captured; see solved_equation_quality");
        }
        const VectorXd rawLinearResidual = J * rawStep + r;
        {
            ScopedPerformanceTimer timer("newton.update_constraints");
            capConfiguredStep(step, J, r);
        }
        if (cfg_.globalContinuityClosure.feasibilityRestoration ==
            "linear_qf_projection") {
            ScopedPerformanceTimer timer(
                "newton.global_continuity_feasibility_restoration");
            VectorXd restoredStep = step;
            if (restoreGlobalContinuityLinearizedStep(
                    restoredStep, J, r, N,
                    electronContactNodes, holeContactNodes)) {
                // The unconstrained rank-2 correction can be large when the
                // two terminal modes are weakly observable. Apply a configured
                // fraction and then reapply the physical QF bounds; this is an
                // inexact restoration that the nonlinear filter must accept.
                step += cfg_.globalContinuityClosure.feasibilityRestorationFactor *
                    (restoredStep - step);
                capConfiguredQfOnly(step);
            }
        }
        limitPopulationStep(step);
        captureVectors("processed.csv",x,r,step);
        VectorXd cappedLinearResidual = J * step + r;
        Real stepNorm = step.norm();
        const Real residualNormBeforeStep = residualNormFn(r);

        const ResidualBlockNormValue currentFilterBlocks =
            ResidualNorm::computeBlocks(r, mesh_.numNodes());
        ResidualBlockNormValue predictedFilterBlocks =
            ResidualNorm::computeBlocks(cappedLinearResidual, mesh_.numNodes());
        const auto closureResidualRatios = [&](const VectorXd& residual) {
            const auto sums = freeCarrierResidualSums(residual);
            return std::array<Real, 2>{
                std::abs(static_cast<Real>(sums[0])) /
                    std::max(activeGlobalElectronScale, Real{1.0e-300}),
                std::abs(static_cast<Real>(sums[1])) /
                    std::max(activeGlobalHoleScale, Real{1.0e-300})};
        };
        const std::array<Real, 2> currentClosureRatios = closureResidualRatios(r);
        std::array<Real, 2> predictedClosureRatios =
            closureResidualRatios(cappedLinearResidual);
        const auto blockFilterAccept = [&](const VectorXd& candidateResidual,
                                           Real alpha) {
            if (cfg_.lineSearchMode != "block_filter" &&
                cfg_.lineSearchMode != "block_closure_filter")
                return false;
            const ResidualBlockNormValue candidateBlocks =
                ResidualNorm::computeBlocks(candidateResidual, mesh_.numNodes());
            const Real current[3] = {
                currentFilterBlocks.psi,
                currentFilterBlocks.phin,
                currentFilterBlocks.phip};
            const Real predicted[3] = {
                predictedFilterBlocks.psi,
                predictedFilterBlocks.phin,
                predictedFilterBlocks.phip};
            const Real candidate[3] = {
                candidateBlocks.psi,
                candidateBlocks.phin,
                candidateBlocks.phip};
            const Real scales[3] = {
                residualScales.psi,
                residualScales.phin,
                residualScales.phip};
            const Real weights[3] = {
                residualWeights.psi,
                residualWeights.phin,
                residualWeights.phip};
            const bool absoluteBlocksEnforced =
                cfg_.blockAbsoluteConvergence.mode == "enforce";
            const Real ceilings[3] = {
                cfg_.blockAbsoluteConvergence.psiResidualCeiling,
                cfg_.blockAbsoluteConvergence.electronResidualCeiling,
                cfg_.blockAbsoluteConvergence.holeResidualCeiling};
            bool carrierImproved = false;
            bool anyImproved = false;
            bool anyActiveImproved = false;
            for (int block = 0; block < 3; ++block) {
                if (weights[block] <= 0.0)
                    continue;
                const Real floor = std::max(
                    cfg_.abstol * scales[block], Real{1.0e-300});
                // An already-qualified block is allowed to move inside its
                // authoritative absolute ceiling while another block is being
                // corrected.  Using only its often sub-ULP initial norm as the
                // envelope makes a Poisson boundary update impossible even
                // though the carrier block remains converged.
                Real envelopeBaseline = absoluteBlocksEnforced
                    ? std::max({scales[block], floor, ceilings[block]})
                    : std::max(scales[block], floor);
                if (cfg_.lineSearchMode == "block_closure_filter")
                    envelopeBaseline = std::max(envelopeBaseline, current[block]);
                const Real envelope = cfg_.residualFilterEnvelopeFactor *
                    envelopeBaseline;
                if (!std::isfinite(candidate[block]) ||
                    candidate[block] > envelope) {
                    return false;
                }
                // A trust-region-limited Newton direction predicts only a
                // fraction of the full residual reduction.  Scale the
                // sufficient-decrease target by that directional linear
                // prediction instead of assuming J*dx = -F after limiting.
                const Real predictedDecrease = std::max(
                    current[block] - predicted[block], Real{0.0});
                const Real required = current[block] -
                    cfg_.residualFilterGamma * alpha * predictedDecrease;
                const bool improved = candidate[block] <= required;
                anyImproved = anyImproved || improved;
                if (block > 0)
                    carrierImproved = carrierImproved || improved;
                const bool active = !absoluteBlocksEnforced ||
                    current[block] > ceilings[block];
                anyActiveImproved = anyActiveImproved || (active && improved);
            }
            if (absoluteBlocksEnforced)
                return anyActiveImproved;
            const bool carrierBlocksEnabled =
                residualWeights.phin > 0.0 || residualWeights.phip > 0.0;
            const bool deviceImproved =
                carrierBlocksEnabled ? carrierImproved : anyImproved;
            if (cfg_.lineSearchMode == "block_filter")
                return deviceImproved;

            const std::array<Real, 2> candidateClosureRatios =
                closureResidualRatios(candidateResidual);
            const Real absoluteClosureTarget =
                cfg_.globalContinuityClosure.absoluteTolerance /
                std::max(cfg_.globalContinuityClosure.sourceFloor, Real{1.0e-300});
            const Real closureTarget = std::max(
                cfg_.globalContinuityClosure.tolerance,
                absoluteClosureTarget);
            bool closureImproved = false;
            for (int carrier = 0; carrier < 2; ++carrier) {
                const Real envelope = cfg_.residualFilterEnvelopeFactor *
                    std::max(currentClosureRatios[carrier], closureTarget);
                if (!std::isfinite(candidateClosureRatios[carrier]) ||
                    candidateClosureRatios[carrier] > envelope) {
                    return false;
                }
                const Real predictedDecrease = std::max(
                    currentClosureRatios[carrier] -
                        predictedClosureRatios[carrier],
                    Real{0.0});
                const Real required = currentClosureRatios[carrier] -
                    cfg_.residualFilterGamma * alpha * predictedDecrease;
                closureImproved = closureImproved ||
                    candidateClosureRatios[carrier] <= required;
            }
            return deviceImproved || closureImproved;
        };

        const auto stableDecrease = [&](const VectorXd& trialResidual, Real) {
            if(cfg_.residualNorm=="l2")return stable_merit::compare(r,trialResidual).accepted;
            return stable_merit::compareBlocks(r,trialResidual,
                {residualScales.psi,residualScales.phin,residualScales.phip},
                {residualWeights.psi,residualWeights.phin,residualWeights.phip}).accepted;
        };
        const bool localRowCorrectionEligible =
            cfg_.lineSearchMode == "block_filter" &&
            cfg_.blockAbsoluteConvergence.mode == "enforce" &&
            cfg_.carrierRowConvergence.mode == "enforce" &&
            blockAbsoluteConvergenceSatisfied(r);
        const NewtonCarrierRowConvergenceEvaluation currentLocalRows =
            localRowCorrectionEligible ? carrierRowEval(x)
                                       : NewtonCarrierRowConvergenceEvaluation{};
        const auto runLineSearch = [&](const VectorXd& trialStep) {
            const VectorXd baseLow = splitRuntime ? splitRuntime->low() : VectorXd{};
            const auto populationFrame=cfg_.populationStepFraction>0.
                ? std::optional<population_step::Frame>(population_step::frame(assembler,x,trialStep,bcs))
                : std::nullopt;
            const auto candidateState = [&](Real alpha) -> VectorXd {
                if(populationFrame) {
                    const VectorXd mapped=population_step::finiteStep(assembler,*populationFrame,trialStep,alpha);
                    if(splitRuntime)return splitRuntime->candidate(x,mapped,1.,baseLow);
                    return x+mapped;
                }
                if (splitRuntime) return splitRuntime->candidate(x, trialStep, alpha, baseLow);
                return x + alpha * trialStep;
            };
            const auto decreaseAccept = [&](const VectorXd& candidateResidual,
                                            Real alpha) {
                // Profiling-only evidence for a full Newton step rejected
                // while the current absolute blocks are already qualified.
                // Observe both gates before deciding; never change acceptance.
                if (activePerformanceProfiler() && localRowCorrectionEligible &&
                    !currentLocalRows.satisfied && alpha == 1.0) {
                    const auto blocks = ResidualNorm::computeBlocks(
                        candidateResidual, mesh_.numNodes());
                    const auto rows = carrierRowEval(candidateState(1.0));
                    const std::string prefix = "newton.local_rows.full_step." +
                        std::to_string(iter) + ".";
                    observePerformanceValue(prefix + "psi", blocks.psi);
                    observePerformanceValue(prefix + "electron", blocks.phin);
                    observePerformanceValue(prefix + "hole", blocks.phip);
                    observePerformanceValue(prefix + "row_ratio", rows.maxRatio);
                    observePerformanceValue(prefix + "current_row_ratio", currentLocalRows.maxRatio);
                    observePerformanceValue(prefix + "block_qualified",
                        blockAbsoluteConvergenceSatisfied(candidateResidual) ? 1.0 : 0.0);
                }
                if (blockFilterAccept(candidateResidual, alpha))
                    return true;
                // Global blocks can already be qualified while a depleted
                // carrier row still needs correction. Preserve every absolute
                // ceiling and require local progress in that otherwise inactive
                // filter; reporting-only rows do not change the search.
                if (!localRowCorrectionEligible || currentLocalRows.satisfied ||
                    !blockAbsoluteConvergenceSatisfied(candidateResidual)) {
                    return false;
                }
                const VectorXd candidate = candidateState(alpha);
                if (!candidate.allFinite())
                    return false;
                const auto candidateRows = carrierRowEval(candidate);
                return std::isfinite(candidateRows.maxRatio) &&
                    candidateRows.maxRatio < currentLocalRows.maxRatio &&
                    candidateRows.maxRatio <=
                        (1.0 - cfg_.residualFilterGamma * alpha) *
                            currentLocalRows.maxRatio;
            };
            auto searched = cfg_.bankRose.enabled ? bankRose.search(
                iter, x, trialStep, r,
                [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
                [&](const VectorXd& candidate, const VectorXd&) {
                    return assembler.hasPositiveFiniteCarriers(candidate);
                }, residualNormFn,
                (splitRuntime||populationFrame) ? BacktrackingLineSearch::CandidateFunction(
                    [&](const VectorXd&, const VectorXd&, Real alpha) { return candidateState(alpha); })
                    : BacktrackingLineSearch::CandidateFunction{},
                cfg_.diagnostics) : lineSearch.search(
                x, trialStep, r,
                [&](const VectorXd& candidate) {
                    return assembler.residual(candidate, bcs);
                },
                [&](const VectorXd& candidate, const VectorXd&) {
                    return assembler.hasPositiveFiniteCarriers(candidate);
                },
                globalClosureLineSearchNorm,
                (cfg_.lineSearchMode == "block_filter" ||
                 cfg_.lineSearchMode == "block_closure_filter")
                    ? BacktrackingLineSearch::DecreaseAcceptFunction(decreaseAccept)
                    : (cfg_.stableMeritComparison
                        ? BacktrackingLineSearch::DecreaseAcceptFunction(stableDecrease)
                        : BacktrackingLineSearch::DecreaseAcceptFunction{}),
                (splitRuntime||populationFrame)?BacktrackingLineSearch::CandidateFunction([&](const VectorXd&,const VectorXd&,Real alpha){return candidateState(alpha);}):BacktrackingLineSearch::CandidateFunction{});
            if(splitRuntime) {
                if(searched.accepted)candidateState(searched.damping);
                else splitRuntime->setLow(baseLow);
            }
            return searched;
        };

        const auto bankRoseBeforeSearch = bankRose;
        const VectorXd auditBaseLow=(splitRuntime && shouldWriteLocalUpdateTrace(iter))?splitRuntime->low():VectorXd{};
        auto ls = runLineSearch(step);
        bool trustRegionContractedThisIteration = false;
        while (!ls.accepted &&
               cfg_.quasiFermiUpdateLimitMode == "uniform_trust_region" &&
               cfg_.quasiFermiTrustRegionShrinkFactor < 1.0 &&
               quasiFermiTrustRegionMultiplier >
                   cfg_.quasiFermiTrustRegionMinMultiplier * (1.0 + 1.0e-12)) {
            const Real contractedMultiplier = std::max(
                cfg_.quasiFermiTrustRegionMinMultiplier,
                quasiFermiTrustRegionMultiplier *
                    cfg_.quasiFermiTrustRegionShrinkFactor);
            if (!(contractedMultiplier < quasiFermiTrustRegionMultiplier))
                break;
            quasiFermiTrustRegionMultiplier = contractedMultiplier;
            trustRegionContractedThisIteration = true;
            step = rawStep;
            {
                ScopedPerformanceTimer timer("newton.update_constraints_retry");
                capConfiguredStep(step, J, r);
            }
            if (cfg_.globalContinuityClosure.feasibilityRestoration ==
                "linear_qf_projection") {
                VectorXd restoredStep = step;
                if (restoreGlobalContinuityLinearizedStep(
                        restoredStep, J, r, N, electronContactNodes,
                        holeContactNodes)) {
                    step +=
                        cfg_.globalContinuityClosure.feasibilityRestorationFactor *
                        (restoredStep - step);
                    capConfiguredQfOnly(step);
                }
            }
            limitPopulationStep(step);
            cappedLinearResidual = J * step + r;
            predictedFilterBlocks = ResidualNorm::computeBlocks(
                cappedLinearResidual, mesh_.numNodes());
            stepNorm = step.norm();
            ls = runLineSearch(step);
        }

        const bool trustRegionBoundaryActive =
            cfg_.quasiFermiUpdateLimitMode == "uniform_trust_region" &&
            stepNorm < rawStep.norm() * (1.0 - 1.0e-12);
        const Real predictedResidualNorm = residualNormFn(cappedLinearResidual);

        // A small damped correction can round to zero while other rows keep
        // making progress, producing a long line-search tail. Repartition the
        // same physical state, then rebuild residuals and retry the direction.
        // This is not an accepted Newton update and never accepts a loose gate.
        bool lostCarrierCorrection = false;
        if (cfg_.quasiFermiRecenterOnStall && quasiFermiRecenters < 8 &&
            cfg_.blockAbsoluteConvergence.mode == "enforce" &&
            currentFilterBlocks.psi <= cfg_.blockAbsoluteConvergence.psiResidualCeiling) {
            const Real alpha = ls.accepted ? ls.damping : 1.0;
            for (int block = 1; block <= 2 && !lostCarrierCorrection; ++block) {
                const auto& boundaries = block == 1 ? bcs.phin : bcs.phip;
                const Real ceiling = block == 1
                    ? cfg_.blockAbsoluteConvergence.electronResidualCeiling
                    : cfg_.blockAbsoluteConvergence.holeResidualCeiling;
                for (int i = 0; i < N; ++i) {
                    const int row = block * N + i;
                    const Real delta = alpha * step(row);
                    if (boundaries.count(static_cast<Index>(i)) == 0 &&
                        std::abs(r(row)) > ceiling / std::sqrt(static_cast<Real>(N)) &&
                        delta != 0.0 && x(row) + delta == x(row)) {
                        lostCarrierCorrection = true;
                        break;
                    }
                }
            }
        }
        if (lostCarrierCorrection) {
            // This search is discarded, not an accepted nonlinear update.
            bankRose = bankRoseBeforeSearch;
            ScopedPerformanceTimer recenterTimer("newton.qf_recenter");
            const bool auditRecenter = cfg_.localUpdateDiagnostics.enabled;
            const VectorXd beforeElectron = auditRecenter ? assembler.electronDensity(x) : VectorXd{};
            const VectorXd beforeHole = auditRecenter ? assembler.holeDensity(x) : VectorXd{};
            const VectorXd beforeResidual = auditRecenter ? r : VectorXd{};
            assembler.recenterQuasiFermiState(x);
            continuityTerms.clear();
            r = assembler.residual(x, bcs);
            if (auditRecenter) {
                const auto densityDifference = [](const VectorXd& before, const VectorXd& after) {
                    return ((after - before).cwiseAbs().array() /
                        before.cwiseAbs().array().max(1.0)).maxCoeff();
                };
                observePerformanceValue("newton.qf_recenter.electron_relative_change",
                    densityDifference(beforeElectron, assembler.electronDensity(x)));
                observePerformanceValue("newton.qf_recenter.hole_relative_change",
                    densityDifference(beforeHole, assembler.holeDensity(x)));
                observePerformanceValue("newton.qf_recenter.residual_absolute_change",
                    (r - beforeResidual).lpNorm<Eigen::Infinity>());
                observePerformanceValue("newton.qf_recenter.jacobian_relative_change",
                    (assembler.assembleJacobian(x, bcs) - J).norm() / std::max(1.0, J.norm()));
            }
            acceptedX = x;
            acceptedR = r;
            activeGlobalEval = globalClosureEval(x);
            ++quasiFermiRecenters;
            incrementPerformanceCounter("newton.qf_recenters");
            incrementPerformanceCounter("newton.qf_recenter_discarded_line_search_trials", ls.attempts);
            // Numeric factorization is performed again on the next iteration;
            // topology and physical bestSolution/history remain valid.
            --iter;
            continue;
        }
        writeLocalUpdateTraceCsv(
            iter, x, r, J, rawStep, step,
            rawLinearResidual, cappedLinearResidual, ls,auditBaseLow);
        experimental::NewtonAuditScope::search(iter,step,ls);
        experimental::NewtonCaptureScope::record([&](auto& capture) {
            if (capture.truncated || capture.iterations.empty()) return;
            auto& item = capture.iterations.back();
            item.protectedDirection = step;
            item.searched = ls;
            item.hasSearch = true;
        });

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
            rejectedTrace.contactMajorityQfDrop =
                contactMajorityQfBranchDrop(acceptedX, acceptedIters);
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
            const bool stalledBlockAbsoluteConverged =
                cfg_.blockAbsoluteConvergence.mode == "enforce" &&
                blockAbsoluteConvergenceSatisfied(acceptedR);
            if (stalledBlockAbsoluteConverged &&
                carrierRowsAcceptConvergence(stalledRowEval) &&
                globalClosureAcceptsConvergence(stalledGlobalEval) &&
                contactMajorityQfBranchAccepts(acceptedX, acceptedIters)) {
                finishConverged(
                    "block_abstol_line_search_stall", acceptedX, acceptedR,
                    acceptedIters, stalledNorm, stalledRowEval,
                    stalledGlobalEval);
                return result;
            }
            if (stalledNorm <= stallResidualFloor &&
                blockAbsoluteConvergenceSatisfied(acceptedR) &&
                carrierRowsAcceptConvergence(stalledRowEval) &&
                globalClosureAcceptsConvergence(stalledGlobalEval) &&
                contactMajorityQfBranchAccepts(acceptedX, acceptedIters) &&
                contactMajorityQfAcceptsFloor(acceptedX, acceptedIters)) {
                finishConverged("stall_residual_floor", acceptedX, acceptedR,
                                acceptedIters, stalledNorm, stalledRowEval,
                                stalledGlobalEval);
                return result;
            }
            const DDSolution stalledSolution = makeSolution(assembler, acceptedX, acceptedIters);
            const Real stalledContactMajorityQfDrop = maxContactMajorityQuasiFermiDrop(stalledSolution);
            Real bestRejectedContactMajorityQfDrop = std::numeric_limits<Real>::infinity();
            if (ls.bestRejectedCandidate) {
                const VectorXd oldLow=splitRuntime?splitRuntime->low():VectorXd{};
                if(splitRuntime)splitRuntime->candidate(x,step,ls.bestRejectedDamping,oldLow);
                const DDSolution bestRejectedSolution = makeSolution(assembler, ls.bestRejectedX, acceptedIters);
                bestRejectedContactMajorityQfDrop = maxContactMajorityQuasiFermiDrop(bestRejectedSolution);
                if(splitRuntime)splitRuntime->setLow(oldLow);
            }

            if (isPoissonLineSearchStall(ls, stalledBlocks, stalledNorm, stalledContactMajorityQfDrop, cfg_) &&
                blockAbsoluteConvergenceSatisfied(acceptedR) &&
                carrierRowsAcceptConvergence(stalledRowEval) &&
                globalClosureAcceptsConvergence(stalledGlobalEval) &&
                contactMajorityQfBranchAccepts(acceptedX, acceptedIters)) {
                finishConverged("poisson_line_search_stall_floor", acceptedX, acceptedR,
                                acceptedIters, stalledNorm, stalledRowEval,
                                stalledGlobalEval);
                return result;
            }
            if (isCarrierRowQualifiedLineSearchStall(
                    ls, stalledBlocks, stalledNorm,
                    stalledContactMajorityQfDrop, stalledRowEval, cfg_) &&
                blockAbsoluteConvergenceSatisfied(acceptedR) &&
                globalClosureAcceptsConvergence(stalledGlobalEval) &&
                contactMajorityQfBranchAccepts(acceptedX, acceptedIters)) {
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
            const bool stalledBranchGuardSatisfied =
                contactMajorityQfBranchAccepts(acceptedX, acceptedIters);
            const std::string rejectedFailureReason =
                (stalledGlobalEval.enforced && !stalledGlobalEval.satisfied)
                    ? std::string("global_continuity_closure_line_search_rejected")
                    : (stalledRowEval.enforced && !stalledRowEval.satisfied)
                    ? std::string("carrier_row_convergence_line_search_rejected")
                    : (!stalledBranchGuardSatisfied)
                    ? std::string("contact_majority_qf_branch_guard_line_search_rejected")
                    : (cfg_.blockAbsoluteConvergence.mode == "enforce" &&
                       !blockAbsoluteConvergenceSatisfied(acceptedR))
                    ? std::string("block_absolute_convergence_line_search_rejected")
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
            if (cfg_.returnBestFeasibleOnFailure && hasBestFeasibleIterate) {
                result.iters = bestFeasibleIters;
                result.finalResidualNorm = bestFeasibleNorm;
                result.finalBlockNorms = blockResidualInfo(
                    bestFeasibleR, mesh_.numNodes());
                result.finalCarrierRowConvergence = bestFeasibleRowEval;
                result.finalGlobalContinuityClosure = bestFeasibleGlobalEval;
                result.solution = makeSolution(
                    assembler, bestFeasibleX, bestFeasibleIters);
            }
            return result;
        }

        auditReplacementRows(iter, x, r, ls.x, ls.residual);
        captureVectors("selected.csv",ls.x,ls.residual,VectorXd(ls.x-x));
        const VectorXd progressDelta = cfg_.progressDirectory.empty() ? VectorXd{} : VectorXd(ls.x-x);
        x = ls.x;
        r = ls.residual;
        acceptedX = x;
        acceptedR = r;
        acceptedIters = iter;

        // Record the norm of the actually applied update (damped step) so that
        // per-iteration metrics are consistent with the accepted solution.
        const Real appliedStepNorm = ls.damping * stepNorm;
        const Real residualNorm = residualNormFn(r);
        updateBestAcceptedIterate(x, r, iter, residualNorm);
        if (cfg_.quasiFermiUpdateLimitMode == "uniform_trust_region" &&
            ls.damping < 1.0 - 1.0e-12 &&
            cfg_.quasiFermiTrustRegionShrinkFactor < 1.0) {
            quasiFermiTrustRegionMultiplier = std::max(
                cfg_.quasiFermiTrustRegionMinMultiplier,
                quasiFermiTrustRegionMultiplier *
                    cfg_.quasiFermiTrustRegionShrinkFactor);
        } else if (!trustRegionContractedThisIteration &&
            trustRegionBoundaryActive && ls.damping >= 1.0 - 1.0e-12 &&
            cfg_.quasiFermiTrustRegionGrowthFactor > 1.0 &&
            quasiFermiTrustRegionMultiplier <
                cfg_.quasiFermiTrustRegionMaxMultiplier) {
            const Real predictedDecrease =
                residualNormBeforeStep - predictedResidualNorm;
            const Real actualDecrease = residualNormBeforeStep - residualNorm;
            if (predictedDecrease >
                    std::numeric_limits<Real>::epsilon() *
                        std::max(residualNormBeforeStep, Real{1.0}) &&
                actualDecrease / predictedDecrease >=
                    cfg_.quasiFermiTrustRegionExpansionThreshold) {
                quasiFermiTrustRegionMultiplier = std::min(
                    cfg_.quasiFermiTrustRegionMaxMultiplier,
                    quasiFermiTrustRegionMultiplier *
                        cfg_.quasiFermiTrustRegionGrowthFactor);
            }
        }
        NewtonIterationInfo info;
        info.iter = iter;
        info.residualNorm = residualNorm;
        info.stepNorm = appliedStepNorm;
        info.dampingFactor = ls.damping;
        info.relativeResidualNorm = ResidualNorm::relative(residualNorm, initialNorm);
        info.contactMajorityQfDrop =
            contactMajorityQfBranchDrop(x, iter);
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
        const NewtonGlobalContinuityClosureEvaluation globalEval =
            globalClosureEval(x);
        activeGlobalEval = globalEval;
        writeProgress(result.history.back(), globalEval, progressDelta,rawStep,step);
        if (carrierRowsAcceptConvergence(
                result.history.back().carrierRowConvergence) &&
            globalClosureAcceptsConvergence(globalEval) &&
            (!hasBestFeasibleIterate || residualNorm < bestFeasibleNorm)) {
            hasBestFeasibleIterate = true;
            bestFeasibleX = x;
            bestFeasibleR = r;
            bestFeasibleIters = iter;
            bestFeasibleNorm = residualNorm;
            bestFeasibleRowEval = result.history.back().carrierRowConvergence;
            bestFeasibleGlobalEval = globalEval;
        }
        if (cfg_.verbose) {
            const NewtonBlockResidualInfo& blocks = result.history.back().blockResiduals;
            long double electronFreeResidualSum = 0.0L;
            long double holeFreeResidualSum = 0.0L;
            for (int i = 0; i < N; ++i) {
                const Index node = static_cast<Index>(i);
                if (bcs.phin.find(node) == bcs.phin.end())
                    electronFreeResidualSum +=
                        static_cast<long double>(r(N + i));
                if (bcs.phip.find(node) == bcs.phip.end())
                    holeFreeResidualSum +=
                        static_cast<long double>(r(2 * N + i));
            }
            const Real electronClosureScale = std::max(
                {std::abs(globalEval.electron.contactFlux),
                 std::abs(globalEval.electron.integratedSource),
                 cfg_.globalContinuityClosure.sourceFloor});
            const Real holeClosureScale = std::max(
                {std::abs(globalEval.hole.contactFlux),
                 std::abs(globalEval.hole.integratedSource),
                 cfg_.globalContinuityClosure.sourceFloor});
            emitVerboseLine(
                "Newton iter " + std::to_string(iter) +
                " residual=" + std::to_string(residualNorm) +
                " step=" + std::to_string(appliedStepNorm) +
                " damping=" + std::to_string(ls.damping) +
                " blocks=(" + std::to_string(blocks.psi) + "," +
                std::to_string(blocks.phin) + "," +
                std::to_string(blocks.phip) + ")" +
                " global_closure=(e=" +
                std::to_string(globalEval.electron.ratio) +
                ",h=" + std::to_string(globalEval.hole.ratio) +
                ",qualified=" +
                std::to_string(globalEval.electron.qualified ? 1 : 0) +
                "/" + std::to_string(globalEval.hole.qualified ? 1 : 0) +
                ",free_residual_ratio=" +
                std::to_string(std::abs(static_cast<Real>(electronFreeResidualSum)) /
                               electronClosureScale) +
                "/" +
                std::to_string(std::abs(static_cast<Real>(holeFreeResidualSum)) /
                               holeClosureScale) +
                ")");
        }

        const Real rel = result.history.back().relativeResidualNorm;
        const bool absoluteConverged = residualNorm <= cfg_.abstol;
        const bool relativeConverged = rel <= cfg_.reltol;
        const bool blockAbsoluteConverged =
            cfg_.blockAbsoluteConvergence.mode == "enforce" &&
            blockAbsoluteConvergenceSatisfied(r);
        const bool primaryToleranceConverged =
            blockAbsoluteConverged ||
            (cfg_.blockAbsoluteConvergence.mode == "off" &&
             (absoluteConverged || relativeConverged));
        const NewtonCarrierRowConvergenceEvaluation& rowEval =
            result.history.back().carrierRowConvergence;
        if (primaryToleranceConverged &&
            carrierRowsAcceptConvergence(rowEval) &&
            globalClosureAcceptsConvergence(globalEval) &&
            contactMajorityQfBranchAccepts(x, iter)) {
            finishConverged(
                blockAbsoluteConverged
                    ? "block_abstol"
                    : (absoluteConverged ? "abstol" : "reltol"),
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
            blockAbsoluteConvergenceSatisfied(r) &&
            contactMajorityQfBranchAccepts(x, iter) &&
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
        if (cfg_.diagnosticStopPolicy == "row_trend_v1" &&
            (!pendingTrendStop.empty() || iter>=cfg_.maxIter)) {
            // All original success checks above have priority; snapshots already saved.
            result.converged=false;
            result.solution=makeSolution(assembler,x,acceptedIters);
            result.iters=acceptedIters;
            result.finalResidualNorm=residualNorm;
            result.finalBlockNorms=result.history.back().blockResiduals;
            result.finalCarrierRowConvergence=rowEval;
            result.finalGlobalContinuityClosure=globalEval;
            result.failureDiagnostics.failureReason=pendingTrendStop.empty()
                ? "iteration_limit" : pendingTrendStop;
            return result;
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
    const Real finalContactMajorityQfDrop =
        maxContactMajorityQuasiFermiDrop(result.solution);
    const bool finalContactMajorityQfAcceptsFloor =
        cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V <= 0.0 ||
        (std::isfinite(finalContactMajorityQfDrop) &&
         finalContactMajorityQfDrop <=
             cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V);
    if (result.finalResidualNorm <= stallResidualFloor &&
        blockAbsoluteConvergenceSatisfied(acceptedR) &&
        carrierRowsAcceptConvergence(finalRowEval) &&
        globalClosureAcceptsConvergence(finalGlobalEval) &&
        contactMajorityQfBranchAccepts(acceptedX, acceptedIters) &&
        finalContactMajorityQfAcceptsFloor) {
        finishConverged("max_iter_stall_residual_floor", acceptedX, acceptedR,
                        acceptedIters, result.finalResidualNorm, finalRowEval,
                        finalGlobalEval);
        return result;
    }

    // Closure is a property of the converged device state, not a required
    // invariant of its Newton path.  When explicitly requested, apply one
    // bounded minimum-norm QF restoration only after the device residual has
    // reached its accepted numerical floor.  Nonlinear backtracking accepts
    // the correction only when closure, carrier rows, and the device floor all
    // pass simultaneously.
    if (cfg_.globalContinuityClosure.feasibilityRestoration ==
            "final_linear_qf_projection" &&
        result.finalResidualNorm <= stallResidualFloor &&
        carrierRowsAcceptConvergence(finalRowEval) &&
        finalGlobalEval.enforced && !finalGlobalEval.satisfied) {
        SparseMatrixd finalJacobian;
        {
            ScopedPerformanceTimer timer(
                "newton.final_global_continuity_restoration_jacobian");
            finalJacobian = (cfg_.jacobian == "finite_difference")
                ? assembler.finiteDifferenceJacobian(
                    acceptedX, bcs, cfg_.finiteDifferenceStep)
                : assembler.assembleJacobian(acceptedX, bcs);
            if (cfg_.globalContinuityClosure.restorationJacobian == "physical") {
                if (cfg_.jacobian != "finite_difference")
                    finalJacobian = assembler.assembleJacobian(acceptedX, bcs, false);
            } else {
                addCarrierRowRegularization(finalJacobian, N, cfg_.carrierRegularizationScale, bcs);
            }
        }
        SparseMatrixd physicalJacobian;
        if (closureAuditEnabled)
            physicalJacobian = cfg_.jacobian == "finite_difference"
                ? assembler.finiteDifferenceJacobian(acceptedX, bcs, cfg_.finiteDifferenceStep)
                : assembler.assembleJacobian(acceptedX, bcs, false);
        auto stepNorms = [](const VectorXd& value) {
            return nlohmann::json{{"l2", value.norm()}, {"linf", value.lpNorm<Eigen::Infinity>()}};
        };
        auto closureJson = [](const NewtonGlobalContinuityClosureEvaluation& value) {
            auto carrier = [](const NewtonGlobalContinuityCarrierClosure& c) {
                return nlohmann::json{{"mismatch", c.mismatch}, {"mode", c.qualified ? "relative" : "absolute"},
                    {"satisfied", c.satisfied}, {"contact_flux", c.contactFlux},
                    {"integrated_source", c.integratedSource}, {"ratio", c.ratio}};
            };
            return nlohmann::json{{"electron", carrier(value.electron)}, {"hole", carrier(value.hole)}};
        };
        auto blocksJson = [&](const VectorXd& value) {
            const auto blocks = ResidualNorm::computeBlocks(value, mesh_.numNodes());
            return nlohmann::json{{"psi", blocks.psi}, {"electron", blocks.phin}, {"hole", blocks.phip}};
        };
        const Real idealAlphaLower = globalContinuityIdealAlphaLowerBound(
            finalGlobalEval, cfg_.globalContinuityClosure);
        VectorXd correction = VectorXd::Zero(acceptedX.size());
        nlohmann::json projection;
        std::vector<CoupledDDCarrierTermDiagnostic> physicalRows;
        if (cfg_.globalContinuityClosure.restorationPhysicalRowsOnly)
            physicalRows = assembler.carrierContinuityEquationTermDiagnostics(acceptedX, {});
        const bool projected = restoreGlobalContinuityLinearizedStep(
            correction, finalJacobian, acceptedR, N, electronContactNodes, holeContactNodes,
            closureAuditEnabled ? &projection : nullptr,
            cfg_.globalContinuityClosure.restorationPhysicalRowsOnly ? &physicalRows : nullptr);
        if (closureAuditEnabled)
            writeClosureAudit({{"event", "final_restoration_trigger"}, {"iteration", acceptedIters},
                {"initial_closure", closureJson(finalGlobalEval)}, {"projection", projection},
                {"raw_step", stepNorms(VectorXd::Zero(acceptedX.size()))},
                {"projected_step", stepNorms(correction)}, {"initial_blocks", blocksJson(acceptedR)},
                {"initial_scalar_residual", result.finalResidualNorm},
                {"initial_carrier_rows_satisfied", carrierRowsAcceptConvergence(finalRowEval)},
                {"ideal_alpha_lower_bound", idealAlphaLower},
                {"carrier_diagonal_floor_enabled", cfg_.carrierDiagonalFloor.enabled},
                {"carrier_diagonal_floor_scale", cfg_.carrierDiagonalFloor.scale},
                {"carrier_regularization_scale", cfg_.carrierRegularizationScale},
                {"rejection_reason", projected ? "none" : "projection_failed"}});
        if (projected) {
            const VectorXd projectedCorrection = correction;
            capConfiguredQfOnly(correction);
            nlohmann::json cappedDofs = nlohmann::json::array();
            Real maxOverLimit = 0;
            if (closureAuditEnabled) {
                for (int row = N; row < 3 * N; ++row) {
                    if (correction(row) == projectedCorrection(row)) continue;
                    const Real over = std::abs(projectedCorrection(row) - correction(row)) * potentialScale;
                    maxOverLimit = std::max(maxOverLimit, over);
                    cappedDofs.push_back({{"node_id", row % N}, {"carrier", row < 2 * N ? "electron" : "hole"},
                        {"over_limit_V", over}});
                }
            }
            const Real correctionNorm = correction.norm();
            Real damping =
                cfg_.globalContinuityClosure.feasibilityRestorationFactor;
            const VectorXd restorationLow=splitRuntime?splitRuntime->low():VectorXd{};
            for (int trial = 0; trial < 13; ++trial, damping *= 0.5) {
                const VectorXd candidateX = splitRuntime
                    ? splitRuntime->candidate(acceptedX,correction,damping,restorationLow)
                    : VectorXd(acceptedX + damping * correction);
                const VectorXd candidateR = assembler.residual(
                    candidateX, bcs);
                const Real candidateNorm = residualNormFn(candidateR);
                // Audit every trial, even those rejected by the scalar gate.
                if (!closureAuditEnabled && (!std::isfinite(candidateNorm) || candidateNorm > stallResidualFloor))
                    continue;
                const NewtonCarrierRowConvergenceEvaluation candidateRowEval = std::isfinite(candidateNorm)
                    ? carrierRowEval(candidateX) : NewtonCarrierRowConvergenceEvaluation{};
                const NewtonGlobalContinuityClosureEvaluation candidateGlobalEval = std::isfinite(candidateNorm)
                    ? globalClosureEval(candidateX) : NewtonGlobalContinuityClosureEvaluation{};
                const std::string rejection = !std::isfinite(candidateNorm) ? "nonfinite_residual"
                    : candidateNorm > stallResidualFloor ? "device_floor"
                    : !blockAbsoluteConvergenceSatisfied(candidateR) ? "block_absolute_convergence"
                    : !carrierRowsAcceptConvergence(candidateRowEval) ? "carrier_rows"
                    : !globalClosureAcceptsConvergence(candidateGlobalEval) ? "global_closure"
                    : !contactMajorityQfBranchAccepts(candidateX, acceptedIters) ? "contact_qf_branch"
                    : !contactMajorityQfAcceptsFloor(candidateX, acceptedIters) ? "contact_qf_floor"
                    : "accepted";
                if (closureAuditEnabled) {
                    const VectorXd applied = candidateX - acceptedX;
                    const VectorXd nominal = damping * correction;
                    std::vector<int> lostRows;
                    for (int row = N; row < 3 * N; ++row)
                        if (nominal(row) != 0.0 && applied(row) == 0.0) lostRows.push_back(row);
                    std::sort(lostRows.begin(), lostRows.end(), [&](int left, int right) {
                        if (std::abs(nominal(left)) != std::abs(nominal(right)))
                            return std::abs(nominal(left)) > std::abs(nominal(right));
                        return left < right;
                    });
                    nlohmann::json largestLost = nlohmann::json::array();
                    for (std::size_t i = 0; i < std::min<std::size_t>(8, lostRows.size()); ++i) {
                        const int row = lostRows[i];
                        largestLost.push_back({{"node_id", row % N},
                            {"carrier", row < 2 * N ? "electron" : "hole"},
                            {"state", acceptedX(row)}, {"nominal_delta", nominal(row)},
                            {"directional_ulp", std::abs(std::nextafter(acceptedX(row),
                                std::copysign(std::numeric_limits<Real>::infinity(), nominal(row))) - acceptedX(row))}});
                    }
                    auto predict = [&](const SparseMatrixd& matrix, const VectorXd& direction) {
                        return globalContinuityLinearizedSums(matrix, acceptedR, direction, N,
                            electronContactNodes, holeContactNodes);
                    };
                    const auto zeroSums = predict(finalJacobian, VectorXd::Zero(acceptedX.size()));
                    const auto physicalSums = predict(physicalJacobian, applied);
                    const std::array<Real, 2> predictedMismatch{
                        finalGlobalEval.electron.mismatch - (physicalSums[0] - zeroSums[0]),
                        finalGlobalEval.hole.mismatch - (physicalSums[1] - zeroSums[1])};
                    writeClosureAudit({{"event", "final_restoration_trial"}, {"iteration", acceptedIters},
                        {"trial", trial}, {"alpha", damping}, {"ideal_alpha_lower_bound", idealAlphaLower},
                        {"ideal_linear_infeasible", damping < idealAlphaLower},
                        {"initial_free_sums", zeroSums}, {"predicted_projected_free_sums", predict(finalJacobian, projectedCorrection)},
                        {"predicted_capped_free_sums", predict(finalJacobian, correction)},
                        {"predicted_nominal_free_sums", predict(finalJacobian, damping * correction)},
                        {"predicted_applied_free_sums", predict(finalJacobian, applied)},
                        {"physical_predicted_applied_free_sums", physicalSums},
                        {"free_sum_inferred_mismatch", predictedMismatch},
                        {"contact_source_vs_free_sum_gap", {
                            candidateGlobalEval.electron.mismatch + freeCarrierResidualSums(candidateR)[0],
                            candidateGlobalEval.hole.mismatch + freeCarrierResidualSums(candidateR)[1]}},
                        {"actual_free_sums", freeCarrierResidualSums(candidateR)},
                        {"actual_closure", closureJson(candidateGlobalEval)},
                        {"scalar_residual", candidateNorm}, {"blocks", blocksJson(candidateR)},
                        {"carrier_rows_satisfied", carrierRowsAcceptConvergence(candidateRowEval)},
                        {"carrier_rows_qualified", candidateRowEval.qualifiedRowCount},
                        {"carrier_rows_max_ratio", candidateRowEval.maxRatio},
                        {"projected_step", stepNorms(projectedCorrection)}, {"capped_step", stepNorms(correction)},
                        {"nominal_step", stepNorms(nominal)}, {"applied_step", stepNorms(applied)},
                        {"electron_nominal_step", stepNorms(nominal.segment(N, N))},
                        {"electron_applied_step", stepNorms(applied.segment(N, N))},
                        {"hole_nominal_step", stepNorms(nominal.segment(2 * N, N))},
                        {"hole_applied_step", stepNorms(applied.segment(2 * N, N))},
                        {"addition_rounding_error", stepNorms(applied - nominal)},
                        {"qf_addition_lost_dofs", lostRows.size()}, {"largest_lost_updates", largestLost},
                        {"capped_dofs", cappedDofs}, {"capped_dof_count", cappedDofs.size()},
                        {"max_over_limit_V", maxOverLimit}, {"poisson_recorrect_applied", false},
                        {"psi_step_before_recorrect", 0.0}, {"psi_step_after_recorrect", 0.0},
                        {"rejection_reason", rejection}});
                }
                if (trial == 0 && closureAuditEnabled && cfg_.globalContinuityClosure.roundoffAudit &&
                    std::isfinite(candidateNorm)) {
                    auto audit = auditGlobalContinuityRoundoff(assembler, candidateX, bcs,
                        cfg_.globalContinuityClosure, mesh_.numEdges());
                    audit["stage"] = "final_trial_alpha1";
                    writeClosureAudit(std::move(audit));
                }
                if (rejection != "accepted") continue;

                NewtonIterationInfo restorationTrace;
                restorationTrace.iter = acceptedIters;
                restorationTrace.event =
                    "final_global_continuity_restoration";
                restorationTrace.residualNorm = candidateNorm;
                restorationTrace.relativeResidualNorm =
                    ResidualNorm::relative(candidateNorm, initialNorm);
                restorationTrace.rawStepNorm = correctionNorm;
                restorationTrace.stepNorm = damping * correctionNorm;
                restorationTrace.dampingFactor = damping;
                restorationTrace.lineSearchAttempts = trial + 1;
                restorationTrace.lineSearchAccepted = true;
                restorationTrace.carrierRowConvergence = candidateRowEval;
                fillNewtonIterationResidualTrace(
                    restorationTrace, candidateX, candidateR,
                    activeRowWeights, mesh_.numNodes(), assembler,
                    cfg_.diagnostics);
                result.trace.push_back(restorationTrace);
                result.history.push_back(std::move(restorationTrace));
                finishConverged(
                    "final_global_continuity_restoration",
                    candidateX, candidateR, acceptedIters, candidateNorm,
                    candidateRowEval, candidateGlobalEval);
                return result;
            }
            if(splitRuntime)splitRuntime->setLow(restorationLow);
        }
    }

    result.converged = false;
    const bool finalBranchGuardSatisfied =
        contactMajorityQfBranchAccepts(acceptedX, acceptedIters);
    const std::string finalFailureReason =
        (finalGlobalEval.enforced && !finalGlobalEval.satisfied)
            ? std::string("global_continuity_closure")
            : (finalRowEval.enforced && !finalRowEval.satisfied)
            ? std::string("carrier_row_convergence")
            : (!finalBranchGuardSatisfied)
            ? std::string("contact_majority_qf_branch_guard")
            : (cfg_.blockAbsoluteConvergence.mode == "enforce" &&
               !blockAbsoluteConvergenceSatisfied(acceptedR))
            ? std::string("block_absolute_convergence")
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
        {},
        finalContactMajorityQfDrop);
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
    if (cfg_.returnBestFeasibleOnFailure && hasBestFeasibleIterate) {
        result.iters = bestFeasibleIters;
        result.finalResidualNorm = bestFeasibleNorm;
        result.finalBlockNorms = blockResidualInfo(
            bestFeasibleR, mesh_.numNodes());
        result.finalCarrierRowConvergence = bestFeasibleRowEval;
        result.finalGlobalContinuityClosure = bestFeasibleGlobalEval;
        result.solution = makeSolution(
            assembler, bestFeasibleX, bestFeasibleIters);
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
