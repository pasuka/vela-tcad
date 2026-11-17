#pragma once

#include "vela/core/Types.h"

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace vela {

/// Piecewise-constant ionization coefficients on one ordered field-line segment.
struct PathIonizationSegment {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real length = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electricField = 0.0;
    Real electronDrivingField = 0.0;
    Real holeDrivingField = 0.0;
    Real electricFieldX = 0.0;
    Real electricFieldY = 0.0;
    Real edgeDirectionX = 0.0; ///< Unit vector from node0 to node1.
    Real edgeDirectionY = 0.0;
    Index cellId = std::numeric_limits<Index>::max();
    Point2 startPointM = Point2::Zero();
    Point2 endPointM = Point2::Zero();
    Real startPotential = 0.0;
    Real endPotential = 0.0;
    bool explicitGeometry = false;
};

/// Piecewise-constant physics used while tracing a field line through one
/// triangle. Vertex order and coordinates must match the mesh cell.
struct CellIonizationSample {
    Index cellId = 0;
    std::array<Index, 3> nodeIds{};
    std::array<Point2, 3> verticesM{};
    std::array<Real, 3> potentials{};
    Point2 electricFieldVPerM = Point2::Zero();
    Real electronDrivingField = 0.0;
    Real holeDrivingField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    /// Reconstructed nodal electric fields. When populated, the tracer uses
    /// their barycentric interpolation for a field continuous across cells.
    std::array<Point2, 3> vertexElectricFieldsVPerM{};
    /// Optional carrier-current direction used for path geometry. Electric
    /// field magnitude still controls seed and stop qualification.
    std::array<Point2, 3> vertexTracingDirections{};
    /// Direction families used by the Sentaurus Eparallel adaptive policy.
    /// Current directions should already contain their configured reliability
    /// fallback; quasi-Fermi families remain separately selectable.
    std::array<Point2, 3> vertexElectronCurrentDirections{};
    std::array<Point2, 3> vertexHoleCurrentDirections{};
    std::array<Point2, 3> vertexElectronQfDirections{};
    std::array<Point2, 3> vertexHoleQfDirections{};
    std::array<Real, 3> vertexElectronQfRelativeMagnitude{};
    std::array<Real, 3> vertexHoleQfRelativeMagnitude{};
    std::array<Real, 3> vertexNetDoping{};
    std::array<bool, 3> vertexTransportBoundary{};
};

/// Sentaurus-style electron- and hole-injection ionization integrals.
struct PathIonizationIntegral {
    Real electron = 0.0;
    Real hole = 0.0;
    Real mean = 0.0;
    Real electronSupportLength = 0.0;
    Real holeSupportLength = 0.0;
};

/// One discrete field path, ordered along E from higher to lower
/// electrostatic potential.
struct IonizationPath {
    Index seedEdgeId = 0;
    Index seedCellId = std::numeric_limits<Index>::max();
    Index seedNodeId = std::numeric_limits<Index>::max();
    /// Reconstructed nodal electric-field magnitude at the path seed.
    Real seedField = 0.0;
    /// Highest bottleneck field connecting this seed to any stronger peak.
    /// This is the merge-tree saddle level of the nodal field graph.
    Real saddleField = 0.0;
    /// Absolute and relative topographic prominence of the seed field.
    Real peakProminence = 0.0;
    Real peakProminenceRatio = 1.0;
    /// Stronger peak reached at saddleField, or invalid for the global peak.
    Index parentPeakNodeId = std::numeric_limits<Index>::max();
    /// Shared identifier for numbered aliases belonging to one physical peak
    /// corridor. Invalid means that every retained row is independently ranked.
    Index physicalPathGroupId = std::numeric_limits<Index>::max();
    Real seedElectronQfRelativeMagnitude = 0.0;
    Real seedHoleQfRelativeMagnitude = 0.0;
    std::vector<Index> nodes;
    std::vector<PathIonizationSegment> segments;
    PathIonizationIntegral integral;
};

struct PathIonizationAnalysis {
    std::vector<IonizationPath> paths; ///< Ordered by decreasing mean integral.
    std::vector<Real> electronNodeValues;
    std::vector<Real> holeNodeValues;
    std::vector<Real> meanNodeValues;
};

enum class ContinuousPathDirection {
    Bidirectional,
    AlongVector,
    OppositeVector,
};

enum class ContinuousPathRetention {
    /// Audit mode retaining every distinct incident-cell trajectory launched
    /// through each nodal local maximum.
    AllSeedTrajectories,
    /// Keep every nodal maximum numbered, while nearby maxima representing
    /// one P1 peak share the strongest traced trajectory. This reproduces
    /// WriteAll's equal-valued paths with separate path numbers.
    NumberedPeakGroups,
    /// Keep one strongest incident-cell trajectory for every nodal local
    /// maximum.
    DistinctLocalMaxima,
    /// Compatibility mode that merges adjacent/two-ring seeds sharing one
    /// high-field corridor.
    CorridorDeduplicated,
};

enum class ContinuousPathSeedMode {
    /// Generate seeds from local maxima of the element field, matching the
    /// element-oriented Sentaurus path search.
    CellLocalMaxima,
    /// Generate seeds from the reconstructed nodal field.
    NodalLocalMaxima,
};

enum class ContinuousPathTracingPolicy {
    /// Use vertexTracingDirections for every seed.
    Configured,
    /// Follow the element-constant electric field obtained directly from the
    /// P1 electrostatic-potential gradient. Reconstructed nodal fields remain
    /// responsible for local-maximum seed discovery.
    CellElectricField,
    /// Follow the continuous nodal P1 electric field with subcell RK4
    /// integration so curved streamlines choose their physical exit edge.
    P1ElectricFieldRungeKutta,
    /// For an interior peak, follow the majority-carrier current. For a peak
    /// on a transport boundary, follow the minority-carrier quasi-Fermi
    /// gradient. This avoids current-to-electric-field fallback truncation in
    /// low-current surface breakdown paths while retaining the resolved
    /// current direction in conducting bulk paths.
    SentaurusEparallelAdaptive,
};

/// Definition used to rank paths. WriteAll control runs establish the
/// arithmetic carrier mean as the exact numeric ordering; the separately
/// exported TDR MeanIonIntegral plateau remains an approximate geometry field.
enum class MeanIonizationDefinition {
    CarrierIntegralArithmetic,
    CarrierAlphaLengthArithmetic,
};

struct PathIonizationIntegrationOptions {
    /// Carrier-specific support cutoffs.  A positive value keeps the
    /// contiguous part of the shared geometric field line around that
    /// carrier's strongest driving field.  The other carrier coefficient is
    /// retained on that interval because Eqs. (469)-(470) are coupled.
    Real electronMinimumDrivingField = 0.0;
    Real holeMinimumDrivingField = 0.0;
    MeanIonizationDefinition meanDefinition =
        MeanIonizationDefinition::CarrierIntegralArithmetic;
};

/// Evaluate Eqs. (469)-(470) of the Sentaurus Device User Guide for an
/// explicitly ordered, piecewise-constant field path.
PathIonizationIntegral integrateIonizationPath(
    const std::vector<PathIonizationSegment>& orderedSegments,
    const PathIonizationIntegrationOptions& options = {});

/// Trace discrete paths through local electric-field maxima and evaluate the
/// Sentaurus-style path integrals.  The graph tracer follows monotone
/// electrostatic-potential edges in both directions from each local maximum.
PathIonizationAnalysis analyzeIonizationPaths(
    std::size_t nodeCount,
    const VectorXd& potential,
    const std::vector<PathIonizationSegment>& edgeSamples,
    std::size_t maxPaths = 0,
    Real minimumElectricField = 0.0,
    const PathIonizationIntegrationOptions& integrationOptions = {});

/// Trace piecewise-linear field lines through triangle interiors. Each cell
/// contributes its best-vertex ionization coefficients over the exact chord
/// cut by the local constant electric-field vector.
PathIonizationAnalysis analyzeContinuousCellIonizationPaths(
    std::size_t nodeCount,
    const std::vector<CellIonizationSample>& cellSamples,
    std::size_t maxPaths = 0,
    Real minimumElectricField = 0.0,
    Real minimumSeedField = 0.0,
    ContinuousPathDirection direction = ContinuousPathDirection::Bidirectional,
    const PathIonizationIntegrationOptions& integrationOptions = {},
    ContinuousPathRetention retention =
        ContinuousPathRetention::DistinctLocalMaxima,
    ContinuousPathSeedMode seedMode =
        ContinuousPathSeedMode::NodalLocalMaxima,
    ContinuousPathTracingPolicy tracingPolicy =
        ContinuousPathTracingPolicy::Configured,
    Real tracingQfRelativeFloor = 0.0,
    std::optional<Real> seedPeakWindow = std::nullopt);

/// Return the nth-largest mean integral (one-based), or zero when unavailable.
Real nthLargestMeanIonizationIntegral(
    const PathIonizationAnalysis& analysis,
    std::size_t oneBasedRank);

/// Return the nth-largest individual carrier integral after flattening the
/// electron and hole values from all paths. This is an explicit alternative
/// ordering for reverse-engineering BreakAtIonIntegral; path-mean ordering
/// remains the validated default.
Real nthLargestCarrierIonizationIntegral(
    const PathIonizationAnalysis& analysis,
    std::size_t oneBasedRank);

} // namespace vela
