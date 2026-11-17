#include "vela/post/PathIonizationIntegral.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <unordered_set>

namespace vela {
namespace {

Real logExponentialAverage(Real exponentIntegral)
{
    if (std::abs(exponentIntegral) < 1.0e-8) {
        const Real d = exponentIntegral;
        return std::log1p(-0.5 * d + d * d / 6.0 - d * d * d / 24.0);
    }
    if (exponentIntegral > 50.0)
        return -std::log(exponentIntegral);
    if (exponentIntegral < -50.0)
        return -exponentIntegral - std::log(-exponentIntegral);
    return std::log(-std::expm1(-exponentIntegral) / exponentIntegral);
}

Real logAdd(Real lhs, Real rhs)
{
    if (!std::isfinite(lhs))
        return rhs;
    if (!std::isfinite(rhs))
        return lhs;
    const Real maximum = std::max(lhs, rhs);
    return maximum + std::log1p(std::exp(-std::abs(lhs - rhs)));
}

Real finiteFromLog(Real logarithm)
{
    if (!std::isfinite(logarithm))
        return 0.0;
    const Real maximumLog = std::log(std::numeric_limits<Real>::max());
    if (logarithm >= maximumLog)
        return std::numeric_limits<Real>::max();
    return std::exp(logarithm);
}

Real edgeStrength(const PathIonizationSegment& edge)
{
    return std::abs(edge.electricField);
}

struct AdjacencyEntry {
    std::size_t sampleIndex = 0;
    Index neighbor = 0;
};

std::vector<std::size_t> traceAwayFromSeed(
    Index start,
    Index previous,
    bool increasingPotential,
    const VectorXd& potential,
    const std::vector<PathIonizationSegment>& samples,
    const std::vector<std::vector<AdjacencyEntry>>& adjacency,
    Real minimumElectricField)
{
    std::vector<std::size_t> traced;
    std::unordered_set<Index> visited{previous, start};
    Index current = start;
    constexpr Real potentialTolerance = 1.0e-14;

    while (current < adjacency.size()) {
        std::size_t best = samples.size();
        Real bestAlignment = -2.0;
        Real bestGradient = -1.0;
        Real bestField = -1.0;
        Index bestNeighbor = current;
        for (const AdjacencyEntry& entry : adjacency[current]) {
            if (entry.sampleIndex >= samples.size() || visited.contains(entry.neighbor))
                continue;
            const Real delta = potential(static_cast<int>(entry.neighbor)) -
                potential(static_cast<int>(current));
            if (increasingPotential ? delta <= potentialTolerance : delta >= -potentialTolerance)
                continue;
            const PathIonizationSegment& sample = samples[entry.sampleIndex];
            if (!(sample.length > 0.0))
                continue;
            const Real gradient = std::abs(delta) / sample.length;
            const Real field = edgeStrength(sample);
            if (field < minimumElectricField)
                continue;
            Real alignment = -2.0;
            const Real vectorField = std::hypot(
                sample.electricFieldX, sample.electricFieldY);
            const Real directionNorm = std::hypot(
                sample.edgeDirectionX, sample.edgeDirectionY);
            if (vectorField > 0.0 && directionNorm > 0.0) {
                const Real orientation = sample.node0 == current ? 1.0 : -1.0;
                const Real moveX = orientation * sample.edgeDirectionX / directionNorm;
                const Real moveY = orientation * sample.edgeDirectionY / directionNorm;
                // E=-grad(psi): increasing-potential tracing follows -E,
                // while decreasing-potential tracing follows +E.
                const Real fieldDirection = increasingPotential ? -1.0 : 1.0;
                alignment = fieldDirection *
                    (moveX * sample.electricFieldX +
                     moveY * sample.electricFieldY) / vectorField;
            }
            const bool vectorChoice = alignment > -1.5 || bestAlignment > -1.5;
            const Index bestEdgeId = best < samples.size()
                ? samples[best].edgeId : std::numeric_limits<Index>::max();
            const bool better =
                (vectorChoice && alignment > bestAlignment + 1.0e-12) ||
                (vectorChoice && std::abs(alignment - bestAlignment) <= 1.0e-12 &&
                 gradient > bestGradient) ||
                (!vectorChoice && gradient > bestGradient) ||
                (!vectorChoice && gradient == bestGradient && field > bestField) ||
                (((!vectorChoice && gradient == bestGradient && field == bestField) ||
                  (vectorChoice && std::abs(alignment - bestAlignment) <= 1.0e-12 &&
                   gradient == bestGradient && field == bestField)) &&
                 sample.edgeId < bestEdgeId);
            if (better) {
                best = entry.sampleIndex;
                bestAlignment = alignment;
                bestGradient = gradient;
                bestField = field;
                bestNeighbor = entry.neighbor;
            }
        }
        if (best >= samples.size())
            break;
        traced.push_back(best);
        visited.insert(bestNeighbor);
        current = bestNeighbor;
    }
    return traced;
}

std::vector<Index> pathSignature(const IonizationPath& path)
{
    std::vector<Index> signature;
    signature.reserve(path.segments.size());
    for (const PathIonizationSegment& segment : path.segments)
        signature.push_back(segment.edgeId);
    return signature;
}

void reversePathOrientation(IonizationPath& path)
{
    std::reverse(path.segments.begin(), path.segments.end());
    for (PathIonizationSegment& segment : path.segments) {
        std::swap(segment.node0, segment.node1);
        std::swap(segment.startPointM, segment.endPointM);
        std::swap(segment.startPotential, segment.endPotential);
    }
    std::reverse(path.nodes.begin(), path.nodes.end());
}

std::vector<Index> pathCellSet(const IonizationPath& path)
{
    std::vector<Index> cells;
    cells.reserve(path.segments.size());
    for (const PathIonizationSegment& segment : path.segments)
        cells.push_back(segment.cellId);
    std::sort(cells.begin(), cells.end());
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
    return cells;
}

bool sameContinuousPathCorridor(
    const IonizationPath& lhs,
    const IonizationPath& rhs)
{
    const Index invalidNode = std::numeric_limits<Index>::max();
    if (lhs.seedNodeId != invalidNode || rhs.seedNodeId != invalidNode) {
        // A continuous vector field has one physical streamline through a
        // nodal local maximum. Incident-cell launches only enumerate the
        // discretization ambiguity at that same point; keep one representative
        // even when the first few crossed cells differ strongly.
        if (lhs.seedNodeId == rhs.seedNodeId)
            return true;
        // Adjacent nodal maxima can still launch the same physical corridor.
        // Fall through to the cell-overlap test so a sub-element shift of the
        // seed does not consume two BreakAtIonIntegral ranks.
    }
    const std::vector<Index> a = pathCellSet(lhs);
    const std::vector<Index> b = pathCellSet(rhs);
    if (a.empty() || b.empty())
        return a == b;
    std::size_t intersection = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            ++intersection;
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return Real(intersection) / Real(std::min(a.size(), b.size())) >= 0.9;
}

struct NodePair {
    Index low = 0;
    Index high = 0;
    bool operator==(const NodePair&) const = default;
};

struct NodePairHash {
    std::size_t operator()(const NodePair& pair) const noexcept
    {
        const std::size_t h0 = std::hash<Index>{}(pair.low);
        const std::size_t h1 = std::hash<Index>{}(pair.high);
        return h0 ^ (h1 + 0x9e3779b97f4a7c15ULL + (h0 << 6U) + (h0 >> 2U));
    }
};

NodePair orderedNodePair(Index a, Index b)
{
    return a < b ? NodePair{a, b} : NodePair{b, a};
}

Real cross2(const Point2& a, const Point2& b)
{
    return a.x() * b.y() - a.y() * b.x();
}

Real interpolateTrianglePotential(
    const CellIonizationSample& sample,
    const Point2& point)
{
    const Point2 v0 = sample.verticesM[1] - sample.verticesM[0];
    const Point2 v1 = sample.verticesM[2] - sample.verticesM[0];
    const Point2 rhs = point - sample.verticesM[0];
    const Real determinant = cross2(v0, v1);
    if (std::abs(determinant) <= 1.0e-300)
        return sample.potentials[0];
    const Real w1 = cross2(rhs, v1) / determinant;
    const Real w2 = cross2(v0, rhs) / determinant;
    const Real w0 = 1.0 - w1 - w2;
    return w0 * sample.potentials[0] +
           w1 * sample.potentials[1] +
           w2 * sample.potentials[2];
}

Point2 interpolateTriangleElectricField(
    const CellIonizationSample& sample,
    const Point2& point)
{
    const bool hasVertexField = std::any_of(
        sample.vertexElectricFieldsVPerM.begin(),
        sample.vertexElectricFieldsVPerM.end(),
        [](const Point2& field) { return field.squaredNorm() > 0.0; });
    if (!hasVertexField)
        return sample.electricFieldVPerM;
    const Point2 v0 = sample.verticesM[1] - sample.verticesM[0];
    const Point2 v1 = sample.verticesM[2] - sample.verticesM[0];
    const Point2 rhs = point - sample.verticesM[0];
    const Real determinant = cross2(v0, v1);
    if (std::abs(determinant) <= 1.0e-300)
        return sample.electricFieldVPerM;
    const Real w1 = cross2(rhs, v1) / determinant;
    const Real w2 = cross2(v0, rhs) / determinant;
    const Real w0 = 1.0 - w1 - w2;
    return w0 * sample.vertexElectricFieldsVPerM[0] +
           w1 * sample.vertexElectricFieldsVPerM[1] +
           w2 * sample.vertexElectricFieldsVPerM[2];
}

enum class TracingDirectionFamily {
    Configured,
    CellElectric,
    ElectronCurrent,
    HoleCurrent,
    ElectronQf,
    HoleQf,
};

const std::array<Point2, 3>& tracingDirections(
    const CellIonizationSample& sample,
    TracingDirectionFamily family)
{
    switch (family) {
    case TracingDirectionFamily::ElectronCurrent:
        return sample.vertexElectronCurrentDirections;
    case TracingDirectionFamily::HoleCurrent:
        return sample.vertexHoleCurrentDirections;
    case TracingDirectionFamily::ElectronQf:
        return sample.vertexElectronQfDirections;
    case TracingDirectionFamily::HoleQf:
        return sample.vertexHoleQfDirections;
    case TracingDirectionFamily::Configured:
    default:
        return sample.vertexTracingDirections;
    }
}

Point2 interpolateTriangleTracingDirection(
    const CellIonizationSample& sample,
    const Point2& point,
    TracingDirectionFamily family)
{
    if (family == TracingDirectionFamily::CellElectric)
        return sample.electricFieldVPerM;
    const auto& directions = tracingDirections(sample, family);
    const bool hasTracingDirection = std::any_of(
        directions.begin(), directions.end(),
        [](const Point2& direction) { return direction.squaredNorm() > 0.0; });
    if (!hasTracingDirection)
        return interpolateTriangleElectricField(sample, point);
    const Point2 v0 = sample.verticesM[1] - sample.verticesM[0];
    const Point2 v1 = sample.verticesM[2] - sample.verticesM[0];
    const Point2 rhs = point - sample.verticesM[0];
    const Real determinant = cross2(v0, v1);
    if (std::abs(determinant) <= 1.0e-300)
        return directions[0];
    const Real w1 = cross2(rhs, v1) / determinant;
    const Real w2 = cross2(v0, rhs) / determinant;
    const Real w0 = 1.0 - w1 - w2;
    return w0 * directions[0] + w1 * directions[1] + w2 * directions[2];
}

struct ContinuousTraceHalf {
    std::vector<PathIonizationSegment> segments;
    std::vector<std::size_t> sampleIndices;
};

ContinuousTraceHalf traceContinuousHalf(
    std::size_t seedIndex,
    Real directionSign,
    const std::vector<CellIonizationSample>& samples,
    const std::vector<std::array<std::size_t, 3>>& neighbors,
    Real minimumElectricField,
    const Point2* initialPoint = nullptr,
    bool initialCellField = false,
    TracingDirectionFamily tracingFamily = TracingDirectionFamily::Configured,
    bool integrateP1Curve = false)
{
    ContinuousTraceHalf trace;
    if (seedIndex >= samples.size())
        return trace;
    const auto noNeighbor = samples.size();
    std::unordered_set<std::size_t> visited;
    std::size_t current = seedIndex;
    std::size_t entryLocalEdge = 3;
    Point2 point = initialPoint != nullptr
        ? *initialPoint
        : (samples[current].verticesM[0] +
           samples[current].verticesM[1] +
           samples[current].verticesM[2]) / 3.0;
    bool firstCell = true;

    while (current < samples.size() && !visited.contains(current)) {
        const CellIonizationSample& sample = samples[current];
        Point2 localElectricField = firstCell && initialCellField
            ? sample.electricFieldVPerM
            : interpolateTriangleElectricField(sample, point);
        Real fieldMagnitude = localElectricField.norm();
        if (fieldMagnitude < minimumElectricField || fieldMagnitude <= 0.0)
            break;
        visited.insert(current);
        Point2 tracingDirection = interpolateTriangleTracingDirection(
            sample, point, tracingFamily);
        if (!(tracingDirection.norm() > 0.0))
            tracingDirection = localElectricField;
        Point2 direction = directionSign * tracingDirection /
            tracingDirection.norm();

        Real bestDistance = std::numeric_limits<Real>::infinity();
        std::size_t exitLocalEdge = 3;
        Point2 exitPoint = point;
        Real cellScale = 0.0;
        for (std::size_t local = 0; local < 3; ++local) {
            cellScale = std::max(
                cellScale,
                (sample.verticesM[(local + 1) % 3] -
                 sample.verticesM[local]).norm());
            if (local == entryLocalEdge)
                continue;
            const Point2 edgeStart = sample.verticesM[local];
            const Point2 edgeVector =
                sample.verticesM[(local + 1) % 3] - edgeStart;
            const Real denominator = cross2(direction, edgeVector);
            if (std::abs(denominator) <= 1.0e-18)
                continue;
            const Point2 offset = edgeStart - point;
            const Real distance = cross2(offset, edgeVector) / denominator;
            const Real fraction = cross2(offset, direction) / denominator;
            const Real distanceTolerance =
                std::max(cellScale * 1.0e-12, Real{1.0e-18});
            if (distance <= distanceTolerance ||
                fraction < -1.0e-10 || fraction > 1.0 + 1.0e-10) {
                continue;
            }
            if (distance < bestDistance) {
                bestDistance = distance;
                exitLocalEdge = local;
                exitPoint = point + distance * direction;
            }
        }
        if ((exitLocalEdge >= 3 || !std::isfinite(bestDistance)) &&
            !integrateP1Curve)
            break;

        // Correct the straight predictor with the affine field at the chord
        // midpoint. This preserves exact edge crossing while following the
        // continuous, linearly interpolated nodal field inside the cell.
        const int correctionCount = integrateP1Curve ||
                (firstCell && initialCellField)
            ? 0 : 4;
        for (int correction = 0; correction < correctionCount; ++correction) {
            const Point2 midpoint = 0.5 * (point + exitPoint);
            const Point2 midpointField =
                interpolateTriangleTracingDirection(
                    sample, midpoint, tracingFamily);
            const Real midpointMagnitude = midpointField.norm();
            if (!(midpointMagnitude > 0.0))
                break;
            const Point2 correctedDirection =
                directionSign * midpointField / midpointMagnitude;
            Real correctedDistance = std::numeric_limits<Real>::infinity();
            std::size_t correctedEdge = 3;
            Point2 correctedExit = exitPoint;
            for (std::size_t local = 0; local < 3; ++local) {
                if (local == entryLocalEdge)
                    continue;
                const Point2 edgeStart = sample.verticesM[local];
                const Point2 edgeVector =
                    sample.verticesM[(local + 1) % 3] - edgeStart;
                const Real denominator = cross2(correctedDirection, edgeVector);
                if (std::abs(denominator) <= 1.0e-18)
                    continue;
                const Point2 offset = edgeStart - point;
                const Real distance =
                    cross2(offset, edgeVector) / denominator;
                const Real fraction =
                    cross2(offset, correctedDirection) / denominator;
                const Real distanceTolerance =
                    std::max(cellScale * 1.0e-12, Real{1.0e-18});
                if (distance <= distanceTolerance || fraction < -1.0e-10 ||
                    fraction > 1.0 + 1.0e-10) {
                    continue;
                }
                if (distance < correctedDistance) {
                    correctedDistance = distance;
                    correctedEdge = local;
                    correctedExit = point + distance * correctedDirection;
                }
            }
            if (correctedEdge >= 3 || !std::isfinite(correctedDistance))
                break;
            direction = correctedDirection;
            bestDistance = correctedDistance;
            exitLocalEdge = correctedEdge;
            exitPoint = correctedExit;
        }

        if (integrateP1Curve) {
            auto barycentric = [&](const Point2& p) {
                const Point2 v0 = sample.verticesM[1] - sample.verticesM[0];
                const Point2 v1 = sample.verticesM[2] - sample.verticesM[0];
                const Point2 rhs = p - sample.verticesM[0];
                const Real determinant = cross2(v0, v1);
                const Real w1 = cross2(rhs, v1) / determinant;
                const Real w2 = cross2(v0, rhs) / determinant;
                return std::array<Real, 3>{1.0 - w1 - w2, w1, w2};
            };
            auto inside = [&](const Point2& p) {
                const auto weights = barycentric(p);
                return std::all_of(
                    weights.begin(), weights.end(),
                    [](Real weight) { return weight >= -1.0e-10; });
            };
            auto unitDirection = [&](const Point2& p) -> Point2 {
                Point2 value = interpolateTriangleTracingDirection(
                    sample, p, tracingFamily);
                if (!(value.norm() > 0.0))
                    value = interpolateTriangleElectricField(sample, p);
                if (!(value.norm() > 0.0))
                    return Point2::Zero();
                return directionSign * value / value.norm();
            };

            Point2 curvePoint = point;
            Real arcLength = 0.0;
            exitLocalEdge = 3;
            const Real step = std::max(cellScale / 64.0, Real{1.0e-18});
            for (int substep = 0; substep < 4096; ++substep) {
                const Point2 k1 = unitDirection(curvePoint);
                if (!(k1.norm() > 0.0))
                    break;
                const Point2 k2 = unitDirection(curvePoint + 0.5 * step * k1);
                const Point2 k3 = unitDirection(curvePoint + 0.5 * step * k2);
                const Point2 k4 = unitDirection(curvePoint + step * k3);
                if (!(k2.norm() > 0.0) || !(k3.norm() > 0.0) ||
                    !(k4.norm() > 0.0)) {
                    break;
                }
                const Point2 proposed = curvePoint +
                    (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                if (inside(proposed)) {
                    arcLength += (proposed - curvePoint).norm();
                    curvePoint = proposed;
                    continue;
                }

                const Point2 delta = proposed - curvePoint;
                Real crossingFraction = std::numeric_limits<Real>::infinity();
                Point2 crossing = curvePoint;
                for (std::size_t local = 0; local < 3; ++local) {
                    if (local == entryLocalEdge)
                        continue;
                    const Point2 edgeStart = sample.verticesM[local];
                    const Point2 edgeVector =
                        sample.verticesM[(local + 1) % 3] - edgeStart;
                    const Real denominator = cross2(delta, edgeVector);
                    if (std::abs(denominator) <= 1.0e-30)
                        continue;
                    const Point2 offset = edgeStart - curvePoint;
                    const Real fraction = cross2(offset, edgeVector) / denominator;
                    const Real edgeFraction = cross2(offset, delta) / denominator;
                    if (fraction <= 1.0e-10 || fraction > 1.0 + 1.0e-10 ||
                        edgeFraction < -1.0e-10 || edgeFraction > 1.0 + 1.0e-10) {
                        continue;
                    }
                    if (fraction < crossingFraction) {
                        crossingFraction = fraction;
                        exitLocalEdge = local;
                        crossing = curvePoint + fraction * delta;
                    }
                }
                if (exitLocalEdge < 3 && std::isfinite(crossingFraction)) {
                    arcLength += (crossing - curvePoint).norm();
                    exitPoint = crossing;
                    direction = unitDirection(crossing);
                }
                break;
            }
            if (exitLocalEdge >= 3 || !(arcLength > 0.0))
                break;
            bestDistance = arcLength;
        }

        const Point2 midpoint = 0.5 * (point + exitPoint);
        localElectricField = firstCell && initialCellField
            ? sample.electricFieldVPerM
            : interpolateTriangleElectricField(sample, midpoint);
        fieldMagnitude = localElectricField.norm();

        PathIonizationSegment segment;
        segment.edgeId = sample.cellId;
        segment.cellId = sample.cellId;
        segment.node0 = sample.nodeIds[0];
        segment.node1 = sample.nodeIds[1];
        segment.length = bestDistance;
        segment.electronAlpha = sample.electronAlpha;
        segment.holeAlpha = sample.holeAlpha;
        segment.electricField = fieldMagnitude;
        segment.electronDrivingField = sample.electronDrivingField;
        segment.holeDrivingField = sample.holeDrivingField;
        segment.electricFieldX = localElectricField.x();
        segment.electricFieldY = localElectricField.y();
        segment.startPointM = point;
        segment.endPointM = exitPoint;
        segment.startPotential = interpolateTrianglePotential(sample, point);
        segment.endPotential = interpolateTrianglePotential(sample, exitPoint);
        segment.explicitGeometry = true;
        trace.segments.push_back(segment);
        trace.sampleIndices.push_back(current);
        firstCell = false;

        const std::size_t next = neighbors[current][exitLocalEdge];
        if (next == noNeighbor || next >= samples.size())
            break;
        std::size_t nextEntry = 3;
        const NodePair crossed = orderedNodePair(
            sample.nodeIds[exitLocalEdge],
            sample.nodeIds[(exitLocalEdge + 1) % 3]);
        for (std::size_t local = 0; local < 3; ++local) {
            const NodePair candidate = orderedNodePair(
                samples[next].nodeIds[local],
                samples[next].nodeIds[(local + 1) % 3]);
            if (candidate == crossed) {
                nextEntry = local;
                break;
            }
        }
        if (nextEntry >= 3)
            break;
        const Real nudge = std::max(cellScale * 1.0e-10, Real{1.0e-18});
        point = exitPoint + nudge * direction;
        current = next;
        entryLocalEdge = nextEntry;
    }
    return trace;
}

} // namespace

namespace {

struct SupportedSegmentRange {
    std::size_t begin = 0;
    std::size_t end = 0;
};

SupportedSegmentRange carrierSupportedRange(
    const std::vector<PathIonizationSegment>& segments,
    Real minimumDrivingField,
    bool electron)
{
    if (segments.empty())
        return {};
    if (minimumDrivingField <= 0.0)
        return {0, segments.size()};

    auto drive = [&](std::size_t index) {
        return std::abs(electron ? segments[index].electronDrivingField
                                 : segments[index].holeDrivingField);
    };
    std::size_t anchor = 0;
    for (std::size_t i = 1; i < segments.size(); ++i) {
        if (drive(i) > drive(anchor))
            anchor = i;
    }
    if (drive(anchor) < minimumDrivingField)
        return {};

    std::size_t begin = anchor;
    std::size_t end = anchor + 1;
    while (begin > 0 && drive(begin - 1) >= minimumDrivingField)
        --begin;
    while (end < segments.size() && drive(end) >= minimumDrivingField)
        ++end;
    return {begin, end};
}

std::pair<Real, Real> integrateCoupledRange(
    const std::vector<PathIonizationSegment>& orderedSegments,
    SupportedSegmentRange range)
{
    if (range.begin >= range.end)
        return {0.0, 0.0};
    Real futureElectronMinusHole = 0.0;
    for (std::size_t i = range.begin; i < range.end; ++i) {
        const PathIonizationSegment& segment = orderedSegments[i];
        if (!std::isfinite(segment.length) || segment.length < 0.0 ||
            !std::isfinite(segment.electronAlpha) || segment.electronAlpha < 0.0 ||
            !std::isfinite(segment.holeAlpha) || segment.holeAlpha < 0.0) {
            throw std::invalid_argument(
                "integrateIonizationPath: lengths and ionization coefficients "
                "must be finite and non-negative.");
        }
        futureElectronMinusHole +=
            (segment.electronAlpha - segment.holeAlpha) * segment.length;
    }

    Real pastHoleMinusElectron = 0.0;
    Real logElectronIntegral = -std::numeric_limits<Real>::infinity();
    Real logHoleIntegral = -std::numeric_limits<Real>::infinity();
    for (std::size_t i = range.begin; i < range.end; ++i) {
        const PathIonizationSegment& segment = orderedSegments[i];
        const Real electronMinusHole =
            (segment.electronAlpha - segment.holeAlpha) * segment.length;
        futureElectronMinusHole -= electronMinusHole;
        if (segment.electronAlpha > 0.0 && segment.length > 0.0) {
            const Real logTerm = std::log(segment.electronAlpha) +
                std::log(segment.length) - futureElectronMinusHole +
                logExponentialAverage(electronMinusHole);
            logElectronIntegral = logAdd(logElectronIntegral, logTerm);
        }

        const Real holeMinusElectron = -electronMinusHole;
        if (segment.holeAlpha > 0.0 && segment.length > 0.0) {
            const Real logTerm = std::log(segment.holeAlpha) +
                std::log(segment.length) - pastHoleMinusElectron +
                logExponentialAverage(holeMinusElectron);
            logHoleIntegral = logAdd(logHoleIntegral, logTerm);
        }
        pastHoleMinusElectron += holeMinusElectron;
    }
    return {finiteFromLog(logElectronIntegral), finiteFromLog(logHoleIntegral)};
}

Real alphaLengthOnRange(
    const std::vector<PathIonizationSegment>& segments,
    SupportedSegmentRange range,
    bool electron)
{
    Real value = 0.0;
    for (std::size_t i = range.begin; i < range.end; ++i) {
        value += (electron ? segments[i].electronAlpha : segments[i].holeAlpha) *
            segments[i].length;
    }
    return value;
}

Real segmentLengthOnRange(
    const std::vector<PathIonizationSegment>& segments,
    SupportedSegmentRange range)
{
    Real value = 0.0;
    for (std::size_t i = range.begin; i < range.end; ++i)
        value += segments[i].length;
    return value;
}

} // namespace

PathIonizationIntegral integrateIonizationPath(
    const std::vector<PathIonizationSegment>& orderedSegments,
    const PathIonizationIntegrationOptions& options)
{
    if (!std::isfinite(options.electronMinimumDrivingField) ||
        options.electronMinimumDrivingField < 0.0 ||
        !std::isfinite(options.holeMinimumDrivingField) ||
        options.holeMinimumDrivingField < 0.0) {
        throw std::invalid_argument(
            "integrateIonizationPath: carrier driving-field thresholds must "
            "be finite and non-negative.");
    }
    // Validate every segment even when a carrier-specific cutoff produces an
    // empty interval.
    for (const PathIonizationSegment& segment : orderedSegments) {
        if (!std::isfinite(segment.length) || segment.length < 0.0 ||
            !std::isfinite(segment.electronAlpha) || segment.electronAlpha < 0.0 ||
            !std::isfinite(segment.holeAlpha) || segment.holeAlpha < 0.0) {
            throw std::invalid_argument(
                "integrateIonizationPath: lengths and ionization coefficients "
                "must be finite and non-negative.");
        }
    }

    const SupportedSegmentRange electronRange = carrierSupportedRange(
        orderedSegments, options.electronMinimumDrivingField, true);
    const SupportedSegmentRange holeRange = carrierSupportedRange(
        orderedSegments, options.holeMinimumDrivingField, false);
    const auto electronPair = integrateCoupledRange(orderedSegments, electronRange);
    const auto holePair = electronRange.begin == holeRange.begin &&
            electronRange.end == holeRange.end
        ? electronPair
        : integrateCoupledRange(orderedSegments, holeRange);

    PathIonizationIntegral result;
    result.electron = electronPair.first;
    result.hole = holePair.second;
    result.electronSupportLength = segmentLengthOnRange(
        orderedSegments, electronRange);
    result.holeSupportLength = segmentLengthOnRange(
        orderedSegments, holeRange);
    if (options.meanDefinition ==
        MeanIonizationDefinition::CarrierAlphaLengthArithmetic) {
        result.mean = 0.5 * alphaLengthOnRange(
            orderedSegments, electronRange, true) +
            0.5 * alphaLengthOnRange(orderedSegments, holeRange, false);
    } else {
        result.mean = 0.5 * result.electron + 0.5 * result.hole;
    }
    return result;
}

PathIonizationAnalysis analyzeIonizationPaths(
    std::size_t nodeCount,
    const VectorXd& potential,
    const std::vector<PathIonizationSegment>& edgeSamples,
    std::size_t maxPaths,
    Real minimumElectricField,
    const PathIonizationIntegrationOptions& integrationOptions)
{
    if (potential.size() != static_cast<int>(nodeCount)) {
        throw std::invalid_argument(
            "analyzeIonizationPaths: potential size must equal node count.");
    }
    if (!std::isfinite(minimumElectricField) || minimumElectricField < 0.0) {
        throw std::invalid_argument(
            "analyzeIonizationPaths: minimum electric field must be finite "
            "and non-negative.");
    }

    PathIonizationAnalysis result;
    result.electronNodeValues.assign(nodeCount, 0.0);
    result.holeNodeValues.assign(nodeCount, 0.0);
    result.meanNodeValues.assign(nodeCount, 0.0);
    if (edgeSamples.empty())
        return result;

    std::vector<std::vector<AdjacencyEntry>> adjacency(nodeCount);
    for (std::size_t i = 0; i < edgeSamples.size(); ++i) {
        const auto& sample = edgeSamples[i];
        if (sample.node0 >= nodeCount || sample.node1 >= nodeCount) {
            throw std::invalid_argument(
                "analyzeIonizationPaths: edge sample references an invalid node.");
        }
        adjacency[sample.node0].push_back({i, sample.node1});
        adjacency[sample.node1].push_back({i, sample.node0});
    }

    std::vector<std::size_t> seedIndices;
    constexpr Real relativeTolerance = 1.0e-12;
    for (std::size_t i = 0; i < edgeSamples.size(); ++i) {
        const auto& sample = edgeSamples[i];
        const Real strength = edgeStrength(sample);
        if (!(strength > 0.0) || strength < minimumElectricField ||
            !(sample.length > 0.0))
            continue;
        bool localMaximum = true;
        for (Index node : {sample.node0, sample.node1}) {
            for (const AdjacencyEntry& entry : adjacency[node]) {
                const Real adjacentStrength = edgeStrength(edgeSamples[entry.sampleIndex]);
                if (adjacentStrength > strength * (1.0 + relativeTolerance) ||
                    (std::abs(adjacentStrength - strength) <=
                         relativeTolerance * std::max(strength, Real{1.0}) &&
                     edgeSamples[entry.sampleIndex].edgeId < sample.edgeId)) {
                    localMaximum = false;
                    break;
                }
            }
            if (!localMaximum)
                break;
        }
        if (localMaximum)
            seedIndices.push_back(i);
    }
    if (seedIndices.empty()) {
        seedIndices.push_back(static_cast<std::size_t>(std::distance(
            edgeSamples.begin(),
            std::max_element(edgeSamples.begin(), edgeSamples.end(),
                [](const auto& a, const auto& b) {
                    return edgeStrength(a) < edgeStrength(b);
                }))));
    }

    std::vector<std::vector<Index>> signatures;
    for (std::size_t seedIndex : seedIndices) {
        const PathIonizationSegment& seed = edgeSamples[seedIndex];
        Index low = seed.node0;
        Index high = seed.node1;
        if (potential(static_cast<int>(low)) > potential(static_cast<int>(high)))
            std::swap(low, high);

        const std::vector<std::size_t> lower = traceAwayFromSeed(
            low, high, false, potential, edgeSamples, adjacency,
            minimumElectricField);
        const std::vector<std::size_t> upper = traceAwayFromSeed(
            high, low, true, potential, edgeSamples, adjacency,
            minimumElectricField);

        IonizationPath path;
        path.seedEdgeId = seed.edgeId;
        Index lowerEnd = low;
        for (std::size_t sampleIndex : lower) {
            const auto& segment = edgeSamples[sampleIndex];
            lowerEnd = segment.node0 == lowerEnd ? segment.node1 : segment.node0;
        }
        path.nodes.push_back(lowerEnd);
        for (auto it = lower.rbegin(); it != lower.rend(); ++it) {
            PathIonizationSegment segment = edgeSamples[*it];
            if (segment.node0 != path.nodes.back())
                std::swap(segment.node0, segment.node1);
            path.segments.push_back(segment);
            path.nodes.push_back(segment.node1);
        }
        PathIonizationSegment orientedSeed = seed;
        if (orientedSeed.node0 != path.nodes.back())
            std::swap(orientedSeed.node0, orientedSeed.node1);
        path.segments.push_back(orientedSeed);
        path.nodes.push_back(orientedSeed.node1);
        for (std::size_t sampleIndex : upper) {
            PathIonizationSegment segment = edgeSamples[sampleIndex];
            if (segment.node0 != path.nodes.back())
                std::swap(segment.node0, segment.node1);
            path.segments.push_back(segment);
            path.nodes.push_back(segment.node1);
        }
        // Sentaurus Eqs. 469-470 use x along the electric-field line. Since
        // E=-grad(psi), orient the ordered quadrature from high to low psi.
        reversePathOrientation(path);
        path.integral = integrateIonizationPath(path.segments, integrationOptions);
        const std::vector<Index> signature = pathSignature(path);
        if (std::find(signatures.begin(), signatures.end(), signature) == signatures.end()) {
            signatures.push_back(signature);
            result.paths.push_back(std::move(path));
        }
    }

    std::sort(result.paths.begin(), result.paths.end(), [](const auto& a, const auto& b) {
        if (a.integral.mean != b.integral.mean)
            return a.integral.mean > b.integral.mean;
        return a.seedEdgeId < b.seedEdgeId;
    });
    if (maxPaths > 0 && result.paths.size() > maxPaths)
        result.paths.resize(maxPaths);

    for (const IonizationPath& path : result.paths) {
        for (Index node : path.nodes) {
            result.electronNodeValues[node] = std::max(
                result.electronNodeValues[node], path.integral.electron);
            result.holeNodeValues[node] = std::max(
                result.holeNodeValues[node], path.integral.hole);
            result.meanNodeValues[node] = std::max(
                result.meanNodeValues[node], path.integral.mean);
        }
    }
    return result;
}

PathIonizationAnalysis analyzeContinuousCellIonizationPaths(
    std::size_t nodeCount,
    const std::vector<CellIonizationSample>& cellSamples,
    std::size_t maxPaths,
    Real minimumElectricField,
    Real minimumSeedField,
    ContinuousPathDirection directionMode,
    const PathIonizationIntegrationOptions& integrationOptions,
    ContinuousPathRetention retention,
    ContinuousPathSeedMode seedMode,
    ContinuousPathTracingPolicy tracingPolicy,
    Real tracingQfRelativeFloor,
    std::optional<Real> seedPeakWindow)
{
    if (!std::isfinite(tracingQfRelativeFloor) ||
        tracingQfRelativeFloor < 0.0 || tracingQfRelativeFloor > 1.0) {
        throw std::invalid_argument(
            "analyzeContinuousCellIonizationPaths: tracing QF relative floor "
            "must be finite and in [0,1].");
    }
    if (!std::isfinite(minimumElectricField) || minimumElectricField < 0.0 ||
        !std::isfinite(minimumSeedField) || minimumSeedField < 0.0) {
        throw std::invalid_argument(
            "analyzeContinuousCellIonizationPaths: field thresholds must be "
            "finite and non-negative.");
    }
    if (seedPeakWindow && (!std::isfinite(*seedPeakWindow) || *seedPeakWindow < 0.0)) {
        throw std::invalid_argument(
            "analyzeContinuousCellIonizationPaths: peak window must be finite and non-negative.");
    }

    PathIonizationAnalysis result;
    result.electronNodeValues.assign(nodeCount, 0.0);
    result.holeNodeValues.assign(nodeCount, 0.0);
    result.meanNodeValues.assign(nodeCount, 0.0);
    if (cellSamples.empty())
        return result;

    const std::size_t noNeighbor = cellSamples.size();
    std::vector<std::array<std::size_t, 3>> neighbors(
        cellSamples.size(), {noNeighbor, noNeighbor, noNeighbor});
    struct EdgeOwner {
        std::size_t sampleIndex = 0;
        std::size_t localEdge = 0;
    };
    std::unordered_map<NodePair, EdgeOwner, NodePairHash> owner;
    std::vector<std::vector<std::size_t>> cellsAtNode(nodeCount);
    std::vector<std::vector<Index>> nodesAdjacent(nodeCount);
    std::vector<Point2> nodeFieldSum(nodeCount, Point2::Zero());
    std::vector<std::size_t> nodeFieldCount(nodeCount, 0);
    for (std::size_t sampleIndex = 0; sampleIndex < cellSamples.size(); ++sampleIndex) {
        const CellIonizationSample& sample = cellSamples[sampleIndex];
        for (std::size_t local = 0; local < 3; ++local) {
            const Index node = sample.nodeIds[local];
            if (node >= nodeCount) {
                throw std::invalid_argument(
                    "analyzeContinuousCellIonizationPaths: cell references an invalid node.");
            }
            cellsAtNode[node].push_back(sampleIndex);
            const Point2 vertexField = sample.vertexElectricFieldsVPerM[local];
            if (vertexField.squaredNorm() > 0.0) {
                nodeFieldSum[node] += vertexField;
                ++nodeFieldCount[node];
            }
        }
        const Real doubleArea = cross2(
            sample.verticesM[1] - sample.verticesM[0],
            sample.verticesM[2] - sample.verticesM[0]);
        if (std::abs(doubleArea) <= 1.0e-300) {
            throw std::invalid_argument(
                "analyzeContinuousCellIonizationPaths: degenerate triangle.");
        }
        for (std::size_t local = 0; local < 3; ++local) {
            const Index node0 = sample.nodeIds[local];
            const Index node1 = sample.nodeIds[(local + 1) % 3];
            nodesAdjacent[node0].push_back(node1);
            nodesAdjacent[node1].push_back(node0);
            const NodePair key = orderedNodePair(
                node0, node1);
            const auto [it, inserted] = owner.emplace(
                key, EdgeOwner{sampleIndex, local});
            if (!inserted) {
                const EdgeOwner other = it->second;
                neighbors[sampleIndex][local] = other.sampleIndex;
                neighbors[other.sampleIndex][other.localEdge] = sampleIndex;
            }
        }
    }

    Real seedThreshold = std::max(minimumElectricField, minimumSeedField);
    if (seedPeakWindow) {
        Real peak = 0.0;
        for (const auto& cell : cellSamples)
            peak = std::max(peak, cell.electricFieldVPerM.norm());
        seedThreshold = std::max(seedThreshold, peak - *seedPeakWindow);
    }
    constexpr Real relativeTolerance = 1.0e-12;
    struct ContinuousSeed {
        std::size_t sampleIndex = 0;
        Index nodeId = std::numeric_limits<Index>::max();
        Point2 pointM = Point2::Zero();
        Real field = 0.0;
        Real saddleField = 0.0;
        Real prominence = 0.0;
        Real prominenceRatio = 1.0;
        Index parentPeakNodeId = std::numeric_limits<Index>::max();
        Real electronQfRelativeMagnitude = 0.0;
        Real holeQfRelativeMagnitude = 0.0;
    };
    std::vector<ContinuousSeed> seeds;
    const bool useNodalField =
        !seedPeakWindow && seedMode == ContinuousPathSeedMode::NodalLocalMaxima &&
        std::any_of(
            nodeFieldCount.begin(), nodeFieldCount.end(),
            [](std::size_t count) { return count > 0; });
    std::vector<Point2> nodeFields(nodeCount, Point2::Zero());
    if (useNodalField) {
        for (Index node = 0; node < nodeCount; ++node) {
            if (nodeFieldCount[node] > 0)
                nodeFields[node] = nodeFieldSum[node] / Real(nodeFieldCount[node]);
        }
        for (Index node = 0; node < nodeCount; ++node) {
            const Real field = nodeFields[node].norm();
            if (!(field >= seedThreshold) || !(field > 0.0) ||
                cellsAtNode[node].empty()) {
                continue;
            }
            bool localMaximum = true;
            for (Index adjacentNode : nodesAdjacent[node]) {
                const Real adjacent = nodeFields[adjacentNode].norm();
                if (adjacent > field * (1.0 + relativeTolerance) ||
                    (std::abs(adjacent - field) <=
                         relativeTolerance * std::max(field, Real{1.0}) &&
                     adjacentNode < node)) {
                    localMaximum = false;
                    break;
                }
            }
            if (!localMaximum)
                continue;
            const std::size_t sampleIndex = *std::max_element(
                cellsAtNode[node].begin(), cellsAtNode[node].end(),
                [&](std::size_t a, std::size_t b) {
                    return cellSamples[a].electricFieldVPerM.norm() <
                           cellSamples[b].electricFieldVPerM.norm();
                });
            Point2 point = Point2::Zero();
            for (std::size_t local = 0; local < 3; ++local) {
                if (cellSamples[sampleIndex].nodeIds[local] == node) {
                    point = cellSamples[sampleIndex].verticesM[local];
                    break;
                }
            }
            ContinuousSeed seed{sampleIndex, node, point, field};
            for (std::size_t local = 0; local < 3; ++local) {
                if (cellSamples[sampleIndex].nodeIds[local] == node) {
                    seed.electronQfRelativeMagnitude =
                        cellSamples[sampleIndex]
                            .vertexElectronQfRelativeMagnitude[local];
                    seed.holeQfRelativeMagnitude =
                        cellSamples[sampleIndex]
                            .vertexHoleQfRelativeMagnitude[local];
                    break;
                }
            }
            seeds.push_back(seed);
        }
    } else {
        for (std::size_t sampleIndex = 0; sampleIndex < cellSamples.size(); ++sampleIndex) {
            const Real field = cellSamples[sampleIndex].electricFieldVPerM.norm();
            if (!(field >= seedThreshold) || !(field > 0.0))
                continue;
            bool localMaximum = true;
            for (std::size_t neighbor : neighbors[sampleIndex]) {
                if (neighbor == noNeighbor)
                    continue;
                const Real adjacent = cellSamples[neighbor].electricFieldVPerM.norm();
                if (adjacent > field * (1.0 + relativeTolerance) ||
                    (std::abs(adjacent - field) <=
                         relativeTolerance * std::max(field, Real{1.0}) &&
                     cellSamples[neighbor].cellId < cellSamples[sampleIndex].cellId)) {
                    localMaximum = false;
                    break;
                }
            }
            if (seedPeakWindow || localMaximum) {
                seeds.push_back({
                    sampleIndex,
                    std::numeric_limits<Index>::max(),
                    Point2::Zero(),
                    field});
            }
        }
    }

    if (useNodalField && seeds.size() > 1) {
        std::unordered_map<Index, Real> peakFields;
        peakFields.reserve(seeds.size());
        for (const ContinuousSeed& seed : seeds)
            peakFields.emplace(seed.nodeId, seed.field);

        struct BottleneckNode {
            Real field = 0.0;
            Index node = 0;
        };
        const auto weakerPriority = [](const BottleneckNode& lhs,
                                      const BottleneckNode& rhs) {
            if (lhs.field != rhs.field)
                return lhs.field < rhs.field;
            return lhs.node > rhs.node;
        };
        for (ContinuousSeed& seed : seeds) {
            std::vector<Real> bottleneck(nodeCount, -1.0);
            std::priority_queue<
                BottleneckNode,
                std::vector<BottleneckNode>,
                decltype(weakerPriority)> frontier(weakerPriority);
            bottleneck[seed.nodeId] = seed.field;
            frontier.push({seed.field, seed.nodeId});

            while (!frontier.empty()) {
                const BottleneckNode current = frontier.top();
                frontier.pop();
                if (current.field < bottleneck[current.node])
                    continue;

                const auto peak = peakFields.find(current.node);
                if (current.node != seed.nodeId && peak != peakFields.end()) {
                    const bool stronger =
                        peak->second > seed.field * (1.0 + relativeTolerance) ||
                        (std::abs(peak->second - seed.field) <=
                             relativeTolerance * std::max(seed.field, Real{1.0}) &&
                         current.node < seed.nodeId);
                    if (stronger) {
                        seed.saddleField = current.field;
                        seed.parentPeakNodeId = current.node;
                        break;
                    }
                }

                for (Index adjacent : nodesAdjacent[current.node]) {
                    const Real candidate = std::min(
                        current.field, nodeFields[adjacent].norm());
                    if (candidate > bottleneck[adjacent]) {
                        bottleneck[adjacent] = candidate;
                        frontier.push({candidate, adjacent});
                    }
                }
            }
            seed.prominence = std::max(seed.field - seed.saddleField, Real{0.0});
            seed.prominenceRatio = seed.field > 0.0
                ? seed.prominence / seed.field
                : 0.0;
        }
    } else {
        for (ContinuousSeed& seed : seeds) {
            seed.prominence = seed.field;
            seed.prominenceRatio = seed.field > 0.0 ? 1.0 : 0.0;
        }
    }

    std::vector<std::vector<Index>> signatures;
    for (const ContinuousSeed& seed : seeds) {
        TracingDirectionFamily tracingFamily = TracingDirectionFamily::Configured;
        const bool integrateP1Curve = tracingPolicy ==
            ContinuousPathTracingPolicy::P1ElectricFieldRungeKutta;
        if (tracingPolicy == ContinuousPathTracingPolicy::CellElectricField) {
            tracingFamily = TracingDirectionFamily::CellElectric;
        } else if (tracingPolicy ==
                ContinuousPathTracingPolicy::SentaurusEparallelAdaptive &&
            seed.nodeId < nodeCount) {
            const CellIonizationSample& seedSample = cellSamples[seed.sampleIndex];
            for (std::size_t local = 0; local < 3; ++local) {
                if (seedSample.nodeIds[local] != seed.nodeId)
                    continue;
                const bool nType = seedSample.vertexNetDoping[local] >= 0.0;
                if (seedSample.vertexTransportBoundary[local]) {
                    const Real qfRelativeMagnitude = nType
                        ? seedSample.vertexHoleQfRelativeMagnitude[local]
                        : seedSample.vertexElectronQfRelativeMagnitude[local];
                    if (qfRelativeMagnitude >= tracingQfRelativeFloor) {
                        tracingFamily = nType
                            ? TracingDirectionFamily::HoleQf
                            : TracingDirectionFamily::ElectronQf;
                    }
                } else {
                    tracingFamily = nType
                        ? TracingDirectionFamily::ElectronCurrent
                        : TracingDirectionFamily::HoleCurrent;
                }
                break;
            }
        }
        auto traceSeedHalves = [&](Real directionSign) {
            std::vector<ContinuousTraceHalf> halves;
            if (seed.nodeId >= nodeCount) {
                halves.push_back(traceContinuousHalf(
                    seed.sampleIndex, directionSign, cellSamples, neighbors,
                    minimumElectricField, nullptr, false, tracingFamily,
                    integrateP1Curve));
                return halves;
            }
            std::vector<std::vector<std::size_t>> halfSignatures;
            for (std::size_t startCell : cellsAtNode[seed.nodeId]) {
                const Point2 centroid =
                    (cellSamples[startCell].verticesM[0] +
                     cellSamples[startCell].verticesM[1] +
                     cellSamples[startCell].verticesM[2]) / 3.0;
                const Point2 startPoint =
                    seed.pointM + 1.0e-8 * (centroid - seed.pointM);
                ContinuousTraceHalf candidate = traceContinuousHalf(
                    startCell, directionSign, cellSamples, neighbors,
                    minimumElectricField, &startPoint, true, tracingFamily,
                    integrateP1Curve);
                if (candidate.segments.empty())
                    continue;
                if (std::find(
                        halfSignatures.begin(), halfSignatures.end(),
                        candidate.sampleIndices) != halfSignatures.end()) {
                    continue;
                }
                halfSignatures.push_back(candidate.sampleIndices);
                halves.push_back(std::move(candidate));
            }
            if (halves.empty())
                halves.emplace_back();
            return halves;
        };
        std::vector<ContinuousTraceHalf> lowers =
            directionMode == ContinuousPathDirection::OppositeVector
                ? std::vector<ContinuousTraceHalf>{ContinuousTraceHalf{}}
                : traceSeedHalves(1.0);
        std::vector<ContinuousTraceHalf> uppers =
            directionMode == ContinuousPathDirection::AlongVector
                ? std::vector<ContinuousTraceHalf>{ContinuousTraceHalf{}}
                : traceSeedHalves(-1.0);
        for (const ContinuousTraceHalf& lower : lowers) {
            for (const ContinuousTraceHalf& upper : uppers) {
                if (lower.segments.empty() && upper.segments.empty())
                    continue;

                IonizationPath path;
                path.seedEdgeId = cellSamples[seed.sampleIndex].cellId;
                path.seedCellId = cellSamples[seed.sampleIndex].cellId;
                path.seedNodeId = seed.nodeId;
                path.seedField = seed.field;
                path.saddleField = seed.saddleField;
                path.peakProminence = seed.prominence;
                path.peakProminenceRatio = seed.prominenceRatio;
                path.parentPeakNodeId = seed.parentPeakNodeId;
                path.seedElectronQfRelativeMagnitude =
                    seed.electronQfRelativeMagnitude;
                path.seedHoleQfRelativeMagnitude =
                    seed.holeQfRelativeMagnitude;
                for (auto it = lower.segments.rbegin();
                     it != lower.segments.rend(); ++it) {
                    PathIonizationSegment segment = *it;
                    std::swap(segment.startPointM, segment.endPointM);
                    std::swap(segment.startPotential, segment.endPotential);
                    path.segments.push_back(std::move(segment));
                }
                path.segments.insert(
                    path.segments.end(), upper.segments.begin(),
                    upper.segments.end());
                if (!path.segments.empty() &&
                    path.segments.front().startPotential <
                        path.segments.back().endPotential) {
                    reversePathOrientation(path);
                }
                for (const PathIonizationSegment& segment : path.segments) {
                    path.nodes.push_back(segment.node0);
                    path.nodes.push_back(segment.node1);
                }
                std::sort(path.nodes.begin(), path.nodes.end());
                path.nodes.erase(
                    std::unique(path.nodes.begin(), path.nodes.end()),
                    path.nodes.end());
                path.integral = integrateIonizationPath(
                    path.segments, integrationOptions);
                const std::vector<Index> signature = pathSignature(path);
                if (retention ==
                    ContinuousPathRetention::AllSeedTrajectories) {
                    // WriteAll retains separate path numbers for distinct
                    // incident-cell launches, even when their final field
                    // line and integrals are numerically identical.
                    result.paths.push_back(std::move(path));
                    continue;
                }
                if (retention == ContinuousPathRetention::CorridorDeduplicated &&
                    std::find(signatures.begin(), signatures.end(), signature) !=
                        signatures.end()) {
                    continue;
                }
                if (retention == ContinuousPathRetention::CorridorDeduplicated)
                    signatures.push_back(signature);
                auto equivalent = std::find_if(
                    result.paths.begin(), result.paths.end(),
                    [&](const IonizationPath& existing) {
                        if (retention ==
                                ContinuousPathRetention::DistinctLocalMaxima ||
                            retention ==
                                ContinuousPathRetention::NumberedPeakGroups) {
                            const Index invalidNode =
                                std::numeric_limits<Index>::max();
                            return path.seedNodeId != invalidNode
                                ? existing.seedNodeId == path.seedNodeId
                                : existing.seedNodeId == invalidNode &&
                                    existing.seedCellId == path.seedCellId;
                        }
                        if (sameContinuousPathCorridor(existing, path))
                            return true;
                        const Index invalidNode =
                            std::numeric_limits<Index>::max();
                        if (existing.seedNodeId == invalidNode ||
                            path.seedNodeId == invalidNode ||
                            existing.seedNodeId >= nodesAdjacent.size()) {
                            return false;
                        }
                        const auto& adjacent = nodesAdjacent[existing.seedNodeId];
                        if (std::find(
                            adjacent.begin(), adjacent.end(), path.seedNodeId) !=
                            adjacent.end()) {
                            return true;
                        }
                        // A P1 field can split one smooth maximum across two
                        // vertices with one intervening mesh node. Treat that
                        // two-edge seed neighborhood as one physical corridor.
                        return std::any_of(
                            adjacent.begin(), adjacent.end(),
                            [&](Index neighbor) {
                                if (neighbor >= nodesAdjacent.size())
                                    return false;
                                const auto& secondRing = nodesAdjacent[neighbor];
                                return std::find(
                                    secondRing.begin(), secondRing.end(),
                                    path.seedNodeId) != secondRing.end();
                            });
                    });
                if (equivalent == result.paths.end()) {
                    result.paths.push_back(std::move(path));
                } else if (path.integral.mean > equivalent->integral.mean) {
                    *equivalent = std::move(path);
                }
            }
        }
    }

    if (retention == ContinuousPathRetention::NumberedPeakGroups &&
        result.paths.size() > 1) {
        const Index invalidNode = std::numeric_limits<Index>::max();
        std::vector<std::size_t> parent(result.paths.size());
        std::iota(parent.begin(), parent.end(), std::size_t{0});
        auto root = [&](std::size_t value) {
            while (parent[value] != value) {
                parent[value] = parent[parent[value]];
                value = parent[value];
            }
            return value;
        };
        auto unite = [&](std::size_t a, std::size_t b) {
            a = root(a);
            b = root(b);
            if (a != b)
                parent[b] = a;
        };
        auto withinTwoNodeEdges = [&](Index a, Index b) {
            if (a == invalidNode || b == invalidNode ||
                a >= nodesAdjacent.size() || b >= nodesAdjacent.size()) {
                return false;
            }
            if (a == b)
                return true;
            const auto& first = nodesAdjacent[a];
            if (std::find(first.begin(), first.end(), b) != first.end())
                return true;
            return std::any_of(first.begin(), first.end(), [&](Index middle) {
                if (middle >= nodesAdjacent.size())
                    return false;
                const auto& second = nodesAdjacent[middle];
                return std::find(second.begin(), second.end(), b) != second.end();
            });
        };
        for (std::size_t i = 0; i < result.paths.size(); ++i) {
            for (std::size_t j = i + 1; j < result.paths.size(); ++j) {
                if (withinTwoNodeEdges(
                        result.paths[i].seedNodeId,
                        result.paths[j].seedNodeId)) {
                    unite(i, j);
                }
            }
        }
        std::vector<std::size_t> strongest(
            result.paths.size(), result.paths.size());
        std::vector<Index> groupId(
            result.paths.size(), std::numeric_limits<Index>::max());
        for (std::size_t i = 0; i < result.paths.size(); ++i) {
            const std::size_t group = root(i);
            groupId[group] = std::min(groupId[group], result.paths[i].seedNodeId);
            if (strongest[group] == result.paths.size() ||
                result.paths[i].integral.mean >
                    result.paths[strongest[group]].integral.mean) {
                strongest[group] = i;
            }
        }
        const std::vector<IonizationPath> original = result.paths;
        for (std::size_t i = 0; i < result.paths.size(); ++i) {
            const std::size_t best = strongest[root(i)];
            const Index seedEdgeId = result.paths[i].seedEdgeId;
            const Index seedCellId = result.paths[i].seedCellId;
            const Index seedNodeId = result.paths[i].seedNodeId;
            const Real seedField = result.paths[i].seedField;
            const Real saddleField = result.paths[i].saddleField;
            const Real peakProminence = result.paths[i].peakProminence;
            const Real peakProminenceRatio = result.paths[i].peakProminenceRatio;
            const Index parentPeakNodeId = result.paths[i].parentPeakNodeId;
            const Real seedElectronQfRelativeMagnitude =
                result.paths[i].seedElectronQfRelativeMagnitude;
            const Real seedHoleQfRelativeMagnitude =
                result.paths[i].seedHoleQfRelativeMagnitude;
            result.paths[i] = original[best];
            result.paths[i].seedEdgeId = seedEdgeId;
            result.paths[i].seedCellId = seedCellId;
            result.paths[i].seedNodeId = seedNodeId;
            result.paths[i].seedField = seedField;
            result.paths[i].saddleField = saddleField;
            result.paths[i].peakProminence = peakProminence;
            result.paths[i].peakProminenceRatio = peakProminenceRatio;
            result.paths[i].parentPeakNodeId = parentPeakNodeId;
            result.paths[i].seedElectronQfRelativeMagnitude =
                seedElectronQfRelativeMagnitude;
            result.paths[i].seedHoleQfRelativeMagnitude =
                seedHoleQfRelativeMagnitude;
            result.paths[i].physicalPathGroupId = groupId[root(i)];
        }
    }

    std::stable_sort(result.paths.begin(), result.paths.end(), [](const auto& a, const auto& b) {
        if (a.integral.mean != b.integral.mean)
            return a.integral.mean > b.integral.mean;
        if (a.seedNodeId != b.seedNodeId)
            return a.seedNodeId < b.seedNodeId;
        return a.seedCellId < b.seedCellId;
    });
    if (maxPaths > 0 && result.paths.size() > maxPaths)
        result.paths.resize(maxPaths);

    for (const IonizationPath& path : result.paths) {
        for (Index node : path.nodes) {
            result.electronNodeValues[node] = std::max(
                result.electronNodeValues[node], path.integral.electron);
            result.holeNodeValues[node] = std::max(
                result.holeNodeValues[node], path.integral.hole);
            result.meanNodeValues[node] = std::max(
                result.meanNodeValues[node], path.integral.mean);
        }
    }
    return result;
}

Real nthLargestMeanIonizationIntegral(
    const PathIonizationAnalysis& analysis,
    std::size_t oneBasedRank)
{
    if (oneBasedRank == 0)
        return 0.0;
    const Index invalid = std::numeric_limits<Index>::max();
    std::unordered_map<Index, Real> grouped;
    std::vector<Real> independent;
    independent.reserve(analysis.paths.size());
    for (std::size_t index = 0; index < analysis.paths.size(); ++index) {
        const IonizationPath& path = analysis.paths[index];
        if (path.physicalPathGroupId == invalid) {
            independent.push_back(path.integral.mean);
        } else {
            auto [it, inserted] = grouped.emplace(
                path.physicalPathGroupId, path.integral.mean);
            if (!inserted)
                it->second = std::max(it->second, path.integral.mean);
        }
    }
    for (const auto& [group, value] : grouped) {
        (void)group;
        independent.push_back(value);
    }
    std::sort(independent.begin(), independent.end(), std::greater<Real>());
    if (oneBasedRank > independent.size())
        return 0.0;
    return independent[oneBasedRank - 1];
}

Real nthLargestCarrierIonizationIntegral(
    const PathIonizationAnalysis& analysis,
    std::size_t oneBasedRank)
{
    if (oneBasedRank == 0)
        return 0.0;
    std::vector<Real> carrierIntegrals;
    carrierIntegrals.reserve(2 * analysis.paths.size());
    for (const IonizationPath& path : analysis.paths) {
        carrierIntegrals.push_back(path.integral.electron);
        carrierIntegrals.push_back(path.integral.hole);
    }
    if (oneBasedRank > carrierIntegrals.size())
        return 0.0;
    std::nth_element(
        carrierIntegrals.begin(),
        carrierIntegrals.begin() + static_cast<std::ptrdiff_t>(oneBasedRank - 1),
        carrierIntegrals.end(), std::greater<Real>{});
    return carrierIntegrals[oneBasedRank - 1];
}

} // namespace vela
