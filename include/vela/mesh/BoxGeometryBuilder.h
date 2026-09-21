#pragma once

#include "vela/core/Types.h"
#include <array>

namespace vela {

class DeviceMesh;
struct Cell;
struct Node;

/**
 * @brief Diagnostics gathered while computing box-method mesh geometry.
 */
struct GeometryBuildReport {
    Index totalCells = 0;
    Index degenerateCells = 0;
    Index negativeCotangentCount = 0;
    Index fallbackCount = 0;
    Index edgeAggregationRepairCount = 0;
    Index aggregatedNegativeEdgeCount = 0;
    Real minAngleDegrees = 0.0;
    Real maxAngleDegrees = 0.0;
    Real minEdgeLength = 0.0;
};

/**
 * @brief Builds 2-D box-method geometric coefficients for Tri3 meshes.
 *
 * Node control-volume areas are computed with a barycentric control volume
 * (one third of each adjacent triangle area). Edge couplings use the
 * cotangent formula with an optional positive barycentric fallback for
 * obtuse/invalid local contributions.
 */
class BoxGeometryBuilder {
public:
    enum class NodeVolumePolicy {
        Barycentric,
        MixedVoronoi,
    };

    enum class EdgeCouplingPolicy {
        LocalPositiveFallback,
        AssembledDelaunayFallback,
    };

    struct Options {
        bool fallbackNegativeCotangent = true;
        bool warnOnNegativeCotangent = false;
        bool requireNonObtuse = false;
        NodeVolumePolicy nodeVolumePolicy = NodeVolumePolicy::Barycentric;
        EdgeCouplingPolicy edgeCouplingPolicy = EdgeCouplingPolicy::LocalPositiveFallback;
    };

    static Real triangleArea(const Node& a, const Node& b, const Node& c);
    /// Exactly the local mixed-Voronoi shares used by buildWithReport.
    static std::array<Real,3> mixedVoronoiVertexVolumes(const DeviceMesh& mesh,const Cell& cell);
    static void build(DeviceMesh& mesh);
    static void build(DeviceMesh& mesh, const Options& options);
    static GeometryBuildReport buildWithReport(DeviceMesh& mesh);
    static GeometryBuildReport buildWithReport(DeviceMesh& mesh, const Options& options);

private:
    static Real triangleArea(const DeviceMesh& mesh, const Cell& cell);
};

} // namespace vela
