#pragma once

#include "vela/core/Types.h"
#include <array>
#include <vector>

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
    Real minAngleDegrees = 0.0;
    Real maxAngleDegrees = 0.0;
    Real minEdgeLength = 0.0;
    Index transferredCellBoxEdges = 0;
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
    enum class PoissonPermittivityPolicy { LegacyAverage, CellMaterial };
    enum class CellBoxPolicy { LegacyPositive, DelaunayTransfer };

    // Dimensionless d/l, ordered along (0,1), (1,2), (2,0).
    // Node identity guards against accidentally reusing a reordered mesh.
    struct PoissonCellEdgeCoefficients {
        std::array<Index, 3> nodeIds{};
        std::array<Real, 3> coefficients{};
    };
    enum class NodeVolumePolicy {
        Barycentric,
        MixedVoronoi,
    };

    struct Options {
        bool fallbackNegativeCotangent = true;
        bool warnOnNegativeCotangent = false;
        bool requireNonObtuse = false;
        NodeVolumePolicy nodeVolumePolicy = NodeVolumePolicy::Barycentric;
        PoissonPermittivityPolicy poissonPermittivityPolicy =
            PoissonPermittivityPolicy::LegacyAverage;
        std::vector<PoissonCellEdgeCoefficients> poissonCellEdgeCoefficients;
        CellBoxPolicy cellBoxPolicy = CellBoxPolicy::LegacyPositive;
    };

    static Real triangleArea(const Node& a, const Node& b, const Node& c);
    static void build(DeviceMesh& mesh);
    static void build(DeviceMesh& mesh, const Options& options);
    static GeometryBuildReport buildWithReport(DeviceMesh& mesh);
    static GeometryBuildReport buildWithReport(DeviceMesh& mesh, const Options& options);

private:
    static Real triangleArea(const DeviceMesh& mesh, const Cell& cell);
};

} // namespace vela
