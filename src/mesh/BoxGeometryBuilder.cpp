#include "vela/mesh/BoxGeometryBuilder.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/mesh/MeshEntity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace vela {

namespace {

constexpr Real kDegenerateTol = 1.0e-30;
constexpr Real kPi = 3.141592653589793238462643383279502884;

std::pair<Index, Index> edgeKey(Index a, Index b)
{
    if (a > b) std::swap(a, b);
    return {a, b};
}

Real angleDegreesAt(const Node& a, const Node& b, const Node& c)
{
    const Real ux = b.x - a.x;
    const Real uy = b.y - a.y;
    const Real vx = c.x - a.x;
    const Real vy = c.y - a.y;
    const Real ul = std::sqrt(ux * ux + uy * uy);
    const Real vl = std::sqrt(vx * vx + vy * vy);
    if (ul < kDegenerateTol || vl < kDegenerateTol)
        return 0.0;

    const Real cosTheta = std::clamp((ux * vx + uy * vy) / (ul * vl), -1.0, 1.0);
    return std::acos(cosTheta) * 180.0 / kPi;
}

Real cotangentAtOppositeVertex(const Node& a, const Node& b, const Node& opp)
{
    const Real ux = a.x - opp.x;
    const Real uy = a.y - opp.y;
    const Real vx = b.x - opp.x;
    const Real vy = b.y - opp.y;

    const Real cross = ux * vy - uy * vx;
    if (std::abs(cross) < kDegenerateTol)
        return 0.0;

    const Real dot = ux * vx + uy * vy;
    return dot / std::abs(cross);
}

void warnNegativeCotangent(Index cellId, Index n0, Index n1, Real cot)
{
    std::cerr << "Vela warning: negative cotangent contribution " << cot
              << " for edge (" << n0 << ", " << n1 << ") in cell "
              << cellId << "; using non-negative barycentric box fallback.\n";
}

} // namespace

Real BoxGeometryBuilder::triangleArea(const Node& a, const Node& b, const Node& c)
{
    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
}

Real BoxGeometryBuilder::triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
        return 0.0;
    return triangleArea(mesh.getNode(cell.node_ids[0]),
                        mesh.getNode(cell.node_ids[1]),
                        mesh.getNode(cell.node_ids[2]));
}

void BoxGeometryBuilder::build(DeviceMesh& mesh)
{
    (void)buildWithReport(mesh, Options{});
}

void BoxGeometryBuilder::build(DeviceMesh& mesh, const Options& options)
{
    (void)buildWithReport(mesh, options);
}

GeometryBuildReport BoxGeometryBuilder::buildWithReport(DeviceMesh& mesh)
{
    return buildWithReport(mesh, Options{});
}

GeometryBuildReport BoxGeometryBuilder::buildWithReport(DeviceMesh& mesh, const Options& options)
{
    if (mesh.edges_.empty() && !mesh.cells_.empty())
        mesh.buildEdgesOnly();

    GeometryBuildReport report;
    report.totalCells = mesh.cells_.size();

    bool hasEdgeLength = false;
    for (const auto& edge : mesh.edges_) {
        if (!hasEdgeLength || edge.length < report.minEdgeLength) {
            report.minEdgeLength = edge.length;
            hasEdgeLength = true;
        }
    }

    for (auto& node : mesh.nodes_)
        node.volume = 0.0;
    for (auto& edge : mesh.edges_)
        edge.couple = 0.0;

    std::map<std::pair<Index, Index>, Index> edgeMap;
    for (Index e = 0; e < mesh.edges_.size(); ++e)
        edgeMap[edgeKey(mesh.edges_[e].n0, mesh.edges_[e].n1)] = e;

    bool hasAngle = false;
    for (const auto& cell : mesh.cells_) {
        if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
            continue;

        const Index ids[3] = {cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        const Node& n0 = mesh.getNode(ids[0]);
        const Node& n1 = mesh.getNode(ids[1]);
        const Node& n2 = mesh.getNode(ids[2]);
        const Real area = triangleArea(n0, n1, n2);
        if (area <= kDegenerateTol) {
            ++report.degenerateCells;
            continue;
        }

        const std::array<Real, 3> angles = {
            angleDegreesAt(n0, n1, n2),
            angleDegreesAt(n1, n2, n0),
            angleDegreesAt(n2, n0, n1),
        };
        for (Real angle : angles) {
            if (!hasAngle) {
                report.minAngleDegrees = angle;
                report.maxAngleDegrees = angle;
                hasAngle = true;
            } else {
                report.minAngleDegrees = std::min(report.minAngleDegrees, angle);
                report.maxAngleDegrees = std::max(report.maxAngleDegrees, angle);
            }
        }

        const Real nodeShare = area / 3.0;
        for (Index id : ids)
            mesh.nodes_.at(id).volume += nodeShare;

        for (int k = 0; k < 3; ++k) {
            const Index a = ids[k];
            const Index b = ids[(k + 1) % 3];
            const Index opp = ids[(k + 2) % 3];

            auto it = edgeMap.find(edgeKey(a, b));
            if (it == edgeMap.end()) {
                std::ostringstream os;
                os << "BoxGeometryBuilder: missing edge (" << a << ", " << b << ")";
                throw std::runtime_error(os.str());
            }

            Edge& edge = mesh.edges_.at(it->second);
            if (edge.length < kDegenerateTol)
                continue;

            const Real cot = cotangentAtOppositeVertex(mesh.getNode(a),
                                                       mesh.getNode(b),
                                                       mesh.getNode(opp));

            Real localCouple = 0.5 * cot * edge.length;
            if (cot < 0.0) {
                ++report.negativeCotangentCount;
                if (options.warnOnNegativeCotangent)
                    warnNegativeCotangent(cell.id, a, b, cot);

                if (options.fallbackNegativeCotangent) {
                    localCouple = area / (3.0 * edge.length);
                    ++report.fallbackCount;
                } else {
                    localCouple = 0.0;
                }
            }

            if (localCouple < 0.0)
                localCouple = 0.0;
            edge.couple += localCouple;
        }
    }

    return report;
}

} // namespace vela
