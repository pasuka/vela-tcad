#include "vela/mesh/BoxGeometryBuilder.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/mesh/MeshEntity.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace vela {

namespace {

constexpr Real kDegenerateTol = 1.0e-30;

std::pair<Index, Index> edgeKey(Index a, Index b)
{
    if (a > b) std::swap(a, b);
    return {a, b};
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
    build(mesh, Options{});
}

void BoxGeometryBuilder::build(DeviceMesh& mesh, const Options& options)
{
    if (mesh.edges_.empty() && !mesh.cells_.empty())
        mesh.buildEdgesOnly();

    for (auto& node : mesh.nodes_)
        node.volume = 0.0;
    for (auto& edge : mesh.edges_)
        edge.couple = 0.0;

    std::map<std::pair<Index, Index>, Index> edgeMap;
    for (Index e = 0; e < mesh.edges_.size(); ++e)
        edgeMap[edgeKey(mesh.edges_[e].n0, mesh.edges_[e].n1)] = e;

    for (const auto& cell : mesh.cells_) {
        if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
            continue;

        const Index ids[3] = {cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        const Real area = triangleArea(mesh, cell);
        if (area <= 0.0)
            continue;

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
                if (options.warnOnNegativeCotangent)
                    warnNegativeCotangent(cell.id, a, b, cot);

                if (options.fallbackNegativeCotangent) {
                    localCouple = area / (3.0 * edge.length);
                } else {
                    localCouple = 0.0;
                }
            }

            if (localCouple < 0.0)
                localCouple = 0.0;
            edge.couple += localCouple;
        }
    }
}

} // namespace vela
