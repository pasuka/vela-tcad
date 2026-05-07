#include "vela/mesh/DeviceMesh.h"
#include "vela/mesh/BoxGeometryBuilder.h"
#include <cmath>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace vela {

// ------------------------------------------------------------------
// Population
// ------------------------------------------------------------------

void DeviceMesh::addNode(const Node& node)
{
    nodes_.push_back(node);
}

void DeviceMesh::addCell(const Cell& cell)
{
    cells_.push_back(cell);
}

void DeviceMesh::addRegion(const Region& region)
{
    regions_.push_back(region);
}

void DeviceMesh::addContact(const Contact& contact)
{
    contacts_.push_back(contact);
}

// ------------------------------------------------------------------
// Edge generation
// ------------------------------------------------------------------

/**
 * Iterate over all triangular cells, collect unique {min,max} node-pairs,
 * then compute their Euclidean lengths from node coordinates.
 * Box geometry is computed by the public buildEdges() wrapper.
 */
void DeviceMesh::buildEdgesOnly()
{
    edges_.clear();

    // Map from sorted node pair → edge index in edges_
    std::map<std::pair<Index,Index>, Index> edgeMap;

    for (const auto& cell : cells_) {
        if (cell.type != CellType::Tri3) continue;
        if (cell.node_ids.size() < 3) continue;

        // Three edges per triangle
        const Index nids[3] = {
            cell.node_ids[0],
            cell.node_ids[1],
            cell.node_ids[2]
        };

        for (int e = 0; e < 3; ++e) {
            Index a = nids[e];
            Index b = nids[(e + 1) % 3];
            if (a > b) std::swap(a, b);

            auto key = std::make_pair(a, b);
            if (edgeMap.find(key) == edgeMap.end()) {
                Edge edge;
                edge.id = edges_.size();
                edge.n0 = a;
                edge.n1 = b;

                const Node& na = getNode(a);
                const Node& nb = getNode(b);
                Real dx = nb.x - na.x;
                Real dy = nb.y - na.y;
                edge.length = std::sqrt(dx*dx + dy*dy);
                edge.couple = 0.0; // Box coupling computed after edge generation

                edgeMap[key] = edge.id;
                edges_.push_back(std::move(edge));
            }
        }
    }

}

void DeviceMesh::buildEdges()
{
    buildEdgesOnly();
    buildBoxGeometry();
}

void DeviceMesh::buildBoxGeometry()
{
    if (edges_.empty() && !cells_.empty())
        buildEdgesOnly();
    BoxGeometryBuilder::build(*this);
}

// ------------------------------------------------------------------
// Accessors
// ------------------------------------------------------------------

const Node& DeviceMesh::getNode(Index id) const
{
    if (id >= nodes_.size())
        throw std::out_of_range("Node id out of range: " + std::to_string(id));
    return nodes_[id];
}

const Edge& DeviceMesh::getEdge(Index id) const
{
    if (id >= edges_.size())
        throw std::out_of_range("Edge id out of range: " + std::to_string(id));
    return edges_[id];
}

const Cell& DeviceMesh::getCell(Index id) const
{
    if (id >= cells_.size())
        throw std::out_of_range("Cell id out of range: " + std::to_string(id));
    return cells_[id];
}

const Region& DeviceMesh::getRegion(Index id) const
{
    if (id >= regions_.size())
        throw std::out_of_range("Region id out of range: " + std::to_string(id));
    return regions_[id];
}

const Contact& DeviceMesh::getContact(Index id) const
{
    if (id >= contacts_.size())
        throw std::out_of_range("Contact id out of range: " + std::to_string(id));
    return contacts_[id];
}

} // namespace vela
