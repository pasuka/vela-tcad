#include "vela/post/ElectricFieldDiagnostics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace vela {
namespace {

constexpr Real kMinGeometry = 1.0e-30;
constexpr Real kMinDet = 1.0e-60;

void validateNodeVector(const char* caller, const DeviceMesh& mesh, const VectorXd& values)
{
    if (values.size() < static_cast<int>(mesh.numNodes())) {
        throw std::invalid_argument(
            std::string(caller) + ": value vector has fewer entries than mesh nodes");
    }
}

CellField2 cellScalarGradient(const DeviceMesh& mesh, const Cell& cell, const VectorXd& value)
{
    CellField2 result;
    result.cellId = cell.id;
    result.regionId = cell.region_id;
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return result;

    const Node& n0 = mesh.getNode(cell.node_ids[0]);
    const Node& n1 = mesh.getNode(cell.node_ids[1]);
    const Node& n2 = mesh.getNode(cell.node_ids[2]);
    const Real dx1 = n1.x - n0.x;
    const Real dy1 = n1.y - n0.y;
    const Real dx2 = n2.x - n0.x;
    const Real dy2 = n2.y - n0.y;
    const Real det = dx1 * dy2 - dy1 * dx2;
    result.area = 0.5 * std::abs(det);
    if (std::abs(det) <= kMinGeometry || !std::isfinite(det))
        return result;

    const Real v0 = value(static_cast<int>(cell.node_ids[0]));
    const Real v1 = value(static_cast<int>(cell.node_ids[1]));
    const Real v2 = value(static_cast<int>(cell.node_ids[2]));
    const Real dv1 = v1 - v0;
    const Real dv2 = v2 - v0;
    const Real gradX = (dv1 * dy2 - dy1 * dv2) / det;
    const Real gradY = (dx1 * dv2 - dv1 * dx2) / det;
    result.vector = Point2{gradX, gradY};
    result.magnitude = result.vector.norm();
    result.valid = std::isfinite(result.magnitude);
    return result;
}

std::vector<std::vector<Index>> buildNodeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeCells(mesh.numNodes());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        for (Index nodeId : cell.node_ids) {
            if (nodeId < nodeCells.size())
                nodeCells[nodeId].push_back(cellId);
        }
    }
    return nodeCells;
}

std::vector<std::set<Index>> buildNodeRegions(const DeviceMesh& mesh)
{
    std::vector<std::set<Index>> nodeRegions(mesh.numNodes());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        for (Index nodeId : cell.node_ids) {
            if (nodeId < nodeRegions.size())
                nodeRegions[nodeId].insert(cell.region_id);
        }
    }
    return nodeRegions;
}

std::vector<std::map<Index, std::set<Index>>> buildRegionNodeNeighbors(const DeviceMesh& mesh)
{
    std::vector<std::map<Index, std::set<Index>>> neighbors(mesh.numNodes());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        if (cell.node_ids.size() < 3)
            continue;
        for (Index a = 0; a < cell.node_ids.size(); ++a) {
            const Index nodeA = cell.node_ids[a];
            if (nodeA >= neighbors.size())
                continue;
            for (Index b = 0; b < cell.node_ids.size(); ++b) {
                if (a == b)
                    continue;
                const Index nodeB = cell.node_ids[b];
                if (nodeB < mesh.numNodes())
                    neighbors[nodeA][cell.region_id].insert(nodeB);
            }
        }
    }
    return neighbors;
}

void setPrimaryFromSamples(NodeField2& node)
{
    node.valid = false;
    node.vector = Point2::Zero();
    node.magnitude = 0.0;
    if (node.regionSamples.empty())
        return;
    const auto it = std::find_if(
        node.regionSamples.begin(), node.regionSamples.end(),
        [](const auto& item) { return item.second.valid; });
    if (it == node.regionSamples.end())
        return;
    node.regionId = it->first;
    node.vector = it->second.vector;
    node.magnitude = it->second.magnitude;
    node.valid = true;
}

std::vector<NodeField2> areaAverageFromCells(const DeviceMesh& mesh,
                                             const std::vector<CellField2>& cellFields)
{
    std::vector<NodeField2> nodes(mesh.numNodes());
    std::vector<std::map<Index, Real>> totalArea(mesh.numNodes());
    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId)
        nodes[nodeId].nodeId = nodeId;

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        if (cellId >= cellFields.size() || !cellFields[cellId].valid || cellFields[cellId].area <= 0.0)
            continue;
        const Cell& cell = mesh.getCell(cellId);
        for (Index nodeId : cell.node_ids) {
            if (nodeId >= nodes.size())
                continue;
            auto& sample = nodes[nodeId].regionSamples[cell.region_id];
            sample.vector += cellFields[cellId].area * cellFields[cellId].vector;
            totalArea[nodeId][cell.region_id] += cellFields[cellId].area;
        }
    }

    for (Index nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        for (auto& item : nodes[nodeId].regionSamples) {
            const Real area = totalArea[nodeId][item.first];
            if (area <= 0.0)
                continue;
            item.second.vector /= area;
            item.second.magnitude = item.second.vector.norm();
            item.second.valid = std::isfinite(item.second.magnitude);
        }
        setPrimaryFromSamples(nodes[nodeId]);
    }
    return nodes;
}

RecoveredField2 leastSquaresRegionSample(const DeviceMesh& mesh,
                                          const VectorXd& potential,
                                          Index nodeId,
                                          Index regionId,
                                          ElectricFieldLeastSquaresWeight weight,
                                          const std::set<Index>& neighbors,
                                          const NodeField2* fallback)
{
    RecoveredField2 result;
    const Node& center = mesh.getNode(nodeId);
    const Real centerValue = potential(static_cast<int>(nodeId));
    Real sxx = 0.0;
    Real sxy = 0.0;
    Real syy = 0.0;
    Real sxv = 0.0;
    Real syv = 0.0;

    for (Index neighborId : neighbors) {
        const Node& neighbor = mesh.getNode(neighborId);
        const Real dx = neighbor.x - center.x;
        const Real dy = neighbor.y - center.y;
        const Real distance = std::hypot(dx, dy);
        if (distance <= kMinGeometry || !std::isfinite(distance))
            continue;
        const Real w = weight == ElectricFieldLeastSquaresWeight::InverseDistanceSquared
            ? 1.0 / (distance * distance)
            : 1.0 / distance;
        const Real dv = potential(static_cast<int>(neighborId)) - centerValue;
        sxx += w * dx * dx;
        sxy += w * dx * dy;
        syy += w * dy * dy;
        sxv += w * dx * dv;
        syv += w * dy * dv;
    }

    const Real det = sxx * syy - sxy * sxy;
    if (std::abs(det) <= kMinDet || !std::isfinite(det)) {
        if (fallback != nullptr) {
            const auto it = fallback->regionSamples.find(regionId);
            if (it != fallback->regionSamples.end())
                return it->second;
        }
        return result;
    }

    const Real gradX = (sxv * syy - syv * sxy) / det;
    const Real gradY = (sxx * syv - sxy * sxv) / det;
    result.vector = Point2{-gradX, -gradY};
    result.magnitude = result.vector.norm();
    result.valid = std::isfinite(result.magnitude);
    return result;
}

bool cellCircumcenter(const DeviceMesh& mesh, const Cell& cell, Point2& center)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return false;

    const Node& n0 = mesh.getNode(cell.node_ids[0]);
    const Node& n1 = mesh.getNode(cell.node_ids[1]);
    const Node& n2 = mesh.getNode(cell.node_ids[2]);
    const Real x0 = n0.x;
    const Real y0 = n0.y;
    const Real x1 = n1.x;
    const Real y1 = n1.y;
    const Real x2 = n2.x;
    const Real y2 = n2.y;
    const Real d = 2.0 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
    if (std::abs(d) <= kMinGeometry || !std::isfinite(d))
        return false;

    const Real r0 = x0 * x0 + y0 * y0;
    const Real r1 = x1 * x1 + y1 * y1;
    const Real r2 = x2 * x2 + y2 * y2;
    center = Point2{
        (r0 * (y1 - y2) + r1 * (y2 - y0) + r2 * (y0 - y1)) / d,
        (r0 * (x2 - x1) + r1 * (x0 - x2) + r2 * (x1 - x0)) / d};
    return std::isfinite(center.x()) && std::isfinite(center.y());
}

std::vector<NodeField2> circumcenterRecoverFromCells(
    const DeviceMesh& mesh,
    const std::vector<CellField2>& cellFields,
    ElectricFieldCircumcenterWeight weight,
    const std::vector<NodeField2>& fallback)
{
    std::vector<NodeField2> nodes(mesh.numNodes());
    std::vector<std::map<Index, Real>> totalWeight(mesh.numNodes());
    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId)
        nodes[nodeId].nodeId = nodeId;

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        if (cellId >= cellFields.size() || !cellFields[cellId].valid || cellFields[cellId].area <= 0.0)
            continue;
        const Cell& cell = mesh.getCell(cellId);
        Point2 circumcenter = Point2::Zero();
        if (!cellCircumcenter(mesh, cell, circumcenter))
            continue;

        for (Index nodeId : cell.node_ids) {
            if (nodeId >= nodes.size())
                continue;
            const Node& node = mesh.getNode(nodeId);
            const Real distance = std::hypot(circumcenter.x() - node.x, circumcenter.y() - node.y);
            if (distance <= kMinGeometry || !std::isfinite(distance))
                continue;
            Real w = 1.0 / distance;
            if (weight == ElectricFieldCircumcenterWeight::AreaOverDistance)
                w *= cellFields[cellId].area;
            auto& sample = nodes[nodeId].regionSamples[cell.region_id];
            sample.vector += w * cellFields[cellId].vector;
            totalWeight[nodeId][cell.region_id] += w;
        }
    }

    for (Index nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        for (auto& item : nodes[nodeId].regionSamples) {
            const Real w = totalWeight[nodeId][item.first];
            if (w <= 0.0)
                continue;
            item.second.vector /= w;
            item.second.magnitude = item.second.vector.norm();
            item.second.valid = std::isfinite(item.second.magnitude);
        }

        if (nodeId < fallback.size()) {
            for (const auto& item : fallback[nodeId].regionSamples) {
                auto existing = nodes[nodeId].regionSamples.find(item.first);
                if (existing == nodes[nodeId].regionSamples.end() || !existing->second.valid)
                    nodes[nodeId].regionSamples[item.first] = item.second;
            }
        }
        setPrimaryFromSamples(nodes[nodeId]);
    }
    return nodes;
}
RecoveredField2 sprRegionSample(const DeviceMesh& mesh,
                                 Index nodeId,
                                 Index regionId,
                                 const std::vector<Index>& nodeCells,
                                 const std::vector<CellField2>& cellFields,
                                 const RecoveredField2& fallback)
{
    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhsX = Eigen::Vector3d::Zero();
    Eigen::Vector3d rhsY = Eigen::Vector3d::Zero();
    std::size_t samples = 0;

    for (Index cellId : nodeCells) {
        if (cellId >= cellFields.size())
            continue;
        const Cell& cell = mesh.getCell(cellId);
        const CellField2& field = cellFields[cellId];
        if (cell.region_id != regionId || !field.valid)
            continue;
        Point2 centroid = Point2::Zero();
        for (Index local = 0; local < 3; ++local) {
            const Node& node = mesh.getNode(cell.node_ids[local]);
            centroid += Point2{node.x, node.y};
        }
        centroid /= 3.0;
        const Eigen::Vector3d row(1.0, centroid.x(), centroid.y());
        normal += row * row.transpose();
        rhsX += row * field.vector.x();
        rhsY += row * field.vector.y();
        ++samples;
    }

    if (samples < 3 || std::abs(normal.determinant()) <= kMinDet || !std::isfinite(normal.determinant()))
        return fallback;

    const Eigen::Vector3d coeffX = normal.ldlt().solve(rhsX);
    const Eigen::Vector3d coeffY = normal.ldlt().solve(rhsY);
    const Node& node = mesh.getNode(nodeId);
    RecoveredField2 result;
    result.vector = Point2{
        coeffX(0) + coeffX(1) * node.x + coeffX(2) * node.y,
        coeffY(0) + coeffY(1) * node.x + coeffY(2) * node.y};
    result.magnitude = result.vector.norm();
    result.valid = std::isfinite(result.magnitude);
    return result.valid ? result : fallback;
}

} // namespace

Real maxEdgeElectricFieldMagnitude(const DeviceMesh& mesh, const VectorXd& potential_V)
{
    validateNodeVector("maxEdgeElectricFieldMagnitude", mesh, potential_V);

    Real maxField = 0.0;
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        if (!(edge.length > 0.0) || !std::isfinite(edge.length))
            continue;
        const int n0 = static_cast<int>(edge.n0);
        const int n1 = static_cast<int>(edge.n1);
        const Real dpsi = potential_V(n1) - potential_V(n0);
        const Real field = std::abs(dpsi) / edge.length;
        if (std::isfinite(field))
            maxField = std::max(maxField, field);
    }
    return maxField;
}

std::vector<CellField2> computeCellElectricField(const DeviceMesh& mesh,
                                                 const VectorXd& potential_V)
{
    validateNodeVector("computeCellElectricField", mesh, potential_V);
    std::vector<CellField2> fields(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        fields[cellId] = cellScalarGradient(mesh, mesh.getCell(cellId), potential_V);
        fields[cellId].vector = -fields[cellId].vector;
        fields[cellId].magnitude = fields[cellId].vector.norm();
        fields[cellId].valid = fields[cellId].valid && std::isfinite(fields[cellId].magnitude);
    }
    return fields;
}

std::vector<CellField2> computeCellGradElectronQuasiFermi(const DeviceMesh& mesh,
                                                          const VectorXd& electronQf_V)
{
    validateNodeVector("computeCellGradElectronQuasiFermi", mesh, electronQf_V);
    std::vector<CellField2> fields(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId)
        fields[cellId] = cellScalarGradient(mesh, mesh.getCell(cellId), electronQf_V);
    return fields;
}

std::vector<CellField2> computeCellGradHoleQuasiFermi(const DeviceMesh& mesh,
                                                      const VectorXd& holeQf_V)
{
    validateNodeVector("computeCellGradHoleQuasiFermi", mesh, holeQf_V);
    std::vector<CellField2> fields(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId)
        fields[cellId] = cellScalarGradient(mesh, mesh.getCell(cellId), holeQf_V);
    return fields;
}

std::vector<EdgeField2> computeEdgeElectricField(const DeviceMesh& mesh,
                                                 const VectorXd& potential_V)
{
    validateNodeVector("computeEdgeElectricField", mesh, potential_V);
    std::vector<EdgeField2> fields(mesh.numEdges());
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        EdgeField2 field;
        field.edgeId = edgeId;
        field.node0 = edge.n0;
        field.node1 = edge.n1;
        if (edge.length > kMinGeometry && std::isfinite(edge.length)) {
            const Node& n0 = mesh.getNode(edge.n0);
            const Node& n1 = mesh.getNode(edge.n1);
            const Real dpsi = potential_V(static_cast<int>(edge.n1)) -
                potential_V(static_cast<int>(edge.n0));
            const Point2 unit{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
            field.projected = -dpsi / edge.length;
            field.vector = field.projected * unit;
            field.magnitude = std::abs(field.projected);
            field.valid = std::isfinite(field.magnitude);
        }
        fields[edgeId] = field;
    }
    return fields;
}

std::vector<NodeField2> computeNodeElectricFieldAreaAverage(const DeviceMesh& mesh,
                                                            const VectorXd& potential_V)
{
    return areaAverageFromCells(mesh, computeCellElectricField(mesh, potential_V));
}

std::vector<NodeField2> computeNodeElectricFieldLeastSquares(
    const DeviceMesh& mesh,
    const VectorXd& potential_V,
    ElectricFieldLeastSquaresWeight weight)
{
    validateNodeVector("computeNodeElectricFieldLeastSquares", mesh, potential_V);
    const auto nodeRegions = buildNodeRegions(mesh);
    const auto neighbors = buildRegionNodeNeighbors(mesh);
    const auto fallback = computeNodeElectricFieldAreaAverage(mesh, potential_V);

    std::vector<NodeField2> nodes(mesh.numNodes());
    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        nodes[nodeId].nodeId = nodeId;
        for (Index regionId : nodeRegions[nodeId]) {
            static const std::set<Index> empty;
            const auto regionIt = neighbors[nodeId].find(regionId);
            const std::set<Index>& regionNeighbors = regionIt == neighbors[nodeId].end()
                ? empty
                : regionIt->second;
            nodes[nodeId].regionSamples[regionId] = leastSquaresRegionSample(
                mesh, potential_V, nodeId, regionId, weight, regionNeighbors, &fallback[nodeId]);
        }
        setPrimaryFromSamples(nodes[nodeId]);
    }
    return nodes;
}

std::vector<NodeField2> computeNodeElectricFieldCircumcenterRecovery(
    const DeviceMesh& mesh,
    const VectorXd& potential_V,
    ElectricFieldCircumcenterWeight weight)
{
    validateNodeVector("computeNodeElectricFieldCircumcenterRecovery", mesh, potential_V);
    const auto cellFields = computeCellElectricField(mesh, potential_V);
    const auto fallback = computeNodeElectricFieldAreaAverage(mesh, potential_V);
    return circumcenterRecoverFromCells(mesh, cellFields, weight, fallback);
}

std::vector<NodeField2> computeNodeElectricFieldSPR(const DeviceMesh& mesh,
                                                    const VectorXd& potential_V)
{
    validateNodeVector("computeNodeElectricFieldSPR", mesh, potential_V);
    const auto nodeRegions = buildNodeRegions(mesh);
    const auto nodeCells = buildNodeCellMap(mesh);
    const auto cellFields = computeCellElectricField(mesh, potential_V);
    const auto lsFallback = computeNodeElectricFieldLeastSquares(
        mesh, potential_V, ElectricFieldLeastSquaresWeight::InverseDistance);

    std::vector<NodeField2> nodes(mesh.numNodes());
    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        nodes[nodeId].nodeId = nodeId;
        for (Index regionId : nodeRegions[nodeId]) {
            RecoveredField2 fallback;
            const auto it = lsFallback[nodeId].regionSamples.find(regionId);
            if (it != lsFallback[nodeId].regionSamples.end())
                fallback = it->second;
            nodes[nodeId].regionSamples[regionId] = sprRegionSample(
                mesh, nodeId, regionId, nodeCells[nodeId], cellFields, fallback);
        }
        setPrimaryFromSamples(nodes[nodeId]);
    }
    return nodes;
}

} // namespace vela
