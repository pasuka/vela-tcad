#include "vela/preprocess/GeometryRegionBuilder.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace vela {
namespace {

std::runtime_error topologyError(const std::string& message)
{
    return std::runtime_error("GeometryRegionBuilder: " + message);
}

std::string pointKey(const Point2D& point)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<Real>::max_digits10)
        << point.x_um << "," << point.y_um;
    return out.str();
}

std::string edgeKey(const Point2D& a, const Point2D& b)
{
    const std::string aKey = pointKey(a);
    const std::string bKey = pointKey(b);
    if (aKey < bKey) {
        return aKey + "->" + bKey;
    }
    return bKey + "->" + aKey;
}

std::vector<Point2D> primitivePolygon(const GeometryPrimitiveIr& primitive)
{
    if (primitive.kind == GeometryPrimitiveKind::Polygon) {
        if (primitive.points.size() < 3) {
            throw topologyError("polygon '" + primitive.name + "' must contain at least 3 points.");
        }
        return primitive.points;
    }

    if (primitive.kind == GeometryPrimitiveKind::Rectangle) {
        if (primitive.points.size() != 2) {
            throw topologyError(
                "rectangle '" + primitive.name + "' must contain exactly 2 corner points.");
        }
        const Point2D p0 = primitive.points[0];
        const Point2D p1 = primitive.points[1];
        if (std::abs(p0.x_um - p1.x_um) == 0.0 || std::abs(p0.y_um - p1.y_um) == 0.0) {
            throw topologyError("rectangle '" + primitive.name + "' has zero area.");
        }
        const Real minX = std::min(p0.x_um, p1.x_um);
        const Real maxX = std::max(p0.x_um, p1.x_um);
        const Real minY = std::min(p0.y_um, p1.y_um);
        const Real maxY = std::max(p0.y_um, p1.y_um);
        return {
            Point2D{minX, minY},
            Point2D{maxX, minY},
            Point2D{maxX, maxY},
            Point2D{minX, maxY},
        };
    }

    throw topologyError("unsupported primitive kind for '" + primitive.name + "'.");
}

} // namespace

GeometryRegionTopology GeometryRegionBuilder::build(const DeviceIr2D& ir) const
{
    GeometryRegionTopology topology;
    std::unordered_map<std::string, Index> regionIdByName;
    std::set<std::string> seenRegionNames;

    for (const auto& region : ir.regions) {
        if (!seenRegionNames.insert(region.name).second) {
            throw topologyError("duplicate region name '" + region.name + "'.");
        }
        const Index regionId = static_cast<Index>(topology.regions.size());
        topology.regions.push_back(RegionTopologyTag{regionId, region.name, region.material});
        regionIdByName.emplace(region.name, regionId);
    }

    for (const auto& primitive : ir.geometry) {
        if (!regionIdByName.contains(primitive.region)) {
            const Index regionId = static_cast<Index>(topology.regions.size());
            topology.regions.push_back(RegionTopologyTag{regionId, primitive.region, primitive.material});
            regionIdByName.emplace(primitive.region, regionId);
        }
    }

    std::map<std::string, std::vector<std::size_t>> edgeToBoundaryIndices;
    for (const auto& primitive : ir.geometry) {
        const auto regionIt = regionIdByName.find(primitive.region);
        if (regionIt == regionIdByName.end()) {
            throw topologyError("primitive '" + primitive.name + "' references unknown region '" +
                                primitive.region + "'.");
        }

        const std::vector<Point2D> polygon = primitivePolygon(primitive);
        for (std::size_t edgeIndex = 0; edgeIndex < polygon.size(); ++edgeIndex) {
            const Point2D& start = polygon[edgeIndex];
            const Point2D& end = polygon[(edgeIndex + 1) % polygon.size()];

            BoundaryTopologyTag boundary;
            boundary.ref = primitive.name + "#edge" + std::to_string(edgeIndex);
            boundary.regionId = regionIt->second;
            boundary.startUm = start;
            boundary.endUm = end;

            topology.boundaries.push_back(boundary);
            edgeToBoundaryIndices[edgeKey(start, end)].push_back(topology.boundaries.size() - 1);
        }
    }

    for (const auto& [_, boundaryIndices] : edgeToBoundaryIndices) {
        if (boundaryIndices.size() > 2) {
            throw topologyError("non-manifold edge shared by more than two regions.");
        }
        if (boundaryIndices.size() != 2) {
            continue;
        }

        auto& a = topology.boundaries[boundaryIndices[0]];
        auto& b = topology.boundaries[boundaryIndices[1]];
        if (a.regionId == b.regionId) {
            throw topologyError("duplicate edge within region '" +
                                topology.regions[a.regionId].name + "'.");
        }
        a.isExternal = false;
        b.isExternal = false;
        a.adjacentRegion = topology.regions[b.regionId].name;
        b.adjacentRegion = topology.regions[a.regionId].name;
    }

    std::unordered_map<std::string, std::size_t> boundaryIndexByRef;
    boundaryIndexByRef.reserve(topology.boundaries.size());
    for (std::size_t i = 0; i < topology.boundaries.size(); ++i) {
        boundaryIndexByRef.emplace(topology.boundaries[i].ref, i);
    }

    for (const auto& contact : ir.contacts) {
        for (const auto& ref : contact.boundaryRefs) {
            const auto boundaryIt = boundaryIndexByRef.find(ref);
            if (boundaryIt == boundaryIndexByRef.end()) {
                throw topologyError("contact '" + contact.name +
                                    "' references unknown boundary '" + ref + "'.");
            }
            auto& boundary = topology.boundaries[boundaryIt->second];
            if (!contact.ownerRegion.empty() &&
                topology.regions[boundary.regionId].name != contact.ownerRegion) {
                throw topologyError("contact '" + contact.name + "' boundary '" + ref +
                                    "' owner region mismatch.");
            }
            if (!boundary.contactName.empty() && boundary.contactName != contact.name) {
                throw topologyError("boundary '" + ref + "' assigned to multiple contacts.");
            }
            boundary.contactName = contact.name;
        }
    }

    return topology;
}

} // namespace vela

