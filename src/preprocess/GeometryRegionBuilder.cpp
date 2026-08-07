#include "vela/preprocess/GeometryRegionBuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace vela {
namespace {

std::runtime_error topologyError(const std::string& message)
{
    return std::runtime_error("GeometryRegionBuilder: " + message);
}

constexpr Real kTopologyToleranceUm = 1e-9;

Real cross(const Point2D& a, const Point2D& b, const Point2D& p)
{
    return (b.x_um - a.x_um) * (p.y_um - a.y_um) -
           (b.y_um - a.y_um) * (p.x_um - a.x_um);
}

Real parameterOnLine(const Point2D& point, const Point2D& start, const Point2D& end)
{
    const Real dx = end.x_um - start.x_um;
    const Real dy = end.y_um - start.y_um;
    const Real length2 = dx * dx + dy * dy;
    if (length2 <= 0.0) {
        return 0.0;
    }
    return ((point.x_um - start.x_um) * dx + (point.y_um - start.y_um) * dy) / length2;
}

Real perpendicularDistance(const Point2D& point, const Point2D& start, const Point2D& end)
{
    const Real dx = end.x_um - start.x_um;
    const Real dy = end.y_um - start.y_um;
    const Real length = std::sqrt(dx * dx + dy * dy);
    if (length <= 0.0) {
        return std::numeric_limits<Real>::infinity();
    }
    return std::abs(cross(start, end, point)) / length;
}

bool collinear(const Point2D& a, const Point2D& b, const Point2D& c, const Point2D& d)
{
    return perpendicularDistance(c, a, b) <= kTopologyToleranceUm &&
           perpendicularDistance(d, a, b) <= kTopologyToleranceUm;
}

bool pointOnSegment(const Point2D& point, const Point2D& start, const Point2D& end)
{
    if (perpendicularDistance(point, start, end) > kTopologyToleranceUm) {
        return false;
    }
    const Real t = parameterOnLine(point, start, end);
    return t >= -kTopologyToleranceUm && t <= 1.0 + kTopologyToleranceUm;
}

bool positiveLengthOverlap(const Point2D& a0,
                           const Point2D& a1,
                           const Point2D& b0,
                           const Point2D& b1)
{
    if (!collinear(a0, a1, b0, b1)) {
        return false;
    }
    const Real b0t = parameterOnLine(b0, a0, a1);
    const Real b1t = parameterOnLine(b1, a0, a1);
    const Real lo = std::max(Real{0.0}, std::min(b0t, b1t));
    const Real hi = std::min(Real{1.0}, std::max(b0t, b1t));
    return hi - lo > kTopologyToleranceUm;
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

    struct OriginalEdge {
        std::string ref;
        Index regionId = 0;
        Point2D start{};
        Point2D end{};
    };
    std::vector<OriginalEdge> originalEdges;
    std::unordered_map<std::string, std::size_t> originalEdgeByRef;

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

            const std::string ref = primitive.name + "#edge" + std::to_string(edgeIndex);
            originalEdgeByRef.emplace(ref, originalEdges.size());
            originalEdges.push_back(OriginalEdge{ref, regionIt->second, start, end});
        }
    }

    std::unordered_map<std::string, std::string> contactByOriginalRef;
    for (const auto& contact : ir.contacts) {
        for (const auto& ref : contact.boundaryRefs) {
            const auto edgeIt = originalEdgeByRef.find(ref);
            if (edgeIt == originalEdgeByRef.end()) {
                throw topologyError("contact '" + contact.name +
                                    "' references unknown boundary '" + ref + "'.");
            }
            const auto& edge = originalEdges[edgeIt->second];
            if (!contact.ownerRegion.empty() &&
                topology.regions[edge.regionId].name != contact.ownerRegion) {
                throw topologyError("contact '" + contact.name + "' boundary '" + ref +
                                    "' owner region mismatch.");
            }
            const auto [it, inserted] = contactByOriginalRef.emplace(ref, contact.name);
            if (!inserted && it->second != contact.name) {
                throw topologyError("boundary '" + ref + "' assigned to multiple contacts.");
            }
        }
    }

    for (std::size_t edgeIndex = 0; edgeIndex < originalEdges.size(); ++edgeIndex) {
        const auto& edge = originalEdges[edgeIndex];
        std::vector<Real> cuts{0.0, 1.0};
        for (const auto& other : originalEdges) {
            if (!positiveLengthOverlap(edge.start, edge.end, other.start, other.end)) {
                continue;
            }
            for (const Point2D endpoint : {other.start, other.end}) {
                const Real t = parameterOnLine(endpoint, edge.start, edge.end);
                if (t > kTopologyToleranceUm && t < 1.0 - kTopologyToleranceUm) {
                    cuts.push_back(t);
                }
            }
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end(), [](Real lhs, Real rhs) {
            return std::abs(lhs - rhs) <= kTopologyToleranceUm;
        }), cuts.end());

        const bool split = cuts.size() > 2;
        for (std::size_t part = 0; part + 1 < cuts.size(); ++part) {
            const Real t0 = cuts[part];
            const Real t1 = cuts[part + 1];
            if (t1 - t0 <= kTopologyToleranceUm) {
                continue;
            }
            const auto interpolate = [&](Real t) {
                return Point2D{
                    edge.start.x_um + t * (edge.end.x_um - edge.start.x_um),
                    edge.start.y_um + t * (edge.end.y_um - edge.start.y_um),
                };
            };
            const Point2D start = interpolate(t0);
            const Point2D end = interpolate(t1);
            const Point2D midpoint = interpolate(0.5 * (t0 + t1));

            std::vector<std::size_t> owners;
            std::set<Index> ownerRegions;
            for (std::size_t ownerIndex = 0; ownerIndex < originalEdges.size(); ++ownerIndex) {
                if (pointOnSegment(midpoint, originalEdges[ownerIndex].start, originalEdges[ownerIndex].end)) {
                    owners.push_back(ownerIndex);
                    if (!ownerRegions.insert(originalEdges[ownerIndex].regionId).second) {
                        throw topologyError(
                            "duplicate overlapping boundary within region '" +
                            topology.regions[originalEdges[ownerIndex].regionId].name + "'.");
                    }
                }
            }
            if (ownerRegions.size() > 2) {
                throw topologyError("non-manifold boundary subsegment shared by more than two regions.");
            }

            BoundaryTopologyTag boundary;
            boundary.ref = split
                ? edge.ref + "#part" + std::to_string(part)
                : edge.ref;
            boundary.regionId = edge.regionId;
            boundary.startUm = start;
            boundary.endUm = end;
            boundary.isExternal = ownerRegions.size() == 1;
            if (!boundary.isExternal) {
                for (const Index ownerRegion : ownerRegions) {
                    if (ownerRegion != edge.regionId) {
                        boundary.adjacentRegion = topology.regions[ownerRegion].name;
                        break;
                    }
                }
            }

            const auto contactIt = contactByOriginalRef.find(edge.ref);
            if (contactIt != contactByOriginalRef.end()) {
                if (!boundary.isExternal) {
                    throw topologyError(
                        "contact '" + contactIt->second + "' boundary '" + edge.ref +
                        "' covers an internal interface subsegment.");
                }
                boundary.contactName = contactIt->second;
            }
            topology.boundaries.push_back(std::move(boundary));
        }
    }

    return topology;
}

} // namespace vela
