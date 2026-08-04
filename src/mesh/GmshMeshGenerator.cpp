#include "vela/mesh/GmshMeshGenerator.h"
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
#include <utility>
#include <vector>

namespace vela {
namespace {

struct TriIndex {
    std::size_t a = 0;
    std::size_t b = 0;
    std::size_t c = 0;
};

struct PointRef {
    Real x = 0.0;
    Real y = 0.0;
    Index globalNodeId = 0;
};

std::runtime_error meshError(const std::string& message)
{
    return std::runtime_error("GmshMeshGenerator: " + message);
}

std::string pointKey(Real x, Real y)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<Real>::max_digits10) << x << "," << y;
    return out.str();
}

std::vector<Point2D> primitivePolygon(const GeometryPrimitiveIr& primitive)
{
    if (primitive.kind == GeometryPrimitiveKind::Polygon) {
        if (primitive.points.size() < 3) {
            throw meshError("polygon '" + primitive.name + "' must have at least 3 points.");
        }
        return primitive.points;
    }

    if (primitive.kind == GeometryPrimitiveKind::Rectangle) {
        if (primitive.points.size() != 2) {
            throw meshError("rectangle '" + primitive.name + "' must have exactly 2 points.");
        }
        const Point2D p0 = primitive.points[0];
        const Point2D p1 = primitive.points[1];
        if (std::abs(p0.x_um - p1.x_um) == 0.0 || std::abs(p0.y_um - p1.y_um) == 0.0) {
            throw meshError("rectangle '" + primitive.name + "' has zero area.");
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

    throw meshError("unsupported primitive kind '" + primitive.name + "'.");
}

Real signedArea2(const Point2D& a, const Point2D& b, const Point2D& c)
{
    return (b.x_um - a.x_um) * (c.y_um - a.y_um) -
           (b.y_um - a.y_um) * (c.x_um - a.x_um);
}

bool isConvexPolygon(const std::vector<Point2D>& poly)
{
    if (poly.size() < 3) {
        return false;
    }
    int sign = 0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        const auto& c = poly[(i + 2) % poly.size()];
        const Real cross = signedArea2(a, b, c);
        if (std::abs(cross) < 1e-15) {
            continue;
        }
        const int current = cross > 0.0 ? 1 : -1;
        if (sign == 0) {
            sign = current;
            continue;
        }
        if (sign != current) {
            return false;
        }
    }
    return sign != 0;
}

bool pointInPolygon(const Point2D& point, const std::vector<Point2D>& poly)
{
    bool inside = false;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const auto& pi = poly[i];
        const auto& pj = poly[j];
        const bool intersect =
            ((pi.y_um > point.y_um) != (pj.y_um > point.y_um)) &&
            (point.x_um < (pj.x_um - pi.x_um) * (point.y_um - pi.y_um) /
                                  (pj.y_um - pi.y_um + 1e-30) +
                              pi.x_um);
        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

Real gaussianGradientMagnitudeCm3PerUm(const DopingProfileIr& profile, const Point2D& point)
{
    if (profile.kind != DopingProfileKind::Gaussian) {
        return 0.0;
    }
    if (profile.gaussianSigmaXUm <= 0.0 || profile.gaussianSigmaYUm <= 0.0) {
        return 0.0;
    }
    const Real dx = point.x_um - profile.gaussianCenterUm.x_um;
    const Real dy = point.y_um - profile.gaussianCenterUm.y_um;
    const Real exponent = -0.5 *
                          ((dx * dx) / (profile.gaussianSigmaXUm * profile.gaussianSigmaXUm) +
                           (dy * dy) / (profile.gaussianSigmaYUm * profile.gaussianSigmaYUm));
    const Real amplitude = profile.gaussianPeak_cm3 * std::exp(exponent);
    const Real gx = -dx / (profile.gaussianSigmaXUm * profile.gaussianSigmaXUm) * amplitude;
    const Real gy = -dy / (profile.gaussianSigmaYUm * profile.gaussianSigmaYUm) * amplitude;
    return std::sqrt(gx * gx + gy * gy);
}

std::pair<Point2D, Point2D> polygonBounds(const std::vector<Point2D>& poly)
{
    Point2D minPoint{poly.front().x_um, poly.front().y_um};
    Point2D maxPoint{poly.front().x_um, poly.front().y_um};
    for (const auto& p : poly) {
        minPoint.x_um = std::min(minPoint.x_um, p.x_um);
        minPoint.y_um = std::min(minPoint.y_um, p.y_um);
        maxPoint.x_um = std::max(maxPoint.x_um, p.x_um);
        maxPoint.y_um = std::max(maxPoint.y_um, p.y_um);
    }
    return {minPoint, maxPoint};
}

bool circumcircleContains(const PointRef& a, const PointRef& b, const PointRef& c, const PointRef& p)
{
    const Real ax = a.x - p.x;
    const Real ay = a.y - p.y;
    const Real bx = b.x - p.x;
    const Real by = b.y - p.y;
    const Real cx = c.x - p.x;
    const Real cy = c.y - p.y;

    const Real det =
        (ax * ax + ay * ay) * (bx * cy - by * cx) -
        (bx * bx + by * by) * (ax * cy - ay * cx) +
        (cx * cx + cy * cy) * (ax * by - ay * bx);

    const Real orient = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if (orient > 0.0) {
        return det > 1e-15;
    }
    return det < -1e-15;
}

std::vector<TriIndex> delaunayTriangulate(const std::vector<PointRef>& inputPoints)
{
    if (inputPoints.size() < 3) {
        return {};
    }

    std::vector<PointRef> points = inputPoints;
    Real minX = inputPoints.front().x;
    Real minY = inputPoints.front().y;
    Real maxX = inputPoints.front().x;
    Real maxY = inputPoints.front().y;
    for (const auto& p : inputPoints) {
        minX = std::min(minX, p.x);
        minY = std::min(minY, p.y);
        maxX = std::max(maxX, p.x);
        maxY = std::max(maxY, p.y);
    }

    const Real dx = maxX - minX;
    const Real dy = maxY - minY;
    const Real delta = std::max(dx, dy);
    const Real midX = 0.5 * (minX + maxX);
    const Real midY = 0.5 * (minY + maxY);

    const std::size_t superA = points.size();
    const std::size_t superB = points.size() + 1;
    const std::size_t superC = points.size() + 2;
    points.push_back(PointRef{midX - 20.0 * delta, midY - delta, 0});
    points.push_back(PointRef{midX, midY + 20.0 * delta, 0});
    points.push_back(PointRef{midX + 20.0 * delta, midY - delta, 0});

    std::vector<TriIndex> triangles;
    triangles.push_back(TriIndex{superA, superB, superC});

    for (std::size_t pIndex = 0; pIndex < inputPoints.size(); ++pIndex) {
        std::vector<TriIndex> badTriangles;
        std::vector<bool> isBad(triangles.size(), false);
        for (std::size_t t = 0; t < triangles.size(); ++t) {
            const auto& tri = triangles[t];
            if (circumcircleContains(points[tri.a], points[tri.b], points[tri.c], points[pIndex])) {
                isBad[t] = true;
                badTriangles.push_back(tri);
            }
        }

        std::map<std::pair<std::size_t, std::size_t>, int> edgeUseCount;
        for (const auto& tri : badTriangles) {
            const std::pair<std::size_t, std::size_t> edges[3] = {
                std::minmax(tri.a, tri.b),
                std::minmax(tri.b, tri.c),
                std::minmax(tri.c, tri.a),
            };
            for (const auto& edge : edges) {
                edgeUseCount[edge] += 1;
            }
        }

        std::vector<TriIndex> kept;
        kept.reserve(triangles.size());
        for (std::size_t t = 0; t < triangles.size(); ++t) {
            if (!isBad[t]) {
                kept.push_back(triangles[t]);
            }
        }
        triangles = std::move(kept);

        for (const auto& [edge, count] : edgeUseCount) {
            if (count != 1) {
                continue;
            }
            triangles.push_back(TriIndex{edge.first, edge.second, pIndex});
        }
    }

    std::vector<TriIndex> result;
    for (const auto& tri : triangles) {
        if (tri.a >= inputPoints.size() || tri.b >= inputPoints.size() || tri.c >= inputPoints.size()) {
            continue;
        }
        result.push_back(tri);
    }
    return result;
}

} // namespace

MeshBundle2D GmshMeshGenerator::generate(const DeviceIr2D& ir, const GmshMeshingOptions& options) const
{
    if (!options.delaunayTriOnly) {
        throw std::invalid_argument("GmshMeshGenerator: phase-1 only supports delaunayTriOnly=true.");
    }

    GeometryRegionBuilder topologyBuilder;
    const GeometryRegionTopology topology = topologyBuilder.build(ir);

    MeshBundle2D bundle;
    std::unordered_map<std::string, Index> nodeByCoord;
    std::vector<Region> meshRegions;
    meshRegions.reserve(topology.regions.size());
    for (const auto& region : topology.regions) {
        meshRegions.push_back(Region{region.id, region.name, region.material, {}});
    }

    auto getOrCreateNode = [&](Real x, Real y) {
        const std::string key = pointKey(x, y);
        const auto it = nodeByCoord.find(key);
        if (it != nodeByCoord.end()) {
            return it->second;
        }
        const Index nodeId = bundle.mesh.numNodes();
        bundle.mesh.addNode(Node{nodeId, x, y, 0.0});
        nodeByCoord.emplace(key, nodeId);
        return nodeId;
    };

    for (const auto& primitive : ir.geometry) {
        auto regionIt = std::find_if(topology.regions.begin(), topology.regions.end(), [&](const auto& item) {
            return item.name == primitive.region;
        });
        if (regionIt == topology.regions.end()) {
            throw meshError("primitive '" + primitive.name + "' maps to unknown region '" +
                            primitive.region + "'.");
        }
        const std::vector<Point2D> polygon = primitivePolygon(primitive);
        if (!isConvexPolygon(polygon)) {
            throw meshError("phase-1 triangulation currently supports convex polygons only ('" +
                            primitive.name + "').");
        }

        std::vector<PointRef> localPoints;
        std::set<Index> localNodeIds;
        localPoints.reserve(polygon.size() + 64);
        auto addLocalPoint = [&](Real x, Real y) {
            const Index nodeId = getOrCreateNode(x, y);
            if (localNodeIds.insert(nodeId).second) {
                localPoints.push_back(PointRef{x, y, nodeId});
            }
        };
        for (const auto& point : polygon) {
            addLocalPoint(point.x_um, point.y_um);
        }

        Real spacing = 0.0;
        if (const auto regionSizeIt = ir.meshControl.regionTargetSizeUm.find(primitive.region);
            regionSizeIt != ir.meshControl.regionTargetSizeUm.end()) {
            spacing = regionSizeIt->second;
        } else if (ir.meshControl.globalTargetSizeUm.has_value()) {
            spacing = *ir.meshControl.globalTargetSizeUm;
        }

        if (spacing > 0.0) {
            const auto [minPoint, maxPoint] = polygonBounds(polygon);
            const Real width = maxPoint.x_um - minPoint.x_um;
            const Real height = maxPoint.y_um - minPoint.y_um;
            if (width > spacing && height > spacing) {
                for (Real y = minPoint.y_um + spacing; y < maxPoint.y_um; y += spacing) {
                    for (Real x = minPoint.x_um + spacing; x < maxPoint.x_um; x += spacing) {
                        const Point2D candidate{x, y};
                        if (pointInPolygon(candidate, polygon)) {
                            addLocalPoint(x, y);
                        }
                    }
                }
            }
        }

        if (ir.meshControl.refineByDopingGradient) {
            Real localMinSize = ir.meshControl.minSizeUm > 0.0 ? ir.meshControl.minSizeUm : spacing * 0.5;
            if (localMinSize <= 0.0) {
                localMinSize = 0.1;
            }
            const Real gradientThreshold =
                ir.meshControl.dopingGradientThresholdCm3PerUm > 0.0
                    ? ir.meshControl.dopingGradientThresholdCm3PerUm
                    : 0.0;
            for (const auto& profile : ir.dopingProfiles) {
                if (profile.targetRegion != primitive.region ||
                    profile.kind != DopingProfileKind::Gaussian) {
                    continue;
                }
                const std::vector<Point2D> probes = {
                    profile.gaussianCenterUm,
                    Point2D{profile.gaussianCenterUm.x_um + profile.gaussianSigmaXUm * 0.5,
                            profile.gaussianCenterUm.y_um},
                    Point2D{profile.gaussianCenterUm.x_um - profile.gaussianSigmaXUm * 0.5,
                            profile.gaussianCenterUm.y_um},
                    Point2D{profile.gaussianCenterUm.x_um,
                            profile.gaussianCenterUm.y_um + profile.gaussianSigmaYUm * 0.5},
                    Point2D{profile.gaussianCenterUm.x_um,
                            profile.gaussianCenterUm.y_um - profile.gaussianSigmaYUm * 0.5},
                    Point2D{profile.gaussianCenterUm.x_um + localMinSize,
                            profile.gaussianCenterUm.y_um + localMinSize},
                    Point2D{profile.gaussianCenterUm.x_um - localMinSize,
                            profile.gaussianCenterUm.y_um - localMinSize},
                };
                for (const auto& probe : probes) {
                    if (pointInPolygon(probe, polygon) &&
                        gaussianGradientMagnitudeCm3PerUm(profile, probe) >= gradientThreshold) {
                        addLocalPoint(probe.x_um, probe.y_um);
                    }
                }
            }
        }

        const auto triangles = delaunayTriangulate(localPoints);
        for (const auto& tri : triangles) {
            const Point2D centroid{
                (localPoints[tri.a].x + localPoints[tri.b].x + localPoints[tri.c].x) / 3.0,
                (localPoints[tri.a].y + localPoints[tri.b].y + localPoints[tri.c].y) / 3.0,
            };
            if (!pointInPolygon(centroid, polygon)) {
                continue;
            }
            const Real area2 = signedArea2(
                Point2D{localPoints[tri.a].x, localPoints[tri.a].y},
                Point2D{localPoints[tri.b].x, localPoints[tri.b].y},
                Point2D{localPoints[tri.c].x, localPoints[tri.c].y});
            if (std::abs(area2) < 1e-14) {
                continue;
            }

            Cell cell;
            cell.id = bundle.mesh.numCells();
            cell.type = CellType::Tri3;
            cell.region_id = regionIt->id;
            cell.node_ids = {
                localPoints[tri.a].globalNodeId,
                localPoints[tri.b].globalNodeId,
                localPoints[tri.c].globalNodeId,
            };
            bundle.mesh.addCell(cell);
            meshRegions[cell.region_id].cell_ids.push_back(cell.id);
        }
    }

    for (const auto& region : meshRegions) {
        bundle.mesh.addRegion(region);
    }

    std::map<std::pair<std::string, Index>, std::set<Index>> contactNodes;
    for (const auto& boundary : topology.boundaries) {
        const std::string n0Key = pointKey(boundary.startUm.x_um, boundary.startUm.y_um);
        const std::string n1Key = pointKey(boundary.endUm.x_um, boundary.endUm.y_um);
        const auto n0It = nodeByCoord.find(n0Key);
        const auto n1It = nodeByCoord.find(n1Key);
        if (n0It == nodeByCoord.end() || n1It == nodeByCoord.end()) {
            throw meshError("boundary '" + boundary.ref + "' cannot be mapped to mesh nodes.");
        }
        bundle.boundaryEdges.push_back(
            BoundaryEdgeTag{n0It->second, n1It->second, boundary.contactName,
                            topology.regions[boundary.regionId].name});

        if (!boundary.contactName.empty()) {
            auto& nodes = contactNodes[{boundary.contactName, boundary.regionId}];
            nodes.insert(n0It->second);
            nodes.insert(n1It->second);
        }
    }

    Index contactId = 0;
    for (const auto& [key, nodes] : contactNodes) {
        Contact contact;
        contact.id = contactId++;
        contact.name = key.first;
        contact.region_id = key.second;
        contact.node_ids.assign(nodes.begin(), nodes.end());
        bundle.mesh.addContact(contact);
    }

    bundle.mesh.buildEdges();
    return bundle;
}

} // namespace vela
