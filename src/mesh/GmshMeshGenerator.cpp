#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/preprocess/GeometryRegionBuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <functional>
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

struct PrimitivePointCloud {
    std::vector<PointRef> points;
    std::set<Index> nodeIds;
};

struct PrimitiveBuildState {
    std::vector<Point2D> polygon;
    Real spacing = 0.0;
};

struct QuantizedPointKey {
    long long x = 0;
    long long y = 0;

    bool operator==(const QuantizedPointKey& other) const
    {
        return x == other.x && y == other.y;
    }
};

struct QuantizedPointKeyHash {
    std::size_t operator()(const QuantizedPointKey& key) const
    {
        return std::hash<long long>{}(key.x) ^ (std::hash<long long>{}(key.y) << 1U);
    }
};

std::runtime_error meshError(const std::string& message)
{
    return std::runtime_error("GmshMeshGenerator: " + message);
}

constexpr Real kCoordToleranceUm = 1e-9;

long long quantizeCoord(Real value)
{
    return static_cast<long long>(std::llround(value / kCoordToleranceUm));
}

QuantizedPointKey pointKey(Real x, Real y)
{
    return QuantizedPointKey{quantizeCoord(x), quantizeCoord(y)};
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
    if (profile.gaussianSigmaXUm <= 0.0) {
        return 0.0;
    }
    const Real dx = point.x_um - profile.gaussianPeakPosUm.x_um;
    const Real exponent = -0.5 *
                          ((dx * dx) / (profile.gaussianSigmaXUm * profile.gaussianSigmaXUm));
    const Real amplitude = profile.gaussianPeak_cm3 * std::exp(exponent);
    const Real gx = -dx / (profile.gaussianSigmaXUm * profile.gaussianSigmaXUm) * amplitude;
    return std::abs(gx);
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

bool pointOnSegment(const Point2D& p, const Point2D& a, const Point2D& b)
{
    const Real dx = b.x_um - a.x_um;
    const Real dy = b.y_um - a.y_um;
    const Real segmentLength = std::sqrt(dx * dx + dy * dy);
    if (segmentLength <= 0.0 ||
        std::abs(signedArea2(a, b, p)) / segmentLength > kCoordToleranceUm) {
        return false;
    }
    const Real minX = std::min(a.x_um, b.x_um) - 1e-10;
    const Real maxX = std::max(a.x_um, b.x_um) + 1e-10;
    const Real minY = std::min(a.y_um, b.y_um) - 1e-10;
    const Real maxY = std::max(a.y_um, b.y_um) + 1e-10;
    return p.x_um >= minX && p.x_um <= maxX && p.y_um >= minY && p.y_um <= maxY;
}

Real segmentParameter(const Point2D& point, const Point2D& start, const Point2D& end)
{
    const Real dx = end.x_um - start.x_um;
    const Real dy = end.y_um - start.y_um;
    const Real length2 = dx * dx + dy * dy;
    if (length2 <= 0.0) {
        return 0.0;
    }
    return ((point.x_um - start.x_um) * dx + (point.y_um - start.y_um) * dy) / length2;
}

bool segmentsColinear(const Point2D& a0, const Point2D& a1, const Point2D& b0, const Point2D& b1)
{
    const Real dx = a1.x_um - a0.x_um;
    const Real dy = a1.y_um - a0.y_um;
    const Real segmentLength = std::sqrt(dx * dx + dy * dy);
    return segmentLength > 0.0 &&
           std::abs(signedArea2(a0, a1, b0)) / segmentLength <= kCoordToleranceUm &&
           std::abs(signedArea2(a0, a1, b1)) / segmentLength <= kCoordToleranceUm;
}

bool overlappingSegmentPoints(const Point2D& a0,
                              const Point2D& a1,
                              const Point2D& b0,
                              const Point2D& b1,
                              Point2D& overlapStart,
                              Point2D& overlapEnd)
{
    if (!segmentsColinear(a0, a1, b0, b1)) {
        return false;
    }

    const Real bax = a1.x_um - a0.x_um;
    const Real bay = a1.y_um - a0.y_um;
    const Real length2 = bax * bax + bay * bay;
    if (length2 <= 0.0) {
        return false;
    }

    const Real t0 = segmentParameter(b0, a0, a1);
    const Real t1 = segmentParameter(b1, a0, a1);
    const Real lo = std::max(Real{0.0}, std::min(t0, t1));
    const Real hi = std::min(Real{1.0}, std::max(t0, t1));
    if (hi - lo <= 1e-12) {
        return false;
    }

    overlapStart = Point2D{a0.x_um + bax * lo, a0.y_um + bay * lo};
    overlapEnd = Point2D{a0.x_um + bax * hi, a0.y_um + bay * hi};
    return true;
}

Real distance(const Point2D& a, const Point2D& b)
{
    const Real dx = a.x_um - b.x_um;
    const Real dy = a.y_um - b.y_um;
    return std::sqrt(dx * dx + dy * dy);
}

void addPointUnique(PrimitivePointCloud& cloud,
                    Real x,
                    Real y,
                    const std::function<PointRef(Real, Real)>& getOrCreateNode)
{
    const PointRef point = getOrCreateNode(x, y);
    if (cloud.nodeIds.insert(point.globalNodeId).second) {
        cloud.points.push_back(point);
    }
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
    std::unordered_map<QuantizedPointKey, Index, QuantizedPointKeyHash> nodeByCoord;
    std::vector<Region> meshRegions;
    meshRegions.reserve(topology.regions.size());
    for (const auto& region : topology.regions) {
        meshRegions.push_back(Region{region.id, region.name, region.material, {}});
    }

    auto getOrCreateNode = [&](Real x, Real y) {
        const QuantizedPointKey key = pointKey(x, y);
        const auto it = nodeByCoord.find(key);
        if (it != nodeByCoord.end()) {
            const auto& existing = bundle.mesh.getNode(it->second);
            return PointRef{existing.x, existing.y, existing.id};
        }
        const Index nodeId = bundle.mesh.numNodes();
        bundle.mesh.addNode(Node{nodeId, x, y, 0.0});
        nodeByCoord.emplace(key, nodeId);
        return PointRef{x, y, nodeId};
    };

    std::vector<PrimitivePointCloud> primitiveClouds(ir.geometry.size());
    std::vector<PrimitiveBuildState> primitiveStates(ir.geometry.size());

    for (std::size_t primitiveIndex = 0; primitiveIndex < ir.geometry.size(); ++primitiveIndex) {
        const auto& primitive = ir.geometry[primitiveIndex];
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
        primitiveStates[primitiveIndex].polygon = polygon;

        auto& cloud = primitiveClouds[primitiveIndex];
        cloud.points.reserve(polygon.size() + 128);
        for (const auto& point : polygon) {
            addPointUnique(cloud, point.x_um, point.y_um, getOrCreateNode);
        }

        Real spacing = 0.0;
        if (const auto regionSizeIt = ir.meshControl.regionTargetSizeUm.find(primitive.region);
            regionSizeIt != ir.meshControl.regionTargetSizeUm.end()) {
            spacing = regionSizeIt->second;
        } else if (ir.meshControl.globalTargetSizeUm.has_value()) {
            spacing = *ir.meshControl.globalTargetSizeUm;
        }
        primitiveStates[primitiveIndex].spacing = spacing;

        if (spacing > 0.0) {
            for (std::size_t e = 0; e < polygon.size(); ++e) {
                const Point2D& a = polygon[e];
                const Point2D& b = polygon[(e + 1) % polygon.size()];
                const int segments =
                    std::max(1, static_cast<int>(std::ceil(distance(a, b) / spacing)));
                for (int i = 1; i < segments; ++i) {
                    const Real t = static_cast<Real>(i) / static_cast<Real>(segments);
                    addPointUnique(cloud,
                                   a.x_um + t * (b.x_um - a.x_um),
                                   a.y_um + t * (b.y_um - a.y_um),
                                   getOrCreateNode);
                }
            }

            const auto [minPoint, maxPoint] = polygonBounds(polygon);
            const Real fineSpacing = spacing / 2.0;
            for (Real x = minPoint.x_um + fineSpacing; x < maxPoint.x_um - 1e-12; x += fineSpacing) {
                for (Real y = minPoint.y_um + fineSpacing; y < maxPoint.y_um - 1e-12; y += fineSpacing) {
                    Point2D probe{x, y};
                    if (pointInPolygon(probe, polygon)) {
                        addPointUnique(cloud, x, y, getOrCreateNode);
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
                const auto [boundsMin, boundsMax] = polygonBounds(polygon);
                const Real midY = 0.5 * (boundsMin.y_um + boundsMax.y_um);
                const std::vector<Point2D> probes = {
                    Point2D{profile.gaussianPeakPosUm.x_um, midY},
                    Point2D{profile.gaussianPeakPosUm.x_um + profile.gaussianSigmaXUm * 0.5,
                            midY},
                    Point2D{profile.gaussianPeakPosUm.x_um - profile.gaussianSigmaXUm * 0.5,
                            midY},
                    Point2D{profile.gaussianPeakPosUm.x_um + localMinSize,
                            midY},
                    Point2D{profile.gaussianPeakPosUm.x_um - localMinSize,
                            midY},
                };
                for (const auto& probe : probes) {
                    const bool isPeakProbe = std::abs(probe.x_um - profile.gaussianPeakPosUm.x_um) < 1e-12;
                    if (pointInPolygon(probe, polygon) &&
                        (isPeakProbe ||
                         gaussianGradientMagnitudeCm3PerUm(profile, probe) >= gradientThreshold)) {
                        addPointUnique(cloud, probe.x_um, probe.y_um, getOrCreateNode);
                    }
                }
            }
        }

    }

    for (std::size_t i = 0; i < primitiveStates.size(); ++i) {
        for (std::size_t j = i + 1; j < primitiveStates.size(); ++j) {
            const auto& leftPoly = primitiveStates[i].polygon;
            const auto& rightPoly = primitiveStates[j].polygon;
            for (std::size_t le = 0; le < leftPoly.size(); ++le) {
                const Point2D& la = leftPoly[le];
                const Point2D& lb = leftPoly[(le + 1) % leftPoly.size()];
                for (std::size_t re = 0; re < rightPoly.size(); ++re) {
                    const Point2D& ra = rightPoly[re];
                    const Point2D& rb = rightPoly[(re + 1) % rightPoly.size()];
                    Point2D overlapStart{};
                    Point2D overlapEnd{};
                    if (!overlappingSegmentPoints(la, lb, ra, rb, overlapStart, overlapEnd)) {
                        continue;
                    }

                    std::vector<Point2D> mergedPoints;
                    const auto collectPoints = [&](std::size_t primitiveIndex, const Point2D& start, const Point2D& end) {
                        for (const auto& point : primitiveClouds[primitiveIndex].points) {
                            const Point2D candidate{point.x, point.y};
                            if (pointOnSegment(candidate, overlapStart, overlapEnd)) {
                                mergedPoints.push_back(candidate);
                            }
                        }
                        mergedPoints.push_back(overlapStart);
                        mergedPoints.push_back(overlapEnd);
                    };
                    collectPoints(i, la, lb);
                    collectPoints(j, ra, rb);

                    std::sort(mergedPoints.begin(), mergedPoints.end(), [&](const Point2D& lhs, const Point2D& rhs) {
                        return segmentParameter(lhs, overlapStart, overlapEnd) <
                               segmentParameter(rhs, overlapStart, overlapEnd);
                    });
                    mergedPoints.erase(std::unique(mergedPoints.begin(), mergedPoints.end(), [](const Point2D& lhs, const Point2D& rhs) {
                        return std::abs(lhs.x_um - rhs.x_um) <= kCoordToleranceUm &&
                               std::abs(lhs.y_um - rhs.y_um) <= kCoordToleranceUm;
                    }), mergedPoints.end());

                    for (const auto& point : mergedPoints) {
                        addPointUnique(primitiveClouds[i], point.x_um, point.y_um, getOrCreateNode);
                        addPointUnique(primitiveClouds[j], point.x_um, point.y_um, getOrCreateNode);
                    }
                }
            }
        }
    }

    for (std::size_t primitiveIndex = 0; primitiveIndex < ir.geometry.size(); ++primitiveIndex) {
        const auto& primitive = ir.geometry[primitiveIndex];
        auto regionIt = std::find_if(topology.regions.begin(), topology.regions.end(), [&](const auto& item) {
            return item.name == primitive.region;
        });
        const auto& cloud = primitiveClouds[primitiveIndex];
        const auto& polygon = primitiveStates[primitiveIndex].polygon;
        const auto triangles = delaunayTriangulate(cloud.points);
        for (const auto& tri : triangles) {
            const Point2D centroid{
                (cloud.points[tri.a].x + cloud.points[tri.b].x + cloud.points[tri.c].x) / 3.0,
                (cloud.points[tri.a].y + cloud.points[tri.b].y + cloud.points[tri.c].y) / 3.0,
            };
            if (!pointInPolygon(centroid, polygon)) {
                continue;
            }
            const Real area2 = signedArea2(
                Point2D{cloud.points[tri.a].x, cloud.points[tri.a].y},
                Point2D{cloud.points[tri.b].x, cloud.points[tri.b].y},
                Point2D{cloud.points[tri.c].x, cloud.points[tri.c].y});
            if (std::abs(area2) < 1e-14) {
                continue;
            }

            Cell cell;
            cell.id = bundle.mesh.numCells();
            cell.type = CellType::Tri3;
            cell.region_id = regionIt->id;
            cell.node_ids = {
                cloud.points[tri.a].globalNodeId,
                cloud.points[tri.b].globalNodeId,
                cloud.points[tri.c].globalNodeId,
            };
            bundle.mesh.addCell(cell);
            meshRegions[cell.region_id].cell_ids.push_back(cell.id);
        }
    }

    for (const auto& region : meshRegions) {
        bundle.mesh.addRegion(region);
    }

    bundle.mesh.buildEdges();
    std::set<std::pair<Index, Index>> meshEdges;
    for (const auto& edge : bundle.mesh.edges()) {
        meshEdges.insert(std::minmax(edge.n0, edge.n1));
    }

    std::unordered_map<std::string, Index> regionIdByName;
    for (const auto& region : topology.regions) {
        regionIdByName.emplace(region.name, region.id);
    }
    for (const auto& interface : topology.boundaries) {
        if (interface.isExternal || interface.adjacentRegion.empty()) {
            continue;
        }
        const auto adjacentIt = regionIdByName.find(interface.adjacentRegion);
        if (adjacentIt == regionIdByName.end()) {
            throw meshError("interface references unknown adjacent region '" +
                            interface.adjacentRegion + "'.");
        }
        if (interface.regionId > adjacentIt->second) {
            continue;
        }

        std::set<Index> interfaceNodes;
        for (const auto& node : bundle.mesh.nodes()) {
            if (pointOnSegment(Point2D{node.x, node.y}, interface.startUm, interface.endUm)) {
                interfaceNodes.insert(node.id);
            }
        }
        std::set<std::pair<Index, Index>> ownerEdges;
        std::set<std::pair<Index, Index>> adjacentEdges;
        std::set<Index> ownerNodes;
        std::set<Index> adjacentNodes;
        for (const auto& cell : bundle.mesh.cells()) {
            const bool owner = cell.region_id == interface.regionId;
            const bool adjacent = cell.region_id == adjacentIt->second;
            if (!owner && !adjacent) {
                continue;
            }
            for (const Index nodeId : cell.node_ids) {
                if (interfaceNodes.contains(nodeId)) {
                    (owner ? ownerNodes : adjacentNodes).insert(nodeId);
                }
            }
            for (int edgeIndex = 0; edgeIndex < 3; ++edgeIndex) {
                const Index a = cell.node_ids[edgeIndex];
                const Index b = cell.node_ids[(edgeIndex + 1) % 3];
                if (interfaceNodes.contains(a) && interfaceNodes.contains(b)) {
                    (owner ? ownerEdges : adjacentEdges).insert(std::minmax(a, b));
                }
            }
        }
        if (ownerNodes != interfaceNodes || adjacentNodes != interfaceNodes ||
            ownerEdges != adjacentEdges || ownerEdges.empty()) {
            throw meshError(
                "non-conformal interface between '" +
                topology.regions[interface.regionId].name + "' and '" +
                interface.adjacentRegion + "' over [(" +
                std::to_string(interface.startUm.x_um) + "," +
                std::to_string(interface.startUm.y_um) + "),(" +
                std::to_string(interface.endUm.x_um) + "," +
                std::to_string(interface.endUm.y_um) + ")].");
        }
    }

    std::map<std::pair<std::string, Index>, std::set<Index>> contactNodes;
    for (const auto& boundary : topology.boundaries) {
        if (!boundary.isExternal) {
            continue;
        }
        std::vector<std::pair<Real, Index>> boundaryNodes;
        for (const auto& node : bundle.mesh.nodes()) {
            const Point2D point{node.x, node.y};
            if (pointOnSegment(point, boundary.startUm, boundary.endUm)) {
                boundaryNodes.push_back({segmentParameter(point, boundary.startUm, boundary.endUm), node.id});
            }
        }
        std::sort(boundaryNodes.begin(), boundaryNodes.end(), [](const auto& lhs, const auto& rhs) {
            if (std::abs(lhs.first - rhs.first) > 1e-12) {
                return lhs.first < rhs.first;
            }
            return lhs.second < rhs.second;
        });
        boundaryNodes.erase(std::unique(boundaryNodes.begin(), boundaryNodes.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.second == rhs.second;
        }), boundaryNodes.end());
        if (boundaryNodes.size() < 2) {
            throw meshError("boundary '" + boundary.ref + "' cannot be mapped to mesh edge chain.");
        }

        for (std::size_t idx = 1; idx < boundaryNodes.size(); ++idx) {
            const Index node0 = boundaryNodes[idx - 1].second;
            const Index node1 = boundaryNodes[idx].second;
            if (!meshEdges.contains(std::minmax(node0, node1))) {
                throw meshError("boundary '" + boundary.ref + "' does not align with generated mesh edges.");
            }
            bundle.boundaryEdges.push_back(
                BoundaryEdgeTag{node0, node1, boundary.contactName, topology.regions[boundary.regionId].name});
        }

        if (!boundary.contactName.empty()) {
            auto& nodes = contactNodes[{boundary.contactName, boundary.regionId}];
            for (const auto& [_, nodeId] : boundaryNodes) {
                nodes.insert(nodeId);
            }
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

    return bundle;
}

} // namespace vela
