#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <set>

#include "vela/io/NeutralMeshReader.h"
#include "vela/io/NeutralMeshWriter.h"
#include "vela/mesh/GmshMeshGenerator.h"

using namespace Catch::Matchers;
using namespace vela;

namespace {

Real edgeLength(const DeviceMesh& mesh, Index a, Index b)
{
    const auto& na = mesh.getNode(a);
    const auto& nb = mesh.getNode(b);
    const Real dx = na.x - nb.x;
    const Real dy = na.y - nb.y;
    return std::sqrt(dx * dx + dy * dy);
}

Real maxCellEdgeLength(const DeviceMesh& mesh)
{
    Real value = 0.0;
    for (const auto& cell : mesh.cells()) {
        value = std::max(value, edgeLength(mesh, cell.node_ids[0], cell.node_ids[1]));
        value = std::max(value, edgeLength(mesh, cell.node_ids[1], cell.node_ids[2]));
        value = std::max(value, edgeLength(mesh, cell.node_ids[2], cell.node_ids[0]));
    }
    return value;
}

std::vector<Index> nodesOnVerticalLine(const DeviceMesh& mesh, Real x)
{
    std::vector<Index> ids;
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        if (std::abs(mesh.getNode(i).x - x) < 1e-12) {
            ids.push_back(i);
        }
    }
    std::sort(ids.begin(), ids.end(), [&](Index lhs, Index rhs) {
        return mesh.getNode(lhs).y < mesh.getNode(rhs).y;
    });
    return ids;
}

bool cellUsesNode(const Cell& cell, Index nodeId)
{
    return cell.node_ids[0] == nodeId || cell.node_ids[1] == nodeId || cell.node_ids[2] == nodeId;
}

bool meshHasEdge(const DeviceMesh& mesh, Index a, Index b)
{
    const auto key = std::minmax(a, b);
    for (const auto& edge : mesh.edges()) {
        if (std::minmax(edge.n0, edge.n1) == key) {
            return true;
        }
    }
    return false;
}

std::set<std::pair<Index, Index>> regionEdgesOnVerticalSegment(
    const DeviceMesh& mesh, Index regionId, Real x, Real minY, Real maxY)
{
    std::set<std::pair<Index, Index>> edges;
    const auto onSegment = [&](Index nodeId) {
        const auto& node = mesh.getNode(nodeId);
        return std::abs(node.x - x) < 1e-10 &&
               node.y >= minY - 1e-10 && node.y <= maxY + 1e-10;
    };
    for (const auto& cell : mesh.cells()) {
        if (cell.region_id != regionId) {
            continue;
        }
        for (int edgeIndex = 0; edgeIndex < 3; ++edgeIndex) {
            const Index a = cell.node_ids[edgeIndex];
            const Index b = cell.node_ids[(edgeIndex + 1) % 3];
            if (onSegment(a) && onSegment(b)) {
                edges.insert(std::minmax(a, b));
            }
        }
    }
    return edges;
}

void checkConformalVerticalInterface(
    const MeshBundle2D& bundle, Real x, Real minY, Real maxY, Index leftRegion, Index rightRegion)
{
    const auto leftEdges =
        regionEdgesOnVerticalSegment(bundle.mesh, leftRegion, x, minY, maxY);
    const auto rightEdges =
        regionEdgesOnVerticalSegment(bundle.mesh, rightRegion, x, minY, maxY);
    REQUIRE_FALSE(leftEdges.empty());
    CHECK(leftEdges == rightEdges);

    std::set<Index> interfaceNodes;
    for (const auto& edge : leftEdges) {
        interfaceNodes.insert(edge.first);
        interfaceNodes.insert(edge.second);
    }
    for (const Index nodeId : interfaceNodes) {
        bool usedByLeft = false;
        bool usedByRight = false;
        for (const auto& cell : bundle.mesh.cells()) {
            if (!cellUsesNode(cell, nodeId)) {
                continue;
            }
            usedByLeft = usedByLeft || cell.region_id == leftRegion;
            usedByRight = usedByRight || cell.region_id == rightRegion;
        }
        CHECK(usedByLeft);
        CHECK(usedByRight);
    }
}

} // namespace

TEST_CASE("gmsh generator builds tri mesh for rectangle region", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);

    GmshMeshGenerator generator;
    const MeshBundle2D bundle = generator.generate(ir);
    CHECK(bundle.mesh.numNodes() == 4);
    CHECK(bundle.mesh.numCells() == 2);
    CHECK(bundle.mesh.numRegions() == 1);
    CHECK(bundle.boundaryEdges.size() == 4);
    CHECK(bundle.mesh.numEdges() == 5);
}

TEST_CASE("gmsh generator tags contacts on boundary refs", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);
    ir.contacts.push_back(ContactIr{"anode", "R.Si", {"rect_0#edge0"}});

    GmshMeshGenerator generator;
    const MeshBundle2D bundle = generator.generate(ir);
    REQUIRE(bundle.mesh.numContacts() == 1);
    CHECK(bundle.mesh.getContact(0).name == "anode");
    CHECK(bundle.mesh.getContact(0).node_ids.size() == 2);
}

TEST_CASE("gmsh generator rejects non-convex polygons in phase-1", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Poly", "Silicon", {0}});
    GeometryPrimitiveIr poly;
    poly.kind = GeometryPrimitiveKind::Polygon;
    poly.name = "poly_0";
    poly.region = "R.Poly";
    poly.material = "Silicon";
    poly.points = {
        {0.0, 0.0},
        {2.0, 0.0},
        {1.0, 0.5},
        {2.0, 1.0},
        {0.0, 1.0},
    };
    ir.geometry.push_back(poly);

    GmshMeshGenerator generator;
    REQUIRE_THROWS_WITH(generator.generate(ir), ContainsSubstring("convex polygons only"));
}

TEST_CASE("gmsh generator applies region size refinement", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D coarse;
    coarse.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{2.0, 1.0}};
    coarse.geometry.push_back(primitive);

    DeviceIr2D refined = coarse;
    refined.meshControl.regionTargetSizeUm["R.Si"] = 0.25;

    GmshMeshGenerator generator;
    const MeshBundle2D coarseMesh = generator.generate(coarse);
    const MeshBundle2D refinedMesh = generator.generate(refined);
    CHECK(refinedMesh.mesh.numCells() > coarseMesh.mesh.numCells());
}

TEST_CASE("gmsh generator enforces target size on 2x1 rectangle", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(primitive);
    ir.meshControl.regionTargetSizeUm["R.Si"] = 0.25;

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    CHECK(maxCellEdgeLength(mesh.mesh) <= Catch::Approx(0.25).margin(0.02));
}

TEST_CASE("gmsh generator enforces target size on thin rectangle", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Thin", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Thin";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{10.0, 0.1}};
    ir.geometry.push_back(primitive);
    ir.meshControl.regionTargetSizeUm["R.Thin"] = 0.5;

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    CHECK(maxCellEdgeLength(mesh.mesh) <= Catch::Approx(0.5).margin(0.05));
}

TEST_CASE("gmsh generator applies gaussian-driven refinement points", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(primitive);

    DopingProfileIr gauss;
    gauss.name = "GP.N";
    gauss.kind = DopingProfileKind::Gaussian;
    gauss.targetRegion = "R.Si";
    gauss.gaussianPeak_cm3 = 1.0e19;
    gauss.gaussianValueAtDepth_cm3 = 1.0e16;
    gauss.gaussianPeakPosUm = Point2D{1.0, 0.0};
    gauss.gaussianSigmaXUm = 0.2;
    ir.dopingProfiles.push_back(gauss);
    ir.meshControl.refineByDopingGradient = true;
    ir.meshControl.minSizeUm = 0.05;
    ir.meshControl.dopingGradientThresholdCm3PerUm = 1.0e18;

    GmshMeshGenerator generator;
    const MeshBundle2D refinedMesh = generator.generate(ir);
    CHECK(refinedMesh.mesh.numNodes() > 4);
    CHECK(refinedMesh.mesh.numCells() > 2);

    DeviceIr2D suppressed = ir;
    suppressed.meshControl.dopingGradientThresholdCm3PerUm = 1.0e25;
    const MeshBundle2D suppressedMesh = generator.generate(suppressed);
    CHECK(refinedMesh.mesh.numNodes() > suppressedMesh.mesh.numNodes());
    CHECK(refinedMesh.mesh.numCells() > suppressedMesh.mesh.numCells());
}

TEST_CASE("gmsh generator refinement probes follow gaussian peak x location", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(primitive);

    DopingProfileIr gauss;
    gauss.name = "GP.Align";
    gauss.kind = DopingProfileKind::Gaussian;
    gauss.targetRegion = "R.Si";
    gauss.gaussianPeak_cm3 = 1.0e19;
    gauss.gaussianValueAtDepth_cm3 = 1.0e16;
    gauss.gaussianPeakPosUm = Point2D{0.4, 0.0};
    gauss.gaussianSigmaXUm = 0.15;
    ir.dopingProfiles.push_back(gauss);
    ir.meshControl.refineByDopingGradient = true;
    ir.meshControl.minSizeUm = 0.05;
    ir.meshControl.dopingGradientThresholdCm3PerUm = 1.0e18;

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    bool foundPeakAlignedNode = false;
    for (Index i = 0; i < mesh.mesh.numNodes(); ++i) {
        const auto& node = mesh.mesh.getNode(i);
        if (std::abs(node.x - 0.4) < 1e-12) {
            foundPeakAlignedNode = true;
            break;
        }
    }
    CHECK(foundPeakAlignedNode);
}

TEST_CASE("gmsh generator builds conformal interface for equal target sizes", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});

    GeometryPrimitiveIr left;
    left.kind = GeometryPrimitiveKind::Rectangle;
    left.name = "rect_0";
    left.region = "R.Left";
    left.material = "Silicon";
    left.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(left);

    GeometryPrimitiveIr right = left;
    right.name = "rect_1";
    right.region = "R.Right";
    right.points = {Point2D{1.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(right);

    ir.meshControl.regionTargetSizeUm["R.Left"] = 0.25;
    ir.meshControl.regionTargetSizeUm["R.Right"] = 0.25;

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    const auto interfaceNodes = nodesOnVerticalLine(mesh.mesh, 1.0);
    CHECK(interfaceNodes.size() >= 5);
    for (Index i = 1; i < interfaceNodes.size(); ++i) {
        CHECK(mesh.mesh.getNode(interfaceNodes[i]).y > mesh.mesh.getNode(interfaceNodes[i - 1]).y);
    }
    for (Index nodeId : interfaceNodes) {
        bool usedByLeft = false;
        bool usedByRight = false;
        for (const auto& cell : mesh.mesh.cells()) {
            if (!cellUsesNode(cell, nodeId)) {
                continue;
            }
            if (cell.region_id == 0) {
                usedByLeft = true;
            }
            if (cell.region_id == 1) {
                usedByRight = true;
            }
        }
        CHECK(usedByLeft);
        CHECK(usedByRight);
    }
}

TEST_CASE("gmsh generator builds conformal interface for mixed target sizes", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});

    GeometryPrimitiveIr left;
    left.kind = GeometryPrimitiveKind::Rectangle;
    left.name = "rect_0";
    left.region = "R.Left";
    left.material = "Silicon";
    left.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(left);

    GeometryPrimitiveIr right = left;
    right.name = "rect_1";
    right.region = "R.Right";
    right.points = {Point2D{1.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(right);

    ir.meshControl.regionTargetSizeUm["R.Left"] = 0.5;
    ir.meshControl.regionTargetSizeUm["R.Right"] = 0.25;

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    const auto interfaceNodes = nodesOnVerticalLine(mesh.mesh, 1.0);
    CHECK(interfaceNodes.size() >= 5);
    for (Index nodeId : interfaceNodes) {
        bool usedByLeft = false;
        bool usedByRight = false;
        for (const auto& cell : mesh.mesh.cells()) {
            if (!cellUsesNode(cell, nodeId)) {
                continue;
            }
            if (cell.region_id == 0) {
                usedByLeft = true;
            }
            if (cell.region_id == 1) {
                usedByRight = true;
            }
        }
        CHECK(usedByLeft);
        CHECK(usedByRight);
    }
}

TEST_CASE("gmsh generator refined contact boundary includes all boundary nodes and round-trips", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);
    ir.contacts.push_back(ContactIr{"anode", "R.Si", {"rect_0#edge0"}});
    ir.meshControl.regionTargetSizeUm["R.Si"] = 0.25;

    GmshMeshGenerator generator;
    const MeshBundle2D bundle = generator.generate(ir);
    REQUIRE(bundle.mesh.numContacts() == 1);
    const auto& contact = bundle.mesh.getContact(0);
    CHECK(contact.node_ids.size() > 2);

    std::set<Index> expectedBoundaryNodes;
    for (const auto& node : bundle.mesh.nodes()) {
        if (std::abs(node.y - 0.0) < 1e-12 && node.x >= -1e-12 && node.x <= 1.0 + 1e-12) {
            expectedBoundaryNodes.insert(node.id);
        }
    }
    CHECK(std::set<Index>(contact.node_ids.begin(), contact.node_ids.end()) == expectedBoundaryNodes);
    for (const auto& edge : bundle.boundaryEdges) {
        CHECK(meshHasEdge(bundle.mesh, edge.node0, edge.node1));
    }

    const auto outDir = std::filesystem::temp_directory_path() / "vela_gmsh_contact_roundtrip";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    NeutralMeshWriter::write(bundle, outDir, nullptr);
    NeutralMeshReader reader;
    const DeviceMesh roundTrip = reader.readDirectory(outDir, UnitScalingConfig{});
    REQUIRE(roundTrip.numContacts() == 1);
    CHECK(std::set<Index>(roundTrip.getContact(0).node_ids.begin(), roundTrip.getContact(0).node_ids.end()) ==
          std::set<Index>(contact.node_ids.begin(), contact.node_ids.end()));
    std::filesystem::remove_all(outDir, ec);
}

TEST_CASE("gmsh generator conforms interface for non-divisible 0.3 and 0.2 targets", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "left", "R.Left", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "right", "R.Right", "Silicon",
        {Point2D{1.0, 0.0}, Point2D{2.0, 1.0}}});
    ir.meshControl.regionTargetSizeUm["R.Left"] = 0.3;
    ir.meshControl.regionTargetSizeUm["R.Right"] = 0.2;

    const MeshBundle2D bundle = GmshMeshGenerator{}.generate(ir);
    checkConformalVerticalInterface(bundle, 1.0, 0.0, 1.0, 0, 1);
}

TEST_CASE("gmsh generator conforms long-short reversed interface without target size", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "left", "R.Left", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 2.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Polygon, "right", "R.Right", "Silicon",
        {Point2D{1.0, 0.0}, Point2D{2.0, 0.0}, Point2D{2.0, 2.0},
         Point2D{1.0, 2.0}, Point2D{1.0, 1.0}}});

    const MeshBundle2D bundle = GmshMeshGenerator{}.generate(ir);
    checkConformalVerticalInterface(bundle, 1.0, 0.0, 2.0, 0, 1);
}

TEST_CASE("gmsh generator conforms only positive partial overlap", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "left", "R.Left", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 2.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "right", "R.Right", "Silicon",
        {Point2D{1.0, 0.5}, Point2D{2.0, 1.5}}});

    const MeshBundle2D bundle = GmshMeshGenerator{}.generate(ir);
    checkConformalVerticalInterface(bundle, 1.0, 0.5, 1.5, 0, 1);
    for (const auto& edge : bundle.boundaryEdges) {
        const auto& a = bundle.mesh.getNode(edge.node0);
        const auto& b = bundle.mesh.getNode(edge.node1);
        const bool liesOnInterface =
            std::abs(a.x - 1.0) < 1e-10 && std::abs(b.x - 1.0) < 1e-10 &&
            std::min(a.y, b.y) >= 0.5 - 1e-10 && std::max(a.y, b.y) <= 1.5 + 1e-10;
        CHECK_FALSE(liesOnInterface);
    }
}

TEST_CASE("gmsh generator T point contact creates no zero-length interface", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.A", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.B", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "a", "R.A", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "b", "R.B", "Silicon",
        {Point2D{1.0, 1.0}, Point2D{2.0, 2.0}}});

    const MeshBundle2D bundle = GmshMeshGenerator{}.generate(ir);
    for (const auto& edge : bundle.boundaryEdges) {
        CHECK(edge.node0 != edge.node1);
    }
}
