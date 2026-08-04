#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/mesh/GmshMeshGenerator.h"

using namespace Catch::Matchers;
using namespace vela;

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
    gauss.gaussianBackground_cm3 = 1.0e16;
    gauss.gaussianCenterUm = Point2D{1.0, 0.5};
    gauss.gaussianSigmaXUm = 0.2;
    gauss.gaussianSigmaYUm = 0.2;
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
