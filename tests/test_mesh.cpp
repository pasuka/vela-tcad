#include <catch2/catch_test_macros.hpp>
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include <cmath>

using namespace vela;

/**
 * Build a unit-square mesh made of two right triangles:
 *
 *   3 ---- 2
 *   |    / |
 *   |   /  |
 *   |  /   |
 *   | /    |
 *   0 ---- 1
 *
 * Nodes:  0=(0,0), 1=(1,0), 2=(1,1), 3=(0,1)
 * Cells:  T0={0,1,2},  T1={0,2,3}
 */
static DeviceMesh makeSquareMesh()
{
    DeviceMesh mesh;

    // Nodes
    Node n0; n0.id=0; n0.x=0.0; n0.y=0.0; mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=1.0; n1.y=0.0; mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=1.0; n2.y=1.0; mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0.0; n3.y=1.0; mesh.addNode(n3);

    // Cells
    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0;
    c0.node_ids = {0,1,2};
    mesh.addCell(c0);

    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=0;
    c1.node_ids = {0,2,3};
    mesh.addCell(c1);

    // Region
    Region r; r.id=0; r.name="Si_body"; r.material="Si";
    r.cell_ids = {0, 1};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

TEST_CASE("DeviceMesh: node and cell counts", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();
    REQUIRE(mesh.numNodes() == 4);
    REQUIRE(mesh.numCells() == 2);
}

TEST_CASE("DeviceMesh: edge count for two triangles sharing a diagonal", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();
    // Two triangles on a unit square share the diagonal (0-2).
    // Unique edges: (0,1),(1,2),(0,2),(2,3),(0,3) -> 5 edges
    REQUIRE(mesh.numEdges() == 5);
}

TEST_CASE("DeviceMesh: all edge lengths are positive", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();
    for (Index i = 0; i < mesh.numEdges(); ++i) {
        REQUIRE(mesh.getEdge(i).length > 0.0);
    }
}

TEST_CASE("DeviceMesh: axis-aligned edge lengths are 1.0", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();

    // Bottom edge (0,1): length == 1.0
    // Left edge (0,3):   length == 1.0
    // Right edge (1,2):  length == 1.0
    // Top edge (2,3):    length == 1.0
    // Diagonal (0,2):    length == sqrt(2)
    const double sqrt2 = std::sqrt(2.0);
    bool foundDiag = false;
    int  axisEdges = 0;
    for (Index i = 0; i < mesh.numEdges(); ++i) {
        Real l = mesh.getEdge(i).length;
        if (std::abs(l - 1.0) < 1e-10)    ++axisEdges;
        if (std::abs(l - sqrt2) < 1e-10)  foundDiag = true;
    }
    REQUIRE(axisEdges == 4);
    REQUIRE(foundDiag);
}

TEST_CASE("DeviceMesh: region name and material", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();
    REQUIRE(mesh.numRegions() == 1);
    REQUIRE(mesh.getRegion(0).name     == "Si_body");
    REQUIRE(mesh.getRegion(0).material == "Si");
}

TEST_CASE("DeviceMesh: getNode out of range throws", "[mesh]")
{
    DeviceMesh mesh = makeSquareMesh();
    REQUIRE_THROWS_AS(mesh.getNode(99), std::out_of_range);
}

TEST_CASE("MaterialDatabase: Si properties", "[material]")
{
    MaterialDatabase db;
    REQUIRE(db.hasMaterial("Si"));
    const Material& si = db.getMaterial("Si");
    REQUIRE(si.eps_r == 11.7);
    REQUIRE(si.ni    == 1.0e16);
    REQUIRE(si.mun   == 0.135);
    REQUIRE(si.mup   == 0.048);
}

TEST_CASE("MaterialDatabase: SiO2 properties", "[material]")
{
    MaterialDatabase db;
    REQUIRE(db.hasMaterial("SiO2"));
    REQUIRE(db.getMaterial("SiO2").eps_r == 3.9);
}

TEST_CASE("MaterialDatabase: unknown material throws", "[material]")
{
    MaterialDatabase db;
    REQUIRE_THROWS_AS(db.getMaterial("GaAs"), std::out_of_range);
}
