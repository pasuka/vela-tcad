#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/mesh/BoxGeometryBuilder.h"
#include "vela/mesh/DeviceMesh.h"

#include <cmath>

using namespace vela;

namespace {

DeviceMesh makeSingleEquilateralTriangle()
{
    DeviceMesh mesh;
    const Real h = std::sqrt(3.0) / 2.0;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5; n2.y = h; mesh.addNode(n2);

    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = {0, 1, 2};
    mesh.addCell(c);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeUnitSquareTwoTriangles()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 1.0; n2.y = 1.0; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = 1.0; mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0, 1};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeObtuseTriangle()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 2.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.2; n2.y = 0.1; mesh.addNode(n2);

    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = {0, 1, 2};
    mesh.addCell(c);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("BoxGeometryBuilder: equilateral triangle area is correct", "[box_geometry]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    const Real expectedArea = std::sqrt(3.0) / 4.0;

    REQUIRE(BoxGeometryBuilder::triangleArea(mesh.getNode(0), mesh.getNode(1), mesh.getNode(2)) ==
            Catch::Approx(expectedArea));

    Real volumeSum = 0.0;
    for (const Node& node : mesh.nodes())
        volumeSum += node.volume;
    REQUIRE(volumeSum == Catch::Approx(expectedArea));
}

TEST_CASE("BoxGeometryBuilder: square node volumes sum to total area", "[box_geometry]")
{
    DeviceMesh mesh = makeUnitSquareTwoTriangles();

    Real volumeSum = 0.0;
    for (const Node& node : mesh.nodes())
        volumeSum += node.volume;

    REQUIRE(volumeSum == Catch::Approx(1.0));
}

TEST_CASE("BoxGeometryBuilder: all edge couplings are non-negative", "[box_geometry]")
{
    DeviceMesh square = makeUnitSquareTwoTriangles();
    for (const Edge& edge : square.edges())
        REQUIRE(edge.couple >= 0.0);

    DeviceMesh obtuse = makeObtuseTriangle();
    for (const Edge& edge : obtuse.edges())
        REQUIRE(edge.couple >= 0.0);
}
