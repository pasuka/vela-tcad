#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/post/ElectricFieldDiagnostics.h"
#include "vela/equation/AssemblerUtils.h"

#include <cmath>

using namespace vela;

namespace {

DeviceMesh makeTriangularMesh()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 3.0, 4.0, 0.0});
    mesh.addNode(Node{2, 0.0, 4.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addRegion(Region{0, "region", "Si", {0}});
    mesh.buildEdges();
    return mesh;
}

VectorXd linearPotential(const DeviceMesh& mesh, Real a, Real b)
{
    VectorXd psi(static_cast<int>(mesh.numNodes()));
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        const Node& node = mesh.getNode(i);
        psi(static_cast<int>(i)) = a * node.x + b * node.y;
    }
    return psi;
}

} // namespace

TEST_CASE("edge electric field diagnostic matches linear potential projections", "[electric_field]")
{
    const DeviceMesh mesh = makeTriangularMesh();
    const VectorXd psi = linearPotential(mesh, 3.0, 4.0);

    // The diagonal edge from (0,0) to (3,4) is aligned with grad(psi)=(3,4),
    // so this edge-based diagnostic reaches |grad psi| = 5 V/m.
    const Real maxField = maxEdgeElectricFieldMagnitude(mesh, psi);

    REQUIRE(std::isfinite(maxField));
    REQUIRE(maxField == Catch::Approx(5.0).epsilon(1.0e-12));
}

TEST_CASE("edge electric field diagnostic reports expected horizontal and vertical projections", "[electric_field]")
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 2.0, 0.0, 0.0});
    mesh.addNode(Node{2, 0.0, 3.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.buildEdges();

    const VectorXd psi = linearPotential(mesh, 4.0, 5.0);
    const Real maxField = maxEdgeElectricFieldMagnitude(mesh, psi);

    REQUIRE(std::isfinite(maxField));
    REQUIRE(maxField == Catch::Approx(5.0).epsilon(1.0e-12));
}

TEST_CASE("edge electric field diagnostic ignores degenerate edges without NaN", "[electric_field]")
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 0.0, 0.0, 0.0});
    mesh.addNode(Node{2, 1.0, 0.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.buildEdges();

    VectorXd psi(3);
    psi << 0.0, 10.0, 2.0;

    const Real maxField = maxEdgeElectricFieldMagnitude(mesh, psi);

    REQUIRE(std::isfinite(maxField));
    REQUIRE(maxField == Catch::Approx(8.0).epsilon(1.0e-12));
}

TEST_CASE("node electric field uses inverse-distance weighted least-squares recovery", "[electric_field]")
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, -1.0, -1.0, 0.0});
    mesh.addNode(Node{2, -1.0, 0.0, 0.0});
    mesh.addNode(Node{3, 0.0, -1.0, 0.0});
    mesh.addNode(Node{4, 1.0, 0.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {1, 0, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 0, {3, 4, 0}});
    mesh.addCell(Cell{2, CellType::Tri3, 0, {0, 1, 3}});
    mesh.buildEdges();

    VectorXd psi(5);
    psi << 0.0, 5.0, 0.0, 2.0, 10.0;

    const std::vector<Real> fields = detail::computeNodeElectricFields(psi, mesh);

    REQUIRE(fields.size() == 5);
    REQUIRE(fields[0] == Catch::Approx(5.97283471086).epsilon(1.0e-12));
}
