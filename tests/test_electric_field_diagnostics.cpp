#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/post/ElectricFieldDiagnostics.h"
#include "vela/equation/AssemblerUtils.h"

#include <cmath>
#include <vector>

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
VectorXd affineNodalValue(const DeviceMesh& mesh, Real a, Real b, Real c)
{
    VectorXd value(static_cast<int>(mesh.numNodes()));
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        const Node& node = mesh.getNode(i);
        value(static_cast<int>(i)) = a * node.x + b * node.y + c;
    }
    return value;
}

DeviceMesh makeFourTrianglePatch()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(Node{2, 0.0, 1.0, 0.0});
    mesh.addNode(Node{3, -1.0, 0.0, 0.0});
    mesh.addNode(Node{4, 0.0, -1.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 0, {0, 2, 3}});
    mesh.addCell(Cell{2, CellType::Tri3, 0, {0, 3, 4}});
    mesh.addCell(Cell{3, CellType::Tri3, 0, {0, 4, 1}});
    mesh.addRegion(Region{0, "silicon", "Si", {0, 1, 2, 3}});
    mesh.addContact(Contact{0, "anode", 0, {1, 2}});
    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeTwoRegionInterfacePatch()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(Node{2, 0.0, 1.0, 0.0});
    mesh.addNode(Node{3, -1.0, 0.0, 0.0});
    mesh.addNode(Node{4, 0.0, -1.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 1, {0, 3, 4}});
    mesh.addRegion(Region{0, "left_si", "Si", {0}});
    mesh.addRegion(Region{1, "right_oxide", "SiO2", {1}});
    mesh.buildEdges();
    return mesh;
}

void requireFieldApprox(const Point2& field, Real expectedX, Real expectedY)
{
    REQUIRE(field.x() == Catch::Approx(expectedX).margin(1.0e-12));
    REQUIRE(field.y() == Catch::Approx(expectedY).margin(1.0e-12));
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

TEST_CASE("Tri3 cell electric field and quasi-Fermi gradients recover affine fields", "[electric_field]")
{
    const DeviceMesh mesh = makeTriangularMesh();
    const VectorXd psi = affineNodalValue(mesh, 2.0, -3.0, 7.0);
    const VectorXd phin = affineNodalValue(mesh, -4.0, 5.0, 1.0);
    const VectorXd phip = affineNodalValue(mesh, 6.0, -8.0, -2.0);

    const auto electric = computeCellElectricField(mesh, psi);
    const auto eGrad = computeCellGradElectronQuasiFermi(mesh, phin);
    const auto hGrad = computeCellGradHoleQuasiFermi(mesh, phip);

    REQUIRE(electric.size() == 1);
    REQUIRE(electric[0].valid);
    requireFieldApprox(electric[0].vector, -2.0, 3.0);
    REQUIRE(electric[0].magnitude == Catch::Approx(std::sqrt(13.0)).epsilon(1.0e-12));

    REQUIRE(eGrad.size() == 1);
    REQUIRE(eGrad[0].valid);
    requireFieldApprox(eGrad[0].vector, -4.0, 5.0);

    REQUIRE(hGrad.size() == 1);
    REQUIRE(hGrad[0].valid);
    requireFieldApprox(hGrad[0].vector, 6.0, -8.0);
}

TEST_CASE("node electric-field recovery methods reproduce affine interior field", "[electric_field]")
{
    const DeviceMesh mesh = makeFourTrianglePatch();
    const VectorXd psi = affineNodalValue(mesh, 3.5, -2.0, 0.25);

    const auto area = computeNodeElectricFieldAreaAverage(mesh, psi);
    const auto ls1d = computeNodeElectricFieldLeastSquares(
        mesh, psi, ElectricFieldLeastSquaresWeight::InverseDistance);
    const auto ls1d2 = computeNodeElectricFieldLeastSquares(
        mesh, psi, ElectricFieldLeastSquaresWeight::InverseDistanceSquared);
    const auto spr = computeNodeElectricFieldSPR(mesh, psi);
    const auto circum1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, psi, ElectricFieldCircumcenterWeight::InverseDistance);
    const auto circumArea1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, psi, ElectricFieldCircumcenterWeight::AreaOverDistance);

    requireFieldApprox(area[0].vector, -3.5, 2.0);
    requireFieldApprox(ls1d[0].vector, -3.5, 2.0);
    requireFieldApprox(ls1d2[0].vector, -3.5, 2.0);
    requireFieldApprox(spr[0].vector, -3.5, 2.0);
    requireFieldApprox(circum1d[0].vector, -3.5, 2.0);
    requireFieldApprox(circumArea1d[0].vector, -3.5, 2.0);
}

TEST_CASE("circumcenter node recovery uses requested distance and area weights", "[electric_field]")
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(Node{2, 0.0, 1.0, 0.0});
    mesh.addNode(Node{3, -2.0, 0.0, 0.0});
    mesh.addNode(Node{4, 0.0, -1.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 0, {0, 3, 4}});
    mesh.addRegion(Region{0, "silicon", "Si", {0, 1}});
    mesh.buildEdges();

    VectorXd psi(5);
    psi << 0.0, -1.0, 0.0, 0.0, -2.0;

    const auto circum1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, psi, ElectricFieldCircumcenterWeight::InverseDistance);
    const auto circumArea1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, psi, ElectricFieldCircumcenterWeight::AreaOverDistance);

    const Real w0 = 1.0 / std::sqrt(0.5);
    const Real w1 = 1.0 / std::sqrt(1.25);
    requireFieldApprox(circum1d[0].vector,
                       w0 / (w0 + w1),
                       -2.0 * w1 / (w0 + w1));

    const Real aw0 = 0.5 / std::sqrt(0.5);
    const Real aw1 = 1.0 / std::sqrt(1.25);
    requireFieldApprox(circumArea1d[0].vector,
                       aw0 / (aw0 + aw1),
                       -2.0 * aw1 / (aw0 + aw1));
}

TEST_CASE("boundary LS and SPR fall back without exceeding area-average affine error", "[electric_field]")
{
    const DeviceMesh mesh = makeFourTrianglePatch();
    const VectorXd psi = affineNodalValue(mesh, -1.25, 4.0, 10.0);

    const auto area = computeNodeElectricFieldAreaAverage(mesh, psi);
    const auto ls1d = computeNodeElectricFieldLeastSquares(
        mesh, psi, ElectricFieldLeastSquaresWeight::InverseDistance);
    const auto spr = computeNodeElectricFieldSPR(mesh, psi);

    const Point2 expected{1.25, -4.0};
    const Real areaError = (area[1].vector - expected).norm();
    const Real lsError = (ls1d[1].vector - expected).norm();
    const Real sprError = (spr[1].vector - expected).norm();

    REQUIRE(lsError <= areaError + 1.0e-12);
    REQUIRE(sprError <= areaError + 1.0e-12);
}

TEST_CASE("region-wise recovery does not cross material interfaces", "[electric_field]")
{
    const DeviceMesh mesh = makeTwoRegionInterfacePatch();
    VectorXd psi(static_cast<int>(mesh.numNodes()));
    psi << 0.0, 2.0, 0.0, 30.0, 10.0;

    const auto area = computeNodeElectricFieldAreaAverage(mesh, psi);
    const auto ls1d = computeNodeElectricFieldLeastSquares(
        mesh, psi, ElectricFieldLeastSquaresWeight::InverseDistance);
    const auto spr = computeNodeElectricFieldSPR(mesh, psi);

    REQUIRE(area[0].regionSamples.size() == 2);
    requireFieldApprox(area[0].regionSamples.at(0).vector, -2.0, -0.0);
    requireFieldApprox(area[0].regionSamples.at(1).vector, 30.0, 10.0);

    REQUIRE(ls1d[0].regionSamples.size() == 2);
    requireFieldApprox(ls1d[0].regionSamples.at(0).vector, -2.0, -0.0);
    requireFieldApprox(ls1d[0].regionSamples.at(1).vector, 30.0, 10.0);

    REQUIRE(spr[0].regionSamples.size() == 2);
    requireFieldApprox(spr[0].regionSamples.at(0).vector, -2.0, -0.0);
    requireFieldApprox(spr[0].regionSamples.at(1).vector, 30.0, 10.0);
}
