#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/solver/LinearSolver.h"
#include "vela/io/VTKWriter.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <vector>

using namespace vela;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Build a simple 2-D p-n junction mesh:
 *
 *   3 -------- 2
 *   |  p-reg  /|
 *   |  (T1)  / |
 *   |       /  |
 *   |      /   |
 *   |     /    |
 *   |    /  T0 |
 *   |   / n-reg|
 *   |  /       |
 *   | /        |
 *   0 -------- 1
 *
 *  Nodes:  0=(0,0), 1=(1e-6,0), 2=(1e-6,1e-6), 3=(0,1e-6)  [1 um square]
 *  Cells:  T0={0,1,2} region 0 (n-Si),  T1={0,2,3} region 1 (p-Si)
 *  Contacts:
 *    cathode (n): nodes 1, 2   bias = 0 V
 *    anode   (p): nodes 0, 3   bias = 0 V
 */
static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;

    const double L = 1.0e-6; // 1 um

    Node n0; n0.id=0; n0.x=0;  n0.y=0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L;  n1.y=0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L;  n2.y=L;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0;  n3.y=L;  mesh.addNode(n3);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0};
    mesh.addRegion(r0);

    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1};
    mesh.addRegion(r1);

    Contact anode;   anode.id=0;   anode.name="anode";   anode.region_id=1;
    anode.node_ids = {0, 3};
    mesh.addContact(anode);

    Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0;
    cathode.node_ids = {1, 2};
    mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

// Build a DopingModel for the p-n mesh:
//   n_region: Nd = 1e23 m^-3, Na = 0
//   p_region: Nd = 0,         Na = 1e23 m^-3
static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        { "n_region", 1e23, 0.0 },
        { "p_region", 0.0,  1e23 }
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("DopingModel: net doping has correct sign per region", "[doping]")
{
    DeviceMesh mesh = makePNMesh();
    DopingModel doping = makePNDoping(mesh);

    REQUIRE(doping.numNodes() == mesh.numNodes());

    // Node 1 belongs only to n_region cell -> positive net doping
    REQUIRE(doping.netDoping(1) > 0.0);

    // Node 3 belongs only to p_region cell -> negative net doping
    REQUIRE(doping.netDoping(3) < 0.0);
}

TEST_CASE("PoissonAssembler: matrix dimensions match node count", "[poisson]")
{
    DeviceMesh      mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    const Index N = mesh.numNodes();
    REQUIRE(asm_.matrix().rows() == static_cast<int>(N));
    REQUIRE(asm_.matrix().cols() == static_cast<int>(N));
    REQUIRE(asm_.rhs().size()    == static_cast<int>(N));
}

TEST_CASE("PoissonAssembler + LinearSolver: solve succeeds, no NaN", "[poisson]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    // Apply Dirichlet: both contacts at 0 V
    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& ct = mesh.getContact(c);
        for (Index nid : ct.node_ids)
            bcs[nid] = 0.0;
    }
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // Solution length must match number of nodes
    REQUIRE(psi.size() == static_cast<int>(mesh.numNodes()));

    // No NaN values
    for (int i = 0; i < psi.size(); ++i)
        REQUIRE_FALSE(std::isnan(psi(i)));
}

TEST_CASE("PoissonAssembler: Dirichlet nodes take prescribed value", "[poisson]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    // Prescribe 0 V on all boundary nodes
    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c)
        for (Index nid : mesh.getContact(c).node_ids)
            bcs[nid] = 0.0;
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // All Dirichlet nodes must return exactly the prescribed value
    for (const auto& [nid, val] : bcs)
        REQUIRE(psi(static_cast<int>(nid)) == Catch::Approx(val).margin(1e-12));
}

TEST_CASE("VTKWriter: writes file with potential field", "[poisson][vtk]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c)
        for (Index nid : mesh.getContact(c).node_ids)
            bcs[nid] = 0.0;
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // Write to a temporary file
    const std::string vtkPath =
        (std::filesystem::temp_directory_path() / "test_poisson_out.vtk").string();
    VTKWriter writer(vtkPath, mesh);
    writer.write();

    std::vector<Real> psiVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        psiVec[i] = psi(static_cast<int>(i));
    writer.addNodeScalar("potential_V", psiVec);

    // File must exist and be non-empty
    REQUIRE(std::filesystem::exists(vtkPath));
    REQUIRE(std::filesystem::file_size(vtkPath) > 0);

    // File must contain the field name
    std::ifstream ifs(vtkPath);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    REQUIRE(content.find("potential_V") != std::string::npos);
}

static DeviceMesh makeMOSCapChargeMesh()
{
    DeviceMesh mesh;

    const double L = 1.0e-6;
    const double H = 1.0e-6;

    Node n0; n0.id=0; n0.x=0; n0.y=0;     mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L; n1.y=0;     mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=0; n2.y=H;     mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=L; n3.y=H;     mesh.addNode(n3);
    Node n4; n4.id=4; n4.x=0; n4.y=2*H;   mesh.addNode(n4);
    Node n5; n5.id=5; n5.x=L; n5.y=2*H;   mesh.addNode(n5);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,3}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=0; c1.node_ids={0,3,2}; mesh.addCell(c1);
    Cell c2; c2.id=2; c2.type=CellType::Tri3; c2.region_id=1; c2.node_ids={2,3,5}; mesh.addCell(c2);
    Cell c3; c3.id=3; c3.type=CellType::Tri3; c3.region_id=1; c3.node_ids={2,5,4}; mesh.addCell(c3);

    Region silicon; silicon.id=0; silicon.name="silicon"; silicon.material="Si"; silicon.cell_ids={0,1};
    mesh.addRegion(silicon);

    Region oxide; oxide.id=1; oxide.name="oxide"; oxide.material="SiO2"; oxide.cell_ids={2,3};
    mesh.addRegion(oxide);

    Contact body; body.id=0; body.name="body"; body.region_id=0; body.node_ids={0,1};
    mesh.addContact(body);

    Contact gate; gate.id=1; gate.name="gate"; gate.region_id=1; gate.node_ids={4,5};
    mesh.addContact(gate);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makeZeroMOSCapDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"silicon", 0.0, 0.0},
        {"oxide", 0.0, 0.0}
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static VectorXd solveMOSCapChargeCase(
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges = {})
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping, fixedCharges, sheetCharges);
    asm_.assemble();

    std::unordered_map<Index, Real> bcs = {
        {0, 0.0}, {1, 0.0}, {4, 0.0}, {5, 0.0}
    };
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    return solver.solve(asm_.matrix(), asm_.rhs());
}

TEST_CASE("PoissonAssembler: zero explicit charge matches legacy RHS", "[poisson][charge]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler legacy(mesh, matdb, doping);
    legacy.assemble();

    PoissonAssembler explicitZero(
        mesh,
        matdb,
        doping,
        {RegionFixedChargeSpec{"oxide", 0.0}},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 0.0}});
    explicitZero.assemble();

    REQUIRE(explicitZero.rhs().size() == legacy.rhs().size());
    for (int i = 0; i < legacy.rhs().size(); ++i)
        REQUIRE(explicitZero.rhs()(i) == Catch::Approx(legacy.rhs()(i)).margin(1e-30));
}

TEST_CASE("PoissonAssembler: fixed charge sign shifts MOS capacitor potential", "[poisson][charge]")
{
    const VectorXd zero = solveMOSCapChargeCase({});
    const VectorXd positive = solveMOSCapChargeCase({RegionFixedChargeSpec{"oxide", 1.0e21}});
    const VectorXd negative = solveMOSCapChargeCase({RegionFixedChargeSpec{"oxide", -1.0e21}});

    const double zeroInterface = 0.5 * (zero(2) + zero(3));
    const double positiveInterface = 0.5 * (positive(2) + positive(3));
    const double negativeInterface = 0.5 * (negative(2) + negative(3));

    REQUIRE(zeroInterface == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(positiveInterface > zeroInterface);
    REQUIRE(negativeInterface < zeroInterface);
    REQUIRE(positiveInterface == Catch::Approx(-negativeInterface).epsilon(1e-12));
}

TEST_CASE("PoissonAssembler: sheet charge is split to shared interface endpoints", "[poisson][charge]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(
        mesh,
        matdb,
        doping,
        {},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 2.0e15}});
    asm_.assemble();

    const Real expectedEndpointCharge = constants::q * 2.0e15 * 1.0e-6 * 0.5;
    REQUIRE(asm_.rhs()(2) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(asm_.rhs()(3) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(asm_.rhs()(0) == Catch::Approx(0.0).margin(1e-30));
    REQUIRE(asm_.rhs()(5) == Catch::Approx(0.0).margin(1e-30));
}
