#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

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
