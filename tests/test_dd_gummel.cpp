#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/core/PhysicalConstants.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <unordered_map>

using namespace vela;

// ---------------------------------------------------------------------------
// Build a simple PN junction mesh (same geometry as test_poisson.cpp)
//
//   3 -------- 2
//   | p-reg  / |
//   | (T1)  /  |
//   |      /   |
//   |  T0 /    |
//   | n-reg    |
//   0 -------- 1
//
//  Nodes: 0=(0,0), 1=(L,0), 2=(L,L), 3=(0,L)   L = 1 um
//  Cells: T0={0,1,2} n-region,  T1={0,2,3} p-region
//  Contacts:
//    cathode (n): nodes 1, 2   V = 0
//    anode   (p): nodes 0, 3   V = 0
// ---------------------------------------------------------------------------

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id=0; n0.x=0;  n0.y=0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L;  n1.y=0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L;  n2.y=L;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0;  n3.y=L;  mesh.addNode(n3);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0;
    c0.node_ids = {0, 1, 2};  mesh.addCell(c0);

    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1;
    c1.node_ids = {0, 2, 3};  mesh.addCell(c1);

    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0};
    mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1};
    mesh.addRegion(r1);

    Contact anode;   anode.id=0;   anode.name="anode";
    anode.region_id=1; anode.node_ids={0,3};
    mesh.addContact(anode);

    Contact cathode; cathode.id=1; cathode.name="cathode";
    cathode.region_id=0; cathode.node_ids={1,2};
    mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        { "n_region", 1.0e23, 0.0   },
        { "p_region", 0.0,    1.0e23 }
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("GummelSolver: equilibrium (0 V bias) does not crash", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter = 30;
    cfg.reltol  = 1.0e-6;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
}

TEST_CASE("GummelSolver: n and p are strictly positive", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(sol.n(i) > 0.0);
        REQUIRE(sol.p(i) > 0.0);
    }
}

TEST_CASE("GummelSolver: no NaN or Inf in any output vector", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.psi (i)));
        REQUIRE(std::isfinite(sol.phin(i)));
        REQUIRE(std::isfinite(sol.phip(i)));
        REQUIRE(std::isfinite(sol.n   (i)));
        REQUIRE(std::isfinite(sol.p   (i)));
    }
}

TEST_CASE("GummelSolver: VTK output is written successfully", "[gummel][vtk]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    const std::string vtkPath =
        (std::filesystem::temp_directory_path() / "test_dd_gummel.vtk").string();

    REQUIRE_NOTHROW(writeDDSolutionVTK(vtkPath, mesh, doping, sol));

    REQUIRE(std::filesystem::exists(vtkPath));
    REQUIRE(std::filesystem::file_size(vtkPath) > 0);

    std::ifstream ifs(vtkPath);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    REQUIRE(content.find("Potential")          != std::string::npos);
    REQUIRE(content.find("Electrons")          != std::string::npos);
    REQUIRE(content.find("Holes")              != std::string::npos);
    REQUIRE(content.find("NetDoping")          != std::string::npos);
    REQUIRE(content.find("ElectronQuasiFermi") != std::string::npos);
    REQUIRE(content.find("HoleQuasiFermi")     != std::string::npos);
}

TEST_CASE("GummelSolver: forward bias converges without crash", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    // Small forward bias: anode = 0.3 V, cathode = 0 V
    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.3},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter    = 50;
    cfg.reltol     = 1.0e-5;
    cfg.dampingPsi = 0.5;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}

TEST_CASE("GummelSolver: abstol can terminate strongly damped high-doping updates", "[gummel]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e24}
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    std::unordered_map<std::string, Real> biases = {
        {"anode", 0.25},
        {"cathode", 0.0}
    };

    GummelConfig noAbs;
    noAbs.maxIter = 2;
    noAbs.reltol = 0.0;
    noAbs.abstol = 0.0;
    noAbs.dampingPsi = 0.05;

    const DDSolution exhausted = runGummel(mesh, matdb, doping, biases, noAbs);
    REQUIRE_FALSE(exhausted.converged);
    REQUIRE(exhausted.iters == noAbs.maxIter);

    GummelConfig withAbs = noAbs;
    withAbs.abstol = 1.0e40;

    const DDSolution converged = runGummel(mesh, matdb, doping, biases, withAbs);
    REQUIRE(converged.converged);
    REQUIRE(converged.iters == 1);

    const GummelConfig parsed = gummelConfigFromJson(nlohmann::json{{"abstol", 1.0e-12}});
    REQUIRE(parsed.abstol == Catch::Approx(1.0e-12));
}
