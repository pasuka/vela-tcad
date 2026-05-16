#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/equation/PoissonAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/LinearSolver.h"

#include <cmath>

using namespace vela;
using Catch::Approx;

namespace {

// Create a simple 2x2 rectangular mesh
// Layout:
//  6 --- 7 --- 8
//  |  \  |  \  |
//  3 --- 4 --- 5
//  |  \  |  \  |
//  0 --- 1 --- 2
//
// Each square is divided into 2 triangles
DeviceMesh makeRectMesh()
{
    DeviceMesh mesh;

    const Real dx = 1.0e-6;
    const Real dy = 1.0e-6;

    // Add nodes (3x3 grid)
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Node n;
            n.id = static_cast<Index>(j * 3 + i);
            n.x = i * dx;
            n.y = j * dy;
            mesh.addNode(n);
        }
    }

    // Add region
    Region r;
    r.id = 0;
    r.name = "silicon";
    r.material = "Si";
    mesh.addRegion(r);

    // Add cells (8 triangles forming 4 squares)
    auto makeCell = [](Index id, Index a, Index b, Index c) {
        Cell cell;
        cell.id = id;
        cell.type = CellType::Tri3;
        cell.region_id = 0;
        cell.node_ids = {a, b, c};
        return cell;
    };

    // Bottom-left square
    mesh.addCell(makeCell(0, 0, 1, 3));
    mesh.addCell(makeCell(1, 1, 4, 3));
    // Bottom-right square
    mesh.addCell(makeCell(2, 1, 2, 4));
    mesh.addCell(makeCell(3, 2, 5, 4));
    // Top-left square
    mesh.addCell(makeCell(4, 3, 4, 6));
    mesh.addCell(makeCell(5, 4, 7, 6));
    // Top-right square
    mesh.addCell(makeCell(6, 4, 5, 7));
    mesh.addCell(makeCell(7, 5, 8, 7));

    mesh.buildEdges();
    return mesh;
}

DopingModel makeUniformDoping(Index numNodes, Real donors, Real acceptors)
{
    DopingModel doping(numNodes);
    for (Index i = 0; i < numNodes; ++i)
        doping.setNodeDoping(i, donors, acceptors);
    return doping;
}

} // namespace

TEST_CASE("PoissonAssembler accepts Neumann boundary specs", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 1e15, 0.0);

    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{0, 1}, 0.0});

    REQUIRE_NOTHROW(PoissonAssembler(mesh, matdb, doping, {}, {}, neumannSpecs));
}

TEST_CASE("PoissonAssembler zero Neumann boundary does not change solution", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 1e15, 0.0);

    // Assemble without Neumann boundary
    PoissonAssembler assembler1(mesh, matdb, doping);
    assembler1.assemble();
    std::unordered_map<Index, Real> bcs{{0, 0.0}, {2, 1.0}};
    assembler1.applyDirichlet(bcs);
    LinearSolver solver;
    VectorXd psi1 = solver.solve(assembler1.matrix(), assembler1.rhs());

    // Assemble with zero Neumann boundary on top edge
    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, 0.0});
    PoissonAssembler assembler2(mesh, matdb, doping, {}, {}, neumannSpecs);
    assembler2.assemble();
    assembler2.applyDirichlet(bcs);
    VectorXd psi2 = solver.solve(assembler2.matrix(), assembler2.rhs());

    // Solutions should be identical
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        REQUIRE(psi1(static_cast<int>(i)) == Approx(psi2(static_cast<int>(i))).margin(1e-12));
    }
}

TEST_CASE("PoissonAssembler non-zero Neumann boundary affects RHS", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 0.0, 0.0); // Zero doping

    // Assemble without Neumann boundary
    PoissonAssembler assembler1(mesh, matdb, doping);
    assembler1.assemble();
    const VectorXd rhs1 = assembler1.rhs();

    // Assemble with non-zero Neumann boundary on top edge
    const Real displacement = 1.0e-8; // C/m^2
    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, displacement});
    PoissonAssembler assembler2(mesh, matdb, doping, {}, {}, neumannSpecs);
    assembler2.assemble();
    const VectorXd rhs2 = assembler2.rhs();

    // RHS should be different for nodes on the boundary
    REQUIRE(rhs2(6) != Approx(rhs1(6)));
    REQUIRE(rhs2(7) != Approx(rhs1(7)));
    REQUIRE(rhs2(8) != Approx(rhs1(8)));

    // RHS should be unchanged for interior nodes
    for (Index i : {1, 3, 4, 5}) {
        REQUIRE(rhs2(static_cast<int>(i)) == Approx(rhs1(static_cast<int>(i))));
    }
}

TEST_CASE("PoissonAssembler Neumann boundary RHS contribution is correct", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 0.0, 0.0);

    const Real displacement = 2.0e-8; // C/m^2
    const Real edgeLength = 1.0e-6;   // m

    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    // Top edge: nodes 6-7-8
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, displacement});

    PoissonAssembler assembler(mesh, matdb, doping, {}, {}, neumannSpecs);
    assembler.assemble();
    const VectorXd rhs = assembler.rhs();

    // Expected contribution per edge: displacement * edgeLength / 2 to each endpoint
    const Real expectedContribution = displacement * edgeLength * 0.5;

    // Node 6: one edge (6-7)
    REQUIRE(rhs(6) == Approx(expectedContribution).margin(1e-20));

    // Node 7: two edges (6-7 and 7-8)
    REQUIRE(rhs(7) == Approx(2.0 * expectedContribution).margin(1e-20));

    // Node 8: one edge (7-8)
    REQUIRE(rhs(8) == Approx(expectedContribution).margin(1e-20));
}

TEST_CASE("PoissonAssembler insulating boundary is equivalent to zero Neumann", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 1e15, 0.0);

    // Zero Neumann
    std::vector<PoissonNeumannBoundarySpec> zeroNeumann;
    zeroNeumann.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, 0.0});
    PoissonAssembler assembler1(mesh, matdb, doping, {}, {}, zeroNeumann);
    assembler1.assemble();

    // Insulating (also zero displacement)
    std::vector<PoissonNeumannBoundarySpec> insulating;
    insulating.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, 0.0});
    PoissonAssembler assembler2(mesh, matdb, doping, {}, {}, insulating);
    assembler2.assemble();

    // RHS should be identical
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        REQUIRE(assembler1.rhs()(static_cast<int>(i)) ==
                Approx(assembler2.rhs()(static_cast<int>(i))).margin(1e-15));
    }
}

TEST_CASE("PoissonAssembler rejects out-of-range Neumann boundary nodes", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 0.0, 0.0);

    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{0, 999}, 0.0}); // Node 999 doesn't exist

    PoissonAssembler assembler(mesh, matdb, doping, {}, {}, neumannSpecs);
    REQUIRE_THROWS_AS(assembler.assemble(), std::out_of_range);
}

TEST_CASE("PoissonAssembler handles multiple Neumann boundaries", "[poisson][neumann][boundary]")
{
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 0.0, 0.0);

    const Real displacement1 = 1.0e-8;
    const Real displacement2 = -0.5e-8;

    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, displacement1}); // Top
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{0, 1, 2}, displacement2}); // Bottom

    PoissonAssembler assembler(mesh, matdb, doping, {}, {}, neumannSpecs);
    REQUIRE_NOTHROW(assembler.assemble());

    // Check RHS contributions
    const Real edgeLength = 1.0e-6;
    const Real top = displacement1 * edgeLength * 0.5;
    const Real bot = displacement2 * edgeLength * 0.5;

    REQUIRE(assembler.rhs()(0) == Approx(bot).margin(1e-20));
    REQUIRE(assembler.rhs()(1) == Approx(2.0 * bot).margin(1e-20));
    REQUIRE(assembler.rhs()(2) == Approx(bot).margin(1e-20));
    REQUIRE(assembler.rhs()(6) == Approx(top).margin(1e-20));
    REQUIRE(assembler.rhs()(7) == Approx(2.0 * top).margin(1e-20));
    REQUIRE(assembler.rhs()(8) == Approx(top).margin(1e-20));
}

TEST_CASE("PoissonAssembler linear potential with Dirichlet+zero Neumann", "[poisson][neumann][boundary]")
{
    // Use a 3x3 mesh, set Dirichlet on left (x=0) at psi=0, right (x=2) at psi=1V
    // With zero doping and zero Neumann on top/bottom, expect linear psi(x) = x/2
    DeviceMesh mesh = makeRectMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeUniformDoping(mesh.numNodes(), 0.0, 0.0);

    // Zero Neumann on top and bottom
    std::vector<PoissonNeumannBoundarySpec> neumannSpecs;
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{0, 1, 2}, 0.0});
    neumannSpecs.push_back(PoissonNeumannBoundarySpec{{6, 7, 8}, 0.0});

    PoissonAssembler assembler(mesh, matdb, doping, {}, {}, neumannSpecs);
    assembler.assemble();

    // Dirichlet: left (nodes 0, 3, 6) at 0V, right (nodes 2, 5, 8) at 1V
    std::unordered_map<Index, Real> bcs;
    bcs[0] = 0.0; bcs[3] = 0.0; bcs[6] = 0.0;
    bcs[2] = 1.0; bcs[5] = 1.0; bcs[8] = 1.0;
    assembler.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(assembler.matrix(), assembler.rhs());

    // Expect linear potential: psi at x=0 is 0, at x=1um is 0.5, at x=2um is 1.0
    REQUIRE(psi(1) == Approx(0.5).margin(1e-9));
    REQUIRE(psi(4) == Approx(0.5).margin(1e-9));
    REQUIRE(psi(7) == Approx(0.5).margin(1e-9));
}
