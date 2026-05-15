#include <catch2/catch_test_macros.hpp>

#include "vela/solver/SolutionValidation.h"

#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>

using namespace vela;

namespace {

DeviceMesh makeContactedMesh()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(Node{2, 1.0, 1.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addRegion(Region{0, "semiconductor", "Si", {0}});
    mesh.addContact(Contact{0, "anode", 0, {0}});
    mesh.addContact(Contact{1, "cathode", 0, {1, 2}});
    mesh.buildEdges();
    return mesh;
}

DDSolution validSolution(Index nodeCount)
{
    const int n = static_cast<int>(nodeCount);
    DDSolution sol;
    sol.psi = VectorXd::Zero(n);
    sol.phin = VectorXd::Zero(n);
    sol.phip = VectorXd::Zero(n);
    sol.n = VectorXd::Constant(n, 1.0e-250);
    sol.p = VectorXd::Constant(n, 2.0e-250);
    sol.converged = true;
    return sol;
}

bool diagnosticsContain(const DDSolutionValidationResult& result, const std::string& needle)
{
    const std::string diagnostics = result.diagnosticsString();
    return diagnostics.find(needle) != std::string::npos;
}

} // namespace

TEST_CASE("DDSolution validation rejects negative carrier density", "[solution_validation]")
{
    const DeviceMesh mesh = makeContactedMesh();
    DDSolution sol = validSolution(mesh.numNodes());
    sol.n(1) = -1.0e-20;

    const auto result = validateDDSolution(sol, mesh, {{"anode", 0.0}, {"cathode", 0.0}});

    REQUIRE_FALSE(result.valid);
    REQUIRE(diagnosticsContain(result, "n[1]"));
}

TEST_CASE("DDSolution validation rejects NaN state entries", "[solution_validation]")
{
    const DeviceMesh mesh = makeContactedMesh();
    DDSolution sol = validSolution(mesh.numNodes());
    sol.psi(2) = std::numeric_limits<Real>::quiet_NaN();

    const auto result = validateDDSolution(sol, mesh, {{"anode", 0.0}, {"cathode", 0.0}});

    REQUIRE_FALSE(result.valid);
    REQUIRE(diagnosticsContain(result, "psi[2]"));
}

TEST_CASE("DDSolution validation accepts legal tiny positive carriers", "[solution_validation]")
{
    const DeviceMesh mesh = makeContactedMesh();
    DDSolution sol = validSolution(mesh.numNodes());

    const auto result = validateDDSolution(sol, mesh, {{"anode", 0.0}, {"cathode", 0.0}});

    REQUIRE(result.valid);
    REQUIRE(result.diagnostics.empty());
    REQUIRE(result.n.min > 0.0);
    REQUIRE(result.p.min > 0.0);
}

TEST_CASE("DDSolution validation tolerates tiny negative carrier roundoff", "[solution_validation]")
{
    const DeviceMesh mesh = makeContactedMesh();
    DDSolution sol = validSolution(mesh.numNodes());
    sol.p(0) = -1.0e-120;

    DDSolutionValidationOptions options;
    options.carrierFloor = 1.0e-100;
    const auto result = validateDDSolution(sol, mesh, {{"anode", 0.0}, {"cathode", 0.0}}, options);

    REQUIRE(result.valid);
}

TEST_CASE("DDSolution validation rejects contact quasi-Fermi mismatch", "[solution_validation]")
{
    const DeviceMesh mesh = makeContactedMesh();
    DDSolution sol = validSolution(mesh.numNodes());
    sol.phin(1) = 0.2;
    sol.phip(1) = 0.2;
    sol.phin(2) = 0.2;
    sol.phip(2) = 0.2;

    const auto result = validateDDSolution(sol, mesh, {{"anode", 0.0}, {"cathode", 0.1}});

    REQUIRE_FALSE(result.valid);
    REQUIRE(diagnosticsContain(result, "cathode"));
    REQUIRE(diagnosticsContain(result, "phin"));
}
