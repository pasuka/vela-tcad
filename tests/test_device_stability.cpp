#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/LinearSolver.h"

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace vela;

namespace {

struct RegionSpec {
    std::string name;
    std::string material;
};

struct StabilityCase {
    std::string name;
    DeviceMesh mesh;
    DopingModel doping;
    std::vector<Index> oxideOnlyNodes;
};

void requireFiniteVector(const VectorXd& values)
{
    for (int i = 0; i < values.size(); ++i)
        REQUIRE(std::isfinite(values(i)));
}

void requireFiniteSparseMatrix(const SparseMatrixd& matrix)
{
    REQUIRE(matrix.rows() > 0);
    REQUIRE(matrix.cols() > 0);
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator it(matrix, outer); it; ++it)
            REQUIRE(std::isfinite(it.value()));
    }
}

std::unordered_map<Index, Real> contactBiases(const DeviceMesh& mesh, Real value = 0.0)
{
    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c) {
        for (Index node : mesh.getContact(c).node_ids)
            bcs[node] = value;
    }
    return bcs;
}

void requireScalarAssemblySystem(const SparseMatrixd& matrix,
                                 const VectorXd& rhs,
                                 Index expectedNodes)
{
    const int n = static_cast<int>(expectedNodes);
    REQUIRE(matrix.rows() == n);
    REQUIRE(matrix.cols() == n);
    REQUIRE(rhs.size() == n);
    requireFiniteSparseMatrix(matrix);
    requireFiniteVector(rhs);
}

void addNode(DeviceMesh& mesh, Index id, Real x, Real y)
{
    Node node;
    node.id = id;
    node.x = x;
    node.y = y;
    mesh.addNode(node);
}

void addTri(DeviceMesh& mesh, Index id, Index region, std::vector<Index> nodes)
{
    Cell cell;
    cell.id = id;
    cell.type = CellType::Tri3;
    cell.region_id = region;
    cell.node_ids = std::move(nodes);
    mesh.addCell(cell);
}

DeviceMesh makeHorizontalStripMesh(const std::vector<RegionSpec>& regions)
{
    DeviceMesh mesh;
    const Real dx = 1.0e-7;
    const Real h = 1.0e-7;
    const Index nx = regions.size() + 1;

    for (Index ix = 0; ix < nx; ++ix)
        addNode(mesh, ix, static_cast<Real>(ix) * dx, 0.0);
    for (Index ix = 0; ix < nx; ++ix)
        addNode(mesh, nx + ix, static_cast<Real>(ix) * dx, h);

    for (Index r = 0; r < regions.size(); ++r) {
        const Index bl = r;
        const Index br = r + 1;
        const Index tl = nx + r;
        const Index tr = nx + r + 1;
        addTri(mesh, 2 * r, r, {bl, br, tr});
        addTri(mesh, 2 * r + 1, r, {bl, tr, tl});
    }

    for (Index r = 0; r < regions.size(); ++r) {
        Region region;
        region.id = r;
        region.name = regions[r].name;
        region.material = regions[r].material;
        region.cell_ids = {2 * r, 2 * r + 1};
        mesh.addRegion(region);
    }

    Contact left;
    left.id = 0;
    left.name = "left";
    left.region_id = 0;
    left.node_ids = {0, nx};
    mesh.addContact(left);

    Contact right;
    right.id = 1;
    right.name = "right";
    right.region_id = regions.size() - 1;
    right.node_ids = {nx - 1, 2 * nx - 1};
    mesh.addContact(right);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeSiSiO2InterfaceMesh()
{
    DeviceMesh mesh;
    const Real L = 1.0e-7;

    addNode(mesh, 0, 0.0, 0.0);
    addNode(mesh, 1, L, 0.0);
    addNode(mesh, 2, 0.0, L);
    addNode(mesh, 3, L, L);
    addNode(mesh, 4, 0.0, 2.0 * L);
    addNode(mesh, 5, L, 2.0 * L);

    addTri(mesh, 0, 0, {0, 1, 3});
    addTri(mesh, 1, 0, {0, 3, 2});
    addTri(mesh, 2, 1, {2, 3, 5});
    addTri(mesh, 3, 1, {2, 5, 4});

    Region silicon;
    silicon.id = 0;
    silicon.name = "silicon";
    silicon.material = "Si";
    silicon.cell_ids = {0, 1};
    mesh.addRegion(silicon);

    Region oxide;
    oxide.id = 1;
    oxide.name = "oxide";
    oxide.material = "SiO2";
    oxide.cell_ids = {2, 3};
    mesh.addRegion(oxide);

    Contact substrate;
    substrate.id = 0;
    substrate.name = "substrate";
    substrate.region_id = 0;
    substrate.node_ids = {0, 1};
    mesh.addContact(substrate);

    Contact gate;
    gate.id = 1;
    gate.name = "gate";
    gate.region_id = 1;
    gate.node_ids = {4, 5};
    mesh.addContact(gate);

    mesh.buildEdges();
    return mesh;
}

VectorXd initialElectrons(const DeviceMesh& mesh, const DopingModel& doping)
{
    VectorXd n(static_cast<int>(mesh.numNodes()));
    for (Index i = 0; i < mesh.numNodes(); ++i)
        n(static_cast<int>(i)) = std::max(1.0e10, doping.donors(i));
    return n;
}

VectorXd initialHoles(const DeviceMesh& mesh, const DopingModel& doping)
{
    VectorXd p(static_cast<int>(mesh.numNodes()));
    for (Index i = 0; i < mesh.numNodes(); ++i)
        p(static_cast<int>(i)) = std::max(1.0e10, doping.acceptors(i));
    return p;
}

StabilityCase makeAbruptPNCase()
{
    DeviceMesh mesh = makeHorizontalStripMesh({{"n_plus", "Si"}, {"p_light", "Si"}});
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_plus", 1.0e24, 0.0},
        {"p_light", 0.0, 1.0e20},
    });
    return {"abrupt PN with 1e4 doping contrast", std::move(mesh), std::move(doping), {}};
}

StabilityCase makeCompensatedInterfaceCase()
{
    DeviceMesh mesh = makeHorizontalStripMesh({{"left_comp", "Si"}, {"right_comp", "Si"}});
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"left_comp", 1.0000e22, 0.9999e22},
        {"right_comp", 0.9999e22, 1.0000e22},
    });
    return {"nearly compensated interface", std::move(mesh), std::move(doping), {}};
}

StabilityCase makeNPlusPBodyNDriftCase()
{
    DeviceMesh mesh = makeHorizontalStripMesh({
        {"n_plus", "Si"},
        {"p_body", "Si"},
        {"n_drift", "Si"},
    });
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_plus", 1.0e24, 0.0},
        {"p_body", 0.0, 5.0e22},
        {"n_drift", 5.0e19, 0.0},
    });
    return {"n+/p-body/n-drift lateral stack", std::move(mesh), std::move(doping), {}};
}

StabilityCase makeSiSiO2InterfaceCase()
{
    DeviceMesh mesh = makeSiSiO2InterfaceMesh();
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"silicon", 0.0, 1.0e21},
        {"oxide", 0.0, 0.0},
    });
    return {"Si/SiO2 shared-edge interface", std::move(mesh), std::move(doping), {4, 5}};
}

void exercisePoissonAssembler(const StabilityCase& testCase, const MaterialDatabase& matdb)
{
    INFO(testCase.name);
    PoissonAssembler assembler(testCase.mesh, matdb, testCase.doping);
    REQUIRE_NOTHROW(assembler.assemble());
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());

    PoissonAssembler solveAssembler(testCase.mesh, matdb, testCase.doping);
    solveAssembler.assemble();
    solveAssembler.applyDirichlet(contactBiases(testCase.mesh));
    requireScalarAssemblySystem(solveAssembler.matrix(), solveAssembler.rhs(), testCase.mesh.numNodes());

    LinearSolver solver;
    VectorXd psi;
    REQUIRE_NOTHROW(psi = solver.solve(solveAssembler.matrix(), solveAssembler.rhs()));
    REQUIRE(psi.size() == static_cast<int>(testCase.mesh.numNodes()));
    requireFiniteVector(psi);
}

void exerciseDDAssembler(const StabilityCase& testCase, const MaterialDatabase& matdb)
{
    INFO(testCase.name);
    const int nNodes = static_cast<int>(testCase.mesh.numNodes());
    const VectorXd psi = VectorXd::Zero(nNodes);
    const VectorXd n = initialElectrons(testCase.mesh, testCase.doping);
    const VectorXd p = initialHoles(testCase.mesh, testCase.doping);

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};

    DDAssembler assembler(testCase.mesh,
                          matdb,
                          testCase.doping,
                          constants::Vt_300,
                          MobilityModelConfig{},
                          recombination);

    REQUIRE_NOTHROW(assembler.assemblePoissonWithCarriers(n, p, psi));
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());

    REQUIRE_NOTHROW(assembler.assembleElectronContinuity(psi, n, p));
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());
    for (Index node : testCase.oxideOnlyNodes) {
        REQUIRE(assembler.matrix().coeff(static_cast<int>(node), static_cast<int>(node)) == Catch::Approx(1.0));
        REQUIRE(assembler.rhs()(static_cast<int>(node)) == Catch::Approx(0.0));
    }

    REQUIRE_NOTHROW(assembler.assembleHoleContinuity(psi, n, p));
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());
    for (Index node : testCase.oxideOnlyNodes) {
        REQUIRE(assembler.matrix().coeff(static_cast<int>(node), static_cast<int>(node)) == Catch::Approx(1.0));
        REQUIRE(assembler.rhs()(static_cast<int>(node)) == Catch::Approx(0.0));
    }
}

void exerciseCoupledDDAssembler(const StabilityCase& testCase, const MaterialDatabase& matdb)
{
    INFO(testCase.name);
    const int nNodes = static_cast<int>(testCase.mesh.numNodes());

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};

    CoupledDDAssembler assembler(testCase.mesh,
                                 matdb,
                                 testCase.doping,
                                 constants::Vt_300,
                                 MobilityModelConfig{},
                                 recombination);

    CoupledDDState state;
    state.psi = VectorXd::Zero(nNodes);
    state.phin = VectorXd::Zero(nNodes);
    state.phip = VectorXd::Zero(nNodes);

    VectorXd x;
    REQUIRE_NOTHROW(x = assembler.pack(state));
    REQUIRE(x.size() == 3 * nNodes);

    CoupledDDBoundaryConditions bcs;
    bcs.psi = contactBiases(testCase.mesh);
    bcs.phin = contactBiases(testCase.mesh);
    bcs.phip = contactBiases(testCase.mesh);

    VectorXd residual;
    REQUIRE_NOTHROW(residual = assembler.residual(x, bcs));
    REQUIRE(residual.size() == 3 * nNodes);
    requireFiniteVector(residual);

    SparseMatrixd jacobian;
    REQUIRE_NOTHROW(jacobian = assembler.assembleJacobian(x, bcs));
    REQUIRE(jacobian.rows() == 3 * nNodes);
    REQUIRE(jacobian.cols() == 3 * nNodes);
    requireFiniteSparseMatrix(jacobian);

    const VectorXd electrons = assembler.electronDensity(x);
    const VectorXd holes = assembler.holeDensity(x);
    REQUIRE(electrons.size() == nNodes);
    REQUIRE(holes.size() == nNodes);
    for (int i = 0; i < nNodes; ++i) {
        REQUIRE(std::isfinite(electrons(i)));
        REQUIRE(std::isfinite(holes(i)));
        REQUIRE(electrons(i) >= 0.0);
        REQUIRE(holes(i) >= 0.0);
    }
}

} // namespace

TEST_CASE("small pathological device assemblies stay finite", "[device_stability][poisson][dd]")
{
    MaterialDatabase matdb;
    std::vector<StabilityCase> cases;
    cases.push_back(makeAbruptPNCase());
    cases.push_back(makeCompensatedInterfaceCase());
    cases.push_back(makeNPlusPBodyNDriftCase());
    cases.push_back(makeSiSiO2InterfaceCase());

    for (const StabilityCase& testCase : cases) {
        exercisePoissonAssembler(testCase, matdb);
        exerciseDDAssembler(testCase, matdb);
        exerciseCoupledDDAssembler(testCase, matdb);
    }
}

TEST_CASE("zero-mobility oxide edges do not create carrier transport rows", "[device_stability][dd]")
{
    const StabilityCase testCase = makeSiSiO2InterfaceCase();
    MaterialDatabase matdb;

    const int nNodes = static_cast<int>(testCase.mesh.numNodes());
    const VectorXd psi = VectorXd::Zero(nNodes);
    const VectorXd n = initialElectrons(testCase.mesh, testCase.doping);
    const VectorXd p = initialHoles(testCase.mesh, testCase.doping);

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};

    DDAssembler assembler(testCase.mesh,
                          matdb,
                          testCase.doping,
                          constants::Vt_300,
                          MobilityModelConfig{},
                          recombination);

    assembler.assembleElectronContinuity(psi, n, p);
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());
    REQUIRE(assembler.matrix().coeff(4, 5) == Catch::Approx(0.0));
    REQUIRE(assembler.matrix().coeff(5, 4) == Catch::Approx(0.0));

    assembler.assembleHoleContinuity(psi, n, p);
    requireScalarAssemblySystem(assembler.matrix(), assembler.rhs(), testCase.mesh.numNodes());
    REQUIRE(assembler.matrix().coeff(4, 5) == Catch::Approx(0.0));
    REQUIRE(assembler.matrix().coeff(5, 4) == Catch::Approx(0.0));
}
