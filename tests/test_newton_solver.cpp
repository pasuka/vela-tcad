#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/SparseLU>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;       n2.y = L;       mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0;     n3.y = L;       mesh.addNode(n3);
    Node n4; n4.id = 4; n4.x = 0.5 * L; n4.y = 0.5 * L; mesh.addNode(n4);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 4}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {1, 2, 4}; mesh.addCell(c1);
    Cell c2; c2.id = 2; c2.type = CellType::Tri3; c2.region_id = 1; c2.node_ids = {2, 3, 4}; mesh.addCell(c2);
    Cell c3; c3.id = 3; c3.type = CellType::Tri3; c3.region_id = 1; c3.node_ids = {3, 0, 4}; mesh.addCell(c3);

    Region r0; r0.id = 0; r0.name = "n_region"; r0.material = "Si"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);
    Region r1; r1.id = 1; r1.name = "p_region"; r1.material = "Si"; r1.cell_ids = {2, 3}; mesh.addRegion(r1);

    Contact cathode; cathode.id = 0; cathode.name = "cathode"; cathode.region_id = 0; cathode.node_ids = {1, 2}; mesh.addContact(cathode);
    Contact anode; anode.id = 1; anode.name = "anode"; anode.region_id = 1; anode.node_ids = {0, 3}; mesh.addContact(anode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e21, 0.0},
        {"p_region", 0.0, 1.0e21},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static DeviceMesh makeOxideMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;   n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = L;   mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);

    Region r0; r0.id = 0; r0.name = "oxide"; r0.material = "SiO2"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);

    mesh.buildEdges();
    return mesh;
}

static DeviceMesh makeContactedOxideMesh()
{
    DeviceMesh mesh = makeOxideMesh();
    Contact gate;
    gate.id = 0;
    gate.name = "gate";
    gate.region_id = 0;
    gate.node_ids = {0, 1, 2, 3};
    mesh.addContact(gate);
    return mesh;
}

static std::unordered_map<std::string, Real> zeroBias()
{
    return {{"anode", 0.0}, {"cathode", 0.0}};
}

static NewtonConfig newtonConfig()
{
    NewtonConfig cfg;
    cfg.maxIter = 10;
    cfg.reltol = 1.0e-7;
    cfg.abstol = 1.0e-20;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.verbose = false;
    return cfg;
}

TEST_CASE("NewtonSolver: PN diode equilibrium converges", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());

    REQUIRE(result.converged);
    REQUIRE(result.iters >= 0);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
}

TEST_CASE("NewtonSolver: Gummel initial guess reduces Newton iterations", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    NewtonResult fromDefault = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    GummelConfig gcfg;
    gcfg.maxIter = 8;
    gcfg.reltol = 1.0e-6;
    gcfg.dampingPsi = 0.5;
    DDSolution gummel = runGummel(mesh, matdb, doping, zeroBias(), gcfg);
    NewtonResult fromGummel = runNewton(mesh, matdb, doping, zeroBias(), gummel, cfg);

    REQUIRE(fromDefault.converged);
    REQUIRE(fromGummel.converged);
    REQUIRE(fromGummel.iters <= fromDefault.iters);
}

TEST_CASE("NewtonSolver: no NaN or Inf and carriers stay positive", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(result.converged);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(result.solution.psi(i)));
        REQUIRE(std::isfinite(result.solution.phin(i)));
        REQUIRE(std::isfinite(result.solution.phip(i)));
        REQUIRE(std::isfinite(result.solution.n(i)));
        REQUIRE(std::isfinite(result.solution.p(i)));
        REQUIRE(result.solution.n(i) > 0.0);
        REQUIRE(result.solution.p(i) > 0.0);
    }
}

TEST_CASE("CoupledDDAssembler: zero-mobility continuity rows are pinned", "[newton][coupled]")
{
    DeviceMesh mesh = makeOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::LinSpaced(N, 0.1, 0.4);
    state.phip = VectorXd::LinSpaced(N, -0.4, -0.1);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    for (Index i = 0; i < mesh.numNodes(); ++i)
        bcs.psi[i] = 0.0;

    const VectorXd r = assembler.residual(x, bcs);
    for (int i = 0; i < N; ++i) {
        REQUIRE(r(N + i) == Catch::Approx(state.phin(i)));
        REQUIRE(r(2 * N + i) == Catch::Approx(state.phip(i)));
    }

    const SparseMatrixd J = assembler.finiteDifferenceJacobian(x, bcs);
    for (int i = 0; i < N; ++i) {
        REQUIRE(J.coeff(N + i, N + i) == Catch::Approx(1.0));
        REQUIRE(J.coeff(2 * N + i, 2 * N + i) == Catch::Approx(1.0));
    }

    Eigen::SparseLU<SparseMatrixd> lu;
    lu.compute(J);
    REQUIRE(lu.info() == Eigen::Success);
}


TEST_CASE("CoupledDDAssembler: analytic pinned rows suppress zero-rate recombination derivatives", "[newton][coupled]")
{
    DeviceMesh mesh = makeOxideMesh();
    MaterialDatabase matdb;
    Material zeroMobilitySemiconductor;
    zeroMobilitySemiconductor.name = "SiO2";
    zeroMobilitySemiconductor.eps_r = 11.7;
    zeroMobilitySemiconductor.ni = 1.0e16;
    zeroMobilitySemiconductor.mun = 0.0;
    zeroMobilitySemiconductor.mup = 0.0;
    matdb.addMaterial(zeroMobilitySemiconductor);

    DopingModel doping(mesh.numNodes());
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::LinSpaced(N, 0.05, 0.2);
    state.phip = state.phin;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    const SparseMatrixd J = assembler.assembleJacobian(x, bcs);
    const Eigen::MatrixXd dense = Eigen::MatrixXd(J);

    for (int i = 0; i < N; ++i) {
        const int electronRow = N + i;
        const int holeRow = 2 * N + i;
        for (int col = 0; col < 3 * N; ++col) {
            REQUIRE(dense(electronRow, col) == Catch::Approx(col == electronRow ? 1.0 : 0.0));
            REQUIRE(dense(holeRow, col) == Catch::Approx(col == holeRow ? 1.0 : 0.0));
        }
    }
}

TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences on small mesh", "[newton][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.04, 0.05);
    state.phin = VectorXd::LinSpaced(N, 0.01, -0.015);
    state.phip = VectorXd::LinSpaced(N, -0.02, 0.012);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 5.0e-5);
}

TEST_CASE("NewtonSolver: defaults to analytic Jacobian", "[newton]")
{
    const NewtonConfig cfg;
    REQUIRE(cfg.jacobian == "analytic");

    const NewtonConfig debugCfg = newtonConfigFromJson(nlohmann::json{{"jacobian", "finite_difference"}});
    REQUIRE(debugCfg.jacobian == "finite_difference");
}



TEST_CASE("NewtonSolver: block residual norm balances mixed equation blocks", "[newton][residual]")
{
    const int N = 2;
    VectorXd initial(3 * N);
    initial << 1.0e-18, -1.0e-18,
               10.0, -10.0,
               5.0, -5.0;
    VectorXd current(3 * N);
    current << 1.0e-18, 0.0,
               1.0, -1.0,
               0.5, -0.5;

    const ResidualBlockNormValue initialBlocks = ResidualNorm::computeBlocks(initial, N);
    const ResidualBlockNormValue currentBlocks = ResidualNorm::computeBlocks(current, N);

    REQUIRE(initialBlocks.psi == Catch::Approx(std::sqrt(2.0) * 1.0e-18));
    REQUIRE(initialBlocks.phin == Catch::Approx(std::sqrt(200.0)));
    REQUIRE(initialBlocks.phip == Catch::Approx(std::sqrt(50.0)));

    const Real balanced = ResidualNorm::normalizedBlockL2(currentBlocks, initialBlocks);
    REQUIRE(balanced == Catch::Approx(std::sqrt(0.5 + 0.01 + 0.01)));

    ResidualBlockWeights continuityOnly;
    continuityOnly.psi = 0.0;
    continuityOnly.phin = 1.0;
    continuityOnly.phip = 4.0;
    const Real weighted = ResidualNorm::normalizedBlockL2(
        currentBlocks, initialBlocks, continuityOnly);
    REQUIRE(weighted == Catch::Approx(std::sqrt(0.01 + 4.0 * 0.01)));
}

TEST_CASE("NewtonSolver: parses block residual norm controls", "[newton][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"residual_norm", "block"},
        {"residual_weights", {{"psi", 0.25}, {"phin", 2.0}, {"phip", 3.0}}},
        {"residual_scales", {{"psi", 1.0e-18}, {"phin", 2.0e4}, {"phip", 3.0e4}}}
    });

    REQUIRE(cfg.residualNorm == "block");
    REQUIRE(cfg.residualWeightPsi == Catch::Approx(0.25));
    REQUIRE(cfg.residualWeightPhin == Catch::Approx(2.0));
    REQUIRE(cfg.residualWeightPhip == Catch::Approx(3.0));
    REQUIRE(cfg.residualScalePsi == Catch::Approx(1.0e-18));
    REQUIRE(cfg.residualScalePhin == Catch::Approx(2.0e4));
    REQUIRE(cfg.residualScalePhip == Catch::Approx(3.0e4));

    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"residual_norm", "unknown"}}),
        std::invalid_argument);
}

TEST_CASE("NewtonSolver: verbose false suppresses failure diagnostics", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::LinSpaced(N, -0.03, 0.04);
    initial.phin = VectorXd::Constant(N, 0.01);
    initial.phip = VectorXd::Constant(N, -0.01);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.verbose = false;

    std::ostringstream capturedStderr;
    std::streambuf* previousStderr = std::cerr.rdbuf(capturedStderr.rdbuf());
    const NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), initial, cfg);
    std::cerr.rdbuf(previousStderr);

    REQUIRE_FALSE(result.converged);
    REQUIRE(capturedStderr.str().empty());
}

TEST_CASE("NewtonSolver: line search rejection returns last accepted state", "[newton][line_search]")
{
    DeviceMesh mesh = makeContactedOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Zero(N);
    initial.phin = VectorXd::Constant(N, 0.1);
    initial.phip = VectorXd::Constant(N, -0.1);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 3;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.verbose = false;
    cfg.lineSearch = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"gate", 0.0}}, initial, cfg);

    REQUIRE_FALSE(result.converged);
    REQUIRE(result.iters == 0);
    REQUIRE(result.history.empty());
    REQUIRE(result.finalResidualNorm == Catch::Approx(result.initialResidualNorm));
    REQUIRE((result.solution.psi - initial.psi).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phin - initial.phin).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phip - initial.phip).norm() == Catch::Approx(0.0));
}
