#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/SparseLU>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
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

static DeviceMesh makePartiallyContactedOxideMesh()
{
    DeviceMesh mesh = makeOxideMesh();
    Contact gate;
    gate.id = 0;
    gate.name = "gate";
    gate.region_id = 0;
    gate.node_ids = {0};
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

static void requireFiniteNewtonSolution(const NewtonResult& result, Index nodeCount);

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

TEST_CASE("NewtonSolver: ohmic contact BC resists compensated-node polarity flips",
          "[newton][contact_bc]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    // Inject an opposite-sign outlier on one anode node to mimic imported
    // compensated/tie ownership artifacts.
    doping.setNodeDoping(0, 8.0e20, 0.0);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(result.converged);

    // Both anode nodes should keep p-side built-in sign (negative psi at 0 V).
    REQUIRE(result.solution.psi(0) < 0.0);
    REQUIRE(result.solution.psi(3) < 0.0);
}

TEST_CASE("NewtonSolver: PN diode equilibrium converges with unit_scaling state",
          "[newton][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    REQUIRE(result.converged);
    REQUIRE(result.iters >= 0);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
    requireFiniteNewtonSolution(result, mesh.numNodes());
}

TEST_CASE("NewtonSolver: high-doping unit-scaled PN cold start reaches near-zero 0V current",
          "[newton][scaling][bgn]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    NewtonConfig cfg;
    cfg.maxIter = 30;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-18;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.verbose = false;
    cfg.maxUpdate = 5.0;
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    cfg.mobility = mobilityModelConfig("caughey_thomas_field");
    cfg.recombination = {"srh", "auger"};
    cfg.bandgapNarrowing = bandgapNarrowingConfig("slotboom");

    const NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    REQUIRE(result.converged);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    ContactCurrent current(mesh, matdb, doping, cfg.mobility, cfg.temperature_K);
    const ContactCurrentResult anode = current.compute(result.solution, "anode");
    REQUIRE(std::abs(anode.totalCurrent) < 1.0e-9);
}

TEST_CASE("ContactCurrent: preserved edge hole QF drop changes reporting only",
          "[contact_current][qf_floor]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    DDSolution solution;
    const int N = static_cast<int>(mesh.numNodes());
    solution.psi = VectorXd::Constant(N, -13.2);
    solution.phin = VectorXd::Constant(N, -12.8);
    solution.phip = VectorXd::Constant(N, -12.8);
    solution.n = VectorXd::Constant(N, 1.0e10);
    solution.p = VectorXd::Constant(N, 1.0e23);

    MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    ContactCurrent current(mesh, matdb, doping, mobility, 300.0);

    const ContactCurrentDetailedResult baseline =
        current.computeDetailed(solution, "anode");
    REQUIRE(std::abs(baseline.totals.holeCurrent) < 1.0e-30);

    REQUIRE_FALSE(baseline.edges.empty());
    const Index edgeId = baseline.edges.front().edgeId;
    ContactCurrentEdgeOverrides overrides;
    const Real ulpAtBias = std::abs(std::nextafter(-12.8, -std::numeric_limits<Real>::infinity()) + 12.8);
    overrides.holeQuasiFermiDropByEdge[edgeId] = -5.0 * ulpAtBias;

    const ContactCurrentDetailedResult reported =
        current.computeDetailed(solution, "anode", overrides);

    REQUIRE(reported.totals.electronCurrent ==
            Catch::Approx(baseline.totals.electronCurrent));
    REQUIRE(std::abs(reported.totals.holeCurrent) > 1.0e-30);
    REQUIRE(reported.totals.totalCurrent ==
            Catch::Approx(reported.totals.electronCurrent - reported.totals.holeCurrent));
    REQUIRE(reported.totals.totalCurrent !=
            Catch::Approx(baseline.totals.totalCurrent));

    const auto changed = std::find_if(
        reported.edges.begin(),
        reported.edges.end(),
        [edgeId](const ContactCurrentEdgeDiagnostic& edge) {
            return edge.edgeId == edgeId;
        });
    REQUIRE(changed != reported.edges.end());
    REQUIRE(changed->holeQfDropOverrideApplied);
    REQUIRE(changed->phip1 - changed->phip0 ==
            Catch::Approx(-5.0 * ulpAtBias));
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

TEST_CASE("NewtonSolver: unit_scaling accepts a physical Gummel warm initial guess",
          "[newton][scaling][warm_start]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    GummelConfig gcfg;
    gcfg.maxIter = 8;
    gcfg.reltol = 1.0e-6;
    gcfg.dampingPsi = 0.5;
    gcfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution gummel = runGummel(mesh, matdb, doping, zeroBias(), gcfg);

    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    cfg.warmStart = true;

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), gummel, cfg);

    REQUIRE(result.converged);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
    requireFiniteNewtonSolution(result, mesh.numNodes());
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

TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences with varying intrinsic density",
          "[newton][coupled][bgn]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e21},
    });

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.025);
    state.phin = VectorXd::LinSpaced(N, 0.006, -0.008);
    state.phip = VectorXd::LinSpaced(N, -0.007, 0.005);
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

    REQUIRE(rel < 1.0e-4);
}

TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences at BV absolute potential scale",
          "[newton][coupled][bgn][bv]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e21},
    });

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi.resize(N);
    state.phin.resize(N);
    state.phip.resize(N);
    state.psi << -13.20365, -13.20340, -13.20315, -13.20390, -13.20355;
    state.phin << -12.79890, -12.79930, -12.79970, -12.80010, -12.79910;
    state.phip << -12.80020, -12.79980, -12.79940, -12.80000, -12.79960;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[1] = state.psi(1);
    bcs.phin[1] = state.phin(1);
    bcs.phip[1] = state.phip(1);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);
    bcs.psi[3] = state.psi(3);
    bcs.phin[3] = state.phin(3);
    bcs.phip[3] = state.phip(3);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-8);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 1.0e-5);
}

TEST_CASE("CoupledDDAssembler: transport Jacobian captures quasi-Fermi high-field mobility",
          "[newton][coupled][mobility][field]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{});

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.025);
    state.phin = VectorXd::LinSpaced(N, 0.7, -0.7);
    state.phip = VectorXd::LinSpaced(N, -0.65, 0.65);
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

    REQUIRE(rel < 1.0e-4);
}
TEST_CASE("CoupledDDAssembler: Slotboom BGN uses total impurity at compensated nodes",
          "[newton][coupled][bgn][doping]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(4, 1.0e23, 1.0e23);

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const Material& si = matdb.getMaterial("Si");
    const SlotboomBandgapNarrowing bgn(bandgapNarrowingConfig("slotboom"));
    const Real expected = effectiveIntrinsicDensity(
        si.ni,
        constants::Vt_300,
        bgn.deltaEg(doping.totalImpurity(4), 0.0, 0.0));

    REQUIRE(doping.netDoping(4) == Catch::Approx(0.0));
    REQUIRE(expected > si.ni);
    REQUIRE(assembler.intrinsicDensity().at(4) == Catch::Approx(expected));
}

TEST_CASE("CoupledDDAssembler: scaled state residual and Jacobian are consistent",
          "[newton][coupled][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    UnitScalingSystem::AutoInputs inputs =
        UnitScalingSystem::autoInputsFrom(mesh, doping, matdb, 1.0e16);
    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        300.0, constants::eps0 * 11.7, inputs);

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = sc.V0();
    scaling.C0 = sc.C0();
    scaling.mu0 = sc.mu0();
    scaling.D0 = sc.D0();
    scaling.L0 = sc.L0();
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{},
        {},
        {},
        scaling);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.04, 0.05) / scaling.V0;
    state.phin = VectorXd::LinSpaced(N, 0.01, -0.015) / scaling.V0;
    state.phip = VectorXd::LinSpaced(N, -0.02, 0.012) / scaling.V0;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const VectorXd r = assembler.residual(x, bcs);
    REQUIRE(r.allFinite());
    REQUIRE(r.norm() < 1.0e8);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());
    REQUIRE(rel < 1.0e-8);
}

TEST_CASE("NewtonSolver: evaluateStep reports one physical Newton correction", "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonStepEvaluation step = solver.evaluateStep(state);

    REQUIRE(step.residual.raw.size() == 3 * N);
    REQUIRE(step.deltaPsi.size() == N);
    REQUIRE(step.deltaPhin.size() == N);
    REQUIRE(step.deltaPhip.size() == N);
    REQUIRE(step.trialSolution.psi.size() == N);
    REQUIRE(step.stepNorm > 0.0);
    REQUIRE(step.rawStepNorm > 0.0);
    REQUIRE(step.trialSolution.psi(0) ==
            Catch::Approx(state.psi(0) + step.deltaPsi(0)));
    REQUIRE(step.trialSolution.phin(1) ==
            Catch::Approx(state.phin(1) + step.deltaPhin(1)));
    REQUIRE(step.trialSolution.phip(2) ==
            Catch::Approx(state.phip(2) + step.deltaPhip(2)));
    REQUIRE(step.trialResidual.blockNorms.combined < step.residual.blockNorms.combined);
}

TEST_CASE("NewtonSolver: evaluateDirectionalDerivative compares analytic and finite-difference Jv",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    DDSolution perturbation;
    perturbation.psi = VectorXd::Zero(N);
    perturbation.phin = VectorXd::Zero(N);
    perturbation.phip = VectorXd::Zero(N);
    perturbation.psi(4) = 0.5e-6;
    perturbation.phin(4) = -0.5e-6;

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonDirectionalDerivativeEvaluation jvp =
        solver.evaluateDirectionalDerivative(state, perturbation);

    REQUIRE(jvp.residual.raw.size() == 3 * N);
    REQUIRE(jvp.analyticJv.size() == 3 * N);
    REQUIRE(jvp.finiteDifferenceJv.size() == 3 * N);
    REQUIRE(jvp.perturbationPsi.size() == N);
    REQUIRE(jvp.perturbationPhin.size() == N);
    REQUIRE(jvp.perturbationPhip.size() == N);
    REQUIRE(jvp.analyticJv.norm() > 0.0);
    REQUIRE(jvp.finiteDifferenceJv.norm() > 0.0);
    REQUIRE(jvp.relativeError < 1.0e-6);
    REQUIRE(jvp.perturbationPsi(4) == Catch::Approx(0.5e-6));
    REQUIRE(jvp.perturbationPhin(4) == Catch::Approx(-0.5e-6));
}

TEST_CASE("NewtonSolver: evaluateJacobianBlockAudit reports finite block rows",
          "[newton][diagnostics][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;
    cfg.impactIonization.model = "van_overstraeten";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "density_gradient";

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::LinSpaced(N, -0.015, 0.015);
    state.phip = VectorXd::LinSpaced(N, 0.012, -0.012);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(state, 1.0e-7);
    const auto hasBlock = [&](const std::string& name) {
        return std::any_of(rows.begin(), rows.end(), [&](const auto& row) {
            return row.block == name &&
                   std::isfinite(row.analyticNorm) &&
                   std::isfinite(row.fdNorm) &&
                   std::isfinite(row.diffNorm) &&
                   std::isfinite(row.relDiff);
        });
    };

    REQUIRE(hasBlock("poisson"));
    REQUIRE(hasBlock("transport"));
    REQUIRE(hasBlock("srh_auger"));
    REQUIRE(hasBlock("sg_avalanche"));
    REQUIRE(hasBlock("dirichlet_or_gauge"));
}

TEST_CASE("NewtonSolver: evaluateJacobianBlockAudit can restrict expensive block rows",
          "[newton][diagnostics][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;
    cfg.impactIonization.model = "van_overstraeten";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "density_gradient";

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::LinSpaced(N, -0.015, 0.015);
    state.phip = VectorXd::LinSpaced(N, 0.012, -0.012);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(std::isfinite(rows.front().analyticNorm));
    REQUIRE(std::isfinite(rows.front().fdNorm));
    REQUIRE(std::isfinite(rows.front().diffNorm));
    REQUIRE(std::isfinite(rows.front().relDiff));
}

TEST_CASE("NewtonSolver: evaluateBlockStep freezes complementary unknown blocks",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonBlockStepEvaluation poisson =
        solver.evaluateBlockStep(state, "poisson_only");
    const NewtonBlockStepEvaluation carriers =
        solver.evaluateBlockStep(state, "carrier_only");

    REQUIRE(poisson.mode == "poisson_only");
    REQUIRE(carriers.mode == "carrier_only");
    REQUIRE(poisson.residual.raw.size() == 3 * N);
    REQUIRE(carriers.residual.raw.size() == 3 * N);
    REQUIRE(poisson.deltaPsi.size() == N);
    REQUIRE(carriers.deltaPhin.size() == N);
    REQUIRE(poisson.deltaPsi.norm() > 0.0);
    REQUIRE(poisson.deltaPhin.norm() == Catch::Approx(0.0));
    REQUIRE(poisson.deltaPhip.norm() == Catch::Approx(0.0));
    REQUIRE(carriers.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE(carriers.deltaPhin.norm() + carriers.deltaPhip.norm() > 0.0);
    REQUIRE(poisson.trialSolution.psi(0) ==
            Catch::Approx(state.psi(0) + poisson.deltaPsi(0)));
    REQUIRE(carriers.trialSolution.phin(1) ==
            Catch::Approx(state.phin(1) + carriers.deltaPhin(1)));
}

TEST_CASE("NewtonSolver: evaluateRegularizedCarrierStep damps carrier-only correction",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxUpdate = 0.0;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonBlockStepEvaluation baseline =
        solver.evaluateBlockStep(state, "carrier_only");
    const NewtonRegularizedCarrierStepEvaluation zero =
        solver.evaluateRegularizedCarrierStep(state, 0.0);
    const NewtonRegularizedCarrierStepEvaluation damped =
        solver.evaluateRegularizedCarrierStep(state, 10.0);

    REQUIRE(zero.regularizationScale == Catch::Approx(0.0));
    REQUIRE(damped.regularizationScale == Catch::Approx(10.0));
    REQUIRE(zero.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE(damped.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE((zero.deltaPhin - baseline.deltaPhin).norm() ==
            Catch::Approx(0.0).margin(1.0e-12));
    REQUIRE((zero.deltaPhip - baseline.deltaPhip).norm() ==
            Catch::Approx(0.0).margin(1.0e-12));
    REQUIRE(damped.rawStepNorm < baseline.rawStepNorm);
}

TEST_CASE("NewtonSolver: evaluateCarrierRowDiagnostics reports carrier row stiffness",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxUpdate = 1.0;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonCarrierRowDiagnosticsEvaluation rows =
        solver.evaluateCarrierRowDiagnostics(state);

    REQUIRE(rows.rows.size() == static_cast<std::size_t>(N));
    REQUIRE(rows.rawCarrierStepNorm > 0.0);
    REQUIRE(rows.cappedCarrierStepNorm > 0.0);
    REQUIRE(rows.rows[0].nodeId == 0);
    REQUIRE(rows.rows[0].electronRowAbsSum >= std::abs(rows.rows[0].electronDiagonal));
    REQUIRE(rows.rows[0].holeRowAbsSum >= std::abs(rows.rows[0].holeDiagonal));
    REQUIRE(rows.rows[0].electronRowL2Norm >= 0.0);
    REQUIRE(rows.rows[0].holeRowL2Norm >= 0.0);
    REQUIRE(rows.rows[0].rawDeltaPhin_V != Catch::Approx(0.0));
    REQUIRE(std::abs(rows.rows[0].cappedDeltaPhin_V) <=
            cfg.maxUpdate * rows.potentialScale + 1.0e-12);
}

TEST_CASE("NewtonSolver: evaluateCarrierTermDiagnostics decomposes continuity residual",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonCarrierTermDiagnosticsEvaluation terms =
        solver.evaluateCarrierTermDiagnostics(state);

    REQUIRE(terms.rows.size() == static_cast<std::size_t>(N));
    const auto& center = terms.rows[4];
    REQUIRE(center.nodeId == 4);
    REQUIRE(center.electronBoundary == Catch::Approx(0.0));
    const Real electronSum = center.electronFlux
        + center.electronRecombination
        + center.electronImpact
        + center.electronGauge
        + center.electronBoundary;
    const Real holeSum = center.holeFlux
        + center.holeRecombination
        + center.holeImpact
        + center.holeGauge
        + center.holeBoundary;
    REQUIRE(electronSum == Catch::Approx(center.electronResidual).margin(1.0e-18));
    REQUIRE(holeSum == Catch::Approx(center.holeResidual).margin(1.0e-18));
}

TEST_CASE("NewtonSolver: defaults to analytic Jacobian", "[newton]")
{
    const NewtonConfig cfg;
    REQUIRE(cfg.jacobian == "analytic");
    REQUIRE_FALSE(cfg.warmStart);
    REQUIRE(cfg.quasiFermiUpdateLimit_V == Catch::Approx(0.0));
    REQUIRE(cfg.carrierRegularizationScale == Catch::Approx(0.0));
    REQUIRE(cfg.contactBoundaryReconstruction == "dominant_signed_contact_mean");
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxation);
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V ==
            Catch::Approx(0.1));
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
            REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationContactSide == "p_contact_only");

    const NewtonConfig debugCfg = newtonConfigFromJson(nlohmann::json{
        {"jacobian", "finite_difference"},
        {"warm_start", true},
        {"quasi_fermi_update_limit_V", 0.0259},
        {"carrier_regularization_scale", 3.0},
        {"contact_boundary_reconstruction", "legacy_node_local"},
    });
    REQUIRE(debugCfg.jacobian == "finite_difference");
    REQUIRE(debugCfg.warmStart);
    REQUIRE(debugCfg.quasiFermiUpdateLimit_V == Catch::Approx(0.0259));
    REQUIRE(debugCfg.carrierRegularizationScale == Catch::Approx(3.0));
    REQUIRE(debugCfg.contactBoundaryReconstruction == "legacy_node_local");
}

TEST_CASE("NewtonSolver: unit_scaling config records scaled mode and preserves analytic Jacobian",
          "[newton][scaling][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(
        nlohmann::json{{"jacobian", "analytic"}},
        UnitScalingConfig{UnitScalingMode::UnitScaling});

    REQUIRE(cfg.inputScaling.isUnitScaling());
    REQUIRE(cfg.jacobian == "analytic");
}

TEST_CASE("NewtonSolver: warm start preserves supplied quasi-Fermi guess", "[newton][warm_start]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Zero(N);
    initial.phin = VectorXd::Zero(N);
    initial.phip = VectorXd::Zero(N);
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.02;
    initial.phip(interiorNode) = -0.015;

    NewtonConfig coldCfg = newtonConfig();
    coldCfg.maxIter = 0;
    coldCfg.reltol = 0.0;
    coldCfg.abstol = 0.0;
    coldCfg.warmStart = false;

    NewtonConfig warmCfg = coldCfg;
    warmCfg.warmStart = true;

    const NewtonResult cold = runNewton(mesh, matdb, doping, zeroBias(), initial, coldCfg);
    const NewtonResult warm = runNewton(mesh, matdb, doping, zeroBias(), initial, warmCfg);

    REQUIRE_FALSE(warm.converged);

    REQUIRE(cold.solution.phin(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(cold.solution.phip(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(warm.solution.phin(interiorNode) ==
            Catch::Approx(initial.phin(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.solution.phip(interiorNode) ==
            Catch::Approx(initial.phip(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.initialResidualNorm != Catch::Approx(cold.initialResidualNorm));
}

TEST_CASE("NewtonSolver: warm start projects contact nodes to the current bias",
          "[newton][warm_start][contacts]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution stale;
    stale.psi = VectorXd::Zero(N);
    stale.phin = VectorXd::Zero(N);
    stale.phip = VectorXd::Zero(N);
    stale.n = VectorXd::Ones(N);
    stale.p = VectorXd::Ones(N);

    const int interiorNode = 4;
    stale.phin(interiorNode) = 0.02;
    stale.phip(interiorNode) = -0.015;

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;

    const Real anodeBias = -0.25;
    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"anode", anodeBias}, {"cathode", 0.0}},
        stale,
        cfg);

    REQUIRE_FALSE(result.converged);
    REQUIRE(result.solution.phin(interiorNode) ==
            Catch::Approx(stale.phin(interiorNode)).margin(1.0e-14));
    REQUIRE(result.solution.phip(interiorNode) ==
            Catch::Approx(stale.phip(interiorNode)).margin(1.0e-14));
    for (Index nid : mesh.getContact(1).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phip(ii) == Catch::Approx(anodeBias).margin(1.0e-14));
    }
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

TEST_CASE("NewtonSolver: evaluates residual for an externally supplied state",
          "[newton][residual]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    NewtonSolver solver(mesh, matdb, doping, zeroBias(), cfg);

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution state;
    state.psi = VectorXd::LinSpaced(n, -0.02, 0.03);
    state.phin = VectorXd::Constant(n, 0.004);
    state.phip = VectorXd::Constant(n, -0.003);
    state.n = VectorXd::Zero(n);
    state.p = VectorXd::Zero(n);

    const NewtonResidualEvaluation residual = solver.evaluateResidual(state);

    REQUIRE(residual.raw.size() == 3 * n);
    REQUIRE(residual.blockNorms.psi > 0.0);
    REQUIRE(residual.blockNorms.phin > 0.0);
    REQUIRE(residual.blockNorms.phip > 0.0);
    REQUIRE(residual.blockNorms.combined == Catch::Approx(residual.raw.norm()));
    REQUIRE(residual.intrinsicDensity.size() == static_cast<std::size_t>(n));
    REQUIRE(residual.potentialScale > 0.0);
}

TEST_CASE("NewtonSolver: parses block residual norm controls", "[newton][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"residual_norm", "block"},
        {"max_update", 0.25},
        {"stall_residual_floor", 2.0e-8},
        {"auger_cn_m6_per_s", 4.0e-43},
        {"auger_cp_m6_per_s", 2.0e-43},
        {"residual_weights", {{"psi", 0.25}, {"phin", 2.0}, {"phip", 3.0}}},
        {"residual_scales", {{"psi", 1.0e-18}, {"phin", 2.0e4}, {"phip", 3.0e4}}}
    });

    REQUIRE(cfg.residualNorm == "block");
    REQUIRE(cfg.maxUpdate == Catch::Approx(0.25));
    REQUIRE(cfg.stallResidualFloor == Catch::Approx(2.0e-8));
    REQUIRE(cfg.residualWeightPsi == Catch::Approx(0.25));
    REQUIRE(cfg.residualWeightPhin == Catch::Approx(2.0));
    REQUIRE(cfg.residualWeightPhip == Catch::Approx(3.0));
    REQUIRE(cfg.residualScalePsi == Catch::Approx(1.0e-18));
    REQUIRE(cfg.residualScalePhin == Catch::Approx(2.0e4));
    REQUIRE(cfg.residualScalePhip == Catch::Approx(3.0e4));
    REQUIRE(cfg.augerCn == Catch::Approx(4.0e-43));
    REQUIRE(cfg.augerCp == Catch::Approx(2.0e-43));

    const NewtonConfig boundaryCfg = newtonConfigFromJson(nlohmann::json{
        {"contact_boundary_minority_electron_relaxation", false},
        {"contact_boundary_minority_electron_relaxation_bias_threshold_V", 0.2},
        {"contact_boundary_minority_electron_relaxation_two_terminal_only", false},
        {"contact_boundary_minority_electron_relaxation_contact_side", "both_contacts"},
        {"contact_boundary_minority_electron_relaxation_strength", 0.5},
    });
    REQUIRE_FALSE(boundaryCfg.contactBoundaryMinorityElectronRelaxation);
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V ==
            Catch::Approx(0.2));
    REQUIRE_FALSE(boundaryCfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationContactSide ==
            "both_contacts");
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationStrength ==
            Catch::Approx(0.5));

    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"residual_norm", "unknown"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"max_update", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_update_limit_V", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"carrier_regularization_scale", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "quasi_fermi_update_limit_V",
            std::numeric_limits<Real>::infinity()
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"contact_boundary_minority_electron_relaxation_bias_threshold_V", -1.0}
        }),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"contact_boundary_reconstruction", "unexpected_mode"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "contact_boundary_minority_electron_relaxation_contact_side",
            "unexpected_side"
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "contact_boundary_minority_electron_relaxation_strength",
            1.5
        }}),
        std::invalid_argument);
}

TEST_CASE("NewtonSolver: rejects disabled residual weights", "[newton][config]")
{
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"residual_weights", {{"psi", 0.0}, {"phin", 0.0}, {"phip", 0.0}}}
        }),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"residual_weights", {{"psi", -1.0}, {"phin", 0.0}, {"phip", 1.0}}}
        }),
        std::invalid_argument);

    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"residual_weights", {{"psi", 0.0}, {"phin", 0.0}, {"phip", 1.0}}}
    });
    REQUIRE(cfg.residualWeightPsi == Catch::Approx(0.0));
    REQUIRE(cfg.residualWeightPhin == Catch::Approx(0.0));
    REQUIRE(cfg.residualWeightPhip == Catch::Approx(1.0));
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
    DeviceMesh mesh = makePartiallyContactedOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Zero(N);
    initial.phin = VectorXd::Constant(N, 0.1);
    initial.phip = VectorXd::Constant(N, -0.1);
    initial.phin(0) = 0.0;
    initial.phip(0) = 0.0;

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 3;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.verbose = false;
    cfg.lineSearch = true;
    cfg.warmStart = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"gate", 0.0}}, initial, cfg);

    REQUIRE_FALSE(result.converged);
    REQUIRE(result.iters == 0);
    REQUIRE(result.history.empty());
    REQUIRE(result.finalResidualNorm == Catch::Approx(result.initialResidualNorm));
    REQUIRE(result.failureDiagnostics.failureReason == "carrier_invalid");
    REQUIRE(result.failureDiagnostics.lineSearchFailureReason == "carrier_invalid");
    REQUIRE(result.failureDiagnostics.blockResiduals.psi >= 0.0);
    REQUIRE_FALSE(result.failureDiagnostics.carrierDiagnostics.positiveFinite);
    REQUIRE_FALSE(result.failureDiagnostics.topPoissonResidualNodes.empty());
    REQUIRE((result.solution.psi - initial.psi).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phin - initial.phin).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phip - initial.phip).norm() == Catch::Approx(0.0));
}

TEST_CASE("NewtonSolver: carrier regularization damps coupled Newton carrier mode",
          "[newton][line_search][regularization]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.18;
    initial.phip(interiorNode) = -0.16;

    NewtonConfig baselineCfg = newtonConfig();
    baselineCfg.maxIter = 1;
    baselineCfg.reltol = 0.0;
    baselineCfg.abstol = 0.0;
    baselineCfg.lineSearch = false;
    baselineCfg.warmStart = true;

    NewtonConfig regularizedCfg = baselineCfg;
    regularizedCfg.carrierRegularizationScale = 10.0;

    const NewtonResult baseline = runNewton(
        mesh, matdb, doping, biases, initial, baselineCfg);
    const NewtonResult regularized = runNewton(
        mesh, matdb, doping, biases, initial, regularizedCfg);

    REQUIRE(baseline.iters == 1);
    REQUIRE(regularized.iters == 1);
    REQUIRE_FALSE(baseline.history.empty());
    REQUIRE_FALSE(regularized.history.empty());

    REQUIRE(regularized.history.front().rawStepNorm < baseline.history.front().rawStepNorm);
    REQUIRE((regularized.solution.psi - initial.psi).norm() > 0.0);
    REQUIRE((regularized.solution.phin - initial.phin).norm() <
            (baseline.solution.phin - initial.phin).norm());
}

TEST_CASE("NewtonSolver: max_update limits a large Newton step before line search",
          "[newton][line_search]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 1;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.verbose = false;
    cfg.lineSearch = false;
    cfg.maxUpdate = 0.05;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.iters > 0);
    REQUIRE_FALSE(result.history.empty());
    REQUIRE(result.history.front().rawStepNorm == Catch::Approx(result.history.front().stepNorm));
    REQUIRE(result.history.front().stepNorm <=
            std::sqrt(static_cast<Real>(3 * mesh.numNodes())) * cfg.maxUpdate);
}

TEST_CASE("NewtonSolver: quasi-Fermi update limit recomputes Poisson correction in unit_scaling",
          "[newton][line_search][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    seedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.18;
    initial.phip(interiorNode) = -0.16;

    NewtonConfig unclampedCfg = newtonConfig();
    unclampedCfg.maxIter = 1;
    unclampedCfg.reltol = 0.0;
    unclampedCfg.abstol = 0.0;
    unclampedCfg.dampingFactor = 1.0;
    unclampedCfg.lineSearch = false;
    unclampedCfg.warmStart = true;
    unclampedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};

    NewtonConfig clampedCfg = unclampedCfg;
    clampedCfg.quasiFermiUpdateLimit_V = 0.01;

    const NewtonResult unclamped = runNewton(
        mesh, matdb, doping, biases, initial, unclampedCfg);
    const NewtonResult clamped = runNewton(
        mesh, matdb, doping, biases, initial, clampedCfg);

    REQUIRE(unclamped.iters == 1);
    REQUIRE(clamped.iters == 1);

    const VectorXd unclampedPsiDelta = unclamped.solution.psi - initial.psi;
    const VectorXd clampedPsiDelta = clamped.solution.psi - initial.psi;
    const VectorXd unclampedPhinDelta = unclamped.solution.phin - initial.phin;
    const VectorXd unclampedPhipDelta = unclamped.solution.phip - initial.phip;
    const VectorXd clampedPhinDelta = clamped.solution.phin - initial.phin;
    const VectorXd clampedPhipDelta = clamped.solution.phip - initial.phip;

    const Real limit = clampedCfg.quasiFermiUpdateLimit_V;
    const Real maxUnclampedQfDelta = std::max(
        unclampedPhinDelta.cwiseAbs().maxCoeff(),
        unclampedPhipDelta.cwiseAbs().maxCoeff());

    REQUIRE(maxUnclampedQfDelta > limit);
    REQUIRE((clampedPsiDelta - unclampedPsiDelta).norm() > 1.0e-12);
    REQUIRE(clampedPhinDelta.cwiseAbs().maxCoeff() <= limit + 1.0e-12);
    REQUIRE(clampedPhipDelta.cwiseAbs().maxCoeff() <= limit + 1.0e-12);
    REQUIRE(clampedPhinDelta(interiorNode) ==
            Catch::Approx((unclampedPhinDelta(interiorNode) > 0.0 ? 1.0 : -1.0) * limit)
                .margin(1.0e-12));
    REQUIRE(clampedPhipDelta(interiorNode) ==
            Catch::Approx((unclampedPhipDelta(interiorNode) > 0.0 ? 1.0 : -1.0) * limit)
                .margin(1.0e-12));
}

TEST_CASE("NewtonSolver: max-iteration exit honors stall residual floor", "[newton][line_search]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.stallResidualFloor = 1.0e9;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.converged);
    REQUIRE(result.finalResidualNorm <= cfg.stallResidualFloor);
}


TEST_CASE("NewtonSolver: optionally records line-search diagnostics in history", "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.diagnostics = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.converged);
    REQUIRE_FALSE(result.history.empty());
    const NewtonIterationInfo& first = result.history.front();
    REQUIRE(first.iter == 1);
    REQUIRE(first.rawStepNorm >= first.stepNorm);
    REQUIRE(first.stepNorm == Catch::Approx(first.dampingFactor * first.rawStepNorm));
    REQUIRE(first.relativeResidualNorm == Catch::Approx(
        ResidualNorm::relative(first.residualNorm, result.initialResidualNorm)));
    REQUIRE(first.lineSearchAccepted);
    REQUIRE(first.lineSearchAttempts >= 1);
    REQUIRE(first.lineSearchHistory.size() == static_cast<std::size_t>(first.lineSearchAttempts));
    REQUIRE(first.lineSearchHistory.back().accepted);
    REQUIRE(first.lineSearchHistory.back().damping == Catch::Approx(first.dampingFactor));
    REQUIRE(first.lineSearchHistory.back().residualNorm == Catch::Approx(first.residualNorm));

    const NewtonConfig parsed = newtonConfigFromJson(nlohmann::json{{"diagnostic_history", true}});
    REQUIRE(parsed.diagnostics);
}

TEST_CASE("NewtonSolver: configured temperature is parsed and passed to initial Gummel guess", "[newton][temperature]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg300 = newtonConfig();
    cfg300.maxIter = 0;
    cfg300.temperature_K = 300.0;
    const NewtonResult result300 = runNewton(mesh, matdb, doping, zeroBias(), cfg300);

    NewtonConfig cfg600 = cfg300;
    cfg600.temperature_K = 600.0;
    const NewtonResult result600 = runNewton(mesh, matdb, doping, zeroBias(), cfg600);

    const Real builtIn300 = result300.solution.psi(1) - result300.solution.psi(0);
    const Real builtIn600 = result600.solution.psi(1) - result600.solution.psi(0);
    REQUIRE(builtIn600 > 0.0);
    REQUIRE(builtIn600 < builtIn300);

    const NewtonConfig parsed = newtonConfigFromJson(nlohmann::json{{"temperature_K", 325.0}});
    REQUIRE(parsed.temperature_K == Catch::Approx(325.0));
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{{"temperature_K", -1.0}}),
                      std::invalid_argument);
}

static void requireFiniteNewtonSolution(const NewtonResult& result, Index nodeCount)
{
    REQUIRE(std::isfinite(result.initialResidualNorm));
    REQUIRE(std::isfinite(result.finalResidualNorm));
    for (Index i = 0; i < nodeCount; ++i) {
        const int ii = static_cast<int>(i);
        REQUIRE(std::isfinite(result.solution.psi(ii)));
        REQUIRE(std::isfinite(result.solution.phin(ii)));
        REQUIRE(std::isfinite(result.solution.phip(ii)));
        REQUIRE(std::isfinite(result.solution.n(ii)));
        REQUIRE(std::isfinite(result.solution.p(ii)));
        REQUIRE(result.solution.n(ii) >= 0.0);
        REQUIRE(result.solution.p(ii) >= 0.0);
    }
}

TEST_CASE("NewtonSolver: high doping gradient reverse bias does not diverge", "[newton][stability]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 1.0e21},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 8;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-18;
    cfg.lineSearch = true;
    cfg.verbose = false;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", -0.10}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.iters <= cfg.maxIter);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm * 1.0e6);
}

TEST_CASE("NewtonSolver: multi-terminal contacted mesh accepts distinct biases", "[newton][contacts]")
{
    DeviceMesh mesh;
    const double L = 1.0e-6;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;   n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = L;   mesh.addNode(n3);
    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);
    Region r0; r0.id = 0; r0.name = "p_region"; r0.material = "Si"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);
    Contact body; body.id = 0; body.name = "body"; body.region_id = 0; body.node_ids = {0}; mesh.addContact(body);
    Contact source; source.id = 1; source.name = "source"; source.region_id = 0; source.node_ids = {1}; mesh.addContact(source);
    Contact gate; gate.id = 2; gate.name = "gate"; gate.region_id = 0; gate.node_ids = {2}; mesh.addContact(gate);
    Contact drain; drain.id = 3; drain.name = "drain"; drain.region_id = 0; drain.node_ids = {3}; mesh.addContact(drain);
    mesh.buildEdges();

    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, std::vector<RegionDopingSpec>{{"p_region", 0.0, 1.0e21}});

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;

    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"body", 0.0}, {"source", 0.02}, {"gate", 0.04}, {"drain", 0.06}},
        cfg);

    REQUIRE(result.iters == 0);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    REQUIRE(result.solution.phin(0) == Catch::Approx(0.0));
    REQUIRE(result.solution.phip(0) == Catch::Approx(0.0));
    REQUIRE(result.solution.phin(1) == Catch::Approx(0.02));
    REQUIRE(result.solution.phip(1) == Catch::Approx(0.02));
    REQUIRE(result.solution.phin(2) == Catch::Approx(0.04));
    REQUIRE(result.solution.phip(2) == Catch::Approx(0.04));
    REQUIRE(result.solution.phin(3) == Catch::Approx(0.06));
    REQUIRE(result.solution.phip(3) == Catch::Approx(0.06));
}

TEST_CASE("NewtonSolver: high-bias ohmic contacts keep quasi-Fermi boundary targets when relaxation is disabled",
          "[newton][contacts]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;
    cfg.contactBoundaryMinorityElectronRelaxation = false;

    const Real anodeBias = 0.828125;
    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"anode", anodeBias}, {"cathode", 0.0}},
        cfg);

    REQUIRE(result.iters == 0);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    for (Index nid : mesh.getContact(1).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phin(ii) == Catch::Approx(anodeBias).margin(1.0e-10));
        REQUIRE(result.solution.phip(ii) == Catch::Approx(anodeBias).margin(1.0e-10));
    }
    for (Index nid : mesh.getContact(0).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phin(ii) == Catch::Approx(0.0).margin(1.0e-10));
        REQUIRE(result.solution.phip(ii) == Catch::Approx(0.0).margin(1.0e-10));
    }
}
