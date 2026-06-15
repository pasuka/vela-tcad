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

#include <cmath>
#include <iostream>
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

TEST_CASE("NewtonSolver: defaults to analytic Jacobian", "[newton]")
{
    const NewtonConfig cfg;
    REQUIRE(cfg.jacobian == "analytic");
    REQUIRE_FALSE(cfg.warmStart);
    REQUIRE(cfg.contactBoundaryReconstruction == "dominant_signed_contact_mean");
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxation);
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V ==
            Catch::Approx(0.1));
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
            REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationContactSide == "p_contact_only");

    const NewtonConfig debugCfg = newtonConfigFromJson(nlohmann::json{
        {"jacobian", "finite_difference"},
        {"warm_start", true},
        {"contact_boundary_reconstruction", "legacy_node_local"},
    });
    REQUIRE(debugCfg.jacobian == "finite_difference");
    REQUIRE(debugCfg.warmStart);
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

    REQUIRE_FALSE(cold.converged);
    REQUIRE_FALSE(warm.converged);

    REQUIRE(cold.solution.phin(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(cold.solution.phip(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(warm.solution.phin(interiorNode) ==
            Catch::Approx(initial.phin(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.solution.phip(interiorNode) ==
            Catch::Approx(initial.phip(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.initialResidualNorm != Catch::Approx(cold.initialResidualNorm));
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
        {"auger_cn_m6_per_s", 4.0e-43},
        {"auger_cp_m6_per_s", 2.0e-43},
        {"residual_weights", {{"psi", 0.25}, {"phin", 2.0}, {"phip", 3.0}}},
        {"residual_scales", {{"psi", 1.0e-18}, {"phin", 2.0e4}, {"phip", 3.0e4}}}
    });

    REQUIRE(cfg.residualNorm == "block");
    REQUIRE(cfg.maxUpdate == Catch::Approx(0.25));
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
    REQUIRE(result.failureDiagnostics.failureReason == "carrier_invalid");
    REQUIRE(result.failureDiagnostics.lineSearchFailureReason == "carrier_invalid");
    REQUIRE(result.failureDiagnostics.blockResiduals.psi >= 0.0);
    REQUIRE_FALSE(result.failureDiagnostics.carrierDiagnostics.positiveFinite);
    REQUIRE_FALSE(result.failureDiagnostics.topPoissonResidualNodes.empty());
    REQUIRE((result.solution.psi - initial.psi).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phin - initial.phin).norm() == Catch::Approx(0.0));
    REQUIRE((result.solution.phip - initial.phip).norm() == Catch::Approx(0.0));
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
