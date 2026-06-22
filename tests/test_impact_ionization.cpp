#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/equation/AssemblerUtils.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;
    Node n0; n0.id=0; n0.x=0; n0.y=0; mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L; n1.y=0; mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L; n2.y=L; mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0; n3.y=L; mesh.addNode(n3);
    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,2}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1; c1.node_ids={0,2,3}; mesh.addCell(c1);
    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0}; mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1}; mesh.addRegion(r1);
    Contact anode; anode.id=0; anode.name="anode"; anode.region_id=1; anode.node_ids={0,3}; mesh.addContact(anode);
    Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0; cathode.node_ids={1,2}; mesh.addContact(cathode);
    mesh.buildEdges();
    return mesh;
}


TEST_CASE("Edge avalanche directional weights follow quasi-Fermi gradient direction",
          "[impact][diagnostic]")
{
    DeviceMesh mesh = makePNMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);

    std::vector<Real> phin(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    std::vector<Real> phip(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Real x = mesh.getNode(node).x;
        phin[node] = -x;
        phip[node] = x;
    }

    bool sawHorizontal = false;
    bool sawVertical = false;
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        const Node& node0 = mesh.getNode(edge.n0);
        const Node& node1 = mesh.getNode(edge.n1);
        const auto weights = detail::edgeAvalancheDirectionalWeights(
            edgeCells,
            mesh,
            edgeId,
            [&](Index node) { return phin[node]; },
            [&](Index node) { return phip[node]; });

        if (std::abs(node1.x - node0.x) > 0.0 && std::abs(node1.y - node0.y) <= 1.0e-30) {
            sawHorizontal = true;
            const Real edgeUnitX = (node1.x - node0.x) / edge.length;
            const Real expectedNode0 = 0.5 + 0.5 * edgeUnitX;
            const Real expectedNode1 = 1.0 - expectedNode0;
            REQUIRE(weights.electronNode0 == Catch::Approx(expectedNode0));
            REQUIRE(weights.electronNode1 == Catch::Approx(expectedNode1));
            REQUIRE(weights.holeNode0 == Catch::Approx(expectedNode0));
            REQUIRE(weights.holeNode1 == Catch::Approx(expectedNode1));
        }
        if (std::abs(node1.x - node0.x) <= 1.0e-30 && std::abs(node1.y - node0.y) > 0.0) {
            sawVertical = true;
            REQUIRE(weights.electronNode0 == Catch::Approx(0.5));
            REQUIRE(weights.electronNode1 == Catch::Approx(0.5));
            REQUIRE(weights.holeNode0 == Catch::Approx(0.5));
            REQUIRE(weights.holeNode1 == Catch::Approx(0.5));
        }
    }

    REQUIRE(sawHorizontal);
    REQUIRE(sawVertical);
}

TEST_CASE("Cached edge avalanche directional weights match direct cell-gradient weights",
          "[impact][diagnostic]")
{
    DeviceMesh mesh = makePNMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);

    std::vector<Real> phin(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    std::vector<Real> phip(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Node& point = mesh.getNode(node);
        phin[node] = 0.03 * point.x - 0.02 * point.y;
        phip[node] = -0.01 * point.x + 0.04 * point.y;
    }

    const auto electronGradients = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return phin[node]; });
    const auto holeGradients = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return phip[node]; });

    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const auto direct = detail::edgeAvalancheDirectionalWeights(
            edgeCells,
            mesh,
            edgeId,
            [&](Index node) { return phin[node]; },
            [&](Index node) { return phip[node]; });
        const auto cached = detail::edgeAvalancheDirectionalWeights(
            edgeCells,
            mesh,
            edgeId,
            electronGradients,
            holeGradients);

        REQUIRE(cached.electronNode0 == Catch::Approx(direct.electronNode0));
        REQUIRE(cached.electronNode1 == Catch::Approx(direct.electronNode1));
        REQUIRE(cached.holeNode0 == Catch::Approx(direct.holeNode0));
        REQUIRE(cached.holeNode1 == Catch::Approx(direct.holeNode1));
    }
}

TEST_CASE("Impact ionization none model is zero", "[impact]")
{
    const auto model = makeImpactIonizationModel(impactIonizationModelConfig("none"));
    REQUIRE(model->electronCoefficient(1.0e8) == Catch::Approx(0.0));
    REQUIRE(model->holeCoefficient(1.0e8) == Catch::Approx(0.0));
    REQUIRE(model->generationRate(1.0e8, 1.0e21, 1.0e21) == Catch::Approx(0.0));
}

TEST_CASE("Selberherr impact ionization grows with electric field", "[impact]")
{
    SelberherrImpactIonization model;
    const Real low = model.electronCoefficient(1.0e7);
    const Real high = model.electronCoefficient(5.0e8);
    REQUIRE(low >= 0.0);
    REQUIRE(high > low);
    REQUIRE(model.generationRate(5.0e8, 1.0e20, 2.0e20) > 0.0);
}

TEST_CASE("Van Overstraeten impact ionization matches Sentaurus 2018 silicon defaults",
          "[impact][van_overstraeten]")
{
    const auto model = makeImpactIonizationModel(
        impactIonizationModelConfig("van_overstraeten"));

    const Real lowField = 2.0e7;  // 2e5 V/cm, below E0.
    const Real highField = 5.0e7; // 5e5 V/cm, above E0.

    const Real expectedElectronLow = 7.03e7 * std::exp(-1.231e8 / lowField);
    const Real expectedHoleLow = 1.582e8 * std::exp(-2.036e8 / lowField);
    const Real expectedElectronHigh = 7.03e7 * std::exp(-1.231e8 / highField);
    const Real expectedHoleHigh = 6.71e7 * std::exp(-1.693e8 / highField);

    REQUIRE(model->electronCoefficient(lowField) ==
            Catch::Approx(expectedElectronLow).epsilon(1.0e-12));
    REQUIRE(model->holeCoefficient(lowField) ==
            Catch::Approx(expectedHoleLow).epsilon(1.0e-12));
    REQUIRE(model->electronCoefficient(highField) ==
            Catch::Approx(expectedElectronHigh).epsilon(1.0e-12));
    REQUIRE(model->holeCoefficient(highField) ==
            Catch::Approx(expectedHoleHigh).epsilon(1.0e-12));
}

TEST_CASE("Gummel reverse bias BV regression runs with impact ionization", "[impact][gummel]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);
    const std::unordered_map<std::string, Real> biases = {{"anode", -1.0}, {"cathode", 0.0}};

    GummelConfig cfg;
    cfg.maxIter = 20;
    cfg.reltol = 1.0e-4;
    cfg.abstol = 1.0e12;
    cfg.dampingPsi = 0.3;
    cfg.mobility = mobilityModelConfig("caughey_thomas_field");
    cfg.impactIonization.model = "selberherr";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.psi(i)));
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}


TEST_CASE("Coupled DD residual includes impact-ionization generation", "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});

    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0;
    impactConfig.carrierVelocity = 1.0;
    CoupledDDAssembler withImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd r1 = withImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    bool sawGeneration = false;
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        const Real electronDelta = r1(phinOffset + i) - r0(phinOffset + i);
        const Real holeDelta = r1(phipOffset + i) - r0(phipOffset + i);
        REQUIRE(electronDelta <= 0.0);
        REQUIRE(holeDelta <= 0.0);
        sawGeneration = sawGeneration || electronDelta < 0.0 || holeDelta < 0.0;
    }
    REQUIRE(sawGeneration);
}

TEST_CASE("Quasi-Fermi avalanche driving force ignores built-in electrostatic field",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig qfImpact;
    qfImpact.model = "selberherr";
    qfImpact.drivingForce = "quasi_fermi_gradient";
    qfImpact.generation = "current_density";
    qfImpact.electronA = 1.0;
    qfImpact.electronB = 1.0;
    qfImpact.holeA = 1.0;
    qfImpact.holeB = 1.0;
    CoupledDDAssembler withQuasiFermiImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        qfImpact);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd r1 = withQuasiFermiImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(r1(phinOffset + i) == Catch::Approx(r0(phinOffset + i)).margin(1.0e-18));
        REQUIRE(r1(phipOffset + i) == Catch::Approx(r0(phipOffset + i)).margin(1.0e-18));
    }
}

TEST_CASE("Quasi-Fermi avalanche interpolation falls back to electric field at low density",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig qfImpact;
    qfImpact.model = "selberherr";
    qfImpact.drivingForce = "quasi_fermi_gradient";
    qfImpact.generation = "current_density";
    qfImpact.drivingForceInterpolation = "quasi_fermi_to_electric_field";
    qfImpact.electronDrivingForceRefDensity = 1.0e30;
    qfImpact.holeDrivingForceRefDensity = 1.0e30;
    qfImpact.electronA = 1.0;
    qfImpact.electronB = 1.0;
    qfImpact.holeA = 1.0;
    qfImpact.holeB = 1.0;
    CoupledDDAssembler withInterpolatedImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        qfImpact);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd r1 = withInterpolatedImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    bool sawGeneration = false;
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        const Real electronDelta = r1(phinOffset + i) - r0(phinOffset + i);
        const Real holeDelta = r1(phipOffset + i) - r0(phipOffset + i);
        REQUIRE(electronDelta <= 0.0);
        REQUIRE(holeDelta <= 0.0);
        sawGeneration = sawGeneration || electronDelta < 0.0 || holeDelta < 0.0;
    }
    REQUIRE(sawGeneration);
}

TEST_CASE("SG edge-current avalanche approximation cancels flat quasi-Fermi current",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig localImpact;
    localImpact.model = "selberherr";
    localImpact.drivingForce = "electric_field";
    localImpact.generation = "current_density";
    localImpact.currentApproximation = "mobility_density_gradient";
    localImpact.electronA = 1.0;
    localImpact.electronB = 1.0;
    localImpact.holeA = 1.0;
    localImpact.holeB = 1.0;
    CoupledDDAssembler withLocalCurrentImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        localImpact);

    ImpactIonizationModelConfig sgImpact = localImpact;
    sgImpact.currentApproximation = "density_gradient";
    CoupledDDAssembler withSgCurrentImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        sgImpact);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd rLocal = withLocalCurrentImpact.residual(x, bcs);
    const VectorXd rSg = withSgCurrentImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    bool sawLocalGeneration = false;
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        sawLocalGeneration = sawLocalGeneration ||
            rLocal(phinOffset + i) < r0(phinOffset + i) ||
            rLocal(phipOffset + i) < r0(phipOffset + i);
        REQUIRE(rSg(phinOffset + i) == Catch::Approx(r0(phinOffset + i)).margin(1.0e-18));
        REQUIRE(rSg(phipOffset + i) == Catch::Approx(r0(phipOffset + i)).margin(1.0e-18));
    }
    REQUIRE(sawLocalGeneration);
}

TEST_CASE("Coupled DD analytic avalanche Jacobian matches carrier finite differences", "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.005);
    state.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.008, 0.006);
    state.psi(1) = 0.1;
    state.psi(2) = 0.1;

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0;
    impactConfig.carrierVelocity = 1.0;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd denseAnalytic = Eigen::MatrixXd(analytic);
    const Eigen::MatrixXd denseFiniteDifference = Eigen::MatrixXd(finiteDifference);

    const int N = static_cast<int>(mesh.numNodes());
    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = 0; row < 3 * N; ++row) {
        for (int col = N; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(denseAnalytic(row, col) - denseFiniteDifference(row, col)));
            maxAbsRef = std::max(maxAbsRef, std::abs(denseFiniteDifference(row, col)));
        }
    }

    REQUIRE(maxAbsDiff / std::max<Real>(1.0, maxAbsRef) < 5.0e-5);
}

TEST_CASE("Coupled DD SG edge-current avalanche Jacobian matches carrier finite differences",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.02, 0.025);
    state.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.006);
    state.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.007, 0.005);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "electric_field";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd denseAnalytic = Eigen::MatrixXd(analytic);
    const Eigen::MatrixXd denseFiniteDifference = Eigen::MatrixXd(finiteDifference);

    const int N = static_cast<int>(mesh.numNodes());
    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = N; row < 3 * N; ++row) {
        for (int col = 0; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(denseAnalytic(row, col) - denseFiniteDifference(row, col)));
            maxAbsRef = std::max(maxAbsRef, std::abs(denseFiniteDifference(row, col)));
        }
    }

    REQUIRE(maxAbsDiff / std::max<Real>(1.0, maxAbsRef) < 5.0e-5);
}

TEST_CASE("Coupled DD SG edge-current avalanche Jacobian captures field-dependent alpha",
          "[impact][newton]")
{
    // Strong-avalanche fixture: quasi-Fermi driving force with a field-sensitive
    // ionization coefficient (B comparable to the driving field). The avalanche
    // source therefore depends strongly on the quasi-Fermi gradients through
    // alpha(F), so the analytic Jacobian must carry the dAlpha/dphin and
    // dAlpha/dphip derivatives. A frozen-alpha Jacobian fails this comparison.
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.05, 0.05);
    state.phin = VectorXd::LinSpaced(N, 0.5, -0.5);
    state.phip = VectorXd::LinSpaced(N, -0.5, 0.5);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0e6;
    impactConfig.electronB = 1.0e6;
    impactConfig.holeA = 1.0e6;
    impactConfig.holeB = 1.0e6;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd denseAnalytic = Eigen::MatrixXd(analytic);
    const Eigen::MatrixXd denseFiniteDifference = Eigen::MatrixXd(finiteDifference);

    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = N; row < 3 * N; ++row) {
        for (int col = 0; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(denseAnalytic(row, col) - denseFiniteDifference(row, col)));
            maxAbsRef = std::max(maxAbsRef, std::abs(denseFiniteDifference(row, col)));
        }
    }

    REQUIRE(maxAbsRef > 0.0);
    REQUIRE(maxAbsDiff / std::max<Real>(1.0, maxAbsRef) < 5.0e-5);
}

TEST_CASE("Coupled DD SG avalanche Jacobian captures low-density driving-force interpolation",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.08, 0.08);
    state.phin = VectorXd::LinSpaced(N, 0.3, -0.3);
    state.phip = VectorXd::LinSpaced(N, -0.3, 0.3);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.drivingForceInterpolation = "quasi_fermi_to_electric_field";
    impactConfig.electronDrivingForceRefDensity = 1.0e20;
    impactConfig.holeDrivingForceRefDensity = 1.0e20;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(analytic - finiteDifference);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(finiteDifference);

    REQUIRE(ref.norm() > 0.0);
    REQUIRE(diff.norm() / std::max<Real>(1.0, ref.norm()) < 5.0e-5);
}

TEST_CASE("SG edge-current avalanche records sum to assembled nodal source",
          "[impact][diagnostic]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    VectorXd psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.02, 0.025);
    VectorXd phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.006);
    VectorXd phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.007, 0.005);
    VectorXd n(static_cast<int>(mesh.numNodes()));
    VectorXd p(static_cast<int>(mesh.numNodes()));
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        n(i) = ni[static_cast<std::size_t>(i)] * std::exp((psi(i) - phin(i)) / Vt);
        p(i) = ni[static_cast<std::size_t>(i)] * std::exp((phip(i) - psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "electric_field";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    const auto impact = makeImpactIonizationModel(impactConfig);

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    const auto nodal = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt);
    const auto components = detail::sgEdgeCurrentAvalancheSourceComponentIntegrals(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt);
    const auto records = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt);

    std::vector<Real> fromRecords(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    std::vector<Real> electronFromRecords(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    std::vector<Real> holeFromRecords(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    Real totalEdgeSource = 0.0;
    Real totalElectronEdgeSource = 0.0;
    Real totalHoleEdgeSource = 0.0;
    for (const auto& record : records) {
        const auto weights = detail::edgeAvalancheDirectionalWeights(
            edgeCells,
            mesh,
            record.edgeId,
            [&](Index node) { return phin(static_cast<int>(node)); },
            [&](Index node) { return phip(static_cast<int>(node)); });
        const Real electronNode0Weight = weights.electronNode0;
        const Real electronNode1Weight = weights.electronNode1;
        const Real holeNode0Weight = weights.holeNode0;
        const Real holeNode1Weight = weights.holeNode1;
        const Real expectedNode0 =
            electronNode0Weight * record.electronSourceIntegral +
            holeNode0Weight * record.holeSourceIntegral;
        const Real expectedNode1 =
            electronNode1Weight * record.electronSourceIntegral +
            holeNode1Weight * record.holeSourceIntegral;

        REQUIRE(record.edgeAreaProxy > 0.0);
        REQUIRE(record.edgeSourceIntegral >= 0.0);
        REQUIRE(record.edgeSourceIntegral ==
                Catch::Approx(record.electronSourceIntegral + record.holeSourceIntegral)
                    .margin(1.0e-18));
        REQUIRE(record.node0SourceIntegral == Catch::Approx(expectedNode0));
        REQUIRE(record.node1SourceIntegral == Catch::Approx(expectedNode1));
        fromRecords[static_cast<std::size_t>(record.node0)] += record.node0SourceIntegral;
        fromRecords[static_cast<std::size_t>(record.node1)] += record.node1SourceIntegral;
        electronFromRecords[static_cast<std::size_t>(record.node0)] +=
            electronNode0Weight * record.electronSourceIntegral;
        electronFromRecords[static_cast<std::size_t>(record.node1)] +=
            electronNode1Weight * record.electronSourceIntegral;
        holeFromRecords[static_cast<std::size_t>(record.node0)] +=
            holeNode0Weight * record.holeSourceIntegral;
        holeFromRecords[static_cast<std::size_t>(record.node1)] +=
            holeNode1Weight * record.holeSourceIntegral;
        totalEdgeSource += record.edgeSourceIntegral;
        totalElectronEdgeSource += record.electronSourceIntegral;
        totalHoleEdgeSource += record.holeSourceIntegral;
    }

    REQUIRE(totalEdgeSource > 0.0);
    REQUIRE(totalElectronEdgeSource > 0.0);
    REQUIRE(totalHoleEdgeSource > 0.0);
    Real totalNodalSource = 0.0;
    Real totalElectronNodalSource = 0.0;
    Real totalHoleNodalSource = 0.0;
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        totalNodalSource += nodal[static_cast<std::size_t>(node)];
        totalElectronNodalSource += components.electron[static_cast<std::size_t>(node)];
        totalHoleNodalSource += components.hole[static_cast<std::size_t>(node)];
        REQUIRE(fromRecords[static_cast<std::size_t>(node)] ==
                Catch::Approx(nodal[static_cast<std::size_t>(node)]).margin(1.0e-18));
        REQUIRE(electronFromRecords[static_cast<std::size_t>(node)] ==
                Catch::Approx(components.electron[static_cast<std::size_t>(node)])
                    .margin(1.0e-18));
        REQUIRE(holeFromRecords[static_cast<std::size_t>(node)] ==
                Catch::Approx(components.hole[static_cast<std::size_t>(node)])
                    .margin(1.0e-18));
        REQUIRE(nodal[static_cast<std::size_t>(node)] ==
                Catch::Approx(components.combined[static_cast<std::size_t>(node)])
                    .margin(1.0e-18));
        REQUIRE(components.combined[static_cast<std::size_t>(node)] ==
                Catch::Approx(components.electron[static_cast<std::size_t>(node)] +
                              components.hole[static_cast<std::size_t>(node)])
                    .margin(1.0e-18));
    }
    REQUIRE(totalEdgeSource == Catch::Approx(totalNodalSource).margin(1.0e-18));
    REQUIRE(totalElectronEdgeSource ==
            Catch::Approx(totalElectronNodalSource).margin(1.0e-18));
    REQUIRE(totalHoleEdgeSource == Catch::Approx(totalHoleNodalSource).margin(1.0e-18));
}

static std::vector<Real> readVtkScalar(const std::filesystem::path& path,
                                       const std::string& name,
                                       std::size_t count)
{
    std::ifstream input(path);
    REQUIRE(input.good());
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream header(line);
        std::string token;
        std::string scalarName;
        header >> token >> scalarName;
        if (token != "SCALARS" || scalarName != name)
            continue;
        REQUIRE(std::getline(input, line));
        std::vector<Real> values;
        while (values.size() < count && std::getline(input, line)) {
            std::istringstream row(line);
            Real value = 0.0;
            while (row >> value)
                values.push_back(value);
        }
        REQUIRE(values.size() >= count);
        values.resize(count);
        return values;
    }
    FAIL("missing VTK scalar " << name);
    return {};
}

TEST_CASE("VTK AvalancheGeneration uses SG edge nodal source over node volume",
          "[impact][diagnostic][vtk]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.02, 0.025);
    sol.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.006);
    sol.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.007, 0.005);
    sol.n.resize(static_cast<int>(mesh.numNodes()));
    sol.p.resize(static_cast<int>(mesh.numNodes()));
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        sol.n(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "electric_field";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    const auto impact = makeImpactIonizationModel(impactConfig);

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto expectedNodal = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);

    const auto vtkPath = std::filesystem::temp_directory_path() /
        "vela_sg_avalanche_generation_volume_policy.vtk";
    writeDDSolutionVTK(
        vtkPath.string(),
        mesh,
        matdb,
        doping,
        sol,
        mobilityConfig,
        recombinationConfig,
        impactConfig,
        BandgapNarrowingConfig{},
        constants::T0);
    const std::vector<Real> avalanche =
        readVtkScalar(vtkPath, "AvalancheGeneration", static_cast<std::size_t>(mesh.numNodes()));

    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Real integral = avalanche[static_cast<std::size_t>(node)] * mesh.getNode(node).volume;
        REQUIRE(integral == Catch::Approx(expectedNodal[static_cast<std::size_t>(node)]).margin(1.0e-18));
    }
    std::error_code removeError;
    std::filesystem::remove(vtkPath, removeError);
}

TEST_CASE("JSON solver config selects impact ionization model", "[impact][json]")
{
    const GummelConfig cfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"electron_A_m_inv", 1.0e6},
            {"source_geometry_scale", 2.0},
        }}
    });
    REQUIRE(cfg.impactIonization.model == "selberherr");
    REQUIRE(cfg.impactIonization.electronA == Catch::Approx(1.0e6));
    REQUIRE(cfg.impactIonization.sourceGeometryScale == Catch::Approx(2.0));

    const NewtonConfig stringCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", "selberherr"}
    });
    REQUIRE(stringCfg.impactIonization.model == "selberherr");

    const NewtonConfig objectCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"electron_A_m_inv", 2.0e6},
            {"electron_B_V_m", 3.0e7},
            {"hole_A_m_inv", 4.0e6},
            {"hole_B_V_m", 5.0e7},
            {"carrier_velocity_m_s", 6.0e4},
        }}
    });
    REQUIRE(objectCfg.impactIonization.model == "selberherr");
    REQUIRE(objectCfg.impactIonization.electronA == Catch::Approx(2.0e6));
    REQUIRE(objectCfg.impactIonization.electronB == Catch::Approx(3.0e7));
    REQUIRE(objectCfg.impactIonization.holeA == Catch::Approx(4.0e6));
    REQUIRE(objectCfg.impactIonization.holeB == Catch::Approx(5.0e7));
    REQUIRE(objectCfg.impactIonization.carrierVelocity == Catch::Approx(6.0e4));

    const NewtonConfig vanOverstraetenCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"electron_a_low_m_inv", 1.0e6},
            {"electron_b_low_V_m", 2.0e7},
            {"hole_a_high_m_inv", 3.0e6},
            {"hole_b_high_V_m", 4.0e7},
            {"switch_field_V_m", 5.0e7},
            {"phonon_energy_eV", 0.063},
            {"temperature_K", 300.0},
        }}
    });
    REQUIRE(vanOverstraetenCfg.impactIonization.model == "van_overstraeten");
    REQUIRE(vanOverstraetenCfg.impactIonization.electronALow == Catch::Approx(1.0e6));
    REQUIRE(vanOverstraetenCfg.impactIonization.electronBLow == Catch::Approx(2.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.holeAHigh == Catch::Approx(3.0e6));
    REQUIRE(vanOverstraetenCfg.impactIonization.holeBHigh == Catch::Approx(4.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.switchField == Catch::Approx(5.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.phononEnergy == Catch::Approx(0.063));
    REQUIRE(vanOverstraetenCfg.impactIonization.temperature_K == Catch::Approx(300.0));

    const NewtonConfig sentaurusCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
        }}
    });
    REQUIRE(sentaurusCfg.impactIonization.model == "van_overstraeten");
    REQUIRE(sentaurusCfg.impactIonization.drivingForce == "quasi_fermi_gradient");
    REQUIRE(sentaurusCfg.impactIonization.generation == "current_density");
    REQUIRE(sentaurusCfg.impactIonization.currentApproximation == "density_gradient");

    const NewtonConfig gradQfCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"current_approximation", "grad_qf"},
        }}
    });
    REQUIRE(gradQfCfg.impactIonization.currentApproximation == "grad_qf");

    const NewtonConfig interpolatedCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"driving_force_interpolation", {
                {"mode", "quasi_fermi_to_electric_field"},
                {"electron_ref_density_m3", 1.0e16},
                {"hole_ref_density_m3", 2.0e16},
            }},
        }}
    });
    REQUIRE(interpolatedCfg.impactIonization.drivingForceInterpolation ==
            "quasi_fermi_to_electric_field");
    REQUIRE(interpolatedCfg.impactIonization.electronDrivingForceRefDensity ==
            Catch::Approx(1.0e16));
    REQUIRE(interpolatedCfg.impactIonization.holeDrivingForceRefDensity ==
            Catch::Approx(2.0e16));

    const NewtonConfig sourceGeometryCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_geometry_scale", 4.0},
            {"quasi_fermi_carrier_truncation", 1.0e-2},
        }}
    });
    REQUIRE(sourceGeometryCfg.impactIonization.sourceGeometryScale == Catch::Approx(4.0));
    REQUIRE(sourceGeometryCfg.impactIonization.quasiFermiCarrierTruncation == Catch::Approx(1.0e-2));

    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "electrostatic"},
        }}
    }), std::invalid_argument);
}

TEST_CASE("Grad-QF avalanche source can rebuild driving field with GSS carrier truncation",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(nodeCount, -0.02, 0.03);
    sol.phin = VectorXd::LinSpaced(nodeCount, 0.9, -0.7);
    sol.phip = VectorXd::LinSpaced(nodeCount, -0.6, 0.8);
    sol.n.resize(nodeCount);
    sol.p.resize(nodeCount);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < nodeCount; ++i) {
        sol.n(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "grad_qf";
    impactConfig.quasiFermiCarrierTruncation = 1.0e-2;
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto impact = makeImpactIonizationModel(impactConfig);

    const auto records = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);

    REQUIRE_FALSE(records.empty());
    const auto& record = records.front();
    const int i = static_cast<int>(record.node0);
    const int j = static_cast<int>(record.node1);
    const auto truncatedElectronQf = [&](int node) {
        const Real carrier = std::max(sol.n(node), impactConfig.quasiFermiCarrierTruncation * ni[static_cast<std::size_t>(node)]);
        return sol.psi(node) - Vt * std::log(carrier / ni[static_cast<std::size_t>(node)]);
    };
    const auto truncatedHoleQf = [&](int node) {
        const Real carrier = std::max(sol.p(node), impactConfig.quasiFermiCarrierTruncation * ni[static_cast<std::size_t>(node)]);
        return sol.psi(node) + Vt * std::log(carrier / ni[static_cast<std::size_t>(node)]);
    };

    REQUIRE(record.electronImpactField ==
            Catch::Approx(std::abs(truncatedElectronQf(j) - truncatedElectronQf(i)) /
                          record.edgeLength));
    REQUIRE(record.holeImpactField ==
            Catch::Approx(std::abs(truncatedHoleQf(j) - truncatedHoleQf(i)) /
                          record.edgeLength));
}

TEST_CASE("Grad-QF avalanche source uses quasi-Fermi field with SG current proxy",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(nodeCount, -0.02, 0.03);
    sol.phin = VectorXd::LinSpaced(nodeCount, 0.015, -0.009);
    sol.phip = VectorXd::LinSpaced(nodeCount, -0.011, 0.007);
    sol.n.resize(nodeCount);
    sol.p.resize(nodeCount);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < nodeCount; ++i) {
        sol.n(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "grad_qf";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto impact = makeImpactIonizationModel(impactConfig);

    const auto records = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);

    REQUIRE_FALSE(records.empty());
    const auto& record = records.front();
    const int i = static_cast<int>(record.node0);
    const int j = static_cast<int>(record.node1);
    const Real electronQfField = std::abs(sol.phin(j) - sol.phin(i)) / record.edgeLength;
    const Real holeQfField = std::abs(sol.phip(j) - sol.phip(i)) / record.edgeLength;
    const Real electronSgFlux = std::abs(sgElectronContinuityFluxFromQuasiFermiVariableNi(
        ni[record.node0],
        ni[record.node1],
        sol.psi(i),
        sol.psi(j),
        sol.phin(i),
        sol.phin(j),
        Vt,
        record.electronMobility * Vt / record.edgeLength));
    const Real holeSgFlux = std::abs(sgHoleContinuityFluxFromQuasiFermiVariableNi(
        ni[record.node0],
        ni[record.node1],
        sol.psi(i),
        sol.psi(j),
        sol.phip(i),
        sol.phip(j),
        Vt,
        record.holeMobility * Vt / record.edgeLength));

    REQUIRE(record.electronImpactField == Catch::Approx(electronQfField));
    REQUIRE(record.holeImpactField == Catch::Approx(holeQfField));
    REQUIRE(record.electronFluxProxy == Catch::Approx(electronSgFlux));
    REQUIRE(record.holeFluxProxy == Catch::Approx(holeSgFlux));
    REQUIRE(record.edgeSourceIntegral == Catch::Approx(
        (record.electronAlpha * electronSgFlux + record.holeAlpha * holeSgFlux)
        * record.edgeAreaProxy));
}
TEST_CASE("SG edge current avalanche source supports diagnostic geometry scale",
          "[impact][diagnostic]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.02, 0.025);
    sol.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.006);
    sol.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.007, 0.005);
    sol.n.resize(static_cast<int>(mesh.numNodes()));
    sol.p.resize(static_cast<int>(mesh.numNodes()));
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        sol.n(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(i)] * std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "electric_field";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    const auto baseImpact = makeImpactIonizationModel(impactConfig);
    const auto base = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *baseImpact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);
    const auto baseNodal = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig,
        *baseImpact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);

    impactConfig.sourceGeometryScale = 4.0;
    const auto scaledImpact = makeImpactIonizationModel(impactConfig);
    const auto scaled = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *scaledImpact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);
    const auto scaledNodal = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig,
        *scaledImpact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        sol.psi,
        sol.phin,
        sol.phip,
        sol.n,
        sol.p,
        ni,
        Vt);

    REQUIRE(scaled.size() == base.size());
    bool sawNonzeroSource = false;
    for (std::size_t i = 0; i < base.size(); ++i) {
        REQUIRE(scaled[i].edgeAreaProxy == Catch::Approx(4.0 * base[i].edgeAreaProxy));
        REQUIRE(scaled[i].edgeSourceIntegral == Catch::Approx(4.0 * base[i].edgeSourceIntegral));
        REQUIRE(scaled[i].node0SourceIntegral == Catch::Approx(4.0 * base[i].node0SourceIntegral));
        REQUIRE(scaled[i].node1SourceIntegral == Catch::Approx(4.0 * base[i].node1SourceIntegral));
        sawNonzeroSource = sawNonzeroSource || base[i].edgeSourceIntegral > 0.0;
    }
    REQUIRE(scaledNodal.size() == baseNodal.size());
    for (std::size_t i = 0; i < baseNodal.size(); ++i)
        REQUIRE(scaledNodal[i] == Catch::Approx(4.0 * baseNodal[i]));
    REQUIRE(sawNonzeroSource);
}
