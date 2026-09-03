#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
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
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh(bool withContacts = true)
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
    if (withContacts) {
        Contact anode; anode.id=0; anode.name="anode"; anode.region_id=1; anode.node_ids={0,3}; mesh.addContact(anode);
        Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0; cathode.node_ids={1,2}; mesh.addContact(cathode);
    }
    mesh.buildEdges();
    return mesh;
}

static DeviceMesh makeContactInteriorMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;
    Node n0; n0.id=0; n0.x=0.0; n0.y=0.0; mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L; n1.y=0.0; mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=0.0; n2.y=L; mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0.5 * L; n3.y=0.5 * L; mesh.addNode(n3);
    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,3}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=0; c1.node_ids={0,3,2}; mesh.addCell(c1);
    Region r0; r0.id=0; r0.name="si"; r0.material="Si"; r0.cell_ids={0,1}; mesh.addRegion(r0);
    Contact anode; anode.id=0; anode.name="anode"; anode.region_id=0; anode.node_ids={0}; mesh.addContact(anode);
    mesh.buildEdges();
    return mesh;
}

static DeviceMesh makeObtuseAvalancheMesh()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 2.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.2; n2.y = 0.1; mesh.addNode(n2);

    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = {0, 1, 2};
    mesh.addCell(c);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

static DeviceMesh makeEparallelVertexStarMesh()
{
    DeviceMesh mesh;
    const std::vector<Point2> points = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0},
        {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
    };
    for (Index id = 0; id < points.size(); ++id) {
        Node node;
        node.id = id;
        node.x = points[id].x();
        node.y = points[id].y();
        mesh.addNode(node);
    }
    const std::vector<std::vector<Index>> triangles = {
        {0, 1, 4}, {0, 4, 3}, {1, 2, 5}, {1, 5, 4},
    };
    for (Index id = 0; id < triangles.size(); ++id) {
        Cell cell;
        cell.id = id;
        cell.type = CellType::Tri3;
        cell.region_id = 0;
        cell.node_ids = triangles[id];
        mesh.addCell(cell);
    }
    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0, 1, 2, 3};
    mesh.addRegion(region);
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

TEST_CASE("Density-gradient SG avalanche source defaults to symmetric edge partition",
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
    VectorXd phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.05, -0.03);
    VectorXd phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.04, 0.035);
    VectorXd n(static_cast<int>(mesh.numNodes()));
    VectorXd p(static_cast<int>(mesh.numNodes()));
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        n(i) = ni[static_cast<std::size_t>(i)] * std::exp((psi(i) - phin(i)) / Vt);
        p(i) = ni[static_cast<std::size_t>(i)] * std::exp((phip(i) - psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
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

    bool sawPositiveSource = false;
    for (const auto& record : records) {
        if (record.electronSourceIntegral > 0.0) {
            sawPositiveSource = true;
            REQUIRE(record.electronNode0SourceIntegral ==
                    Catch::Approx(0.5 * record.electronSourceIntegral).margin(1.0e-18));
            REQUIRE(record.electronNode1SourceIntegral ==
                    Catch::Approx(0.5 * record.electronSourceIntegral).margin(1.0e-18));
        }
        if (record.holeSourceIntegral > 0.0) {
            sawPositiveSource = true;
            REQUIRE(record.holeNode0SourceIntegral ==
                    Catch::Approx(0.5 * record.holeSourceIntegral).margin(1.0e-18));
            REQUIRE(record.holeNode1SourceIntegral ==
                    Catch::Approx(0.5 * record.holeSourceIntegral).margin(1.0e-18));
        }
        REQUIRE(record.node0SourceIntegral ==
                Catch::Approx(0.5 * record.edgeSourceIntegral).margin(1.0e-18));
        REQUIRE(record.node1SourceIntegral ==
                Catch::Approx(0.5 * record.edgeSourceIntegral).margin(1.0e-18));
    }
    REQUIRE(sawPositiveSource);
}

TEST_CASE("Fermi-Dirac SG avalanche post-processing matches generalized transport flux",
          "[impact][diagnostic][fermi_dirac]")
{
    DeviceMesh mesh = makePNMesh(false);
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e24, 0.0},
        {"p_region", 0.0, 5.0e24},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);
    const Real Vt = constants::Vt_300;
    const CarrierStatisticsConfig statistics{"fermi_dirac"};
    const std::vector<Real> ni(mesh.numNodes(), 1.0e16);
    const std::vector<Real> Nc =
        detail::buildNodeDensityOfStates(mesh, matdb, constants::T0, true);
    const std::vector<Real> Nv =
        detail::buildNodeDensityOfStates(mesh, matdb, constants::T0, false);

    VectorXd psi(static_cast<int>(mesh.numNodes()));
    VectorXd phin(static_cast<int>(mesh.numNodes()));
    VectorXd phip(static_cast<int>(mesh.numNodes()));
    VectorXd n(static_cast<int>(mesh.numNodes()));
    VectorXd p(static_cast<int>(mesh.numNodes()));
    for (int row = 0; row < static_cast<int>(mesh.numNodes()); ++row) {
        psi(row) = 0.58 + 0.035 * static_cast<Real>(row);
        phin(row) = 0.01 * static_cast<Real>(row);
        phip(row) = psi(row) + 0.20 - 0.005 * static_cast<Real>(row);
        n(row) = electronDensity(
            ni[static_cast<std::size_t>(row)], Nc[static_cast<std::size_t>(row)],
            psi(row), phin(row), Vt, statistics);
        p(row) = holeDensity(
            ni[static_cast<std::size_t>(row)], Nv[static_cast<std::size_t>(row)],
            psi(row), phip(row), Vt, statistics);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
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
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    const auto records = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh,
        doping, cellMaterials, psi, phin, phip, n, p, ni, Vt, 1.0, true,
        statistics, Nc, Nv);

    bool sawDegenerateFactor = false;
    for (const auto& record : records) {
        const Edge& edge = mesh.getEdge(record.edgeId);
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real eta0 = (psi(i) - phin(i)) / Vt
            + std::log(ni[edge.n0] / Nc[edge.n0]);
        const Real eta1 = (psi(j) - phin(j)) / Vt
            + std::log(ni[edge.n1] / Nc[edge.n1]);
        const Real driftPotential = psi(j) - psi(i) + Vt * std::log(
            (ni[edge.n1] / Nc[edge.n1]) / (ni[edge.n0] / Nc[edge.n0]));
        const Real factor = sgGeneralizedEinsteinFactor(n(i), n(j), eta0, eta1);
        const Real expected = sgElectronFermiDiracContinuityFlux(
            n(i), n(j), eta0, eta1, driftPotential, phin(i), phin(j), Vt,
            record.electronMobility * Vt / edge.length);

        REQUIRE(record.electronSgUsesFermiDirac);
        REQUIRE(record.holeSgUsesFermiDirac);
        REQUIRE(record.electronSgGeneralizedEinsteinFactor ==
                Catch::Approx(factor).epsilon(1.0e-13));
        REQUIRE(record.electronSgGeneralizedBernoulliArgument == Catch::Approx(
            driftPotential / (Vt * factor)).epsilon(1.0e-13));
        REQUIRE(record.electronSgProductionSignedFluxNative ==
                Catch::Approx(expected).epsilon(1.0e-12));
        sawDegenerateFactor = sawDegenerateFactor || factor > 1.01;
    }
    REQUIRE(sawDegenerateFactor);
}

TEST_CASE("Fermi-Dirac avalanche carrier diagnostics reproduce the coupled residual",
          "[impact][diagnostic][fermi_dirac][newton]")
{
    DeviceMesh mesh = makePNMesh(false);
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e24, 0.0},
        {"p_region", 0.0, 5.0e24},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);
    const Real Vt = constants::Vt_300;

    CoupledDDState state;
    const int nodeCount = static_cast<int>(mesh.numNodes());
    state.psi.resize(nodeCount);
    state.phin.resize(nodeCount);
    state.phip.resize(nodeCount);
    for (int node = 0; node < nodeCount; ++node) {
        state.psi(node) = 0.58 + 0.035 * static_cast<Real>(node);
        state.phin(node) = 0.01 * static_cast<Real>(node);
        state.phip(node) = state.psi(node) + 0.20
            - 0.005 * static_cast<Real>(node);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    const CarrierStatisticsConfig statistics{"fermi_dirac"};
    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig,
        {},
        {},
        {},
        {},
        statistics);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd residual = assembler.residual(x, bcs);
    const auto terms =
        assembler.carrierContinuityEquationTermDiagnostics(x, bcs);

    bool sawImpactSource = false;
    for (int node = 0; node < nodeCount; ++node) {
        const auto& term = terms[static_cast<std::size_t>(node)];
        const Real electronExpected = residual(nodeCount + node);
        const Real holeExpected = residual(2 * nodeCount + node);
        const Real electronScale = std::max<Real>(1.0, std::abs(electronExpected));
        const Real holeScale = std::max<Real>(1.0, std::abs(holeExpected));
        CHECK(std::abs(term.electronResidual - electronExpected) <=
              1.0e-12 * electronScale);
        CHECK(std::abs(term.holeResidual - holeExpected) <=
              1.0e-12 * holeScale);
        sawImpactSource = sawImpactSource ||
            std::abs(term.impactCombinedSource) > 0.0;
    }
    REQUIRE(sawImpactSource);
}

TEST_CASE("Edge-source partition is directional only for grad-QF or explicit switch",
          "[impact][diagnostic]")
{
    ImpactIonizationModelConfig densityGradient;
    densityGradient.generation = "current_density";
    densityGradient.currentApproximation = "density_gradient";
    REQUIRE_FALSE(detail::usesDirectionalEdgeAvalancheSourcePartition(densityGradient));

    ImpactIonizationModelConfig explicitDirectional = densityGradient;
    explicitDirectional.edgeSourcePartition = "qf_gradient";
    REQUIRE(detail::usesDirectionalEdgeAvalancheSourcePartition(explicitDirectional));

    ImpactIonizationModelConfig gradQf = densityGradient;
    gradQf.currentApproximation = "grad_qf";
    REQUIRE(detail::usesDirectionalEdgeAvalancheSourcePartition(gradQf));
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

TEST_CASE("Raw Van Overstraeten diagnostic bypasses minimum field and RefDens damping",
          "[impact][van_overstraeten][diagnostic]")
{
    ImpactIonizationModelConfig cutoffConfig =
        impactIonizationModelConfig("van_overstraeten");
    cutoffConfig.minimumField = 1.0e9;
    const Real field = 2.0e7;
    const auto cutoffModel = makeImpactIonizationModel(cutoffConfig);
    REQUIRE(cutoffModel->electronCoefficient(field) == Catch::Approx(0.0));

    ImpactIonizationModelConfig rawConfig = cutoffConfig;
    rawConfig.debugRawVanOverstraeten = true;
    const auto rawModel = makeImpactIonizationModel(rawConfig);
    const Real expectedRaw = 7.03e7 * std::exp(-1.231e8 / field);
    REQUIRE(rawModel->electronCoefficient(field) ==
            Catch::Approx(expectedRaw).epsilon(1.0e-12));

    rawConfig.drivingForce = "quasi_fermi_gradient";
    rawConfig.drivingForceInterpolation = "quasi_fermi_to_electric_field";
    rawConfig.electronDrivingForceRefDensity = 1.0e30;
    rawConfig.holeDrivingForceRefDensity = 1.0e30;

    REQUIRE(detail::electronAvalancheDrivingField(
                rawConfig,
                field,
                1.0e5,
                1.0e10) == Catch::Approx(field));
    REQUIRE(detail::holeAvalancheDrivingField(
                rawConfig,
                field,
                1.0e5,
                1.0e10) == Catch::Approx(field));
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

TEST_CASE("Coupled DD assembles E2 band-to-band generation into both carrier rows",
          "[btbt][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 0.8;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    RecombinationModelConfig disabled = recombinationModelConfig({"none"});
    RecombinationModelConfig enabled = disabled;
    enabled.bandToBand.model = "e2";
    enabled.bandToBand.prefactorA_SI = 2.0e10;
    enabled.bandToBand.exponentialB_V_per_m = 1.0e6;
    enabled.bandToBand.jacobian = "potential_finite_difference";

    CoupledDDAssembler baseline(mesh, matdb, doping, 0.025852,
                                mobilityConfig, disabled);
    CoupledDDAssembler withBtbt(mesh, matdb, doping, 0.025852,
                                mobilityConfig, enabled);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd x = baseline.pack(state);
    const VectorXd delta = withBtbt.residual(x, bcs) - baseline.residual(x, bcs);
    const int nodeCount = static_cast<int>(mesh.numNodes());

    bool sawGeneration = false;
    for (int node = 0; node < nodeCount; ++node) {
        REQUIRE(delta(nodeCount + node) ==
                Catch::Approx(delta(2 * nodeCount + node)).epsilon(1.0e-13));
        REQUIRE(delta(nodeCount + node) <= 0.0);
        sawGeneration = sawGeneration || delta(nodeCount + node) < 0.0;
    }
    REQUIRE(sawGeneration);

    const Eigen::MatrixXd analyticDelta = Eigen::MatrixXd(
        withBtbt.assembleJacobian(x, bcs) - baseline.assembleJacobian(x, bcs));
    const Eigen::MatrixXd finiteDifferenceDelta = Eigen::MatrixXd(
        withBtbt.finiteDifferenceJacobian(x, bcs, 1.0e-7) -
        baseline.finiteDifferenceJacobian(x, bcs, 1.0e-7));
    const Real denominator = std::max<Real>(1.0, finiteDifferenceDelta.norm());
    REQUIRE((analyticDelta - finiteDifferenceDelta).norm() / denominator < 2.0e-5);
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

TEST_CASE("Postprocess-only avalanche observes source without solver feedback",
          "[impact][postprocess_only]")
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
    state.phin = VectorXd::LinSpaced(
        static_cast<int>(mesh.numNodes()), 0.01, -0.005);
    state.phip = VectorXd::LinSpaced(
        static_cast<int>(mesh.numNodes()), -0.008, 0.006);
    state.psi(1) = 0.1;
    state.psi(2) = 0.1;

    ImpactIonizationModelConfig selfConsistent;
    selfConsistent.model = "selberherr";
    selfConsistent.electronA = 1.0;
    selfConsistent.electronB = 1.0;
    selfConsistent.holeA = 1.0;
    selfConsistent.holeB = 1.0;
    selfConsistent.carrierVelocity = 1.0;

    ImpactIonizationModelConfig postprocess = selfConsistent;
    postprocess.couplingMode = "postprocess_only";

    const auto makeAssembler = [&](const ImpactIonizationModelConfig& impact) {
        return std::make_unique<CoupledDDAssembler>(
            mesh,
            matdb,
            doping,
            Vt,
            mobilityModelConfig("constant"),
            recombinationModelConfig({"none"}),
            BandgapNarrowingConfig{},
            impact);
    };
    const auto disabled = makeAssembler(ImpactIonizationModelConfig{});
    const auto observed = makeAssembler(postprocess);
    const auto coupled = makeAssembler(selfConsistent);

    const VectorXd x = disabled->pack(state);
    const CoupledDDBoundaryConditions bcs;
    REQUIRE(observed->residual(x, bcs).isApprox(disabled->residual(x, bcs), 0.0));
    REQUIRE(
        Eigen::MatrixXd(observed->assembleJacobian(x, bcs))
            .isApprox(Eigen::MatrixXd(disabled->assembleJacobian(x, bcs)), 0.0));

    const auto observedTerms = observed->carrierContinuityTermDiagnostics(x, bcs);
    const auto coupledTerms = coupled->carrierContinuityTermDiagnostics(x, bcs);
    const auto disabledEquationTerms =
        disabled->carrierContinuityEquationTermDiagnostics(x, bcs);
    const auto observedEquationTerms =
        observed->carrierContinuityEquationTermDiagnostics(x, bcs);
    const auto coupledEquationTerms =
        coupled->carrierContinuityEquationTermDiagnostics(x, bcs);
    bool sawSource = false;
    for (std::size_t node = 0; node < observedTerms.size(); ++node) {
        REQUIRE(
            observedTerms[node].electronImpact ==
            Catch::Approx(coupledTerms[node].electronImpact));
        REQUIRE(
            observedTerms[node].holeImpact ==
            Catch::Approx(coupledTerms[node].holeImpact));
        sawSource = sawSource || observedTerms[node].electronImpact != 0.0 ||
            observedTerms[node].holeImpact != 0.0;
        REQUIRE(observedEquationTerms[node].electronResidual ==
                Catch::Approx(disabledEquationTerms[node].electronResidual));
        REQUIRE(observedEquationTerms[node].holeResidual ==
                Catch::Approx(disabledEquationTerms[node].holeResidual));
        REQUIRE(observedEquationTerms[node].electronImpact == 0.0);
        REQUIRE(observedEquationTerms[node].holeImpact == 0.0);
        REQUIRE(coupledEquationTerms[node].electronImpact ==
                Catch::Approx(coupledTerms[node].electronImpact));
        REQUIRE(coupledEquationTerms[node].holeImpact ==
                Catch::Approx(coupledTerms[node].holeImpact));
    }
    REQUIRE(sawSource);
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

TEST_CASE("Coupled DD element-edge GSS Laux avalanche Jacobian matches finite differences",
          "[impact][newton][element_edge_gss_laux]")
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
    state.psi =
        (VectorXd(4) << 0.0, -5.0, -20.0, -2.0).finished();
    state.phin =
        (VectorXd(4) << 0.02, -5.05, -20.10, -2.08).finished();
    state.phip =
        (VectorXd(4) << -0.15, -5.35, -20.82, -2.60).finished();

    ImpactIonizationModelConfig impactConfig =
        impactIonizationModelConfig("van_overstraeten");
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    MobilityModelConfig constantMobilityConfig =
        mobilityModelConfig("caughey_thomas_field");
    constantMobilityConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    const auto mobility = makeMobilityModel(constantMobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const Real temperature_K = Vt * constants::q / constants::kb;
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        constantMobilityConfig,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);
    CoupledDDAssembler baselineAssembler(
        mesh,
        matdb,
        doping,
        Vt,
        constantMobilityConfig,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{});

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const Eigen::MatrixXd analytic =
        Eigen::MatrixXd(assembler.assembleJacobian(x, bcs)) -
        Eigen::MatrixXd(baselineAssembler.assembleJacobian(x, bcs));

    // Independently scatter the shared local-AD derivatives. This fails if the
    // assembler regresses to its former cell-local finite-difference surrogate,
    // even when both external and internal FD steps happen to be 1e-7.
    const int nodeCount = static_cast<int>(mesh.numNodes());
    Eigen::MatrixXd scatteredLocalAd = Eigen::MatrixXd::Zero(x.size(), x.size());
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    const auto& intrinsicDensity = assembler.intrinsicDensity();
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        std::array<detail::Tri3LocalForwardDual, 3> psiAd{};
        std::array<detail::Tri3LocalForwardDual, 3> phinAd{};
        std::array<detail::Tri3LocalForwardDual, 3> phipAd{};
        std::array<detail::Tri3LocalForwardDual, 3> nAd{};
        std::array<detail::Tri3LocalForwardDual, 3> pAd{};
        std::array<Real, 3> niLocal{};
        for (std::size_t local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            psiAd[local] = detail::Tri3LocalForwardDual::variable(
                state.psi(node), local);
            phinAd[local] = detail::Tri3LocalForwardDual::variable(
                state.phin(node), 3 + local);
            phipAd[local] = detail::Tri3LocalForwardDual::variable(
                state.phip(node), 6 + local);
            niLocal[local] = intrinsicDensity[static_cast<std::size_t>(node)];
            nAd[local] = detail::Tri3LocalForwardDual(niLocal[local]) *
                detail::localAdLimitedExp(
                    (psiAd[local] - phinAd[local]) /
                    detail::Tri3LocalForwardDual(Vt));
            pAd[local] = detail::Tri3LocalForwardDual(niLocal[local]) *
                detail::localAdLimitedExp(
                    (phipAd[local] - psiAd[local]) /
                    detail::Tri3LocalForwardDual(Vt));
        }
        const auto sourceAd =
            detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<
                detail::Tri3LocalForwardDual>(
                impactConfig, constantMobilityConfig, *mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh, doping,
                cellMaterials, cellId, psiAd, phinAd, phipAd, nAd, pAd,
                niLocal, Vt, 1.0).combined;
        for (std::size_t localRow = 0; localRow < 3; ++localRow) {
            const int rowNode = static_cast<int>(cell.node_ids[localRow]);
            for (std::size_t localDof = 0;
                 localDof < detail::Tri3LocalPotentialDofCount; ++localDof) {
                const int block = static_cast<int>(localDof / 3);
                const int localColumn = static_cast<int>(localDof % 3);
                const int columnNode = static_cast<int>(
                    cell.node_ids[static_cast<std::size_t>(localColumn)]);
                const int column = block == 0 ? columnNode
                    : (block == 1 ? nodeCount + columnNode
                                  : 2 * nodeCount + columnNode);
                const Real derivative =
                    sourceAd[localRow].derivative[localDof];
                scatteredLocalAd(nodeCount + rowNode, column) -= derivative;
                scatteredLocalAd(2 * nodeCount + rowNode, column) -= derivative;
            }
        }
    }
    const Real scatterReference = std::max<Real>(1.0, analytic.norm());
    REQUIRE((analytic - scatteredLocalAd).norm() / scatterReference <= 1.0e-12);

    const Eigen::MatrixXd directSourceJacobian =
        Eigen::MatrixXd(assembler.impactIonizationSourceJacobian(x, bcs));
    REQUIRE((analytic - directSourceJacobian).norm() /
            std::max<Real>(1.0, analytic.norm()) <= 1.0e-12);
    for (const Real step : {1.0e-14, 3.0e-15, 1.0e-15}) {
        const Eigen::MatrixXd branchResolved =
            Eigen::MatrixXd(
                assembler
                    .impactIonizationSourceBranchResolvedFiniteDifferenceJacobian(
                        x, bcs, step));
        const Real scale =
            std::max(directSourceJacobian.norm(), branchResolved.norm());
        CAPTURE(step, directSourceJacobian.norm(), branchResolved.norm());
        REQUIRE(scale > 0.0);
        REQUIRE((directSourceJacobian - branchResolved).norm() / scale <=
                1.0e-8);
    }

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = Vt;
    scaling.C0 = 2.0;
    scaling.mu0 = 1.0;
    scaling.D0 = 3.0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.fieldFromCoordinateDeltaFactor = 1.0;
    CoupledDDAssembler scaledAssembler(
        mesh, matdb, doping, Vt, constantMobilityConfig,
        recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
        impactConfig, {}, {}, scaling);
    CoupledDDAssembler scaledBaselineAssembler(
        mesh, matdb, doping, Vt, constantMobilityConfig,
        recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{}, {}, {}, scaling);
    CoupledDDState scaledState = state;
    scaledState.psi /= scaling.V0;
    scaledState.phin /= scaling.V0;
    scaledState.phip /= scaling.V0;
    const VectorXd scaledX = scaledAssembler.pack(scaledState);
    const Eigen::MatrixXd scaledAnalytic =
        Eigen::MatrixXd(scaledAssembler.assembleJacobian(scaledX, bcs)) -
        Eigen::MatrixXd(scaledBaselineAssembler.assembleJacobian(scaledX, bcs));
    const Real sourceIntegralFactor =
        scaling.unitSystem.continuitySourceIntegralFactor();
    const Real scaledScatterFactor =
        sourceIntegralFactor * scaling.V0 / (scaling.C0 * scaling.D0);
    const Eigen::MatrixXd scaledScatteredLocalAd =
        scatteredLocalAd * scaledScatterFactor;
    const Real scaledScatterReference =
        std::max<Real>(1.0, scaledAnalytic.norm());
    REQUIRE((scaledAnalytic - scaledScatteredLocalAd).norm() /
            scaledScatterReference <= 1.0e-12);

    const auto sourceResidual = [&](const VectorXd& values) {
        const CoupledDDState replayState = assembler.unpack(values);
        const auto components =
            detail::currentDensityAvalancheSourceComponentIntegrals(
                impactConfig,
                *impact,
                constantMobilityConfig,
                *mobility,
                edgeCells,
                mesh,
                doping,
                cellMaterials,
                replayState.psi,
                replayState.phin,
                replayState.phip,
                assembler.electronDensity(values),
                assembler.holeDensity(values),
                assembler.intrinsicDensity(),
                Vt);
        const int nodeCount = static_cast<int>(mesh.numNodes());
        VectorXd replayResidual = VectorXd::Zero(values.size());
        for (int node = 0; node < nodeCount; ++node) {
            const Real source =
                components.combined[static_cast<std::size_t>(node)];
            replayResidual(nodeCount + node) = -source;
            replayResidual(2 * nodeCount + node) = -source;
        }
        return replayResidual;
    };
    const auto diagnosticTerms =
        assembler.carrierContinuityTermDiagnostics(x, bcs);
    VectorXd diagnosticSourceResidual = VectorXd::Zero(x.size());
    for (int node = 0; node < static_cast<int>(mesh.numNodes()); ++node) {
        const Real source =
            diagnosticTerms[static_cast<std::size_t>(node)].impactCombinedSource;
        diagnosticSourceResidual(static_cast<int>(mesh.numNodes()) + node) = -source;
        diagnosticSourceResidual(2 * static_cast<int>(mesh.numNodes()) + node) = -source;
    }
    const VectorXd replaySourceResidual = sourceResidual(x);
    REQUIRE((diagnosticSourceResidual - replaySourceResidual).norm() /
            std::max<Real>(1.0, diagnosticSourceResidual.norm()) <= 1.0e-12);

    // Keep the strict assembled gate at the validated 1e-7 scale. Independent
    // three-step convergence is covered by the shared local-source AD test.
    constexpr Real relativeStep = 1.0e-7;
    const int M = static_cast<int>(x.size());
    Eigen::MatrixXd finiteDifference = Eigen::MatrixXd::Zero(M, M);
    for (int col = 0; col < M; ++col) {
        const Real step =
            relativeStep * std::max<Real>(1.0, std::abs(x(col)));
        VectorXd plus = x;
        VectorXd minus = x;
        plus(col) += step;
        minus(col) -= step;
        const VectorXd plusSource = sourceResidual(plus);
        const VectorXd minusSource = sourceResidual(minus);
        finiteDifference.col(col) =
            (plusSource - minusSource) / (2.0 * step);
    }

    const int N = static_cast<int>(mesh.numNodes());
    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = N; row < 3 * N; ++row) {
        for (int col = 0; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(analytic(row, col) - finiteDifference(row, col)));
            maxAbsRef =
                std::max(maxAbsRef, std::abs(finiteDifference(row, col)));
        }
    }
    REQUIRE(maxAbsRef > 0.0);
    CAPTURE(maxAbsDiff, maxAbsRef, maxAbsDiff / maxAbsRef);
    REQUIRE(maxAbsDiff / maxAbsRef <= 1.0e-8);

    constexpr Real nearZeroAbsoluteTolerance = 1.0e-12;
    ImpactIonizationModelConfig nearZeroImpactConfig =
        impactIonizationModelConfig("van_overstraeten");
    nearZeroImpactConfig.drivingForce = "quasi_fermi_gradient";
    nearZeroImpactConfig.generation = "current_density";
    nearZeroImpactConfig.currentApproximation = "element_edge_sg_gss_laux";
    nearZeroImpactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    nearZeroImpactConfig.sourceMappingMode = "element_vertex_box_measure";
    CoupledDDAssembler nearZeroAssembler(
        mesh, matdb, doping, Vt, constantMobilityConfig,
        recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
        nearZeroImpactConfig);
    CoupledDDState zeroSourceState = state;
    zeroSourceState.psi.setConstant(0.0);
    zeroSourceState.phin.setConstant(0.0);
    zeroSourceState.phip.setConstant(0.0);
    const VectorXd zeroSourceX = nearZeroAssembler.pack(zeroSourceState);
    const auto zeroSourceTerms =
        nearZeroAssembler.carrierContinuityTermDiagnostics(zeroSourceX, bcs);
    for (const auto& term : zeroSourceTerms)
        REQUIRE(std::abs(term.impactCombinedSource) <= nearZeroAbsoluteTolerance);
    const Eigen::MatrixXd zeroSourceAnalytic =
        Eigen::MatrixXd(nearZeroAssembler.assembleJacobian(zeroSourceX, bcs)) -
        Eigen::MatrixXd(baselineAssembler.assembleJacobian(zeroSourceX, bcs));
    const Eigen::MatrixXd zeroSourceDirect =
        Eigen::MatrixXd(
            nearZeroAssembler.impactIonizationSourceJacobian(
                zeroSourceX, bcs));
    const Eigen::MatrixXd zeroSourceBranchResolved =
        Eigen::MatrixXd(
            nearZeroAssembler
                .impactIonizationSourceBranchResolvedFiniteDifferenceJacobian(
                    zeroSourceX, bcs, 1.0e-15));
    REQUIRE((zeroSourceDirect - zeroSourceBranchResolved)
                .cwiseAbs()
                .maxCoeff() <= nearZeroAbsoluteTolerance);
    Eigen::MatrixXd zeroSourceFiniteDifference = Eigen::MatrixXd::Zero(M, M);
    for (int col = 0; col < M; ++col) {
        const Real step =
            relativeStep * std::max<Real>(1.0, std::abs(zeroSourceX(col)));
        VectorXd plus = zeroSourceX;
        VectorXd minus = zeroSourceX;
        plus(col) += step;
        minus(col) -= step;
        const VectorXd plusSource =
            nearZeroAssembler.residual(plus, bcs) -
            baselineAssembler.residual(plus, bcs);
        const VectorXd minusSource =
            nearZeroAssembler.residual(minus, bcs) -
            baselineAssembler.residual(minus, bcs);
        zeroSourceFiniteDifference.col(col) =
            (plusSource - minusSource) / (2.0 * step);
    }
    REQUIRE((zeroSourceAnalytic - zeroSourceFiniteDifference)
                .cwiseAbs()
                .maxCoeff() <= nearZeroAbsoluteTolerance);
}

TEST_CASE("Scaled source Jacobian perturbation matches scaled coordinates",
          "[impact][newton][element_edge_gss_laux][scaling]")
{
    constexpr Real relativeStep = 1.0e-7;
    constexpr Real potentialScale = 0.025852;
    REQUIRE(detail::physicalPotentialCentralDifferenceStep(
                0.0, potentialScale, relativeStep) ==
            Catch::Approx(relativeStep * potentialScale));
    REQUIRE(detail::physicalPotentialCentralDifferenceStep(
                20.0, potentialScale, relativeStep) ==
            Catch::Approx(relativeStep * 20.0));
}

TEST_CASE("Coupled DD psi-gradient avalanche Jacobian matches carrier finite differences",
          "[impact][newton][psi_gradient_proxy]")
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
    state.psi = VectorXd::LinSpaced(N, -0.01, 0.01);
    state.phin = VectorXd::LinSpaced(N, 0.4, -0.4);
    state.phip = VectorXd::LinSpaced(N, -0.3, 0.3);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "psi_gradient_proxy";
    impactConfig.electronA = 1.0e6;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0e6;
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
    const Eigen::MatrixXd analytic =
        Eigen::MatrixXd(assembler.assembleJacobian(x, bcs));
    const Eigen::MatrixXd finiteDifference =
        Eigen::MatrixXd(assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7));

    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = N; row < 3 * N; ++row) {
        for (int col = 0; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(analytic(row, col) - finiteDifference(row, col)));
            maxAbsRef = std::max(maxAbsRef, std::abs(finiteDifference(row, col)));
        }
    }

    REQUIRE(maxAbsRef > 0.0);
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

        REQUIRE(record.edgeAreaProxy > 0.0);
        REQUIRE(record.edgeSourceIntegral >= 0.0);
        REQUIRE(record.edgeSourceIntegral ==
                Catch::Approx(record.electronSourceIntegral + record.holeSourceIntegral)
                    .margin(1.0e-18));
        REQUIRE(record.node0SourceIntegral ==
                Catch::Approx(record.electronNode0SourceIntegral +
                              record.holeNode0SourceIntegral).margin(1.0e-18));
        REQUIRE(record.node1SourceIntegral ==
                Catch::Approx(record.electronNode1SourceIntegral +
                              record.holeNode1SourceIntegral).margin(1.0e-18));
        fromRecords[static_cast<std::size_t>(record.node0)] += record.node0SourceIntegral;
        fromRecords[static_cast<std::size_t>(record.node1)] += record.node1SourceIntegral;
        electronFromRecords[static_cast<std::size_t>(record.node0)] +=
            record.electronNode0SourceIntegral;
        electronFromRecords[static_cast<std::size_t>(record.node1)] +=
            record.electronNode1SourceIntegral;
        holeFromRecords[static_cast<std::size_t>(record.node0)] +=
            record.holeNode0SourceIntegral;
        holeFromRecords[static_cast<std::size_t>(record.node1)] +=
            record.holeNode1SourceIntegral;
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

TEST_CASE("Avalanche source mapping modes preserve total source while changing node support",
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
    VectorXd phin(static_cast<int>(mesh.numNodes()));
    VectorXd phip(static_cast<int>(mesh.numNodes()));
    phin << 0.020, -0.005, -0.018, 0.004;
    phip << -0.010, 0.006, 0.018, -0.002;
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

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto impact = makeImpactIonizationModel(impactConfig);

    const auto nodeMapped = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt);
    impactConfig.sourceMappingMode = "cell_F_cell_alpha_cell_G_to_node";
    const auto cellMapped = detail::sgEdgeCurrentAvalancheSourceIntegrals(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt);

    REQUIRE(nodeMapped.size() == cellMapped.size());
    Real nodeTotal = 0.0;
    Real cellTotal = 0.0;
    Real l1Difference = 0.0;
    for (std::size_t i = 0; i < nodeMapped.size(); ++i) {
        nodeTotal += nodeMapped[i];
        cellTotal += cellMapped[i];
        l1Difference += std::abs(nodeMapped[i] - cellMapped[i]);
    }
    REQUIRE(nodeTotal > 0.0);
    REQUIRE(cellTotal == Catch::Approx(nodeTotal).epsilon(1.0e-12));
    REQUIRE(l1Difference > 0.0);
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

TEST_CASE("VTK triangle GSS source projection uses triangle records",
          "[impact][diagnostic][vtk][triangle_gss]")
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
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_reconstructed";
    impactConfig.cellReconstructedMidpointDensity = "gss_logistic";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "triangle_gss_gradqf_truncated";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto expectedNodal = detail::currentDensityAvalancheSourceIntegrals(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh, doping,
        cellMaterials, sol.psi, sol.phin, sol.phip, sol.n, sol.p, ni, Vt);
    const auto records = detail::triangleGssAvalancheSourceRecords(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh, doping,
        cellMaterials, sol.psi, sol.phin, sol.phip, sol.n, sol.p, ni, Vt);
    std::vector<Real> expectedElectronIon(mesh.numNodes(), 0.0);
    std::vector<Real> expectedHoleIon(mesh.numNodes(), 0.0);
    for (const auto& record : records) {
        const Real halfLength = 0.5 * record.edgeLength;
        expectedElectronIon[record.node0] += record.electronAlpha * halfLength;
        expectedElectronIon[record.node1] += record.electronAlpha * halfLength;
        expectedHoleIon[record.node0] += record.holeAlpha * halfLength;
        expectedHoleIon[record.node1] += record.holeAlpha * halfLength;
    }

    const auto vtkPath = std::filesystem::temp_directory_path() /
        "vela_triangle_gss_avalanche_generation.vtk";
    writeDDSolutionVTK(
        vtkPath.string(), mesh, matdb, doping, sol, mobilityConfig, recombinationConfig,
        impactConfig, BandgapNarrowingConfig{}, constants::T0);
    const auto avalanche = readVtkScalar(
        vtkPath, "AvalancheGeneration", static_cast<std::size_t>(mesh.numNodes()));
    const auto electronIon = readVtkScalar(
        vtkPath, "ElectronIonIntegral", static_cast<std::size_t>(mesh.numNodes()));
    const auto holeIon = readVtkScalar(
        vtkPath, "HoleIonIntegral", static_cast<std::size_t>(mesh.numNodes()));
    const auto meanIon = readVtkScalar(
        vtkPath, "MeanIonIntegral", static_cast<std::size_t>(mesh.numNodes()));

    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const std::size_t i = static_cast<std::size_t>(node);
        REQUIRE(avalanche[i] * mesh.getNode(node).volume ==
                Catch::Approx(expectedNodal[i]).margin(1.0e-18));
        REQUIRE(electronIon[i] == Catch::Approx(expectedElectronIon[i]).margin(1.0e-18));
        REQUIRE(holeIon[i] == Catch::Approx(expectedHoleIon[i]).margin(1.0e-18));
        REQUIRE(meanIon[i] == Catch::Approx(0.5 * (expectedElectronIon[i] + expectedHoleIon[i]))
                              .margin(1.0e-18));
    }
    std::error_code removeError;
    std::filesystem::remove(vtkPath, removeError);
}
TEST_CASE("VTK exports direct avalanche velocity alpha ion integral and impact drive scalars",
          "[impact][diagnostic][vtk]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -2.0, 2.0);
    sol.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.4, -0.4);
    sol.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.3, 0.3);
    sol.n = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 1.0e21);
    sol.p = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 2.0e21);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "grad_qf";
    impactConfig.electronA = 2.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 3.0;
    impactConfig.holeB = 1.0e-30;

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    mobilityConfig.highFieldDrivingForce = "electric_field";
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});

    const auto vtkPath = std::filesystem::temp_directory_path() /
        "vela_avalanche_direct_scalars.vtk";
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

    const auto electronVelocity =
        readVtkScalar(vtkPath, "ElectronVelocity", static_cast<std::size_t>(mesh.numNodes()));
    const auto holeVelocity =
        readVtkScalar(vtkPath, "HoleVelocity", static_cast<std::size_t>(mesh.numNodes()));
    const auto electronAlpha =
        readVtkScalar(vtkPath, "ElectronAlphaAvalanche", static_cast<std::size_t>(mesh.numNodes()));
    const auto holeAlpha =
        readVtkScalar(vtkPath, "HoleAlphaAvalanche", static_cast<std::size_t>(mesh.numNodes()));
    const auto electronImpactDrive =
        readVtkScalar(vtkPath, "ElectronImpactIonizationDrive", static_cast<std::size_t>(mesh.numNodes()));
    const auto holeImpactDrive =
        readVtkScalar(vtkPath, "HoleImpactIonizationDrive", static_cast<std::size_t>(mesh.numNodes()));
    const auto electronMobilityDrive =
        readVtkScalar(vtkPath, "ElectronHighFieldDrive", static_cast<std::size_t>(mesh.numNodes()));
    const auto electronIon =
        readVtkScalar(vtkPath, "ElectronIonIntegral", static_cast<std::size_t>(mesh.numNodes()));
    const auto holeIon =
        readVtkScalar(vtkPath, "HoleIonIntegral", static_cast<std::size_t>(mesh.numNodes()));
    const auto meanIon =
        readVtkScalar(vtkPath, "MeanIonIntegral", static_cast<std::size_t>(mesh.numNodes()));
    const auto localElectronAlphaLength = readVtkScalar(
        vtkPath,
        "LocalElectronAlphaLengthProxy",
        static_cast<std::size_t>(mesh.numNodes()));
    const auto localHoleAlphaLength = readVtkScalar(
        vtkPath,
        "LocalHoleAlphaLengthProxy",
        static_cast<std::size_t>(mesh.numNodes()));
    const auto localMeanAlphaLength = readVtkScalar(
        vtkPath,
        "LocalMeanAlphaLengthProxy",
        static_cast<std::size_t>(mesh.numNodes()));
    const auto nodeReconstructedElectronMobility = readVtkScalar(
        vtkPath,
        "NodeReconstructedElectronMobility",
        static_cast<std::size_t>(mesh.numNodes()));
    const auto legacyElectronMobility = readVtkScalar(
        vtkPath,
        "ElectronMobility",
        static_cast<std::size_t>(mesh.numNodes()));

    REQUIRE(*std::max_element(electronVelocity.begin(), electronVelocity.end()) > 0.0);
    REQUIRE(*std::max_element(holeVelocity.begin(), holeVelocity.end()) > 0.0);
    REQUIRE(*std::max_element(electronAlpha.begin(), electronAlpha.end()) > 0.0);
    REQUIRE(*std::max_element(holeAlpha.begin(), holeAlpha.end()) > 0.0);
    REQUIRE(*std::max_element(electronImpactDrive.begin(), electronImpactDrive.end()) > 0.0);
    REQUIRE(*std::max_element(holeImpactDrive.begin(), holeImpactDrive.end()) > 0.0);
    REQUIRE(electronImpactDrive.size() == static_cast<std::size_t>(mesh.numNodes()));
    REQUIRE(holeImpactDrive.size() == static_cast<std::size_t>(mesh.numNodes()));
    REQUIRE(electronImpactDrive != electronMobilityDrive);
    REQUIRE(electronIon.size() == static_cast<std::size_t>(mesh.numNodes()));
    REQUIRE(holeIon.size() == static_cast<std::size_t>(mesh.numNodes()));
    REQUIRE(meanIon.size() == static_cast<std::size_t>(mesh.numNodes()));
    for (std::size_t i = 0; i < meanIon.size(); ++i) {
        REQUIRE(meanIon[i] == Catch::Approx(0.5 * (electronIon[i] + holeIon[i])));
        REQUIRE(localElectronAlphaLength[i] ==
                Catch::Approx(electronIon[i]).margin(0.0));
        REQUIRE(localHoleAlphaLength[i] ==
                Catch::Approx(holeIon[i]).margin(0.0));
        REQUIRE(localMeanAlphaLength[i] ==
                Catch::Approx(meanIon[i]).margin(0.0));
        REQUIRE(nodeReconstructedElectronMobility[i] ==
                Catch::Approx(legacyElectronMobility[i]).margin(0.0));
    }

    std::error_code removeError;
    std::filesystem::remove(vtkPath, removeError);
}

TEST_CASE("VTK exports dual-face nodal currents from supplied production SG fluxes",
          "[impact][diagnostic][vtk][dual_face]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    DDSolution sol;
    sol.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    sol.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    sol.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    sol.n = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 1.0e21);
    sol.p = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 2.0e21);

    std::vector<CoupledDDEdgeFluxDiagnostic> edgeFluxes;
    edgeFluxes.reserve(static_cast<std::size_t>(mesh.numEdges()));
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        CoupledDDEdgeFluxDiagnostic diagnostic;
        diagnostic.edgeId = edgeId;
        diagnostic.couple_m = edge.couple;
        diagnostic.electronParticleLineFlux_per_m_s =
            -2.0 * edge.couple * 1.0e4 / constants::q;
        diagnostic.holeParticleLineFlux_per_m_s =
            3.0 * edge.couple * 1.0e4 / constants::q;
        edgeFluxes.push_back(diagnostic);
    }

    const auto vtkPath = std::filesystem::temp_directory_path() /
        "vela_dual_face_sg_current_vectors.vtk";
    writeDDSolutionVTK(
        vtkPath.string(),
        mesh,
        matdb,
        doping,
        sol,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        ImpactIonizationModelConfig{},
        BandgapNarrowingConfig{},
        constants::T0,
        UnitScalingConfig{},
        CarrierStatisticsConfig{},
        &edgeFluxes);

    std::ifstream input(vtkPath);
    REQUIRE(input.good());
    const std::string vtkText(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>());
    REQUIRE(vtkText.find("VECTORS DualFaceSgElectronCurrentDensityVector") !=
            std::string::npos);
    REQUIRE(vtkText.find("VECTORS DualFaceSgHoleCurrentDensityVector") !=
            std::string::npos);
    REQUIRE(vtkText.find("VECTORS DualFaceSgTotalCurrentDensityVector") !=
            std::string::npos);

    std::error_code removeError;
    std::filesystem::remove(vtkPath, removeError);
}

TEST_CASE("Genius-style avalanche source volume truncates obtuse Tri3 edge support",
          "[impact][diagnostic]")
{
    DeviceMesh mesh = makeObtuseAvalancheMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, std::vector<RegionDopingSpec>{{"body", 1.0e21, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    DDSolution sol;
    sol.psi = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.0, 0.2);
    sol.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.01, 0.02);
    sol.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.03, -0.01);
    sol.n = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 1.0e21);
    sol.p = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 2.0e21);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);

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
        0.025852);

    std::unordered_map<std::string, Real> edgeAreaByNodes;
    for (const auto& record : records) {
        edgeAreaByNodes[std::to_string(record.node0) + "-" +
                        std::to_string(record.node1)] = record.edgeAreaProxy;
    }

    REQUIRE(edgeAreaByNodes.size() == 3);
    REQUIRE(edgeAreaByNodes.at("0-1") == Catch::Approx(0.0).margin(1.0e-18));
    REQUIRE(edgeAreaByNodes.at("1-2") == Catch::Approx(0.045138888888888895));
    REQUIRE(edgeAreaByNodes.at("0-2") == Catch::Approx(0.006250000000000002));
}

TEST_CASE("Conservative Genius avalanche support closes transport triangle area",
          "[impact][diagnostic][source-volume]")
{
    const DeviceMesh mesh = makeObtuseAvalancheMesh();
    const MaterialDatabase matdb;
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig config;
    config.sourceVolumePolicy = "genius_conservative";
    Real sourceArea = 0.0;
    for (Index edge = 0; edge < mesh.numEdges(); ++edge) {
        sourceArea += detail::avalancheSourceEdgeArea(
            config, edgeCells, mesh, edge, &cellMaterials);
    }

    REQUIRE(sourceArea == Catch::Approx(0.1).epsilon(1.0e-13));
    REQUIRE_THROWS_AS(
        detail::avalancheSourceEdgeArea(config, edgeCells, mesh, 0),
        std::invalid_argument);

    const NewtonConfig parsed = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "genius_conservative"},
        }}
    });
    REQUIRE(parsed.impactIonization.sourceVolumePolicy == "genius_conservative");
}

TEST_CASE("Nodal Eparallel P1 avalanche mapping conserves assembled source",
          "[impact][diagnostic][nodal-source]")
{
    const DeviceMesh mesh = makePNMesh(false);
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, std::vector<RegionDopingSpec>{
            {"n_region", 1.0e21, 0.0}, {"p_region", 1.0e21, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);
    VectorXd psi(4), phin(4), phip(4);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Real x = mesh.getNode(node).x;
        psi(static_cast<int>(node)) = -1.0e5 * x;
        phin(static_cast<int>(node)) = -5.0e4 * x;
        phip(static_cast<int>(node)) = 0.0;
    }
    const VectorXd n = VectorXd::Constant(4, 1.0e21);
    const VectorXd p = VectorXd::Constant(4, 1.0e18);
    const std::vector<Real> ni(4, 1.0e16);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.couplingMode = "postprocess_only";
    impactConfig.drivingForce = "eparallel";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "nodal_vector_current_reconstructed";
    impactConfig.eparallelFieldRecovery = "nodal_vertex_star";
    impactConfig.sourceMappingMode = "nodal_eparallel_p1";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    const auto impact = makeImpactIonizationModel(impactConfig);
    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);

    const auto records = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh,
        doping, cellMaterials, psi, phin, phip, n, p, ni, 0.025852);
    const auto assembled = detail::currentDensityAvalancheSourceIntegrals(
        impactConfig, *impact, mobilityConfig, *mobility, edgeCells, mesh,
        doping, cellMaterials, psi, phin, phip, n, p, ni, 0.025852);

    Real recordTotal = 0.0;
    std::vector<Real> recordNodeSource(mesh.numNodes(), 0.0);
    for (const auto& record : records) {
        recordTotal += record.edgeSourceIntegral;
        recordNodeSource[record.node0] += record.node0SourceIntegral;
        recordNodeSource[record.node1] += record.node1SourceIntegral;
        REQUIRE(record.electronNode0ImpactField == Catch::Approx(
            detail::sentaurusEparallelAvalancheDrivingField(
                record.node0ElectricFieldVector,
                record.electronNode0CurrentVector)));
        REQUIRE(record.electronNode1ImpactField == Catch::Approx(
            detail::sentaurusEparallelAvalancheDrivingField(
                record.node1ElectricFieldVector,
                record.electronNode1CurrentVector)));
        REQUIRE(record.holeNode0ImpactField == Catch::Approx(
            detail::sentaurusEparallelAvalancheDrivingField(
                record.node0ElectricFieldVector,
                record.holeNode0CurrentVector)));
        REQUIRE(record.holeNode1ImpactField == Catch::Approx(
            detail::sentaurusEparallelAvalancheDrivingField(
                record.node1ElectricFieldVector,
                record.holeNode1CurrentVector)));
        REQUIRE(record.electronNode0Alpha == Catch::Approx(
            impact->electronCoefficient(record.electronNode0ImpactField)));
        REQUIRE(record.electronNode1Alpha == Catch::Approx(
            impact->electronCoefficient(record.electronNode1ImpactField)));
        REQUIRE(record.holeNode0Alpha == Catch::Approx(
            impact->holeCoefficient(record.holeNode0ImpactField)));
        REQUIRE(record.holeNode1Alpha == Catch::Approx(
            impact->holeCoefficient(record.holeNode1ImpactField)));
        REQUIRE(record.node0SourceIntegral == Catch::Approx(
            record.electronNode0SourceIntegral + record.holeNode0SourceIntegral));
        REQUIRE(record.node1SourceIntegral == Catch::Approx(
            record.electronNode1SourceIntegral + record.holeNode1SourceIntegral));
        REQUIRE(record.edgeSourceIntegral == Catch::Approx(
            record.node0SourceIntegral + record.node1SourceIntegral));
    }
    const Real assembledTotal =
        std::accumulate(assembled.begin(), assembled.end(), 0.0);
    REQUIRE(recordTotal > 0.0);
    REQUIRE(assembledTotal == Catch::Approx(recordTotal).epsilon(1.0e-12));
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        REQUIRE(assembled[node] == Catch::Approx(recordNodeSource[node])
            .epsilon(1.0e-12));
    }

    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"coupling_mode", "postprocess_only"},
            {"driving_force", "eparallel"},
            {"generation", "current_density"},
            {"current_approximation", "nodal_vector_current_reconstructed"},
            {"source_mapping_mode", "nodal_eparallel_p1"},
        }}
    }), std::invalid_argument);

    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"coupling_mode", "self_consistent"},
            {"driving_force", "eparallel"},
            {"generation", "current_density"},
            {"current_approximation", "nodal_vector_current_reconstructed"},
            {"source_mapping_mode", "nodal_eparallel_p1"},
        }}
    }), std::invalid_argument);
}

TEST_CASE("JSON solver config selects impact ionization model", "[impact][json]")
{
    REQUIRE(ImpactIonizationModelConfig{}.sourceVolumePolicy == "genius_truncated");
    REQUIRE(ImpactIonizationModelConfig{}.edgeSourcePartition == "symmetric");
    REQUIRE(ImpactIonizationModelConfig{}.couplingMode == "self_consistent");
    REQUIRE(ImpactIonizationModelConfig{}.currentApproximation ==
            "mobility_density_gradient");
    REQUIRE(ImpactIonizationModelConfig{}.cellReconstructedMidpointDensity ==
            "bernoulli");
    REQUIRE(ImpactIonizationModelConfig{}.sourceMappingMode ==
            "node_F_node_alpha_node_G");
    REQUIRE(ImpactIonizationModelConfig{}.eparallelFieldRecovery ==
            "edge_adjacent_cells");

    const GummelConfig cfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"electron_A_m_inv", 1.0e6},
            {"source_geometry_scale", 2.0},
        }}
    });
    REQUIRE(cfg.impactIonization.model == "selberherr");
    REQUIRE(cfg.impactIonization.couplingMode == "self_consistent");
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

    const NewtonConfig postprocessCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"coupling_mode", "postprocess_only"},
        }}
    });
    REQUIRE(postprocessCfg.impactIonization.couplingMode == "postprocess_only");
    const GummelConfig gummelPostprocessCfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"coupling_mode", "postprocess_only"},
        }}
    });
    REQUIRE(gummelPostprocessCfg.impactIonization.couplingMode == "postprocess_only");
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"impact_ionization", {
                {"model", "van_overstraeten"},
                {"coupling_mode", "unknown"},
            }}
        }),
        std::invalid_argument);

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

    const NewtonConfig rawVanOverstraetenCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"debug_raw_vanoverstraeten", true},
        }}
    });
    REQUIRE(rawVanOverstraetenCfg.impactIonization.debugRawVanOverstraeten);

    const NewtonConfig scaledVanOverstraetenCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"A_scale", 2.0},
        }}
    });
    REQUIRE(scaledVanOverstraetenCfg.impactIonization.aScale == Catch::Approx(2.0));

    const GummelConfig gummelScaledVanOverstraetenCfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"A_scale", 4.0},
        }}
    });
    REQUIRE(gummelScaledVanOverstraetenCfg.impactIonization.aScale == Catch::Approx(4.0));
    const NewtonConfig bScaledVanOverstraetenCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"A_scale", 2.0},
            {"B_scale", 0.95},
        }}
    });
    REQUIRE(bScaledVanOverstraetenCfg.impactIonization.aScale == Catch::Approx(2.0));
    REQUIRE(bScaledVanOverstraetenCfg.impactIonization.bScale == Catch::Approx(0.95));

    const GummelConfig gummelBScaledVanOverstraetenCfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"B_scale", 1.10},
        }}
    });
    REQUIRE(gummelBScaledVanOverstraetenCfg.impactIonization.bScale == Catch::Approx(1.10));

    const NewtonConfig sentaurusFitSetCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"parameter_set", "sentaurus_fit_A_B_switch"},
        }}
    });
    REQUIRE(sentaurusFitSetCfg.impactIonization.parameterSet ==
            "sentaurus_fit_A_B_switch");

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
            {"quasi_fermi_gradient_discretization", "cell_gradient"},
        }}
    });
    REQUIRE(gradQfCfg.impactIonization.currentApproximation == "grad_qf");
    REQUIRE(gradQfCfg.impactIonization.quasiFermiGradientDiscretization == "cell_gradient");

    const NewtonConfig dualFaceCurrentMagnitudeCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"current_approximation", "grad_qf"},
            {"current_magnitude_mode", "dual_face_vector_mag"},
        }}
    });
    REQUIRE(dualFaceCurrentMagnitudeCfg.impactIonization.currentMagnitudeMode == "dual_face_vector_mag");

    const NewtonConfig eparallelCellVectorCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "eparallel"},
            {"generation", "current_density"},
            {"current_approximation", "cell_vector_current_reconstructed"},
        }}
    });
    REQUIRE(eparallelCellVectorCfg.impactIonization.drivingForce == "eparallel");
    REQUIRE(eparallelCellVectorCfg.impactIonization.currentApproximation ==
            "cell_vector_current_reconstructed");

    const NewtonConfig eparallelNodalVectorCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "eparallel"},
            {"generation", "current_density"},
            {"current_approximation", "nodal_vector_current_reconstructed"},
            {"eparallel_field_recovery", "nodal_vertex_star"},
        }}
    });
    REQUIRE(eparallelNodalVectorCfg.impactIonization.currentApproximation ==
            "nodal_vector_current_reconstructed");
    REQUIRE(eparallelNodalVectorCfg.impactIonization.eparallelFieldRecovery ==
            "nodal_vertex_star");

    const GummelConfig gummelEparallelNodalVectorCfg = gummelConfigFromJson(
        nlohmann::json{
            {"impact_ionization", {
                {"model", "van_overstraeten"},
                {"driving_force", "eparallel"},
                {"generation", "current_density"},
                {"current_approximation", "nodal_vector_current_reconstructed"},
                {"eparallel_field_recovery", "nodal_vertex_star"},
            }}
        });
    REQUIRE(gummelEparallelNodalVectorCfg.impactIonization.eparallelFieldRecovery ==
            "nodal_vertex_star");

    const NewtonConfig arithmeticMidpointCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"current_approximation", "cell_reconstructed"},
            {"cell_reconstructed_midpoint_density", "arithmetic"},
        }}
    });
    REQUIRE(arithmeticMidpointCfg.impactIonization.cellReconstructedMidpointDensity == "arithmetic");

    const GummelConfig gummelGradQfCfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"quasi_fermi_gradient_discretization", "cell_gradient"},
        }}
    });
    REQUIRE(gummelGradQfCfg.impactIonization.quasiFermiGradientDiscretization ==
            "cell_gradient");

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
            {"edge_source_partition", "qf_gradient"},
            {"quasi_fermi_carrier_truncation", 1.0e-2},
        }}
    });
    REQUIRE(sourceGeometryCfg.impactIonization.sourceGeometryScale == Catch::Approx(4.0));
    REQUIRE(sourceGeometryCfg.impactIonization.sourceVolumePolicy == "genius_truncated");
    REQUIRE(sourceGeometryCfg.impactIonization.edgeSourcePartition == "qf_gradient");

    auto geniusVolumePolicyCfg = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "genius_truncated"},
        }}
    });
    REQUIRE(geniusVolumePolicyCfg.impactIonization.sourceVolumePolicy == "genius_truncated");

    auto volumePolicyCfg = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "edge_box"},
        }}
    });
    REQUIRE(volumePolicyCfg.impactIonization.sourceVolumePolicy == "edge_box");

    const GummelConfig gummelVolumePolicyCfg = gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "edge_box"},
        }}
    });
    REQUIRE(gummelVolumePolicyCfg.impactIonization.sourceVolumePolicy == "edge_box");

    const NewtonConfig cellSourceMappingCfg = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "grad_qf"},
            {"source_mapping_mode", "cell_F_cell_alpha_cell_G_to_node"},
        }}
    });
    REQUIRE(cellSourceMappingCfg.impactIonization.sourceMappingMode ==
            "cell_F_cell_alpha_cell_G_to_node");

    const GummelConfig edgeSourceMappingCfg = gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "grad_qf"},
            {"source_mapping_mode", "edge_F_edge_alpha_edge_G_to_node"},
        }}
    });
    REQUIRE(edgeSourceMappingCfg.impactIonization.sourceMappingMode ==
            "edge_F_edge_alpha_edge_G_to_node");

    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "grad_qf"},
            {"source_mapping_mode", "cell_F_cell_alpha_cell_G_integral_only"},
        }}
    }), std::invalid_argument);

    const NewtonConfig sourceVolumeFactorCfg = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_factor", 0.75},
        }}
    });
    REQUIRE(sourceVolumeFactorCfg.impactIonization.sourceVolumeFactor == Catch::Approx(0.75));

    const GummelConfig gummelSourceVolumeFactorCfg = gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_factor", 0.75},
        }}
    });
    REQUIRE(gummelSourceVolumeFactorCfg.impactIonization.sourceVolumeFactor == Catch::Approx(0.75));

    const NewtonConfig charonLikeCfg = newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "effective_field_parallel_j"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"minimum_field_V_m", 5.0e6},
        }}
    });
    REQUIRE(charonLikeCfg.impactIonization.drivingForce == "effective_field_parallel_j");
    REQUIRE(charonLikeCfg.impactIonization.minimumField == Catch::Approx(5.0e6));

    const GummelConfig gummelCharonLikeCfg = gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "grad_potential_parallel_j"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"minimum_field_V_m", 5.0e6},
        }}
    });
    REQUIRE(gummelCharonLikeCfg.impactIonization.drivingForce == "grad_potential_parallel_j");
    REQUIRE(gummelCharonLikeCfg.impactIonization.minimumField == Catch::Approx(5.0e6));

    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "unsupported"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_policy", "unsupported"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "eparallel"},
            {"generation", "current_density"},
            {"current_approximation", "nodal_vector_current_reconstructed"},
            {"eparallel_field_recovery", "unsupported"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "electric_field"},
            {"eparallel_field_recovery", "nodal_vertex_star"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"edge_source_partition", "unsupported"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"edge_source_partition", "unsupported"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_factor", 1.25},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(gummelConfigFromJson({
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"source_volume_factor", 0.25},
        }}
    }), std::invalid_argument);
    REQUIRE(sourceGeometryCfg.impactIonization.quasiFermiCarrierTruncation == Catch::Approx(1.0e-2));

    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "electrostatic"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"minimum_field_V_m", -1.0},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"A_scale", 0.0},
        }}
    }), std::invalid_argument);    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"B_scale", 0.0},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"B_scale", 0.95},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"driving_force", "quasi_fermi_gradient"},
            {"quasi_fermi_gradient_discretization", "cell_average"},
        }}
    }), std::invalid_argument);
    REQUIRE_THROWS_AS(gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"driving_force", "electric_field"},
            {"quasi_fermi_gradient_discretization", "cell_gradient"},
        }}
    }), std::invalid_argument);
}

TEST_CASE("VanOverstraeten Sentaurus fit parameter sets are explicit overrides",
          "[impact][van_overstraeten]")
{
    const Real fieldBelowDefaultSwitch = 3.0e7; // 3e5 V/cm.

    ImpactIonizationModelConfig defaultConfig = impactIonizationModelConfig("van_overstraeten");
    defaultConfig.temperature_K = 300.0;
    defaultConfig.referenceTemperature_K = 300.0;
    const auto defaultModel = makeImpactIonizationModel(defaultConfig);
    const Real expectedDefaultElectron =
        7.03e7 * std::exp(-1.231e8 / fieldBelowDefaultSwitch);
    REQUIRE(defaultModel->electronCoefficient(fieldBelowDefaultSwitch) ==
            Catch::Approx(expectedDefaultElectron));

    ImpactIonizationModelConfig aOnlyConfig = defaultConfig;
    aOnlyConfig.parameterSet = "sentaurus_fit_A_only";
    const auto aOnlyModel = makeImpactIonizationModel(aOnlyConfig);
    const Real expectedAOnlyElectron =
        2.35990376332e9 * std::exp(-1.231e8 / fieldBelowDefaultSwitch);
    REQUIRE(aOnlyModel->electronCoefficient(fieldBelowDefaultSwitch) ==
            Catch::Approx(expectedAOnlyElectron));
    REQUIRE(aOnlyModel->holeCoefficient(fieldBelowDefaultSwitch) >
            defaultModel->holeCoefficient(fieldBelowDefaultSwitch));

    ImpactIonizationModelConfig abSwitchConfig = defaultConfig;
    abSwitchConfig.parameterSet = "sentaurus_fit_A_B_switch";
    const auto abSwitchModel = makeImpactIonizationModel(abSwitchConfig);
    const Real expectedSwitchElectron =
        6.78391642452e9 * std::exp(-1.21718982697e8 / fieldBelowDefaultSwitch);
    const Real expectedSwitchHole =
        1.41230834668e10 * std::exp(-1.99067614831e8 / fieldBelowDefaultSwitch);
    REQUIRE(abSwitchModel->electronCoefficient(fieldBelowDefaultSwitch) ==
            Catch::Approx(expectedSwitchElectron));
    REQUIRE(abSwitchModel->holeCoefficient(fieldBelowDefaultSwitch) ==
            Catch::Approx(expectedSwitchHole));
}

TEST_CASE("VanOverstraeten B_scale only multiplies B critical fields",
          "[impact][van_overstraeten]")
{
    const Real lowField = 3.0e7;
    const Real highField = 5.0e7;

    ImpactIonizationModelConfig defaultConfig = impactIonizationModelConfig("van_overstraeten");
    defaultConfig.temperature_K = 300.0;
    defaultConfig.referenceTemperature_K = 300.0;
    const auto defaultModel = makeImpactIonizationModel(defaultConfig);

    ImpactIonizationModelConfig scaledConfig = defaultConfig;
    scaledConfig.bScale = 0.90;
    const auto scaledModel = makeImpactIonizationModel(scaledConfig);

    const Real gamma = 1.0;
    const Real expectedElectronLow = gamma * defaultConfig.electronALow *
        std::exp(-defaultConfig.electronBLow * scaledConfig.bScale * gamma / lowField);
    const Real expectedHoleLow = gamma * defaultConfig.holeALow *
        std::exp(-defaultConfig.holeBLow * scaledConfig.bScale * gamma / lowField);
    const Real expectedElectronHigh = gamma * defaultConfig.electronAHigh *
        std::exp(-defaultConfig.electronBHigh * scaledConfig.bScale * gamma / highField);
    const Real expectedHoleHigh = gamma * defaultConfig.holeAHigh *
        std::exp(-defaultConfig.holeBHigh * scaledConfig.bScale * gamma / highField);

    REQUIRE(scaledModel->electronCoefficient(lowField) == Catch::Approx(expectedElectronLow));
    REQUIRE(scaledModel->holeCoefficient(lowField) == Catch::Approx(expectedHoleLow));
    REQUIRE(scaledModel->electronCoefficient(highField) == Catch::Approx(expectedElectronHigh));
    REQUIRE(scaledModel->holeCoefficient(highField) == Catch::Approx(expectedHoleHigh));
    REQUIRE(scaledModel->electronCoefficient(lowField) > defaultModel->electronCoefficient(lowField));
    REQUIRE(scaledModel->holeCoefficient(lowField) > defaultModel->holeCoefficient(lowField));
}
TEST_CASE("VanOverstraeten A_scale only multiplies A prefactors",
          "[impact][van_overstraeten]")
{
    const Real lowField = 3.0e7;
    const Real highField = 5.0e7;

    ImpactIonizationModelConfig defaultConfig = impactIonizationModelConfig("van_overstraeten");
    defaultConfig.temperature_K = 300.0;
    defaultConfig.referenceTemperature_K = 300.0;
    const auto defaultModel = makeImpactIonizationModel(defaultConfig);

    ImpactIonizationModelConfig scaledConfig = defaultConfig;
    scaledConfig.aScale = 2.0;
    const auto scaledModel = makeImpactIonizationModel(scaledConfig);

    REQUIRE(scaledModel->electronCoefficient(lowField) ==
            Catch::Approx(2.0 * defaultModel->electronCoefficient(lowField)));
    REQUIRE(scaledModel->holeCoefficient(lowField) ==
            Catch::Approx(2.0 * defaultModel->holeCoefficient(lowField)));
    REQUIRE(scaledModel->electronCoefficient(highField) ==
            Catch::Approx(2.0 * defaultModel->electronCoefficient(highField)));
    REQUIRE(scaledModel->holeCoefficient(highField) ==
            Catch::Approx(2.0 * defaultModel->holeCoefficient(highField)));
}

TEST_CASE("Van Overstraeten impact ionization supports Charon-style minimum field cutoff",
          "[impact]")
{
    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.minimumField = 5.0e6;
    impactConfig.electronALow = 10.0;
    impactConfig.electronAHigh = 10.0;
    impactConfig.electronBLow = 1.0;
    impactConfig.electronBHigh = 1.0;
    impactConfig.holeALow = 20.0;
    impactConfig.holeAHigh = 20.0;
    impactConfig.holeBLow = 1.0;
    impactConfig.holeBHigh = 1.0;
    impactConfig.switchField = 1.0e7;

    const auto impact = makeImpactIonizationModel(impactConfig);

    REQUIRE(impact->electronCoefficient(4.999e6) == 0.0);
    REQUIRE(impact->holeCoefficient(4.999e6) == 0.0);
    REQUIRE(impact->electronCoefficient(5.001e6) > 0.0);
    REQUIRE(impact->holeCoefficient(5.001e6) > 0.0);
}

TEST_CASE("Current-aligned avalanche driving field keeps only field parallel to carrier current",
          "[impact][diagnostic]")
{
    REQUIRE(detail::parallelCurrentAvalancheDrivingField(2.0e6, 3.0) ==
            Catch::Approx(2.0e6));
    REQUIRE(detail::parallelCurrentAvalancheDrivingField(2.0e6, -3.0) == 0.0);
    REQUIRE(detail::parallelCurrentAvalancheDrivingField(-2.0e6, -3.0) ==
            Catch::Approx(2.0e6));
    REQUIRE(detail::parallelCurrentAvalancheDrivingField(-2.0e6, 3.0) == 0.0);
    REQUIRE(detail::parallelCurrentAvalancheDrivingField(2.0e6, 0.0) == 0.0);
}

TEST_CASE("Sentaurus Eparallel uses independent conventional carrier-current vectors",
          "[impact][eparallel]")
{
    const Point2 electricField{3.0e6, 4.0e6};
    const Point2 electronParticleFlux{-6.0, -8.0};
    const Point2 holeParticleFlux{-8.0, 6.0};

    // Jn = -q*Fn, whereas Jp = q*Fp.
    REQUIRE(detail::sentaurusEparallelAvalancheDrivingField(
                electricField, -electronParticleFlux)
            == Catch::Approx(5.0e6));
    REQUIRE(detail::sentaurusEparallelAvalancheDrivingField(
                electricField, holeParticleFlux)
            == 0.0);
    REQUIRE(detail::sentaurusEparallelAvalancheDrivingField(
                electricField, Point2::Zero())
            == 0.0);
}

TEST_CASE("Nodal SG least-squares reconstruction preserves a constant vector field",
          "[impact][eparallel]")
{
    DeviceMesh mesh = makePNMesh(false);
    const Point2 expected{3.0, -4.0};
    std::vector<Real> signedFlux(
        static_cast<std::size_t>(mesh.numEdges()), 0.0);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        const Node& node0 = mesh.getNode(edge.n0);
        const Node& node1 = mesh.getNode(edge.n1);
        const Point2 tangent{
            (node1.x - node0.x) / edge.length,
            (node1.y - node0.y) / edge.length};
        signedFlux[static_cast<std::size_t>(edgeId)] = expected.dot(tangent);
    }
    const auto nodeEdges = detail::buildNodeEdgeMap(mesh);
    const auto flux = [&](Index edgeId) {
        return signedFlux[static_cast<std::size_t>(edgeId)];
    };
    const auto active = [](Index) { return true; };
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Point2 recovered = detail::nodalLeastSquaresCurrentVector(
            node, nodeEdges, mesh, flux, active);
        REQUIRE(recovered.x() == Catch::Approx(expected.x()).margin(1.0e-12));
        REQUIRE(recovered.y() == Catch::Approx(expected.y()).margin(1.0e-12));
    }
}

TEST_CASE("Dual-face SG nodal reconstruction preserves a constant vector field",
          "[impact][eparallel][dual_face]")
{
    DeviceMesh mesh = makePNMesh(false);
    const Point2 expected{-2.5, 6.0};
    std::vector<Real> signedFlux(
        static_cast<std::size_t>(mesh.numEdges()), 0.0);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        const Node& node0 = mesh.getNode(edge.n0);
        const Node& node1 = mesh.getNode(edge.n1);
        const Point2 tangent{
            (node1.x - node0.x) / edge.length,
            (node1.y - node0.y) / edge.length};
        signedFlux[static_cast<std::size_t>(edgeId)] = expected.dot(tangent);
    }
    const auto nodeEdges = detail::buildNodeEdgeMap(mesh);
    const auto flux = [&](Index edgeId) {
        return signedFlux[static_cast<std::size_t>(edgeId)];
    };
    const auto active = [&](Index edgeId) {
        return mesh.getEdge(edgeId).couple > 0.0;
    };
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Point2 recovered = detail::nodalDualFaceSgCurrentVector(
            node, nodeEdges, mesh, flux, active);
        REQUIRE(recovered.x() == Catch::Approx(expected.x()).margin(1.0e-12));
        REQUIRE(recovered.y() == Catch::Approx(expected.y()).margin(1.0e-12));
    }
}

TEST_CASE("Eparallel nodal vertex-star recovery widens the electric-field stencil",
          "[impact][eparallel]")
{
    DeviceMesh mesh = makeEparallelVertexStarMesh();
    MaterialDatabase matdb;
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto nodeCells = detail::buildNodeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(
        mesh, matdb, constants::T0);
    VectorXd psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    psi(5) = 1.0;
    const auto cellCache = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return psi(static_cast<int>(node)); });
    const auto nodalCache = detail::computeTransportNodalScalarGradientCache(
        mesh, nodeCells, cellMaterials, cellCache);

    const auto target = std::find_if(
        mesh.edges().begin(), mesh.edges().end(), [](const Edge& edge) {
            return (edge.n0 == 1 && edge.n1 == 4) ||
                   (edge.n0 == 4 && edge.n1 == 1);
        });
    REQUIRE(target != mesh.edges().end());

    ImpactIonizationModelConfig config;
    bool valid = false;
    const Point2 legacy = detail::eparallelElectricGradientForEdge(
        config, edgeCells, mesh, target->id, cellCache, nodalCache, valid);
    REQUIRE(valid);
    REQUIRE(legacy.x() == Catch::Approx(0.5));
    REQUIRE(legacy.y() == Catch::Approx(0.0).margin(1.0e-15));

    config.eparallelFieldRecovery = "nodal_vertex_star";
    const Point2 vertexStar = detail::eparallelElectricGradientForEdge(
        config, edgeCells, mesh, target->id, cellCache, nodalCache, valid);
    REQUIRE(valid);
    REQUIRE(vertexStar.x() == Catch::Approx(1.0 / 3.0));
    REQUIRE(vertexStar.y() == Catch::Approx(1.0 / 6.0));
}

TEST_CASE("Grad-QF avalanche source can rebuild driving field with GSS carrier truncation",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makePNMesh(false);
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

TEST_CASE("Grad-QF avalanche source falls back to electric field on contact-adjacent edges",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makeContactInteriorMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"si", 1.0e21, 0.0},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::Zero(nodeCount);
    sol.phin = VectorXd::Zero(nodeCount);
    sol.phip = VectorXd::Zero(nodeCount);
    sol.psi(0) = 0.0;
    sol.psi(3) = 0.1;
    sol.phin(0) = 0.0;
    sol.phin(3) = 2.0;
    sol.phip(0) = 0.0;
    sol.phip(3) = -2.0;
    sol.n = VectorXd::Constant(nodeCount, 1.0e21);
    sol.p = VectorXd::Constant(nodeCount, 1.0e15);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
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

    const auto contactEdgeRecord = std::find_if(
        records.begin(),
        records.end(),
        [](const auto& record) {
            return (record.node0 == 0 && record.node1 == 3) ||
                   (record.node0 == 3 && record.node1 == 0);
        });
    REQUIRE(contactEdgeRecord != records.end());
    const Real electricField = std::abs((sol.psi(3) - sol.psi(0)) /
                                        contactEdgeRecord->edgeLength);
    const Real electronQfField = std::abs((sol.phin(3) - sol.phin(0)) /
                                          contactEdgeRecord->edgeLength);

    REQUIRE(electronQfField > 10.0 * electricField);
    REQUIRE(contactEdgeRecord->electronImpactField == Catch::Approx(electricField));
    REQUIRE(contactEdgeRecord->holeImpactField == Catch::Approx(electricField));

    const auto contactCellInteriorRecord = std::find_if(
        records.begin(),
        records.end(),
        [](const auto& record) {
            return (record.node0 == 1 && record.node1 == 3) ||
                   (record.node0 == 3 && record.node1 == 1);
        });
    REQUIRE(contactCellInteriorRecord != records.end());
    const Real interiorElectricField = std::abs((sol.psi(3) - sol.psi(1)) /
                                                contactCellInteriorRecord->edgeLength);
    const Real interiorElectronQfField = std::abs((sol.phin(3) - sol.phin(1)) /
                                                  contactCellInteriorRecord->edgeLength);

    REQUIRE(interiorElectronQfField > 10.0 * interiorElectricField);
    REQUIRE(contactCellInteriorRecord->electronImpactField == Catch::Approx(interiorElectricField));
    REQUIRE(contactCellInteriorRecord->holeImpactField == Catch::Approx(interiorElectricField));
}

TEST_CASE("High-field driving helper falls back to electric field in contact elements",
          "[impact][mobility]")
{
    DeviceMesh mesh = makeContactInteriorMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto contactNodes = detail::contactNodeMask(mesh);

    const auto contactCellInteriorEdge = std::find_if(
        mesh.edges().begin(),
        mesh.edges().end(),
        [](const Edge& edge) {
            return (edge.n0 == 1 && edge.n1 == 3) ||
                   (edge.n0 == 3 && edge.n1 == 1);
        });
    REQUIRE(contactCellInteriorEdge != mesh.edges().end());

    const Real qfField = 2.0e6;
    const Real electricField = 1.0e5;
    const Real selected = detail::edgeHighFieldDrivingField(
        true,
        qfField,
        electricField,
        edgeCells,
        mesh,
        contactCellInteriorEdge->id,
        contactNodes);

    REQUIRE(selected == Catch::Approx(electricField));
}

TEST_CASE("Grad-QF avalanche source uses quasi-Fermi field with SG current proxy",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makePNMesh(false);
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

    const auto defaultRecords = detail::sgEdgeCurrentAvalancheSourceRecords(
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
    REQUIRE_FALSE(defaultRecords.empty());
    REQUIRE_FALSE(defaultRecords.front().electronSgDiagnosticsCollected);
    REQUIRE(defaultRecords.front().electronSgFluxDecomposition.highPrecisionReferenceFlux == 0.0);


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
        Vt,
        1.0,
        true);

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
    REQUIRE(record.electronSgProductionSignedFluxNative ==
            Catch::Approx(record.electronRawSignedFluxProxy));
    REQUIRE(record.electronSgFluxDecomposition.reconstructedFlux ==
            Catch::Approx(record.electronSgProductionSignedFluxNative));
    REQUIRE(record.electronSgFluxDecomposition.coef ==
            Catch::Approx(record.electronMobility * Vt / record.edgeLength));
    REQUIRE(record.electronSgDiagnosticsCollected);
    REQUIRE(record.electronSgFluxDecomposition.includeNiGradientDrift);
    REQUIRE(std::isfinite(record.electronSgFluxDecomposition.stableFactorizedFlux));
    REQUIRE(std::isfinite(record.electronSgFluxDecomposition.highPrecisionReferenceFlux));
    REQUIRE(std::isfinite(record.electronSgFluxDecomposition.cancellationCondition));
    REQUIRE(std::isfinite(record.electronSgReconstructionRelativeError));
    REQUIRE(std::isfinite(record.electronSgProductionVsReferenceRelativeError));
    REQUIRE(std::isfinite(record.electronSgStableVsReferenceRelativeError));
    REQUIRE(record.electronSgReconstructionRelativeError <= 1.0e-12);
    REQUIRE(record.electronSgProductionVsReferenceRelativeError <= 1.0e-6);
    REQUIRE(record.electronSgStableVsReferenceRelativeError <= 1.0e-6);

    REQUIRE(record.edgeSourceIntegral == Catch::Approx(
        (record.electronAlpha * electronSgFlux + record.holeAlpha * holeSgFlux)
        * record.edgeAreaProxy));
}

TEST_CASE("Genius-style Grad-QF avalanche source uses cell-gradient magnitude",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makePNMesh(false);
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::Zero(nodeCount);
    sol.phin.resize(nodeCount);
    sol.phip.resize(nodeCount);
    sol.n.resize(nodeCount);
    sol.p.resize(nodeCount);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const int i = static_cast<int>(node);
        const Node& point = mesh.getNode(node);
        const Real qf = point.x + 2.0 * point.y;
        sol.phin(i) = qf;
        sol.phip(i) = qf;
        sol.n(i) = ni[static_cast<std::size_t>(node)] *
            std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(node)] *
            std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
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
    const Real expectedCellGradient = std::sqrt(5.0);
    bool sawDifferentEdgeDifference = false;
    for (const auto& record : records) {
        REQUIRE(record.electronImpactField == Catch::Approx(expectedCellGradient));
        REQUIRE(record.holeImpactField == Catch::Approx(expectedCellGradient));
        const Edge& edge = mesh.getEdge(record.edgeId);
        const Real edgeDifference =
            std::abs((sol.phin(static_cast<int>(edge.n1)) -
                      sol.phin(static_cast<int>(edge.n0))) / edge.length);
        if (std::abs(edgeDifference - expectedCellGradient) > 1.0e-9)
            sawDifferentEdgeDifference = true;
    }
    REQUIRE(sawDifferentEdgeDifference);

    constexpr Real fieldFactor = 1.0e4;
    const auto scaledRecords = detail::sgEdgeCurrentAvalancheSourceRecords(
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
        Vt,
        fieldFactor);

    REQUIRE(scaledRecords.size() == records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
        REQUIRE(scaledRecords[i].electronImpactField ==
                Catch::Approx(records[i].electronImpactField * fieldFactor));
        REQUIRE(scaledRecords[i].holeImpactField ==
                Catch::Approx(records[i].holeImpactField * fieldFactor));
        REQUIRE(scaledRecords[i].electronRawFluxProxy ==
                Catch::Approx(records[i].electronRawFluxProxy * fieldFactor));
        REQUIRE(scaledRecords[i].holeRawFluxProxy ==
                Catch::Approx(records[i].holeRawFluxProxy * fieldFactor));
    }
}

TEST_CASE("Genius-style Grad-QF contact edges keep cell-gradient drive",
          "[impact][grad_qf]")
{
    DeviceMesh mesh = makeContactInteriorMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"si", 5.0e22, 0.0},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    const Real Vt = 0.025852;
    DDSolution sol;
    sol.psi = VectorXd::Zero(nodeCount);
    sol.phin.resize(nodeCount);
    sol.phip.resize(nodeCount);
    sol.n.resize(nodeCount);
    sol.p.resize(nodeCount);
    const std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const int i = static_cast<int>(node);
        const Node& point = mesh.getNode(node);
        const Real qf = point.x + 2.0 * point.y;
        sol.phin(i) = qf;
        sol.phip(i) = qf;
        sol.n(i) = ni[static_cast<std::size_t>(node)] *
            std::exp((sol.psi(i) - sol.phin(i)) / Vt);
        sol.p(i) = ni[static_cast<std::size_t>(node)] *
            std::exp((sol.phip(i) - sol.psi(i)) / Vt);
    }

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
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
    const auto impact = makeImpactIonizationModel(impactConfig);
    const auto contactNodes = detail::contactNodeMask(mesh);

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
    const Real expectedCellGradient = std::sqrt(5.0);
    bool checkedContactEdge = false;
    for (const auto& record : records) {
        if (!contactNodes[record.node0] && !contactNodes[record.node1])
            continue;
        checkedContactEdge = true;
        REQUIRE(record.electronImpactField == Catch::Approx(expectedCellGradient));
        REQUIRE(record.holeImpactField == Catch::Approx(expectedCellGradient));
        REQUIRE(record.electricField == Catch::Approx(0.0));
    }
    REQUIRE(checkedContactEdge);
}

TEST_CASE("SG edge current avalanche source supports diagnostic volume policy",
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

    impactConfig.sourceVolumeFactor = 0.75;
    const auto interpolatedImpact = makeImpactIonizationModel(impactConfig);
    const auto interpolated = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *interpolatedImpact,
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

    impactConfig.sourceVolumeFactor = 0.0;
    impactConfig.sourceVolumePolicy = "edge_box";
    const auto fullEdgeImpact = makeImpactIonizationModel(impactConfig);
    const auto fullEdge = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *fullEdgeImpact,
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

    REQUIRE(fullEdge.size() == base.size());
    REQUIRE(interpolated.size() == base.size());
    bool sawNonzeroSource = false;
    for (std::size_t i = 0; i < base.size(); ++i) {
        REQUIRE(interpolated[i].edgeAreaProxy == Catch::Approx(1.5 * base[i].edgeAreaProxy));
        REQUIRE(interpolated[i].edgeSourceIntegral == Catch::Approx(1.5 * base[i].edgeSourceIntegral));
        REQUIRE(interpolated[i].node0SourceIntegral == Catch::Approx(1.5 * base[i].node0SourceIntegral));
        REQUIRE(interpolated[i].node1SourceIntegral == Catch::Approx(1.5 * base[i].node1SourceIntegral));
        REQUIRE(fullEdge[i].edgeAreaProxy == Catch::Approx(2.0 * base[i].edgeAreaProxy));
        REQUIRE(fullEdge[i].edgeSourceIntegral == Catch::Approx(2.0 * base[i].edgeSourceIntegral));
        REQUIRE(fullEdge[i].node0SourceIntegral == Catch::Approx(2.0 * base[i].node0SourceIntegral));
        REQUIRE(fullEdge[i].node1SourceIntegral == Catch::Approx(2.0 * base[i].node1SourceIntegral));
        sawNonzeroSource = sawNonzeroSource || base[i].edgeSourceIntegral > 0.0;
    }
    REQUIRE(sawNonzeroSource);
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

TEST_CASE("SG edge avalanche source records apply coordinate field scaling",
          "[impact][scaling][diagnostic]")
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
    const auto impact = makeImpactIonizationModel(impactConfig);

    const auto unscaled = detail::sgEdgeCurrentAvalancheSourceRecords(
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
    const Real fieldFactor = 1.0e6;
    const auto scaled = detail::sgEdgeCurrentAvalancheSourceRecords(
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
        Vt,
        fieldFactor);

    REQUIRE(scaled.size() == unscaled.size());
    bool checkedEdge = false;
    for (std::size_t i = 0; i < unscaled.size(); ++i) {
        if (unscaled[i].electricField <= 0.0 ||
            unscaled[i].electronRawFluxProxy <= 0.0 ||
            unscaled[i].holeRawFluxProxy <= 0.0) {
            continue;
        }
        REQUIRE(scaled[i].electricField ==
                Catch::Approx(unscaled[i].electricField * fieldFactor));
        REQUIRE(scaled[i].electronImpactField ==
                Catch::Approx(unscaled[i].electronImpactField * fieldFactor));
        REQUIRE(scaled[i].holeImpactField ==
                Catch::Approx(unscaled[i].holeImpactField * fieldFactor));
        REQUIRE(scaled[i].electronRawFluxProxy ==
                Catch::Approx(unscaled[i].electronRawFluxProxy * fieldFactor));
        REQUIRE(scaled[i].holeRawFluxProxy ==
                Catch::Approx(unscaled[i].holeRawFluxProxy * fieldFactor));
        checkedEdge = true;
        break;
    }
    REQUIRE(checkedEdge);
}
TEST_CASE("VanOverstraeten Sentaurus fit parameter sets follow active internal units",
          "[impact][scaling][van_overstraeten]")
{
    ImpactIonizationModelConfig legacyConfig = impactIonizationModelConfig("van_overstraeten");
    legacyConfig.parameterSet = "sentaurus_fit_A_B_switch";
    legacyConfig = applyImpactIonizationParameterSet(legacyConfig);

    ImpactIonizationModelConfig tcadConfig = impactIonizationModelConfig(
        "van_overstraeten", UnitScalingConfig{UnitScalingMode::UnitScaling});
    tcadConfig.parameterSet = "sentaurus_fit_A_B_switch";
    tcadConfig = applyImpactIonizationParameterSet(tcadConfig);

    REQUIRE(legacyConfig.electronALow == Catch::Approx(2.35990376332e9));
    REQUIRE(legacyConfig.electronBLow == Catch::Approx(6.68288073314e7));
    REQUIRE(legacyConfig.switchField == Catch::Approx(2.5e7));

    REQUIRE(tcadConfig.electronALow == Catch::Approx(2.35990376332e7));
    REQUIRE(tcadConfig.electronBLow == Catch::Approx(6.68288073314e5));
    REQUIRE(tcadConfig.switchField == Catch::Approx(2.5e5));
}

TEST_CASE("Triangle GSS production source Jacobian supports mobility doping bases",
          "[impact][newton][triangle_gss][mobility_basis]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);
    doping.setNodeDoping(1, 5.0e22, 5.0e22);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = (VectorXd(4) << 0.0, -5.0, -20.0, -2.0).finished();
    state.phin = (VectorXd(4) << 0.02, -5.05, -20.10, -2.08).finished();
    state.phip = (VectorXd(4) << -0.15, -5.35, -20.82, -2.60).finished();

    ImpactIonizationModelConfig impactConfig =
        impactIonizationModelConfig("van_overstraeten");
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_reconstructed";
    impactConfig.currentMagnitudeMode = "edge_scalar_abs";
    impactConfig.cellReconstructedMidpointDensity = "gss_logistic";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceVolumePolicy = "genius_truncated";
    impactConfig.sourceMappingMode = "triangle_gss_gradqf_truncated";
    impactConfig.quasiFermiCarrierTruncation =
        GENERATE(0.0, 1.0e-2);

    const auto basis = GENERATE(
        std::string("net_doping"),
        std::string("cell_reconstructed_total_impurity"));
    CAPTURE(basis, impactConfig.quasiFermiCarrierTruncation);
    MobilityModelConfig mobilityConfig =
        mobilityModelConfig("caughey_thomas_field");
    mobilityConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    mobilityConfig.dopingConcentrationBasis = basis;
    mobilityConfig.jacobianFieldDerivatives = false;

    CoupledDDAssembler assembler(
        mesh, matdb, doping, Vt, mobilityConfig,
        recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
        impactConfig);
    CoupledDDAssembler baselineAssembler(
        mesh, matdb, doping, Vt, mobilityConfig,
        recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{});

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const Eigen::MatrixXd analyticSource =
        Eigen::MatrixXd(assembler.assembleJacobian(x, bcs)) -
        Eigen::MatrixXd(baselineAssembler.assembleJacobian(x, bcs));
    const Eigen::MatrixXd finiteDifferenceSource =
        Eigen::MatrixXd(assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7)) -
        Eigen::MatrixXd(
            baselineAssembler.finiteDifferenceJacobian(x, bcs, 1.0e-7));

    const int N = static_cast<int>(mesh.numNodes());
    const auto analyticCarrier = analyticSource.block(N, 0, 2 * N, 3 * N);
    const auto finiteDifferenceCarrier =
        finiteDifferenceSource.block(N, 0, 2 * N, 3 * N);
    const Real reference = std::max<Real>(1.0, finiteDifferenceCarrier.norm());
    const Real relativeDifference =
        (analyticCarrier - finiteDifferenceCarrier).norm() / reference;
    CAPTURE(analyticCarrier.norm(), finiteDifferenceCarrier.norm(),
            relativeDifference);
    REQUIRE(analyticCarrier.norm() > 0.0);
    REQUIRE(relativeDifference < 2.0e-4);
}
