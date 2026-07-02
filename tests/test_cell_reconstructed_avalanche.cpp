#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace vela;

static DeviceMesh makeSingleCellMesh()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0; n2.y = 1.0e-6; mesh.addNode(n2);

    Cell cell; cell.id = 0; cell.type = CellType::Tri3; cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);

    Region region; region.id = 0; region.name = "si"; region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);

    mesh.buildEdges();
    return mesh;
}

TEST_CASE("Cell reconstructed avalanche support uses local current density magnitude",
          "[impact][cell_reconstructed]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_reconstructed";
    impactConfig.sourceVolumePolicy = "edge_box";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impactConfig, "test"));
    REQUIRE(detail::usesEdgeCurrentAvalancheSource(impactConfig));

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);

    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    psi << 0.0, -0.20, 0.0;
    phin << 0.0, -0.40, 0.0;
    phip << 0.0, 0.30, 0.0;
    n << 1.0e20, 3.0e20, 1.0e20;
    p << 2.0e20, 6.0e20, 2.0e20;

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
        constants::kb * constants::T0 / constants::q);

    const auto it = std::find_if(records.begin(), records.end(), [](const auto& record) {
        return record.node0 == 0 && record.node1 == 1;
    });
    REQUIRE(it != records.end());

    const Real Vt = constants::kb * constants::T0 / constants::q;
    auto aux2 = [](Real x) {
        return x >= 0.0 ? 1.0 / (1.0 + std::exp(x))
                        : std::exp(x) / (1.0 + std::exp(x));
    };
    const Real electronArg = (psi(0) - psi(1)) / (2.0 * Vt);
    const Real electronCarrier =
        n(0) * aux2(electronArg) + n(1) * aux2(-electronArg);
    const Real holeArg = (psi(1) - psi(0)) / (2.0 * Vt);
    const Real holeCarrier =
        p(0) * aux2(holeArg) + p(1) * aux2(-holeArg);
    const Real expectedElectronFlux =
        it->electronMobility * electronCarrier * std::abs((phin(1) - phin(0)) / it->edgeLength);
    const Real expectedHoleFlux =
        it->holeMobility * holeCarrier * std::abs((phip(1) - phip(0)) / it->edgeLength);

    REQUIRE(it->electronImpactField == Catch::Approx(std::abs((phin(1) - phin(0)) / it->edgeLength)));
    REQUIRE(it->holeImpactField == Catch::Approx(std::abs((phip(1) - phip(0)) / it->edgeLength)));
    REQUIRE(it->electronFluxProxy == Catch::Approx(expectedElectronFlux).epsilon(1.0e-12));
    REQUIRE(it->holeFluxProxy == Catch::Approx(expectedHoleFlux).epsilon(1.0e-12));
}

TEST_CASE("Cell-current reconstructed avalanche support uses cell-smoothed SG current magnitude",
          "[impact][cell_current_reconstructed]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_current_reconstructed";
    impactConfig.sourceVolumePolicy = "edge_box";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impactConfig, "test"));
    REQUIRE(detail::usesEdgeCurrentAvalancheSource(impactConfig));

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);

    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;

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
        constants::kb * constants::T0 / constants::q);

    std::vector<Real> rawElectron(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawHole(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    for (const auto& record : records) {
        rawElectron[static_cast<std::size_t>(record.edgeId)] = record.electronRawFluxProxy;
        rawHole[static_cast<std::size_t>(record.edgeId)] = record.holeRawFluxProxy;
    }
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);

    const auto it = std::find_if(records.begin(), records.end(), [](const auto& record) {
        return record.node0 == 0 && record.node1 == 1;
    });
    REQUIRE(it != records.end());
    const Real expectedElectron = detail::cellSmoothedEdgeFluxMagnitude(
        it->edgeId, rawElectron, edgeCells, cellEdges);
    const Real expectedHole = detail::cellSmoothedEdgeFluxMagnitude(
        it->edgeId, rawHole, edgeCells, cellEdges);

    REQUIRE(it->electronReconstructedFluxProxy == Catch::Approx(expectedElectron).epsilon(1.0e-12));
    REQUIRE(it->holeReconstructedFluxProxy == Catch::Approx(expectedHole).epsilon(1.0e-12));
    REQUIRE(it->electronFluxProxy == Catch::Approx(expectedElectron).epsilon(1.0e-12));
    REQUIRE(it->holeFluxProxy == Catch::Approx(expectedHole).epsilon(1.0e-12));
    REQUIRE(it->electronFinalOverRawFluxProxy ==
            Catch::Approx(expectedElectron / it->electronRawFluxProxy).epsilon(1.0e-12));
    REQUIRE(it->holeFinalOverRawFluxProxy ==
            Catch::Approx(expectedHole / it->holeRawFluxProxy).epsilon(1.0e-12));
}

TEST_CASE("Conserved total-current avalanche support uses |F_n + F_p| on both carriers",
          "[impact][conserved_total_current]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "conserved_total_current";
    impactConfig.sourceVolumePolicy = "edge_box";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impactConfig, "test"));
    REQUIRE(detail::usesEdgeCurrentAvalancheSource(impactConfig));
    REQUIRE(detail::usesConservedTotalCurrentAvalancheCurrent(impactConfig));

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);

    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;

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
        constants::kb * constants::T0 / constants::q);

    const auto it = std::find_if(records.begin(), records.end(), [](const auto& record) {
        return record.node0 == 0 && record.node1 == 1;
    });
    REQUIRE(it != records.end());

    // Both carriers must see the same conserved total-current magnitude
    // |F_n + F_p| built from the signed SG continuity fluxes.
    const Real expectedConserved = std::abs(
        it->electronRawSignedFluxProxy + it->holeRawSignedFluxProxy);
    REQUIRE(expectedConserved > 0.0);
    REQUIRE(it->electronFluxProxy == Catch::Approx(expectedConserved).epsilon(1.0e-12));
    REQUIRE(it->holeFluxProxy == Catch::Approx(expectedConserved).epsilon(1.0e-12));
    // And it must differ from the per-carrier local-density fluxes in general.
    REQUIRE(it->electronFluxProxy != Catch::Approx(it->electronRawFluxProxy).epsilon(1.0e-12));
}

TEST_CASE("Cell-vector current reconstruction recovers a constant edge-projected current",
          "[impact][cell_vector_current_reconstructed]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);

    std::vector<Real> signedFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    const Real jx = 3.0;
    const Real jy = 4.0;
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Node& n0 = mesh.getNode(edge.n0);
        const Node& n1 = mesh.getNode(edge.n1);
        signedFlux[static_cast<std::size_t>(e)] =
            jx * (n1.x - n0.x) / edge.length +
            jy * (n1.y - n0.y) / edge.length;
    }

    const auto it = std::find_if(mesh.edges().begin(), mesh.edges().end(), [](const Edge& edge) {
        return edge.n0 == 0 && edge.n1 == 1;
    });
    REQUIRE(it != mesh.edges().end());

    const Real reconstructed = detail::cellVectorReconstructedEdgeFluxMagnitude(
        it->id, signedFlux, edgeCells, cellEdges, mesh);

    REQUIRE(reconstructed == Catch::Approx(5.0).epsilon(1.0e-12));
}


TEST_CASE("Median-dual face-normal current reconstruction recovers a constant vector current",
          "[impact][dual_face_vector_mag]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);

    std::vector<Real> signedFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    const Real jx = 3.0;
    const Real jy = 4.0;
    const Cell& cell = mesh.getCell(0);
    for (int k = 0; k < 3; ++k) {
        const Index a = cell.node_ids[static_cast<std::size_t>(k)];
        const Index b = cell.node_ids[static_cast<std::size_t>((k + 1) % 3)];
        const Index edgeId = detail::edgeIdForNodePair(mesh, cellEdges[0], a, b);
        const Point2 normal = detail::medianDualFaceNormal(mesh, cell, a, b);
        const Edge& edge = mesh.getEdge(edgeId);
        const Real orientation = (edge.n0 == a && edge.n1 == b) ? 1.0 : -1.0;
        signedFlux[static_cast<std::size_t>(edgeId)] = orientation * (jx * normal.x() + jy * normal.y());
    }

    const auto it = std::find_if(mesh.edges().begin(), mesh.edges().end(), [](const Edge& edge) {
        return edge.n0 == 0 && edge.n1 == 1;
    });
    REQUIRE(it != mesh.edges().end());

    const Real reconstructed = detail::medianDualFaceVectorReconstructedEdgeFluxMagnitude(
        it->id, signedFlux, edgeCells, cellEdges, mesh);

    REQUIRE(reconstructed == Catch::Approx(5.0).epsilon(1.0e-12));
}

TEST_CASE("Cell-vector current reconstructed avalanche support uses vector SG current magnitude",
          "[impact][cell_vector_current_reconstructed]")
{
    DeviceMesh mesh = makeSingleCellMesh();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_vector_current_reconstructed";
    impactConfig.sourceVolumePolicy = "edge_box";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impactConfig, "test"));
    REQUIRE(detail::usesEdgeCurrentAvalancheSource(impactConfig));

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);

    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;

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
        constants::kb * constants::T0 / constants::q);

    std::vector<Real> signedElectron(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> signedHole(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    for (const auto& record : records) {
        signedElectron[static_cast<std::size_t>(record.edgeId)] = record.electronRawSignedFluxProxy;
        signedHole[static_cast<std::size_t>(record.edgeId)] = record.holeRawSignedFluxProxy;
    }
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);

    const auto it = std::find_if(records.begin(), records.end(), [](const auto& record) {
        return record.node0 == 0 && record.node1 == 1;
    });
    REQUIRE(it != records.end());
    const Real expectedElectron = detail::cellVectorReconstructedEdgeFluxMagnitude(
        it->edgeId, signedElectron, edgeCells, cellEdges, mesh);
    const Real expectedHole = detail::cellVectorReconstructedEdgeFluxMagnitude(
        it->edgeId, signedHole, edgeCells, cellEdges, mesh);

    REQUIRE(it->electronReconstructedFluxProxy == Catch::Approx(expectedElectron).epsilon(1.0e-12));
    REQUIRE(it->holeReconstructedFluxProxy == Catch::Approx(expectedHole).epsilon(1.0e-12));
    REQUIRE(it->electronFluxProxy == Catch::Approx(expectedElectron).epsilon(1.0e-12));
    REQUIRE(it->holeFluxProxy == Catch::Approx(expectedHole).epsilon(1.0e-12));
    REQUIRE(it->electronFinalOverRawFluxProxy ==
            Catch::Approx(expectedElectron / it->electronRawFluxProxy).epsilon(1.0e-12));
    REQUIRE(it->holeFinalOverRawFluxProxy ==
            Catch::Approx(expectedHole / it->holeRawFluxProxy).epsilon(1.0e-12));
}
