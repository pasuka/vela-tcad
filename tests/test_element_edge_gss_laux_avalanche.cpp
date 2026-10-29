#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"

#include <array>
#include <vector>

using namespace vela;

namespace {

DeviceMesh makeSingleRightTriangle()
{
    DeviceMesh mesh;
    Node n0;
    n0.id = 0;
    n0.x = 0.0;
    n0.y = 0.0;
    mesh.addNode(n0);
    Node n1;
    n1.id = 1;
    n1.x = 1.0e-6;
    n1.y = 0.0;
    mesh.addNode(n1);
    Node n2;
    n2.id = 2;
    n2.x = 0.0;
    n2.y = 0.5e-6;
    mesh.addNode(n2);

    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);

    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);

    mesh.buildEdges();
    return mesh;
}
DeviceMesh makeSingleAcuteScaleneTriangle(bool reverse)
{
    DeviceMesh mesh;
    Node n0;
    n0.id = 0;
    n0.x = 0.0;
    n0.y = 0.0;
    mesh.addNode(n0);
    Node n1;
    n1.id = 1;
    n1.x = 1.0e-6;
    n1.y = 0.0;
    mesh.addNode(n1);
    Node n2;
    n2.id = 2;
    n2.x = 0.2e-6;
    n2.y = 0.8e-6;
    mesh.addNode(n2);

    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = reverse ? std::vector<Index>{0, 2, 1}
                            : std::vector<Index>{0, 1, 2};
    mesh.addCell(cell);

    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeSingleObtuseTriangle(bool reverse)
{
    DeviceMesh mesh;
    Node n0;
    n0.id = 0;
    n0.x = 0.0;
    n0.y = 0.0;
    mesh.addNode(n0);
    Node n1;
    n1.id = 1;
    n1.x = 2.0e-6;
    n1.y = 0.0;
    mesh.addNode(n1);
    Node n2;
    n2.id = 2;
    n2.x = 0.2e-6;
    n2.y = 0.1e-6;
    mesh.addNode(n2);

    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = reverse ? std::vector<Index>{0, 2, 1}
                            : std::vector<Index>{0, 1, 2};
    mesh.addCell(cell);

    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeContactAndInteriorTriangles()
{
    DeviceMesh mesh;
    const std::array<Point2, 4> points = {
        Point2{0.0, 0.0}, Point2{1.0e-6, 0.0},
        Point2{0.0, 1.0e-6}, Point2{1.0e-6, 1.0e-6}};
    for (Index id = 0; id < points.size(); ++id) {
        Node node;
        node.id = id;
        node.x = points[id].x();
        node.y = points[id].y();
        mesh.addNode(node);
    }

    Cell contactCell;
    contactCell.id = 0;
    contactCell.type = CellType::Tri3;
    contactCell.region_id = 0;
    contactCell.node_ids = {0, 1, 2};
    mesh.addCell(contactCell);

    Cell interiorCell;
    interiorCell.id = 1;
    interiorCell.type = CellType::Tri3;
    interiorCell.region_id = 0;
    interiorCell.node_ids = {1, 3, 2};
    mesh.addCell(interiorCell);

    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0, 1};
    mesh.addRegion(region);

    Contact contact;
    contact.id = 0;
    contact.name = "contact";
    contact.region_id = 0;
    contact.node_ids = {0, 1};
    mesh.addContact(contact);

    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("Element-edge GSS Laux mode has an explicit canonical contract",
          "[impact][element_edge_gss_laux]")
{
    ImpactIonizationModelConfig config;
    config.model = "van_overstraeten";
    config.drivingForce = "quasi_fermi_gradient";
    config.generation = "current_density";
    config.currentApproximation = "element_edge_sg_gss_laux";
    config.quasiFermiGradientDiscretization = "cell_gradient";
    config.sourceMappingMode = "element_vertex_box_measure";

    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(config, "test"));
    REQUIRE(detail::usesElementEdgeGssLauxAvalancheSource(config));
    config.drivingForce = "electric_field";
    config.quasiFermiGradientDiscretization = "edge_difference";
    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(config, "test"));

    config.sourceMappingMode = "node_F_node_alpha_node_G";
    REQUIRE_THROWS_WITH(
        detail::validateImpactIonizationDrivingForce(config, "test"),
        Catch::Matchers::ContainsSubstring(
            "element_edge_sg_gss_laux requires the canonical element-box configuration"));

    config.sourceMappingMode = "element_vertex_box_measure";
    config.drivingForce = "quasi_fermi_gradient";
    config.quasiFermiGradientDiscretization = "cell_gradient";
    config.contactElectricFieldFallbackMode = "face_normal";
    REQUIRE_THROWS_WITH(
        detail::validateImpactIonizationDrivingForce(config, "test"),
        Catch::Matchers::ContainsSubstring(
            "require scope 'contact_boundary_face'"));
    config.contactElectricFieldFallbackScope = "contact_boundary_face";
    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(config, "test"));
    config.contactElectricFieldFallbackMode = "unsupported";
    REQUIRE_THROWS_WITH(
        detail::validateImpactIonizationDrivingForce(config, "test"),
        Catch::Matchers::ContainsSubstring(
            "contact_electric_field_fallback_mode must be"));
}

TEST_CASE("A zero-box diagonal is inactive on a right triangle",
          "[impact][element_edge_gss_laux][geometry][right_triangle]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    const Cell& cell = mesh.getCell(0);
    const Point2 expected{3.0, -4.0};
    std::array<Real, 3> signedEdgeCurrent{};

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 =
            cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        const Point2 delta =
            detail::meshPoint(mesh, node1) - detail::meshPoint(mesh, node0);
        const Point2 tangent = delta / delta.norm();
        signedEdgeCurrent[static_cast<std::size_t>(localEdge)] =
            tangent.dot(expected);
    }

    const auto partialVolumes =
        detail::tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    REQUIRE(partialVolumes[1] == Catch::Approx(0.0).margin(1.0e-30));
    REQUIRE(signedEdgeCurrent[1] != Catch::Approx(0.0));

    const Point2 reconstructed = detail::gssLauxTri3CurrentVector(
        mesh, cell, signedEdgeCurrent);
    REQUIRE(reconstructed.x() ==
            Catch::Approx(expected.x()).epsilon(1.0e-13));
    REQUIRE(reconstructed.y() ==
            Catch::Approx(expected.y()).epsilon(1.0e-13));

    auto perturbedDiagonalCurrent = signedEdgeCurrent;
    perturbedDiagonalCurrent[1] += 123.0;
    const Point2 perturbed = detail::gssLauxTri3CurrentVector(
        mesh, cell, perturbedDiagonalCurrent);
    REQUIRE(perturbed.x() ==
            Catch::Approx(reconstructed.x()).epsilon(1.0e-13));
    REQUIRE(perturbed.y() ==
            Catch::Approx(reconstructed.y()).epsilon(1.0e-13));

    const auto vertexMeasures =
        detail::tri3ElementVertexBoxMeasures(mesh, cell);
    const Real measureSum =
        vertexMeasures[0] + vertexMeasures[1] + vertexMeasures[2];
    REQUIRE(measureSum == Catch::Approx(0.25e-12).epsilon(1.0e-13));
}

TEST_CASE("GSS Laux reconstruction is orientation invariant on a scalene triangle",
          "[impact][element_edge_gss_laux][geometry][general_mesh]")
{
    const Point2 expected{-2.75, 6.5};
    for (const bool reverse : {false, true}) {
        DeviceMesh mesh = makeSingleAcuteScaleneTriangle(reverse);
        const Cell& cell = mesh.getCell(0);
        std::array<Real, 3> signedEdgeCurrent{};
        for (int localEdge = 0; localEdge < 3; ++localEdge) {
            const Index node0 =
                cell.node_ids[static_cast<std::size_t>(localEdge)];
            const Index node1 =
                cell.node_ids[
                    static_cast<std::size_t>((localEdge + 1) % 3)];
            const Point2 delta =
                detail::meshPoint(mesh, node1) -
                detail::meshPoint(mesh, node0);
            signedEdgeCurrent[static_cast<std::size_t>(localEdge)] =
                (delta / delta.norm()).dot(expected);
        }

        const Point2 reconstructed = detail::gssLauxTri3CurrentVector(
            mesh, cell, signedEdgeCurrent);
        REQUIRE(reconstructed.x() ==
                Catch::Approx(expected.x()).epsilon(1.0e-13));
        REQUIRE(reconstructed.y() ==
                Catch::Approx(expected.y()).epsilon(1.0e-13));

        for (int localEdge = 0; localEdge < 3; ++localEdge) {
            auto perturbedEdgeCurrent = signedEdgeCurrent;
            perturbedEdgeCurrent[static_cast<std::size_t>(localEdge)] +=
                123.0;
            const Point2 perturbed = detail::gssLauxTri3CurrentVector(
                mesh, cell, perturbedEdgeCurrent);
            REQUIRE((perturbed - reconstructed).norm() >
                    1.0e-6);
        }

        const auto partialVolumes =
            detail::tri3ElementEdgeBoxPartialVolumes(mesh, cell);
        REQUIRE(partialVolumes[0] > 0.0);
        REQUIRE(partialVolumes[1] > 0.0);
        REQUIRE(partialVolumes[2] > 0.0);

        const auto vertexMeasures =
            detail::tri3ElementVertexBoxMeasures(mesh, cell);
        const Real measureSum =
            vertexMeasures[0] + vertexMeasures[1] + vertexMeasures[2];
        REQUIRE(measureSum ==
                Catch::Approx(0.4e-12).epsilon(1.0e-13));
        REQUIRE(partialVolumes[0] + partialVolumes[1] +
                    partialVolumes[2] ==
                Catch::Approx(0.4e-12).epsilon(1.0e-13));
    }
}

TEST_CASE("Obtuse element-edge box support is nonnegative and area conservative",
          "[impact][element_edge_gss_laux][geometry][general_mesh][obtuse]")
{
    constexpr Real expectedArea = 0.1e-12;
    for (const bool reverse : {false, true}) {
        DeviceMesh mesh = makeSingleObtuseTriangle(reverse);
        const Cell& cell = mesh.getCell(0);
        const auto partialVolumes =
            detail::tri3ElementEdgeBoxPartialVolumes(mesh, cell);
        const auto vertexMeasures =
            detail::tri3ElementVertexBoxMeasures(mesh, cell);

        for (const Real partialVolume : partialVolumes)
            REQUIRE(partialVolume >= 0.0);
        for (const Real vertexMeasure : vertexMeasures)
            REQUIRE(vertexMeasure >= 0.0);

        const Real partialVolumeSum =
            partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
        const Real vertexMeasureSum =
            vertexMeasures[0] + vertexMeasures[1] + vertexMeasures[2];
        REQUIRE(partialVolumeSum ==
                Catch::Approx(expectedArea).epsilon(1.0e-12));
        REQUIRE(vertexMeasureSum ==
                Catch::Approx(expectedArea).epsilon(1.0e-12));
    }
}

TEST_CASE("Element-edge QFP driver remains global in contact and interior cells",
          "[impact][element_edge_gss_laux][driver][contact][interior]")
{
    DeviceMesh mesh = makeContactAndInteriorTriangles();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    const auto impact = makeImpactIonizationModel(impactConfig);
    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);

    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Node& point = mesh.getNode(node);
        psi(static_cast<int>(node)) = point.x * 1.0e6;
        phin(static_cast<int>(node)) =
            point.x * 3.0e6 + point.y * 4.0e6;
        phip(static_cast<int>(node)) =
            point.x * -2.0e6 + point.y * 5.0e6;
    }
    const std::vector<Real> ni(
        static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    const Real Vt = constants::kb * constants::T0 / constants::q;
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const int index = static_cast<int>(node);
        n(index) = ni[node] * std::exp((psi(index) - phin(index)) / Vt);
        p(index) = ni[node] * std::exp((phip(index) - psi(index)) / Vt);
    }

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const auto record =
            detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                impactConfig, *impact, mobilityConfig, *mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh, doping,
                cellMaterials, cellId, psi, phin, phip, n, p, ni, Vt);
        REQUIRE(record.electronImpactField ==
                Catch::Approx(5.0e6).epsilon(1.0e-13));
        REQUIRE(record.holeImpactField ==
                Catch::Approx(std::sqrt(29.0) * 1.0e6).epsilon(1.0e-13));
        REQUIRE(record.electronImpactField != Catch::Approx(1.0e6));
        REQUIRE(record.holeImpactField != Catch::Approx(1.0e6));
    }

    impactConfig.contactElectricFieldFallback = true;
    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(impactConfig, "test"));
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const auto record =
            detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                impactConfig, *impact, mobilityConfig, *mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh, doping,
                cellMaterials, cellId, psi, phin, phip, n, p, ni, Vt);
        const Real expectedElectron = 1.0e6;
        const Real expectedHole = 1.0e6;
        REQUIRE(record.electronImpactField ==
                Catch::Approx(expectedElectron).epsilon(1.0e-13));
        REQUIRE(record.holeImpactField ==
                Catch::Approx(expectedHole).epsilon(1.0e-13));
    }

    impactConfig.contactElectricFieldFallbackScope = "contact_boundary_face";
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const auto record =
            detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                impactConfig, *impact, mobilityConfig, *mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh, doping,
                cellMaterials, cellId, psi, phin, phip, n, p, ni, Vt);
        const Real expectedElectron = cellId == 0 ? 1.0e6 : 5.0e6;
        const Real expectedHole = cellId == 0
            ? 1.0e6 : std::sqrt(29.0) * 1.0e6;
        REQUIRE(record.electronImpactField ==
                Catch::Approx(expectedElectron).epsilon(1.0e-13));
        REQUIRE(record.holeImpactField ==
                Catch::Approx(expectedHole).epsilon(1.0e-13));
    }

    struct ModeExpectation {
        const char* mode;
        Real electronField;
        Real holeField;
    };
    const std::array<ModeExpectation, 3> expectations{{
        {"face_normal", 0.0, 0.0},
        {"one_sided", 1.0e6 / std::sqrt(5.0),
         1.0e6 / std::sqrt(5.0)},
        {"distance_weighted_blend", 5.0e6 / 3.0,
         std::sqrt(29.0) * 1.0e6 / 3.0}}};
    for (const auto& expectation : expectations) {
        impactConfig.contactElectricFieldFallbackMode = expectation.mode;
        const auto contactRecord =
            detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                impactConfig, *impact, mobilityConfig, *mobility,
                cellEdges.at(0), mesh, doping, cellMaterials, 0,
                psi, phin, phip, n, p, ni, Vt);
        REQUIRE(contactRecord.electronImpactField ==
                Catch::Approx(expectation.electronField).epsilon(1.0e-13));
        REQUIRE(contactRecord.holeImpactField ==
                Catch::Approx(expectation.holeField).epsilon(1.0e-13));
        const auto interiorRecord =
            detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                impactConfig, *impact, mobilityConfig, *mobility,
                cellEdges.at(1), mesh, doping, cellMaterials, 1,
                psi, phin, phip, n, p, ni, Vt);
        REQUIRE(interiorRecord.electronImpactField ==
                Catch::Approx(5.0e6).epsilon(1.0e-13));
        REQUIRE(interiorRecord.holeImpactField ==
                Catch::Approx(std::sqrt(29.0) * 1.0e6).epsilon(1.0e-13));
    }
}

TEST_CASE("Element-edge GSS Laux records use exact SG currents and box measures",
          "[impact][element_edge_gss_laux][source]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    const auto impact = makeImpactIonizationModel(impactConfig);

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;
    const std::vector<Real> ni(
        static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    const Real Vt = constants::kb * constants::T0 / constants::q;

    const auto record =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, mobilityConfig, *mobility, cellEdges.at(0),
            mesh, doping, cellMaterials, 0, psi, phin, phip, n, p, ni, Vt);

    REQUIRE(record.cellId == 0);
    REQUIRE(record.vertexMeasures[0] + record.vertexMeasures[1] +
                record.vertexMeasures[2] ==
            Catch::Approx(0.25e-12).epsilon(1.0e-13));
    REQUIRE(record.edgePartialVolumes[1] ==
            Catch::Approx(0.0).margin(1.0e-30));
    REQUIRE(record.electronSignedEdgeFlux[1] != Catch::Approx(0.0));
    REQUIRE(record.holeSignedEdgeFlux[1] != Catch::Approx(0.0));

    const Real totalMeasure = 0.25e-12;
    REQUIRE(record.electronCurrentVector.norm() > 0.0);
    REQUIRE(record.holeCurrentVector.norm() > 0.0);
    REQUIRE(record.electronSourceIntegrals[0] +
                record.electronSourceIntegrals[1] +
                record.electronSourceIntegrals[2] ==
            Catch::Approx(record.electronAlpha *
                          record.electronCurrentVector.norm() *
                          totalMeasure)
                .epsilon(1.0e-13));
    REQUIRE(record.holeSourceIntegrals[0] +
                record.holeSourceIntegrals[1] +
                record.holeSourceIntegrals[2] ==
            Catch::Approx(record.holeAlpha *
                          record.holeCurrentVector.norm() *
                          totalMeasure)
                .epsilon(1.0e-13));

    const auto components =
        detail::currentDensityAvalancheSourceComponentIntegrals(
            impactConfig, *impact, mobilityConfig, *mobility, edgeCells,
            mesh, doping, cellMaterials, psi, phin, phip, n, p, ni, Vt);
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        const Index node = mesh.getCell(0).node_ids[localNode];
        REQUIRE(components.electron[node] ==
                Catch::Approx(record.electronSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.hole[node] ==
                Catch::Approx(record.holeSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.combined[node] ==
                Catch::Approx(record.combinedSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.electron[node] + components.hole[node] ==
                Catch::Approx(components.combined[node])
                    .epsilon(1.0e-13));
    }
    impactConfig.drivingForce = "electric_field";
    const auto electricRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, mobilityConfig, *mobility, cellEdges.at(0),
            mesh, doping, cellMaterials, 0, psi, phin, phip, n, p, ni, Vt);
    const Real expectedElectricField = std::hypot(2.0e5, 1.6e5);
    REQUIRE(electricRecord.electronImpactField ==
            Catch::Approx(expectedElectricField).epsilon(1.0e-13));
    REQUIRE(electricRecord.holeImpactField ==
            Catch::Approx(expectedElectricField).epsilon(1.0e-13));
    REQUIRE(electricRecord.electronImpactField !=
            Catch::Approx(record.electronImpactField));

    MobilityModelConfig fieldMobilityConfig =
        mobilityModelConfig("masetti_field");
    fieldMobilityConfig.highFieldDrivingForce = "electric_field";
    const auto fieldMobility = makeMobilityModel(fieldMobilityConfig);
    const auto electricMobilityRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, fieldMobilityConfig, *fieldMobility,
            cellEdges.at(0), mesh, doping, cellMaterials, 0, psi, phin, phip,
            n, p, ni, Vt);

    fieldMobilityConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    const auto qfMobilityRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, fieldMobilityConfig, *fieldMobility,
            cellEdges.at(0), mesh, doping, cellMaterials, 0, psi, phin, phip,
            n, p, ni, Vt);

    REQUIRE(electricMobilityRecord.electronMobilities[0] >
            qfMobilityRecord.electronMobilities[0]);
    REQUIRE(electricMobilityRecord.holeMobilities[0] >
            qfMobilityRecord.holeMobilities[0]);
    REQUIRE(electricMobilityRecord.electronSignedEdgeFlux[0] !=
            Catch::Approx(qfMobilityRecord.electronSignedEdgeFlux[0]));
    REQUIRE(electricMobilityRecord.holeSignedEdgeFlux[0] !=
            Catch::Approx(qfMobilityRecord.holeSignedEdgeFlux[0]));
}

TEST_CASE("Element-edge avalanche composes with Enhanced Lombardi mobility",
          "[impact][element_edge_gss_laux][lombardi]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto materials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    const auto impact = makeImpactIonizationModel(impactConfig);

    MobilityModelConfig bulkConfig = mobilityModelConfig("masetti_field");
    bulkConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    MobilityModelConfig lombardiConfig = bulkConfig;
    lombardiConfig.model = "masetti_field_lombardi";
    lombardiConfig.surface.surfaceRegion = "si";
    const auto bulkMobility = makeMobilityModel(bulkConfig);
    const auto lombardiMobility = makeMobilityModel(lombardiConfig);

    VectorXd psi(3), phin(3), phip(3), n(3), p(3);
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;
    const std::vector<Real> ni(3, 1.0e16);
    const Real Vt = constants::kb * constants::T0 / constants::q;

    const auto bulkRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, bulkConfig, *bulkMobility, cellEdges.at(0),
            mesh, doping, materials, 0, psi, phin, phip, n, p, ni, Vt);
    const auto lombardiRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, lombardiConfig, *lombardiMobility,
            cellEdges.at(0), mesh, doping, materials, 0, psi, phin, phip,
            n, p, ni, Vt);

    REQUIRE(lombardiRecord.electronMobilities[0] <
            bulkRecord.electronMobilities[0]);
    REQUIRE(lombardiRecord.holeMobilities[0] <
            bulkRecord.holeMobilities[0]);
    REQUIRE(lombardiRecord.combinedSourceIntegrals[0] > 0.0);
    REQUIRE(lombardiRecord.combinedSourceIntegrals[0] <
            bulkRecord.combinedSourceIntegrals[0]);
}

TEST_CASE("Contact-face fallback driver modes have consistent local AD derivatives",
          "[impact][element_edge_gss_laux][ad][contact_face]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    Contact contact;
    contact.id = 0;
    contact.name = "contact";
    contact.region_id = 0;
    contact.node_ids = {0, 1};
    mesh.addContact(contact);
    const Cell& cell = mesh.getCell(0);

    ImpactIonizationModelConfig config;
    config.contactElectricFieldFallbackScope = "contact_boundary_face";
    constexpr Real fieldFactor = 100.0;
    const std::array<Real, 3> basePsi{0.0, -0.10, 0.04};
    const std::array<Real, 3> basePhin{0.01, 0.01, 0.20};

    for (const std::string mode : {
             "face_normal", "one_sided", "distance_weighted_blend"}) {
        config.contactElectricFieldFallbackMode = mode;
        std::array<detail::Tri3LocalForwardDual, 3> psi{};
        std::array<detail::Tri3LocalForwardDual, 3> phin{};
        for (std::size_t local = 0; local < 3; ++local) {
            psi[local] = detail::Tri3LocalForwardDual::variable(
                basePsi[local], local);
            phin[local] = detail::Tri3LocalForwardDual::variable(
                basePhin[local], 3 + local);
        }
        const auto electricGradient =
            detail::localAdTri3Gradient(mesh, cell, psi);
        const auto qfGradient =
            detail::localAdTri3Gradient(mesh, cell, phin);
        const auto qfField = detail::localAdNorm2(
            qfGradient[0], qfGradient[1]) *
            detail::Tri3LocalForwardDual(fieldFactor);
        const auto result =
            detail::localAdContactElectricFallbackImpactField(
                config, mesh, cell, electricGradient, qfField,
                psi, fieldFactor);

        const auto evaluate = [&](const std::array<Real, 3>& psiValues,
                                  const std::array<Real, 3>& phinValues) {
            const auto electric =
                detail::localAdTri3Gradient(mesh, cell, psiValues);
            const auto qf =
                detail::localAdTri3Gradient(mesh, cell, phinValues);
            const Real qfMagnitude =
                detail::localAdNorm2(qf[0], qf[1]) * fieldFactor;
            const Point2 electricPoint{electric[0], electric[1]};
            return detail::contactElectricFallbackImpactField(
                config, mesh, cell, electricPoint, qfMagnitude,
                [&](Index node) {
                    return psiValues[static_cast<std::size_t>(node)];
                }, fieldFactor);
        };
        const Real base = evaluate(basePsi, basePhin);
        CAPTURE(mode, base, result.value);
        REQUIRE(result.value == Catch::Approx(base).epsilon(1.0e-13));
        for (std::size_t dof = 0; dof < 6; ++dof) {
            auto plusPsi = basePsi;
            auto minusPsi = basePsi;
            auto plusPhin = basePhin;
            auto minusPhin = basePhin;
            const std::size_t local = dof % 3;
            const Real step = 1.0e-7;
            if (dof < 3) {
                plusPsi[local] += step;
                minusPsi[local] -= step;
            } else {
                plusPhin[local] += step;
                minusPhin[local] -= step;
            }
            const Real finiteDifference =
                (evaluate(plusPsi, plusPhin) -
                 evaluate(minusPsi, minusPhin)) / (2.0 * step);
            CAPTURE(mode, dof, finiteDifference, result.derivative[dof]);
            REQUIRE(result.derivative[dof] ==
                    Catch::Approx(finiteDifference).epsilon(2.0e-8)
                        .margin(1.0e-4));
        }
    }
}

TEST_CASE("Local forward AD source derivatives converge independently",
          "[impact][element_edge_gss_laux][ad]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    Contact contact;
    contact.id = 0;
    contact.name = "contact";
    contact.region_id = 0;
    contact.node_ids = {0, 1};
    mesh.addContact(contact);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig =
        impactIonizationModelConfig("van_overstraeten");
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    impactConfig.contactElectricFieldFallback = true;
    impactConfig.contactElectricFieldFallbackScope = "contact_boundary_face";
    MobilityModelConfig mobilityConfig = mobilityModelConfig("masetti_field");
    mobilityConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    mobilityConfig.highFieldGradientDiscretization = GENERATE(
        std::string("edge_projection"), std::string("transport_cell_vector"));
    CAPTURE(mobilityConfig.highFieldGradientDiscretization);
    const auto mobility = makeMobilityModel(mobilityConfig);

    constexpr Real Vt = 0.025852;
    constexpr Real fieldFactor = 100.0;
    const std::array<Real, 3> basePsi{0.0, -0.10, 0.04};
    // Edge 0 is exactly flat in both QFPs; the third vertex keeps the cell
    // QFP-gradient drive active and verifies that AD does not lose the
    // derivative at the production value short-circuit.
    const std::array<Real, 3> basePhin{0.01, 0.01, 0.20};
    const std::array<Real, 3> basePhip{-0.01, -0.01, -0.20};
    const std::array<Real, 3> intrinsicDensity{1.0e16, 1.0e16, 1.0e16};

    const auto evaluateReal = [&](const std::array<Real, 3>& psi,
                                  const std::array<Real, 3>& phin,
                                  const std::array<Real, 3>& phip) {
        std::array<Real, 3> n{};
        std::array<Real, 3> p{};
        for (std::size_t local = 0; local < 3; ++local) {
            n[local] = intrinsicDensity[local] *
                std::exp((psi[local] - phin[local]) / Vt);
            p[local] = intrinsicDensity[local] *
                std::exp((phip[local] - psi[local]) / Vt);
        }
        return detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<Real>(
            impactConfig, mobilityConfig, *mobility, cellEdges.at(0), mesh,
            doping, cellMaterials, 0, psi, phin, phip, n, p,
            intrinsicDensity, Vt, fieldFactor).combined;
    };

    std::array<detail::Tri3LocalForwardDual, 3> psi{};
    std::array<detail::Tri3LocalForwardDual, 3> phin{};
    std::array<detail::Tri3LocalForwardDual, 3> phip{};
    std::array<detail::Tri3LocalForwardDual, 3> n{};
    std::array<detail::Tri3LocalForwardDual, 3> p{};
    for (std::size_t local = 0; local < 3; ++local) {
        psi[local] = detail::Tri3LocalForwardDual::variable(basePsi[local], local);
        phin[local] =
            detail::Tri3LocalForwardDual::variable(basePhin[local], 3 + local);
        phip[local] =
            detail::Tri3LocalForwardDual::variable(basePhip[local], 6 + local);
        n[local] = detail::Tri3LocalForwardDual(intrinsicDensity[local]) *
            detail::localAdLimitedExp(
                (psi[local] - phin[local]) /
                detail::Tri3LocalForwardDual(Vt));
        p[local] = detail::Tri3LocalForwardDual(intrinsicDensity[local]) *
            detail::localAdLimitedExp(
                (phip[local] - psi[local]) /
                detail::Tri3LocalForwardDual(Vt));
    }
    const auto adSource =
        detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<
            detail::Tri3LocalForwardDual>(
            impactConfig, mobilityConfig, *mobility, cellEdges.at(0), mesh,
            doping, cellMaterials, 0, psi, phin, phip, n, p,
            intrinsicDensity, Vt, fieldFactor).combined;

    const auto baseSource = evaluateReal(basePsi, basePhin, basePhip);
    VectorXd recordPsi(3), recordPhin(3), recordPhip(3), recordN(3), recordP(3);
    for (std::size_t local = 0; local < 3; ++local) {
        recordPsi(local) = basePsi[local];
        recordPhin(local) = basePhin[local];
        recordPhip(local) = basePhip[local];
        recordN(local) = intrinsicDensity[local] * std::exp((basePsi[local] - basePhin[local]) / Vt);
        recordP(local) = intrinsicDensity[local] * std::exp((basePhip[local] - basePsi[local]) / Vt);
    }
    const auto impact = makeImpactIonizationModel(impactConfig);
    const auto record = detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
        impactConfig, *impact, mobilityConfig, *mobility, cellEdges.at(0), mesh,
        doping, cellMaterials, 0, recordPsi, recordPhin, recordPhip, recordN,
        recordP, std::vector<Real>(intrinsicDensity.begin(), intrinsicDensity.end()), Vt, fieldFactor);
    for (std::size_t localRow = 0; localRow < 3; ++localRow) {
        REQUIRE(adSource[localRow].value ==
                Catch::Approx(baseSource[localRow]).epsilon(1.0e-13));
        // Independently reconstruct generation from the selected mobility
        // drive and GSS/Laux current, rather than comparing two AD calls.
        const Real generation = record.electronAlpha * record.electronCurrentVector.norm() +
                                record.holeAlpha * record.holeCurrentVector.norm();
        REQUIRE(baseSource[localRow] == Catch::Approx(
            generation * record.vertexMeasures[localRow]).epsilon(1.0e-12));
    }

    const auto errorAt = [&](Real relativeStep) {
        Real maxDifference = 0.0;
        Real maxReference = 0.0;
        for (std::size_t localDof = 0;
             localDof < detail::Tri3LocalPotentialDofCount; ++localDof) {
            std::array<Real, 3> plusPsi = basePsi;
            std::array<Real, 3> minusPsi = basePsi;
            std::array<Real, 3> plusPhin = basePhin;
            std::array<Real, 3> minusPhin = basePhin;
            std::array<Real, 3> plusPhip = basePhip;
            std::array<Real, 3> minusPhip = basePhip;
            const std::size_t block = localDof / 3;
            const std::size_t local = localDof % 3;
            const Real value = block == 0 ? basePsi[local]
                : (block == 1 ? basePhin[local] : basePhip[local]);
            const Real step = relativeStep * std::max<Real>(1.0, std::abs(value));
            auto* plus = block == 0 ? &plusPsi : (block == 1 ? &plusPhin : &plusPhip);
            auto* minus = block == 0 ? &minusPsi : (block == 1 ? &minusPhin : &minusPhip);
            (*plus)[local] += step;
            (*minus)[local] -= step;
            const auto plusSource = evaluateReal(plusPsi, plusPhin, plusPhip);
            const auto minusSource = evaluateReal(minusPsi, minusPhin, minusPhip);
            for (std::size_t localRow = 0; localRow < 3; ++localRow) {
                const Real finiteDifference =
                    (plusSource[localRow] - minusSource[localRow]) /
                    (2.0 * step);
                const Real analytic = adSource[localRow].derivative[localDof];
                maxDifference = std::max(
                    maxDifference, std::abs(analytic - finiteDifference));
                maxReference = std::max(maxReference, std::abs(finiteDifference));
            }
        }
        return maxDifference / maxReference;
    };

    const Real coarseError = errorAt(1.0e-6);
    const Real mediumError = errorAt(3.0e-7);
    const Real fineError = errorAt(1.0e-7);
    CAPTURE(coarseError, mediumError, fineError);
    REQUIRE(coarseError > mediumError);
    REQUIRE(mediumError > fineError);
    REQUIRE(fineError <= 1.0e-8);
}
