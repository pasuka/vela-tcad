#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using namespace vela;

namespace {

nlohmann::json readJson(const std::filesystem::path& path)
{
    std::ifstream input(path);
    REQUIRE(input.is_open());
    nlohmann::json json;
    input >> json;
    return json;
}

DopingModel dopingFromDeck(const DeviceMesh& mesh, const nlohmann::json& cfg)
{
    std::vector<RegionDopingSpec> specs;
    for (const auto& entry : cfg.at("doping")) {
        specs.push_back({
            entry.at("region").get<std::string>(),
            entry.at("donors").get<Real>(),
            entry.at("acceptors").get<Real>(),
        });
    }
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

bool vectorIsFinite(const VectorXd& values)
{
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values(i)))
            return false;
    }
    return true;
}

std::vector<std::vector<Index>> buildNodeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeCells(mesh.numNodes());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        for (Index nodeId : mesh.getCell(cellId).node_ids)
            nodeCells.at(nodeId).push_back(cellId);
    }
    return nodeCells;
}

std::vector<Index> contactNodesOnlyInRegion(const DeviceMesh& mesh,
                                             const std::string& contactName,
                                             const std::string& regionName)
{
    const auto nodeCells = buildNodeCellMap(mesh);
    std::vector<Index> nodes;

    for (const Contact& contact : mesh.contacts()) {
        if (contact.name != contactName)
            continue;

        for (Index nodeId : contact.node_ids) {
            bool hasCell = false;
            bool onlyRegion = true;
            for (Index cellId : nodeCells.at(nodeId)) {
                hasCell = true;
                const Region& region = mesh.getRegion(mesh.getCell(cellId).region_id);
                if (region.name != regionName) {
                    onlyRegion = false;
                    break;
                }
            }
            if (hasCell && onlyRegion)
                nodes.push_back(nodeId);
        }
        break;
    }

    return nodes;
}

Index firstSemiconductorOxideInterfaceEdge(const DeviceMesh& mesh,
                                           const MaterialDatabase& matdb)
{
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        bool hasTransportCell = false;
        bool hasInsulatingCell = false;
        for (Index cellId : edgeCells.at(edgeId)) {
            const Region& region = mesh.getRegion(mesh.getCell(cellId).region_id);
            const Material material = matdb.getMaterial(region.material);
            hasTransportCell =
                hasTransportCell || material.mun > 0.0 || material.mup > 0.0;
            hasInsulatingCell =
                hasInsulatingCell || (material.mun <= 0.0 && material.mup <= 0.0);
        }
        if (hasTransportCell && hasInsulatingCell)
            return edgeId;
    }

    FAIL("expected the MOS mesh to contain a semiconductor/oxide interface edge");
    return 0;
}

std::vector<std::filesystem::path> mosExampleDirs()
{
    const std::filesystem::path examplesRoot =
        std::filesystem::path(VELA_SOURCE_DIR) / "examples";
    return {
        examplesRoot / "nmos2d_mos_dd",
        examplesRoot / "pmos2d_mos_dd",
    };
}

DeviceMesh makeHorizontalInterfaceMesh()
{
    DeviceMesh mesh;
    const Real L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;  mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;  mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5 * L; n2.y = -L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.5 * L; n3.y = L;    mesh.addNode(n3);

    Cell si; si.id = 0; si.type = CellType::Tri3; si.region_id = 0;
    si.node_ids = {0, 1, 2}; mesh.addCell(si);
    Cell oxide; oxide.id = 1; oxide.type = CellType::Tri3; oxide.region_id = 1;
    oxide.node_ids = {0, 3, 1}; mesh.addCell(oxide);

    Region channel; channel.id = 0; channel.name = "channel"; channel.material = "Si";
    channel.cell_ids = {0}; mesh.addRegion(channel);
    Region gateOxide; gateOxide.id = 1; gateOxide.name = "gate_oxide"; gateOxide.material = "SiO2";
    gateOxide.cell_ids = {1}; mesh.addRegion(gateOxide);

    mesh.buildEdges();
    return mesh;
}

Index findEdgeByNodes(const DeviceMesh& mesh, Index a, Index b)
{
    if (b < a)
        std::swap(a, b);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        Index e0 = edge.n0;
        Index e1 = edge.n1;
        if (e1 < e0)
            std::swap(e0, e1);
        if (e0 == a && e1 == b)
            return edgeId;
    }
    FAIL("expected edge in test mesh");
    return 0;
}

} // namespace


TEST_CASE("surface mobility uses the reconstructed normal interface field",
          "[mos_mixed][mobility][surface]")
{
    DeviceMesh mesh = makeHorizontalInterfaceMesh();
    MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 0.0, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const Real temperature_K = constants::Vt_300 * constants::q / constants::kb;
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);

    MobilityModelConfig config = mobilityModelConfig("caughey_thomas_surface");
    config.surface.thetaElectron = 2.0e-6;
    config.surface.surfaceRegion = "channel";
    config.surface.surfaceInterface = {"channel", "gate_oxide"};
    const auto mobility = makeMobilityModel(config);

    VectorXd psi(4);
    psi << 0.0, 0.0, 0.0, 1.0;
    const Index interfaceEdge = findEdgeByNodes(mesh, 0, 1);

    const Real tangentialOnly = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config);
    const Real normalLimited = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi);

    REQUIRE(tangentialOnly > 0.0);
    REQUIRE(normalLimited > 0.0);
    REQUIRE(normalLimited < tangentialOnly);
}

TEST_CASE("mixed Si/SiO2 MOS edge mobility preserves semiconductor interface transport",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            const auto edgeCells = detail::buildEdgeCellMap(mesh);
            const Real temperature_K = constants::Vt_300 * constants::q / constants::kb;
            const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);
            const ConstantMobility mobility;
            const Index interfaceEdge = firstSemiconductorOxideInterfaceEdge(mesh, matdb);

            const Real mun = detail::edgeMobility(edgeCells,
                                                  mesh,
                                                  doping,
                                                  mobility,
                                                  cellMaterials,
                                                  interfaceEdge,
                                                  CarrierType::Electron,
                                                  0.0);
            const Real mup = detail::edgeMobility(edgeCells,
                                                  mesh,
                                                  doping,
                                                  mobility,
                                                  cellMaterials,
                                                  interfaceEdge,
                                                  CarrierType::Hole,
                                                  0.0);

            REQUIRE(mun > 0.0);
            REQUIRE(mup > 0.0);
        }
    }
}

TEST_CASE("mixed Si/SiO2 MOS DD scalar assembly keeps oxide carrier rows finite",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            DDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

            const int nNodes = static_cast<int>(mesh.numNodes());
            VectorXd psi = VectorXd::Zero(nNodes);
            VectorXd n = VectorXd::Constant(nNodes, 1.0e16);
            VectorXd p = VectorXd::Constant(nNodes, 1.0e16);

            assembler.assembleElectronContinuity(psi, n, p);
            REQUIRE(vectorIsFinite(assembler.rhs()));
            REQUIRE(assembler.matrix().rows() == nNodes);
            REQUIRE(assembler.matrix().cols() == nNodes);

            const std::vector<Index> oxideGateOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideGateOnlyNodes.empty());
            for (Index node : oxideGateOnlyNodes) {
                INFO("electron oxide node " << node);
                const int row = static_cast<int>(node);
                REQUIRE(assembler.matrix().coeff(row, row) == 1.0);
                REQUIRE(assembler.rhs()(row) == 0.0);
            }

            assembler.assembleHoleContinuity(psi, n, p);
            REQUIRE(vectorIsFinite(assembler.rhs()));
            for (Index node : oxideGateOnlyNodes) {
                INFO("hole oxide node " << node);
                const int row = static_cast<int>(node);
                REQUIRE(assembler.matrix().coeff(row, row) == 1.0);
                REQUIRE(assembler.rhs()(row) == 0.0);
            }
        }
    }
}

TEST_CASE("mixed Si/SiO2 MOS coupled DD residual and Jacobian are finite",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

            const int nNodes = static_cast<int>(mesh.numNodes());
            CoupledDDState state;
            state.psi = VectorXd::Zero(nNodes);
            state.phin = VectorXd::Constant(nNodes, 0.25);
            state.phip = VectorXd::Constant(nNodes, -0.25);
            const VectorXd x = assembler.pack(state);

            CoupledDDBoundaryConditions bcs;
            const VectorXd residual = assembler.residual(x, bcs);
            REQUIRE(residual.size() == 3 * nNodes);
            REQUIRE(vectorIsFinite(residual));

            const SparseMatrixd jacobian = assembler.assembleJacobian(x, bcs);
            REQUIRE(jacobian.rows() == 3 * nNodes);
            REQUIRE(jacobian.cols() == 3 * nNodes);
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator it(jacobian, outer); it; ++it)
                    REQUIRE(std::isfinite(it.value()));
            }

            const std::vector<Index> oxideGateOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideGateOnlyNodes.empty());
            for (Index node : oxideGateOnlyNodes) {
                INFO("coupled oxide node " << node);
                const int electronRow = nNodes + static_cast<int>(node);
                const int holeRow = 2 * nNodes + static_cast<int>(node);
                REQUIRE(residual(electronRow) == state.phin(static_cast<int>(node)));
                REQUIRE(residual(holeRow) == state.phip(static_cast<int>(node)));
                REQUIRE(jacobian.coeff(electronRow, electronRow) == 1.0);
                REQUIRE(jacobian.coeff(holeRow, holeRow) == 1.0);
            }
        }
    }
}
