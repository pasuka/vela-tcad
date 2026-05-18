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

} // namespace

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
