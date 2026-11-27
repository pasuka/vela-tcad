#include "vela/equation/ExtendedAveragePoissonGeometry.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/core/ArithmeticStageInput.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/mesh/AverageBox.h"
#include <fstream>

namespace vela {

ExtendedPermittivities readExtendedPermittivities(
    const std::string& materialsFile, const MaterialDatabase& database)
{
    std::ifstream input(materialsFile);
    if (!input)
        throw std::invalid_argument("Cannot open extended AverageBox materials");
    const auto root = nlohmann::json::parse(input);
    ExtendedPermittivities result;
    for (const auto& entry : root.at("materials")) {
        if (!entry.contains("decimal_eps_r"))
            continue;
        const auto name = entry.at("name").get<std::string>();
        const auto value = decimalCompatibilityNumber(entry.at("decimal_eps_r"),
            name + ".eps_r", database.getMaterial(name).eps_r);
        if (!(value >= 1 && value <= 1000) || !result.emplace(name, value).second)
            throw std::invalid_argument("Invalid extended permittivity binding: " + name);
    }
    return result;
}

ExtendedAveragePoissonGeometry buildExtendedAveragePoissonGeometry(
    const DeviceMesh& mesh, const MaterialDatabase& database, long double temperature_K)
{
    return buildExtendedAveragePoissonGeometry(mesh, database, temperature_K, {});
}

ExtendedAveragePoissonGeometry buildExtendedAveragePoissonGeometry(
    const DeviceMesh& mesh, const MaterialDatabase& database, long double temperature_K,
    const ExtendedPermittivities& permittivities)
{
    if (mesh.poissonPermittivityPolicy() != BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial ||
        !(temperature_K > 0) || !std::isfinite(temperature_K))
        throw std::invalid_argument("Extended AverageBox requires cell-material Poisson geometry");
    const auto faces = buildArithmeticAverageCellBox(mesh);
    const auto adjacency = detail::buildEdgeCellMap(mesh);
    const auto materials = detail::buildCellMaterials(mesh, database, static_cast<Real>(temperature_K));
    const auto eps0 = constants::extendedVacuumPermittivity();
    const auto permittivity = [&](Index cell) {
        const auto& material = materials.at(cell);
        const auto it = permittivities.find(material.name);
        return eps0 * (it == permittivities.end() ? static_cast<long double>(material.eps_r) : it->second);
    };
    ExtendedAveragePoissonGeometry result{std::vector<long double>(mesh.numEdges()),
        std::vector<long double>(mesh.numNodes()), eps0};
    std::vector<bool> transport(mesh.numNodes());
    for (const auto& cell : mesh.cells()) {
        const auto& material = materials.at(cell.id);
        result.permittivityReference_F_per_m = std::max(result.permittivityReference_F_per_m,
                                                        permittivity(cell.id));
        if (!detail::isTransportMaterial(material))
            continue;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto i = cell.node_ids[k], j = cell.node_ids[(k + 1) % 3];
            const auto [ix, iy] = mesh.arithmeticCoordinates(i);
            const auto [jx, jy] = mesh.arithmeticCoordinates(j);
            const auto dx = jx - ix, dy = jy - iy;
            const auto measure = faces[cell.id][k] * (dx * dx + dy * dy) / 4;
            result.chargeVolumes[i] += measure;
            result.chargeVolumes[j] += measure;
            transport[i] = transport[j] = true;
        }
    }
    for (const auto& edge : mesh.edges()) {
        for (const auto cellId : adjacency.at(edge.id)) {
            const auto& cell = mesh.getCell(cellId);
            const auto epsilon = permittivity(cellId);
            if (!(epsilon > 0) || !std::isfinite(epsilon))
                throw std::invalid_argument("Invalid permittivity for extended AverageBox");
            for (std::size_t k = 0; k < 3; ++k) {
                const auto i = cell.node_ids[k], j = cell.node_ids[(k + 1) % 3];
                if ((i == edge.n0 && j == edge.n1) || (j == edge.n0 && i == edge.n1))
                    result.edgeStiffness[edge.id] += epsilon * faces[cellId][k];
            }
        }
        if (!std::isfinite(result.edgeStiffness[edge.id]))
            throw std::invalid_argument("Nonfinite extended AverageBox stiffness");
    }
    for (Index node = 0; node < mesh.numNodes(); ++node)
        if (transport[node] && (!(result.chargeVolumes[node] > 0) ||
                                !std::isfinite(result.chargeVolumes[node])))
            throw std::invalid_argument("Invalid extended AverageBox transport node volume");
    return result;
}

} // namespace vela
