#include "vela/equation/AssemblerUtils.h"
#include "vela/io/MeshReader.h"
#include "vela/simulation/ConfigParsing.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>

int main(int argc, char** argv)
{
    using namespace vela;
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: vela_tri3_box_audit mesh.json output.csv [simulation.json]\n";
        return 2;
    }
    try {
        const nlohmann::json defaults = {{"format_version", 2}, {"mesh_geometry", {
            {"poisson_permittivity_policy", "cell_material"},
            {"cell_box_policy", "delaunay_transfer"}}}};
        nlohmann::json config = defaults;
        if (argc == 4) {
            std::ifstream input(argv[3]);
            input >> config;
        }
        const auto scaling = parseUnitScalingConfig(config);
        JsonMeshReader reader;
        auto mesh = reader.read(argv[1], scaling);
        mesh.buildBoxGeometry(parseBoxGeometryOptions(config));
        const Real areaToUm2 = scaling.unitSystem().areaM2PerInternal() * 1e12;
        std::vector<std::string> nodeContact(mesh.numNodes());
        for (const auto& contact : mesh.contacts())
            for (Index node : contact.node_ids)
                nodeContact.at(node) = contact.name;
        std::ofstream out(argv[2]);
        if (!out)
            throw std::runtime_error("Cannot write geometry audit");
        out << std::setprecision(17)
            << "cell,region,local_vertex,node,next_node,contact,next_contact,g_edge,area_um2,"
               "coefficient_um2,truncated_um2,coefficient_normalized_um2\n";
        for (Index id = 0; id < mesh.numCells(); ++id) {
            const auto& cell = mesh.getCell(id);
            if (cell.type != CellType::Tri3)
                throw std::invalid_argument("Tri3 box audit requires triangular cells");
            const auto raw = detail::cellBoxNodeMeasures(mesh, id);
            const auto truncated = detail::tri3ElementVertexBoxMeasures(mesh, cell);
            std::array<Real, 3> normalized;
            try {
                normalized = detail::conservativeCoefficientBoxMeasures(mesh, id);
            } catch (const std::invalid_argument&) {
                normalized.fill(std::numeric_limits<Real>::quiet_NaN());
            }
            const Real area = .5 * std::abs(detail::triangleSignedDoubleArea(
                detail::meshPoint(mesh, cell.node_ids[0]), detail::meshPoint(mesh, cell.node_ids[1]),
                detail::meshPoint(mesh, cell.node_ids[2])));
            const auto& coefficients = mesh.poissonCellEdgeCoefficients(id);
            for (std::size_t k = 0; k < 3; ++k) {
                const Index node = cell.node_ids[k];
                const Index next = cell.node_ids[(k + 1) % 3];
                out << id << ',' << cell.region_id << ',' << k << ',' << node << ',' << next << ','
                    << nodeContact[node] << ',' << nodeContact[next] << ',' << coefficients[k] << ','
                    << area * areaToUm2 << ',' << raw[k] * areaToUm2 << ',' << truncated[k] * areaToUm2 << ','
                    << normalized[k] * areaToUm2 << '\n';
            }
        }
        return out ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
