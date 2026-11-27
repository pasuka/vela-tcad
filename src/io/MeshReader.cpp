#include "vela/io/MeshReader.h"
#include "vela/core/ArithmeticDecimal.h"
#include <algorithm>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

std::runtime_error meshError(const std::string& filename, const std::string& message)
{
    return std::runtime_error("Mesh file '" + filename + "': " + message);
}

template <typename EntityRange>
void validateSequentialIds(const EntityRange& entities,
                           const std::string& entity_name,
                           const std::string& filename)
{
    for (vela::Index expected = 0; expected < entities.size(); ++expected) {
        const auto actual = entities[expected].id;
        if (actual != expected) {
            throw meshError(filename,
                            entity_name + " id " + std::to_string(actual) +
                            " must be " + std::to_string(expected) +
                            " for zero-based contiguous ids");
        }
    }
}

void validateMesh(const vela::DeviceMesh& mesh, const std::string& filename)
{
    validateSequentialIds(mesh.nodes(), "node", filename);
    validateSequentialIds(mesh.cells(), "triangle", filename);
    validateSequentialIds(mesh.regions(), "region", filename);
    validateSequentialIds(mesh.contacts(), "contact", filename);

    for (const auto& cell : mesh.cells()) {
        if (cell.node_ids.size() != 3) {
            throw meshError(filename,
                            "triangle id " + std::to_string(cell.id) +
                            " must have exactly 3 node ids");
        }

        for (const auto node_id : cell.node_ids) {
            if (node_id >= mesh.numNodes()) {
                throw meshError(filename,
                                "triangle id " + std::to_string(cell.id) +
                                " references missing node id " + std::to_string(node_id));
            }
        }

        if (cell.region_id >= mesh.numRegions()) {
            throw meshError(filename,
                            "triangle id " + std::to_string(cell.id) +
                            " references missing region id " + std::to_string(cell.region_id));
        }
    }

    for (const auto& region : mesh.regions()) {
        for (const auto cell_id : region.cell_ids) {
            if (cell_id >= mesh.numCells()) {
                throw meshError(filename,
                                "region id " + std::to_string(region.id) +
                                " references missing cell id " + std::to_string(cell_id));
            }

            const auto& cell = mesh.getCell(cell_id);
            if (cell.region_id != region.id) {
                throw meshError(filename,
                                "region id " + std::to_string(region.id) +
                                " references cell id " + std::to_string(cell_id) +
                                " whose region id is " + std::to_string(cell.region_id));
            }
        }
    }

    for (const auto& contact : mesh.contacts()) {
        if (contact.region_id >= mesh.numRegions()) {
            throw meshError(filename,
                            "contact id " + std::to_string(contact.id) +
                            " references missing region id " + std::to_string(contact.region_id));
        }

        for (const auto node_id : contact.node_ids) {
            if (node_id >= mesh.numNodes()) {
                throw meshError(filename,
                                "contact id " + std::to_string(contact.id) +
                                " references missing node id " + std::to_string(node_id));
            }
        }
        for (const auto& edge : contact.edge_node_ids) {
            for (const auto node_id : edge) {
                if (node_id >= mesh.numNodes()) {
                    throw meshError(filename,
                                    "contact id " + std::to_string(contact.id) +
                                    " edge references missing node id " +
                                    std::to_string(node_id));
                }
                if (std::find(contact.node_ids.begin(), contact.node_ids.end(), node_id) ==
                    contact.node_ids.end()) {
                    throw meshError(filename,
                                    "contact id " + std::to_string(contact.id) +
                                    " edge endpoint is absent from node_ids");
                }
            }
        }
    }
}

} // namespace

namespace vela {

DeviceMesh JsonMeshReader::read(const std::string& filename)
{
    return read(filename, UnitScalingConfig{});
}

DeviceMesh JsonMeshReader::read(const std::string& filename, UnitScalingConfig scaling)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open mesh file: " + filename);

    nlohmann::json j;
    DeviceMesh mesh;

    try {
        ifs >> j;

        // ------------------------------------------------------------------
        // Nodes
        // ------------------------------------------------------------------
        for (const auto& jn : j.at("nodes")) {
            Node n;
            n.id = jn.at("id").get<Index>();
            const auto literals = jn.find("coordinate_literals");
            if (literals == jn.end()) {
                n.x = scaling.lengthToInternal(jn.at("x").get<Real>());
                n.y = scaling.lengthToInternal(jn.at("y").get<Real>());
                mesh.addNode(n);
                continue;
            }
            const auto coordinate = [&](const char* key, std::size_t component) {
                const Real compatibility = jn.at(key).get<Real>();
                if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
                    throw std::invalid_argument("coordinate_literals requires an extended80 build");
                if (!literals->is_array() || literals->size() != 2)
                    throw std::invalid_argument("coordinate_literals requires two decimal operands");
                const auto literal = literals->at(component).get<std::string>();
                if (parseArithmeticDecimal<double>(literal) != compatibility)
                    throw std::invalid_argument("coordinate_literals differs from the compatibility coordinate");
                return parseArithmeticDecimal(literal) * scaling.lengthToInternal(1.0);
            };
            const DeviceMesh::ArithmeticPoint coordinates{coordinate("x", 0), coordinate("y", 1)};
            n.x = static_cast<Real>(coordinates[0]);
            n.y = static_cast<Real>(coordinates[1]);
            mesh.addNode(n, coordinates);
        }

        // ------------------------------------------------------------------
        // Triangular cells
        // ------------------------------------------------------------------
        for (const auto& jc : j.at("triangles")) {
            Cell c;
            c.id        = jc.at("id").get<Index>();
            c.type      = CellType::Tri3;
            c.region_id = jc.at("region_id").get<Index>();
            c.node_ids  = jc.at("node_ids").get<std::vector<Index>>();
            mesh.addCell(c);
        }

        // ------------------------------------------------------------------
        // Regions
        // ------------------------------------------------------------------
        for (const auto& jr : j.at("regions")) {
            Region r;
            r.id       = jr.at("id").get<Index>();
            r.name     = jr.at("name").get<std::string>();
            r.material = jr.at("material").get<std::string>();
            r.cell_ids = jr.at("cell_ids").get<std::vector<Index>>();
            mesh.addRegion(r);
        }

        // ------------------------------------------------------------------
        // Contacts
        // ------------------------------------------------------------------
        for (const auto& jct : j.at("contacts")) {
            Contact ct;
            ct.id        = jct.at("id").get<Index>();
            ct.name      = jct.at("name").get<std::string>();
            ct.region_id = jct.at("region_id").get<Index>();
            ct.node_ids  = jct.at("node_ids").get<std::vector<Index>>();
            if (jct.contains("edge_node_ids")) {
                ct.edge_node_ids =
                    jct.at("edge_node_ids").get<std::vector<std::array<Index, 2>>>();
            }
            mesh.addContact(ct);
        }
    } catch (const nlohmann::json::exception& e) {
        throw meshError(filename, std::string("invalid JSON mesh: ") + e.what());
    }

    validateMesh(mesh, filename);

    // Build edges from the cells just loaded
    mesh.buildEdges();

    return mesh;
}

} // namespace vela
