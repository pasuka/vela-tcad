#include "vela/io/MeshReader.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

namespace vela {

DeviceMesh JsonMeshReader::read(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open mesh file: " + filename);

    nlohmann::json j;
    ifs >> j;

    DeviceMesh mesh;

    // ------------------------------------------------------------------
    // Nodes
    // ------------------------------------------------------------------
    for (const auto& jn : j.at("nodes")) {
        Node n;
        n.id = jn.at("id").get<Index>();
        n.x  = jn.at("x").get<Real>();
        n.y  = jn.at("y").get<Real>();
        mesh.addNode(n);
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
        mesh.addContact(ct);
    }

    // Build edges from the cells just loaded
    mesh.buildEdges();

    return mesh;
}

} // namespace vela
