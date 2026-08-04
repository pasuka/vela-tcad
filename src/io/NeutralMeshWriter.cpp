#include "vela/io/NeutralMeshWriter.h"

#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace vela {

void NeutralMeshWriter::write(const MeshBundle2D& meshBundle,
                              const std::filesystem::path& outputDirectory,
                              const DopingModel* doping)
{
    std::filesystem::create_directories(outputDirectory);

    {
        std::ofstream out(outputDirectory / "nodes.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open nodes.csv.");
        }
        out.precision(std::numeric_limits<double>::max_digits10);
        out << "id,x_um,y_um\n";
        for (const auto& node : meshBundle.mesh.nodes()) {
            out << node.id << "," << node.x << "," << node.y << "\n";
        }
    }

    {
        std::ofstream out(outputDirectory / "elements.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open elements.csv.");
        }
        out << "id,node0,node1,node2,region_id,region_name,material\n";
        for (const auto& cell : meshBundle.mesh.cells()) {
            const auto& region = meshBundle.mesh.getRegion(cell.region_id);
            out << cell.id << "," << cell.node_ids[0] << "," << cell.node_ids[1] << ","
                << cell.node_ids[2] << "," << region.id << "," << region.name << ","
                << region.material << "\n";
        }
    }

    {
        std::ofstream out(outputDirectory / "regions.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open regions.csv.");
        }
        out << "region_id,region_name,material\n";
        for (const auto& region : meshBundle.mesh.regions()) {
            out << region.id << "," << region.name << "," << region.material << "\n";
        }
    }

    {
        std::ofstream out(outputDirectory / "boundaries.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open boundaries.csv.");
        }
        out << "id,node0,node1,contact_name,owner_region\n";
        for (std::size_t i = 0; i < meshBundle.boundaryEdges.size(); ++i) {
            const auto& edge = meshBundle.boundaryEdges[i];
            out << i << "," << edge.node0 << "," << edge.node1 << ","
                << edge.contactName << "," << edge.ownerRegion << "\n";
        }
    }

    {
        std::ofstream out(outputDirectory / "contacts.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open contacts.csv.");
        }
        out << "id,name,region_id,node_ids\n";
        for (const auto& contact : meshBundle.mesh.contacts()) {
            out << contact.id << "," << contact.name << "," << contact.region_id << ",";
            for (std::size_t i = 0; i < contact.node_ids.size(); ++i) {
                if (i > 0) {
                    out << ";";
                }
                out << contact.node_ids[i];
            }
            out << "\n";
        }
    }

    if (doping != nullptr) {
        if (doping->numNodes() != meshBundle.mesh.numNodes()) {
            throw std::invalid_argument("NeutralMeshWriter: doping node count does not match mesh.");
        }
        std::ofstream out(outputDirectory / "doping.csv");
        if (!out.is_open()) {
            throw std::runtime_error("NeutralMeshWriter: cannot open doping.csv.");
        }
        out.precision(std::numeric_limits<double>::max_digits10);
        out << "node_id,donors_cm3,acceptors_cm3\n";
        for (Index i = 0; i < meshBundle.mesh.numNodes(); ++i) {
            out << i << "," << doping->donors(i) << "," << doping->acceptors(i) << "\n";
        }
    }
}

} // namespace vela

