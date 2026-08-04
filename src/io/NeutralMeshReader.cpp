#include "vela/io/NeutralMeshReader.h"

#include "vela/io/CsvUtils.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {
namespace {

std::unordered_map<std::string, std::size_t> headerIndex(const std::vector<std::string>& header)
{
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < header.size(); ++i) {
        index[header[i]] = i;
    }
    return index;
}

std::size_t requireColumn(const std::unordered_map<std::string, std::size_t>& index,
                          const std::string& name,
                          const std::string& fileLabel)
{
    const auto it = index.find(name);
    if (it == index.end()) {
        throw std::runtime_error("NeutralMeshReader: missing column '" + name + "' in " + fileLabel + ".");
    }
    return it->second;
}

} // namespace

DeviceMesh NeutralMeshReader::readDirectory(const std::filesystem::path& directory,
                                            UnitScalingConfig scaling) const
{
    DeviceMesh mesh;
    std::vector<Region> parsedRegions;

    const std::filesystem::path nodesPath = directory / "nodes.csv";
    const std::filesystem::path elementsPath = directory / "elements.csv";
    const std::filesystem::path regionsPath = directory / "regions.csv";
    const std::filesystem::path contactsPath = directory / "contacts.csv";

    {
        std::ifstream in(nodesPath);
        if (!in.is_open()) {
            throw std::runtime_error("NeutralMeshReader: cannot open " + nodesPath.string());
        }
        std::string line;
        if (!std::getline(in, line)) {
            throw std::runtime_error("NeutralMeshReader: nodes.csv is empty.");
        }
        const auto header = splitCsvLine(line);
        const auto idx = headerIndex(header);
        const std::size_t idCol = requireColumn(idx, "id", "nodes.csv");
        const std::size_t xCol = requireColumn(idx, "x_um", "nodes.csv");
        const std::size_t yCol = requireColumn(idx, "y_um", "nodes.csv");

        while (std::getline(in, line)) {
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            Node node;
            node.id = static_cast<Index>(std::stoul(row.at(idCol)));
            node.x = scaling.lengthToInternal(std::stod(row.at(xCol)));
            node.y = scaling.lengthToInternal(std::stod(row.at(yCol)));
            mesh.addNode(node);
        }
    }

    {
        std::ifstream in(regionsPath);
        if (!in.is_open()) {
            throw std::runtime_error("NeutralMeshReader: cannot open " + regionsPath.string());
        }
        std::string line;
        if (!std::getline(in, line)) {
            throw std::runtime_error("NeutralMeshReader: regions.csv is empty.");
        }
        const auto header = splitCsvLine(line);
        const auto idx = headerIndex(header);
        const std::size_t idCol = requireColumn(idx, "region_id", "regions.csv");
        const std::size_t nameCol = requireColumn(idx, "region_name", "regions.csv");
        const std::size_t materialCol = requireColumn(idx, "material", "regions.csv");

        while (std::getline(in, line)) {
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            Region region;
            region.id = static_cast<Index>(std::stoul(row.at(idCol)));
            region.name = row.at(nameCol);
            region.material = row.at(materialCol);
            parsedRegions.push_back(region);
        }
    }

    {
        std::ifstream in(elementsPath);
        if (!in.is_open()) {
            throw std::runtime_error("NeutralMeshReader: cannot open " + elementsPath.string());
        }
        std::string line;
        if (!std::getline(in, line)) {
            throw std::runtime_error("NeutralMeshReader: elements.csv is empty.");
        }
        const auto header = splitCsvLine(line);
        const auto idx = headerIndex(header);
        const std::size_t idCol = requireColumn(idx, "id", "elements.csv");
        const std::size_t n0Col = requireColumn(idx, "node0", "elements.csv");
        const std::size_t n1Col = requireColumn(idx, "node1", "elements.csv");
        const std::size_t n2Col = requireColumn(idx, "node2", "elements.csv");
        const std::size_t ridCol = requireColumn(idx, "region_id", "elements.csv");

        while (std::getline(in, line)) {
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            Cell cell;
            cell.id = static_cast<Index>(std::stoul(row.at(idCol)));
            cell.type = CellType::Tri3;
            cell.region_id = static_cast<Index>(std::stoul(row.at(ridCol)));
            cell.node_ids = {
                static_cast<Index>(std::stoul(row.at(n0Col))),
                static_cast<Index>(std::stoul(row.at(n1Col))),
                static_cast<Index>(std::stoul(row.at(n2Col))),
            };
            mesh.addCell(cell);
        }
    }

    for (Index c = 0; c < mesh.numCells(); ++c) {
        const Cell& cell = mesh.getCell(c);
        if (cell.region_id >= parsedRegions.size()) {
            throw std::runtime_error("NeutralMeshReader: cell region id out of range.");
        }
        parsedRegions[cell.region_id].cell_ids.push_back(cell.id);
    }

    for (const auto& region : parsedRegions) {
        mesh.addRegion(region);
    }

    if (std::filesystem::exists(contactsPath)) {
        std::ifstream in(contactsPath);
        if (!in.is_open()) {
            throw std::runtime_error("NeutralMeshReader: cannot open " + contactsPath.string());
        }
        std::string line;
        if (std::getline(in, line)) {
            const auto header = splitCsvLine(line);
            const auto idx = headerIndex(header);
            const std::size_t idCol = requireColumn(idx, "id", "contacts.csv");
            const std::size_t nameCol = requireColumn(idx, "name", "contacts.csv");
            const std::size_t regionCol = requireColumn(idx, "region_id", "contacts.csv");
            const std::size_t nodesCol = requireColumn(idx, "node_ids", "contacts.csv");
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                const auto row = splitCsvLine(line);
                Contact contact;
                contact.id = static_cast<Index>(std::stoul(row.at(idCol)));
                contact.name = row.at(nameCol);
                contact.region_id = static_cast<Index>(std::stoul(row.at(regionCol)));
                std::string nodes = row.at(nodesCol);
                std::size_t start = 0;
                while (start < nodes.size()) {
                    const std::size_t sep = nodes.find(';', start);
                    const std::string token = sep == std::string::npos
                        ? nodes.substr(start)
                        : nodes.substr(start, sep - start);
                    if (!token.empty()) {
                        contact.node_ids.push_back(static_cast<Index>(std::stoul(token)));
                    }
                    if (sep == std::string::npos) break;
                    start = sep + 1;
                }
                mesh.addContact(contact);
            }
        }
    }

    mesh.buildEdges();
    return mesh;
}

DopingModel NeutralMeshReader::readDopingCsv(const std::filesystem::path& path,
                                             Index nodeCount,
                                             UnitScalingConfig scaling) const
{
    DopingModel model(nodeCount);
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("NeutralMeshReader: cannot open " + path.string());
    }
    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("NeutralMeshReader: doping csv is empty.");
    }
    const auto header = splitCsvLine(line);
    const auto idx = headerIndex(header);
    const std::size_t nodeCol = requireColumn(idx, "node_id", "doping.csv");
    const std::size_t donorCol = requireColumn(idx, "donors_cm3", "doping.csv");
    const std::size_t acceptorCol = requireColumn(idx, "acceptors_cm3", "doping.csv");
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto row = splitCsvLine(line);
        const Index nodeId = static_cast<Index>(std::stoul(row.at(nodeCol)));
        const Real donors = scaling.concentrationToInternal(std::stod(row.at(donorCol)));
        const Real acceptors = scaling.concentrationToInternal(std::stod(row.at(acceptorCol)));
        model.setNodeDoping(nodeId, donors, acceptors);
    }
    return model;
}

} // namespace vela
