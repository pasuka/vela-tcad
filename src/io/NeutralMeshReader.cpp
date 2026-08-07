#include "vela/io/NeutralMeshReader.h"

#include "vela/io/CsvUtils.h"

#include <fstream>
#include <limits>
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

std::runtime_error rowError(const std::string& fileLabel,
                            std::size_t rowNumber,
                            const std::string& fieldLabel,
                            const std::string& message)
{
    return std::runtime_error("NeutralMeshReader: " + message + " for " + fieldLabel +
                              " in " + fileLabel + " at row " + std::to_string(rowNumber) + ".");
}

Index parseStrictIndex(const std::string& value,
                       const std::string& fileLabel,
                       const std::string& fieldLabel,
                       std::size_t rowNumber)
{
    if (value.empty()) {
        throw rowError(fileLabel, rowNumber, fieldLabel, "invalid");
    }
    for (char c : value) {
        if (c < '0' || c > '9') {
            throw rowError(fileLabel, rowNumber, fieldLabel, "invalid");
        }
    }
    Index parsed = 0;
    for (char c : value) {
        if (parsed > (std::numeric_limits<Index>::max() - static_cast<Index>(c - '0')) / 10) {
            throw rowError(fileLabel, rowNumber, fieldLabel, "invalid");
        }
        parsed = parsed * 10 + static_cast<Index>(c - '0');
    }
    return parsed;
}

Real parseStrictReal(const std::string& value,
                     const std::string& fileLabel,
                     const std::string& fieldLabel,
                     std::size_t rowNumber)
{
    try {
        std::size_t pos = 0;
        const Real parsed = std::stod(value, &pos);
        if (pos != value.size() || !std::isfinite(parsed)) {
            throw rowError(fileLabel, rowNumber, fieldLabel, "invalid");
        }
        return parsed;
    } catch (const std::exception&) {
        throw rowError(fileLabel, rowNumber, fieldLabel, "invalid");
    }
}

} // namespace

DeviceMesh NeutralMeshReader::readDirectory(const std::filesystem::path& directory,
                                            UnitScalingConfig scaling) const
{
    DeviceMesh mesh;
    std::vector<Region> parsedRegions;
    std::vector<std::size_t> nodeIds;
    std::vector<std::size_t> regionIds;
    std::vector<std::size_t> cellIds;
    std::vector<std::size_t> contactIds;

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

        std::size_t rowNumber = 1;
        while (std::getline(in, line)) {
            ++rowNumber;
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            if (row.size() <= std::max({idCol, xCol, yCol})) {
                throw rowError("nodes.csv", rowNumber, "row", "malformed");
            }
            Node node;
            node.id = parseStrictIndex(row.at(idCol), "nodes.csv", "node id", rowNumber);
            if (node.id != mesh.numNodes()) {
                throw rowError("nodes.csv", rowNumber, "node id", "out-of-order");
            }
            node.x = scaling.lengthToInternal(parseStrictReal(row.at(xCol), "nodes.csv", "x_um", rowNumber));
            node.y = scaling.lengthToInternal(parseStrictReal(row.at(yCol), "nodes.csv", "y_um", rowNumber));
            if (!std::isfinite(node.x) || !std::isfinite(node.y)) {
                throw rowError("nodes.csv", rowNumber, "coordinates", "invalid");
            }
            mesh.addNode(node);
            nodeIds.push_back(node.id);
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

        std::size_t rowNumber = 1;
        while (std::getline(in, line)) {
            ++rowNumber;
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            if (row.size() <= std::max({idCol, nameCol, materialCol})) {
                throw rowError("regions.csv", rowNumber, "row", "malformed");
            }
            Region region;
            region.id = parseStrictIndex(row.at(idCol), "regions.csv", "region id", rowNumber);
            if (region.id != parsedRegions.size()) {
                throw rowError("regions.csv", rowNumber, "region id", "out-of-order");
            }
            region.name = row.at(nameCol);
            region.material = row.at(materialCol);
            parsedRegions.push_back(region);
            regionIds.push_back(region.id);
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

        std::size_t rowNumber = 1;
        while (std::getline(in, line)) {
            ++rowNumber;
            if (line.empty()) continue;
            const auto row = splitCsvLine(line);
            if (row.size() <= std::max({idCol, n0Col, n1Col, n2Col, ridCol})) {
                throw rowError("elements.csv", rowNumber, "row", "malformed");
            }
            Cell cell;
            cell.id = parseStrictIndex(row.at(idCol), "elements.csv", "cell id", rowNumber);
            if (cell.id != mesh.numCells()) {
                throw rowError("elements.csv", rowNumber, "cell id", "out-of-order");
            }
            cell.type = CellType::Tri3;
            cell.region_id = parseStrictIndex(row.at(ridCol), "elements.csv", "region id", rowNumber);
            cell.node_ids = {
                parseStrictIndex(row.at(n0Col), "elements.csv", "node0", rowNumber),
                parseStrictIndex(row.at(n1Col), "elements.csv", "node1", rowNumber),
                parseStrictIndex(row.at(n2Col), "elements.csv", "node2", rowNumber),
            };
            mesh.addCell(cell);
            cellIds.push_back(cell.id);
        }
    }

    for (Index c = 0; c < mesh.numCells(); ++c) {
        const Cell& cell = mesh.getCell(c);
        if (cell.region_id >= parsedRegions.size()) {
            throw std::runtime_error("NeutralMeshReader: cell region id out of range.");
        }
        for (const auto nodeId : cell.node_ids) {
            if (nodeId >= mesh.numNodes()) {
                throw std::runtime_error("NeutralMeshReader: cell node id out of range.");
            }
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
            std::size_t rowNumber = 1;
            while (std::getline(in, line)) {
                ++rowNumber;
                if (line.empty()) continue;
                const auto row = splitCsvLine(line);
                if (row.size() <= std::max({idCol, nameCol, regionCol, nodesCol})) {
                    throw rowError("contacts.csv", rowNumber, "row", "malformed");
                }
                Contact contact;
                contact.id = parseStrictIndex(row.at(idCol), "contacts.csv", "contact id", rowNumber);
                if (contact.id != mesh.numContacts()) {
                    throw rowError("contacts.csv", rowNumber, "contact id", "out-of-order");
                }
                contact.name = row.at(nameCol);
                contact.region_id = parseStrictIndex(row.at(regionCol), "contacts.csv", "region id", rowNumber);
                std::string nodes = row.at(nodesCol);
                std::size_t start = 0;
                while (start < nodes.size()) {
                    const std::size_t sep = nodes.find(';', start);
                    const std::string token = sep == std::string::npos
                        ? nodes.substr(start)
                        : nodes.substr(start, sep - start);
                    if (!token.empty()) {
                        contact.node_ids.push_back(parseStrictIndex(token, "contacts.csv", "node id", rowNumber));
                    }
                    if (sep == std::string::npos) break;
                    start = sep + 1;
                }
                mesh.addContact(contact);
                contactIds.push_back(contact.id);
            }
        }
    }

    for (const auto& contact : mesh.contacts()) {
        if (contact.region_id >= mesh.numRegions()) {
            throw std::runtime_error("NeutralMeshReader: contact region id out of range.");
        }
        for (const auto nodeId : contact.node_ids) {
            if (nodeId >= mesh.numNodes()) {
                throw std::runtime_error("NeutralMeshReader: contact node id out of range.");
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
    return readNamedDopingCsv(path, nodeCount, scaling, "doping.csv");
}

DopingModel NeutralMeshReader::readNamedDopingCsv(const std::filesystem::path& path,
                                                  Index nodeCount,
                                                  UnitScalingConfig scaling,
                                                  const std::string& fileLabel) const
{
    DopingModel model(nodeCount);
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("NeutralMeshReader: cannot open " + path.string());
    }
    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("NeutralMeshReader: " + fileLabel + " is empty.");
    }
    const auto header = splitCsvLine(line);
    const auto idx = headerIndex(header);
    const std::size_t nodeCol = requireColumn(idx, "node_id", fileLabel);
    const std::size_t donorCol = requireColumn(idx, "donors_cm3", fileLabel);
    const std::size_t acceptorCol = requireColumn(idx, "acceptors_cm3", fileLabel);
    std::vector<bool> seen(static_cast<std::size_t>(nodeCount), false);
    std::size_t rowNumber = 1;
    while (std::getline(in, line)) {
        ++rowNumber;
        if (line.empty()) continue;
        const auto row = splitCsvLine(line);
        if (row.size() <= std::max({nodeCol, donorCol, acceptorCol})) {
            throw rowError(fileLabel, rowNumber, "row", "malformed");
        }
        const Index nodeId = parseStrictIndex(row.at(nodeCol), fileLabel, "node id", rowNumber);
        const Real donors = scaling.concentrationToInternal(parseStrictReal(row.at(donorCol), fileLabel, "donors_cm3", rowNumber));
        const Real acceptors = scaling.concentrationToInternal(parseStrictReal(row.at(acceptorCol), fileLabel, "acceptors_cm3", rowNumber));
        if (!std::isfinite(donors) || !std::isfinite(acceptors)) {
            throw rowError(fileLabel, rowNumber, "doping", "invalid");
        }
        if (donors < 0.0 || acceptors < 0.0) {
            throw rowError(fileLabel, rowNumber, "doping", "invalid");
        }
        if (nodeId >= static_cast<Index>(nodeCount)) {
            throw rowError(fileLabel, rowNumber, "node id", "invalid");
        }
        if (seen[static_cast<std::size_t>(nodeId)]) {
            throw rowError(fileLabel, rowNumber, "node id", "duplicate");
        }
        seen[static_cast<std::size_t>(nodeId)] = true;
        model.setNodeDoping(nodeId, donors, acceptors);
    }
    for (Index nodeId = 0; nodeId < nodeCount; ++nodeId) {
        if (!seen[static_cast<std::size_t>(nodeId)]) {
            throw std::runtime_error("NeutralMeshReader: missing node id in " + fileLabel + ".");
        }
    }
    return model;
}

} // namespace vela
