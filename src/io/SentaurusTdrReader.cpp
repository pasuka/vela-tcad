#include "vela/io/SentaurusTdrReader.h"

#include <hdf5.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace vela {
namespace {

struct H5Object {
    hid_t id = -1;
    herr_t (*close)(hid_t) = nullptr;

    H5Object() = default;
    H5Object(hid_t objectId, herr_t (*closeFn)(hid_t)) : id(objectId), close(closeFn) {}
    H5Object(const H5Object&) = delete;
    H5Object& operator=(const H5Object&) = delete;
    H5Object(H5Object&& other) noexcept : id(other.id), close(other.close)
    {
        other.id = -1;
        other.close = nullptr;
    }
    H5Object& operator=(H5Object&& other) noexcept
    {
        if (this != &other) {
            reset();
            id = other.id;
            close = other.close;
            other.id = -1;
            other.close = nullptr;
        }
        return *this;
    }
    ~H5Object() { reset(); }

    void reset()
    {
        if (id >= 0 && close != nullptr) {
            close(id);
        }
        id = -1;
        close = nullptr;
    }

    explicit operator bool() const { return id >= 0; }
};

struct VertexRecord {
    double x;
    double y;
};

bool exists(hid_t parent, const std::string& path)
{
    return H5Lexists(parent, path.c_str(), H5P_DEFAULT) > 0;
}

std::string readStringAttribute(hid_t object, const std::string& name, const std::string& fallback = {})
{
    if (H5Aexists(object, name.c_str()) <= 0) {
        return fallback;
    }
    H5Object attr(H5Aopen(object, name.c_str(), H5P_DEFAULT), H5Aclose);
    if (!attr) {
        return fallback;
    }
    H5Object type(H5Aget_type(attr.id), H5Tclose);
    if (!type) {
        return fallback;
    }

    if (H5Tis_variable_str(type.id) > 0) {
        char* value = nullptr;
        if (H5Aread(attr.id, type.id, &value) < 0 || value == nullptr) {
            return fallback;
        }
        std::string result(value);
        H5free_memory(value);
        return result;
    }

    const std::size_t size = H5Tget_size(type.id);
    std::string buffer(size, '\0');
    if (H5Aread(attr.id, type.id, buffer.data()) < 0) {
        return fallback;
    }
    const auto end = std::find(buffer.begin(), buffer.end(), '\0');
    buffer.erase(end, buffer.end());
    return buffer;
}

int readIntAttribute(hid_t object, const std::string& name, int fallback = 0)
{
    if (H5Aexists(object, name.c_str()) <= 0) {
        return fallback;
    }
    H5Object attr(H5Aopen(object, name.c_str(), H5P_DEFAULT), H5Aclose);
    int value = fallback;
    H5Aread(attr.id, H5T_NATIVE_INT, &value);
    return value;
}

std::size_t readSizeAttribute(hid_t object, const std::string& name, std::size_t fallback = 0)
{
    if (H5Aexists(object, name.c_str()) <= 0) {
        return fallback;
    }
    H5Object attr(H5Aopen(object, name.c_str(), H5P_DEFAULT), H5Aclose);
    unsigned long long value = static_cast<unsigned long long>(fallback);
    H5Aread(attr.id, H5T_NATIVE_ULLONG, &value);
    return static_cast<std::size_t>(value);
}

std::vector<std::string> childNames(hid_t group)
{
    std::vector<std::string> names;
    auto callback = [](hid_t, const char* name, const H5L_info2_t*, void* data) -> herr_t {
        auto* out = static_cast<std::vector<std::string>*>(data);
        out->emplace_back(name);
        return 0;
    };
    hsize_t index = 0;
    H5Literate2(group, H5_INDEX_NAME, H5_ITER_INC, &index, callback, &names);
    return names;
}

bool startsWith(std::string_view text, std::string_view prefix)
{
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

int trailingIndex(const std::string& name, const std::string& prefix)
{
    if (!startsWith(name, prefix)) {
        return -1;
    }
    return std::stoi(name.substr(prefix.size()));
}

SentaurusTdrRegionType regionTypeFromInt(int value)
{
    if (value == 0) {
        return SentaurusTdrRegionType::Material;
    }
    if (value == 1) {
        return SentaurusTdrRegionType::Contact;
    }
    if (value == 2) {
        return SentaurusTdrRegionType::Interface;
    }
    return SentaurusTdrRegionType::Other;
}

std::vector<int> readIntVector(hid_t group, const std::string& name)
{
    H5Object dataset(H5Dopen2(group, name.c_str(), H5P_DEFAULT), H5Dclose);
    if (!dataset) {
        return {};
    }
    H5Object space(H5Dget_space(dataset.id), H5Sclose);
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space.id, dims, nullptr);
    std::vector<int> values(static_cast<std::size_t>(dims[0]));
    if (!values.empty()) {
        H5Dread(dataset.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    }
    return values;
}

std::vector<double> readDoubleVector(hid_t groupOrDataset, const std::string& name)
{
    H5Object dataset(H5Dopen2(groupOrDataset, name.c_str(), H5P_DEFAULT), H5Dclose);
    if (!dataset) {
        return {};
    }
    H5Object space(H5Dget_space(dataset.id), H5Sclose);
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space.id, dims, nullptr);
    std::vector<double> values(static_cast<std::size_t>(dims[0]));
    if (!values.empty()) {
        H5Dread(dataset.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    }
    return values;
}

std::vector<SentaurusTdrVertex> readVertices(hid_t geometry)
{
    H5Object dataset(H5Dopen2(geometry, "vertex", H5P_DEFAULT), H5Dclose);
    if (!dataset) {
        throw std::runtime_error("TDR geometry is missing vertex dataset");
    }
    H5Object space(H5Dget_space(dataset.id), H5Sclose);
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space.id, dims, nullptr);

    H5Object type(H5Tcreate(H5T_COMPOUND, sizeof(VertexRecord)), H5Tclose);
    H5Tinsert(type.id, "x", HOFFSET(VertexRecord, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(type.id, "y", HOFFSET(VertexRecord, y), H5T_NATIVE_DOUBLE);

    std::vector<VertexRecord> raw(static_cast<std::size_t>(dims[0]));
    if (!raw.empty()) {
        H5Dread(dataset.id, type.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
    }

    std::vector<SentaurusTdrVertex> vertices;
    vertices.reserve(raw.size());
    for (const auto& item : raw) {
        vertices.push_back({item.x, item.y});
    }
    return vertices;
}

void parseElements(SentaurusTdrRegion& region, const std::vector<int>& raw)
{
    std::size_t i = 0;
    while (i < raw.size()) {
        const int kind = raw[i++];
        if (kind == 2 && i + 2 < raw.size()) {
            region.triangles.push_back({
                static_cast<std::size_t>(raw[i]),
                static_cast<std::size_t>(raw[i + 1]),
                static_cast<std::size_t>(raw[i + 2]),
            });
            i += 3;
        } else if (kind == 1 && i + 1 < raw.size()) {
            region.edges.push_back({
                static_cast<std::size_t>(raw[i]),
                static_cast<std::size_t>(raw[i + 1]),
            });
            i += 2;
        } else {
            throw std::runtime_error("unsupported Sentaurus element encoding");
        }
    }
}

std::vector<SentaurusTdrRegion> readRegions(hid_t geometry)
{
    std::vector<SentaurusTdrRegion> regions;
    for (const auto& name : childNames(geometry)) {
        const int index = trailingIndex(name, "region_");
        if (index < 0) {
            continue;
        }
        H5Object group(H5Gopen2(geometry, name.c_str(), H5P_DEFAULT), H5Gclose);
        SentaurusTdrRegion region;
        region.index = index;
        region.name = readStringAttribute(group.id, "name", name);
        region.material = readStringAttribute(group.id, "material", "");
        region.type = regionTypeFromInt(readIntAttribute(group.id, "type", 99));
        parseElements(region, readIntVector(group.id, "elements_0"));
        regions.push_back(std::move(region));
    }
    std::sort(regions.begin(), regions.end(), [](const auto& a, const auto& b) {
        return a.index < b.index;
    });
    return regions;
}

std::string readUnitAttribute(hid_t object)
{
    if (H5Aexists(object, "unit") > 0) {
        return readStringAttribute(object, "unit");
    }
    if (H5Aexists(object, "unit:name") > 0) {
        return readStringAttribute(object, "unit:name");
    }
    return {};
}

SentaurusTdrField readFieldGroup(hid_t group, const std::string& groupName)
{
    H5Object values(H5Dopen2(group, "values", H5P_DEFAULT), H5Dclose);
    if (!values) {
        throw std::runtime_error("state dataset is missing values");
    }

    SentaurusTdrField field;
    field.index = trailingIndex(groupName, "dataset_");
    field.name = readStringAttribute(values.id, "name", readStringAttribute(group, "name", groupName));
    field.region_index = readIntAttribute(values.id, "region", readIntAttribute(group, "region", -1));
    field.value_count = readSizeAttribute(values.id, "number of values",
                                          readSizeAttribute(group, "number of values", 0));
    field.unit = readUnitAttribute(values.id);
    if (field.unit.empty()) {
        field.unit = readUnitAttribute(group);
    }
    field.values = readDoubleVector(group, "values");
    if (field.value_count == 0) {
        field.value_count = field.values.size();
    }
    if (field.value_count > 0 && field.values.size() % field.value_count == 0) {
        field.component_count = field.values.size() / field.value_count;
    }
    return field;
}

std::vector<SentaurusTdrField> readFields(hid_t geometry)
{
    if (!exists(geometry, "state_0")) {
        return {};
    }
    H5Object state(H5Gopen2(geometry, "state_0", H5P_DEFAULT), H5Gclose);
    std::vector<SentaurusTdrField> fields;
    for (const auto& name : childNames(state.id)) {
        const int index = trailingIndex(name, "dataset_");
        if (index < 0) {
            continue;
        }
        H5Object group(H5Gopen2(state.id, name.c_str(), H5P_DEFAULT), H5Gclose);
        fields.push_back(readFieldGroup(group.id, name));
    }
    std::sort(fields.begin(), fields.end(), [](const auto& a, const auto& b) {
        return a.index < b.index;
    });
    return fields;
}

std::vector<std::size_t> regionNodeOrder(const SentaurusTdrRegion& region)
{
    std::vector<std::size_t> order;
    std::set<std::size_t> seen;
    auto add = [&](std::size_t node) {
        if (seen.insert(node).second) {
            order.push_back(node);
        }
    };
    for (const auto& tri : region.triangles) {
        add(tri[0]);
        add(tri[1]);
        add(tri[2]);
    }
    for (const auto& edge : region.edges) {
        add(edge[0]);
        add(edge[1]);
    }
    return order;
}

std::string materialName(const std::string& material)
{
    if (material == "Silicon") {
        return "Si";
    }
    if (material == "Oxide") {
        return "SiO2";
    }
    return material;
}

std::string sanitizeFilename(std::string value)
{
    for (char& ch : value) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '-') {
            ch = '_';
        }
    }
    return value.empty() ? "field" : value;
}

const SentaurusTdrRegion* findRegion(const SentaurusTdrInventory& inventory, int index)
{
    for (const auto& region : inventory.regions) {
        if (region.index == index) {
            return &region;
        }
    }
    return nullptr;
}

const SentaurusTdrRegion* bestMaterialRegionForContact(const SentaurusTdrInventory& inventory,
                                                       const SentaurusTdrRegion& contact)
{
    std::set<std::size_t> contactNodes;
    for (const auto& edge : contact.edges) {
        contactNodes.insert(edge[0]);
        contactNodes.insert(edge[1]);
    }

    const SentaurusTdrRegion* best = nullptr;
    std::size_t bestScore = 0;
    for (const auto& region : inventory.regions) {
        if (region.type != SentaurusTdrRegionType::Material) {
            continue;
        }
        const auto nodes = regionNodeOrder(region);
        std::size_t score = 0;
        for (const auto node : nodes) {
            if (contactNodes.contains(node)) {
                ++score;
            }
        }
        if (best == nullptr || score > bestScore) {
            best = &region;
            bestScore = score;
        }
    }
    return best;
}

void writeFieldCsv(const std::filesystem::path& path,
                   const SentaurusTdrField& field,
                   const std::vector<std::size_t>& nodes)
{
    std::ofstream out(path);
    out << "node_id";
    for (std::size_t component = 0; component < field.component_count; ++component) {
        out << ",component" << component;
    }
    out << "\n";
    const std::size_t rows = std::min(nodes.size(), field.value_count);
    for (std::size_t row = 0; row < rows; ++row) {
        out << nodes[row];
        for (std::size_t component = 0; component < field.component_count; ++component) {
            out << "," << field.values[row * field.component_count + component];
        }
        out << "\n";
    }
}

} // namespace

const SentaurusTdrField* SentaurusTdrInventory::findField(const std::string& name, int regionIndex) const
{
    for (const auto& field : fields) {
        if (field.name == name && field.region_index == regionIndex) {
            return &field;
        }
    }
    return nullptr;
}

SentaurusTdrInventory SentaurusTdrReader::readInventory(const std::string& filename) const
{
    H5Object file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!file) {
        throw std::runtime_error("failed to open Sentaurus TDR/HDF5 file: " + filename);
    }
    H5Object geometry(H5Gopen2(file.id, "/collection/geometry_0", H5P_DEFAULT), H5Gclose);
    if (!geometry) {
        throw std::runtime_error("TDR file is missing /collection/geometry_0");
    }

    SentaurusTdrInventory inventory;
    inventory.vertices = readVertices(geometry.id);
    inventory.regions = readRegions(geometry.id);
    inventory.fields = readFields(geometry.id);
    return inventory;
}

void SentaurusTdrReader::exportNeutral(const std::string& filename, const std::string& outputDirectory) const
{
    const SentaurusTdrInventory inventory = readInventory(filename);
    const std::filesystem::path outDir(outputDirectory);
    std::filesystem::create_directories(outDir);
    std::filesystem::create_directories(outDir / "fields");

    {
        std::ofstream out(outDir / "nodes.csv");
        out << "id,x_um,y_um\n";
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            out << i << "," << inventory.vertices[i].x << "," << inventory.vertices[i].y << "\n";
        }
    }

    std::map<int, std::string> materialRegionNames;
    {
        std::ofstream out(outDir / "elements.csv");
        out << "id,node0,node1,node2,region,material\n";
        std::size_t cellId = 0;
        for (const auto& region : inventory.regions) {
            if (region.type != SentaurusTdrRegionType::Material) {
                continue;
            }
            materialRegionNames[region.index] = region.name;
            for (const auto& tri : region.triangles) {
                out << cellId++ << "," << tri[0] << "," << tri[1] << "," << tri[2]
                    << "," << region.name << "," << materialName(region.material) << "\n";
            }
        }
    }

    {
        std::ofstream out(outDir / "contacts.csv");
        out << "name,node_ids,region\n";
        for (const auto& region : inventory.regions) {
            if (region.type != SentaurusTdrRegionType::Contact) {
                continue;
            }
            std::set<std::size_t> nodes;
            for (const auto& edge : region.edges) {
                nodes.insert(edge[0]);
                nodes.insert(edge[1]);
            }
            const auto* owner = bestMaterialRegionForContact(inventory, region);
            out << region.name << ",";
            bool first = true;
            for (const auto node : nodes) {
                if (!first) {
                    out << ";";
                }
                out << node;
                first = false;
            }
            out << "," << (owner != nullptr ? owner->name : "") << "\n";
        }
    }

    std::vector<double> donors(inventory.vertices.size(), 0.0);
    std::vector<double> acceptors(inventory.vertices.size(), 0.0);
    for (const auto& field : inventory.fields) {
        if (field.name != "DonorConcentration" && field.name != "AcceptorConcentration") {
            continue;
        }
        const auto* region = findRegion(inventory, field.region_index);
        if (region == nullptr) {
            continue;
        }
        const auto nodes = regionNodeOrder(*region);
        const std::size_t rows = std::min(nodes.size(), field.value_count);
        for (std::size_t row = 0; row < rows; ++row) {
            if (nodes[row] >= inventory.vertices.size() || field.component_count == 0) {
                continue;
            }
            if (field.name == "DonorConcentration") {
                donors[nodes[row]] = field.values[row * field.component_count];
            } else {
                acceptors[nodes[row]] = field.values[row * field.component_count];
            }
        }
    }
    {
        std::ofstream out(outDir / "doping.csv");
        out << "node_id,donors_cm3,acceptors_cm3\n";
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            out << i << "," << donors[i] << "," << acceptors[i] << "\n";
        }
    }

    for (const auto& field : inventory.fields) {
        const auto* region = findRegion(inventory, field.region_index);
        if (region == nullptr || field.values.empty()) {
            continue;
        }
        const auto nodes = regionNodeOrder(*region);
        const auto fieldPath = outDir / "fields" /
            (sanitizeFilename(field.name) + "_region" + std::to_string(field.region_index) + ".csv");
        writeFieldCsv(fieldPath, field, nodes);
    }

    nlohmann::json metadata;
    metadata["source"] = filename;
    metadata["vertex_count"] = inventory.vertices.size();
    metadata["region_count"] = inventory.regions.size();
    metadata["dataset_count"] = inventory.fields.size();
    metadata["regions"] = nlohmann::json::array();
    for (const auto& region : inventory.regions) {
        metadata["regions"].push_back({
            {"index", region.index},
            {"name", region.name},
            {"material", region.material},
            {"type", static_cast<int>(region.type == SentaurusTdrRegionType::Other ? 99 : region.index >= 0 ? static_cast<int>(region.type) : 99)},
            {"triangles", region.triangles.size()},
            {"edges", region.edges.size()},
        });
    }
    metadata["fields"] = nlohmann::json::array();
    for (const auto& field : inventory.fields) {
        nlohmann::json fieldJson = {
            {"index", field.index},
            {"name", field.name},
            {"region", field.region_index},
            {"unit", field.unit},
            {"values", field.value_count},
            {"components", field.component_count},
        };
        if (field.value_count <= 4 && field.values.size() <= 16) {
            fieldJson["raw_values"] = field.values;
        }
        metadata["fields"].push_back(std::move(fieldJson));
    }
    std::ofstream metaOut(outDir / "metadata.json");
    metaOut << metadata.dump(2) << "\n";
}

} // namespace vela
