#include "vela/io/SentaurusTdrReader.h"

#include <hdf5.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
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
    // Keep the callback and iterator on the configured HDF5 API version.
    // HDF5 1.10 has no H5L_info2_t or H5Literate2.
    auto callback = [](hid_t, const char* name, const H5L_info_t*, void* data) -> herr_t {
        auto* out = static_cast<std::vector<std::string>*>(data);
        out->emplace_back(name);
        return 0;
    };
    hsize_t index = 0;
    H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, &index, callback, &names);
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
        if (kind == 0 && i < raw.size()) {
            region.points.push_back(static_cast<std::size_t>(raw[i]));
            i += 1;
        } else if (kind == 2 && i + 2 < raw.size()) {
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
        region.material_type = readIntAttribute(group.id, "material type", -1);
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
    field.location_type = readIntAttribute(
        values.id, "location type", readIntAttribute(group, "location type", -1));
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
    // Sentaurus stores region-supported nodal datasets in ascending global
    // vertex-id order, not in the order in which vertices first occur in the
    // region element stream.  The latter can be highly permuted on production
    // meshes and silently binds otherwise coherent state fields to the wrong
    // geometry nodes.
    std::set<std::size_t> nodes;
    for (const auto& tri : region.triangles) {
        nodes.insert(tri[0]);
        nodes.insert(tri[1]);
        nodes.insert(tri[2]);
    }
    for (const auto& edge : region.edges) {
        nodes.insert(edge[0]);
        nodes.insert(edge[1]);
    }
    return {nodes.begin(), nodes.end()};
}

std::vector<std::size_t> globalVertexOrder(std::size_t vertexCount)
{
    std::vector<std::size_t> order(vertexCount);
    std::iota(order.begin(), order.end(), std::size_t{0});
    return order;
}

bool usesGlobalVertexOrder(const SentaurusTdrInventory& inventory,
                           const SentaurusTdrField& field,
                           const SentaurusTdrRegion& region)
{
    return region.type != SentaurusTdrRegionType::Contact &&
        field.value_count == inventory.vertices.size();
}

bool usesRegionCellOrder(const SentaurusTdrField& field,
                         const SentaurusTdrRegion& region)
{
    constexpr int elementLocationType = 3;
    return field.location_type == elementLocationType &&
        region.type == SentaurusTdrRegionType::Material &&
        field.value_count == region.triangles.size();
}

std::vector<std::size_t> regionCellOrder(const SentaurusTdrInventory& inventory,
                                         const SentaurusTdrRegion& target)
{
    std::size_t offset = 0;
    for (const auto& region : inventory.regions) {
        if (region.type != SentaurusTdrRegionType::Material) {
            continue;
        }
        if (region.index == target.index) {
            std::vector<std::size_t> order(region.triangles.size());
            std::iota(order.begin(), order.end(), offset);
            return order;
        }
        offset += region.triangles.size();
    }
    return {};
}

std::vector<std::size_t> fieldNodeOrder(const SentaurusTdrInventory& inventory,
                                        const SentaurusTdrField& field,
                                        const SentaurusTdrRegion& region)
{
    if (usesRegionCellOrder(field, region)) {
        return {};
    }
    if (usesGlobalVertexOrder(inventory, field, region)) {
        return globalVertexOrder(inventory.vertices.size());
    }
    return regionNodeOrder(region);
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

bool isDonorConcentrationField(const std::string& name)
{
    return name == "DonorConcentration" ||
        name == "PActive" ||
        name == "AsActive" ||
        name == "PhosphorusActiveConcentration" ||
        name == "ArsenicActiveConcentration" ||
        name == "AntimonyActiveConcentration" ||
        name == "PActive" ||
        name == "AsActive" ||
        name == "SbActive";
}

bool isAggregateDonorConcentrationField(const std::string& name)
{
    return name == "DonorConcentration";
}

bool isAcceptorConcentrationField(const std::string& name)
{
    return name == "AcceptorConcentration" ||
        name == "BActive" ||
        name == "BoronActiveConcentration" ||
        name == "AluminumActiveConcentration" ||
        name == "IndiumActiveConcentration" ||
        name == "BActive" ||
        name == "AlActive" ||
        name == "InActive";
}

bool isAggregateAcceptorConcentrationField(const std::string& name)
{
    return name == "AcceptorConcentration";
}

bool isDopableMaterial(const std::string& material)
{
    return material == "Silicon" || material == "PolySilicon" ||
        material == "Germanium" || material == "SiliconGermanium";
}

bool isSemiconductorMaterialRegion(const SentaurusTdrRegion& region)
{
    return region.type == SentaurusTdrRegionType::Material &&
        (region.material_type != 0 || isDopableMaterial(region.material));
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
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
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

void writeCellFieldCsv(const std::filesystem::path& path,
                       const SentaurusTdrField& field,
                       const std::vector<std::size_t>& cells)
{
    std::ofstream out(path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "cell_id";
    for (std::size_t component = 0; component < field.component_count; ++component) {
        out << ",component" << component;
    }
    out << "\n";
    const std::size_t rows = std::min(cells.size(), field.value_count);
    for (std::size_t row = 0; row < rows; ++row) {
        out << cells[row];
        for (std::size_t component = 0; component < field.component_count; ++component) {
            out << "," << field.values[row * field.component_count + component];
        }
        out << "\n";
    }
}

nlohmann::json fieldManifestEntry(const SentaurusTdrInventory& inventory,
                                  const SentaurusTdrField& field)
{
    nlohmann::json entry = {
        {"index", field.index},
        {"name", field.name},
        {"region", field.region_index},
        {"location_type", field.location_type},
        {"support_kind", "unknown"},
        {"unit", field.unit},
        {"values", field.value_count},
        {"raw_value_count", field.values.size()},
        {"components", field.component_count},
        {"global_node_mapping", "none"},
        {"mapping_status", "unmapped"},
        {"warnings", nlohmann::json::array()},
    };

    const auto* region = findRegion(inventory, field.region_index);
    if (region == nullptr) {
        entry["warnings"].push_back(
            "field region " + std::to_string(field.region_index) + " is not present in the geometry");
        return entry;
    }

    entry["region_name"] = region->name;
    entry["region_type"] = static_cast<int>(region->type);
    const auto nodes = regionNodeOrder(*region);
    entry["region_node_count"] = nodes.size();
    entry["global_vertex_count"] = inventory.vertices.size();
    entry["region_cell_count"] = region->triangles.size();
    entry["support_kind"] = "node";
    entry["csv_file"] =
        sanitizeFilename(field.name) + "_region" + std::to_string(field.region_index) + ".csv";

    if (usesRegionCellOrder(field, *region)) {
        entry["support_kind"] = "cell";
        entry["global_node_mapping"] = "region_cell_order";
        entry["mapping_status"] = "complete";
        entry["csv_file"] =
            sanitizeFilename(field.name) + "_region" + std::to_string(field.region_index) +
            "_cells.csv";
        return entry;
    }

    if (region->type == SentaurusTdrRegionType::Contact && field.value_count == 1) {
        entry["support_kind"] = "contact";
        entry["global_node_mapping"] = "contact_scalar";
        entry["mapping_status"] = "scalar";
        return entry;
    }

    if (usesGlobalVertexOrder(inventory, field, *region)) {
        entry["global_node_mapping"] = "global_vertex_order";
        entry["mapping_status"] = "complete";
        return entry;
    }

    entry["global_node_mapping"] = "region_node_order";
    entry["region_node_ordering"] = "ascending_global_vertex_id";
    if (field.value_count == nodes.size()) {
        entry["mapping_status"] = "complete";
        return entry;
    }

    if (field.value_count == 0) {
        entry["mapping_status"] = "empty";
    } else {
        entry["mapping_status"] = "partial";
    }
    entry["warnings"].push_back(
        "value_count " + std::to_string(field.value_count) +
        " does not match region node count " + std::to_string(nodes.size()) +
        "; CSV export uses the overlapping prefix only");
    return entry;
}

nlohmann::json fieldManifest(const SentaurusTdrInventory& inventory)
{
    nlohmann::json manifest;
    manifest["fields"] = nlohmann::json::array();
    for (const auto& field : inventory.fields) {
        manifest["fields"].push_back(fieldManifestEntry(inventory, field));
    }
    return manifest;
}

double coordinateScaleToMicrometers(const std::string& unit)
{
    if (unit.empty() || unit == "um") {
        return 1.0;
    }
    if (unit == "cm") {
        return 1.0e4;
    }
    if (unit == "mm") {
        return 1.0e3;
    }
    if (unit == "m") {
        return 1.0e6;
    }
    if (unit == "nm") {
        return 1.0e-3;
    }
    throw std::runtime_error(
        "unsupported Sentaurus TDR coordinate unit for micrometer export: " + unit);
}

bool containsString(const std::vector<std::string>& values, const std::string& value)
{
    return std::find(values.begin(), values.end(), value) != values.end();
}

std::set<std::size_t> materialNodes(const SentaurusTdrRegion& region)
{
    const auto ordered = regionNodeOrder(region);
    return {ordered.begin(), ordered.end()};
}

bool fieldCoversRegion(const SentaurusTdrInventory& inventory,
                       const SentaurusTdrField& field,
                       const SentaurusTdrRegion& region)
{
    return field.component_count == 1 && field.values.size() == field.value_count &&
        (field.value_count == inventory.vertices.size() ||
         field.value_count == materialNodes(region).size());
}

bool fieldValuesAreFiniteAndNonnegative(const SentaurusTdrField& field)
{
    return std::all_of(field.values.begin(), field.values.end(), [](double value) {
        return std::isfinite(value) && value >= 0.0;
    });
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
    inventory.coordinate_unit = readUnitAttribute(geometry.id);
    inventory.vertices = readVertices(geometry.id);
    inventory.regions = readRegions(geometry.id);
    inventory.fields = readFields(geometry.id);
    return inventory;
}

SentaurusTdrQualificationReport SentaurusTdrReader::qualify(
    const SentaurusTdrInventory& inventory,
    const SentaurusTdrQualificationContract& contract) const
{
    SentaurusTdrQualificationReport report;
    auto addCheck = [&](std::string code, bool passed, std::string message) {
        report.checks.push_back({std::move(code), passed, std::move(message)});
    };

    bool coordinateUnitConvertible = false;
    try {
        (void)coordinateScaleToMicrometers(inventory.coordinate_unit);
        coordinateUnitConvertible = !inventory.coordinate_unit.empty();
    } catch (const std::runtime_error&) {
        coordinateUnitConvertible = false;
    }
    const bool coordinateUnitAccepted = coordinateUnitConvertible &&
        !inventory.coordinate_unit.empty() &&
        containsString(contract.acceptedCoordinateUnits, inventory.coordinate_unit);
    addCheck(
        "geometry.coordinate_unit",
        coordinateUnitAccepted,
        coordinateUnitAccepted
            ? "coordinate unit '" + inventory.coordinate_unit + "' is explicitly accepted"
            : "coordinate unit '" + inventory.coordinate_unit + "' is missing or not accepted");

    const bool finiteVertices = !inventory.vertices.empty() &&
        std::all_of(inventory.vertices.begin(), inventory.vertices.end(), [](const auto& vertex) {
            return std::isfinite(vertex.x) && std::isfinite(vertex.y);
        });
    addCheck(
        "geometry.vertices",
        finiteVertices,
        finiteVertices ? "all vertices are finite" : "geometry has no vertices or non-finite coordinates");

    std::vector<std::string> materialProblems;
    std::vector<const SentaurusTdrRegion*> materialRegions;
    for (const auto& region : inventory.regions) {
        if (region.type != SentaurusTdrRegionType::Material) {
            continue;
        }
        materialRegions.push_back(&region);
        if (region.name.empty() || region.material.empty()) {
            materialProblems.push_back("region " + std::to_string(region.index) + " has no name or material");
        }
        if (!contract.allowedMaterials.empty() &&
            !containsString(contract.allowedMaterials, region.material)) {
            materialProblems.push_back(
                "region '" + region.name + "' uses unaccepted material '" + region.material + "'");
        }
        if (region.material_type >= 0) {
            const int expectedMaterialType =
                containsString(contract.semiconductorMaterials, region.material) ? 1 : 0;
            if (region.material_type != expectedMaterialType) {
                materialProblems.push_back(
                    "region '" + region.name + "' has material type " +
                    std::to_string(region.material_type) + " but contract expects " +
                    std::to_string(expectedMaterialType));
            }
        }
        if (region.triangles.empty()) {
            materialProblems.push_back("region '" + region.name + "' has no triangles");
        }
        for (const auto& triangle : region.triangles) {
            if (triangle[0] >= inventory.vertices.size() ||
                triangle[1] >= inventory.vertices.size() ||
                triangle[2] >= inventory.vertices.size()) {
                materialProblems.push_back("region '" + region.name + "' has an out-of-range triangle");
                break;
            }
            const auto& a = inventory.vertices[triangle[0]];
            const auto& b = inventory.vertices[triangle[1]];
            const auto& c = inventory.vertices[triangle[2]];
            const double twiceArea =
                (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (!std::isfinite(twiceArea) || twiceArea == 0.0) {
                materialProblems.push_back("region '" + region.name + "' has a degenerate triangle");
                break;
            }
        }
    }
    const bool materialsValid = !materialRegions.empty() && materialProblems.empty();
    addCheck(
        "geometry.material_regions",
        materialsValid,
        materialsValid
            ? std::to_string(materialRegions.size()) + " material regions are valid"
            : (materialProblems.empty() ? "geometry has no material regions" : materialProblems.front()));

    std::set<std::string> contactNames;
    std::vector<std::string> contactProblems;
    for (const auto& region : inventory.regions) {
        if (region.type != SentaurusTdrRegionType::Contact) {
            continue;
        }
        const bool inserted = contactNames.insert(region.name).second;
        if (!inserted) {
            contactProblems.push_back("duplicate contact '" + region.name + "'");
        }
        if (region.name.empty() || region.edges.empty()) {
            contactProblems.push_back(
                "contact region " + std::to_string(region.index) + " has no name or edges");
            continue;
        }
        std::set<std::size_t> nodes;
        bool indicesValid = true;
        for (const auto& edge : region.edges) {
            indicesValid = indicesValid && edge[0] < inventory.vertices.size() &&
                edge[1] < inventory.vertices.size();
            nodes.insert(edge[0]);
            nodes.insert(edge[1]);
        }
        if (!indicesValid) {
            contactProblems.push_back("contact '" + region.name + "' has an out-of-range edge");
            continue;
        }
        bool touchesMaterial = false;
        for (const auto* materialRegion : materialRegions) {
            const auto supportedNodes = materialNodes(*materialRegion);
            if (std::any_of(nodes.begin(), nodes.end(), [&](std::size_t node) {
                    return supportedNodes.contains(node);
                })) {
                touchesMaterial = true;
                break;
            }
        }
        if (!touchesMaterial) {
            contactProblems.push_back("contact '" + region.name + "' does not touch a material region");
        }
    }
    std::set<std::string> requiredContacts(
        contract.requiredContactNames.begin(), contract.requiredContactNames.end());
    std::vector<std::string> missingContacts;
    std::set_difference(
        requiredContacts.begin(), requiredContacts.end(),
        contactNames.begin(), contactNames.end(),
        std::back_inserter(missingContacts));
    std::vector<std::string> extraContacts;
    if (contract.requireExactContactSet) {
        std::set_difference(
            contactNames.begin(), contactNames.end(),
            requiredContacts.begin(), requiredContacts.end(),
            std::back_inserter(extraContacts));
    }
    const bool contactsValid =
        contactProblems.empty() && missingContacts.empty() && extraContacts.empty();
    std::string contactMessage = std::to_string(contactNames.size()) + " contact regions are valid";
    if (!missingContacts.empty()) {
        contactMessage = "missing required contact '" + missingContacts.front() + "'";
    } else if (!extraContacts.empty()) {
        contactMessage = "unexpected contact '" + extraContacts.front() + "'";
    } else if (!contactProblems.empty()) {
        contactMessage = contactProblems.front();
    }
    addCheck("geometry.contacts", contactsValid, contactMessage);

    std::vector<std::string> dopingCoverageProblems;
    std::vector<std::string> dopingUnitProblems;
    std::size_t semiconductorRegionCount = 0;
    for (const auto* region : materialRegions) {
        if (!containsString(contract.semiconductorMaterials, region->material)) {
            continue;
        }
        ++semiconductorRegionCount;
        bool donorCoverage = false;
        bool acceptorCoverage = false;
        for (const auto& field : inventory.fields) {
            if (field.region_index != region->index ||
                (!isDonorConcentrationField(field.name) &&
                 !isAcceptorConcentrationField(field.name))) {
                continue;
            }
            if (!containsString(contract.acceptedDopingUnits, field.unit)) {
                dopingUnitProblems.push_back(
                    "field '" + field.name + "' has unaccepted unit '" + field.unit + "'");
                continue;
            }
            if (!fieldCoversRegion(inventory, field, *region) ||
                !fieldValuesAreFiniteAndNonnegative(field)) {
                dopingCoverageProblems.push_back(
                    "field '" + field.name + "' does not completely cover region '" +
                    region->name + "' with finite nonnegative scalar values");
                continue;
            }
            donorCoverage = donorCoverage || isDonorConcentrationField(field.name);
            acceptorCoverage = acceptorCoverage || isAcceptorConcentrationField(field.name);
        }
        if (contract.requireCompleteSemiconductorDoping &&
            (!donorCoverage || !acceptorCoverage)) {
            dopingCoverageProblems.push_back(
                "semiconductor region '" + region->name +
                "' lacks complete donor and acceptor coverage");
        }
    }
    if (semiconductorRegionCount == 0) {
        dopingCoverageProblems.push_back("no configured semiconductor material region is present");
    }
    addCheck(
        "doping.units",
        dopingUnitProblems.empty(),
        dopingUnitProblems.empty() ? "all selected doping fields use accepted units"
                                   : dopingUnitProblems.front());
    addCheck(
        "doping.coverage",
        dopingCoverageProblems.empty(),
        dopingCoverageProblems.empty()
            ? std::to_string(semiconductorRegionCount) +
                " semiconductor regions have complete donor and acceptor coverage"
            : dopingCoverageProblems.front());

    report.passed = std::all_of(report.checks.begin(), report.checks.end(), [](const auto& check) {
        return check.passed;
    });
    return report;
}

void SentaurusTdrReader::exportNeutral(const std::string& filename, const std::string& outputDirectory) const
{
    exportNeutral(filename, outputDirectory, SentaurusTdrExportOptions{});
}

void SentaurusTdrReader::exportNeutral(const std::string& filename,
                                       const std::string& outputDirectory,
                                       const SentaurusTdrExportOptions& options) const
{
    const SentaurusTdrInventory inventory = readInventory(filename);
    if (!options.coordinateUnit.empty() && options.coordinateUnit != "um" &&
        options.coordinateUnit != "cm") {
        throw std::invalid_argument("SentaurusTdrReader: coordinateUnit override must be 'um' or 'cm'.");
    }
    const std::string sourceCoordinateUnit = options.coordinateUnit.empty()
        ? inventory.coordinate_unit : options.coordinateUnit;
    const double coordinateScale = coordinateScaleToMicrometers(sourceCoordinateUnit);
    const double coordinateToUm = coordinateScale;
    const std::filesystem::path outDir(outputDirectory);
    std::filesystem::create_directories(outDir);
    std::filesystem::create_directories(outDir / "fields");

    {
        std::ofstream out(outDir / "nodes.csv");
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        out << "id,x_um,y_um\n";
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            out << i << "," << inventory.vertices[i].x * coordinateScale << ","
                << inventory.vertices[i].y * coordinateScale << "\n";
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
        out << "name,node_ids,edge_node_ids,region\n";
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
            out << ",";
            first = true;
            for (const auto& edge : region.edges) {
                if (!first)
                    out << ";";
                out << edge[0] << "-" << edge[1];
                first = false;
            }
            out << "," << (owner != nullptr ? owner->name : "") << "\n";
        }
    }

    std::vector<double> donors(inventory.vertices.size(), 0.0);
    std::vector<double> acceptors(inventory.vertices.size(), 0.0);
    std::vector<std::vector<double>> signedDopingByNode(inventory.vertices.size());
    std::vector<bool> signedDopingTie(inventory.vertices.size(), false);
    std::vector<std::string> resolutionSource(inventory.vertices.size(), "reported");
    std::set<int> aggregateDonorRegions;
    std::set<int> aggregateAcceptorRegions;
    for (const auto& field : inventory.fields) {
        const auto* region = findRegion(inventory, field.region_index);
        if (region == nullptr || !isSemiconductorMaterialRegion(*region)) {
            continue;
        }
        if (isAggregateDonorConcentrationField(field.name)) {
            aggregateDonorRegions.insert(field.region_index);
        } else if (isAggregateAcceptorConcentrationField(field.name)) {
            aggregateAcceptorRegions.insert(field.region_index);
        }
    }
    for (const auto& field : inventory.fields) {
        const bool donorField = isDonorConcentrationField(field.name);
        const bool acceptorField = isAcceptorConcentrationField(field.name);
        const auto* region = findRegion(inventory, field.region_index);
        if (region == nullptr || !isSemiconductorMaterialRegion(*region)) {
            continue;
        }
        const auto nodes = fieldNodeOrder(inventory, field, *region);
        const std::size_t rows = std::min(nodes.size(), field.value_count);
        if (field.name == "DopingConcentration" || field.name == "NetActive") {
            for (std::size_t row = 0; row < rows; ++row) {
                if (nodes[row] >= inventory.vertices.size() || field.component_count == 0) {
                    continue;
                }
                signedDopingByNode[nodes[row]].push_back(
                    field.values[row * field.component_count]);
            }
        }
        if (!donorField && !acceptorField) {
            continue;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (nodes[row] >= inventory.vertices.size() || field.component_count == 0) {
                continue;
            }
            if (donorField) {
                if (isAggregateDonorConcentrationField(field.name)) {
                    donors[nodes[row]] = field.values[row * field.component_count];
                } else if (!aggregateDonorRegions.contains(field.region_index)) {
                    donors[nodes[row]] += field.values[row * field.component_count];
                }
            } else {
                if (isAggregateAcceptorConcentrationField(field.name)) {
                    acceptors[nodes[row]] = field.values[row * field.component_count];
                } else if (!aggregateAcceptorRegions.contains(field.region_index)) {
                    acceptors[nodes[row]] += field.values[row * field.component_count];
                }
            }
        }
    }

    const std::vector<double> originalDonors = donors;
    const std::vector<double> originalAcceptors = acceptors;
    auto isCompensated = [](double donor, double acceptor) {
        const double scale = std::max({std::abs(donor), std::abs(acceptor), 1.0});
        return donor > 0.0 && acceptor > 0.0 &&
            std::abs(donor - acceptor) <= 1.0e-6 * scale;
    };
    std::vector<std::set<std::size_t>> nodeAdjacency(inventory.vertices.size());
    auto addAdjacency = [&](std::size_t a, std::size_t b) {
        if (a >= inventory.vertices.size() || b >= inventory.vertices.size() || a == b) {
            return;
        }
        nodeAdjacency[a].insert(b);
        nodeAdjacency[b].insert(a);
    };
    for (const auto& region : inventory.regions) {
        if (region.type != SentaurusTdrRegionType::Material) {
            continue;
        }
        for (const auto& tri : region.triangles) {
            addAdjacency(tri[0], tri[1]);
            addAdjacency(tri[1], tri[2]);
            addAdjacency(tri[2], tri[0]);
        }
        for (const auto& edge : region.edges) {
            addAdjacency(edge[0], edge[1]);
        }
    }
    auto dominantNeighbourSignedDoping = [&](std::size_t node) {
        double positive = 0.0;
        double negative = 0.0;
        for (const std::size_t neighbour : nodeAdjacency[node]) {
            const double donor = originalDonors[neighbour];
            const double acceptor = originalAcceptors[neighbour];
            const double net = donor - acceptor;
            const double scale = std::max({std::abs(donor), std::abs(acceptor), 1.0});
            if (std::abs(net) <= 1.0e-6 * scale) {
                continue;
            }
            if (net > 0.0) {
                positive += net;
            } else {
                negative += -net;
            }
        }
        const double magnitude = std::max(originalDonors[node], originalAcceptors[node]);
        if (magnitude <= 0.0) {
            return 0.0;
        }
        if (positive > negative) {
            return magnitude;
        }
        if (negative > positive) {
            return -magnitude;
        }
        return 0.0;
    };

    auto dominantSignedAggregateDoping = [&](std::size_t node) {
        double positive = 0.0;
        double negative = 0.0;
        double firstNonzero = 0.0;
        for (const double value : signedDopingByNode[node]) {
            if (firstNonzero == 0.0 && value != 0.0) {
                firstNonzero = value;
            }
            if (value > positive) {
                positive = value;
            } else if (-value > negative) {
                negative = -value;
            }
        }
        const double scale = std::max({positive, negative, 1.0});
        if (positive > negative + 1.0e-6 * scale) {
            return positive;
        }
        if (negative > positive + 1.0e-6 * scale) {
            return -negative;
        }
        if (positive > 0.0 && negative > 0.0) {
            signedDopingTie[node] = true;
            return firstNonzero;
        }
        return 0.0;
    };

    if (options.compensatedDopingPolicy == "dominant_signed_region") {
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            if (!isCompensated(donors[i], acceptors[i])) {
                continue;
            }
            const bool hasSignedDoping = !signedDopingByNode[i].empty();
            double signedPick = dominantSignedAggregateDoping(i);
            if (signedPick != 0.0) {
                resolutionSource[i] = signedDopingTie[i]
                    ? "signed_aggregate_tie_first_region"
                    : "signed_aggregate_doping";
            }
            if (signedPick == 0.0 && hasSignedDoping && !signedDopingTie[i]) {
                resolutionSource[i] = "signed_aggregate_zero";
                continue;
            }
            if (signedPick == 0.0) {
                signedPick = dominantNeighbourSignedDoping(i);
                if (signedPick != 0.0) {
                    resolutionSource[i] = "neighbour_region_sign";
                } else if (signedDopingTie[i]) {
                    resolutionSource[i] = "unresolved_tie";
                }
            }
            if (signedPick > 0.0) {
                donors[i] = std::abs(signedPick);
                acceptors[i] = 0.0;
            } else if (signedPick < 0.0) {
                donors[i] = 0.0;
                acceptors[i] = std::abs(signedPick);
            }
        }
    } else if (options.compensatedDopingPolicy != "reported") {
        throw std::invalid_argument(
            "SentaurusTdrExportOptions: compensatedDopingPolicy must be 'reported' or "
            "'dominant_signed_region'.");
    }

    {
        std::ofstream out(outDir / "doping.csv");
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        out << "node_id,donors_cm3,acceptors_cm3\n";
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            out << i << "," << donors[i] << "," << acceptors[i] << "\n";
        }
    }
    {
        nlohmann::json dopingMetadata;
        dopingMetadata["schema"] = "vela.sentaurus_tdr.doping_metadata.v1";
        dopingMetadata["compensated_nodes"] = {
            {"policy", options.compensatedDopingPolicy},
            {"count", 0},
            {"nodes", nlohmann::json::array()},
        };
        for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
            const double donor = originalDonors[i];
            const double acceptor = originalAcceptors[i];
            if (isCompensated(donor, acceptor)) {
                const bool resolved =
                    donors[i] != originalDonors[i] || acceptors[i] != originalAcceptors[i];
                dopingMetadata["compensated_nodes"]["nodes"].push_back({
                    {"node_id", i},
                    {"donors_cm3", donor},
                    {"acceptors_cm3", acceptor},
                    {"policy", options.compensatedDopingPolicy},
                    {"resolved", resolved},
                    {"resolved_donors_cm3", donors[i]},
                    {"resolved_acceptors_cm3", acceptors[i]},
                    {"resolution_source", resolutionSource[i]},
                    {"signed_doping_tie", signedDopingTie[i]},
                });
            }
        }
        dopingMetadata["compensated_nodes"]["count"] =
            dopingMetadata["compensated_nodes"]["nodes"].size();

        std::ofstream out(outDir / "doping_metadata.json");
        out << dopingMetadata.dump(2) << "\n";
    }

    for (const auto& field : inventory.fields) {
        const auto* region = findRegion(inventory, field.region_index);
        if (region == nullptr || field.values.empty()) {
            continue;
        }
        if (usesRegionCellOrder(field, *region)) {
            const auto cells = regionCellOrder(inventory, *region);
            const auto fieldPath = outDir / "fields" /
                (sanitizeFilename(field.name) + "_region" +
                 std::to_string(field.region_index) + "_cells.csv");
            writeCellFieldCsv(fieldPath, field, cells);
            continue;
        }
        const auto nodes = fieldNodeOrder(inventory, field, *region);
        const auto fieldPath = outDir / "fields" /
            (sanitizeFilename(field.name) + "_region" + std::to_string(field.region_index) + ".csv");
        writeFieldCsv(fieldPath, field, nodes);
    }

    nlohmann::json metadata;
    metadata["source"] = filename;
    metadata["source_coordinate_unit"] = sourceCoordinateUnit;
    metadata["export_coordinate_unit"] = "um";
    metadata["coordinate_scale_to_um"] = coordinateScale;
    metadata["source_coordinate_unit"] = sourceCoordinateUnit;
    metadata["exported_coordinate_unit"] = "um";
    metadata["coordinate_to_um_scale"] = coordinateToUm;
    metadata["vertex_count"] = inventory.vertices.size();
    metadata["region_count"] = inventory.regions.size();
    metadata["dataset_count"] = inventory.fields.size();
    metadata["regions"] = nlohmann::json::array();
    for (const auto& region : inventory.regions) {
        metadata["regions"].push_back({
            {"index", region.index},
            {"name", region.name},
            {"material", region.material},
            {"material_type", region.material_type},
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
            {"location_type", field.location_type},
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

    std::ofstream manifestOut(outDir / "field_manifest.json");
    manifestOut << fieldManifest(inventory).dump(2) << "\n";
}

} // namespace vela
