#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/io/SentaurusTdrReader.h"

#include <hdf5.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace vela;

namespace {

struct Vertex2D {
    double x;
    double y;
};

struct SyntheticField {
    std::string name;
    std::vector<double> values;
};

void writeStringAttribute(hid_t object, const char* name, const std::string& value)
{
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, value.size());
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(object, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, value.data());
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(type);
}

void writeIntAttribute(hid_t object, const char* name, int value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(object, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

void writeSizeAttribute(hid_t object, const char* name, hsize_t value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(object, name, H5T_NATIVE_HSIZE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_HSIZE, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

hid_t createGroup(hid_t parent, const char* name)
{
    return H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void writeIntDataset(hid_t group, const char* name, const std::vector<int>& values)
{
    hsize_t dims[] = {values.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate2(group, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    H5Dclose(dataset);
    H5Sclose(space);
}

void writeDoubleDatasetWithAttrs(hid_t group,
                                 const char* datasetName,
                                 const std::vector<double>& values,
                                 const std::string& name,
                                 int region,
                                 hsize_t numberOfValues,
                                 const std::string& unit)
{
    hsize_t dims[] = {values.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate2(group, datasetName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    writeStringAttribute(dataset, "name", name);
    writeIntAttribute(dataset, "region", region);
    writeSizeAttribute(dataset, "number of values", numberOfValues);
    writeStringAttribute(dataset, "unit", unit);
    H5Dclose(dataset);
    H5Sclose(space);
}

std::filesystem::path writeSyntheticTdr(const std::vector<SyntheticField>& dopingFields = {
    {"PhosphorusActiveConcentration", {1.0e17, 2.0e17, 3.0e17, 4.0e17}},
    {"BoronActiveConcentration", {0.0, 0.0, 1.0e16, 1.0e16}},
})
{
    const auto path = std::filesystem::temp_directory_path() / "vela_synthetic_sentaurus.tdr";
    std::filesystem::remove(path);

    hid_t file = H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t collection = createGroup(file, "collection");
    hid_t geometry = createGroup(collection, "geometry_0");

    const std::vector<Vertex2D> vertices = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
    };
    hsize_t vertexDims[] = {vertices.size()};
    hid_t vertexSpace = H5Screate_simple(1, vertexDims, nullptr);
    hid_t vertexType = H5Tcreate(H5T_COMPOUND, sizeof(Vertex2D));
    H5Tinsert(vertexType, "x", HOFFSET(Vertex2D, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(vertexType, "y", HOFFSET(Vertex2D, y), H5T_NATIVE_DOUBLE);
    hid_t vertexDataset = H5Dcreate2(geometry, "vertex", vertexType, vertexSpace,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vertexDataset, vertexType, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices.data());
    H5Dclose(vertexDataset);
    H5Tclose(vertexType);
    H5Sclose(vertexSpace);

    hid_t silicon = createGroup(geometry, "region_0");
    writeStringAttribute(silicon, "name", "Silicon_1");
    writeStringAttribute(silicon, "material", "Silicon");
    writeIntAttribute(silicon, "type", 0);
    writeSizeAttribute(silicon, "number of elements", 2);
    writeIntDataset(silicon, "elements_0", {2, 0, 1, 2, 2, 0, 2, 3});
    H5Gclose(silicon);

    hid_t drain = createGroup(geometry, "region_1");
    writeStringAttribute(drain, "name", "drain");
    writeIntAttribute(drain, "type", 1);
    writeSizeAttribute(drain, "number of elements", 2);
    writeIntDataset(drain, "elements_0", {1, 1, 2, 1, 2, 3});
    H5Gclose(drain);

    hid_t iface = createGroup(geometry, "region_2");
    writeStringAttribute(iface, "name", "Oxide_1+Silicon_1");
    writeIntAttribute(iface, "type", 2);
    writeSizeAttribute(iface, "number of elements", 3);
    writeIntDataset(iface, "elements_0", {0, 0, 0, 3, 1, 2, 3});
    H5Gclose(iface);

    hid_t state = createGroup(geometry, "state_0");
    int datasetIndex = 0;
    for (const auto& field : dopingFields) {
        const std::string datasetName = "dataset_" + std::to_string(datasetIndex++);
        hid_t dataset = createGroup(state, datasetName.c_str());
        writeDoubleDatasetWithAttrs(dataset, "values", field.values, field.name, 0, 4, "cm^-3");
        H5Gclose(dataset);
    }
    hid_t d2 = createGroup(state, ("dataset_" + std::to_string(datasetIndex++)).c_str());
    writeDoubleDatasetWithAttrs(d2, "values", {0.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
                                "ElectricField", 0, 4, "V*cm^-1");
    H5Gclose(d2);
    hid_t d3 = createGroup(state, ("dataset_" + std::to_string(datasetIndex++)).c_str());
    writeDoubleDatasetWithAttrs(d3, "values", {5.0}, "ContactExternalVoltage", 1, 1, "V");
    H5Gclose(d3);
    hid_t d4 = createGroup(state, ("dataset_" + std::to_string(datasetIndex++)).c_str());
    writeDoubleDatasetWithAttrs(d4, "values", {2.5e-3}, "ContactCurrentFlux", 1, 1, "A");
    H5Gclose(d4);
    hid_t d5 = createGroup(state, ("dataset_" + std::to_string(datasetIndex++)).c_str());
    writeDoubleDatasetWithAttrs(d5, "values", {7.0, 8.0, 9.0}, "MismatchedScalar", 0, 3, "1");
    H5Gclose(d5);
    H5Gclose(state);

    H5Gclose(geometry);
    H5Gclose(collection);
    H5Fclose(file);
    return path;
}

std::string readFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    return std::string(std::istreambuf_iterator<char>(input), {});
}

std::string exportSyntheticDopingCsv(const std::vector<SyntheticField>& dopingFields)
{
    const auto path = writeSyntheticTdr(dopingFields);
    const auto outDir = std::filesystem::temp_directory_path() / "vela_synthetic_sentaurus_export";
    std::filesystem::remove_all(outDir);

    SentaurusTdrReader reader;
    reader.exportNeutral(path.string(), outDir.string());
    return readFile(outDir / "doping.csv");
}

} // namespace

TEST_CASE("SentaurusTdrReader reads mesh regions contacts and state datasets", "[sentaurus][tdr]")
{
    const auto path = writeSyntheticTdr();

    SentaurusTdrReader reader;
    const SentaurusTdrInventory inventory = reader.readInventory(path.string());

    REQUIRE(inventory.vertices.size() == 4);
    REQUIRE(inventory.regions.size() == 3);
    REQUIRE(inventory.fields.size() == 6);

    const auto& silicon = inventory.regions.at(0);
    REQUIRE(silicon.name == "Silicon_1");
    REQUIRE(silicon.material == "Silicon");
    REQUIRE(silicon.type == SentaurusTdrRegionType::Material);
    REQUIRE(silicon.triangles.size() == 2);
    REQUIRE(silicon.triangles.at(0) == std::array<std::size_t, 3>{0, 1, 2});

    const auto& contact = inventory.regions.at(1);
    REQUIRE(contact.type == SentaurusTdrRegionType::Contact);
    REQUIRE(contact.edges.size() == 2);

    const auto& interface = inventory.regions.at(2);
    REQUIRE(interface.type == SentaurusTdrRegionType::Interface);
    REQUIRE(interface.edges.size() == 1);
    REQUIRE(interface.points.size() == 2);

    const auto* electricField = inventory.findField("ElectricField", 0);
    REQUIRE(electricField != nullptr);
    REQUIRE(electricField->component_count == 2);
    REQUIRE(electricField->values.size() == 8);

    const auto* contactVoltage = inventory.findField("ContactExternalVoltage", 1);
    REQUIRE(contactVoltage != nullptr);
    REQUIRE(contactVoltage->values.at(0) == Catch::Approx(5.0));
}

TEST_CASE("SentaurusTdrReader exports neutral reference TCAD CSV files", "[sentaurus][tdr]")
{
    const auto path = writeSyntheticTdr();
    const auto outDir = std::filesystem::temp_directory_path() / "vela_synthetic_sentaurus_export";
    std::filesystem::remove_all(outDir);

    SentaurusTdrReader reader;
    reader.exportNeutral(path.string(), outDir.string());

    REQUIRE(std::filesystem::is_regular_file(outDir / "nodes.csv"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "elements.csv"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "contacts.csv"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "doping.csv"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "fields" / "ElectricField_region0.csv"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "metadata.json"));
    REQUIRE(std::filesystem::is_regular_file(outDir / "field_manifest.json"));

    const std::string contacts = readFile(outDir / "contacts.csv");
    REQUIRE(contacts.find("drain,1;2;3,Silicon_1") != std::string::npos);

    const std::string doping = readFile(outDir / "doping.csv");
    REQUIRE(doping.find("0,1e+17,0") != std::string::npos);
    REQUIRE(doping.find("3,4e+17,1e+16") != std::string::npos);

    const auto dopingMetadata = nlohmann::json::parse(readFile(outDir / "doping_metadata.json"));
    REQUIRE(dopingMetadata["compensated_nodes"]["count"].get<int>() == 0);

    const std::string metadata = readFile(outDir / "metadata.json");
    REQUIRE(metadata.find("\"vertex_count\": 4") != std::string::npos);
    REQUIRE(metadata.find("\"dataset_count\": 6") != std::string::npos);

    const auto manifest = nlohmann::json::parse(readFile(outDir / "field_manifest.json"));
    REQUIRE(manifest["fields"].size() == 6);
    const auto electric = std::find_if(
        manifest["fields"].begin(), manifest["fields"].end(),
        [](const nlohmann::json& field) {
            return field["name"] == "ElectricField";
        });
    REQUIRE(electric != manifest["fields"].end());
    REQUIRE((*electric)["unit"] == "V*cm^-1");
    REQUIRE((*electric)["components"] == 2);
    REQUIRE((*electric)["global_node_mapping"] == "region_node_order");
    REQUIRE((*electric)["mapping_status"] == "complete");

    const auto mismatch = std::find_if(
        manifest["fields"].begin(), manifest["fields"].end(),
        [](const nlohmann::json& field) {
            return field["name"] == "MismatchedScalar";
        });
    REQUIRE(mismatch != manifest["fields"].end());
    REQUIRE((*mismatch)["mapping_status"] == "partial");
    REQUIRE((*mismatch)["warnings"][0].get<std::string>().find("value_count 3 does not match region node count 4")
            != std::string::npos);
}

TEST_CASE("SentaurusTdrReader reports compensated dopant nodes", "[sentaurus][tdr]")
{
    const auto path = writeSyntheticTdr({
        {"PhosphorusActiveConcentration", {1.0e17, 1.0e17, 3.0e17, 4.0e17}},
        {"BoronActiveConcentration", {0.0, 1.0e17, 1.0e16, 1.0e16}},
    });
    const auto outDir = std::filesystem::temp_directory_path() /
        "vela_synthetic_sentaurus_compensated_export";
    std::filesystem::remove_all(outDir);

    SentaurusTdrReader reader;
    reader.exportNeutral(path.string(), outDir.string());

    const auto metadata = nlohmann::json::parse(readFile(outDir / "doping_metadata.json"));
    REQUIRE(metadata["compensated_nodes"]["count"].get<int>() == 1);
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["node_id"].get<int>() == 1);
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["donors_cm3"].get<double>() ==
            Catch::Approx(1.0e17));
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["acceptors_cm3"].get<double>() ==
            Catch::Approx(1.0e17));
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["policy"].get<std::string>() == "reported");
}

TEST_CASE("SentaurusTdrReader resolves compensated dopant nodes from signed aggregate doping",
          "[sentaurus][tdr]")
{
    const auto path = writeSyntheticTdr({
        {"DopingConcentration", {1.0e17, 1.0e17, 3.0e17, 4.0e17}},
        {"PhosphorusActiveConcentration", {1.0e17, 1.0e17, 3.0e17, 4.0e17}},
        {"BoronActiveConcentration", {0.0, 1.0e17, 1.0e16, 1.0e16}},
    });
    const auto outDir = std::filesystem::temp_directory_path() /
        "vela_synthetic_sentaurus_resolved_export";
    std::filesystem::remove_all(outDir);

    SentaurusTdrExportOptions options;
    options.compensatedDopingPolicy = "dominant_signed_region";
    SentaurusTdrReader reader;
    reader.exportNeutral(path.string(), outDir.string(), options);

    const std::string doping = readFile(outDir / "doping.csv");
    REQUIRE(doping.find("1,1e+17,0") != std::string::npos);
    const auto metadata = nlohmann::json::parse(readFile(outDir / "doping_metadata.json"));
    REQUIRE(metadata["compensated_nodes"]["count"].get<int>() == 1);
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["policy"].get<std::string>() ==
            "dominant_signed_region");
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["resolved"].get<bool>());
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["resolved_donors_cm3"].get<double>() ==
            Catch::Approx(1.0e17));
    REQUIRE(metadata["compensated_nodes"]["nodes"][0]["resolved_acceptors_cm3"].get<double>() ==
            Catch::Approx(0.0));
}

TEST_CASE("SentaurusTdrReader exports legacy aggregate dopant fields", "[sentaurus][tdr]")
{
    const std::string doping = exportSyntheticDopingCsv({
        {"DonorConcentration", {1.0e17, 2.0e17, 3.0e17, 4.0e17}},
        {"AcceptorConcentration", {0.0, 0.0, 1.0e16, 1.0e16}},
    });

    REQUIRE(doping.find("0,1e+17,0") != std::string::npos);
    REQUIRE(doping.find("3,4e+17,1e+16") != std::string::npos);
}

TEST_CASE("SentaurusTdrReader sums active dopant species when no aggregate field is present", "[sentaurus][tdr]")
{
    const std::string doping = exportSyntheticDopingCsv({
        {"PhosphorusActiveConcentration", {1.0e17, 2.0e17, 3.0e17, 4.0e17}},
        {"ArsenicActiveConcentration", {5.0e16, 6.0e16, 7.0e16, 8.0e16}},
        {"BoronActiveConcentration", {1.0e16, 2.0e16, 3.0e16, 4.0e16}},
        {"AluminumActiveConcentration", {2.0e16, 3.0e16, 4.0e16, 5.0e16}},
    });

    REQUIRE(doping.find("0,1.5e+17,3e+16") != std::string::npos);
    REQUIRE(doping.find("3,4.8e+17,9e+16") != std::string::npos);
}

TEST_CASE("SentaurusTdrReader prefers aggregate dopant totals over active species for the same region", "[sentaurus][tdr]")
{
    const std::string doping = exportSyntheticDopingCsv({
        {"DonorConcentration", {1.0e17, 2.0e17, 3.0e17, 4.0e17}},
        {"PhosphorusActiveConcentration", {9.0e17, 9.0e17, 9.0e17, 9.0e17}},
        {"AcceptorConcentration", {0.0, 0.0, 1.0e16, 1.0e16}},
        {"BoronActiveConcentration", {9.0e16, 9.0e16, 9.0e16, 9.0e16}},
    });

    REQUIRE(doping.find("0,1e+17,0") != std::string::npos);
    REQUIRE(doping.find("3,4e+17,1e+16") != std::string::npos);
}
