#include "vela/io/SentaurusTdrReader.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace vela;

namespace {

int regionTypeCode(SentaurusTdrRegionType type)
{
    switch (type) {
    case SentaurusTdrRegionType::Material:
        return 0;
    case SentaurusTdrRegionType::Contact:
        return 1;
    case SentaurusTdrRegionType::Interface:
        return 2;
    case SentaurusTdrRegionType::Other:
        return 99;
    }
    return 99;
}

nlohmann::json inventoryJson(const SentaurusTdrInventory& inventory,
                             bool includeFieldValues = false)
{
    nlohmann::json data;
    data["coordinate_unit"] = inventory.coordinate_unit;
    data["vertex_count"] = inventory.vertices.size();
    data["region_count"] = inventory.regions.size();
    data["dataset_count"] = inventory.fields.size();
    data["regions"] = nlohmann::json::array();
    for (const auto& region : inventory.regions) {
        data["regions"].push_back({
            {"index", region.index},
            {"name", region.name},
            {"material", region.material},
            {"material_type", region.material_type},
            {"type", regionTypeCode(region.type)},
            {"triangles", region.triangles.size()},
            {"edges", region.edges.size()},
            {"points", region.points.size()},
        });
    }
    data["geometry"] = {
        {"vertices", nlohmann::json::array()},
        {"regions", nlohmann::json::array()},
    };
    for (const auto& vertex : inventory.vertices) {
        data["geometry"]["vertices"].push_back({vertex.x, vertex.y});
    }
    for (const auto& region : inventory.regions) {
        nlohmann::json regionJson = {
            {"index", region.index},
            {"name", region.name},
            {"material", region.material},
            {"material_type", region.material_type},
            {"type", regionTypeCode(region.type)},
            {"triangles", nlohmann::json::array()},
            {"edges", nlohmann::json::array()},
            {"points", region.points},
        };
        for (const auto& triangle : region.triangles) {
            regionJson["triangles"].push_back({triangle[0], triangle[1], triangle[2]});
        }
        for (const auto& edge : region.edges) {
            regionJson["edges"].push_back({edge[0], edge[1]});
        }
        data["geometry"]["regions"].push_back(std::move(regionJson));
    }
    data["fields"] = nlohmann::json::array();
    for (const auto& field : inventory.fields) {
        nlohmann::json fieldJson = {
            {"index", field.index},
            {"name", field.name},
            {"region", field.region_index},
            {"unit", field.unit},
            {"values", field.value_count},
            {"components", field.component_count},
        };
        if (includeFieldValues || (field.value_count <= 4 && field.values.size() <= 16)) {
            fieldJson["raw_values"] = field.values;
        }
        data["fields"].push_back(std::move(fieldJson));
    }
    return data;
}

SentaurusTdrQualificationContract qualificationContractFromJson(
    const nlohmann::json& document)
{
    const nlohmann::json& data = document.contains("input_tdr_qualification")
        ? document.at("input_tdr_qualification")
        : document;
    if (data.value("schema", std::string{}) != "vela.sentaurus_tdr.qualification.v1") {
        throw std::runtime_error(
            "qualification contract schema must be vela.sentaurus_tdr.qualification.v1");
    }
    if (data.value("export_coordinate_unit", std::string{"um"}) != "um") {
        throw std::runtime_error(
            "qualification contract export_coordinate_unit must be um");
    }

    SentaurusTdrQualificationContract contract;
    contract.requiredContactNames =
        data.at("required_contacts").get<std::vector<std::string>>();
    contract.requireExactContactSet = data.value("exact_contact_set", true);
    contract.allowedMaterials =
        data.at("allowed_materials").get<std::vector<std::string>>();
    contract.semiconductorMaterials =
        data.at("semiconductor_materials").get<std::vector<std::string>>();
    contract.acceptedCoordinateUnits =
        data.at("accepted_coordinate_units").get<std::vector<std::string>>();
    contract.acceptedDopingUnits =
        data.at("accepted_doping_units").get<std::vector<std::string>>();
    contract.requireCompleteSemiconductorDoping =
        data.value("require_complete_semiconductor_doping", true);
    return contract;
}

nlohmann::json qualificationJson(
    const SentaurusTdrQualificationReport& report,
    const SentaurusTdrInventory& inventory,
    const std::string& source,
    const std::string& contractPath)
{
    nlohmann::json data = {
        {"schema", "vela.sentaurus_tdr.qualification_report.v1"},
        {"source", source},
        {"contract", contractPath},
        {"passed", report.passed},
        {"inventory", {
            {"coordinate_unit", inventory.coordinate_unit},
            {"vertex_count", inventory.vertices.size()},
            {"region_count", inventory.regions.size()},
            {"dataset_count", inventory.fields.size()},
        }},
        {"checks", nlohmann::json::array()},
    };
    for (const auto& check : report.checks) {
        data["checks"].push_back({
            {"code", check.code},
            {"passed", check.passed},
            {"message", check.message},
        });
    }
    return data;
}

void writeJsonFile(const std::string& path, const nlohmann::json& data)
{
    const std::filesystem::path output(path);
    if (!output.parent_path().empty()) {
        std::filesystem::create_directories(output.parent_path());
    }
    std::ofstream out(output);
    if (!out.is_open()) {
        throw std::runtime_error("cannot open JSON output: " + path);
    }
    out << data.dump(2) << "\n";
}

void usage()
{
    std::cerr
        << "Usage: sentaurus_import --tdr FILE [--inventory-json FILE] [--export-dir DIR] "
           "[--field-values-json FILE] "
           "[--qualification-contract FILE --qualification-report FILE] "
           "[--coordinate-unit um|cm] "
           "[--compensated-doping-policy reported|dominant_signed_region]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::string tdrPath;
        std::string inventoryPath;
        std::string fieldValuesPath;
        std::string exportDir;
        std::string qualificationContractPath;
        std::string qualificationReportPath;
        SentaurusTdrExportOptions exportOptions;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto requireValue = [&](const char* option) -> std::string {
                if (i + 1 >= argc) {
                    throw std::runtime_error(std::string("missing value for ") + option);
                }
                return argv[++i];
            };
            if (arg == "--tdr") {
                tdrPath = requireValue("--tdr");
            } else if (arg == "--inventory-json") {
                inventoryPath = requireValue("--inventory-json");
            } else if (arg == "--field-values-json") {
                fieldValuesPath = requireValue("--field-values-json");
            } else if (arg == "--export-dir") {
                exportDir = requireValue("--export-dir");
            } else if (arg == "--qualification-contract") {
                qualificationContractPath = requireValue("--qualification-contract");
            } else if (arg == "--qualification-report") {
                qualificationReportPath = requireValue("--qualification-report");
            } else if (arg == "--compensated-doping-policy") {
                exportOptions.compensatedDopingPolicy = requireValue("--compensated-doping-policy");
            } else if (arg == "--coordinate-unit") {
                exportOptions.coordinateUnit = requireValue("--coordinate-unit");
            } else if (arg == "--help" || arg == "-h") {
                usage();
                return 0;
            } else {
                throw std::runtime_error("unknown argument: " + arg);
            }
        }
        if (tdrPath.empty()) {
            usage();
            return 2;
        }
        if (qualificationContractPath.empty() != qualificationReportPath.empty()) {
            throw std::runtime_error(
                "--qualification-contract and --qualification-report must be provided together");
        }

        SentaurusTdrReader reader;
        const auto inventory = reader.readInventory(tdrPath);
        const auto json = inventoryJson(inventory);
        if (!inventoryPath.empty()) {
            writeJsonFile(inventoryPath, json);
        } else {
            std::cout << json.dump(2) << "\n";
        }
        if (!qualificationContractPath.empty()) {
            std::ifstream contractInput(qualificationContractPath);
            if (!contractInput.is_open()) {
                throw std::runtime_error(
                    "cannot open qualification contract: " + qualificationContractPath);
            }
            nlohmann::json contractDocument;
            contractInput >> contractDocument;
            const auto contract = qualificationContractFromJson(contractDocument);
            const auto report = reader.qualify(inventory, contract);
            writeJsonFile(
                qualificationReportPath,
                qualificationJson(
                    report, inventory, tdrPath, qualificationContractPath));
            if (!report.passed) {
                throw std::runtime_error(
                    "TDR qualification failed; see " + qualificationReportPath);
            }
        }
        if (!exportDir.empty()) {
            reader.exportNeutral(tdrPath, exportDir, exportOptions);
        }
        if (!fieldValuesPath.empty()) {
            writeJsonFile(fieldValuesPath, inventoryJson(inventory, true));
        }
    } catch (const std::exception& ex) {
        std::cerr << "sentaurus_import: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
