#include "vela/io/SentaurusTdrReader.h"

#include <nlohmann/json.hpp>

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

nlohmann::json inventoryJson(const SentaurusTdrInventory& inventory)
{
    nlohmann::json data;
    data["vertex_count"] = inventory.vertices.size();
    data["region_count"] = inventory.regions.size();
    data["dataset_count"] = inventory.fields.size();
    data["regions"] = nlohmann::json::array();
    for (const auto& region : inventory.regions) {
        data["regions"].push_back({
            {"index", region.index},
            {"name", region.name},
            {"material", region.material},
            {"type", regionTypeCode(region.type)},
            {"triangles", region.triangles.size()},
            {"edges", region.edges.size()},
            {"points", region.points.size()},
        });
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
        if (field.value_count <= 4 && field.values.size() <= 16) {
            fieldJson["raw_values"] = field.values;
        }
        data["fields"].push_back(std::move(fieldJson));
    }
    return data;
}

void usage()
{
    std::cerr
        << "Usage: sentaurus_import --tdr FILE [--inventory-json FILE] [--export-dir DIR] "
           "[--compensated-doping-policy reported|dominant_signed_region]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::string tdrPath;
        std::string inventoryPath;
        std::string exportDir;
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
            } else if (arg == "--export-dir") {
                exportDir = requireValue("--export-dir");
            } else if (arg == "--compensated-doping-policy") {
                exportOptions.compensatedDopingPolicy = requireValue("--compensated-doping-policy");
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

        SentaurusTdrReader reader;
        if (!exportDir.empty()) {
            reader.exportNeutral(tdrPath, exportDir, exportOptions);
        }
        const auto inventory = reader.readInventory(tdrPath);
        const auto json = inventoryJson(inventory);
        if (!inventoryPath.empty()) {
            std::ofstream out(inventoryPath);
            out << json.dump(2) << "\n";
        } else {
            std::cout << json.dump(2) << "\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "sentaurus_import: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
