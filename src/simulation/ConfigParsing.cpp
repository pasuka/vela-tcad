#include "vela/simulation/ConfigParsing.h"
#include "vela/io/CsvUtils.h"
#include "vela/mesh/DeviceMesh.h"

#include <nlohmann/json.hpp>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace vela {
namespace {

void appendConfigFixedChargeSpec(
    std::vector<RegionFixedChargeSpec>& specs,
    std::unordered_map<std::string, std::string>& sourcesByRegion,
    std::string region,
    Real fixedCharge,
    std::string source)
{
    const auto [_, inserted] = sourcesByRegion.emplace(region, source);
    if (!inserted) {
        throw std::runtime_error(
            "ConfigParsing: duplicate fixed_charge_m3 for region '" + region +
            "' from " + sourcesByRegion.at(region) + " and " + source +
            ". Specify fixed charge for each region only once.");
    }

    specs.push_back(RegionFixedChargeSpec{std::move(region), fixedCharge});
}

} // namespace

std::vector<RegionFixedChargeSpec> parseRegionFixedChargeSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling)
{
    std::vector<RegionFixedChargeSpec> specs;
    std::unordered_map<std::string, std::string> sourcesByRegion;

    if (cfg.contains("doping")) {
        for (const auto& entry : cfg.at("doping")) {
            if (!entry.contains("fixed_charge_m3")) continue;
            appendConfigFixedChargeSpec(
                specs,
                sourcesByRegion,
                entry.at("region").get<std::string>(),
                scaling.concentrationToInternal(entry.at("fixed_charge_m3").get<Real>()),
                "doping entry");
        }
    }

    if (cfg.contains("regions")) {
        for (const auto& entry : cfg.at("regions")) {
            if (!entry.contains("fixed_charge_m3")) continue;
            appendConfigFixedChargeSpec(
                specs,
                sourcesByRegion,
                entry.at("name").get<std::string>(),
                scaling.concentrationToInternal(entry.at("fixed_charge_m3").get<Real>()),
                "regions entry");
        }
    }

    return specs;
}

std::vector<RegionDopingSpec> parseDopingSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling)
{
    std::vector<RegionDopingSpec> specs;
    for (const auto& entry : cfg.at("doping")) {
        RegionDopingSpec spec;
        spec.region = entry.at("region").get<std::string>();
        spec.donors = scaling.concentrationToInternal(entry.at("donors").get<Real>());
        spec.acceptors = scaling.concentrationToInternal(entry.at("acceptors").get<Real>());
        specs.push_back(std::move(spec));
    }
    return specs;
}

std::vector<InterfaceSheetChargeSpec> parseInterfaceSheetChargeSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling)
{
    std::vector<InterfaceSheetChargeSpec> specs;
    if (!cfg.contains("interfaces"))
        return specs;

    for (const auto& entry : cfg.at("interfaces")) {
        const bool hasTrapOccupancy = entry.contains("trap_occupancy");
        const bool hasTrapDensity = entry.contains("trap_density_m2");
        if (hasTrapOccupancy && !hasTrapDensity) {
            throw std::runtime_error(
                "ConfigParsing: trap_occupancy requires trap_density_m2.");
        }

        const Real trapOccupancy = entry.value("trap_occupancy", 0.0);
        if (hasTrapOccupancy && (trapOccupancy < 0.0 || trapOccupancy > 1.0)) {
            throw std::runtime_error(
                "ConfigParsing: trap_occupancy must be in [0, 1].");
        }

        if (!entry.contains("sheet_charge_m2") &&
            !entry.contains("fixed_charge_m2") &&
            !hasTrapDensity) continue;

        const Real sheetCharge = scaling.sheetDensityToInternal(
            entry.value("sheet_charge_m2", 0.0));
        const Real fixedCharge = scaling.sheetDensityToInternal(
            entry.value("fixed_charge_m2", 0.0));
        const Real trapDensity = scaling.sheetDensityToInternal(
            entry.value("trap_density_m2", 0.0));

        if (entry.contains("regions")) {
            const auto regions = entry.at("regions").get<std::vector<std::string>>();
            if (regions.size() != 2)
                throw std::runtime_error(
                    "ConfigParsing: interface regions must contain exactly two names.");
            specs.push_back(InterfaceSheetChargeSpec{
                regions[0], regions[1], sheetCharge, fixedCharge, trapDensity, trapOccupancy});
        } else {
            specs.push_back(InterfaceSheetChargeSpec{
                entry.at("region0").get<std::string>(),
                entry.at("region1").get<std::string>(),
                sheetCharge,
                fixedCharge,
                trapDensity,
                trapOccupancy});
        }
    }

    return specs;
}

BoxGeometryBuilder::Options parseBoxGeometryOptions(const nlohmann::json& cfg)
{
    BoxGeometryBuilder::Options options;
    if (!cfg.contains("mesh_geometry"))
        return options;

    const auto& geometry = cfg.at("mesh_geometry");
    if (!geometry.is_object())
        throw std::runtime_error("ConfigParsing: mesh_geometry must be an object.");

    const std::string policy = geometry.value("node_volume_policy", "barycentric");
    if (policy == "barycentric") {
        options.nodeVolumePolicy = BoxGeometryBuilder::NodeVolumePolicy::Barycentric;
    } else if (policy == "mixed_voronoi") {
        options.nodeVolumePolicy = BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    } else {
        throw std::runtime_error(
            "ConfigParsing: mesh_geometry.node_volume_policy must be "
            "'barycentric' or 'mixed_voronoi'.");
    }

    if (geometry.contains("require_non_obtuse") &&
        !geometry.at("require_non_obtuse").is_boolean()) {
        throw std::runtime_error(
            "ConfigParsing: mesh_geometry.require_non_obtuse must be boolean.");
    }
    options.requireNonObtuse = geometry.value("require_non_obtuse", false);

    if (geometry.contains("fallback_negative_cotangent") &&
        !geometry.at("fallback_negative_cotangent").is_boolean()) {
        throw std::runtime_error(
            "ConfigParsing: mesh_geometry.fallback_negative_cotangent must be boolean.");
    }
    options.fallbackNegativeCotangent =
        geometry.value("fallback_negative_cotangent", true);

    return options;
}

CarrierTransportCoupleProfileReport applyCarrierTransportCoupleProfile(
    DeviceMesh& mesh,
    const nlohmann::json& cfg,
    const std::filesystem::path& configDirectory,
    UnitScalingConfig scaling)
{
    CarrierTransportCoupleProfileReport report;
    if (!cfg.contains("mesh_geometry"))
        return report;

    const auto& geometry = cfg.at("mesh_geometry");
    if (!geometry.is_object())
        throw std::runtime_error("ConfigParsing: mesh_geometry must be an object.");
    report.profile = geometry.value(
        "carrier_transport_couple_profile", "mesh_default");
    const bool hasFile = geometry.contains("external_averagebox_couples_file");
    const bool hasCount = geometry.contains("external_averagebox_expected_edges");
    if (report.profile == "mesh_default") {
        if (hasFile || hasCount) {
            throw std::runtime_error(
                "ConfigParsing: external AverageBox fields require "
                "mesh_geometry.carrier_transport_couple_profile="
                "'templates_ldmos_external_averagebox'.");
        }
        return report;
    }
    if (report.profile != "templates_ldmos_external_averagebox") {
        throw std::runtime_error(
            "ConfigParsing: mesh_geometry.carrier_transport_couple_profile must be "
            "'mesh_default' or the template-private diagnostic "
            "'templates_ldmos_external_averagebox'.");
    }
    if (geometry.value("node_volume_policy", "barycentric") != "barycentric") {
        throw std::runtime_error(
            "ConfigParsing: templates_ldmos_external_averagebox requires "
            "mesh_geometry.node_volume_policy='barycentric'.");
    }
    if (!hasFile || !geometry.at("external_averagebox_couples_file").is_string()) {
        throw std::runtime_error(
            "ConfigParsing: templates_ldmos_external_averagebox requires a string "
            "external_averagebox_couples_file.");
    }
    if (!hasCount ||
        (!geometry.at("external_averagebox_expected_edges").is_number_unsigned() &&
         !geometry.at("external_averagebox_expected_edges").is_number_integer())) {
        throw std::runtime_error(
            "ConfigParsing: templates_ldmos_external_averagebox requires an integer "
            "external_averagebox_expected_edges.");
    }
    const auto expectedWide =
        geometry.at("external_averagebox_expected_edges").get<std::int64_t>();
    if (expectedWide <= 0) {
        throw std::runtime_error(
            "ConfigParsing: external_averagebox_expected_edges must be positive.");
    }
    const Index expected = static_cast<Index>(expectedWide);

    std::filesystem::path source =
        geometry.at("external_averagebox_couples_file").get<std::string>();
    if (source.is_relative())
        source = configDirectory / source;
    source = std::filesystem::weakly_canonical(source);
    report.sourceFile = source.string();
    std::ifstream input(source);
    if (!input.is_open()) {
        throw std::runtime_error(
            "ConfigParsing: cannot open external AverageBox couples file: " +
            source.string());
    }

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("ConfigParsing: external AverageBox CSV is empty.");
    const auto header = splitCsvLine(line);
    if (header != std::vector<std::string>{"node0", "node1", "couple_m"}) {
        throw std::runtime_error(
            "ConfigParsing: external AverageBox CSV header must be exactly "
            "node0,node1,couple_m.");
    }

    std::map<std::pair<Index, Index>, Index> edgeByNodes;
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        edgeByNodes[{std::min(edge.n0, edge.n1), std::max(edge.n0, edge.n1)}] = edgeId;
    }
    std::map<std::pair<Index, Index>, Real> overrides;
    Index lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (trimCsvToken(line).empty())
            continue;
        const auto fields = splitCsvLine(line);
        if (fields.size() != 3) {
            throw std::runtime_error(
                "ConfigParsing: external AverageBox CSV row " +
                std::to_string(lineNumber) + " must contain three fields.");
        }
        Index node0 = 0;
        Index node1 = 0;
        Real coupleM = 0.0;
        try {
            node0 = static_cast<Index>(std::stoull(fields[0]));
            node1 = static_cast<Index>(std::stoull(fields[1]));
            coupleM = std::stod(fields[2]);
        } catch (const std::exception&) {
            throw std::runtime_error(
                "ConfigParsing: invalid external AverageBox CSV value at row " +
                std::to_string(lineNumber) + ".");
        }
        if (node0 == node1 || !std::isfinite(coupleM) || coupleM < 0.0) {
            throw std::runtime_error(
                "ConfigParsing: external AverageBox row " +
                std::to_string(lineNumber) +
                " requires distinct nodes and a finite non-negative couple_m.");
        }
        const auto key = std::minmax(node0, node1);
        const std::pair<Index, Index> pair{key.first, key.second};
        if (!edgeByNodes.contains(pair)) {
            throw std::runtime_error(
                "ConfigParsing: external AverageBox row " +
                std::to_string(lineNumber) + " references a non-mesh edge.");
        }
        const Real internal = scaling.unitSystem().metersToInternalLength(coupleM);
        if (!overrides.emplace(pair, internal).second) {
            throw std::runtime_error(
                "ConfigParsing: duplicate external AverageBox edge at row " +
                std::to_string(lineNumber) + ".");
        }
    }
    if (overrides.size() != static_cast<std::size_t>(expected)) {
        throw std::runtime_error(
            "ConfigParsing: external AverageBox edge count mismatch: expected " +
            std::to_string(expected) + ", read " +
            std::to_string(overrides.size()) + ".");
    }

    bool first = true;
    for (const auto& [nodes, couple] : overrides) {
        mesh.setTransportCouple(edgeByNodes.at(nodes), couple);
        ++report.records;
        if (couple == 0.0)
            ++report.zeroCouples;
        if (first) {
            report.minimumCouple = couple;
            report.maximumCouple = couple;
            first = false;
        } else {
            report.minimumCouple = std::min(report.minimumCouple, couple);
            report.maximumCouple = std::max(report.maximumCouple, couple);
        }
    }
    return report;
}

} // namespace vela
