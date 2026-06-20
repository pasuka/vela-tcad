#include "vela/simulation/ConfigParsing.h"

#include <nlohmann/json.hpp>
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
                scaling.concentrationToSI(entry.at("fixed_charge_m3").get<Real>()),
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
                scaling.concentrationToSI(entry.at("fixed_charge_m3").get<Real>()),
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
        spec.donors = scaling.concentrationToSI(entry.at("donors").get<Real>());
        spec.acceptors = scaling.concentrationToSI(entry.at("acceptors").get<Real>());
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

        const Real sheetCharge = scaling.sheetDensityToSI(
            entry.value("sheet_charge_m2", 0.0));
        const Real fixedCharge = scaling.sheetDensityToSI(
            entry.value("fixed_charge_m2", 0.0));
        const Real trapDensity = scaling.sheetDensityToSI(
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

    return options;
}

} // namespace vela
