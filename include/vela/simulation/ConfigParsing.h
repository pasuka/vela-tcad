#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/mesh/BoxGeometryBuilder.h"
#include "vela/physics/DopingModel.h"
#include <nlohmann/json_fwd.hpp>
#include <filesystem>
#include <string>
#include <vector>

namespace vela {

std::vector<RegionDopingSpec> parseDopingSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

std::vector<RegionFixedChargeSpec> parseRegionFixedChargeSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

std::vector<InterfaceSheetChargeSpec> parseInterfaceSheetChargeSpecs(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

BoxGeometryBuilder::Options parseBoxGeometryOptions(const nlohmann::json& cfg);

struct CarrierTransportCoupleProfileReport {
    std::string profile = "mesh_default";
    std::string sourceFile;
    Index records = 0;
    Index zeroCouples = 0;
    Real minimumCouple = 0.0;
    Real maximumCouple = 0.0;
};

/// Apply an explicit carrier-transport-only edge-coupling profile. The
/// default profile performs no I/O and leaves every edge on mesh geometry.
CarrierTransportCoupleProfileReport applyCarrierTransportCoupleProfile(
    DeviceMesh& mesh,
    const nlohmann::json& cfg,
    const std::filesystem::path& configDirectory,
    UnitScalingConfig scaling = {});

struct PoissonCoupleProfileReport {
    std::string profile = "mesh_default";
    std::string sourceFile;
    Index records = 0;
    Index zeroCouples = 0;
    Real minimumCouple = 0.0;
    Real maximumCouple = 0.0;
};

/// Apply the explicit, Templates/LDMOS-only region-local AverageBox
/// electrostatic couple oracle. The default leaves Poisson geometry intact.
PoissonCoupleProfileReport applyPoissonCoupleProfile(
    DeviceMesh& mesh,
    const nlohmann::json& cfg,
    const std::filesystem::path& configDirectory,
    UnitScalingConfig scaling = {});

} // namespace vela
