#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/physics/DopingModel.h"
#include <nlohmann/json_fwd.hpp>
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

} // namespace vela
