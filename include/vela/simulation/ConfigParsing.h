#pragma once

#include "vela/equation/ChargeSpec.h"
#include <nlohmann/json_fwd.hpp>
#include <vector>

namespace vela {

std::vector<RegionFixedChargeSpec> parseRegionFixedChargeSpecs(
    const nlohmann::json& cfg);

std::vector<InterfaceSheetChargeSpec> parseInterfaceSheetChargeSpecs(
    const nlohmann::json& cfg);

} // namespace vela
