#include "vela/core/UnitScaling.h"

#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

namespace vela {

Real UnitScalingConfig::lengthToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e-6 : value;
}

Real UnitScalingConfig::concentrationToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e6 : value;
}

Real UnitScalingConfig::sheetDensityToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e4 : value;
}

Real UnitScalingConfig::mobilityToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e-4 : value;
}

Real UnitScalingConfig::electricFieldToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e2 : value;
}

Real UnitScalingConfig::inverseLengthToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e2 : value;
}

Real UnitScalingConfig::surfaceFieldCoefficientToSI(Real value) const
{
    return isUnitScaling() ? value * 1.0e-2 : value;
}

UnitScalingConfig parseUnitScalingConfig(const nlohmann::json& cfg)
{
    if (!cfg.contains("scaling"))
        return {};

    const auto& scaling = cfg.at("scaling");
    if (!scaling.is_object()) {
        throw std::invalid_argument(
            "scaling must be an object with mode set to 'unit_scaling'.");
    }
    if (scaling.contains("system")) {
        throw std::invalid_argument(
            "scaling.system is not supported; use scaling.mode = 'unit_scaling'.");
    }
    if (!scaling.contains("mode")) {
        throw std::invalid_argument(
            "scaling.mode is required when scaling is present; supported value is 'unit_scaling'.");
    }

    const std::string mode = scaling.at("mode").get<std::string>();
    if (mode == "unit_scaling")
        return UnitScalingConfig{UnitScalingMode::UnitScaling};

    throw std::invalid_argument(
        "Unsupported scaling.mode '" + mode + "'. Supported value is 'unit_scaling'. "
        "Omit scaling to keep legacy SI input behavior.");
}

} // namespace vela
