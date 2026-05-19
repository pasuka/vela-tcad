#pragma once

#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>

namespace vela {

enum class UnitScalingMode {
    LegacySI,
    UnitScaling,
};

struct UnitScalingConfig {
    UnitScalingMode mode = UnitScalingMode::LegacySI;

    bool isUnitScaling() const { return mode == UnitScalingMode::UnitScaling; }

    Real lengthToSI(Real value) const;
    Real concentrationToSI(Real value) const;
    Real sheetDensityToSI(Real value) const;
    Real mobilityToSI(Real value) const;
    Real electricFieldToSI(Real value) const;
    Real inverseLengthToSI(Real value) const;
    Real surfaceFieldCoefficientToSI(Real value) const;
};

UnitScalingConfig parseUnitScalingConfig(const nlohmann::json& cfg);

} // namespace vela
