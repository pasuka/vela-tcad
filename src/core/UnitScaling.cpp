#include "vela/core/UnitScaling.h"

#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

namespace vela {

PhysicalUnitSystem::PhysicalUnitSystem(Real lengthMPerInternal,
                                       Real concentrationM3PerInternal,
                                       Real sheetDensityM2PerInternal,
                                       Real mobilityM2PerVSPerInternal,
                                       Real velocityMPerSPerInternal,
                                       Real electricFieldVPerMPerInternal,
                                       Real inverseLengthMInvPerInternal,
                                       Real surfaceFieldCoefficientMPerVPerInternal,
                                       Real currentDensityAM2PerInternal)
    : lengthMPerInternal_(lengthMPerInternal),
      concentrationM3PerInternal_(concentrationM3PerInternal),
      sheetDensityM2PerInternal_(sheetDensityM2PerInternal),
      mobilityM2PerVSPerInternal_(mobilityM2PerVSPerInternal),
      velocityMPerSPerInternal_(velocityMPerSPerInternal),
      electricFieldVPerMPerInternal_(electricFieldVPerMPerInternal),
      inverseLengthMInvPerInternal_(inverseLengthMInvPerInternal),
      surfaceFieldCoefficientMPerVPerInternal_(surfaceFieldCoefficientMPerVPerInternal),
      currentDensityAM2PerInternal_(currentDensityAM2PerInternal)
{}

PhysicalUnitSystem PhysicalUnitSystem::legacySI()
{
    return PhysicalUnitSystem(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
}

PhysicalUnitSystem PhysicalUnitSystem::tcadInternal()
{
    return PhysicalUnitSystem(1.0e-6, 1.0e6, 1.0e4, 1.0e-4, 1.0e-2, 1.0e2, 1.0e2, 1.0e-2, 1.0e4);
}

Real PhysicalUnitSystem::chargeAreaFactor() const
{
    return concentrationM3PerInternal_ * areaM2PerInternal();
}

Real PhysicalUnitSystem::chargeLineFactor() const
{
    return sheetDensityM2PerInternal_ * lengthMPerInternal_;
}

Real PhysicalUnitSystem::chargeVolumeFactor() const
{
    return concentrationM3PerInternal_ * volumeM3PerInternal();
}

Real PhysicalUnitSystem::chargeSheetFactor() const
{
    return sheetDensityM2PerInternal_ * areaM2PerInternal();
}

Real PhysicalUnitSystem::fieldFromCoordinateDeltaFactor() const
{
    return (1.0 / lengthMPerInternal_) / electricFieldVPerMPerInternal_;
}

Real PhysicalUnitSystem::continuitySourceIntegralFactor() const
{
    // Match a volumetric particle source integrated over a 2-D control area
    // to the native SG edge-flux line integral.  The latter is converted to
    // physical current per device depth by
    //   currentDensityAM2PerInternal * lengthMPerInternal,
    // with q applied by the caller.  The previous area/mobility expression
    // omitted the concentration and field/current-density unit scales; in the
    // TCAD cm/um convention it understated every continuity source by 1e4.
    return concentrationM3PerInternal_ * areaM2PerInternal()
         / (currentDensityAM2PerInternal_ * lengthMPerInternal_);
}

Real PhysicalUnitSystem::currentPerInternalDepthFactor() const
{
    return 1.0e-6;
}

Real PhysicalUnitSystem::internalCurrentPerDeviceDepthToAPerUm(Real value) const
{
    return value * currentPerInternalDepthFactor();
}

UnitScalingConfig::UnitScalingConfig(UnitScalingMode modeValue)
    : mode(modeValue),
      physicalUnitSystem(modeValue == UnitScalingMode::UnitScaling
                             ? PhysicalUnitSystem::tcadInternal()
                             : PhysicalUnitSystem::legacySI())
{}

Real UnitScalingConfig::lengthToSI(Real value) const
{
    return unitSystem().internalLengthToMeters(value);
}

Real UnitScalingConfig::concentrationToSI(Real value) const
{
    return unitSystem().internalConcentrationToM3(value);
}

Real UnitScalingConfig::sheetDensityToSI(Real value) const
{
    return unitSystem().internalSheetDensityToM2(value);
}

Real UnitScalingConfig::mobilityToSI(Real value) const
{
    return unitSystem().internalMobilityToM2PerVS(value);
}

Real UnitScalingConfig::velocityToSI(Real value) const
{
    return unitSystem().internalVelocityToMPerS(value);
}

Real UnitScalingConfig::electricFieldToSI(Real value) const
{
    return unitSystem().internalElectricFieldToVPerM(value);
}

Real UnitScalingConfig::inverseLengthToSI(Real value) const
{
    return unitSystem().internalInverseLengthToMInv(value);
}

Real UnitScalingConfig::surfaceFieldCoefficientToSI(Real value) const
{
    return unitSystem().internalSurfaceFieldCoefficientToMPerV(value);
}

int parseDeckFormatVersion(const nlohmann::json& cfg)
{
    if (!cfg.is_object() || !cfg.contains("format_version"))
        return 0;

    const auto& value = cfg.at("format_version");
    if (!value.is_number_integer()) {
        throw std::invalid_argument(
            "format_version must be the integer 2.");
    }

    const int version = value.get<int>();
    if (version == 2)
        return 2;

    throw std::invalid_argument(
        "Unsupported format_version " + std::to_string(version) +
        ". The only supported deck format version is 2, in which every value is "
        "expressed in the TCAD internal units (um, cm^-3, cm^2/(V s), V/cm). "
        "Run the deck migration tool to convert an older deck.");
}

UnitScalingConfig parseUnitScalingConfig(const nlohmann::json& cfg)
{
    if (parseDeckFormatVersion(cfg) == 2) {
        if (cfg.contains("scaling")) {
            throw std::invalid_argument(
                "scaling is not accepted with format_version 2; the deck already "
                "uses the TCAD internal units. Move "
                "characteristic_length_um, reference_concentration_cm3 and "
                "reference_mobility_cm2_V_s to solver.normalization.");
        }
        return UnitScalingConfig{UnitScalingMode::UnitScaling};
    }

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
