#pragma once

#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>

namespace vela {

enum class UnitScalingMode {
    LegacySI,
    UnitScaling,
};

class PhysicalUnitSystem {
public:
    static PhysicalUnitSystem legacySI();
    static PhysicalUnitSystem tcadInternal();

    Real lengthMPerInternal() const { return lengthMPerInternal_; }
    Real areaM2PerInternal() const { return lengthMPerInternal_ * lengthMPerInternal_; }
    Real volumeM3PerInternal() const { return areaM2PerInternal() * lengthMPerInternal_; }
    Real concentrationM3PerInternal() const { return concentrationM3PerInternal_; }
    Real sheetDensityM2PerInternal() const { return sheetDensityM2PerInternal_; }
    Real mobilityM2PerVSPerInternal() const { return mobilityM2PerVSPerInternal_; }
    Real velocityMPerSPerInternal() const { return velocityMPerSPerInternal_; }
    Real electricFieldVPerMPerInternal() const { return electricFieldVPerMPerInternal_; }
    Real inverseLengthMInvPerInternal() const { return inverseLengthMInvPerInternal_; }
    Real surfaceFieldCoefficientMPerVPerInternal() const
    {
        return surfaceFieldCoefficientMPerVPerInternal_;
    }
    Real currentDensityAM2PerInternal() const { return currentDensityAM2PerInternal_; }
    /**
     * Auger coefficients multiply a carrier density by an excess carrier
     * product, so their unit scale is the inverse square of the concentration
     * scale: SI [m^6/s] per internal [cm^6/s] is 1e-12 in the TCAD system.
     */
    Real augerCoefficientM6SPerInternal() const
    {
        return 1.0 / (concentrationM3PerInternal_ * concentrationM3PerInternal_);
    }

    /**
     * Converts a volumetric concentration integrated over a 2-D control area
     * into the SI per-device-depth measure used by the Poisson rows.
     */
    Real chargeAreaFactor() const;
    /**
     * Converts a sheet density integrated along a 2-D interface line into the
     * SI per-device-depth measure used by the Poisson rows.
     */
    Real chargeLineFactor() const;
    /// Three-dimensional volume conversion retained for genuinely 3-D uses.
    Real chargeVolumeFactor() const;
    /// Three-dimensional surface conversion retained for genuinely 3-D uses.
    Real chargeSheetFactor() const;
    Real fieldFromCoordinateDeltaFactor() const;
    /**
     * Converts a volumetric rate integrated over 2-D area into the native
     * Scharfetter-Gummel particle-flux line-integral units used by continuity
     * rows.  This includes concentration, area, current-density, and device
     * depth unit conversions; the elementary charge is applied separately.
     */
    Real continuitySourceIntegralFactor() const;
    Real currentPerInternalDepthFactor() const;

    Real internalLengthToMeters(Real value) const { return value * lengthMPerInternal_; }
    Real metersToInternalLength(Real value) const { return value / lengthMPerInternal_; }
    Real internalConcentrationToM3(Real value) const
    {
        return value * concentrationM3PerInternal_;
    }
    Real m3ToInternalConcentration(Real value) const
    {
        return value / concentrationM3PerInternal_;
    }
    Real internalSheetDensityToM2(Real value) const { return value * sheetDensityM2PerInternal_; }
    Real m2ToInternalSheetDensity(Real value) const { return value / sheetDensityM2PerInternal_; }
    Real internalMobilityToM2PerVS(Real value) const
    {
        return value * mobilityM2PerVSPerInternal_;
    }
    Real m2PerVSToInternalMobility(Real value) const
    {
        return value / mobilityM2PerVSPerInternal_;
    }
    Real internalVelocityToMPerS(Real value) const
    {
        return value * velocityMPerSPerInternal_;
    }
    Real mPerSToInternalVelocity(Real value) const
    {
        return value / velocityMPerSPerInternal_;
    }
    Real internalElectricFieldToVPerM(Real value) const
    {
        return value * electricFieldVPerMPerInternal_;
    }
    Real vPerMToInternalElectricField(Real value) const
    {
        return value / electricFieldVPerMPerInternal_;
    }
    Real internalInverseLengthToMInv(Real value) const
    {
        return value * inverseLengthMInvPerInternal_;
    }
    Real mInvToInternalInverseLength(Real value) const
    {
        return value / inverseLengthMInvPerInternal_;
    }
    Real internalSurfaceFieldCoefficientToMPerV(Real value) const
    {
        return value * surfaceFieldCoefficientMPerVPerInternal_;
    }
    Real mPerVToInternalSurfaceFieldCoefficient(Real value) const
    {
        return value / surfaceFieldCoefficientMPerVPerInternal_;
    }
    Real internalCurrentDensityToAPerM2(Real value) const
    {
        return value * currentDensityAM2PerInternal_;
    }
    /**
     * Native continuity particle flux uses the same area scale as native
     * current density because the elementary charge is already in SI units.
     */
    Real internalContinuityParticleFluxToPerM2PerS(Real value) const
    {
        return value * currentDensityAM2PerInternal_;
    }

    Real aPerM2ToInternalCurrentDensity(Real value) const
    {
        return value / currentDensityAM2PerInternal_;
    }
    Real internalAugerCoefficientToM6PerS(Real value) const
    {
        return value * augerCoefficientM6SPerInternal();
    }
    Real m6PerSToInternalAugerCoefficient(Real value) const
    {
        return value / augerCoefficientM6SPerInternal();
    }
    Real internalCurrentPerDeviceDepthToAPerUm(Real value) const;

private:
    PhysicalUnitSystem(Real lengthMPerInternal,
                       Real concentrationM3PerInternal,
                       Real sheetDensityM2PerInternal,
                       Real mobilityM2PerVSPerInternal,
                       Real velocityMPerSPerInternal,
                       Real electricFieldVPerMPerInternal,
                       Real inverseLengthMInvPerInternal,
                       Real surfaceFieldCoefficientMPerVPerInternal,
                       Real currentDensityAM2PerInternal);

    Real lengthMPerInternal_ = 1.0;
    Real concentrationM3PerInternal_ = 1.0;
    Real sheetDensityM2PerInternal_ = 1.0;
    Real mobilityM2PerVSPerInternal_ = 1.0;
    Real velocityMPerSPerInternal_ = 1.0;
    Real electricFieldVPerMPerInternal_ = 1.0;
    Real inverseLengthMInvPerInternal_ = 1.0;
    Real surfaceFieldCoefficientMPerVPerInternal_ = 1.0;
    Real currentDensityAM2PerInternal_ = 1.0;
};

struct UnitScalingConfig {
    UnitScalingMode mode = UnitScalingMode::LegacySI;
    PhysicalUnitSystem physicalUnitSystem = PhysicalUnitSystem::legacySI();

    UnitScalingConfig() = default;
    explicit UnitScalingConfig(UnitScalingMode mode);

    bool isUnitScaling() const { return mode == UnitScalingMode::UnitScaling; }
    const PhysicalUnitSystem& unitSystem() const { return physicalUnitSystem; }

    Real lengthToInternal(Real value) const { return value; }
    Real concentrationToInternal(Real value) const { return value; }
    Real sheetDensityToInternal(Real value) const { return value; }
    Real mobilityToInternal(Real value) const { return value; }
    Real velocityToInternal(Real value) const { return value; }
    Real electricFieldToInternal(Real value) const { return value; }
    Real inverseLengthToInternal(Real value) const { return value; }
    Real surfaceFieldCoefficientToInternal(Real value) const { return value; }
    Real augerCoefficientToInternal(Real value) const { return value; }

    Real lengthToSI(Real value) const;
    Real concentrationToSI(Real value) const;
    Real sheetDensityToSI(Real value) const;
    Real mobilityToSI(Real value) const;
    Real velocityToSI(Real value) const;
    Real electricFieldToSI(Real value) const;
    Real inverseLengthToSI(Real value) const;
    Real surfaceFieldCoefficientToSI(Real value) const;
};

// Returns the declared deck format version, or 0 when `format_version` is
// absent. Version 2 is the single-unit-system deck: values are expressed in
// the TCAD internal units and the `scaling` block is not accepted. Version 0
// is the transitional legacy form that still selects its unit system through
// `scaling.mode`.
int parseDeckFormatVersion(const nlohmann::json& cfg);

UnitScalingConfig parseUnitScalingConfig(const nlohmann::json& cfg);

} // namespace vela
