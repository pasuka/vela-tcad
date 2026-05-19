#pragma once

#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>
#include <optional>

namespace vela {

class DeviceMesh;
class DopingModel;
class MaterialDatabase;

struct UnitScalingReferenceConfig {
    // Nullopt means "auto".
    std::optional<Real> characteristicLength_m;
    std::optional<Real> referenceConcentration_m3;
    std::optional<Real> referenceMobility_m2_V_s;
};

UnitScalingReferenceConfig parseUnitScalingReferenceConfig(const nlohmann::json& cfg);

class UnitScalingSystem {
public:
    struct AutoInputs {
        Real maxAbsNetDoping_m3 = 0.0;
        Real niFloor_m3 = 0.0;
        Real meshMaxLength_m = 0.0;
        Real maxMobility_m2_V_s = 0.0;
    };

    UnitScalingSystem(Real temperature_K,
                      Real epsRef_F_per_m,
                      Real concentrationScale_m3,
                      Real lengthScale_m,
                      Real mobilityScale_m2_V_s);

    static UnitScalingSystem fromInputs(Real temperature_K,
                                        Real epsRef_F_per_m,
                                        const AutoInputs& inputs,
                                        const UnitScalingReferenceConfig& refs = {});

    static AutoInputs autoInputsFrom(const DeviceMesh& mesh,
                                     const DopingModel& doping,
                                     const MaterialDatabase& materials,
                                     Real niFloor_m3);

    Real V0() const { return V0_; }
    Real C0() const { return C0_; }
    Real L0() const { return L0_; }
    Real mu0() const { return mu0_; }
    Real D0() const { return D0_; }
    Real lambda2() const { return lambda2_; }
    Real J0() const { return J0_; }
    Real R0() const { return R0_; }
    Real E0() const { return E0_; }
    Real rho0() const { return rho0_; }

    Real scalePotential(Real value) const { return value / V0_; }
    Real scaleLength(Real value) const { return value / L0_; }
    Real scaleConcentration(Real value) const { return value / C0_; }
    Real scaleElectricField(Real value) const { return value / E0_; }
    Real scaleCurrentDensity(Real value) const { return value / J0_; }

    Real unscalePotential(Real value) const { return value * V0_; }
    Real unscaleLength(Real value) const { return value * L0_; }
    Real unscaleConcentration(Real value) const { return value * C0_; }
    Real unscaleElectricField(Real value) const { return value * E0_; }
    Real unscaleCurrentDensity(Real value) const { return value * J0_; }

private:
    Real V0_ = 0.0;
    Real C0_ = 0.0;
    Real L0_ = 0.0;
    Real mu0_ = 0.0;
    Real D0_ = 0.0;
    Real lambda2_ = 0.0;
    Real J0_ = 0.0;
    Real R0_ = 0.0;
    Real E0_ = 0.0;
    Real rho0_ = 0.0;
};

} // namespace vela
