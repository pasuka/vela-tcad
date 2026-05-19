#include "vela/core/UnitScalingSystem.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"

#include <algorithm>
#include <cmath>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

namespace vela {

namespace {

Real requirePositive(Real value, const char* name)
{
    if (value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive.");
    }
    return value;
}

std::optional<Real> parseAutoOrPositiveNumber(const nlohmann::json& scaling,
                                              const char* key,
                                              Real (UnitScalingConfig::*toSI)(Real) const)
{
    if (!scaling.contains(key)) {
        return std::nullopt;
    }

    const nlohmann::json& value = scaling.at(key);
    if (value.is_string()) {
        const std::string text = value.get<std::string>();
        if (text == "auto") {
            return std::nullopt;
        }
        throw std::invalid_argument(std::string("scaling.") + key +
                                    " must be 'auto' or a positive number.");
    }

    if (!value.is_number()) {
        throw std::invalid_argument(std::string("scaling.") + key +
                                    " must be 'auto' or a positive number.");
    }

    const UnitScalingConfig unit{UnitScalingMode::UnitScaling};
    const Real numericValue = (unit.*toSI)(value.get<Real>());
    if (numericValue <= 0.0) {
        throw std::invalid_argument(std::string("scaling.") + key + " must be positive.");
    }
    return numericValue;
}

} // namespace

UnitScalingReferenceConfig parseUnitScalingReferenceConfig(const nlohmann::json& cfg)
{
    UnitScalingReferenceConfig refs;

    if (!cfg.contains("scaling")) {
        return refs;
    }

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
    if (mode != "unit_scaling") {
        throw std::invalid_argument(
            "Unsupported scaling.mode '" + mode + "'. Supported value is 'unit_scaling'. "
            "Omit scaling to keep legacy SI input behavior.");
    }

    refs.characteristicLength_m =
        parseAutoOrPositiveNumber(scaling,
                                  "characteristic_length_um",
                                  &UnitScalingConfig::lengthToSI);
    refs.referenceConcentration_m3 =
        parseAutoOrPositiveNumber(scaling,
                                  "reference_concentration_cm3",
                                  &UnitScalingConfig::concentrationToSI);
    refs.referenceMobility_m2_V_s =
        parseAutoOrPositiveNumber(scaling,
                                  "reference_mobility_cm2_V_s",
                                  &UnitScalingConfig::mobilityToSI);

    return refs;
}

UnitScalingSystem::UnitScalingSystem(Real temperature_K,
                                     Real epsRef_F_per_m,
                                     Real concentrationScale_m3,
                                     Real lengthScale_m,
                                     Real mobilityScale_m2_V_s)
{
    using namespace constants;

    const Real T = requirePositive(temperature_K, "temperature");
    const Real epsRef = requirePositive(epsRef_F_per_m, "eps_ref");
    C0_ = requirePositive(concentrationScale_m3, "reference concentration");
    L0_ = requirePositive(lengthScale_m, "characteristic length");
    mu0_ = requirePositive(mobilityScale_m2_V_s, "reference mobility");

    V0_ = kb * T / q;
    D0_ = mu0_ * V0_;
    E0_ = V0_ / L0_;
    rho0_ = q * C0_;
    lambda2_ = epsRef * V0_ / (q * C0_ * L0_ * L0_);
    J0_ = q * D0_ * C0_ / L0_;
    R0_ = D0_ * C0_ / (L0_ * L0_);
}

UnitScalingSystem UnitScalingSystem::fromInputs(Real temperature_K,
                                                Real epsRef_F_per_m,
                                                const AutoInputs& inputs,
                                                const UnitScalingReferenceConfig& refs)
{
    const Real autoC0 = std::max(std::abs(inputs.maxAbsNetDoping_m3),
                                 requirePositive(inputs.niFloor_m3, "ni_floor"));
    const Real autoL0 = requirePositive(inputs.meshMaxLength_m, "mesh max length");
    const Real autoMu0 = requirePositive(inputs.maxMobility_m2_V_s, "max mobility");

    const Real C0 = refs.referenceConcentration_m3.value_or(autoC0);
    const Real L0 = refs.characteristicLength_m.value_or(autoL0);
    const Real mu0 = refs.referenceMobility_m2_V_s.value_or(autoMu0);

    return UnitScalingSystem(temperature_K, epsRef_F_per_m, C0, L0, mu0);
}

UnitScalingSystem::AutoInputs UnitScalingSystem::autoInputsFrom(const DeviceMesh& mesh,
                                                                const DopingModel& doping,
                                                                const MaterialDatabase& materials,
                                                                Real niFloor_m3)
{
    AutoInputs inputs;

    if (mesh.numNodes() == 0) {
        throw std::invalid_argument("UnitScalingSystem requires a non-empty mesh.");
    }

    inputs.niFloor_m3 = requirePositive(niFloor_m3, "ni_floor");

    Real minX = mesh.getNode(0).x;
    Real maxX = mesh.getNode(0).x;
    Real minY = mesh.getNode(0).y;
    Real maxY = mesh.getNode(0).y;
    for (Index i = 1; i < mesh.numNodes(); ++i) {
        const Node& node = mesh.getNode(i);
        minX = std::min(minX, node.x);
        maxX = std::max(maxX, node.x);
        minY = std::min(minY, node.y);
        maxY = std::max(maxY, node.y);
    }
    inputs.meshMaxLength_m = std::max(maxX - minX, maxY - minY);

    Real maxAbsNet = 0.0;
    for (Index i = 0; i < doping.numNodes(); ++i) {
        maxAbsNet = std::max(maxAbsNet, std::abs(doping.netDoping(i)));
    }
    inputs.maxAbsNetDoping_m3 = maxAbsNet;

    Real maxMobility = 0.0;
    for (const Region& region : mesh.regions()) {
        const Material& material = materials.getMaterial(region.material);
        maxMobility = std::max(maxMobility, std::abs(material.mun));
        maxMobility = std::max(maxMobility, std::abs(material.mup));
    }
    inputs.maxMobility_m2_V_s = maxMobility;

    return inputs;
}

} // namespace vela
