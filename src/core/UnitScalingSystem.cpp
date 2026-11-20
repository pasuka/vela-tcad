#include "vela/core/UnitScalingSystem.h"
#include "vela/core/ArithmeticStageInput.h"

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
#include <type_traits>

namespace vela {

namespace {

Real requirePositive(Real value, const char* name)
{
    if (value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive.");
    }
    return value;
}

std::optional<Real> parseAutoOrPositiveNumber(const nlohmann::json& block,
                                              const std::string& blockName,
                                              const char* key,
                                              Real (UnitScalingConfig::*toInternal)(Real) const)
{
    if (!block.contains(key)) {
        return std::nullopt;
    }

    const nlohmann::json& value = block.at(key);
    if (value.is_string()) {
        const std::string text = value.get<std::string>();
        if (text == "auto") {
            return std::nullopt;
        }
        throw std::invalid_argument(blockName + "." + key +
                                    " must be 'auto' or a positive number.");
    }

    if (!value.is_number()) {
        throw std::invalid_argument(blockName + "." + key +
                                    " must be 'auto' or a positive number.");
    }

    const UnitScalingConfig unit{UnitScalingMode::UnitScaling};
    const Real numericValue = (unit.*toInternal)(value.get<Real>());
    if (numericValue <= 0.0) {
        throw std::invalid_argument(blockName + "." + key + " must be positive.");
    }
    return numericValue;
}

UnitScalingReferenceConfig parseReferenceBlock(const nlohmann::json& block,
                                               const std::string& blockName,
                                               UnitScalingReferenceConfig refs)
{
    refs.characteristicLength_m =
        parseAutoOrPositiveNumber(block, blockName,
                                  "characteristic_length_um",
                                  &UnitScalingConfig::lengthToInternal);
    refs.referenceConcentration_m3 =
        parseAutoOrPositiveNumber(block, blockName,
                                  "reference_concentration_cm3",
                                  &UnitScalingConfig::concentrationToInternal);
    refs.referenceMobility_m2_V_s =
        parseAutoOrPositiveNumber(block, blockName,
                                  "reference_mobility_cm2_V_s",
                                  &UnitScalingConfig::mobilityToInternal);
    return refs;
}

ArithmeticReal boundTemperature(Real compatibility, const UnitScalingReferenceConfig& refs)
{
    const auto& source = *refs.decimalTemperature_K;
    if (!std::isfinite(source.kelvin) || source.kelvin <= 0 ||
        source.compatibilityKelvin != compatibility)
        throw std::invalid_argument("Decimal temperature differs from the active solver temperature");
    return source.kelvin;
}

template<class Scalar>
UnitScalingValues<Scalar> deriveScales(
    Scalar temperature, Scalar epsRef, Scalar concentration, Scalar length, Scalar mobility,
    const PhysicalUnitSystem& units)
{
    Scalar lengthFactor = units.lengthMPerInternal();
    Scalar mobilityFactor = units.mobilityM2PerVSPerInternal();
    if constexpr (std::is_same_v<Scalar, long double>) {
        if (units.lengthMPerInternal() == 1.0e-6) {
            lengthFactor = 1.0e-6L;
            mobilityFactor = 1.0e-4L;
        }
    }
    const auto q = constants::arithmeticConstants<Scalar>().q;
    const Scalar concentrationSI = concentration * units.concentrationM3PerInternal();
    const Scalar lengthSI = length * lengthFactor;
    const Scalar mobilitySI = mobility * mobilityFactor;
    const Scalar v = constants::thermalVoltage(temperature);
    const Scalar dSI = mobilitySI * v;
    const Scalar rhoSI = q * concentrationSI;
    const Scalar currentSI = q * dSI * concentrationSI / lengthSI;
    const Scalar rateSI = dSI * concentrationSI / (lengthSI * lengthSI);
    return {v, concentration, length, mobility, dSI / mobilityFactor,
        epsRef * v / (q * concentrationSI * lengthSI * lengthSI),
        currentSI / units.currentDensityAM2PerInternal(),
        rateSI / units.concentrationM3PerInternal(),
        (v / lengthSI) / units.electricFieldVPerMPerInternal(),
        rhoSI / units.concentrationM3PerInternal()};
}

} // namespace

UnitScalingReferenceConfig parseUnitScalingReferenceConfig(const nlohmann::json& cfg)
{
    UnitScalingReferenceConfig refs;
    if (cfg.contains("decimal_stage_controls") && cfg.at("decimal_stage_controls").contains("solver")) {
        const auto compatibility = cfg.at("solver").value("temperature_K", constants::T0);
        refs.decimalTemperature_K = UnitScalingReferenceConfig::TemperatureInput{
            stageControl(cfg, "solver", "temperature_K", compatibility), compatibility
        };
        (void)boundTemperature(compatibility, refs);
    }

    const bool isVersion2 = parseDeckFormatVersion(cfg) == 2;
    const bool hasNormalization = cfg.contains("solver") && cfg.at("solver").is_object() &&
        cfg.at("solver").contains("normalization");

    if (isVersion2) {
        if (cfg.contains("scaling")) {
            throw std::invalid_argument(
                "scaling is not accepted with format_version 2; move "
                "characteristic_length_um, reference_concentration_cm3 and "
                "reference_mobility_cm2_V_s to solver.normalization.");
        }
        if (!hasNormalization) {
            return refs;
        }
        const auto& normalization = cfg.at("solver").at("normalization");
        if (!normalization.is_object()) {
            throw std::invalid_argument("solver.normalization must be an object.");
        }
        return parseReferenceBlock(normalization, "solver.normalization", refs);
    }

    if (hasNormalization) {
        throw std::invalid_argument(
            "solver.normalization requires format_version 2.");
    }

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

    return parseReferenceBlock(scaling, "scaling", refs);
}

Real thermalVoltageForInputs(Real temperature_K, const UnitScalingReferenceConfig& refs)
{
    if (refs.decimalTemperature_K)
        return static_cast<Real>(arithmeticThermalVoltageForInputs(temperature_K, refs));
    return constants::thermalVoltage(requirePositive(temperature_K, "temperature"));
}

ArithmeticReal arithmeticThermalVoltageForInputs(
    Real temperature_K, const UnitScalingReferenceConfig& refs)
{
    const auto temperature = refs.decimalTemperature_K
        ? boundTemperature(temperature_K, refs)
        : static_cast<ArithmeticReal>(requirePositive(temperature_K, "temperature"));
    return constants::thermalVoltage(temperature);
}

UnitScalingSystem::UnitScalingSystem(Real temperature_K,
                                     Real epsRef_F_per_m,
                                     Real concentrationScale_m3,
                                     Real lengthScale_m,
                                     Real mobilityScale_m2_V_s,
                                     PhysicalUnitSystem unitSystem)
    : UnitScalingSystem(temperature_K, epsRef_F_per_m, concentrationScale_m3, lengthScale_m,
                        mobilityScale_m2_V_s, unitSystem, {})
{}

UnitScalingSystem::UnitScalingSystem(Real temperature_K,
                                     Real epsRef_F_per_m,
                                     Real concentrationScale_m3,
                                     Real lengthScale_m,
                                     Real mobilityScale_m2_V_s,
                                     const PhysicalUnitSystem& unitSystem,
                                     const UnitScalingReferenceConfig& refs)
{
    const Real T = requirePositive(temperature_K, "temperature");
    const Real epsRef = requirePositive(epsRef_F_per_m, "eps_ref");
    C0_ = requirePositive(concentrationScale_m3, "reference concentration");
    L0_ = requirePositive(lengthScale_m, "characteristic length");
    mu0_ = requirePositive(mobilityScale_m2_V_s, "reference mobility");
    const auto wideTemperature = refs.decimalTemperature_K
        ? boundTemperature(T, refs) : static_cast<ArithmeticReal>(T);
    arithmetic_ = deriveScales(wideTemperature, static_cast<ArithmeticReal>(epsRef),
        static_cast<ArithmeticReal>(C0_), static_cast<ArithmeticReal>(L0_),
        static_cast<ArithmeticReal>(mu0_), unitSystem);
    const auto project = [this](const auto& values) {
        V0_ = static_cast<Real>(values.V0);
        D0_ = static_cast<Real>(values.D0);
        E0_ = static_cast<Real>(values.E0);
        rho0_ = static_cast<Real>(values.rho0);
        lambda2_ = static_cast<Real>(values.lambda2);
        J0_ = static_cast<Real>(values.J0);
        R0_ = static_cast<Real>(values.R0);
    };
    if (refs.decimalTemperature_K || std::is_same_v<ArithmeticReal, Real>)
        project(arithmetic_);
    else
        project(deriveScales(T, epsRef, C0_, L0_, mu0_, unitSystem));
}
UnitScalingSystem UnitScalingSystem::fromInputs(Real temperature_K,
                                                Real epsRef_F_per_m,
                                                const AutoInputs& inputs,
                                                const UnitScalingReferenceConfig& refs,
                                                PhysicalUnitSystem unitSystem){
    const Real autoC0 = std::max(std::abs(inputs.maxAbsNetDoping_m3),
                                 requirePositive(inputs.niFloor_m3, "ni_floor"));
    const Real autoL0 = requirePositive(inputs.meshMaxLength_m, "mesh max length");
    const Real autoMu0 = requirePositive(inputs.maxMobility_m2_V_s, "max mobility");

    const Real C0 = refs.referenceConcentration_m3.value_or(autoC0);
    const Real L0 = refs.characteristicLength_m.value_or(autoL0);
    const Real mu0 = refs.referenceMobility_m2_V_s.value_or(autoMu0);

    UnitScalingSystem system(temperature_K, epsRef_F_per_m, C0, L0, mu0, unitSystem, refs);
    if constexpr (compiledArithmeticPrecision == ArithmeticPrecision::Extended80) {
        const auto wideConcentration = refs.referenceConcentration_m3 || !inputs.arithmeticMaxAbsNetDoping
            ? static_cast<ArithmeticReal>(C0)
            : std::max(*inputs.arithmeticMaxAbsNetDoping, static_cast<ArithmeticReal>(inputs.niFloor_m3));
        const auto wideLength = refs.characteristicLength_m || !inputs.arithmeticMeshMaxLength
            ? static_cast<ArithmeticReal>(L0) : *inputs.arithmeticMeshMaxLength;
        const auto wideTemperature = refs.decimalTemperature_K
            ? boundTemperature(temperature_K, refs) : static_cast<ArithmeticReal>(temperature_K);
        system.arithmetic_ = deriveScales(wideTemperature, static_cast<ArithmeticReal>(epsRef_F_per_m),
            wideConcentration, wideLength, static_cast<ArithmeticReal>(mu0), unitSystem);
    }
    return system;
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
    if constexpr (compiledArithmeticPrecision == ArithmeticPrecision::Extended80) {
        auto low = mesh.arithmeticCoordinates(0);
        auto high = low;
        for (Index i = 1; i < mesh.numNodes(); ++i) {
            const auto point = mesh.arithmeticCoordinates(i);
            for (std::size_t axis = 0; axis < point.size(); ++axis) {
                low[axis] = std::min(low[axis], point[axis]);
                high[axis] = std::max(high[axis], point[axis]);
            }
        }
        inputs.arithmeticMeshMaxLength = std::max(high[0] - low[0], high[1] - low[1]);
    }

    Real maxAbsNet = 0.0;
    for (Index i = 0; i < doping.numNodes(); ++i) {
        maxAbsNet = std::max(maxAbsNet, std::abs(doping.netDoping(i)));
    }
    inputs.maxAbsNetDoping_m3 = maxAbsNet;
    if constexpr (compiledArithmeticPrecision == ArithmeticPrecision::Extended80) {
        ArithmeticReal wideNet = 0;
        for (Index i = 0; i < doping.numNodes(); ++i)
            wideNet = std::max(wideNet, std::abs(doping.netDopingForAssembly(i)));
        inputs.arithmeticMaxAbsNetDoping = wideNet;
    }

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
