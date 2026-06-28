#include "vela/physics/ImpactIonizationModel.h"
#include "vela/core/PhysicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

namespace {

constexpr Real electronFitLowA_cm_inv = 2.35990376332e7;
constexpr Real electronFitLowB_V_per_cm = 6.68288073314e5;
constexpr Real electronFitHighA_cm_inv = 6.78391642452e7;
constexpr Real electronFitHighB_V_per_cm = 1.21718982697e6;
constexpr Real holeFitLowA_cm_inv = 3.90747663281e7;
constexpr Real holeFitLowB_V_per_cm = 1.10514810627e6;
constexpr Real holeFitHighA_cm_inv = 1.41230834668e8;
constexpr Real holeFitHighB_V_per_cm = 1.99067614831e6;
constexpr Real sentaurusFitSwitchField_V_per_cm = 2.5e5;

Real inverseCmToInverseM(Real value) { return value * 100.0; }
Real fieldVPerCmToVPerM(Real value) { return value * 100.0; }

void applySentaurusFitA(ImpactIonizationModelConfig& config)
{
    config.electronALow = inverseCmToInverseM(electronFitLowA_cm_inv);
    config.electronAHigh = inverseCmToInverseM(electronFitHighA_cm_inv);
    config.holeALow = inverseCmToInverseM(holeFitLowA_cm_inv);
    config.holeAHigh = inverseCmToInverseM(holeFitHighA_cm_inv);
}

void applySentaurusFitB(ImpactIonizationModelConfig& config)
{
    config.electronBLow = fieldVPerCmToVPerM(electronFitLowB_V_per_cm);
    config.electronBHigh = fieldVPerCmToVPerM(electronFitHighB_V_per_cm);
    config.holeBLow = fieldVPerCmToVPerM(holeFitLowB_V_per_cm);
    config.holeBHigh = fieldVPerCmToVPerM(holeFitHighB_V_per_cm);
}

void applyVanOverstraetenAScale(ImpactIonizationModelConfig& config)
{
    if (!std::isfinite(config.aScale) || config.aScale <= 0.0) {
        throw std::invalid_argument(
            "ImpactIonizationModelConfig: A_scale must be positive and finite.");
    }
    config.electronALow *= config.aScale;
    config.electronAHigh *= config.aScale;
    config.holeALow *= config.aScale;
    config.holeAHigh *= config.aScale;
}

} // namespace

Real NoImpactIonization::electronCoefficient(Real) const { return 0.0; }
Real NoImpactIonization::holeCoefficient(Real) const { return 0.0; }
Real NoImpactIonization::generationRate(Real, Real, Real) const { return 0.0; }

SelberherrImpactIonization::SelberherrImpactIonization(ImpactIonizationModelConfig config)
    : config_(std::move(config))
{
    if (config_.electronA < 0.0 || config_.holeA < 0.0 ||
        config_.electronB <= 0.0 || config_.holeB <= 0.0 ||
        config_.carrierVelocity < 0.0 || config_.minimumField < 0.0) {
        throw std::invalid_argument(
            "SelberherrImpactIonization: prefactors/velocity/minimum field must be non-negative and critical fields positive.");
    }
}

Real SelberherrImpactIonization::coefficient(Real electricField,
                                             Real prefactor,
                                             Real criticalField)
{
    const Real field = std::abs(electricField);
    if (field <= 0.0 || prefactor <= 0.0)
        return 0.0;
    const Real exponent = std::clamp(-criticalField / field, -700.0, 0.0);
    return prefactor * std::exp(exponent);
}

Real SelberherrImpactIonization::electronCoefficient(Real electricField) const
{
    if (std::abs(electricField) < config_.minimumField)
        return 0.0;
    return coefficient(electricField, config_.electronA, config_.electronB);
}

Real SelberherrImpactIonization::holeCoefficient(Real electricField) const
{
    if (std::abs(electricField) < config_.minimumField)
        return 0.0;
    return coefficient(electricField, config_.holeA, config_.holeB);
}

Real SelberherrImpactIonization::generationRate(Real electricField, Real n, Real p) const
{
    if (config_.carrierVelocity <= 0.0)
        return 0.0;
    return config_.carrierVelocity *
           (electronCoefficient(electricField) * std::max(n, 0.0) +
            holeCoefficient(electricField) * std::max(p, 0.0));
}

VanOverstraetenImpactIonization::VanOverstraetenImpactIonization(
    ImpactIonizationModelConfig config)
    : config_(applyImpactIonizationParameterSet(std::move(config)))
{
    const Real prefactors[] = {
        config_.electronALow,
        config_.electronAHigh,
        config_.holeALow,
        config_.holeAHigh,
    };
    for (Real prefactor : prefactors) {
        if (prefactor < 0.0)
            throw std::invalid_argument(
                "VanOverstraetenImpactIonization: prefactors must be non-negative.");
    }
    const Real criticalFields[] = {
        config_.electronBLow,
        config_.electronBHigh,
        config_.holeBLow,
        config_.holeBHigh,
        config_.switchField,
    };
    if (config_.minimumField < 0.0)
        throw std::invalid_argument(
            "VanOverstraetenImpactIonization: minimum field must be non-negative.");
    for (Real field : criticalFields) {
        if (field <= 0.0)
            throw std::invalid_argument(
                "VanOverstraetenImpactIonization: critical/switch fields must be positive.");
    }
    if (config_.carrierVelocity < 0.0 || config_.phononEnergy <= 0.0 ||
        config_.referenceTemperature_K <= 0.0 || config_.temperature_K <= 0.0) {
        throw std::invalid_argument(
            "VanOverstraetenImpactIonization: velocity, phonon energy, and temperatures must be valid.");
    }
}

Real VanOverstraetenImpactIonization::gamma() const
{
    constexpr Real kBoltzmann_eV_per_K = constants::kb / constants::q;
    const Real refArg =
        config_.phononEnergy / (2.0 * kBoltzmann_eV_per_K * config_.referenceTemperature_K);
    const Real arg =
        config_.phononEnergy / (2.0 * kBoltzmann_eV_per_K * config_.temperature_K);
    const Real denominator = std::tanh(arg);
    if (std::abs(denominator) <= 0.0)
        return 1.0;
    return std::tanh(refArg) / denominator;
}

Real VanOverstraetenImpactIonization::coefficient(Real electricField,
                                                  Real switchField,
                                                  Real lowPrefactor,
                                                  Real highPrefactor,
                                                  Real lowCriticalField,
                                                  Real highCriticalField,
                                                  Real gamma)
{
    const Real field = std::abs(electricField);
    if (field <= 0.0 || gamma <= 0.0)
        return 0.0;
    const bool lowField = field < switchField;
    const Real prefactor = lowField ? lowPrefactor : highPrefactor;
    const Real criticalField = lowField ? lowCriticalField : highCriticalField;
    if (prefactor <= 0.0)
        return 0.0;
    const Real exponent = std::clamp(-criticalField * gamma / field, -700.0, 0.0);
    return gamma * prefactor * std::exp(exponent);
}

Real VanOverstraetenImpactIonization::electronCoefficient(Real electricField) const
{
    if (!config_.debugRawVanOverstraeten && std::abs(electricField) < config_.minimumField)
        return 0.0;
    return coefficient(
        electricField,
        config_.switchField,
        config_.electronALow,
        config_.electronAHigh,
        config_.electronBLow,
        config_.electronBHigh,
        gamma());
}

Real VanOverstraetenImpactIonization::holeCoefficient(Real electricField) const
{
    if (!config_.debugRawVanOverstraeten && std::abs(electricField) < config_.minimumField)
        return 0.0;
    return coefficient(
        electricField,
        config_.switchField,
        config_.holeALow,
        config_.holeAHigh,
        config_.holeBLow,
        config_.holeBHigh,
        gamma());
}

Real VanOverstraetenImpactIonization::generationRate(Real electricField, Real n, Real p) const
{
    if (config_.carrierVelocity <= 0.0)
        return 0.0;
    return config_.carrierVelocity *
           (electronCoefficient(electricField) * std::max(n, 0.0) +
            holeCoefficient(electricField) * std::max(p, 0.0));
}

ImpactIonizationModelConfig impactIonizationModelConfig(std::string modelName)
{
    ImpactIonizationModelConfig config;
    config.model = std::move(modelName);
    return config;
}

ImpactIonizationModelConfig applyImpactIonizationParameterSet(
    ImpactIonizationModelConfig config)
{
    if (config.parameterSet == "default") {
        applyVanOverstraetenAScale(config);
        return config;
    }
    if (config.model != "van_overstraeten") {
        throw std::invalid_argument(
            "ImpactIonizationModelConfig: parameter_set requires model 'van_overstraeten'.");
    }
    if (config.parameterSet == "sentaurus_fit_A_only") {
        applySentaurusFitA(config);
        applyVanOverstraetenAScale(config);
        return config;
    }
    if (config.parameterSet == "sentaurus_fit_A_B") {
        applySentaurusFitA(config);
        applySentaurusFitB(config);
        applyVanOverstraetenAScale(config);
        return config;
    }
    if (config.parameterSet == "sentaurus_fit_A_B_switch") {
        applySentaurusFitA(config);
        applySentaurusFitB(config);
        config.switchField = fieldVPerCmToVPerM(sentaurusFitSwitchField_V_per_cm);
        applyVanOverstraetenAScale(config);
        return config;
    }
    throw std::invalid_argument(
        "ImpactIonizationModelConfig: unsupported parameter_set '" +
        config.parameterSet + "'.");
}

std::unique_ptr<ImpactIonizationModel> makeImpactIonizationModel(
    const ImpactIonizationModelConfig& config)
{
    if (config.parameterSet != "default" && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            "makeImpactIonizationModel: non-default parameter_set requires model 'van_overstraeten'.");
    }
    if (config.model == "none")
        return std::make_unique<NoImpactIonization>();
    if (config.model == "selberherr")
        return std::make_unique<SelberherrImpactIonization>(config);
    if (config.model == "van_overstraeten")
        return std::make_unique<VanOverstraetenImpactIonization>(config);
    throw std::invalid_argument(
        "makeImpactIonizationModel: unknown impact ionization model '" + config.model + "'.");
}

} // namespace vela
