#include "vela/physics/ImpactIonizationModel.h"
#include "vela/core/PhysicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

Real NoImpactIonization::electronCoefficient(Real) const { return 0.0; }
Real NoImpactIonization::holeCoefficient(Real) const { return 0.0; }
Real NoImpactIonization::generationRate(Real, Real, Real) const { return 0.0; }

SelberherrImpactIonization::SelberherrImpactIonization(ImpactIonizationModelConfig config)
    : config_(std::move(config))
{
    if (config_.electronA < 0.0 || config_.holeA < 0.0 ||
        config_.electronB <= 0.0 || config_.holeB <= 0.0 ||
        config_.carrierVelocity < 0.0) {
        throw std::invalid_argument(
            "SelberherrImpactIonization: prefactors/velocity must be non-negative and critical fields positive.");
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
    return coefficient(electricField, config_.electronA, config_.electronB);
}

Real SelberherrImpactIonization::holeCoefficient(Real electricField) const
{
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
    : config_(std::move(config))
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

std::unique_ptr<ImpactIonizationModel> makeImpactIonizationModel(
    const ImpactIonizationModelConfig& config)
{
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
