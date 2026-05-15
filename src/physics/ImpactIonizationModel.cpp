#include "vela/physics/ImpactIonizationModel.h"
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
    throw std::invalid_argument(
        "makeImpactIonizationModel: unknown impact ionization model '" + config.model + "'.");
}

} // namespace vela
