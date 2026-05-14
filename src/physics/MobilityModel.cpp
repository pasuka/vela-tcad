#include "vela/physics/MobilityModel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

Real ConstantMobility::electronMobility(const Material& material,
                                        Real,
                                        Real,
                                        Real,
                                        Real) const
{
    return material.mun;
}

Real ConstantMobility::holeMobility(const Material& material,
                                    Real,
                                    Real,
                                    Real,
                                    Real) const
{
    return material.mup;
}

DopingDependentMobility::DopingDependentMobility(MobilityModelConfig config)
    : config_(std::move(config))
{}

Real DopingDependentMobility::electronMobility(const Material& material,
                                               Real netDoping,
                                               Real,
                                               Real,
                                               Real electricField) const
{
    const Real lowField = caugheyThomas(material.mun, netDoping, config_.electronCT);
    if (config_.model == "caughey_thomas_field")
        return fieldLimit(lowField, electricField, config_.electronField);
    return lowField;
}

Real DopingDependentMobility::holeMobility(const Material& material,
                                           Real netDoping,
                                           Real,
                                           Real,
                                           Real electricField) const
{
    const Real lowField = caugheyThomas(material.mup, netDoping, config_.holeCT);
    if (config_.model == "caughey_thomas_field")
        return fieldLimit(lowField, electricField, config_.holeField);
    return lowField;
}

Real DopingDependentMobility::caugheyThomas(
    Real muMax,
    Real netDoping,
    const CaugheyThomasParameters& params)
{
    if (muMax <= 0.0)
        return 0.0;
    if (params.nRef <= 0.0 || params.alpha <= 0.0)
        throw std::invalid_argument(
            "DopingDependentMobility: Caughey-Thomas nRef and alpha must be positive.");

    const Real muMin = std::clamp(params.muMin, 0.0, muMax);
    const Real normalizedDoping = std::abs(netDoping) / params.nRef;
    const Real rolloff = std::pow(normalizedDoping, params.alpha);
    return muMin + (muMax - muMin) / (1.0 + rolloff);
}

Real DopingDependentMobility::fieldLimit(Real lowFieldMobility,
                                         Real electricField,
                                         const FieldMobilityParameters& params)
{
    if (lowFieldMobility <= 0.0)
        return 0.0;
    if (params.saturationVelocity <= 0.0 || params.beta <= 0.0)
        throw std::invalid_argument(
            "DopingDependentMobility: field saturation velocity and beta must be positive.");
    const Real field = std::abs(electricField);
    if (field <= 0.0)
        return lowFieldMobility;
    const Real ratio = lowFieldMobility * field / params.saturationVelocity;
    return lowFieldMobility / std::pow(1.0 + std::pow(ratio, params.beta), 1.0 / params.beta);
}

MobilityModelConfig mobilityModelConfig(std::string modelName)
{
    MobilityModelConfig config;
    config.model = std::move(modelName);
    return config;
}

std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config)
{
    if (config.model == "constant")
        return std::make_unique<ConstantMobility>();
    if (config.model == "caughey_thomas" || config.model == "caughey_thomas_field")
        return std::make_unique<DopingDependentMobility>(config);

    throw std::invalid_argument(
        "makeMobilityModel: unknown mobility model '" + config.model + "'.");
}

} // namespace vela
