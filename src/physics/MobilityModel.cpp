#include "vela/physics/MobilityModel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

Real ConstantMobility::electronMobility(const Material& material,
                                        Real,
                                        Real,
                                        Real) const
{
    return material.mun;
}

Real ConstantMobility::holeMobility(const Material& material,
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
                                               Real) const
{
    return caugheyThomas(material.mun, netDoping, config_.electronCT);
}

Real DopingDependentMobility::holeMobility(const Material& material,
                                           Real netDoping,
                                           Real,
                                           Real) const
{
    return caugheyThomas(material.mup, netDoping, config_.holeCT);
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
    if (config.model == "caughey_thomas")
        return std::make_unique<DopingDependentMobility>(config);

    throw std::invalid_argument(
        "makeMobilityModel: unknown mobility model '" + config.model + "'.");
}

} // namespace vela
