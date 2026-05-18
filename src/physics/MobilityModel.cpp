#include "vela/physics/MobilityModel.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace vela {

namespace {

void parseCaugheyThomas(const nlohmann::json& json,
                        CaugheyThomasParameters& params,
                        const char* prefix)
{
    params.muMin = json.value((std::string(prefix) + "_mu_min_m2_V_s").c_str(), params.muMin);
    params.nRef = json.value((std::string(prefix) + "_nref_m3").c_str(), params.nRef);
    params.alpha = json.value((std::string(prefix) + "_alpha").c_str(), params.alpha);
}

void parseField(const nlohmann::json& json,
                FieldMobilityParameters& params,
                const char* prefix)
{
    params.saturationVelocity = json.value(
        (std::string(prefix) + "_saturation_velocity_m_s").c_str(),
        params.saturationVelocity);
    params.beta = json.value((std::string(prefix) + "_field_beta").c_str(), params.beta);
}

} // namespace

Real ConstantMobility::electronMobility(const Material& material,
                                        Real,
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
                                               Real electricField,
                                               Real surfaceNormalField) const
{
    Real mobility = caugheyThomas(material.mun, netDoping, config_.electronCT);
    if (config_.model == "caughey_thomas_field" ||
        config_.model == "caughey_thomas_field_surface")
        mobility = fieldLimit(mobility, electricField, config_.electronField);
    if (isSurfaceMobilityModel(config_))
        mobility = surfaceLimit(
            mobility, surfaceNormalField, config_.surface.thetaElectron, config_.surface);
    return mobility;
}

Real DopingDependentMobility::holeMobility(const Material& material,
                                           Real netDoping,
                                           Real,
                                           Real,
                                           Real electricField,
                                           Real surfaceNormalField) const
{
    Real mobility = caugheyThomas(material.mup, netDoping, config_.holeCT);
    if (config_.model == "caughey_thomas_field" ||
        config_.model == "caughey_thomas_field_surface")
        mobility = fieldLimit(mobility, electricField, config_.holeField);
    if (isSurfaceMobilityModel(config_))
        mobility = surfaceLimit(
            mobility, surfaceNormalField, config_.surface.thetaHole, config_.surface);
    return mobility;
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

Real DopingDependentMobility::surfaceLimit(Real bulkMobility,
                                           Real surfaceNormalField,
                                           Real theta,
                                           const SurfaceMobilityParameters& params)
{
    if (bulkMobility <= 0.0)
        return 0.0;
    if (theta < 0.0 || params.beta <= 0.0 || params.referenceField < 0.0 ||
        params.minFactor < 0.0 || params.maxFactor <= 0.0 ||
        params.minFactor > params.maxFactor)
        throw std::invalid_argument(
            "DopingDependentMobility: surface mobility parameters must be nonnegative "
            "with beta > 0 and min_factor <= max_factor.");
    if (theta == 0.0 || !std::isfinite(surfaceNormalField))
        return bulkMobility;

    const Real field = std::max<Real>(0.0, std::abs(surfaceNormalField) - params.referenceField);
    if (field <= 0.0)
        return bulkMobility;

    const Real thetaField = theta * field;
    Real factor = 1.0 / std::pow(1.0 + std::pow(thetaField, params.beta), 1.0 / params.beta);
    factor = std::clamp(factor, params.minFactor, params.maxFactor);
    return bulkMobility * factor;
}

MobilityModelConfig mobilityModelConfig(std::string modelName)
{
    MobilityModelConfig config;
    config.model = std::move(modelName);
    return config;
}

MobilityModelConfig mobilityModelConfigFromJson(const nlohmann::json& value)
{
    if (value.is_null())
        return {};
    if (value.is_string())
        return mobilityModelConfig(value.get<std::string>());
    if (!value.is_object())
        throw std::invalid_argument("mobility config must be a string or object.");

    MobilityModelConfig config;
    config.model = value.value("model", config.model);

    parseCaugheyThomas(value, config.electronCT, "electron");
    parseCaugheyThomas(value, config.holeCT, "hole");
    parseField(value, config.electronField, "electron");
    parseField(value, config.holeField, "hole");

    if (value.contains("surface")) {
        const auto& surface = value.at("surface");
        if (!surface.is_object())
            throw std::invalid_argument("mobility.surface must be an object.");
        config.surface.thetaElectron = surface.value(
            "theta_electron_m_per_V", config.surface.thetaElectron);
        config.surface.thetaHole = surface.value(
            "theta_hole_m_per_V", config.surface.thetaHole);
        config.surface.beta = surface.value("beta", config.surface.beta);
        config.surface.referenceField = surface.value(
            "reference_field_V_per_m", config.surface.referenceField);
        config.surface.minFactor = surface.value("min_factor", config.surface.minFactor);
        config.surface.maxFactor = surface.value("max_factor", config.surface.maxFactor);
        config.surface.surfaceRegion = surface.value(
            "surface_region", config.surface.surfaceRegion);
        if (surface.contains("surface_interface") && surface.contains("interface"))
            throw std::invalid_argument(
                "mobility.surface must not specify both surface_interface and interface.");
        if (surface.contains("surface_interface"))
            config.surface.surfaceInterface =
                surface.at("surface_interface").get<std::vector<std::string>>();
        else if (surface.contains("interface"))
            config.surface.surfaceInterface =
                surface.at("interface").get<std::vector<std::string>>();
    }

    return config;
}

bool isSurfaceMobilityModel(const MobilityModelConfig& config)
{
    return config.model == "caughey_thomas_surface" ||
           config.model == "caughey_thomas_field_surface";
}

bool surfaceMobilityAppliesToRegionPair(const MobilityModelConfig& config,
                                        const std::string& regionName,
                                        const std::vector<std::string>& adjacentRegionNames)
{
    if (!isSurfaceMobilityModel(config))
        return false;
    if (!config.surface.surfaceRegion.empty() &&
        config.surface.surfaceRegion != regionName)
        return false;
    if (config.surface.surfaceInterface.empty())
        return true;
    if (config.surface.surfaceInterface.size() != 2)
        throw std::invalid_argument(
            "surface mobility surface_interface must contain exactly two region names.");

    const std::string& a = config.surface.surfaceInterface[0];
    const std::string& b = config.surface.surfaceInterface[1];
    if (regionName != a && regionName != b)
        return false;
    const std::string& other = (regionName == a) ? b : a;
    return std::find(adjacentRegionNames.begin(), adjacentRegionNames.end(), other) !=
           adjacentRegionNames.end();
}

std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config)
{
    if (config.model == "constant")
        return std::make_unique<ConstantMobility>();
    if (config.model == "caughey_thomas" ||
        config.model == "caughey_thomas_field" ||
        isSurfaceMobilityModel(config))
        return std::make_unique<DopingDependentMobility>(config);

    throw std::invalid_argument(
        "makeMobilityModel: unknown mobility model '" + config.model + "'.");
}

} // namespace vela
