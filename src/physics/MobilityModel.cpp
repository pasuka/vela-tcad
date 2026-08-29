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
                        const char* prefix,
                        UnitScalingConfig scaling)
{
    const std::string muMinKey = std::string(prefix) + "_mu_min_m2_V_s";
    if (json.contains(muMinKey))
        params.muMin = scaling.mobilityToInternal(json.at(muMinKey).get<Real>());
    const std::string nRefKey = std::string(prefix) + "_nref_m3";
    if (json.contains(nRefKey))
        params.nRef = scaling.concentrationToInternal(json.at(nRefKey).get<Real>());
    params.alpha = json.value((std::string(prefix) + "_alpha").c_str(), params.alpha);
}

void parseMasetti(const nlohmann::json& json,
                  MasettiParameters& params,
                  const char* prefix,
                  UnitScalingConfig scaling)
{
    const std::string base = std::string(prefix) + "_";
    const std::pair<const char*, Real*> mobilityFields[] = {
        {"mu_const_m2_V_s", &params.muConst},
        {"mumin1_m2_V_s", &params.muMin1},
        {"mumin2_m2_V_s", &params.muMin2},
        {"mu1_m2_V_s", &params.mu1},
    };
    for (const auto& [name, target] : mobilityFields) {
        const std::string key = base + name;
        if (json.contains(key))
            *target = scaling.mobilityToInternal(json.at(key).get<Real>());
    }

    const std::pair<const char*, Real*> concentrationFields[] = {
        {"pc_m3", &params.pc},
        {"cr_m3", &params.cr},
        {"cs_m3", &params.cs},
    };
    for (const auto& [name, target] : concentrationFields) {
        const std::string key = base + name;
        if (json.contains(key))
            *target = scaling.concentrationToInternal(json.at(key).get<Real>());
    }

    params.alpha = json.value((base + "masetti_alpha").c_str(), params.alpha);
    params.beta = json.value((base + "masetti_beta").c_str(), params.beta);
}

bool isMasettiModel(const std::string& model)
{
    return model == "masetti" || model == "masetti_field" ||
           model == "masetti_surface" ||
           model == "masetti_field_surface" ||
           model == "masetti_lombardi" ||
           model == "masetti_field_lombardi";
}

bool isFieldMobilityModel(const std::string& model)
{
    return model == "constant_field" ||
           model == "caughey_thomas_field" ||
           model == "caughey_thomas_field_surface" ||
           model == "masetti_field" ||
           model == "masetti_field_surface" ||
           model == "masetti_field_lombardi";
}

bool isLombardiModel(const std::string& model)
{
    return model == "masetti_lombardi" ||
           model == "masetti_field_lombardi";
}

void parseField(const nlohmann::json& json,
                FieldMobilityParameters& params,
                const char* prefix,
                UnitScalingConfig scaling)
{
    const std::string velocityKey =
        std::string(prefix) + "_saturation_velocity_m_s";
    if (json.contains(velocityKey)) {
        params.saturationVelocity = scaling.velocityToInternal(
            json.at(velocityKey).get<Real>());
    }
    params.beta = json.value((std::string(prefix) + "_field_beta").c_str(), params.beta);
}

void validateHighFieldDrivingForce(const std::string& value)
{
    if (value != "electric_field" && value != "quasi_fermi_gradient")
        throw std::invalid_argument(
            "mobility.high_field_driving_force must be 'electric_field' or "
            "'quasi_fermi_gradient'.");
}

void validateHighFieldGradientDiscretization(const std::string& value)
{
    if (value != "edge_projection" && value != "transport_cell_vector") {
        throw std::invalid_argument(
            "mobility.high_field_gradient_discretization must be "
            "'edge_projection' or 'transport_cell_vector'.");
    }
}

void validateCarrierCurrentDiscretization(const std::string& value)
{
    if (value != "scharfetter_gummel_edge" &&
        value != "element_qf_gradient") {
        throw std::invalid_argument(
            "mobility.carrier_current_discretization must be "
            "'scharfetter_gummel_edge' or 'element_qf_gradient'.");
    }
}

void validateDopingConcentrationBasis(const std::string& value)
{
    if (value != "net_doping" && value != "total_impurity" &&
        value != "cell_reconstructed_total_impurity") {
        throw std::invalid_argument(
            "mobility.doping_concentration_basis must be 'net_doping', "
            "'total_impurity', or 'cell_reconstructed_total_impurity'.");
    }
}
void convertMobilityDefaultsToInternal(MobilityModelConfig& config,
                                       UnitScalingConfig scaling)
{
    const PhysicalUnitSystem& units = scaling.unitSystem();
    config.internalFieldToVPerM = units.electricFieldVPerMPerInternal();
    config.internalConcentrationToM3 = units.concentrationM3PerInternal();
    config.internalMobilityToM2PerVS = units.mobilityM2PerVSPerInternal();
    config.internalLengthToM = units.lengthMPerInternal();
    config.surface.coordinateFieldFactor =
        scaling.unitSystem().fieldFromCoordinateDeltaFactor();
    if (!scaling.isUnitScaling())
        return;

    auto convertCT = [&](CaugheyThomasParameters& params) {
        params.muMin = units.m2PerVSToInternalMobility(params.muMin);
        params.nRef = units.m3ToInternalConcentration(params.nRef);
    };
    auto convertMasetti = [&](MasettiParameters& params) {
        params.muConst = units.m2PerVSToInternalMobility(params.muConst);
        params.muMin1 = units.m2PerVSToInternalMobility(params.muMin1);
        params.muMin2 = units.m2PerVSToInternalMobility(params.muMin2);
        params.mu1 = units.m2PerVSToInternalMobility(params.mu1);
        params.pc = units.m3ToInternalConcentration(params.pc);
        params.cr = units.m3ToInternalConcentration(params.cr);
        params.cs = units.m3ToInternalConcentration(params.cs);
    };

    convertCT(config.electronCT);
    convertCT(config.holeCT);
    convertMasetti(config.electronMasetti);
    convertMasetti(config.holeMasetti);
    config.electronField.saturationVelocity =
        units.mPerSToInternalVelocity(
            config.electronField.saturationVelocity);
    config.holeField.saturationVelocity =
        units.mPerSToInternalVelocity(
            config.holeField.saturationVelocity);
    config.surface.thetaElectron =
        units.mPerVToInternalSurfaceFieldCoefficient(config.surface.thetaElectron);
    config.surface.thetaHole =
        units.mPerVToInternalSurfaceFieldCoefficient(config.surface.thetaHole);
    config.surface.referenceField =
        units.vPerMToInternalElectricField(config.surface.referenceField);
}

} // namespace

Real ConstantMobility::electronMobility(const Material& material,
                                        Real,
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
                                               Real n,
                                               Real p,
                                               Real electricField,
                                               Real surfaceNormalField,
                                               Real surfaceDistance) const
{
    Real mobility = config_.model == "constant_field"
        ? material.mun
        : (isMasettiModel(config_.model)
            ? masetti(netDoping, config_.electronMasetti)
            : caugheyThomas(material.mun, netDoping, config_.electronCT));
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.electronField);
    if (isSurfaceMobilityModel(config_))
        mobility = isLombardiModel(config_.model)
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance, CarrierType::Electron,
                            config_.electronLombardi)
            : surfaceLimit(mobility, surfaceNormalField,
                           config_.surface.thetaElectron, config_.surface);
    return mobility;
}

Real DopingDependentMobility::holeMobility(const Material& material,
                                           Real netDoping,
                                           Real n,
                                           Real p,
                                           Real electricField,
                                           Real surfaceNormalField,
                                           Real surfaceDistance) const
{
    Real mobility = config_.model == "constant_field"
        ? material.mup
        : (isMasettiModel(config_.model)
            ? masetti(netDoping, config_.holeMasetti)
            : caugheyThomas(material.mup, netDoping, config_.holeCT));
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.holeField);
    if (isSurfaceMobilityModel(config_))
        mobility = isLombardiModel(config_.model)
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance, CarrierType::Hole,
                            config_.holeLombardi)
            : surfaceLimit(mobility, surfaceNormalField,
                           config_.surface.thetaHole, config_.surface);
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

Real DopingDependentMobility::masetti(Real netDoping,
                                      const MasettiParameters& params)
{
    if (params.muConst <= 0.0)
        return 0.0;
    if (params.cr <= 0.0 || params.cs <= 0.0 || params.alpha <= 0.0 ||
        params.beta <= 0.0)
        throw std::invalid_argument(
            "DopingDependentMobility: Masetti cr, cs, alpha, and beta must be positive.");

    const Real doping = std::abs(netDoping);
    if (doping <= 0.0)
        return params.muConst;

    const Real exponential =
        params.muMin1 * std::exp(-std::max<Real>(0.0, params.pc) / doping);
    const Real rolloff = (params.muConst - params.muMin2) /
        (1.0 + std::pow(doping / params.cr, params.alpha));
    const Real highDopingCorrection = params.mu1 /
        (1.0 + std::pow(params.cs / doping, params.beta));
    const Real mobility = exponential + rolloff - highDopingCorrection;
    return std::max<Real>(0.0, mobility);
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

Real DopingDependentMobility::lombardiLimit(
    Real bulkMobility,
    Real netDoping,
    Real n,
    Real p,
    Real surfaceNormalField,
    Real surfaceDistance,
    CarrierType carrier,
    const LombardiParameters& params) const
{
    if (bulkMobility <= 0.0 || !std::isfinite(surfaceNormalField) ||
        !std::isfinite(surfaceDistance))
        return bulkMobility;
    if (params.B <= 0.0 || params.C < 0.0 || params.N0 <= 0.0 ||
        params.N1 <= 0.0 || params.N2 < 0.0 || params.delta <= 0.0 ||
        params.eta <= 0.0 || params.criticalLength <= 0.0 ||
        params.nu <= 0.0 || params.acousticFactor < 0.0 ||
        params.roughnessFactor < 0.0) {
        throw std::invalid_argument(
            "DopingDependentMobility: invalid Enhanced Lombardi parameters.");
    }

    const Real field = std::abs(surfaceNormalField) * config_.internalFieldToVPerM;
    if (field <= 0.0)
        return bulkMobility;
    const Real distance = std::max<Real>(0.0, surfaceDistance) *
        config_.internalLengthToM;
    const Real damping = std::exp(-distance / params.criticalLength);
    if (damping <= std::numeric_limits<Real>::min())
        return bulkMobility;

    const Real doping = std::abs(netDoping) *
        config_.internalConcentrationToM3;
    const Real electronDensity = std::max<Real>(0.0, n) *
        config_.internalConcentrationToM3;
    const Real holeDensity = std::max<Real>(0.0, p) *
        config_.internalConcentrationToM3;
    const Real acousticMobility = params.B / field +
        params.C * std::pow((doping + params.N2) / params.N0, params.lambda) /
        std::cbrt(field); // T=300 K for the reference BVmethods deck.

    const Real sameCarrier = carrier == CarrierType::Electron
        ? electronDensity : holeDensity;
    const Real otherCarrier = carrier == CarrierType::Electron
        ? holeDensity : electronDensity;
    const Real exponent = params.A +
        params.alpha * (sameCarrier + params.aOther * otherCarrier) /
        std::pow(doping + params.N1, params.nu);
    const Real referenceField = 100.0; // 1 V/cm in V/m.
    const Real inverseRoughness =
        std::pow(field / referenceField, exponent) / params.delta +
        std::pow(field, 3.0) / params.eta;
    if (!(acousticMobility > 0.0) || !(inverseRoughness > 0.0))
        return bulkMobility;

    const Real bulkSI = bulkMobility * config_.internalMobilityToM2PerVS;
    const Real inverseTotal = 1.0 / bulkSI +
        damping * params.acousticFactor / acousticMobility +
        damping * params.roughnessFactor * inverseRoughness;
    return (1.0 / inverseTotal) / config_.internalMobilityToM2PerVS;
}

MobilityModelConfig mobilityModelConfig(std::string modelName)
{
    MobilityModelConfig config;
    config.model = std::move(modelName);
    validateHighFieldDrivingForce(config.highFieldDrivingForce);
    validateHighFieldGradientDiscretization(
        config.highFieldGradientDiscretization);
    validateCarrierCurrentDiscretization(
        config.carrierCurrentDiscretization);
    validateDopingConcentrationBasis(config.dopingConcentrationBasis);
    return config;
}

MobilityModelConfig mobilityModelConfigFromJson(
    const nlohmann::json& value,
    UnitScalingConfig scaling)
{
    if (value.is_null())
        return {};
    if (value.is_string()) {
        MobilityModelConfig config;
        convertMobilityDefaultsToInternal(config, scaling);
        config.model = value.get<std::string>();
        validateHighFieldDrivingForce(config.highFieldDrivingForce);
        validateHighFieldGradientDiscretization(
            config.highFieldGradientDiscretization);
        validateCarrierCurrentDiscretization(
            config.carrierCurrentDiscretization);
        validateDopingConcentrationBasis(config.dopingConcentrationBasis);
        return config;
    }
    if (!value.is_object())
        throw std::invalid_argument("mobility config must be a string or object.");

    MobilityModelConfig config;
    convertMobilityDefaultsToInternal(config, scaling);
    config.model = value.value("model", config.model);
    config.highFieldDrivingForce = value.value(
        "high_field_driving_force", config.highFieldDrivingForce);
    config.highFieldGradientDiscretization = value.value(
        "high_field_gradient_discretization",
        config.highFieldGradientDiscretization);
    config.carrierCurrentDiscretization = value.value(
        "carrier_current_discretization",
        config.carrierCurrentDiscretization);
    config.dopingConcentrationBasis = value.value(
        "doping_concentration_basis", config.dopingConcentrationBasis);
    config.jacobianFieldDerivatives = value.value(
        "jacobian_field_derivatives", config.jacobianFieldDerivatives);
    validateHighFieldDrivingForce(config.highFieldDrivingForce);
    validateHighFieldGradientDiscretization(
        config.highFieldGradientDiscretization);
    validateCarrierCurrentDiscretization(
        config.carrierCurrentDiscretization);
    validateDopingConcentrationBasis(config.dopingConcentrationBasis);

    parseCaugheyThomas(value, config.electronCT, "electron", scaling);
    parseCaugheyThomas(value, config.holeCT, "hole", scaling);
    parseMasetti(value, config.electronMasetti, "electron", scaling);
    parseMasetti(value, config.holeMasetti, "hole", scaling);
    parseField(value, config.electronField, "electron", scaling);
    parseField(value, config.holeField, "hole", scaling);

    if (value.contains("surface")) {
        const auto& surface = value.at("surface");
        if (!surface.is_object())
            throw std::invalid_argument("mobility.surface must be an object.");
        if (surface.contains("theta_electron_m_per_V")) {
            config.surface.thetaElectron = scaling.surfaceFieldCoefficientToInternal(
                surface.at("theta_electron_m_per_V").get<Real>());
        }
        if (surface.contains("theta_hole_m_per_V")) {
            config.surface.thetaHole = scaling.surfaceFieldCoefficientToInternal(
                surface.at("theta_hole_m_per_V").get<Real>());
        }
        config.surface.beta = surface.value("beta", config.surface.beta);
        if (surface.contains("reference_field_V_per_m")) {
            config.surface.referenceField = scaling.electricFieldToInternal(
                surface.at("reference_field_V_per_m").get<Real>());
        }
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
           config.model == "caughey_thomas_field_surface" ||
           config.model == "masetti_surface" ||
           config.model == "masetti_field_surface" ||
           isLombardiModel(config.model);
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
    if (config.model == "constant_field" ||
        config.model == "caughey_thomas" ||
        config.model == "caughey_thomas_field" ||
        isMasettiModel(config.model) ||
        isSurfaceMobilityModel(config))
        return std::make_unique<DopingDependentMobility>(config);

    throw std::invalid_argument(
        "makeMobilityModel: unknown mobility model '" + config.model + "'.");
}

} // namespace vela
