#include "vela/physics/MobilityModel.h"
#include "vela/equation/IalTransport.h"
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

void convertPhuMobDefaultsToInternal(PhuMobParameters& params,
                                     const PhysicalUnitSystem& units)
{
    auto convertCarrier = [&](PhuMobCarrierParameters& carrier) {
        carrier.muMax = units.m2PerVSToInternalMobility(carrier.muMax);
        carrier.muMin = units.m2PerVSToInternalMobility(carrier.muMin);
        carrier.nRef = units.m3ToInternalConcentration(carrier.nRef);
    };
    convertCarrier(params.electronArsenic);
    convertCarrier(params.electronPhosphorus);
    convertCarrier(params.holeBoron);
    params.donorClusterReference = units.m3ToInternalConcentration(
        params.donorClusterReference);
    params.acceptorClusterReference = units.m3ToInternalConcentration(
        params.acceptorClusterReference);
    params.internalConcentrationToCm3 =
        units.concentrationM3PerInternal() * 1.0e-6;
    params.internalMobilityToCm2PerVS =
        units.mobilityM2PerVSPerInternal() * 1.0e4;
}

void parsePhuMobCarrier(const nlohmann::json& json,
                        PhuMobCarrierParameters& params,
                        UnitScalingConfig scaling)
{
    if (!json.is_object())
        throw std::invalid_argument("mobility.phumob carrier parameters must be an object.");
    if (json.contains("mu_max_m2_V_s")) {
        params.muMax = scaling.mobilityToInternal(
            json.at("mu_max_m2_V_s").get<Real>());
    }
    if (json.contains("mu_min_m2_V_s")) {
        params.muMin = scaling.mobilityToInternal(
            json.at("mu_min_m2_V_s").get<Real>());
    }
    if (json.contains("n_ref_m3")) {
        params.nRef = scaling.concentrationToInternal(
            json.at("n_ref_m3").get<Real>());
    }
    params.theta = json.value("theta", params.theta);
    params.alpha = json.value("alpha", params.alpha);
}

void parsePhuMob(const nlohmann::json& json,
                 PhuMobParameters& params,
                 UnitScalingConfig scaling)
{
    if (!json.is_object())
        throw std::invalid_argument("mobility.phumob must be an object.");

    const std::string donorSpecies = json.value("donor_species", "arsenic");
    if (donorSpecies == "arsenic")
        params.donorSpecies = PhuMobDonorSpecies::Arsenic;
    else if (donorSpecies == "phosphorus")
        params.donorSpecies = PhuMobDonorSpecies::Phosphorus;
    else
        throw std::invalid_argument(
            "mobility.phumob.donor_species must be 'arsenic' or 'phosphorus'.");

    if (json.contains("electron_arsenic")) {
        parsePhuMobCarrier(
            json.at("electron_arsenic"), params.electronArsenic, scaling);
    }
    if (json.contains("electron_phosphorus")) {
        parsePhuMobCarrier(
            json.at("electron_phosphorus"), params.electronPhosphorus, scaling);
    }
    if (json.contains("hole_boron")) {
        parsePhuMobCarrier(json.at("hole_boron"), params.holeBoron, scaling);
    }

    const std::pair<const char*, Real*> concentrationFields[] = {
        {"donor_cluster_reference_m3", &params.donorClusterReference},
        {"acceptor_cluster_reference_m3", &params.acceptorClusterReference},
    };
    for (const auto& [name, target] : concentrationFields) {
        if (json.contains(name))
            *target = scaling.concentrationToInternal(json.at(name).get<Real>());
    }

    const std::pair<const char*, Real*> dimensionlessFields[] = {
        {"donor_cluster_coefficient", &params.donorClusterCoefficient},
        {"acceptor_cluster_coefficient", &params.acceptorClusterCoefficient},
        {"electron_mass_ratio", &params.electronMassRatio},
        {"hole_mass_ratio", &params.holeMassRatio},
        {"conwell_weisskopf_factor", &params.conwellWeisskopfFactor},
        {"brooks_herring_factor", &params.brooksHerringFactor},
        {"electron_hole_scattering_factor", &params.electronHoleScatteringFactor},
        {"hole_electron_scattering_factor", &params.holeElectronScatteringFactor},
        {"g_a", &params.gA},
        {"g_b", &params.gB},
        {"g_c", &params.gC},
        {"g_alpha", &params.gAlpha},
        {"g_alpha_prime", &params.gAlphaPrime},
        {"g_beta", &params.gBeta},
        {"g_gamma", &params.gGamma},
    };
    for (const auto& [name, target] : dimensionlessFields)
        *target = json.value(name, *target);
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
           model == "constant_field_lombardi" ||
           model == "caughey_thomas_field" ||
           model == "caughey_thomas_field_surface" ||
           model == "masetti_field" ||
           model == "masetti_field_surface" ||
           model == "masetti_field_lombardi" ||
           model == "phumob_field" ||
           model == "phumob_field_lombardi";
}

bool isLombardiModel(const std::string& model)
{
    return model == "constant_lombardi" ||
           model == "constant_field_lombardi" ||
           model == "masetti_lombardi" ||
           model == "masetti_field_lombardi" ||
           model == "phumob_lombardi" ||
           model == "phumob_field_lombardi";
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
    if (value != "edge_projection" && value != "transport_cell_vector" && value != "element_vertex_partial_layer") {
        throw std::invalid_argument(
            "mobility.high_field_gradient_discretization must be "
            "'edge_projection' or 'transport_cell_vector'.");
    }
}

void validateContactElectricFieldFallback(const MobilityModelConfig& config)
{
    if (config.contactElectricFieldFallbackScope != "contact_node_cell") {
        throw std::invalid_argument(
            "mobility.contact_electric_field_fallback_scope must be "
            "'contact_node_cell'.");
    }
    if (config.contactElectricFieldFallbackMode != "cell_gradient_magnitude") {
        throw std::invalid_argument(
            "mobility.contact_electric_field_fallback_mode must be "
            "'cell_gradient_magnitude'.");
    }
    if (config.contactElectricFieldFallback &&
        config.highFieldDrivingForce != "quasi_fermi_gradient") {
        throw std::invalid_argument(
            "mobility.contact_electric_field_fallback requires "
            "high_field_driving_force='quasi_fermi_gradient'.");
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
    convertPhuMobDefaultsToInternal(config.phuMob, units);
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

Real MobilityModel::electronMobilityWithIonizedImpurities(
    const Material& material,
    Real donors,
    Real acceptors,
    Real n,
    Real p,
    Real electricField,
    Real surfaceNormalField,
    Real surfaceDistance) const
{
    return electronMobility(
        material, donors - acceptors, n, p, electricField,
        surfaceNormalField, surfaceDistance);
}

Real MobilityModel::holeMobilityWithIonizedImpurities(
    const Material& material,
    Real donors,
    Real acceptors,
    Real n,
    Real p,
    Real electricField,
    Real surfaceNormalField,
    Real surfaceDistance) const
{
    return holeMobility(
        material, donors - acceptors, n, p, electricField,
        surfaceNormalField, surfaceDistance);
}

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
    if (isPhuMobModel(config_)) {
        return electronMobilityWithIonizedImpurities(
            material, std::max<Real>(netDoping, 0.0),
            std::max<Real>(-netDoping, 0.0), n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    Real mobility = config_.model.rfind("constant_", 0) == 0
        ? material.mun
        : isMasettiModel(config_.model)
            ? masetti(netDoping, config_.electronMasetti)
            : caugheyThomas(material.mun, netDoping, config_.electronCT);
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.electronField);
    if (isSurfaceMobilityModel(config_))
        mobility = isLombardiModel(config_.model)
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance,
                            material.temperature_K.value_or(300.0),
                            CarrierType::Electron,
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
    if (isPhuMobModel(config_)) {
        return holeMobilityWithIonizedImpurities(
            material, std::max<Real>(netDoping, 0.0),
            std::max<Real>(-netDoping, 0.0), n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    Real mobility = config_.model.rfind("constant_", 0) == 0
        ? material.mup
        : isMasettiModel(config_.model)
            ? masetti(netDoping, config_.holeMasetti)
            : caugheyThomas(material.mup, netDoping, config_.holeCT);
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.holeField);
    if (isSurfaceMobilityModel(config_))
        mobility = isLombardiModel(config_.model)
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance,
                            material.temperature_K.value_or(300.0),
                            CarrierType::Hole,
                            config_.holeLombardi)
            : surfaceLimit(mobility, surfaceNormalField,
                           config_.surface.thetaHole, config_.surface);
    return mobility;
}

Real DopingDependentMobility::electronMobilityWithIonizedImpurities(
    const Material& material,
    Real donors,
    Real acceptors,
    Real n,
    Real p,
    Real electricField,
    Real surfaceNormalField,
    Real surfaceDistance) const
{
    if (!isPhuMobModel(config_)) {
        return MobilityModel::electronMobilityWithIonizedImpurities(
            material, donors, acceptors, n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    if (material.mun <= 0.0)
        return 0.0;
    const PhuMobScalarState state{
        donors,
        acceptors,
        n,
        p,
        material.temperature_K.value_or(300.0),
    };
    Real mobility = evaluatePhuMobScalar(
        CarrierType::Electron, state, config_.phuMob).mobility;
    // EnormalDependence is part of the low-field mobility.  The high-field
    // saturation law therefore consumes the PhuMob + Enormal result.
    if (isLombardiModel(config_.model)) {
        mobility = lombardiLimit(
            mobility, donors + acceptors, n, p, surfaceNormalField,
            surfaceDistance, state.temperature_K, CarrierType::Electron,
            config_.electronLombardi);
    }
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.electronField);
    return mobility;
}

Real DopingDependentMobility::holeMobilityWithIonizedImpurities(
    const Material& material,
    Real donors,
    Real acceptors,
    Real n,
    Real p,
    Real electricField,
    Real surfaceNormalField,
    Real surfaceDistance) const
{
    if (!isPhuMobModel(config_)) {
        return MobilityModel::holeMobilityWithIonizedImpurities(
            material, donors, acceptors, n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    if (material.mup <= 0.0)
        return 0.0;
    const PhuMobScalarState state{
        donors,
        acceptors,
        n,
        p,
        material.temperature_K.value_or(300.0),
    };
    Real mobility = evaluatePhuMobScalar(
        CarrierType::Hole, state, config_.phuMob).mobility;
    if (isLombardiModel(config_.model)) {
        mobility = lombardiLimit(
            mobility, donors + acceptors, n, p, surfaceNormalField,
            surfaceDistance, state.temperature_K, CarrierType::Hole,
            config_.holeLombardi);
    }
    if (isFieldMobilityModel(config_.model))
        mobility = fieldLimit(mobility, electricField, config_.holeField);
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
    Real latticeTemperature_K,
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
        params.roughnessFactor < 0.0 || latticeTemperature_K <= 0.0 ||
        !std::isfinite(latticeTemperature_K)) {
        throw std::invalid_argument(
            "DopingDependentMobility: invalid Enhanced Lombardi parameters.");
    }

    const Real field = std::abs(surfaceNormalField) * config_.internalFieldToVPerM;
    if (field <= 0.0)
        return bulkMobility;
    const Real distance = std::max<Real>(0.0, surfaceDistance) *
        config_.internalLengthToM;
    // Sentaurus documents l_crit > 100 cm as the no-distance-damping
    // sentinel.  Parameters are stored in metres here.
    const Real damping = params.criticalLength > 1.0
        ? 1.0 : std::exp(-distance / params.criticalLength);
    if (damping <= std::numeric_limits<Real>::min())
        return bulkMobility;

    const Real doping = std::abs(netDoping) *
        config_.internalConcentrationToM3;
    const Real electronDensity = std::max<Real>(0.0, n) *
        config_.internalConcentrationToM3;
    const Real holeDensity = std::max<Real>(0.0, p) *
        config_.internalConcentrationToM3;
    const Real temperatureRatio = latticeTemperature_K / 300.0;
    const Real acousticMobility = params.B / field +
        params.C * std::pow(temperatureRatio, -params.k) *
        std::pow((doping + params.N2) / params.N0, params.lambda) /
        std::cbrt(field);

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
    validateContactElectricFieldFallback(config);
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
        validateContactElectricFieldFallback(config);
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
    config.edgeAveraging = value.value("edge_averaging", "legacy");
    if (config.edgeAveraging != "legacy" && config.edgeAveraging != "element_box" &&
        config.edgeAveraging != "element_box_phumob")
        throw std::invalid_argument("mobility.edge_averaging must be legacy, element_box or element_box_phumob.");
    if (config.edgeAveraging == "element_box" && config.model != "masetti" && config.model != "constant")
        throw std::invalid_argument("element_box mobility supports constant and Masetti without field or carrier dependence.");
    if (config.edgeAveraging == "element_box_phumob" && config.model != "phumob" && config.model != "phumob_lombardi" && config.model != "phumob_field_lombardi")
        throw std::invalid_argument("element_box_phumob is an explicit plain PhuMob candidate; field and surface models are unsupported.");
    config.highFieldDrivingForce = value.value(
        "high_field_driving_force", config.highFieldDrivingForce);
    config.highFieldGradientDiscretization = value.value(
        "high_field_gradient_discretization",
        config.highFieldGradientDiscretization);
    config.contactElectricFieldFallback = value.value(
        "contact_electric_field_fallback",
        config.contactElectricFieldFallback);
    config.contactElectricFieldFallbackScope = value.value(
        "contact_electric_field_fallback_scope",
        config.contactElectricFieldFallbackScope);
    config.contactElectricFieldFallbackMode = value.value(
        "contact_electric_field_fallback_mode",
        config.contactElectricFieldFallbackMode);
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
    validateContactElectricFieldFallback(config);
    validateCarrierCurrentDiscretization(
        config.carrierCurrentDiscretization);
    validateDopingConcentrationBasis(config.dopingConcentrationBasis);

    parseCaugheyThomas(value, config.electronCT, "electron", scaling);
    parseCaugheyThomas(value, config.holeCT, "hole", scaling);
    parseMasetti(value, config.electronMasetti, "electron", scaling);
    parseMasetti(value, config.holeMasetti, "hole", scaling);
    if (value.contains("phumob"))
        parsePhuMob(value.at("phumob"), config.phuMob, scaling);
    parseField(value, config.electronField, "electron", scaling);
    parseField(value, config.holeField, "hole", scaling);

    if (value.contains("surface")) {
        const auto& surface = value.at("surface");
        if (!surface.is_object())
            throw std::invalid_argument("mobility.surface must be an object.");
        config.surface.discretization=surface.value("discretization",config.surface.discretization);
        if(config.surface.discretization!="legacy_cell_centroid" && config.surface.discretization!="element_distance_gradient")
            throw std::invalid_argument("Unknown surface mobility discretization");
        const Real ac=surface.value("acoustic_factor",1.),sr=surface.value("roughness_factor",1.);
        if(!std::isfinite(ac)||!std::isfinite(sr)||ac<0||sr<0)
            throw std::invalid_argument("Surface scattering factors must be finite and nonnegative");
        config.electronLombardi.acousticFactor=config.holeLombardi.acousticFactor=ac;
        config.electronLombardi.roughnessFactor=config.holeLombardi.roughnessFactor=sr;
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

    if (config.model == "ialmob") {
        if (config.carrierCurrentDiscretization != "scharfetter_gummel_edge" ||
            config.highFieldDrivingForce != "quasi_fermi_gradient")
            throw std::invalid_argument("IALMob requires SG edge transport and QF high-field drive");
        config.ialmob = ialTransportOptionsFromJson(value.at("ialmob"));
    } else if (value.contains("ialmob")) {
        throw std::invalid_argument("ialmob block requires model=ialmob");
    }
    const bool hfs=config.highFieldGradientDiscretization=="element_vertex_partial_layer";
    if(hfs != (config.model=="phumob_field_lombardi" && config.edgeAveraging=="element_box_phumob") ||
       (hfs && (config.highFieldDrivingForce!="quasi_fermi_gradient" || !config.jacobianFieldDerivatives)))
        throw std::invalid_argument("Element HFS requires its explicit PhuMob/Enormal box profile and complete QF field derivatives");
    const bool elementSurface=config.surface.discretization=="element_distance_gradient";
    if(elementSurface != (config.edgeAveraging=="element_box_phumob" && (config.model=="phumob_lombardi" || config.model=="phumob_field_lombardi")))
        throw std::invalid_argument("Element-box PhuMob Lombardi requires explicit element_distance_gradient surface discretization");
    if(elementSurface && (!config.surface.surfaceRegion.empty() || !config.surface.surfaceInterface.empty()))
        throw std::invalid_argument("Element-distance Lombardi currently qualifies the default semiconductor-insulator interface selector only");
    return config;
}

bool isSurfaceMobilityModel(const MobilityModelConfig& config)
{
    return config.model == "constant_lombardi" ||
           config.model == "constant_field_lombardi" ||
           config.model == "caughey_thomas_surface" ||
           config.model == "caughey_thomas_field_surface" ||
           config.model == "masetti_surface" ||
           config.model == "masetti_field_surface" ||
           isLombardiModel(config.model);
}

bool isPhuMobModel(const MobilityModelConfig& config)
{
    return config.model == "phumob" || config.model == "phumob_field" ||
           config.model == "phumob_lombardi" ||
           config.model == "phumob_field_lombardi";
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

namespace {
class IalCellOnlyMobility final : public MobilityModel {
public:
    Real electronMobility(const Material&,Real,Real,Real,Real,Real,Real) const override {
        throw std::logic_error("IALMob requires the live element transport interface");
    }
    Real holeMobility(const Material&,Real,Real,Real,Real,Real,Real) const override {
        throw std::logic_error("IALMob requires the live element transport interface");
    }
};
}

std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config)
{
    if (config.model == "ialmob") {
        if (!config.ialmob) throw std::invalid_argument("IALMob requires an explicit configuration block");
        return std::make_unique<IalCellOnlyMobility>();
    }
    if (config.model == "constant")
        return std::make_unique<ConstantMobility>();
    if (config.model == "constant_field" ||
        config.model == "constant_lombardi" ||
        config.model == "constant_field_lombardi" ||
        config.model == "caughey_thomas" ||
        config.model == "caughey_thomas_field" ||
        isMasettiModel(config.model) ||
        isPhuMobModel(config) ||
        isSurfaceMobilityModel(config))
        return std::make_unique<DopingDependentMobility>(config);

    throw std::invalid_argument(
        "makeMobilityModel: unknown mobility model '" + config.model + "'.");
}

namespace {

void validatePhuMobCarrier(const PhuMobCarrierParameters& carrier,
                           const char* name)
{
    if (!std::isfinite(carrier.muMax) || !std::isfinite(carrier.muMin) ||
        !std::isfinite(carrier.theta) || !std::isfinite(carrier.nRef) ||
        !std::isfinite(carrier.alpha) || carrier.muMax <= carrier.muMin ||
        carrier.muMin < 0.0 || carrier.theta < 0.0 || carrier.nRef <= 0.0 ||
        carrier.alpha <= 0.0) {
        throw std::invalid_argument(
            std::string("PhuMob: invalid ") + name + " carrier parameters.");
    }
}

void validatePhuMob(const PhuMobParameters& params,
                    const PhuMobScalarState& state)
{
    validatePhuMobCarrier(params.electronArsenic, "electron-arsenic");
    validatePhuMobCarrier(params.electronPhosphorus, "electron-phosphorus");
    validatePhuMobCarrier(params.holeBoron, "hole-boron");
    const Real values[] = {
        params.donorClusterReference, params.acceptorClusterReference,
        params.donorClusterCoefficient, params.acceptorClusterCoefficient,
        params.electronMassRatio, params.holeMassRatio,
        params.conwellWeisskopfFactor, params.brooksHerringFactor,
        params.electronHoleScatteringFactor,
        params.holeElectronScatteringFactor, params.gA, params.gB, params.gC,
        params.gAlpha, params.gAlphaPrime, params.gBeta, params.gGamma,
        params.internalConcentrationToCm3, params.internalMobilityToCm2PerVS,
    };
    for (const Real value : values) {
        if (!std::isfinite(value))
            throw std::invalid_argument("PhuMob: parameters must be finite.");
    }
    if (params.donorClusterReference <= 0.0 ||
        params.acceptorClusterReference <= 0.0 ||
        params.donorClusterCoefficient < 0.0 ||
        params.acceptorClusterCoefficient < 0.0 ||
        params.electronMassRatio <= 0.0 || params.holeMassRatio <= 0.0 ||
        params.conwellWeisskopfFactor <= 0.0 ||
        params.brooksHerringFactor <= 0.0 ||
        params.electronHoleScatteringFactor < 0.0 ||
        params.holeElectronScatteringFactor < 0.0 || params.gA <= 0.0 ||
        params.gB <= 0.0 || params.gC <= 0.0 || params.gAlpha <= 0.0 ||
        params.gAlphaPrime <= 0.0 || params.gBeta <= 0.0 ||
        params.gGamma <= 0.0 || params.internalConcentrationToCm3 <= 0.0 ||
        params.internalMobilityToCm2PerVS <= 0.0) {
        throw std::invalid_argument("PhuMob: parameters must be physically positive.");
    }

    const Real stateValues[] = {
        state.donors, state.acceptors, state.electrons, state.holes,
        state.temperature_K,
    };
    for (const Real value : stateValues) {
        if (!std::isfinite(value))
            throw std::invalid_argument("PhuMob: scalar state must be finite.");
    }
    if (state.donors < 0.0 || state.acceptors < 0.0 ||
        state.electrons < 0.0 || state.holes < 0.0 ||
        state.temperature_K <= 0.0) {
        throw std::invalid_argument(
            "PhuMob: concentrations must be nonnegative and temperature positive.");
    }
}

Real clusteredConcentration(Real concentration, Real reference, Real coefficient)
{
    if (concentration <= 0.0)
        return 0.0;
    const Real ratio = concentration / reference;
    const Real ratioSquared = ratio * ratio;
    const Real correction = ratioSquared /
        (1.0 + coefficient * ratioSquared);
    return concentration * (1.0 + correction);
}

Real screeningF(Real screening, Real massRatio, Real otherMassRatio)
{
    if (!std::isfinite(screening))
        return 0.7643;
    const Real power = std::pow(screening, 0.6478);
    const Real massQuotient = massRatio / otherMassRatio;
    return (0.7643 * power + 2.2999 + 6.5502 * massQuotient) /
        (power + 2.3670 - 0.8552 * massQuotient);
}

Real rawScreeningG(Real screening, Real massRatio, Real temperature_K,
                   const PhuMobParameters& params)
{
    if (!std::isfinite(screening))
        return 1.0;
    const Real temperatureRatio = temperature_K / 300.0;
    const Real firstScale = std::pow(
        temperatureRatio / massRatio, params.gAlpha);
    const Real secondScale = std::pow(
        massRatio / temperatureRatio, params.gAlphaPrime);
    return 1.0 - params.gA * std::pow(
        params.gB + screening * firstScale, -params.gBeta) +
        params.gC * std::pow(screening * secondScale, -params.gGamma);
}

bool screeningGDerivativeIsPositive(Real logScreening, Real massRatio,
                                    Real temperature_K,
                                    const PhuMobParameters& params)
{
    const Real screening = std::exp(logScreening);
    const Real temperatureRatio = temperature_K / 300.0;
    const Real firstScale = std::pow(
        temperatureRatio / massRatio, params.gAlpha);
    const Real secondScale = std::pow(
        massRatio / temperatureRatio, params.gAlphaPrime);
    const Real logPositive = std::log(
        params.gA * params.gBeta * firstScale) -
        (params.gBeta + 1.0) * std::log(
            params.gB + screening * firstScale);
    const Real logNegative = std::log(params.gC * params.gGamma) -
        params.gGamma * std::log(secondScale) -
        (params.gGamma + 1.0) * logScreening;
    return logPositive > logNegative;
}

Real screeningG(Real screening, Real massRatio, Real temperature_K,
                const PhuMobParameters& params, Real* derivative)
{
    *derivative = 0.0;
    if (!std::isfinite(screening))
        return 1.0;
    const Real logFloor = std::log(std::numeric_limits<Real>::min());
    Real low = std::max(logFloor, -80.0);
    Real high = 0.0;
    while (!screeningGDerivativeIsPositive(
               high, massRatio, temperature_K, params) && high < 80.0) {
        high += 4.0;
    }
    for (int iteration = 0; iteration < 80; ++iteration) {
        const Real middle = 0.5 * (low + high);
        if (screeningGDerivativeIsPositive(
                middle, massRatio, temperature_K, params)) {
            high = middle;
        } else {
            low = middle;
        }
    }
    const Real minimumScreening = std::exp(0.5 * (low + high));
    const Real minimum = rawScreeningG(
        minimumScreening, massRatio, temperature_K, params);
    const Real minimumFloor = std::abs(minimum);
    if (screening < minimumScreening)
        return minimumFloor;
    const Real raw = rawScreeningG(screening, massRatio, temperature_K, params);
    if (raw > minimumFloor) {
        const Real firstScale = std::pow(temperature_K / (300.0 * massRatio), params.gAlpha);
        const Real secondScale = std::pow(300.0 * massRatio / temperature_K, params.gAlphaPrime);
        *derivative = params.gA * params.gBeta * firstScale *
            std::pow(params.gB + screening * firstScale, -params.gBeta - 1.0) -
            params.gC * params.gGamma * std::pow(screening * secondScale, -params.gGamma) / screening;
    }
    return std::max(raw, minimumFloor);
}

} // namespace

PhuMobScalarResult evaluatePhuMobScalar(
    CarrierType carrier,
    const PhuMobScalarState& state,
    const PhuMobParameters& params)
{
    validatePhuMob(params, state);
    const auto& carrierParams = carrier == CarrierType::Electron
        ? (params.donorSpecies == PhuMobDonorSpecies::Arsenic
               ? params.electronArsenic
               : params.electronPhosphorus)
        : params.holeBoron;
    const Real massRatio = carrier == CarrierType::Electron
        ? params.electronMassRatio : params.holeMassRatio;
    const Real otherMassRatio = carrier == CarrierType::Electron
        ? params.holeMassRatio : params.electronMassRatio;

    const Real concentrationScale = params.internalConcentrationToCm3;
    const Real mobilityScale = params.internalMobilityToCm2PerVS;
    const Real donors = state.donors * concentrationScale;
    const Real acceptors = state.acceptors * concentrationScale;
    const Real electrons = state.electrons * concentrationScale;
    const Real holes = state.holes * concentrationScale;
    const Real donorReference = params.donorClusterReference * concentrationScale;
    const Real acceptorReference =
        params.acceptorClusterReference * concentrationScale;
    const Real donorCluster = clusteredConcentration(
        donors, donorReference, params.donorClusterCoefficient);
    const Real acceptorCluster = clusteredConcentration(
        acceptors, acceptorReference, params.acceptorClusterCoefficient);

    const Real muMax = carrierParams.muMax * mobilityScale;
    const Real muMin = carrierParams.muMin * mobilityScale;
    const Real nRef = carrierParams.nRef * concentrationScale;
    const Real temperatureRatio = state.temperature_K / 300.0;
    const Real latticeMobility = muMax * std::pow(
        temperatureRatio, -carrierParams.theta);

    const Real otherCarrier = carrier == CarrierType::Electron
        ? holes : electrons;
    const Real scatteringConcentration =
        donorCluster + acceptorCluster + otherCarrier;

    PhuMobScalarResult result;
    result.latticeMobility = latticeMobility / mobilityScale;
    result.scatteringConcentration = scatteringConcentration / concentrationScale;
    if (scatteringConcentration <= 0.0) {
        result.mobility = result.latticeMobility;
        result.scatteringMobility = std::numeric_limits<Real>::infinity();
        result.effectiveScatteringConcentration = 0.0;
        result.screeningParameter = std::numeric_limits<Real>::infinity();
        result.screeningF = 0.7643;
        result.screeningG = 1.0;
        return result;
    }

    const Real freeCarriers = electrons + holes;
    const Real screeningDenominator =
        params.conwellWeisskopfFactor *
            std::pow(scatteringConcentration, 2.0 / 3.0) / 3.97e13 +
        params.brooksHerringFactor * freeCarriers /
            (1.36e20 * massRatio);
    const Real screening = screeningDenominator > 0.0
        ? temperatureRatio * temperatureRatio / screeningDenominator
        : std::numeric_limits<Real>::infinity();
    const Real functionF = screeningF(screening, massRatio, otherMassRatio);
    const Real functionG = screeningG(
        screening, massRatio, state.temperature_K, params, &result.screeningGDerivative);

    const Real effectiveScatteringConcentration = carrier == CarrierType::Electron
        ? donorCluster + functionG * acceptorCluster +
            params.electronHoleScatteringFactor * holes / functionF
        : acceptorCluster + functionG * donorCluster +
            params.holeElectronScatteringFactor * electrons / functionF;
    result.effectiveScatteringConcentration =
        effectiveScatteringConcentration / concentrationScale;
    result.screeningParameter = screening;
    result.screeningF = functionF;
    result.screeningG = functionG;

    if (effectiveScatteringConcentration <= 0.0) {
        result.mobility = result.latticeMobility;
        result.scatteringMobility = std::numeric_limits<Real>::infinity();
        return result;
    }

    const Real muN = muMax * muMax / (muMax - muMin) *
        std::pow(temperatureRatio, 3.0 * carrierParams.alpha - 1.5);
    const Real muC = muMax * muMin / (muMax - muMin) *
        std::sqrt(1.0 / temperatureRatio);
    const Real scatteringMobility =
        muN * (scatteringConcentration / effectiveScatteringConcentration) *
            std::pow(nRef / scatteringConcentration, carrierParams.alpha) +
        muC * (freeCarriers / effectiveScatteringConcentration);
    const Real bulkMobility = 1.0 /
        (1.0 / latticeMobility + 1.0 / scatteringMobility);
    result.scatteringMobility = scatteringMobility / mobilityScale;
    result.mobility = bulkMobility / mobilityScale;
    return result;
}

PhuMobLogDensityDerivatives evaluatePhuMobLogDensityDerivatives(
    CarrierType carrier, const PhuMobScalarState& state, const PhuMobParameters& params)
{
    const auto value = evaluatePhuMobScalar(carrier, state, params);
    PhuMobLogDensityDerivatives result;
    if (!std::isfinite(value.screeningParameter) ||
        value.effectiveScatteringConcentration <= 0.0)
        return result;
    const bool electron = carrier == CarrierType::Electron;
    const auto& cp = electron
        ? (params.donorSpecies == PhuMobDonorSpecies::Arsenic
            ? params.electronArsenic : params.electronPhosphorus)
        : params.holeBoron;
    const Real mass = electron ? params.electronMassRatio : params.holeMassRatio;
    const Real otherMass = electron ? params.holeMassRatio : params.electronMassRatio;
    const Real cs = params.internalConcentrationToCm3;
    const Real ms = params.internalMobilityToCm2PerVS;
    const Real n = state.electrons * cs, p = state.holes * cs;
    const Real nsc = value.scatteringConcentration * cs;
    const Real neff = value.effectiveScatteringConcentration * cs;
    const Real P = value.screeningParameter, F = value.screeningF;
    const Real temperatureRatio = state.temperature_K / 300.0;
    const Real aa = params.conwellWeisskopfFactor / 3.97e13;
    const Real bb = params.brooksHerringFactor / (1.36e20 * mass);
    const Real ratio = mass / otherMass;
    const Real denominator = std::pow(P, .6478) + 2.3670 - .8552 * ratio;
    // Cancel the common P^0.6478 terms algebraically, not numerically.
    const Real dFdP = .6478 * std::pow(P, -.3522) *
        (.7643 * (2.3670 - .8552 * ratio) - 2.2999 - 6.5502 * ratio) /
        (denominator * denominator);
    const Real oppositeImpurity = electron
        ? clusteredConcentration(state.acceptors * cs,
              params.acceptorClusterReference * cs, params.acceptorClusterCoefficient)
        : clusteredConcentration(state.donors * cs,
              params.donorClusterReference * cs, params.donorClusterCoefficient);
    const Real oppositePopulation = electron ? p : n;
    const Real scatteringFactor = electron ? params.electronHoleScatteringFactor
                                          : params.holeElectronScatteringFactor;
    const Real muMax = cp.muMax * ms, muMin = cp.muMin * ms;
    const Real muN = muMax * muMax / (muMax - muMin) *
        std::pow(temperatureRatio, 3.0 * cp.alpha - 1.5);
    const Real muC = muMax * muMin / (muMax - muMin) / std::sqrt(temperatureRatio);
    const Real first = muN * nsc * std::pow(cp.nRef * cs / nsc, cp.alpha);
    const Real scattering = value.scatteringMobility * ms;
    const Real lattice = value.latticeMobility * ms;
    const Real latticeWeight = lattice / (lattice + scattering);
    auto partial = [&](Real population, bool opposite) {
        const Real dnsc = opposite ? population : 0.0;
        const Real dP = -P * P / (temperatureRatio * temperatureRatio) *
            (aa * (2.0 / 3.0) * std::pow(nsc, -1.0 / 3.0) * dnsc + bb * population);
        const Real dneff = (value.screeningGDerivative * oppositeImpurity -
            scatteringFactor * oppositePopulation * dFdP / (F * F)) * dP +
            scatteringFactor * dnsc / F;
        const Real dscattering = (first * (1.0 - cp.alpha) / nsc * dnsc +
            muC * population - scattering * dneff) / neff;
        return latticeWeight * latticeWeight * dscattering / ms;
    };
    result.electrons = partial(n, !electron);
    result.holes = partial(p, electron);
    return result;
}

} // namespace vela
