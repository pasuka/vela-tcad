#include "vela/physics/MobilityModel.h"
#include "vela/physics/ExtendedPhuMobEvaluation.h"
#include "vela/physics/SlotIALMobility.h"
#include "vela/equation/IalTransport.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <array>
#include <cstdlib>
#include <string_view>

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
    if (value != "edge_projection" && value != "transport_cell_vector" &&
        value != "element_vertex_partial_layer" && value != "element_box_vertex") {
        throw std::invalid_argument("Unknown mobility.high_field_gradient_discretization: " + value);
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
                                        Real,
                                        MobilityDopants) const
{
    return material.mun;
}

Real ConstantMobility::holeMobility(const Material& material,
                                    Real,
                                    Real,
                                    Real,
                                    Real,
                                    Real,
                                    Real,
                                        MobilityDopants) const
{
    return material.mup;
}

DopingDependentMobility::DopingDependentMobility(MobilityModelConfig config)
    : config_(std::move(config))
    , bulkLaw_(isPhuMobModel(config_) ? BulkLaw::PhuMob :
               config_.model.rfind("constant_", 0) == 0 ? BulkLaw::Constant :
               isMasettiModel(config_.model) ? BulkLaw::Masetti : BulkLaw::CaugheyThomas)
    , surfaceLaw_(isLombardiModel(config_.model) ? SurfaceLaw::Lombardi :
                  isSurfaceMobilityModel(config_) ? SurfaceLaw::Generic : SurfaceLaw::None)
    , fieldMobility_(isFieldMobilityModel(config_.model))
{}

Real DopingDependentMobility::electronMobility(const Material& material,
                                               Real netDoping,
                                               Real n,
                                               Real p,
                                               Real electricField,
                                               Real surfaceNormalField,
                                               Real surfaceDistance,
                                               MobilityDopants) const
{
    if (bulkLaw_ == BulkLaw::PhuMob) {
        return electronMobilityWithIonizedImpurities(
            material, std::max<Real>(netDoping, 0.0),
            std::max<Real>(-netDoping, 0.0), n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    Real mobility = bulkLaw_ == BulkLaw::Constant
        ? material.mun
        : bulkLaw_ == BulkLaw::Masetti
            ? masetti(netDoping, config_.electronMasetti)
            : caugheyThomas(material.mun, netDoping, config_.electronCT);
    // Enormal contributes to the low-field mobility consumed by Canali HFS.
    // Applying the surface resistance after HFS changes both the constitutive
    // value and its field derivatives (unlike the PhuMob path below).
    if (surfaceLaw_ != SurfaceLaw::None)
        mobility = surfaceLaw_ == SurfaceLaw::Lombardi
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance,
                            material.temperature_K.value_or(300.0),
                            CarrierType::Electron,
                            config_.electronLombardi)
            : surfaceLimit(mobility, surfaceNormalField,
                           config_.surface.thetaElectron, config_.surface);
    if (fieldMobility_)
        mobility = fieldLimit(mobility, electricField, config_.electronField);
    return mobility;
}

Real DopingDependentMobility::holeMobility(const Material& material,
                                           Real netDoping,
                                           Real n,
                                           Real p,
                                           Real electricField,
                                           Real surfaceNormalField,
                                           Real surfaceDistance,
                                               MobilityDopants) const
{
    if (bulkLaw_ == BulkLaw::PhuMob) {
        return holeMobilityWithIonizedImpurities(
            material, std::max<Real>(netDoping, 0.0),
            std::max<Real>(-netDoping, 0.0), n, p, electricField,
            surfaceNormalField, surfaceDistance);
    }
    Real mobility = bulkLaw_ == BulkLaw::Constant
        ? material.mup
        : bulkLaw_ == BulkLaw::Masetti
            ? masetti(netDoping, config_.holeMasetti)
            : caugheyThomas(material.mup, netDoping, config_.holeCT);
    // Use the same low-field surface -> high-field ordering for holes.
    if (surfaceLaw_ != SurfaceLaw::None)
        mobility = surfaceLaw_ == SurfaceLaw::Lombardi
            ? lombardiLimit(mobility, netDoping, n, p, surfaceNormalField,
                            surfaceDistance,
                            material.temperature_K.value_or(300.0),
                            CarrierType::Hole,
                            config_.holeLombardi)
            : surfaceLimit(mobility, surfaceNormalField,
                           config_.surface.thetaHole, config_.surface);
    if (fieldMobility_)
        mobility = fieldLimit(mobility, electricField, config_.holeField);
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
    if (bulkLaw_ != BulkLaw::PhuMob) {
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
    if (surfaceLaw_ == SurfaceLaw::Lombardi) {
        mobility = lombardiLimit(
            mobility, donors + acceptors, n, p, surfaceNormalField,
            surfaceDistance, state.temperature_K, CarrierType::Electron,
            config_.electronLombardi);
    }
    if (fieldMobility_)
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
    if (bulkLaw_ != BulkLaw::PhuMob) {
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
    if (surfaceLaw_ == SurfaceLaw::Lombardi) {
        mobility = lombardiLimit(
            mobility, donors + acceptors, n, p, surfaceNormalField,
            surfaceDistance, state.temperature_K, CarrierType::Hole,
            config_.holeLombardi);
    }
    if (fieldMobility_)
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

namespace {
bool validLombardiParameters(const LombardiParameters& params, Real temperature)
{
    return !(params.B <= 0.0 || params.C < 0.0 || params.N0 <= 0.0 ||
        params.N1 <= 0.0 || params.N2 < 0.0 || params.delta <= 0.0 ||
        params.eta <= 0.0 || params.criticalLength <= 0.0 ||
        params.nu <= 0.0 || params.acousticFactor < 0.0 ||
        params.roughnessFactor < 0.0 || temperature <= 0.0 ||
        !std::isfinite(temperature));
}

Real lombardiDamping(Real distance, Real criticalLength)
{
    // A critical length greater than 100 cm disables distance damping.
    return criticalLength > 1.0 ? 1.0 : std::exp(-distance / criticalLength);
}

Real lombardiAcousticCoefficient(Real doping, Real temperature, const LombardiParameters& p)
{
    return p.C * std::pow(temperature / 300.0, -p.k) *
        std::pow((doping + p.N2) / p.N0, p.lambda);
}
}

DopingDependentMobility::PreparedMasettiLombardi
DopingDependentMobility::prepareMasettiLombardi(
    const Material& material, Real impurity, CarrierType carrier) const
{
    const bool electron = carrier == CarrierType::Electron;
    const auto& params = electron ? config_.electronLombardi : config_.holeLombardi;
    const Real temperature = material.temperature_K.value_or(300.0);
    PreparedMasettiLombardi out;
    out.bulkMobility = masetti(impurity, electron ? config_.electronMasetti : config_.holeMasetti);
    out.validSurfaceParameters = validLombardiParameters(params, temperature);
    if (out.validSurfaceParameters && out.bulkMobility > 0.0) {
        const Real doping = std::abs(impurity) * config_.internalConcentrationToM3;
        out.acousticCoefficient = lombardiAcousticCoefficient(doping, temperature, params);
        // This vertex-box path passes zero for both carrier populations.
        out.roughnessExponent = params.A +
            params.alpha * (0.0 + params.aOther * 0.0) /
            std::pow(doping + params.N1, params.nu);
    }
    return out;
}

Real DopingDependentMobility::evaluatePreparedMasettiLombardi(
    const PreparedMasettiLombardi& prepared, CarrierType carrier, Real electricField,
    Real surfaceNormalField, Real surfaceDistance) const
{
    const bool electron = carrier == CarrierType::Electron;
    Real mobility = prepared.bulkMobility;
    if (mobility > 0.0 && std::isfinite(surfaceNormalField) && std::isfinite(surfaceDistance)) {
        if (!prepared.validSurfaceParameters)
            throw std::invalid_argument("DopingDependentMobility: invalid Enhanced Lombardi parameters.");
        const auto& params = electron ? config_.electronLombardi : config_.holeLombardi;
        const Real field = std::abs(surfaceNormalField) * config_.internalFieldToVPerM;
        if (!(field <= 0.0)) {
            const Real distance = std::max(0.0, surfaceDistance) * config_.internalLengthToM;
            const Real damping = lombardiDamping(distance, params.criticalLength);
            if (!(damping <= std::numeric_limits<Real>::min()))
                mobility = lombardiMobilityFromField(mobility, field, damping,
                    prepared.acousticCoefficient, prepared.roughnessExponent, params);
        }
    }
    return fieldMobility_ ? fieldLimit(mobility, electricField,
        electron ? config_.electronField : config_.holeField) : mobility;
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
    if (!validLombardiParameters(params, latticeTemperature_K)) {
        throw std::invalid_argument(
            "DopingDependentMobility: invalid Enhanced Lombardi parameters.");
    }

    const Real field = std::abs(surfaceNormalField) * config_.internalFieldToVPerM;
    if (field <= 0.0)
        return bulkMobility;
    const Real distance = std::max<Real>(0.0, surfaceDistance) *
        config_.internalLengthToM;
    const Real damping = lombardiDamping(distance, params.criticalLength);
    if (damping <= std::numeric_limits<Real>::min())
        return bulkMobility;

    const Real doping = std::abs(netDoping) *
        config_.internalConcentrationToM3;
    const Real electronDensity = std::max<Real>(0.0, n) *
        config_.internalConcentrationToM3;
    const Real holeDensity = std::max<Real>(0.0, p) *
        config_.internalConcentrationToM3;
    const Real sameCarrier = carrier == CarrierType::Electron
        ? electronDensity : holeDensity;
    const Real otherCarrier = carrier == CarrierType::Electron
        ? holeDensity : electronDensity;
    const Real exponent = params.A +
        params.alpha * (sameCarrier + params.aOther * otherCarrier) /
        std::pow(doping + params.N1, params.nu);
    return lombardiMobilityFromField(bulkMobility, field, damping,
        lombardiAcousticCoefficient(doping, latticeTemperature_K, params), exponent, params);
}

Real DopingDependentMobility::lombardiMobilityFromField(
    Real bulkMobility, Real field, Real damping, Real acousticCoefficient,
    Real exponent, const LombardiParameters& params) const
{
    const Real acousticMobility = params.B / field + acousticCoefficient / std::cbrt(field);
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
    if (isIALMobilityModel(config)) config.dopingConcentrationBasis = "total_impurity";
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
        if (isIALMobilityModel(config)) config.dopingConcentrationBasis = "total_impurity";
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
    if (isIALMobilityModel(config)) config.dopingConcentrationBasis = "total_impurity";
    config.edgeAveraging = value.value("edge_averaging", config.edgeAveraging);
    config.prepareBoxMasettiLombardi =
        value.value("prepare_box_masetti_lombardi", config.prepareBoxMasettiLombardi);
    if (value.contains("constant_materials")) {
        config.constantMaterials = value.at("constant_materials").get<std::vector<std::string>>();
        if (config.edgeAveraging != "element_box_phumob")
            throw std::invalid_argument("constant_materials requires element_box_phumob");
    }
    if (config.edgeAveraging != "auto" && config.edgeAveraging != "legacy" &&
        config.edgeAveraging != "element_box" && config.edgeAveraging != "element_box_phumob" &&
        config.edgeAveraging != "element_box_masetti")
        throw std::invalid_argument("Unknown mobility.edge_averaging: " + config.edgeAveraging);
    if (config.edgeAveraging == "element_box_masetti" && config.model != "masetti_field_lombardi")
        throw std::invalid_argument("element_box_masetti requires masetti_field_lombardi.");
    if (config.edgeAveraging == "element_box" && config.model != "masetti" &&
        config.model != "constant" && config.model != "masetti_field")
        throw std::invalid_argument("element_box mobility supports constant and Masetti with optional element high field.");
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
    config.ialmobTransportFeedbackDerivatives=value.value(
        "ialmob_transport_feedback_derivatives",true);
    if(!config.ialmobTransportFeedbackDerivatives && config.model!="ialmob")
        throw std::invalid_argument("ialmob_transport_feedback_derivatives=false requires ialmob");
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
    if (isIALMobilityModel(config) && value.contains("ialmob"))
        parseIALMobilityParameters(value.at("ialmob"), config.slotIalMobility);

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
        config.ialmob = ialTransportOptionsFromJson(value.at("ialmob"), &value, scaling);
    } else if (value.contains("ialmob") && !isIALMobilityModel(config)) {
        throw std::invalid_argument("ialmob block requires model=ialmob");
    }
    const bool masettiElement = config.model == "masetti_field" && config.edgeAveraging == "element_box";
    if ((config.highFieldGradientDiscretization == "element_box_vertex") != masettiElement ||
        (masettiElement && (config.highFieldDrivingForce != "quasi_fermi_gradient" ||
                           !config.jacobianFieldDerivatives ||
                           config.carrierCurrentDiscretization != "scharfetter_gummel_edge")))
        throw std::invalid_argument("element_box_vertex requires Masetti field mobility, element_box, SG edge currents and full GradQF derivatives.");
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

bool usesVertexSampledMasetti(const MobilityModelConfig& config)
{
    return isMasettiModel(config.model) && !mobilityDependsOnCarrierDensity(config) &&
        (config.edgeAveraging == "element_box_masetti" ||
         (config.edgeAveraging == "auto" &&
          config.dopingConcentrationBasis != "cell_reconstructed_total_impurity"));
}

bool isSurfaceMobilityModel(const MobilityModelConfig& config)
{
    return config.model == "constant_lombardi" ||
           config.model == "constant_field_lombardi" ||
           config.model == "caughey_thomas_surface" ||
           config.model == "caughey_thomas_field_surface" ||
           config.model == "masetti_surface" ||
           config.model == "masetti_field_surface" ||
           isLombardiModel(config.model) || isIALMobilityModel(config);
}

bool isIALMobilityModel(const MobilityModelConfig& config)
{
    return config.model == "slot_ialmob" || config.model == "ialmob_field";
}

bool mobilityDependsOnCarrierDensity(const MobilityModelConfig& config)
{
    return isIALMobilityModel(config) || (isLombardiModel(config.model) &&
        ((config.electronLombardi.alpha != 0.0 && config.electronLombardi.roughnessFactor != 0.0) ||
         (config.holeLombardi.alpha != 0.0 && config.holeLombardi.roughnessFactor != 0.0)));
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
    Real electronMobility(const Material&,Real,Real,Real,Real,Real,Real,MobilityDopants) const override {
        throw std::logic_error("IALMob requires the live element transport interface");
    }
    Real holeMobility(const Material&,Real,Real,Real,Real,Real,Real,MobilityDopants) const override {
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
    if (isIALMobilityModel(config)) return std::make_unique<IALMobility>(config);
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

template<class Scalar>
void validatePhuMobCarrier(const BasicPhuMobCarrierParameters<Scalar>& carrier,
                           const char* name)
{
    if (!std::isfinite(carrier.muMax) || !std::isfinite(carrier.muMin) ||
        !std::isfinite(carrier.theta) || !std::isfinite(carrier.nRef) ||
        !std::isfinite(carrier.alpha) || carrier.muMax <= carrier.muMin ||
        carrier.muMin < 0 || carrier.theta < 0 || carrier.nRef <= 0 ||
        carrier.alpha <= 0) {
        throw std::invalid_argument(
            std::string("PhuMob: invalid ") + name + " carrier parameters.");
    }
}

template<class Scalar>
void validatePhuMob(const BasicPhuMobParameters<Scalar>& params,
                    const BasicPhuMobScalarState<Scalar>& state)
{
    validatePhuMobCarrier(params.electronArsenic, "electron-arsenic");
    validatePhuMobCarrier(params.electronPhosphorus, "electron-phosphorus");
    validatePhuMobCarrier(params.holeBoron, "hole-boron");
    const Scalar values[] = {
        params.donorClusterReference, params.acceptorClusterReference,
        params.donorClusterCoefficient, params.acceptorClusterCoefficient,
        params.electronMassRatio, params.holeMassRatio,
        params.conwellWeisskopfFactor, params.brooksHerringFactor,
        params.electronHoleScatteringFactor,
        params.holeElectronScatteringFactor, params.gA, params.gB, params.gC,
        params.gAlpha, params.gAlphaPrime, params.gBeta, params.gGamma,
        params.internalConcentrationToCm3, params.internalMobilityToCm2PerVS,
    };
    for (const Scalar value : values) {
        if (!std::isfinite(value))
            throw std::invalid_argument("PhuMob: parameters must be finite.");
    }
    if (params.donorClusterReference <= 0 ||
        params.acceptorClusterReference <= 0 ||
        params.donorClusterCoefficient < 0 ||
        params.acceptorClusterCoefficient < 0 ||
        params.electronMassRatio <= 0 || params.holeMassRatio <= 0 ||
        params.conwellWeisskopfFactor <= 0 ||
        params.brooksHerringFactor <= 0 ||
        params.electronHoleScatteringFactor < 0 ||
        params.holeElectronScatteringFactor < 0 || params.gA <= 0 ||
        params.gB <= 0 || params.gC <= 0 || params.gAlpha <= 0 ||
        params.gAlphaPrime <= 0 || params.gBeta <= 0 ||
        params.gGamma <= 0 || params.internalConcentrationToCm3 <= 0 ||
        params.internalMobilityToCm2PerVS <= 0) {
        throw std::invalid_argument("PhuMob: parameters must be physically positive.");
    }

    const Scalar stateValues[] = {
        state.donors, state.acceptors, state.electrons, state.holes,
        state.temperature_K,
    };
    for (const Scalar value : stateValues) {
        if (!std::isfinite(value))
            throw std::invalid_argument("PhuMob: scalar state must be finite.");
    }
    if (state.donors < 0 || state.acceptors < 0 ||
        state.electrons < 0 || state.holes < 0 ||
        state.temperature_K <= 0) {
        throw std::invalid_argument(
            "PhuMob: concentrations must be nonnegative and temperature positive.");
    }
}

template<class Scalar>
Scalar clusteredConcentration(Scalar concentration, Scalar reference, Scalar coefficient)
{
    if (concentration <= 0)
        return 0;
    const Scalar ratio = concentration / reference;
    const Scalar ratioSquared = ratio * ratio;
    const Scalar correction = ratioSquared /
        (1 + coefficient * ratioSquared);
    return concentration * (1 + correction);
}

template<class Scalar>
Scalar screeningF(Scalar screening, Scalar massRatio, Scalar otherMassRatio)
{
    if (!std::isfinite(screening))
        return Scalar(0.7643L);
    const Scalar power = std::pow(screening, Scalar(0.6478L));
    const Scalar massQuotient = massRatio / otherMassRatio;
    return (Scalar(0.7643L) * power + Scalar(2.2999L) + Scalar(6.5502L) * massQuotient) /
        (power + Scalar(2.3670L) - Scalar(0.8552L) * massQuotient);
}

template<class Scalar>
Scalar rawScreeningG(Scalar screening, Scalar massRatio, Scalar temperature_K,
                   const BasicPhuMobParameters<Scalar>& params)
{
    if (!std::isfinite(screening))
        return 1;
    const Scalar temperatureRatio = temperature_K / 300;
    const Scalar firstScale = std::pow(
        temperatureRatio / massRatio, params.gAlpha);
    const Scalar secondScale = std::pow(
        massRatio / temperatureRatio, params.gAlphaPrime);
    return 1 - params.gA * std::pow(
        params.gB + screening * firstScale, -params.gBeta) +
        params.gC * std::pow(screening * secondScale, -params.gGamma);
}

template<class Scalar>
bool screeningGDerivativeIsPositive(Scalar logScreening, Scalar massRatio,
                                    Scalar temperature_K,
                                    const BasicPhuMobParameters<Scalar>& params)
{
    const Scalar screening = std::exp(logScreening);
    const Scalar temperatureRatio = temperature_K / 300;
    const Scalar firstScale = std::pow(
        temperatureRatio / massRatio, params.gAlpha);
    const Scalar secondScale = std::pow(
        massRatio / temperatureRatio, params.gAlphaPrime);
    const Scalar logPositive = std::log(
        params.gA * params.gBeta * firstScale) -
        (params.gBeta + 1) * std::log(
            params.gB + screening * firstScale);
    const Scalar logNegative = std::log(params.gC * params.gGamma) -
        params.gGamma * std::log(secondScale) -
        (params.gGamma + 1) * logScreening;
    return logPositive > logNegative;
}

// The clamp minimum depends on temperature, mass and all seven G constants,
// never on carrier/dopant populations. Keep exact (not quantized) keys and a
// bounded thread-local cache: no shared mutation and no stale state across
// devices or temperature sweeps. The miss path preserves the original search.
template<class Scalar>
std::array<Scalar,2> screeningGMinimum(Scalar massRatio, Scalar temperature_K,
                                   const BasicPhuMobParameters<Scalar>& params)
{
    using Key=std::array<Scalar,9>;
    struct Entry { Key key{}; std::array<Scalar,2> value{}; bool valid=false; };
    static thread_local std::array<Entry,8> cache{};
    static thread_local std::size_t next=0;
    // Diagnostic rollback for paired numerical/performance verification.
    static const bool enabled=[](){const auto* flag=std::getenv("VELA_PHUMOB_G_MINIMUM_CACHE");return !flag || std::string_view(flag)!="0";}();
    const Key key{massRatio,temperature_K,params.gA,params.gB,params.gC,
                  params.gAlpha,params.gAlphaPrime,params.gBeta,params.gGamma};
    if(enabled) for(const auto& entry:cache)
        if(entry.valid && entry.key==key) return entry.value;
    const Scalar logFloor = std::log(std::numeric_limits<Scalar>::min());
    Scalar low = std::max(logFloor, Scalar(-80));
    Scalar high = 0;
    while (!screeningGDerivativeIsPositive(
               high, massRatio, temperature_K, params) && high < 80) {
        high += 4;
    }
    for (int iteration = 0; iteration < 80; ++iteration) {
        const Scalar middle = Scalar(0.5L) * (low + high);
        if (screeningGDerivativeIsPositive(
                middle, massRatio, temperature_K, params)) {
            high = middle;
        } else {
            low = middle;
        }
    }
    const Scalar minimumScreening = std::exp(Scalar(0.5L) * (low + high));
    const Scalar minimum = rawScreeningG(
        minimumScreening, massRatio, temperature_K, params);
    const Scalar minimumFloor = std::abs(minimum);
    const std::array<Scalar,2> result{minimumScreening,minimumFloor};
    if(enabled){cache[next]={key,result,true};next=(next+1)%cache.size();}
    return result;
}

template<class Scalar>
Scalar screeningG(Scalar screening, Scalar massRatio, Scalar temperature_K,
                const BasicPhuMobParameters<Scalar>& params, Scalar* derivative)
{
    *derivative = 0;
    if (!std::isfinite(screening)) return 1;
    const auto minimum=screeningGMinimum(massRatio,temperature_K,params);
    const Scalar minimumScreening=minimum[0],minimumFloor=minimum[1];
    if (screening < minimumScreening)
        return minimumFloor;
    const Scalar raw = rawScreeningG(screening, massRatio, temperature_K, params);
    if (raw > minimumFloor) {
        const Scalar firstScale = std::pow(temperature_K / (300 * massRatio), params.gAlpha);
        const Scalar secondScale = std::pow(300 * massRatio / temperature_K, params.gAlphaPrime);
        *derivative = params.gA * params.gBeta * firstScale *
            std::pow(params.gB + screening * firstScale, -params.gBeta - 1) -
            params.gC * params.gGamma * std::pow(screening * secondScale, -params.gGamma) / screening;
    }
    return std::max(raw, minimumFloor);
}

template<class Scalar>
BasicPhuMobScalarResult<Scalar> evaluatePhuMobScalarImpl(
    CarrierType carrier,
    const BasicPhuMobScalarState<Scalar>& state,
    const BasicPhuMobParameters<Scalar>& params)
{
    validatePhuMob(params, state);
    const auto& carrierParams = carrier == CarrierType::Electron
        ? (params.donorSpecies == PhuMobDonorSpecies::Arsenic
               ? params.electronArsenic
               : params.electronPhosphorus)
        : params.holeBoron;
    const Scalar massRatio = carrier == CarrierType::Electron
        ? params.electronMassRatio : params.holeMassRatio;
    const Scalar otherMassRatio = carrier == CarrierType::Electron
        ? params.holeMassRatio : params.electronMassRatio;

    const Scalar concentrationScale = params.internalConcentrationToCm3;
    const Scalar mobilityScale = params.internalMobilityToCm2PerVS;
    const Scalar donors = state.donors * concentrationScale;
    const Scalar acceptors = state.acceptors * concentrationScale;
    const Scalar electrons = state.electrons * concentrationScale;
    const Scalar holes = state.holes * concentrationScale;
    const Scalar donorReference = params.donorClusterReference * concentrationScale;
    const Scalar acceptorReference =
        params.acceptorClusterReference * concentrationScale;
    const Scalar donorCluster = clusteredConcentration(
        donors, donorReference, params.donorClusterCoefficient);
    const Scalar acceptorCluster = clusteredConcentration(
        acceptors, acceptorReference, params.acceptorClusterCoefficient);

    const Scalar muMax = carrierParams.muMax * mobilityScale;
    const Scalar muMin = carrierParams.muMin * mobilityScale;
    const Scalar nRef = carrierParams.nRef * concentrationScale;
    const Scalar temperatureRatio = state.temperature_K / 300;
    const Scalar latticeMobility = muMax * std::pow(
        temperatureRatio, -carrierParams.theta);

    const Scalar otherCarrier = carrier == CarrierType::Electron
        ? holes : electrons;
    const Scalar scatteringConcentration =
        donorCluster + acceptorCluster + otherCarrier;

    BasicPhuMobScalarResult<Scalar> result;
    result.latticeMobility = latticeMobility / mobilityScale;
    result.scatteringConcentration = scatteringConcentration / concentrationScale;
    if (scatteringConcentration <= 0) {
        result.mobility = result.latticeMobility;
        result.scatteringMobility = std::numeric_limits<Scalar>::infinity();
        result.effectiveScatteringConcentration = 0;
        result.screeningParameter = std::numeric_limits<Scalar>::infinity();
        result.screeningF = Scalar(0.7643L);
        result.screeningG = 1;
        return result;
    }

    const Scalar freeCarriers = electrons + holes;
    const Scalar screeningDenominator =
        params.conwellWeisskopfFactor *
            std::pow(scatteringConcentration, Scalar(2) / 3) / Scalar(3.97e13L) +
        params.brooksHerringFactor * freeCarriers /
            (Scalar(1.36e20L) * massRatio);
    const Scalar screening = screeningDenominator > 0
        ? temperatureRatio * temperatureRatio / screeningDenominator
        : std::numeric_limits<Scalar>::infinity();
    const Scalar functionF = screeningF(screening, massRatio, otherMassRatio);
    const Scalar functionG = screeningG(
        screening, massRatio, state.temperature_K, params, &result.screeningGDerivative);

    const Scalar effectiveScatteringConcentration = carrier == CarrierType::Electron
        ? donorCluster + functionG * acceptorCluster +
            params.electronHoleScatteringFactor * holes / functionF
        : acceptorCluster + functionG * donorCluster +
            params.holeElectronScatteringFactor * electrons / functionF;
    result.effectiveScatteringConcentration =
        effectiveScatteringConcentration / concentrationScale;
    result.screeningParameter = screening;
    result.screeningF = functionF;
    result.screeningG = functionG;

    if (effectiveScatteringConcentration <= 0) {
        result.mobility = result.latticeMobility;
        result.scatteringMobility = std::numeric_limits<Scalar>::infinity();
        return result;
    }

    const Scalar muN = muMax * muMax / (muMax - muMin) *
        std::pow(temperatureRatio, 3 * carrierParams.alpha - Scalar(1.5L));
    const Scalar muC = muMax * muMin / (muMax - muMin) *
        std::sqrt(1 / temperatureRatio);
    const Scalar scatteringMobility =
        muN * (scatteringConcentration / effectiveScatteringConcentration) *
            std::pow(nRef / scatteringConcentration, carrierParams.alpha) +
        muC * (freeCarriers / effectiveScatteringConcentration);
    const Scalar bulkMobility = 1 /
        (1 / latticeMobility + 1 / scatteringMobility);
    result.scatteringMobility = scatteringMobility / mobilityScale;
    result.mobility = bulkMobility / mobilityScale;
    return result;
}

template<class Scalar>
BasicPhuMobLogDensityDerivatives<Scalar> evaluatePhuMobLogDensityDerivativesImpl(
    CarrierType carrier, const BasicPhuMobScalarState<Scalar>& state,
    const BasicPhuMobParameters<Scalar>& params, const BasicPhuMobScalarResult<Scalar>& value)
{
    BasicPhuMobLogDensityDerivatives<Scalar> result;
    if (!std::isfinite(value.screeningParameter) ||
        value.effectiveScatteringConcentration <= 0)
        return result;
    const bool electron = carrier == CarrierType::Electron;
    const auto& cp = electron
        ? (params.donorSpecies == PhuMobDonorSpecies::Arsenic
            ? params.electronArsenic : params.electronPhosphorus)
        : params.holeBoron;
    const Scalar mass = electron ? params.electronMassRatio : params.holeMassRatio;
    const Scalar otherMass = electron ? params.holeMassRatio : params.electronMassRatio;
    const Scalar cs = params.internalConcentrationToCm3;
    const Scalar ms = params.internalMobilityToCm2PerVS;
    const Scalar n = state.electrons * cs, p = state.holes * cs;
    const Scalar nsc = value.scatteringConcentration * cs;
    const Scalar neff = value.effectiveScatteringConcentration * cs;
    const Scalar P = value.screeningParameter, F = value.screeningF;
    const Scalar temperatureRatio = state.temperature_K / 300;
    const Scalar aa = params.conwellWeisskopfFactor / Scalar(3.97e13L);
    const Scalar bb = params.brooksHerringFactor / (Scalar(1.36e20L) * mass);
    const Scalar ratio = mass / otherMass;
    const Scalar denominator = std::pow(P, Scalar(.6478L)) + Scalar(2.3670L) - Scalar(.8552L) * ratio;
    // Cancel the common P^0.6478 terms algebraically, not numerically.
    const Scalar dFdP = Scalar(.6478L) * std::pow(P, -Scalar(.3522L)) *
        (Scalar(.7643L) * (Scalar(2.3670L) - Scalar(.8552L) * ratio) - Scalar(2.2999L) - Scalar(6.5502L) * ratio) /
        (denominator * denominator);
    const Scalar oppositeImpurity = electron
        ? clusteredConcentration(state.acceptors * cs,
              params.acceptorClusterReference * cs, params.acceptorClusterCoefficient)
        : clusteredConcentration(state.donors * cs,
              params.donorClusterReference * cs, params.donorClusterCoefficient);
    const Scalar oppositePopulation = electron ? p : n;
    const Scalar scatteringFactor = electron ? params.electronHoleScatteringFactor
                                          : params.holeElectronScatteringFactor;
    const Scalar muMax = cp.muMax * ms, muMin = cp.muMin * ms;
    const Scalar muN = muMax * muMax / (muMax - muMin) *
        std::pow(temperatureRatio, 3 * cp.alpha - Scalar(1.5L));
    const Scalar muC = muMax * muMin / (muMax - muMin) / std::sqrt(temperatureRatio);
    const Scalar first = muN * nsc * std::pow(cp.nRef * cs / nsc, cp.alpha);
    const Scalar scattering = value.scatteringMobility * ms;
    const Scalar lattice = value.latticeMobility * ms;
    const Scalar latticeWeight = lattice / (lattice + scattering);
    auto partial = [&](Scalar population, bool opposite) {
        const Scalar dnsc = opposite ? population : 0;
        const Scalar dP = -P * P / (temperatureRatio * temperatureRatio) *
            (aa * (Scalar(2) / 3) * std::pow(nsc, Scalar(-1) / 3) * dnsc + bb * population);
        const Scalar dneff = (value.screeningGDerivative * oppositeImpurity -
            scatteringFactor * oppositePopulation * dFdP / (F * F)) * dP +
            scatteringFactor * dnsc / F;
        const Scalar dscattering = (first * (1 - cp.alpha) / nsc * dnsc +
            muC * population - scattering * dneff) / neff;
        return latticeWeight * latticeWeight * dscattering / ms;
    };
    result.electrons = partial(n, !electron);
    result.holes = partial(p, electron);
    return result;
}

} // namespace

PhuMobScalarResult evaluatePhuMobScalar(
    CarrierType carrier, const PhuMobScalarState& state, const PhuMobParameters& params)
{
    return evaluatePhuMobScalarImpl(carrier, state, params);
}

PhuMobLogDensityDerivatives evaluatePhuMobLogDensityDerivatives(
    CarrierType carrier, const PhuMobScalarState& state, const PhuMobParameters& params)
{
    return evaluatePhuMobLogDensityDerivativesImpl(carrier, state, params,
        evaluatePhuMobScalarImpl(carrier, state, params));
}

ExtendedPhuMobScalarResult evaluateExtendedPhuMobScalar(
    CarrierType carrier, const ExtendedPhuMobScalarState& state, const ExtendedPhuMobParameters& params)
{
    return evaluatePhuMobScalarImpl(carrier, state, params);
}

ExtendedPhuMobLogDensityDerivatives evaluateExtendedPhuMobLogDensityDerivatives(
    CarrierType carrier, const ExtendedPhuMobScalarState& state, const ExtendedPhuMobParameters& params)
{
    return evaluateExtendedPhuMobWithLogDensityDerivatives(carrier, state, params).slopes;
}

ExtendedPhuMobEvaluation evaluateExtendedPhuMobWithLogDensityDerivatives(
    CarrierType carrier, const ExtendedPhuMobScalarState& state, const ExtendedPhuMobParameters& params)
{
    const auto value = evaluatePhuMobScalarImpl(carrier, state, params);
    return {value.mobility, evaluatePhuMobLogDensityDerivativesImpl(carrier, state, params, value)};
}

} // namespace vela
