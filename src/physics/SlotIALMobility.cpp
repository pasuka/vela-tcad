#include "vela/physics/SlotIALMobility.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace vela {
namespace {
constexpr Real infinity = std::numeric_limits<Real>::infinity();
using CarrierParameters = IALMobilityCarrierParameters;
using ParameterMember = std::pair<std::string_view, Real CarrierParameters::*>;

// Explicit source-model units avoid ambiguous fractional-power SI conversions.
constexpr ParameterMember carrierMembers[] = {
#define IAL_PARAMETER(name) {#name, &CarrierParameters::name}
    IAL_PARAMETER(mumax), IAL_PARAMETER(mumin), IAL_PARAMETER(theta),
    IAL_PARAMETER(n_ref), IAL_PARAMETER(alpha), IAL_PARAMETER(ndop_ref),
    IAL_PARAMETER(nsc_ref), IAL_PARAMETER(S), IAL_PARAMETER(p),
    IAL_PARAMETER(l_crit_c), IAL_PARAMETER(B), IAL_PARAMETER(C),
    IAL_PARAMETER(lambda), IAL_PARAMETER(k), IAL_PARAMETER(alpha_ph2d_A),
    IAL_PARAMETER(alpha_ph2d_D), IAL_PARAMETER(lambda_ph2d_A),
    IAL_PARAMETER(lambda_ph2d_D), IAL_PARAMETER(delta), IAL_PARAMETER(lambda_sr),
    IAL_PARAMETER(A), IAL_PARAMETER(alpha_sr), IAL_PARAMETER(nu),
    IAL_PARAMETER(eta), IAL_PARAMETER(N1), IAL_PARAMETER(N2),
    IAL_PARAMETER(alpha_sr_A), IAL_PARAMETER(alpha_sr_D),
    IAL_PARAMETER(lambda_sr_A), IAL_PARAMETER(lambda_sr_D), IAL_PARAMETER(l_crit),
    IAL_PARAMETER(D1_inv), IAL_PARAMETER(D2_inv), IAL_PARAMETER(nu0_inv),
    IAL_PARAMETER(nu1_inv), IAL_PARAMETER(nu2_inv), IAL_PARAMETER(alpha1_inv),
    IAL_PARAMETER(alpha2_inv), IAL_PARAMETER(D1_acc), IAL_PARAMETER(D2_acc),
    IAL_PARAMETER(nu0_acc), IAL_PARAMETER(nu1_acc), IAL_PARAMETER(nu2_acc),
    IAL_PARAMETER(alpha1_acc), IAL_PARAMETER(alpha2_acc), IAL_PARAMETER(EnormMinimum)
#undef IAL_PARAMETER
};

void validate(const IALMobilityParameters& parameters)
{
    for (const auto* value : {&parameters.electron, &parameters.hole}) {
        for (const auto& [name, member] : carrierMembers)
            if (!std::isfinite(value->*member))
                throw std::invalid_argument("IALMob: non-finite parameter " + std::string(name));
        if (!(value->mumax > value->mumin && value->mumin > 0.0) ||
            value->n_ref <= 0.0 || value->ndop_ref <= 0.0 || value->nsc_ref <= 0.0 ||
            value->alpha <= 0.0 || value->S < 0.0 || value->l_crit_c <= 0.0 ||
            value->l_crit <= 0.0 || value->B <= 0.0 || value->C < 0.0 ||
            value->delta <= 0.0 || value->eta <= 0.0 || value->N1 <= 0.0 ||
            value->N2 <= 0.0 || value->D1_inv < 0.0 || value->D2_inv <= 0.0 ||
            value->D1_acc < 0.0 || value->D2_acc <= 0.0 || value->EnormMinimum < 0.0 ||
            value->alpha_sr < 0.0 || value->nu0_inv < 0.0 || value->nu0_acc < 0.0 ||
            value->alpha_ph2d_A <= 0.0 || value->alpha_ph2d_D <= 0.0 ||
            value->alpha_sr_A <= 0.0 || value->alpha_sr_D <= 0.0)
            throw std::invalid_argument("IALMob: invalid mobility coefficients");
    }
    for (Real value : {parameters.nref_D, parameters.nref_A, parameters.cref_D,
                       parameters.cref_A, parameters.me_over_m0, parameters.mh_over_m0})
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("IALMob: clustering and mass parameters must be positive");
    if (parameters.phonon_combination < 0 || parameters.phonon_combination > 2)
        throw std::invalid_argument("IALMob: phonon_combination must be 0, 1, or 2");
}

Real rawG(Real screening, Real relativeTemperature, Real mass)
{
    if (std::isinf(screening)) return 1.0;
    return 1.0 - 0.89233 / std::pow(0.41372 + screening *
            std::pow(relativeTemperature / mass, 0.28227), 0.19778) +
        0.005978 / std::pow(screening *
            std::pow(mass / relativeTemperature, 0.72169), 1.80618);
}

Real screeningMinimum(Real relativeTemperature, Real mass)
{
    // Minimize in log(P), independent of n/p and fixed throughout a solve.
    Real left = -30.0, right = 30.0;
    constexpr Real fraction = 0.38196601125010515;
    for (int iteration = 0; iteration < 90; ++iteration) {
        const Real first = left + fraction * (right - left);
        const Real second = right - fraction * (right - left);
        if (rawG(std::exp(first), relativeTemperature, mass) <
            rawG(std::exp(second), relativeTemperature, mass)) right = second;
        else left = first;
    }
    return std::exp(0.5 * (left + right));
}

Real clustered(Real concentration, Real reference, Real coefficient)
{
    if (concentration == 0.0) return 0.0;
    const Real ratio = reference / concentration;
    return concentration * (1.0 + 1.0 / (coefficient + ratio * ratio));
}

Real coulomb2DPart(Real concentration, Real carriers, Real temperature,
                   const CarrierParameters& parameters, bool inversion)
{
    if (concentration == 0.0) return infinity;
    const Real doping = concentration / parameters.ndop_ref;
    const Real normalizedCarriers = carriers / parameters.nsc_ref;
    const Real first = (inversion ? parameters.D1_inv : parameters.D1_acc) *
        std::pow(temperature, inversion ? parameters.alpha1_inv : parameters.alpha1_acc) *
        std::pow(normalizedCarriers, inversion ? parameters.nu0_inv : parameters.nu0_acc) /
        std::pow(doping, inversion ? parameters.nu1_inv : parameters.nu1_acc);
    const Real second = (inversion ? parameters.D2_inv : parameters.D2_acc) *
        std::pow(temperature, inversion ? parameters.alpha2_inv : parameters.alpha2_acc) /
        std::pow(doping, inversion ? parameters.nu2_inv : parameters.nu2_acc);
    return std::hypot(first, second);
}

Real inverse(Real value) { return 1.0 / value; }

Real weightedMobility(Real first, Real second, Real secondWeight)
{
    if (secondWeight == 0.0) return first;
    if (secondWeight == 1.0) return second;
    return (1.0 - secondWeight) * first + secondWeight * second;
}
} // namespace

void parseIALMobilityParameters(const nlohmann::json& document, IALMobilityParameters& parameters)
{
    if (!document.is_object()) throw std::invalid_argument("mobility.ialmob must be an object");
    for (const auto& [key, value] : document.items()) {
        if (key == "electron" || key == "hole") {
            if (!value.is_object()) throw std::invalid_argument("IALMob carrier parameters must be an object");
            auto& target = key == "electron" ? parameters.electron : parameters.hole;
            for (const auto& [name, coefficient] : value.items()) {
                const auto found = std::find_if(std::begin(carrierMembers), std::end(carrierMembers),
                    [&](const auto& member) { return member.first == name; });
                if (found == std::end(carrierMembers))
                    throw std::invalid_argument("IALMob: unsupported carrier parameter '" + name + "'");
                target.*(found->second) = coefficient.get<Real>();
            }
        } else if (key == "full_phumob") parameters.full_phumob = value.get<bool>();
        else if (key == "clustering_everywhere") parameters.clustering_everywhere = value.get<bool>();
        else if (key == "phonon_combination") parameters.phonon_combination = value.get<int>();
        else if (key == "nref_D") parameters.nref_D = value.get<Real>();
        else if (key == "nref_A") parameters.nref_A = value.get<Real>();
        else if (key == "cref_D") parameters.cref_D = value.get<Real>();
        else if (key == "cref_A") parameters.cref_A = value.get<Real>();
        else if (key == "me_over_m0") parameters.me_over_m0 = value.get<Real>();
        else if (key == "mh_over_m0") parameters.mh_over_m0 = value.get<Real>();
        else throw std::invalid_argument("IALMob: unsupported option '" + key + "'");
    }
    validate(parameters);
}

IALMobility::IALMobility(MobilityModelConfig config) : config_(std::move(config))
{
    validate(config_.ialmob);
    if (config_.dopingConcentrationBasis != "total_impurity")
        throw std::invalid_argument("IALMob requires total_impurity with separate nodal donors and acceptors");
    electronScreeningMinimum_ = screeningMinimum(1.0, config_.ialmob.me_over_m0);
    holeScreeningMinimum_ = screeningMinimum(1.0, config_.ialmob.mh_over_m0);
}

IALMobilityComponents IALMobility::components(const Material& material, CarrierType carrier,
    MobilityDopants dopants, Real n, Real p, Real normal, Real distance) const
{
    for (Real value : {dopants.donors, dopants.acceptors, n, p})
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument("IALMob requires finite nonnegative n, p, donors, and acceptors");
    const Real temperature = material.temperature_K.value_or(300.0);
    if (!std::isfinite(temperature) || temperature < 50.0)
        throw std::invalid_argument("IALMob: temperatures below 50 K require unsupported low-temperature corrections");
    const Real relativeTemperature = temperature / 300.0;
    const bool electron = carrier == CarrierType::Electron;
    const auto& model = config_.ialmob;
    const auto& parameters = electron ? model.electron : model.hole;
    const Real concentrationFactor = config_.internalConcentrationToM3 * 1.0e-6;
    const Real donors = dopants.donors * concentrationFactor;
    const Real acceptors = dopants.acceptors * concentrationFactor;
    const Real carriers = (electron ? n : p) * concentrationFactor;
    const Real other = model.full_phumob ? (electron ? p : n) * concentrationFactor : 0.0;
    const Real mass = electron ? model.me_over_m0 : model.mh_over_m0;
    const Real otherMass = electron ? model.mh_over_m0 : model.me_over_m0;
    const Real clusteredDonors = clustered(donors, model.nref_D, model.cref_D);
    const Real clusteredAcceptors = clustered(acceptors, model.nref_A, model.cref_A);
    const Real totalScatterers = clusteredDonors + clusteredAcceptors + other;
    IALMobilityComponents result;
    // Eq. 306. The concentration reference includes cm^-2, hence Nsc^(2/3).
    const Real screeningDenominator = 2.459 * std::pow(totalScatterers, 2.0 / 3.0) / 3.97e13 +
        3.828 * (carriers + other) / (1.36e20 * mass);
    result.screening = screeningDenominator > 0.0 ?
        relativeTemperature * relativeTemperature / screeningDenominator : infinity;
    const Real power = std::pow(result.screening, 0.6478);
    result.F = std::isinf(power) ? 0.7643 :
        (0.7643 * power + 2.2999 + 6.5502 * mass / otherMass) /
        (power + 2.3670 - 0.8552 * mass / otherMass);
    Real minimum = electron ? electronScreeningMinimum_ : holeScreeningMinimum_;
    if (temperature != 300.0) minimum = screeningMinimum(relativeTemperature, mass);
    const Real screeningForG = model.full_phumob ? std::max(result.screening, minimum) : result.screening;
    result.G = rawG(screeningForG, relativeTemperature, mass);
    const Real accumulationStar = electron ? clusteredDonors : clusteredAcceptors;
    const Real inversionStar = electron ? clusteredAcceptors : clusteredDonors;
    const Real effectiveScatterers = accumulationStar + result.G * inversionStar + other / result.F;
    if (effectiveScatterers == 0.0) result.coulomb3D = infinity;
    else {
        const Real muN = parameters.mumax * parameters.mumax / (parameters.mumax - parameters.mumin) *
            std::pow(relativeTemperature, 3.0 * parameters.alpha - 1.5);
        const Real muC = parameters.mumax * parameters.mumin / (parameters.mumax - parameters.mumin) /
            std::sqrt(relativeTemperature);
        result.coulomb3D = muN * totalScatterers / effectiveScatterers *
            std::pow(parameters.n_ref / totalScatterers, parameters.alpha) +
            muC * (carriers + other) / effectiveScatterers;
    }
    const Real rawDonors = model.clustering_everywhere ? clusteredDonors : donors;
    const Real rawAcceptors = model.clustering_everywhere ? clusteredAcceptors : acceptors;
    const Real inversion = electron ? rawAcceptors : rawDonors;
    const Real accumulation = electron ? rawDonors : rawAcceptors;
    const Real muInv = coulomb2DPart(inversion, carriers, relativeTemperature, parameters, true);
    const Real muAcc = coulomb2DPart(accumulation, carriers, relativeTemperature, parameters, false) * result.G;
    result.coulomb2D = inverse(inverse(muInv) + inverse(muAcc));
    const bool hasInterface = std::isfinite(normal) && std::isfinite(distance);
    const Real field = hasInterface ? std::abs(normal) * config_.internalFieldToVPerM * 0.01 : 0.0;
    if (hasInterface) {
        const Real distanceCm = std::max(0.0, distance) * config_.internalLengthToM * 100.0;
        result.interfaceDamping = std::exp(-distanceCm / parameters.l_crit);
        result.coulombDamping = std::exp(-distanceCm / parameters.l_crit_c);
    }
    // No LayerThickness/ThinLayer: S_t=0 and erf((t+t1)/tcoulomb)=1.
    const Real mixingArgument = parameters.S * std::pow(field, 2.0 / 3.0) / temperature - parameters.p;
    if (mixingArgument >= 0.0) {
        const Real exponential = std::exp(-mixingArgument);
        result.bulkWeight = exponential / (1.0 + exponential);
    } else result.bulkWeight = 1.0 / (1.0 + std::exp(mixingArgument));
    result.coulomb = weightedMobility(result.coulomb3D, result.coulomb2D,
        result.coulombDamping * (1.0 - result.bulkWeight));
    result.phonon3D = parameters.mumax * std::pow(relativeTemperature, -parameters.theta);
    result.phonon2D = infinity;
    result.roughness = infinity;
    if (field > 0.0 && field >= parameters.EnormMinimum && result.interfaceDamping > 0.0) {
        const Real acousticDoping = parameters.alpha_ph2d_A *
            std::pow(rawAcceptors + parameters.N2 / 2.0, parameters.lambda_ph2d_A) +
            parameters.alpha_ph2d_D * std::pow(rawDonors + parameters.N2 / 2.0, parameters.lambda_ph2d_D);
        result.phonon2D = parameters.B / field + parameters.C * std::pow(acousticDoping, parameters.lambda) /
            (std::cbrt(field) * std::pow(relativeTemperature, parameters.k));
        const Real roughnessDoping = parameters.alpha_sr_A *
            std::pow(rawAcceptors + parameters.N2 / 2.0, parameters.lambda_sr_A) +
            parameters.alpha_sr_D * std::pow(rawDonors + parameters.N2 / 2.0, parameters.lambda_sr_D);
        // The roughness exponent depends on both carriers even with -FullPhuMob.
        const Real exponent = parameters.A + parameters.alpha_sr * (n + p) * concentrationFactor /
            std::pow(rawDonors + rawAcceptors + parameters.N1, parameters.nu);
        const Real inverseRoughness = (std::exp(exponent * std::log(field) - std::log(parameters.delta)) +
            field * field * field / parameters.eta) / std::pow(roughnessDoping, parameters.lambda_sr);
        result.roughness = inverse(inverseRoughness);
    }
    const Real inverse2D = result.interfaceDamping / result.phonon2D;
    if (model.phonon_combination == 0)
        result.phonon = inverse(std::max(inverse2D, inverse(result.phonon3D)));
    else {
        const Real bulkFactor = model.phonon_combination == 2 && hasInterface ? result.bulkWeight : 1.0;
        result.phonon = inverse(inverse2D + bulkFactor / result.phonon3D);
    }
    result.lowField = inverse(inverse(result.coulomb) + inverse(result.phonon) +
        result.interfaceDamping / result.roughness);
    if (!std::isfinite(result.lowField) || result.lowField < 0.0)
        throw std::runtime_error("IALMob produced a non-finite or negative mobility");
    return result;
}

Real IALMobility::evaluate(const Material& material, CarrierType carrier, MobilityDopants dopants,
    Real n, Real p, Real field, Real normal, Real distance) const
{
    Real mobility = components(material, carrier, dopants, n, p, normal, distance).lowField *
        1.0e-4 / config_.internalMobilityToM2PerVS;
    if (config_.model == "ialmob_field") {
        const auto& parameters = carrier == CarrierType::Electron ? config_.electronField : config_.holeField;
        if (!std::isfinite(field) || parameters.saturationVelocity <= 0.0 || parameters.beta <= 0.0)
            throw std::invalid_argument("IALMob: invalid high-field drive or saturation parameters");
        const Real ratio = mobility * std::abs(field) / parameters.saturationVelocity;
        mobility /= std::pow(1.0 + std::pow(ratio, parameters.beta), 1.0 / parameters.beta);
    }
    return mobility;
}

Real IALMobility::electronMobility(const Material& material, Real, Real n, Real p,
    Real field, Real normal, Real distance, MobilityDopants dopants) const
{
    return evaluate(material, CarrierType::Electron, dopants, n, p, field, normal, distance);
}
Real IALMobility::holeMobility(const Material& material, Real, Real n, Real p,
    Real field, Real normal, Real distance, MobilityDopants dopants) const
{
    return evaluate(material, CarrierType::Hole, dopants, n, p, field, normal, distance);
}
} // namespace vela
