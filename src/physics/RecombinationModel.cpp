#include "vela/physics/RecombinationModel.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

namespace {

constexpr Real kMaxAugerCarrier = 1.0e30;

Real limitedAugerCarrier(Real value)
{
    return std::clamp(value, -kMaxAugerCarrier, kMaxAugerCarrier);
}

Real limitedExcessProduct(Real n, Real p, Real ni)
{
    const Real limitedN = limitedAugerCarrier(n);
    const Real limitedP = limitedAugerCarrier(p);
    return limitedN * limitedP - ni * ni;
}

Real limitedExcessValue(Real value)
{
    constexpr Real limit = kMaxAugerCarrier * kMaxAugerCarrier;
    return std::clamp(value, -limit, limit);
}

Real dopingDependentLifetime(Real concentration,
                             const SRHLifetimeParameters& parameters)
{
    const Real doping = std::abs(concentration);
    const Real ratio = doping / parameters.referenceDoping;
    return parameters.tauMin +
        (parameters.tauMax - parameters.tauMin) /
            (1.0 + std::pow(ratio, parameters.gamma));
}

void validateLifetimeParameters(const SRHLifetimeParameters& parameters,
                                const char* carrier)
{
    if (!std::isfinite(parameters.tauMin) || parameters.tauMin < 0.0 ||
        !std::isfinite(parameters.tauMax) || parameters.tauMax <= 0.0 ||
        parameters.tauMin > parameters.tauMax ||
        !std::isfinite(parameters.referenceDoping) ||
        parameters.referenceDoping <= 0.0 ||
        !std::isfinite(parameters.gamma) || parameters.gamma <= 0.0) {
        throw std::invalid_argument(
            std::string("RecombinationModel: invalid ") + carrier +
            " SRH doping-dependent lifetime parameters.");
    }
}

} // namespace

RecombinationModel::RecombinationModel(RecombinationModelConfig config)
    : config_(std::move(config))
    , bandToBand_(config_.bandToBand)
{
    for (const std::string& mechanism : config_.mechanisms) {
        if (mechanism == "srh") {
            srhEnabled_ = true;
        } else if (mechanism == "auger") {
            augerEnabled_ = true;
        } else if (mechanism == "none") {
            // Explicitly accepted for experiments that disable recombination.
        } else {
            throw std::invalid_argument(
                "RecombinationModel: unknown recombination mechanism '" + mechanism + "'.");
        }
    }

    if (srhEnabled_ && (config_.taun <= 0.0 || config_.taup <= 0.0))
        throw std::invalid_argument("RecombinationModel: SRH lifetimes must be positive.");
    if (config_.srhDopingDependence.enabled) {
        if (config_.srhDopingDependence.concentrationBasis != "total_impurity" &&
            config_.srhDopingDependence.concentrationBasis != "net_doping") {
            throw std::invalid_argument(
                "RecombinationModel: SRH doping concentration_basis must be "
                "'total_impurity' or 'net_doping'.");
        }
        validateLifetimeParameters(config_.srhDopingDependence.electron, "electron");
        validateLifetimeParameters(config_.srhDopingDependence.hole, "hole");
        if (config_.srhDopingDependence.temperatureDependence &&
            (!std::isfinite(config_.srhDopingDependence.temperature_K) ||
             config_.srhDopingDependence.temperature_K <= 0.0 ||
             !std::isfinite(config_.srhDopingDependence.referenceTemperature_K) ||
             config_.srhDopingDependence.referenceTemperature_K <= 0.0 ||
             !std::isfinite(config_.srhDopingDependence.electronTemperatureExponent) ||
             !std::isfinite(config_.srhDopingDependence.holeTemperatureExponent))) {
            throw std::invalid_argument(
                "RecombinationModel: invalid SRH temperature-dependence parameters.");
        }
    }
    if (augerEnabled_ && (config_.augerCn < 0.0 || config_.augerCp < 0.0))
        throw std::invalid_argument("RecombinationModel: Auger coefficients cannot be negative.");
}

Real RecombinationModel::electronLifetime(Real dopingConcentration) const
{
    const Real base = config_.srhDopingDependence.enabled
        ? dopingDependentLifetime(
            dopingConcentration, config_.srhDopingDependence.electron)
        : config_.taun;
    if (!config_.srhDopingDependence.temperatureDependence)
        return base;
    return base * std::pow(
        config_.srhDopingDependence.temperature_K /
            config_.srhDopingDependence.referenceTemperature_K,
        config_.srhDopingDependence.electronTemperatureExponent);
}

Real RecombinationModel::holeLifetime(Real dopingConcentration) const
{
    const Real base = config_.srhDopingDependence.enabled
        ? dopingDependentLifetime(
            dopingConcentration, config_.srhDopingDependence.hole)
        : config_.taup;
    if (!config_.srhDopingDependence.temperatureDependence)
        return base;
    return base * std::pow(
        config_.srhDopingDependence.temperature_K /
            config_.srhDopingDependence.referenceTemperature_K,
        config_.srhDopingDependence.holeTemperatureExponent);
}

Real RecombinationModel::srhDopingConcentration(
    Real donors, Real acceptors) const
{
    if (config_.srhDopingDependence.concentrationBasis == "net_doping")
        return std::abs(donors - acceptors);
    return std::abs(donors) + std::abs(acceptors);
}

Real RecombinationModel::srhDenominator(
    Real n, Real p, Real ni, Real dopingConcentration) const
{
    return holeLifetime(dopingConcentration) * (n + ni) +
           electronLifetime(dopingConcentration) * (p + ni);
}

Real RecombinationModel::srhRate(
    Real n, Real p, Real ni, Real dopingConcentration) const
{
    return srhRateFromExcessProduct(
        n * p - ni * ni, n, p, ni, dopingConcentration);
}

Real RecombinationModel::srhRateFromExcessProduct(Real excessProduct,
                                                  Real n,
                                                  Real p,
                                                  Real ni,
                                                  Real dopingConcentration) const
{
    if (!srhEnabled_)
        return 0.0;
    const Real den = srhDenominator(n, p, ni, dopingConcentration);
    if (std::abs(den) < 1.0e-100)
        return 0.0;
    return limitedExcessValue(excessProduct) / den;
}

Real RecombinationModel::augerRate(Real n, Real p, Real ni) const
{
    return augerRateFromExcessProduct(limitedExcessProduct(n, p, ni), n, p);
}

Real RecombinationModel::augerRateFromExcessProduct(Real excessProduct,
                                                    Real n,
                                                    Real p) const
{
    if (!augerEnabled_)
        return 0.0;
    const Real limitedN = limitedAugerCarrier(n);
    const Real limitedP = limitedAugerCarrier(p);
    return (config_.augerCn * limitedN + config_.augerCp * limitedP) *
           limitedExcessValue(excessProduct);
}

Real RecombinationModel::totalRate(
    Real n, Real p, Real ni, Real dopingConcentration) const
{
    const Real excessProduct = limitedExcessProduct(n, p, ni);
    return totalRateFromExcessProduct(
        excessProduct, n, p, ni, dopingConcentration);
}

Real RecombinationModel::totalRateFromExcessProduct(Real excessProduct,
                                                    Real n,
                                                    Real p,
                                                    Real ni,
                                                    Real dopingConcentration) const
{
    return srhRateFromExcessProduct(
               excessProduct, n, p, ni, dopingConcentration)
         + augerRateFromExcessProduct(excessProduct, n, p);
}

RecombinationRateDerivatives RecombinationModel::totalRateDerivativesFromExcessProduct(
    Real excessProduct,
    Real n,
    Real p,
    Real ni,
    Real dopingConcentration) const
{
    RecombinationRateDerivatives derivatives;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni, dopingConcentration);
        if (std::abs(den) >= 1.0e-100) {
            derivatives.dRateDExcess += 1.0 / den;
            const Real invDen2 = 1.0 / (den * den);
            const Real limitedExcess = limitedExcessValue(excessProduct);
            derivatives.dRateDn -= limitedExcess *
                holeLifetime(dopingConcentration) * invDen2;
            derivatives.dRateDp -= limitedExcess *
                electronLifetime(dopingConcentration) * invDen2;
        }
    }

    if (augerEnabled_) {
        const Real limitedN = limitedAugerCarrier(n);
        const Real limitedP = limitedAugerCarrier(p);
        excessProduct = limitedExcessValue(excessProduct);
        const Real prefactor = config_.augerCn * limitedN + config_.augerCp * limitedP;
        derivatives.dRateDExcess += prefactor;
        derivatives.dRateDn += config_.augerCn * excessProduct;
        derivatives.dRateDp += config_.augerCp * excessProduct;
    }

    return derivatives;
}

RecombinationLinearization RecombinationModel::electronLinearization(
    Real n,
    Real p,
    Real ni,
    Real dopingConcentration) const
{
    RecombinationLinearization linearization;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni, dopingConcentration);
        if (den > 1.0e-100) {
            linearization.diagonal += p / den;
            linearization.rhs += ni * ni / den;
        }
    }

    if (augerEnabled_) {
        const Real limitedN = limitedAugerCarrier(n);
        const Real limitedP = limitedAugerCarrier(p);
        const Real excessProduct = limitedExcessProduct(n, p, ni);
        const Real prefactor = config_.augerCn * limitedN + config_.augerCp * limitedP;
        const Real rate = prefactor * excessProduct;
        const Real derivative = config_.augerCn * excessProduct + prefactor * limitedP;
        const Real positiveDerivative = std::max(derivative, 0.0);
        linearization.diagonal += positiveDerivative;
        linearization.rhs += positiveDerivative * limitedN - rate;
    }

    return linearization;
}

RecombinationLinearization RecombinationModel::holeLinearization(
    Real n,
    Real p,
    Real ni,
    Real dopingConcentration) const
{
    RecombinationLinearization linearization;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni, dopingConcentration);
        if (den > 1.0e-100) {
            linearization.diagonal += n / den;
            linearization.rhs += ni * ni / den;
        }
    }

    if (augerEnabled_) {
        const Real limitedN = limitedAugerCarrier(n);
        const Real limitedP = limitedAugerCarrier(p);
        const Real excessProduct = limitedExcessProduct(n, p, ni);
        const Real prefactor = config_.augerCn * limitedN + config_.augerCp * limitedP;
        const Real rate = prefactor * excessProduct;
        const Real derivative = config_.augerCp * excessProduct + prefactor * limitedN;
        const Real positiveDerivative = std::max(derivative, 0.0);
        linearization.diagonal += positiveDerivative;
        linearization.rhs += positiveDerivative * limitedP - rate;
    }

    return linearization;
}

RecombinationModelConfig recombinationModelConfig(
    std::vector<std::string> mechanisms,
    Real taun,
    Real taup,
    SRHDopingDependenceConfig srhDopingDependence)
{
    RecombinationModelConfig config;
    config.mechanisms = std::move(mechanisms);
    config.taun = taun;
    config.taup = taup;
    config.srhDopingDependence = std::move(srhDopingDependence);
    return config;
}

SRHDopingDependenceConfig srhDopingDependenceConfigFromJson(
    const nlohmann::json& value,
    UnitScalingConfig scaling)
{
    if (!value.is_object())
        throw std::invalid_argument(
            "srh_doping_dependence must be an object.");

    SRHDopingDependenceConfig config;
    // The compiled reference concentrations are SI literals, so express them
    // in the active internal concentration unit before any deck override.
    config.electron.referenceDoping =
        scaling.unitSystem().m3ToInternalConcentration(config.electron.referenceDoping);
    config.hole.referenceDoping =
        scaling.unitSystem().m3ToInternalConcentration(config.hole.referenceDoping);
    config.enabled = value.value("enabled", true);
    config.concentrationBasis = value.value(
        "concentration_basis", config.concentrationBasis);
    config.temperatureDependence = value.value(
        "temperature_dependence", config.temperatureDependence);
    config.temperature_K = value.value("temperature_K", config.temperature_K);
    config.referenceTemperature_K = value.value(
        "reference_temperature_K", config.referenceTemperature_K);
    config.electronTemperatureExponent = value.value(
        "electron_temperature_exponent", config.electronTemperatureExponent);
    config.holeTemperatureExponent = value.value(
        "hole_temperature_exponent", config.holeTemperatureExponent);

    auto parseCarrier = [&](const char* name, SRHLifetimeParameters& parameters) {
        if (!value.contains(name))
            return;
        const auto& carrier = value.at(name);
        if (!carrier.is_object())
            throw std::invalid_argument(
                std::string("srh_doping_dependence.") + name +
                " must be an object.");
        parameters.tauMin = carrier.value("tau_min_s", parameters.tauMin);
        parameters.tauMax = carrier.value("tau_max_s", parameters.tauMax);
        if (carrier.contains("reference_doping_m3")) {
            parameters.referenceDoping = scaling.concentrationToInternal(
                carrier.at("reference_doping_m3").get<Real>());
        }
        parameters.gamma = carrier.value("gamma", parameters.gamma);
    };
    parseCarrier("electron", config.electron);
    parseCarrier("hole", config.hole);

    RecombinationModelConfig validation;
    validation.mechanisms = {"srh"};
    validation.srhDopingDependence = config;
    (void)RecombinationModel(validation);
    return config;
}

} // namespace vela
