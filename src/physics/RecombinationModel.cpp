#include "vela/physics/RecombinationModel.h"
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

} // namespace

RecombinationModel::RecombinationModel(RecombinationModelConfig config)
    : config_(std::move(config))
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
    if (augerEnabled_ && (config_.augerCn < 0.0 || config_.augerCp < 0.0))
        throw std::invalid_argument("RecombinationModel: Auger coefficients cannot be negative.");
}

Real RecombinationModel::srhDenominator(Real n, Real p, Real ni) const
{
    return config_.taup * (n + ni) + config_.taun * (p + ni);
}

Real RecombinationModel::srhRate(Real n, Real p, Real ni) const
{
    return srhRateFromExcessProduct(n * p - ni * ni, n, p, ni);
}

Real RecombinationModel::srhRateFromExcessProduct(Real excessProduct,
                                                  Real n,
                                                  Real p,
                                                  Real ni) const
{
    if (!srhEnabled_)
        return 0.0;
    const Real den = srhDenominator(n, p, ni);
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

Real RecombinationModel::totalRate(Real n, Real p, Real ni) const
{
    const Real excessProduct = limitedExcessProduct(n, p, ni);
    return totalRateFromExcessProduct(excessProduct, n, p, ni);
}

Real RecombinationModel::totalRateFromExcessProduct(Real excessProduct,
                                                    Real n,
                                                    Real p,
                                                    Real ni) const
{
    return srhRateFromExcessProduct(excessProduct, n, p, ni)
         + augerRateFromExcessProduct(excessProduct, n, p);
}

RecombinationRateDerivatives RecombinationModel::totalRateDerivativesFromExcessProduct(
    Real excessProduct,
    Real n,
    Real p,
    Real ni) const
{
    RecombinationRateDerivatives derivatives;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni);
        if (std::abs(den) >= 1.0e-100) {
            derivatives.dRateDExcess += 1.0 / den;
            const Real invDen2 = 1.0 / (den * den);
            const Real limitedExcess = limitedExcessValue(excessProduct);
            derivatives.dRateDn -= limitedExcess * config_.taup * invDen2;
            derivatives.dRateDp -= limitedExcess * config_.taun * invDen2;
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
    Real ni) const
{
    RecombinationLinearization linearization;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni);
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
    Real ni) const
{
    RecombinationLinearization linearization;

    if (srhEnabled_) {
        const Real den = srhDenominator(n, p, ni);
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
    Real taup)
{
    RecombinationModelConfig config;
    config.mechanisms = std::move(mechanisms);
    config.taun = taun;
    config.taup = taup;
    return config;
}

} // namespace vela
