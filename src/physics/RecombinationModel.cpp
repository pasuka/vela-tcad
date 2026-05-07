#include "vela/physics/RecombinationModel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

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
    return excessProduct / den;
}

Real RecombinationModel::augerRate(Real n, Real p, Real ni) const
{
    return augerRateFromExcessProduct(n * p - ni * ni, n, p);
}

Real RecombinationModel::augerRateFromExcessProduct(Real excessProduct,
                                                    Real n,
                                                    Real p) const
{
    if (!augerEnabled_)
        return 0.0;
    return (config_.augerCn * n + config_.augerCp * p) * excessProduct;
}

Real RecombinationModel::totalRate(Real n, Real p, Real ni) const
{
    const Real excessProduct = n * p - ni * ni;
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
        const Real excessProduct = n * p - ni * ni;
        const Real prefactor = config_.augerCn * n + config_.augerCp * p;
        const Real rate = prefactor * excessProduct;
        const Real derivative = config_.augerCn * excessProduct + prefactor * p;
        const Real positiveDerivative = std::max(derivative, 0.0);
        linearization.diagonal += positiveDerivative;
        linearization.rhs += positiveDerivative * n - rate;
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
        const Real excessProduct = n * p - ni * ni;
        const Real prefactor = config_.augerCn * n + config_.augerCp * p;
        const Real rate = prefactor * excessProduct;
        const Real derivative = config_.augerCp * excessProduct + prefactor * n;
        const Real positiveDerivative = std::max(derivative, 0.0);
        linearization.diagonal += positiveDerivative;
        linearization.rhs += positiveDerivative * p - rate;
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
