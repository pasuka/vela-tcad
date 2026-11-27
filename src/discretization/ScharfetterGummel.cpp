#include "vela/core/ArithmeticPrecision.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/SgFluxArithmetic.h"
#include "vela/physics/CarrierStatistics.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>


namespace vela {

namespace {

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
}

constexpr Real ProductionExponentClamp = 500.0;

Real logBernoulli(Real value)
{
    if (std::abs(value) < 1.0e-10) {
        const Real series = 1.0 - value * 0.5 + value * value / 12.0;
        return std::log(series);
    }
    if (value > 500.0)
        return std::log(value) - value;
    if (value < -500.0)
        return std::log(-value);
    return std::log(value / std::expm1(value));
}

Real signedSaturatedExp(Real logMagnitude, Real sign)
{
    // Invalid arithmetic is not a finite overflow and must remain observable.
    if (std::isnan(logMagnitude) || !std::isfinite(sign))
        return std::numeric_limits<Real>::quiet_NaN();
    if (sign == 0.0)
        return 0.0;
    const Real logMaximum = std::log(std::numeric_limits<Real>::max());
    if (logMagnitude >= logMaximum)
        return std::copysign(std::numeric_limits<Real>::max(), sign);
    const Real magnitude = std::exp(logMagnitude);
    if (!std::isfinite(magnitude))
        return std::copysign(std::numeric_limits<Real>::max(), sign);
    return std::copysign(magnitude, sign);
}

Real stableBernoulliDensityDifferenceFlux(Real ni1,
                                          Real exponent1,
                                          Real eta,
                                          Real coef,
                                          Real logLeftOverRight)
{
    if (coef == 0.0 || logLeftOverRight == 0.0)
        return 0.0;

    const Real logAbsRelativeDifference = logLeftOverRight > 50.0
        ? logLeftOverRight + std::log1p(-std::exp(-logLeftOverRight))
        : std::log(std::abs(std::expm1(logLeftOverRight)));
    const Real logRight = logBernoulli(eta) + std::log(ni1) + exponent1;
    const Real logMagnitude =
        std::log(std::abs(coef)) + logRight + logAbsRelativeDifference;
    const Real sign = (coef > 0.0 ? 1.0 : -1.0)
        * (logLeftOverRight > 0.0 ? 1.0 : -1.0);
    return signedSaturatedExp(logMagnitude, sign);
}

Real finiteCancellationCondition(Real leftTerm,
                                 Real rightTerm,
                                 Real signedDifference)
{
    const Real termScale = std::abs(leftTerm) + std::abs(rightTerm);
    if (termScale == 0.0)
        return 0.0;
    if (signedDifference == 0.0)
        return std::numeric_limits<Real>::max();

    const Real condition = termScale / std::abs(signedDifference);
    return std::isfinite(condition)
        ? condition
        : std::numeric_limits<Real>::max();
}
vela::ArithmeticReal longDoubleBernoulli(vela::ArithmeticReal value)
{
    if (std::abs(value) < VELA_ARITHMETIC_LITERAL(1.0e-12))
        return VELA_ARITHMETIC_LITERAL(1.0) - value * VELA_ARITHMETIC_LITERAL(0.5) + value * value / VELA_ARITHMETIC_LITERAL(12.0);
    if (value > VELA_ARITHMETIC_LITERAL(50.0)) {
        const vela::ArithmeticReal expNegative = std::exp(-value);
        return value * expNegative / (VELA_ARITHMETIC_LITERAL(1.0) - expNegative);
    }
    if (value < -VELA_ARITHMETIC_LITERAL(50.0)) {
        const vela::ArithmeticReal expValue = std::exp(value);
        return -value / (VELA_ARITHMETIC_LITERAL(1.0) - expValue);
    }
    return value / std::expm1(value);
}

ArithmeticReal arithmeticBernoulliDerivative(ArithmeticReal x)
{
    if (std::abs(x) < VELA_ARITHMETIC_LITERAL(0.01)) {
        const ArithmeticReal x2 = x * x;
        return -VELA_ARITHMETIC_LITERAL(0.5) + x *
            (ArithmeticReal(1) / 6 - x2 * (ArithmeticReal(1) / 180 -
                x2 * (ArithmeticReal(1) / 5040 - x2 / 151200)));
    }
    if (x > 0) {
        const ArithmeticReal e = std::exp(-x);
        const ArithmeticReal d = -std::expm1(-x);
        return e * (d - x) / (d * d);
    }
    const ArithmeticReal e = std::exp(x);
    const ArithmeticReal d = std::expm1(x);
    return (d - x * e) / (d * d);
}

Real finiteRealFromLongDouble(vela::ArithmeticReal value)
{
    const vela::ArithmeticReal maximum =
        static_cast<vela::ArithmeticReal>(std::numeric_limits<Real>::max());
    if (std::isnan(value))
        return 0.0;
    if (value >= maximum)
        return std::numeric_limits<Real>::max();
    if (value <= -maximum)
        return -std::numeric_limits<Real>::max();
    return static_cast<Real>(value);
}

double positiveFrexp(double value, int& exponent)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t));
    static_assert(std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits == 53);
    constexpr std::uint64_t fractionMask = (std::uint64_t{1} << 52) - 1;
    constexpr std::uint64_t halfExponent = std::uint64_t{1022} << 52;
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto storedExponent = static_cast<int>(bits >> 52);
    if (storedExponent == 0)
        return std::frexp(value, &exponent);
    exponent = storedExponent - 1022;
    return std::bit_cast<double>((bits & fractionMask) | halfExponent);
}



} // namespace

SGEdgeWeights sgEdgeWeights(Real dpsi, Real Vt)
{
    const Real u = dpsi / Vt;
    return SGEdgeWeights{bernoulli(u), bernoulli(-u)};
}

Real sgElectronContinuityFlux(Real n0, Real n1, Real dpsi, Real Vt, Real coef)
{
    const SGEdgeWeights w = sgEdgeWeights(dpsi, Vt);
    return detail::sgWeightedPopulationDifference(coef, n0, n1, w.b_minus, w.b_plus);
}

Real sgHoleContinuityFlux(Real p0, Real p1, Real dpsi, Real Vt, Real coef)
{
    const SGEdgeWeights w = sgEdgeWeights(dpsi, Vt);
    return detail::sgWeightedPopulationDifference(coef, p0, p1, w.b_plus, w.b_minus);
}

ArithmeticReal arithmeticPopulationLogRatio(ArithmeticReal numerator, ArithmeticReal denominator)
{
    const vela::ArithmeticReal low = std::min(numerator, denominator);
    const vela::ArithmeticReal high = std::max(numerator, denominator);
    const vela::ArithmeticReal relativeDifference = (high - low) / low;
    const vela::ArithmeticReal magnitude = std::isfinite(relativeDifference)
        ? std::log1p(relativeDifference) : std::log(high) - std::log(low);
    return numerator < denominator ? -magnitude : magnitude;
}

ArithmeticReal sgPopulationLogRatio(Real numerator, Real denominator)
{
    return arithmeticPopulationLogRatio(numerator, denominator);
}

numerics::SplitCoordinate sgPopulationLogRatioCompensated(Real numerator, Real denominator)
{
    using numerics::SplitCoordinate;
    if (!(numerator > 0. && denominator > 0. && std::isfinite(numerator) && std::isfinite(denominator)))
        throw std::domain_error("Expected finite positive SG populations");
    if (numerator == denominator) return {};
    // The positive ratio reduces to m * 2^e with m close to one. Evaluate
    // log(m) = 2*atanh((m-1)/(m+1)) without rounding the ratio first.
    int aExponent, bExponent;
    double a = positiveFrexp(numerator, aExponent);
    const double b = positiveFrexp(denominator, bExponent);
    int exponent = aExponent - bExponent;
    constexpr double sqrtTwo = 0x1.6a09e667f3bcdp+0;
    if (a > b * sqrtTwo) {
        a *= .5;
        ++exponent;
    } else if (b > a * sqrtTwo) {
        a *= 2.;
        --exponent;
    }
    const auto sum = SplitCoordinate::sum(a, b);
    const double difference = a - b;  // Exact for the reduced ratio.
    const double high = difference / sum.hi;
    const auto t = SplitCoordinate::sum(high, (std::fma(-high, sum.hi, difference) - high * sum.lo) / sum.hi);
    const auto t2 = t.multiplied(t);
    // |t| <= 0.172. Rounding the terms after t^5 contributes less than
    // 1e-21. The first omitted term after t^25 is smaller still.
    constexpr SplitCoordinate oneFifth{0x1.999999999999ap-3, -0x1.999999999999ap-57};
    constexpr SplitCoordinate oneThird{0x1.5555555555555p-2, 0x1.5555555555555p-56};
    constexpr SplitCoordinate logTwo{0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56};
    auto lowTerms = oneFifth;
    if (std::abs(t.hi) > 1e-3) {
        double polynomial = 1. / 25.;
        for (int k = 11; k >= 3; --k)
            polynomial = 1. / (2 * k + 1) + t2.hi * polynomial;
        lowTerms = lowTerms.added(t2.scaled(polynomial));
    }
    const auto series = oneThird.added(t2.multiplied(lowTerms));
    const auto logMantissa = t.added(t.multiplied(t2).multiplied(series)).scaled(2.);
    return exponent == 0 ? logMantissa : logMantissa.added(logTwo.scaled(exponent));
}

namespace {
ArithmeticReal populationFluxFromLogImbalance(
    ArithmeticReal population0, ArithmeticReal population1, ArithmeticReal eta,
    ArithmeticReal coef, ArithmeticReal correction0, ArithmeticReal correction1, ArithmeticReal logRatio)
{
    if (logRatio >= 0.)
        return coef * longDoubleBernoulli(-eta) * population0
            * (ArithmeticReal(1.) + correction0) * (-std::expm1(-logRatio));
    return coef * longDoubleBernoulli(eta) * population1
        * (ArithmeticReal(1.) + correction1) * std::expm1(logRatio);
}
} // namespace

SGPopulationEdgeDifferential sgBoltzmannPopulationEdge(
    ArithmeticReal psi0, ArithmeticReal psi1, ArithmeticReal population0,
    ArithmeticReal population1, ArithmeticReal ni0, ArithmeticReal ni1,
    ArithmeticReal signedPotentialScale, ArithmeticReal Vt, ArithmeticReal coefficient,
    bool includeNiGradientDrift)
{
    if (!(Vt > 0) || !std::isfinite(Vt))
        throw std::domain_error("SG population edge requires finite positive thermal voltage");
    // Keep the coordinate difference before any conversion to volts.
    auto eta = (psi1 - psi0) * signedPotentialScale / Vt;
    if (includeNiGradientDrift && ni0 > 0 && ni1 > 0)
        eta += arithmeticPopulationLogRatio(ni1, ni0);
    const auto left = coefficient * longDoubleBernoulli(-eta);
    const auto right = -coefficient * longDoubleBernoulli(eta);
    const auto flux = population0 > 0 && population1 > 0
        ? populationFluxFromLogImbalance(population0, population1, eta, coefficient, 0, 0,
                                         eta + arithmeticPopulationLogRatio(population0, population1))
        : left * population0 + right * population1;
    const auto potential1 = -coefficient * signedPotentialScale / Vt *
        (arithmeticBernoulliDerivative(-eta) * population0 +
         arithmeticBernoulliDerivative(eta) * population1);
    return {flux, left, right, potential1};
}

SGPopulationEdgeDifferential sgFermiPopulationEdge(
    ArithmeticReal psi0, ArithmeticReal psi1, ArithmeticReal population0,
    ArithmeticReal population1, ArithmeticReal intrinsic0, ArithmeticReal intrinsic1,
    ArithmeticReal dos0, ArithmeticReal dos1, ArithmeticReal signedPotentialScale,
    ArithmeticReal Vt, ArithmeticReal coefficient, CarrierStatisticsModel model)
{
    for (const auto value : {population0, population1, intrinsic0, intrinsic1, dos0, dos1, Vt})
        if (!(value > 0) || !std::isfinite(value))
            throw std::domain_error("Effective-band SG requires positive finite populations, material and scale");
    const auto ratio0 = static_cast<long double>(population0) / dos0;
    const auto ratio1 = static_cast<long double>(population1) / dos1;
    const auto correction0 = extendedFermiDegeneracy(ratio0, model);
    const auto correction1 = extendedFermiDegeneracy(ratio1, model);
    const auto delta = static_cast<ArithmeticReal>(correction0.logCorrection - correction1.logCorrection);
    const auto potentialFactor = signedPotentialScale / Vt;
    const auto classical = (psi1 - psi0) * potentialFactor +
        arithmeticPopulationLogRatio(intrinsic1, intrinsic0);
    const auto drift = classical + delta;
    const auto imbalance = (classical + arithmeticPopulationLogRatio(population0, population1)) + delta;
    const auto flux = populationFluxFromLogImbalance(population0, population1, drift, coefficient, 0, 0, imbalance);
    const auto driftDerivative = -coefficient * (arithmeticBernoulliDerivative(-drift) * population0 +
                                                 arithmeticBernoulliDerivative(drift) * population1);
    // At fixed electrostatic potential the effective band itself depends on
    // concentration. A ratio that underflows still has the dilute derivative.
    const auto slope = [](long double ratio, long double logSlope, ArithmeticReal density, ArithmeticReal dos) {
        return static_cast<ArithmeticReal>(ratio >= std::numeric_limits<long double>::min()
            ? logSlope / density : 1 / (2 * std::sqrt(2.L) * dos));
    };
    return {flux,
            coefficient * longDoubleBernoulli(-drift) + driftDerivative * slope(ratio0, correction0.logSlope, population0, dos0),
            -coefficient * longDoubleBernoulli(drift) - driftDerivative * slope(ratio1, correction1.logSlope, population1, dos1),
            driftDerivative * potentialFactor};
}

numerics::SplitCoordinate sgPopulationEdgeDrift(
    Real psi0, Real psi1, Real signedPotentialScale, Real Vt,
    Real ni0, Real ni1, bool includeNiGradientDrift)
{
    auto eta = numerics::SplitCoordinate::sum(psi1, -psi0)
        .scaled(signedPotentialScale).divided(Vt);
    if (includeNiGradientDrift && ni0 > 0. && ni1 > 0.)
        eta = eta.added(sgPopulationLogRatioCompensated(ni1, ni0));
    return eta;
}

Real sgBoltzmannPopulationContinuityFlux(
    Real population0, Real population1, Real ni0, Real ni1,
    vela::ArithmeticReal signedDpsi, Real Vt, Real coef, bool includeNiGradientDrift,
    Real relativeCorrection0, Real relativeCorrection1, vela::ArithmeticReal signedDpsiCorrection)
{
    if constexpr (compiledArithmeticPrecision == ArithmeticPrecision::Double) {
        if (std::isfinite(signedDpsi) && Vt > 0. && std::isfinite(signedDpsi / Vt)) {
            auto eta = numerics::SplitCoordinate{static_cast<Real>(signedDpsi), 0.}.divided(Vt);
            if (includeNiGradientDrift && ni0 > 0. && ni1 > 0.)
                eta = eta.added(sgPopulationLogRatioCompensated(ni1, ni0));
            return sgBoltzmannPopulationContinuityFluxFromDrift(
                population0, population1, eta, Vt, coef, relativeCorrection0, relativeCorrection1,
                static_cast<Real>(signedDpsiCorrection));
        }
    }
    vela::ArithmeticReal eta = static_cast<vela::ArithmeticReal>(signedDpsi) / Vt;
    if (includeNiGradientDrift && ni0 > 0.0 && ni1 > 0.0)
        eta += sgPopulationLogRatio(ni1, ni0);
    return sgBoltzmannPopulationContinuityFluxFromDrift(
        population0, population1, eta, Vt, coef,
        relativeCorrection0, relativeCorrection1, signedDpsiCorrection);
}

Real sgBoltzmannPopulationContinuityFluxFromDrift(
    Real population0, Real population1, vela::ArithmeticReal eta, Real Vt, Real coef,
    Real relativeCorrection0, Real relativeCorrection1, vela::ArithmeticReal signedDpsiCorrection)
{
    if (population0 <= 0.0 || population1 <= 0.0) {
        return static_cast<Real>(coef * (longDoubleBernoulli(-eta) * population0
            - longDoubleBernoulli(eta) * population1));
    }

    // B(-eta) = exp(eta) * B(eta). Factor out the larger term and
    // evaluate the small relative imbalance without subtracting two fluxes.
    const vela::ArithmeticReal deltaEta = signedDpsiCorrection / Vt;
    const vela::ArithmeticReal logLeftOverRight = eta + sgPopulationLogRatio(population0, population1)
        + (deltaEta + std::log1p(static_cast<vela::ArithmeticReal>(relativeCorrection0))
           - std::log1p(static_cast<vela::ArithmeticReal>(relativeCorrection1)));
    eta += deltaEta;
    return populationFluxFromLogImbalance(
        population0, population1, eta, coef, relativeCorrection0, relativeCorrection1, logLeftOverRight);
}

Real sgBoltzmannPopulationContinuityFluxFromDrift(
    Real population0, Real population1, numerics::SplitCoordinate eta, Real Vt, Real coef,
    Real relativeCorrection0, Real relativeCorrection1, Real signedDpsiCorrection)
{
    if (population0 <= 0. || population1 <= 0.)
        return sgBoltzmannPopulationContinuityFluxFromDrift(
            population0, population1, ArithmeticReal(eta.hi), Vt, coef,
            relativeCorrection0, relativeCorrection1, signedDpsiCorrection);
    const auto correction = numerics::SplitCoordinate{signedDpsiCorrection, 0.}.divided(Vt);
    const double logCorrection0 = std::log1p(relativeCorrection0);
    const double logCorrection1 = std::log1p(relativeCorrection1);
    const auto imbalance = eta.added(sgPopulationLogRatioCompensated(population0, population1))
        .added(correction).added(numerics::SplitCoordinate::sum(logCorrection0, -logCorrection1));
    return populationFluxFromLogImbalance(population0, population1, eta.added(correction).hi,
                                          coef, relativeCorrection0, relativeCorrection1, imbalance.hi);
}

Real sgElectronContinuityFluxFromQuasiFermi(Real ni0,
                                            Real psi1,
                                            Real phin0,
                                            Real phin1,
                                            Real dpsi,
                                            Real Vt,
                                            Real coef)
{
    if (phin0 == phin1)
        return 0.0;

    const Real psi0 = psi1 - dpsi;
    return sgElectronContinuityFlux(
        ni0 * limitedExp((psi0 - phin0) / Vt),
        ni0 * limitedExp((psi1 - phin1) / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgElectronContinuityFluxFromQuasiFermiStable(Real ni0,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phin0,
                                                  Real phin1,
                                                  Real Vt,
                                                  Real coef)
{
    if (phin0 == phin1)
        return 0.0;

    // Equal-ni is the homogeneous limit of the factorized SG law. Retain
    // expm1(delta phin/Vt): subtracting exp(-phin0/Vt) and exp(-phin1/Vt)
    // loses the tiny majority-carrier increments in low-current PN sweeps.
    return sgElectronContinuityFluxFromQuasiFermiVariableNi(
        ni0, ni0, psi0, psi1, phin0, phin1, Vt, coef);
}

Real sgElectronContinuityFluxFromQuasiFermiFactors(Real ni0,
                                                   Real expPsi1,
                                                   Real expNegPhin0,
                                                   Real expNegPhin1,
                                                   Real dpsi,
                                                   Real Vt,
                                                   Real coef)
{
    const Real Bu = bernoulli(dpsi / Vt);
    return coef * Bu * ni0 * expPsi1 * (expNegPhin0 - expNegPhin1);
}

Real sgElectronContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                      Real ni1,
                                                      Real psi0,
                                                      Real psi1,
                                                      Real phin0,
                                                      Real phin1,
                                                      Real Vt,
                                                      Real coef,
                                                      bool includeNiGradientDrift)
{
    if (phin0 == phin1)
        return 0.0;
    if (ni0 <= 0.0 || ni1 <= 0.0)
        return sgElectronContinuityFlux(
            ni0 * limitedExp((psi0 - phin0) / Vt),
            ni1 * limitedExp((psi1 - phin1) / Vt),
            psi1 - psi0,
            Vt,
            coef);

    const Real eta = (psi1 - psi0) / Vt
        + (includeNiGradientDrift ? std::log(ni1 / ni0) : 0.0);
    const Real endpointExponent0 = (psi0 - phin0) / Vt;
    const Real endpointExponent1 = (psi1 - phin1) / Vt;
    const Real clampedExponent0 = std::clamp(
        endpointExponent0, -ProductionExponentClamp, ProductionExponentClamp);
    const Real clampedExponent1 = std::clamp(
        endpointExponent1, -ProductionExponentClamp, ProductionExponentClamp);
    Real logLeftOverRight = (phin1 - phin0) / Vt;
    if (!includeNiGradientDrift)
        logLeftOverRight += std::log(ni0 / ni1);
    logLeftOverRight +=
        (clampedExponent0 - endpointExponent0)
        - (clampedExponent1 - endpointExponent1);
    return stableBernoulliDensityDifferenceFlux(
        ni1, clampedExponent1, eta, coef, logLeftOverRight);
}

Real sgElectronBoltzmannContinuityFlux(Real ni0,
                                       Real ni1,
                                       Real psi0,
                                       Real psi1,
                                       Real phin0,
                                       Real phin1,
                                       Real Vt,
                                       Real coef,
                                       SGBoltzmannFluxPolicy policy)
{
    if (policy.includeNiGradientDrift ||
        (policy.compensatedEqualNi && ni0 == ni1)) {
        return sgElectronContinuityFluxFromQuasiFermiVariableNi(
            ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef,
            policy.includeNiGradientDrift);
    }
    if (ni0 == ni1) {
        return sgElectronContinuityFluxFromQuasiFermiStable(
            ni0, psi0, psi1, phin0, phin1, Vt, coef);
    }
    return sgElectronContinuityFlux(
        ni0 * limitedExp((psi0 - phin0) / Vt),
        ni1 * limitedExp((psi1 - phin1) / Vt),
        psi1 - psi0, Vt, coef);
}

SgElectronVariableNiFluxDecomposition
sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
    Real ni0,
    Real ni1,
    Real psi0,
    Real psi1,
    Real phin0,
    Real phin1,
    Real Vt,
    Real coef,
    bool includeNiGradientDrift)
{
    SgElectronVariableNiFluxDecomposition result;
    result.ni0 = ni0;
    result.ni1 = ni1;
    result.psi0 = psi0;
    result.psi1 = psi1;
    result.phin0 = phin0;
    result.phin1 = phin1;
    result.coef = coef;
    result.includeNiGradientDrift = includeNiGradientDrift;
    result.flatQuasiFermiShortCircuit = phin0 == phin1;

    const Real endpointExponent0 = (psi0 - phin0) / Vt;
    const Real endpointExponent1 = (psi1 - phin1) / Vt;
    result.node0ExponentClampedLow =
        endpointExponent0 < -ProductionExponentClamp;
    result.node0ExponentClampedHigh =
        endpointExponent0 > ProductionExponentClamp;
    result.node1ExponentClampedLow =
        endpointExponent1 < -ProductionExponentClamp;
    result.node1ExponentClampedHigh =
        endpointExponent1 > ProductionExponentClamp;
    const Real clampedExponent0 = std::clamp(
        endpointExponent0, -ProductionExponentClamp, ProductionExponentClamp);
    const Real clampedExponent1 = std::clamp(
        endpointExponent1, -ProductionExponentClamp, ProductionExponentClamp);

    result.n0 = ni0 * std::exp(clampedExponent0);
    result.n1 = ni1 * std::exp(clampedExponent1);
    result.eta = (psi1 - psi0) / Vt;
    if (ni0 > 0.0 && ni1 > 0.0 && includeNiGradientDrift)
        result.eta += std::log(ni1 / ni0);
    result.bernoulliMinusEta = bernoulli(-result.eta);
    result.bernoulliEta = bernoulli(result.eta);
    result.leftTerm = result.bernoulliMinusEta * result.n0;
    result.rightTerm = result.bernoulliEta * result.n1;
    result.signedDifference = result.flatQuasiFermiShortCircuit
        ? 0.0
        : result.leftTerm - result.rightTerm;
    result.reconstructedFlux = result.coef * result.signedDifference;
    result.cancellationCondition = finiteCancellationCondition(
        result.leftTerm, result.rightTerm, result.signedDifference);

    const vela::ArithmeticReal ni0Long = static_cast<vela::ArithmeticReal>(ni0);
    const vela::ArithmeticReal ni1Long = static_cast<vela::ArithmeticReal>(ni1);
    const vela::ArithmeticReal psi0Long = static_cast<vela::ArithmeticReal>(psi0);
    const vela::ArithmeticReal psi1Long = static_cast<vela::ArithmeticReal>(psi1);
    const vela::ArithmeticReal phin0Long = static_cast<vela::ArithmeticReal>(phin0);
    const vela::ArithmeticReal phin1Long = static_cast<vela::ArithmeticReal>(phin1);
    const vela::ArithmeticReal VtLong = static_cast<vela::ArithmeticReal>(Vt);
    const vela::ArithmeticReal coefLong = static_cast<vela::ArithmeticReal>(coef);
    const vela::ArithmeticReal endpointExponent0Long =
        (psi0Long - phin0Long) / VtLong;
    const vela::ArithmeticReal endpointExponent1Long =
        (psi1Long - phin1Long) / VtLong;
    const vela::ArithmeticReal clampedExponent0Long = std::clamp(
        endpointExponent0Long, -VELA_ARITHMETIC_LITERAL(500.0), VELA_ARITHMETIC_LITERAL(500.0));
    const vela::ArithmeticReal clampedExponent1Long = std::clamp(
        endpointExponent1Long, -VELA_ARITHMETIC_LITERAL(500.0), VELA_ARITHMETIC_LITERAL(500.0));
    const vela::ArithmeticReal n0Long = ni0Long * std::exp(clampedExponent0Long);
    const vela::ArithmeticReal n1Long = ni1Long * std::exp(clampedExponent1Long);
    vela::ArithmeticReal etaLong = (psi1Long - psi0Long) / VtLong;
    if (ni0 > 0.0 && ni1 > 0.0 && includeNiGradientDrift)
        etaLong += std::log(ni1Long / ni0Long);
    const vela::ArithmeticReal leftTermLong =
        longDoubleBernoulli(-etaLong) * n0Long;
    const vela::ArithmeticReal rightTermLong =
        longDoubleBernoulli(etaLong) * n1Long;
    const vela::ArithmeticReal referenceTermScaleLong =
        std::abs(coefLong)
        * (std::abs(leftTermLong) + std::abs(rightTermLong));
    vela::ArithmeticReal referenceFluxLong = VELA_ARITHMETIC_LITERAL(0.0);
    if (!result.flatQuasiFermiShortCircuit) {
        if (ni0 > 0.0 && ni1 > 0.0) {
            vela::ArithmeticReal logLeftOverRightLong =
                (phin1Long - phin0Long) / VtLong;
            if (!includeNiGradientDrift)
                logLeftOverRightLong += std::log(ni0Long / ni1Long);
            logLeftOverRightLong +=
                (clampedExponent0Long - endpointExponent0Long)
                - (clampedExponent1Long - endpointExponent1Long);
            referenceFluxLong = coefLong * rightTermLong
                * std::expm1(logLeftOverRightLong);
        } else {
            referenceFluxLong = coefLong * (leftTermLong - rightTermLong);
        }
    }
    result.highPrecisionReferenceFlux =
        finiteRealFromLongDouble(referenceFluxLong);
    result.highPrecisionReferenceTermScale =
        finiteRealFromLongDouble(referenceTermScaleLong);


    if (result.flatQuasiFermiShortCircuit) {
        result.stableFactorizedFlux = 0.0;
    } else if (ni0 > 0.0 && ni1 > 0.0
               && std::isfinite(clampedExponent0)
               && std::isfinite(clampedExponent1)
               && std::isfinite(result.eta)
               && std::isfinite(coef)) {
        Real logLeftOverRight = (phin1 - phin0) / Vt;
        if (!includeNiGradientDrift)
            logLeftOverRight += std::log(ni0 / ni1);
        logLeftOverRight +=
            (clampedExponent0 - endpointExponent0)
            - (clampedExponent1 - endpointExponent1);
        result.logLeftOverRight = logLeftOverRight;
        result.stableFactorizedFlux = stableBernoulliDensityDifferenceFlux(
            ni1, clampedExponent1, result.eta, coef, logLeftOverRight);
    } else {
        result.stableFactorizedFlux = result.reconstructedFlux;
    }

    return result;
}


Real sgHoleContinuityFluxFromQuasiFermi(Real ni0,
                                        Real psi0,
                                        Real phip0,
                                        Real phip1,
                                        Real dpsi,
                                        Real Vt,
                                        Real coef)
{
    if (phip0 == phip1)
        return 0.0;

    const Real psi1 = psi0 + dpsi;
    return sgHoleContinuityFlux(
        ni0 * limitedExp((phip0 - psi0) / Vt),
        ni0 * limitedExp((phip1 - psi1) / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgHoleContinuityFluxFromQuasiFermiStable(Real ni0,
                                              Real psi0,
                                              Real psi1,
                                              Real phip0,
                                              Real phip1,
                                              Real Vt,
                                              Real coef)
{
    if (phip0 == phip1)
        return 0.0;

    return sgHoleContinuityFluxFromQuasiFermiVariableNi(
        ni0, ni0, psi0, psi1, phip0, phip1, Vt, coef);
}

Real sgHoleContinuityFluxFromQuasiFermiFactors(Real ni0,
                                               Real expNegPsi0,
                                               Real expPhip0,
                                               Real expPhip1,
                                               Real dpsi,
                                               Real Vt,
                                               Real coef)
{
    const Real Bu = bernoulli(dpsi / Vt);
    return coef * Bu * ni0 * expNegPsi0 * (expPhip0 - expPhip1);
}

Real sgHoleContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                  Real ni1,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phip0,
                                                  Real phip1,
                                                  Real Vt,
                                                  Real coef,
                                                  bool includeNiGradientDrift)
{
    if (phip0 == phip1)
        return 0.0;
    if (ni0 <= 0.0 || ni1 <= 0.0)
        return sgHoleContinuityFlux(
            ni0 * limitedExp((phip0 - psi0) / Vt),
            ni1 * limitedExp((phip1 - psi1) / Vt),
            psi1 - psi0,
            Vt,
            coef);

    const Real eta = (psi1 - psi0) / Vt
        + (includeNiGradientDrift ? std::log(ni0 / ni1) : 0.0);
    const Real endpointExponent0 = (phip0 - psi0) / Vt;
    const Real endpointExponent1 = (phip1 - psi1) / Vt;
    const Real clampedExponent0 = std::clamp(
        endpointExponent0, -ProductionExponentClamp, ProductionExponentClamp);
    const Real clampedExponent1 = std::clamp(
        endpointExponent1, -ProductionExponentClamp, ProductionExponentClamp);
    Real logLeftOverRight = (phip0 - phip1) / Vt;
    if (!includeNiGradientDrift)
        logLeftOverRight += std::log(ni0 / ni1);
    logLeftOverRight +=
        (clampedExponent0 - endpointExponent0)
        - (clampedExponent1 - endpointExponent1);
    return stableBernoulliDensityDifferenceFlux(
        ni1, clampedExponent1, -eta, coef, logLeftOverRight);
}

Real sgHoleBoltzmannContinuityFlux(Real ni0,
                                   Real ni1,
                                   Real psi0,
                                   Real psi1,
                                   Real phip0,
                                   Real phip1,
                                   Real Vt,
                                   Real coef,
                                   SGBoltzmannFluxPolicy policy)
{
    if (policy.includeNiGradientDrift ||
        (policy.compensatedEqualNi && ni0 == ni1)) {
        return sgHoleContinuityFluxFromQuasiFermiVariableNi(
            ni0, ni1, psi0, psi1, phip0, phip1, Vt, coef,
            policy.includeNiGradientDrift);
    }
    if (ni0 == ni1) {
        return sgHoleContinuityFluxFromQuasiFermiStable(
            ni0, psi0, psi1, phip0, phip1, Vt, coef);
    }
    return sgHoleContinuityFlux(
        ni0 * limitedExp((phip0 - psi0) / Vt),
        ni1 * limitedExp((phip1 - psi1) / Vt),
        psi1 - psi0, Vt, coef);
}

Real sgGeneralizedEinsteinFactor(Real density0,
                                 Real density1,
                                 Real eta0,
                                 Real eta1)
{
    if (!(density0 > 0.0) || !(density1 > 0.0))
        return 1.0;
    const Real logDensityDifference = std::log(density1 / density0);
    const Real etaDifference = eta1 - eta0;
    if (std::abs(logDensityDifference) > 1.0e-10) {
        const Real factor = etaDifference / logDensityDifference;
        if (factor > 0.0 && std::isfinite(factor))
            return factor;
    }

    // Secant limit is the generalized Einstein factor F_{1/2}/F'_{1/2}.
    const Real etaMidpoint = 0.5 * (eta0 + eta1);
    const Real derivative = fermiDiracHalfDerivative(etaMidpoint);
    return derivative > 0.0
        ? std::max<Real>(1.0, fermiDiracHalf(etaMidpoint) / derivative)
        : 1.0;
}

namespace {

Real stableGeneralizedFermiDifference(Real rightDensity,
                                      Real rightBernoulliArgument,
                                      Real logarithmicLeftOverRight,
                                      Real coefficient)
{
    if (!std::isfinite(rightDensity) || !std::isfinite(rightBernoulliArgument) ||
        !std::isfinite(logarithmicLeftOverRight) || !std::isfinite(coefficient))
        return std::numeric_limits<Real>::quiet_NaN();
    if (!(rightDensity > 0.0) || coefficient == 0.0 ||
        logarithmicLeftOverRight == 0.0)
        return 0.0;

    const Real logAbsRelativeDifference = logarithmicLeftOverRight > 50.0
        ? logarithmicLeftOverRight
            + std::log1p(-std::exp(-logarithmicLeftOverRight))
        : std::log(std::abs(std::expm1(logarithmicLeftOverRight)));
    const Real logMagnitude = std::log(std::abs(coefficient))
        + logBernoulli(rightBernoulliArgument)
        + std::log(rightDensity)
        + logAbsRelativeDifference;
    const Real sign = (coefficient > 0.0 ? 1.0 : -1.0)
        * (logarithmicLeftOverRight > 0.0 ? 1.0 : -1.0);
    return signedSaturatedExp(logMagnitude, sign);
}

// Fold the Fermi density correction into the Bernoulli argument while
// preserving the quasi-Fermi driving drop and exact zero at equilibrium.
Real effectivePotentialFermiFlux(Real left, Real right, Real reducedDrop, Real coef)
{
    return stableGeneralizedFermiDifference(
        right, sgEffectiveBandArgument(left, right, reducedDrop), reducedDrop, coef);
}

vela::ArithmeticReal effectivePotentialFermiFluxExtended(
    Real left, Real right, vela::ArithmeticReal reducedDrop, Real coef)
{
    if (!(left > 0.0) || !(right > 0.0))
        return std::numeric_limits<ArithmeticReal>::quiet_NaN();
    const auto argument = std::log(static_cast<vela::ArithmeticReal>(right))
        - std::log(static_cast<vela::ArithmeticReal>(left)) + reducedDrop;
    return detail::sgFactorizedPopulationDifference(
        static_cast<ArithmeticReal>(coef) * longDoubleBernoulli(argument) *
            static_cast<ArithmeticReal>(right),
        reducedDrop, [](ArithmeticReal x) { return std::expm1(x); });
}

} // namespace

Real sgEffectiveBandArgument(Real density0, Real density1, Real signedQfDropOverVt)
{
    if (!(density0 > 0.0) || !(density1 > 0.0))
        return std::numeric_limits<Real>::quiet_NaN();
    const Real ratio = density1 / density0;
    const Real logRatio = ratio > 0.5 && ratio < 2.0
        ? std::log1p((density1 - density0) / density0)
        : std::log(density1) - std::log(density0);
    return logRatio + signedQfDropOverVt;
}

Real sgElectronFermiDiracContinuityFlux(
    Real n0, Real n1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef, FermiSgFormulation formulation)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    if (formulation == FermiSgFormulation::EffectivePotential)
        return effectivePotentialFermiFlux(n0, n1, (quasiFermi1 - quasiFermi0) / Vt, coef);
    const Real factor = sgGeneralizedEinsteinFactor(n0, n1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight = (quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        n1, argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return detail::sgWeightedPopulationDifference(
        coef * factor, n0, n1, bernoulli(-argument), bernoulli(argument));
}

Real sgElectronFermiDiracQuantumContinuityFlux(
    Real qmDensity0, Real qmDensity1,
    Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef, FermiSgFormulation formulation)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    if (formulation == FermiSgFormulation::EffectivePotential)
        return effectivePotentialFermiFlux(qmDensity0, qmDensity1, (quasiFermi1 - quasiFermi0) / Vt, coef);
    const Real factor = sgGeneralizedEinsteinFactor(
        fermiDensity0, fermiDensity1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight =
        (quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        qmDensity1, argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return detail::sgWeightedPopulationDifference(
        coef * factor, qmDensity0, qmDensity1, bernoulli(-argument), bernoulli(argument));
}

Real sgHoleFermiDiracContinuityFlux(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef, FermiSgFormulation formulation)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    if (formulation == FermiSgFormulation::EffectivePotential)
        return effectivePotentialFermiFlux(p0, p1, (quasiFermi0 - quasiFermi1) / Vt, coef);
    const Real factor = sgGeneralizedEinsteinFactor(p0, p1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight = -(quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        p1, -argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return detail::sgWeightedPopulationDifference(
        coef * factor, p0, p1, bernoulli(argument), bernoulli(-argument));
}

vela::ArithmeticReal sgElectronFermiDiracQuantumContinuityFluxLongDouble(
    Real qmDensity0, Real qmDensity1, Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    vela::ArithmeticReal quasiFermiDrop, Real Vt, Real coef, FermiSgFormulation formulation)
{
    if (!(qmDensity1 > 0.0) || !(Vt > 0.0) || coef == 0.0 ||
        quasiFermiDrop == VELA_ARITHMETIC_LITERAL(0.0)) {
        return VELA_ARITHMETIC_LITERAL(0.0);
    }
    if (formulation == FermiSgFormulation::EffectivePotential)
        return effectivePotentialFermiFluxExtended(qmDensity0, qmDensity1, quasiFermiDrop / static_cast<vela::ArithmeticReal>(Vt), coef);
    const vela::ArithmeticReal factor = static_cast<vela::ArithmeticReal>(
        sgGeneralizedEinsteinFactor(
            fermiDensity0, fermiDensity1, eta0, eta1));
    const vela::ArithmeticReal thermalVoltage = static_cast<vela::ArithmeticReal>(Vt);
    const vela::ArithmeticReal argument = static_cast<vela::ArithmeticReal>(driftPotential) /
        (thermalVoltage * factor);
    const vela::ArithmeticReal logarithmicLeftOverRight =
        quasiFermiDrop / (thermalVoltage * factor);
    return detail::sgFactorizedPopulationDifference(
        static_cast<ArithmeticReal>(coef) * factor * longDoubleBernoulli(argument) *
            static_cast<ArithmeticReal>(qmDensity1),
        logarithmicLeftOverRight, [](ArithmeticReal x) { return std::expm1(x); });
}

vela::ArithmeticReal sgHoleFermiDiracContinuityFluxLongDouble(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    vela::ArithmeticReal quasiFermiDrop, Real Vt, Real coef, FermiSgFormulation formulation)
{
    if (!(p1 > 0.0) || !(Vt > 0.0) || coef == 0.0 ||
        quasiFermiDrop == VELA_ARITHMETIC_LITERAL(0.0)) {
        return VELA_ARITHMETIC_LITERAL(0.0);
    }
    if (formulation == FermiSgFormulation::EffectivePotential)
        return effectivePotentialFermiFluxExtended(p0, p1, -quasiFermiDrop / static_cast<vela::ArithmeticReal>(Vt), coef);
    const vela::ArithmeticReal factor = static_cast<vela::ArithmeticReal>(
        sgGeneralizedEinsteinFactor(p0, p1, eta0, eta1));
    const vela::ArithmeticReal thermalVoltage = static_cast<vela::ArithmeticReal>(Vt);
    const vela::ArithmeticReal argument = static_cast<vela::ArithmeticReal>(driftPotential) /
        (thermalVoltage * factor);
    const vela::ArithmeticReal logarithmicLeftOverRight =
        -quasiFermiDrop / (thermalVoltage * factor);
    return detail::sgFactorizedPopulationDifference(
        static_cast<ArithmeticReal>(coef) * factor * longDoubleBernoulli(-argument) *
            static_cast<ArithmeticReal>(p1),
        logarithmicLeftOverRight, [](ArithmeticReal x) { return std::expm1(x); });
}

double sgElectronFlux(double n0, double n1, double dpsi, double Vt,
                      double mu, double h)
{
    return -sgElectronContinuityFlux(n0, n1, dpsi, Vt, mu * Vt / h);
}

double sgHoleFlux(double p0, double p1, double dpsi, double Vt,
                  double mu, double h)
{
    return -sgHoleContinuityFlux(p0, p1, dpsi, Vt, mu * Vt / h);
}

} // namespace vela
