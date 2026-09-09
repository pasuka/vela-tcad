#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/physics/CarrierStatistics.h"

#include <algorithm>
#include <cmath>
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
long double longDoubleBernoulli(long double value)
{
    if (std::abs(value) < 1.0e-12L)
        return 1.0L - value * 0.5L + value * value / 12.0L;
    if (value > 50.0L) {
        const long double expNegative = std::exp(-value);
        return value * expNegative / (1.0L - expNegative);
    }
    if (value < -50.0L) {
        const long double expValue = std::exp(value);
        return -value / (1.0L - expValue);
    }
    return value / std::expm1(value);
}

Real finiteRealFromLongDouble(long double value)
{
    const long double maximum =
        static_cast<long double>(std::numeric_limits<Real>::max());
    if (std::isnan(value))
        return 0.0;
    if (value >= maximum)
        return std::numeric_limits<Real>::max();
    if (value <= -maximum)
        return -std::numeric_limits<Real>::max();
    return static_cast<Real>(value);
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
    return coef * (w.b_minus * n0 - w.b_plus * n1);
}

Real sgHoleContinuityFlux(Real p0, Real p1, Real dpsi, Real Vt, Real coef)
{
    const SGEdgeWeights w = sgEdgeWeights(dpsi, Vt);
    return coef * (w.b_plus * p0 - w.b_minus * p1);
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

    const long double ni0Long = static_cast<long double>(ni0);
    const long double ni1Long = static_cast<long double>(ni1);
    const long double psi0Long = static_cast<long double>(psi0);
    const long double psi1Long = static_cast<long double>(psi1);
    const long double phin0Long = static_cast<long double>(phin0);
    const long double phin1Long = static_cast<long double>(phin1);
    const long double VtLong = static_cast<long double>(Vt);
    const long double coefLong = static_cast<long double>(coef);
    const long double endpointExponent0Long =
        (psi0Long - phin0Long) / VtLong;
    const long double endpointExponent1Long =
        (psi1Long - phin1Long) / VtLong;
    const long double clampedExponent0Long = std::clamp(
        endpointExponent0Long, -500.0L, 500.0L);
    const long double clampedExponent1Long = std::clamp(
        endpointExponent1Long, -500.0L, 500.0L);
    const long double n0Long = ni0Long * std::exp(clampedExponent0Long);
    const long double n1Long = ni1Long * std::exp(clampedExponent1Long);
    long double etaLong = (psi1Long - psi0Long) / VtLong;
    if (ni0 > 0.0 && ni1 > 0.0 && includeNiGradientDrift)
        etaLong += std::log(ni1Long / ni0Long);
    const long double leftTermLong =
        longDoubleBernoulli(-etaLong) * n0Long;
    const long double rightTermLong =
        longDoubleBernoulli(etaLong) * n1Long;
    const long double referenceTermScaleLong =
        std::abs(coefLong)
        * (std::abs(leftTermLong) + std::abs(rightTermLong));
    long double referenceFluxLong = 0.0L;
    if (!result.flatQuasiFermiShortCircuit) {
        if (ni0 > 0.0 && ni1 > 0.0) {
            long double logLeftOverRightLong =
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
    if (!(rightDensity > 0.0) || coefficient == 0.0 ||
        logarithmicLeftOverRight == 0.0)
        return 0.0;
    if (!std::isfinite(logarithmicLeftOverRight) || !std::isfinite(coefficient))
        return std::numeric_limits<Real>::quiet_NaN();

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

} // namespace

Real sgElectronFermiDiracContinuityFlux(
    Real n0, Real n1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    const Real factor = sgGeneralizedEinsteinFactor(n0, n1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight = (quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        n1, argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return coef * factor *
        (bernoulli(-argument) * n0 - bernoulli(argument) * n1);
}

Real sgElectronFermiDiracQuantumContinuityFlux(
    Real qmDensity0, Real qmDensity1,
    Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    const Real factor = sgGeneralizedEinsteinFactor(
        fermiDensity0, fermiDensity1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight =
        (quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        qmDensity1, argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return coef * factor * (
        bernoulli(-argument) * qmDensity0 -
        bernoulli(argument) * qmDensity1);
}

Real sgHoleFermiDiracContinuityFlux(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef)
{
    if (quasiFermi0 == quasiFermi1)
        return 0.0;
    if (!(Vt > 0.0))
        return 0.0;
    const Real factor = sgGeneralizedEinsteinFactor(p0, p1, eta0, eta1);
    const Real argument = driftPotential / (Vt * factor);
    const Real logLeftOverRight = -(quasiFermi1 - quasiFermi0) / (Vt * factor);
    const Real stable = stableGeneralizedFermiDifference(
        p1, -argument, logLeftOverRight, coef * factor);
    if (std::isfinite(stable))
        return stable;
    return coef * factor *
        (bernoulli(argument) * p0 - bernoulli(-argument) * p1);
}

long double sgElectronFermiDiracQuantumContinuityFluxLongDouble(
    Real qmDensity1, Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    long double quasiFermiDrop, Real Vt, Real coef)
{
    if (!(qmDensity1 > 0.0) || !(Vt > 0.0) || coef == 0.0 ||
        quasiFermiDrop == 0.0L) {
        return 0.0L;
    }
    const long double factor = static_cast<long double>(
        sgGeneralizedEinsteinFactor(
            fermiDensity0, fermiDensity1, eta0, eta1));
    const long double thermalVoltage = static_cast<long double>(Vt);
    const long double argument = static_cast<long double>(driftPotential) /
        (thermalVoltage * factor);
    const long double logarithmicLeftOverRight =
        quasiFermiDrop / (thermalVoltage * factor);
    return static_cast<long double>(coef) * factor *
        longDoubleBernoulli(argument) *
        static_cast<long double>(qmDensity1) *
        std::expm1(logarithmicLeftOverRight);
}

long double sgHoleFermiDiracContinuityFluxLongDouble(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    long double quasiFermiDrop, Real Vt, Real coef)
{
    if (!(p1 > 0.0) || !(Vt > 0.0) || coef == 0.0 ||
        quasiFermiDrop == 0.0L) {
        return 0.0L;
    }
    const long double factor = static_cast<long double>(
        sgGeneralizedEinsteinFactor(p0, p1, eta0, eta1));
    const long double thermalVoltage = static_cast<long double>(Vt);
    const long double argument = static_cast<long double>(driftPotential) /
        (thermalVoltage * factor);
    const long double logarithmicLeftOverRight =
        -quasiFermiDrop / (thermalVoltage * factor);
    return static_cast<long double>(coef) * factor *
        longDoubleBernoulli(-argument) *
        static_cast<long double>(p1) *
        std::expm1(logarithmicLeftOverRight);
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
