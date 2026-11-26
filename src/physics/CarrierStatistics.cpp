#include "vela/core/ArithmeticPrecision.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/experimental/JoyceDixonStatistics.h"
#include "vela/core/PhysicsCallCounters.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iterator>
#include <stdexcept>
#include <span>

namespace vela {

namespace {

constexpr Real SqrtPi = 1.7724538509055160273;
#include "FermiHalfCoefficients.inc"

template<int Order>
Real evaluateFermiHalf(Real eta)
{
    constexpr bool Derivative=Order>0;
    if (std::isnan(eta)) return eta;
    if (eta == -std::numeric_limits<Real>::infinity()) return 0.0;
    if (eta == std::numeric_limits<Real>::infinity())
        return Order==2?0.:std::numeric_limits<Real>::infinity();
    if (eta < -8.0) {
        // Convergent fugacity series: sum (-1)^(k+1) exp(k eta)/k^(3/2).
        // Differentiating term by term changes the denominator to sqrt(k).
        const Real z = std::exp(eta);
        Real sum = 0.0;
        for (int k = 6; k >= 1; --k) {
            const Real denominator = Order==2 ? 1./std::sqrt(Real(k)) : Derivative ? std::sqrt(Real(k))
                                               : Real(k) * std::sqrt(Real(k));
            sum = (k % 2 ? 1.0 : -1.0) / denominator + z * sum;
        }
        return z * sum;
    }
    if (eta >= 64.0) {
        const Real reciprocal = 1.0 / eta;
        const Real u = reciprocal * reciprocal;
        Real sum = 0.0;
        for (int k = 6; k >= 0; --k)
            sum = sum * u + FermiSommerfeld[k] * (Derivative ? 1.5 - 2.0*k : 1.0) * (Order==2 ? .5-2.*k : 1.);
        const Real prefactor = (4.0 / (3.0 * SqrtPi)) * std::sqrt(eta);
        if constexpr (Order==2) return prefactor * sum / eta;
        if constexpr (Derivative) return prefactor * sum;
        return std::min(std::numeric_limits<Real>::max(), (prefactor * sum) * eta);
    }
    std::size_t segment = 0;
    while (eta > FermiBounds[segment + 1]) ++segment;
    const Real scale = 2.0 / (FermiBounds[segment + 1] - FermiBounds[segment]);
    const Real t = (eta - FermiBounds[segment]) * scale - 1.0;
    const auto coefficients = FermiCoefficients[segment];
    Real b1 = 0.0, b2 = 0.0, d1 = 0.0, d2 = 0.0, c1=0., c2=0.;
    for (std::size_t k = coefficients.size() - 1; k > 0; --k) {
        if constexpr (Order==2) {const Real second=4.*d1+2.*t*c1-c2;c2=c1;c1=second;}
        if constexpr (Derivative) {
            const Real derivative = 2.0*b1 + 2.0*t*d1 - d2;
            d2 = d1; d1 = derivative;
        }
        const Real b = coefficients[k] + 2.0*t*b1 - b2;
        b2 = b1; b1 = b;
    }
    if constexpr (Order==2) return (2.*d1+t*c1-c2)*scale*scale;
    if constexpr (Derivative) return (b1 + t*d1 - d2) * scale;
    return coefficients[0] + t*b1 - b2;
}

#include "ExtendedFermiHalfCoefficients.inc"

using ExtendedJoyceDixon = experimental::BasicJoyceDixonStatistics<long double>;

template<int Order>
long double evaluateExtendedFermiHalf(long double eta)
{
    constexpr long double sqrtPi = 1.77245385090551602729816748334114518L;
    constexpr auto infinity = std::numeric_limits<long double>::infinity();
    if (std::isnan(eta)) return eta;
    if (eta == -infinity) return 0;
    if (eta == infinity) return Order == 2 ? 0 : infinity;
    if (eta < WideFermiBounds[0]) {
        const auto z = std::exp(eta);
        if (z < std::numeric_limits<long double>::epsilon() / 2) return z;
        long double sum = 0;
        for (int k = 14; k >= 1; --k) {
            const auto root = std::sqrt(static_cast<long double>(k));
            const auto coefficient = Order == 2 ? root : Order == 1 ? 1 / root : 1 / (k * root);
            sum = (k % 2 ? coefficient : -coefficient) + z * sum;
        }
        return z * sum;
    }
    if (eta >= WideFermiBounds[std::size(WideFermiBounds) - 1]) {
        const auto reciprocal = 1 / eta;
        const auto u = reciprocal * reciprocal;
        long double sum = 0;
        for (int k = static_cast<int>(std::size(WideFermiSommerfeld)) - 1; k >= 0; --k) {
            auto term = WideFermiSommerfeld[k];
            if constexpr (Order >= 1) term *= 1.5L - 2 * k;
            if constexpr (Order == 2) term *= .5L - 2 * k;
            sum = sum * u + term;
        }
        const auto prefactor = (4 / (3 * sqrtPi)) * std::sqrt(eta);
        if constexpr (Order == 2) return prefactor * sum / eta;
        if constexpr (Order == 1) return prefactor * sum;
        return (prefactor * sum) * eta;
    }
    std::size_t segment = 0;
    while (eta > WideFermiBounds[segment + 1]) ++segment;
    const auto scale = 2 / (WideFermiBounds[segment + 1] - WideFermiBounds[segment]);
    const auto t = (eta - WideFermiBounds[segment]) * scale - 1;
    const auto& coefficients = WideFermiCoefficients[segment];
    long double b1 = 0, b2 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;
    for (std::size_t k = std::size(coefficients) - 1; k > 0; --k) {
        if constexpr (Order == 2) {
            const auto second = 4*d1 + 2*t*c1 - c2;
            c2 = c1;
            c1 = second;
        }
        if constexpr (Order >= 1) {
            const auto derivative = 2*b1 + 2*t*d1 - d2;
            d2 = d1;
            d1 = derivative;
        }
        const auto b = coefficients[k] + 2*t*b1 - b2;
        b2 = b1;
        b1 = b;
    }
    if constexpr (Order == 2) return (2*d1 + t*c1 - c2) * scale * scale;
    if constexpr (Order == 1) return (b1 + t*d1 - d2) * scale;
    return coefficients[0] + t*b1 - b2;
}

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
}

// At eta < -700 the fugacity correction is far below binary64 precision,
// while exp(eta) may be subnormal or zero even when DOS*exp(eta) is normal.
// Combine exponents first only in that tail; retain ordinary arithmetic and
// the full Fermi integral elsewhere. The same expression is its eta derivative.
template<bool Derivative = false>
Real scaledFermiHalf(Real dos, Real eta,
                     CarrierStatisticsModel model = CarrierStatisticsModel::FermiDirac)
{
    if (model == CarrierStatisticsModel::FermiDiracJoyceDixon) {
        if constexpr (Derivative)
            return experimental::JoyceDixonStatistics::densityDerivativeEta(eta, dos);
        return experimental::JoyceDixonStatistics::density(eta, dos);
    }
    if (eta < -700.0) return std::exp(std::log(dos) + eta);
    if constexpr (Derivative) return dos * fermiDiracHalfDerivative(eta);
    return dos * fermiDiracHalf(eta);
}

void validateDensityArguments(Real ni, Real densityOfStates, Real Vt,
                              CarrierStatisticsModel model)
{
    if (!(Vt > 0.0) || !std::isfinite(Vt))
        throw std::invalid_argument("carrier statistics: thermal voltage must be positive and finite.");
    if (ni < 0.0 || !std::isfinite(ni))
        throw std::invalid_argument("carrier statistics: intrinsic density must be nonnegative and finite.");
    if (usesFermiDirac(model) && ni > 0.0 &&
        (!(densityOfStates > 0.0) || !std::isfinite(densityOfStates))) {
        throw std::invalid_argument(
            "carrier statistics: Fermi-Dirac statistics require a positive finite density of states.");
    }
}

Real boltzmannElectronEquilibrium(Real netDoping, Real ni)
{
    const Real half = 0.5 * netDoping;
    const Real root = std::hypot(half, ni);
    if (netDoping >= 0.0)
        return half + root;
    const Real p = root - half;
    return p > 0.0 ? ni * ni / p : 0.0;
}

} // namespace

CarrierStatisticsModel carrierStatisticsModel(const CarrierStatisticsConfig& config)
{
    if (config.transport != "generalized_einstein" && config.transport != "effective_band_sg")
        throw std::invalid_argument("carrier_statistics.transport must be 'generalized_einstein' or 'effective_band_sg'.");
    if (config.transport == "effective_band_sg" && config.model != "fermi_dirac" &&
        config.model != "fermi_dirac_joyce_dixon")
        throw std::invalid_argument("effective_band_sg transport requires Fermi-Dirac statistics.");
    if (config.model == "boltzmann")
        return CarrierStatisticsModel::Boltzmann;
    if (config.model == "fermi_dirac")
        return CarrierStatisticsModel::FermiDirac;
    if (config.model == "fermi_dirac_joyce_dixon")
        return CarrierStatisticsModel::FermiDiracJoyceDixon;
    throw std::invalid_argument(
        "CarrierStatisticsConfig.model must be 'boltzmann', 'fermi_dirac', or 'fermi_dirac_joyce_dixon'.");
}

bool usesFermiDirac(const CarrierStatisticsConfig& config)
{
    return usesFermiDirac(carrierStatisticsModel(config));
}

Real fermiDiracHalf(Real eta)
{
    ++physicsCallCounters.fermiDiracHalf;
    return evaluateFermiHalf<false>(eta);
}

Real scaledFermiDiracHalf(Real densityOfStates, Real eta)
{
    return scaledFermiHalf(densityOfStates, eta);
}

Real scaledFermiDiracHalf(Real densityOfStates, Real eta, CarrierStatisticsModel model)
{
    return scaledFermiHalf(densityOfStates, eta, model);
}

Real fermiDiracHalfDerivative(Real eta)
{
    ++physicsCallCounters.fermiDiracHalfDerivative;
    return evaluateFermiHalf<true>(eta);
}

Real fermiDiracHalfSecondDerivative(Real eta) { return evaluateFermiHalf<2>(eta); }

Real inverseFermiDiracHalf(Real value)
{
    ++physicsCallCounters.inverseFermiDiracHalf;
    if (!(value > 0.0) || !std::isfinite(value)) {
        if (value == std::numeric_limits<Real>::infinity())
            return std::numeric_limits<Real>::infinity();
        throw std::invalid_argument(
            "inverseFermiDiracHalf: argument must be positive and finite.");
    }

    // For tiny positive values log(value) already resolves the inverse to
    // floating-point precision; do not stop by an absolute density tolerance.
    if (value < 1.0e-12) return std::log(value);
    Real lower = std::log(value);
    const Real degenerateEstimate = std::pow(value, 2.0 / 3.0)
        * std::pow(0.75 * SqrtPi, 2.0 / 3.0);
    Real upper = std::max<Real>(2.0, 2.0 * degenerateEstimate);
    Real eta = value < 0.5 ? lower : degenerateEstimate;
    eta = std::clamp(eta, lower, upper);
    for (int iteration = 0; iteration < 80; ++iteration) {
        const Real evaluated = fermiDiracHalf(eta);
        const Real function = evaluated - value;
        if (std::abs(function / value) <= 4.0e-14) break;
        if (function > 0.0) upper = eta;
        else lower = eta;
        const Real derivative = fermiDiracHalfDerivative(eta);
        const Real trial = derivative > 0.0 ? eta - function / derivative
                                            : lower + 0.5 * (upper - lower);
        const Real next = trial > lower && trial < upper && std::isfinite(trial)
            ? trial : lower + 0.5 * (upper - lower);
        if (next == eta || upper - lower <= 4.0 * std::numeric_limits<Real>::epsilon()
                * std::max<Real>(1.0, std::abs(eta))) break;
        eta = next;
    }
    return eta;
}

long double extendedFermiDiracHalf(long double eta) { return evaluateExtendedFermiHalf<0>(eta); }
long double extendedFermiDiracHalfDerivative(long double eta) { return evaluateExtendedFermiHalf<1>(eta); }
long double extendedFermiDiracHalfSecondDerivative(long double eta) { return evaluateExtendedFermiHalf<2>(eta); }

std::pair<long double, long double> extendedFermiPopulation(
    long double densityOfStates, long double intrinsic, long double reducedPotential,
    CarrierStatisticsModel model)
{
    const auto ratio = intrinsic / densityOfStates;
    const auto eta = reducedPotential + (ratio >= std::numeric_limits<long double>::min() && std::isfinite(ratio)
        ? std::log(ratio) : std::log(intrinsic) - std::log(densityOfStates));
    static const auto logMin = std::log(std::numeric_limits<long double>::min()) + 8;
    if (eta < logMin) {
        const auto value = std::exp(std::log(densityOfStates) + eta);
        return {value, value};
    }
    if (model == CarrierStatisticsModel::FermiDiracJoyceDixon) {
        const auto law = ExtendedJoyceDixon::evaluate(eta);
        return {densityOfStates * law.value, densityOfStates * law.first};
    }
    return {densityOfStates * extendedFermiDiracHalf(eta),
            densityOfStates * extendedFermiDiracHalfDerivative(eta)};
}

long double extendedInverseFermiDiracHalf(long double value)
{
    if (value == std::numeric_limits<long double>::infinity()) return value;
    if (!(value > 0) || !std::isfinite(value))
        throw std::invalid_argument("extendedInverseFermiDiracHalf: argument must be positive and finite.");

    const auto logarithm = std::log(value);
    constexpr auto epsilon = std::numeric_limits<long double>::epsilon();
    if (value < epsilon) return logarithm;

    auto lower = logarithm;
    // The Maxwell-Boltzmann limit is a lower bound. A Sommerfeld
    // estimate brackets the upper end even when the ratio is very large.
    auto estimate = std::pow(value, 2.L / 3) * 1.20899396551235231651L;
    auto upper = std::max(2.L, 2 * estimate);
    auto eta = value < .5L ? lower : estimate;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const auto residual = extendedFermiDiracHalf(eta) - value;
        if (std::abs(residual / value) <= 4 * epsilon) return eta;
        if (residual > 0) upper = eta;
        else lower = eta;
        const auto midpoint = lower + (upper - lower) / 2;
        const auto trial = eta - residual / extendedFermiDiracHalfDerivative(eta);
        if (trial == eta || midpoint == lower || midpoint == upper)
            return eta;
        eta = trial > lower && trial < upper ? trial : midpoint;
    }
    throw std::runtime_error("extendedInverseFermiDiracHalf: failed to converge.");
}

ExtendedFermiDegeneracy extendedFermiDegeneracy(long double densityOverDOS, CarrierStatisticsModel model)
{
    if (densityOverDOS == 0) return {};
    if (model == CarrierStatisticsModel::FermiDiracJoyceDixon) {
        return {-ExtendedJoyceDixon::primitiveLogGamma(densityOverDOS),
                -densityOverDOS * ExtendedJoyceDixon::primitiveLogGammaSlope(densityOverDOS)};
    }
    const auto eta = extendedInverseFermiDiracHalf(densityOverDOS);
    if (!std::isfinite(eta))
        throw std::invalid_argument("Fermi degeneracy correction requires a finite ratio.");
    if (densityOverDOS < std::numeric_limits<long double>::epsilon()) {
        const auto dilute = densityOverDOS / (2 * std::sqrt(2.L));
        return {dilute, dilute};
    }
    if (eta >= WideFermiBounds[0]) {
        const auto slope = extendedFermiDiracHalfDerivative(eta);
        return {eta - std::log(densityOverDOS), (densityOverDOS - slope) / slope};
    }

    // F = z(1 + correction) and F' = z(1 + derivativeCorrection).
    // Work with the two small corrections so the dilute limit does not
    // subtract nearly equal full logarithms or slopes.
    const auto z = std::exp(eta);
    long double correction = 0, derivativeCorrection = 0, difference = 0;
    for (int k = 14; k >= 2; --k) {
        const auto coefficient = (k % 2 ? 1.L : -1.L) / (k * std::sqrt(static_cast<long double>(k)));
        correction = coefficient + z * correction;
        derivativeCorrection = k * coefficient + z * derivativeCorrection;
        difference = (1 - k) * coefficient + z * difference;
    }
    return {-std::log1p(z * correction), z * difference / (1 + z * derivativeCorrection)};
}

Real fermiDiracHalf(Real eta, CarrierStatisticsModel model) {
    return model == CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::evaluate(eta).value : fermiDiracHalf(eta);
}
Real fermiDiracHalfDerivative(Real eta, CarrierStatisticsModel model) {
    return model == CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::evaluate(eta).first : fermiDiracHalfDerivative(eta);
}
Real fermiDiracHalfSecondDerivative(Real eta, CarrierStatisticsModel model) {
    return model == CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::evaluate(eta).second : fermiDiracHalfSecondDerivative(eta);
}
Real inverseFermiDiracHalf(Real value, CarrierStatisticsModel model) {
    return model == CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::inverse(value) : inverseFermiDiracHalf(value);
}

double electronDensity(double ni, double psi, double phin, double Vt)
{
    const double arg = std::clamp((psi - phin) / Vt, -500.0, 500.0);
    return ni * std::exp(arg);
}

double holeDensity(double ni, double psi, double phip, double Vt)
{
    const double arg = std::clamp((phip - psi) / Vt, -500.0, 500.0);
    return ni * std::exp(arg);
}

Real electronDensity(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                     CarrierStatisticsModel model)
{
    validateDensityArguments(ni, Nc, Vt, model);
    if (ni == 0.0)
        return 0.0;
    if (!usesFermiDirac(model))
        return ni * limitedExp((psi - phin) / Vt);
    const Real eta = (psi - phin) / Vt + std::log(ni / Nc);
    return model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::density(eta,Nc) : scaledFermiHalf(Nc, eta);
}

Real holeDensity(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                 CarrierStatisticsModel model)
{
    validateDensityArguments(ni, Nv, Vt, model);
    if (ni == 0.0)
        return 0.0;
    if (!usesFermiDirac(model))
        return ni * limitedExp((phip - psi) / Vt);
    const Real eta = (phip - psi) / Vt + std::log(ni / Nv);
    return model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::density(eta,Nv) : scaledFermiHalf(Nv, eta);
}

Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  CarrierStatisticsModel model)
{
    const Real n = electronDensity(ni, Nc, psi, phin, Vt, model);
    if (!usesFermiDirac(model) || n == 0.0)
        return n;
    const Real eta = (psi - phin) / Vt + std::log(ni / Nc);
    return model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::densityDerivativeEta(eta,Nc) : scaledFermiHalf<true>(Nc, eta);
}

Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              CarrierStatisticsModel model)
{
    const Real p = holeDensity(ni, Nv, psi, phip, Vt, model);
    if (!usesFermiDirac(model) || p == 0.0)
        return p;
    const Real eta = (phip - psi) / Vt + std::log(ni / Nv);
    return model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::densityDerivativeEta(eta,Nv) : scaledFermiHalf<true>(Nv, eta);
}

Real electronDensity(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                     const CarrierStatisticsConfig& config)
{
    return electronDensity(ni, Nc, psi, phin, Vt, carrierStatisticsModel(config));
}

Real holeDensity(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                 const CarrierStatisticsConfig& config)
{
    return holeDensity(ni, Nv, psi, phip, Vt, carrierStatisticsModel(config));
}

Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  const CarrierStatisticsConfig& config)
{
    return electronDensityDerivativeEta(
        ni, Nc, psi, phin, Vt, carrierStatisticsModel(config));
}

Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              const CarrierStatisticsConfig& config)
{
    return holeDensityDerivativeEta(
        ni, Nv, psi, phip, Vt, carrierStatisticsModel(config));
}

Real electronQuasiFermiPotential(Real ni, Real Nc, Real psi, Real n, Real Vt,
                                 const CarrierStatisticsConfig& config)
{
    const CarrierStatisticsModel model = carrierStatisticsModel(config);
    validateDensityArguments(ni, Nc, Vt, model);
    if (!(n > 0.0) || ni == 0.0)
        throw std::invalid_argument("electronQuasiFermiPotential: densities must be positive.");
    if (!usesFermiDirac(model))
        return psi - Vt * std::log(n / ni);
    const Real eta=model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::inverseDensity(n,Nc) : inverseFermiDiracHalf(n / Nc, model);
    return psi - Vt * (eta - std::log(ni / Nc));
}

Real holeQuasiFermiPotential(Real ni, Real Nv, Real psi, Real p, Real Vt,
                             const CarrierStatisticsConfig& config)
{
    const CarrierStatisticsModel model = carrierStatisticsModel(config);
    validateDensityArguments(ni, Nv, Vt, model);
    if (!(p > 0.0) || ni == 0.0)
        throw std::invalid_argument("holeQuasiFermiPotential: densities must be positive.");
    if (!usesFermiDirac(model))
        return psi + Vt * std::log(p / ni);
    const Real eta=model==CarrierStatisticsModel::FermiDiracJoyceDixon
        ? experimental::JoyceDixonStatistics::inverseDensity(p,Nv) : inverseFermiDiracHalf(p / Nv, model);
    return psi + Vt * (eta - std::log(ni / Nv));
}

EquilibriumCarrierState equilibriumCarrierState(
    Real netDoping, Real ni, Real Nc, Real Nv, Real Vt,
    CarrierStatisticsModel model)
{
    validateDensityArguments(ni, Nc, Vt, model);
    validateDensityArguments(ni, Nv, Vt, model);
    if (ni == 0.0)
        return {};

    if (!usesFermiDirac(model)) {
        const Real n = boltzmannElectronEquilibrium(netDoping, ni);
        const Real p = n > 0.0 ? ni * ni / n : 0.0;
        const Real potential = n > 0.0 ? Vt * std::log(n / ni) : 0.0;
        return {potential, n, p};
    }

    ++physicsCallCounters.equilibriumStateSolves;

    // The arguments are validated once above and ni != 0 here, so the
    // per-call validation and the log(ni/Nc), log(ni/Nv) evaluations inside
    // electronDensity()/holeDensity() are hoisted out of the scalar solve.
    // Ordinary-range arithmetic and ordering remain unchanged; the extreme
    // dilute tail shares the scaled evaluation used by the public functions.
    const Real logNiNc = std::log(ni / Nc);
    const Real logNiNv = std::log(ni / Nv);
    const auto electronEta = [&](Real potential) {
        return (potential - 0.0) / Vt + logNiNc;
    };
    const auto holeEta = [&](Real potential) {
        return (0.0 - potential) / Vt + logNiNv;
    };
    const auto charge = [&](Real potential) {
        return scaledFermiHalf(Nc, electronEta(potential), model)
            - scaledFermiHalf(Nv, holeEta(potential), model)
            - netDoping;
    };

    Real lower = -0.25;
    Real upper = 0.25;
    Real lowerCharge = charge(lower);
    while (lowerCharge > 0.0 && lower > -100.0) {
        lower *= 2.0;
        lowerCharge = charge(lower);
    }
    Real upperCharge = charge(upper);
    while (upperCharge < 0.0 && upper < 100.0) {
        upper *= 2.0;
        upperCharge = charge(upper);
    }
    if (!(lowerCharge <= 0.0 && upperCharge >= 0.0))
        throw std::runtime_error("equilibriumCarrierState: failed to bracket charge-neutral potential.");

    Real potential = 0.5 * (lower + upper);
    Real n = 0.0;
    Real p = 0.0;
    bool haveFinalDensities = false;
    for (int iteration = 0; iteration < 100; ++iteration) {
        ++physicsCallCounters.equilibriumStateIterations;
        const Real etaN = electronEta(potential);
        const Real etaP = holeEta(potential);
        n = scaledFermiHalf(Nc, etaN, model);
        p = scaledFermiHalf(Nv, etaP, model);
        const Real residual = n - p - netDoping;
        if (residual > 0.0)
            upper = potential;
        else
            lower = potential;
        if (std::abs(upper - lower) <= 2.0e-14 * std::max<Real>(1.0, std::abs(potential))) {
            haveFinalDensities = true;
            break;
        }

        // Same values as electronDensityDerivativeEta()/holeDensityDerivativeEta()
        // at this potential, reusing the densities evaluated for the residual.
        const Real dn = n == 0.0 ? n : scaledFermiHalf<true>(Nc, etaN, model);
        const Real dp = p == 0.0 ? p : scaledFermiHalf<true>(Nv, etaP, model);
        const Real derivative = (dn + dp) / Vt;
        const Real newton = derivative > 0.0 ? potential - residual / derivative
                                             : 0.5 * (lower + upper);
        potential = (newton > lower && newton < upper && std::isfinite(newton))
            ? newton : 0.5 * (lower + upper);
    }
    if (!haveFinalDensities) {
        n = scaledFermiHalf(Nc, electronEta(potential), model);
        p = scaledFermiHalf(Nv, holeEta(potential), model);
    }
    return {potential, n, p};
}

EquilibriumCarrierState equilibriumCarrierState(
    Real netDoping, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config)
{
    return equilibriumCarrierState(
        netDoping, ni, Nc, Nv, Vt, carrierStatisticsModel(config));
}

Real equilibriumCarrierProduct(
    Real n, Real p, Real ni, Real Nc, Real Nv, Real Vt,
    CarrierStatisticsModel model)
{
    if (!usesFermiDirac(model))
        return ni * ni;
    const EquilibriumCarrierState state = equilibriumCarrierState(
        n - p, ni, Nc, Nv, Vt, model);
    return state.n * state.p;
}

Real equilibriumCarrierProduct(
    Real n, Real p, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config)
{
    return equilibriumCarrierProduct(
        n, p, ni, Nc, Nv, Vt, carrierStatisticsModel(config));
}

namespace {

Real fermiLogDegeneracyFactor(Real density, Real densityOfStates, CarrierStatisticsModel model)
{
    if (!(density > 0.0) || !(densityOfStates > 0.0))return 0.0;
    const Real ratio=density/densityOfStates;
    if (!(ratio > 0.0) || !std::isfinite(ratio))return 0.0;
    if(model==CarrierStatisticsModel::FermiDiracJoyceDixon)
        return experimental::JoyceDixonStatistics::primitiveLogGamma(ratio);
    return std::log(ratio)-inverseFermiDiracHalf(ratio,model);
}

} // namespace

GeneralizedSrhCarrierState generalizedSrhCarrierState(
    Real n, Real p, Real ni, Real Nc, Real Nv,
    Real quasiFermiSplitting_V, Real Vt,
    CarrierStatisticsModel model)
{
    if (!(Vt > 0.0) || !std::isfinite(Vt))
        throw std::invalid_argument(
            "generalizedSrhCarrierState: Vt must be finite and positive.");

    GeneralizedSrhCarrierState state;
    const Real logN=usesFermiDirac(model)?fermiLogDegeneracyFactor(n,Nc,model):0.;
    const Real logP=usesFermiDirac(model)?fermiLogDegeneracyFactor(p,Nv,model):0.;
    state.electronDegeneracy=std::exp(logN);
    state.holeDegeneracy=std::exp(logP);
    state.logEquilibriumProduct=logN+logP+2*std::log(std::abs(ni));
    state.normalizedSplitting=quasiFermiSplitting_V/Vt;
    state.equilibriumProduct=state.electronDegeneracy*state.holeDegeneracy*ni*ni;
    if(!std::isfinite(state.equilibriumProduct) || state.equilibriumProduct==0.)
        state.equilibriumProduct=static_cast<Real>(std::exp((vela::ArithmeticReal)state.logEquilibriumProduct));
    const Real split=state.normalizedSplitting;
    // Keep the accurate expm1 path near equilibrium. Outside that range,
    // combine logs before exponentiation instead of clipping independent factors.
    if(std::abs(split)<=500. && std::isfinite(state.equilibriumProduct) &&
       state.equilibriumProduct>=std::numeric_limits<Real>::min())
        state.excessProduct=state.equilibriumProduct*std::expm1(split);
    else if(ni==0. || split==0.)state.excessProduct=0.;
    else {
        const vela::ArithmeticReal magnitude=std::exp((vela::ArithmeticReal)state.logEquilibriumProduct+
            (split>0.?split:0.))*-std::expm1(-std::abs((vela::ArithmeticReal)split));
        state.excessProduct=std::copysign(static_cast<Real>(magnitude),split);
    }
    return state;
}

Real generalizedSrhExcessDirectionalDerivative(
    const GeneralizedSrhCarrierState& state,Real dLogGammaN,
    Real dLogGammaP,Real dNormalizedSplitting)
{
    using Wide=vela::ArithmeticReal;
    const Wide dLogProduct=Wide(dLogGammaN)+dLogGammaP;
    const auto scaledExp=[](Wide logarithm,Wide direction) {
        if(direction==0.)return Wide(0.);
        return std::copysign(std::exp(logarithm+std::log(std::abs(direction))),direction);
    };
    const Wide logQ=state.logEquilibriumProduct,split=state.normalizedSplitting;
    const Wide direction=split>=1.
        ? scaledExp(logQ+split,dLogProduct+dNormalizedSplitting)-scaledExp(logQ,dLogProduct)
        : (dLogProduct==0.?Wide(0.):Wide(state.excessProduct)*dLogProduct)+
          scaledExp(logQ+split,dNormalizedSplitting);
    return static_cast<Real>(direction);
}

GeneralizedSrhCarrierState generalizedSrhCarrierState(
    Real n, Real p, Real ni, Real Nc, Real Nv,
    Real quasiFermiSplitting_V, Real Vt,
    const CarrierStatisticsConfig& config)
{
    return generalizedSrhCarrierState(
        n, p, ni, Nc, Nv, quasiFermiSplitting_V, Vt,
        carrierStatisticsModel(config));
}


double intrinsicDensity(const Material& material, double temperature_K)
{
    return material.atTemperature(temperature_K).ni;
}

} // namespace vela
