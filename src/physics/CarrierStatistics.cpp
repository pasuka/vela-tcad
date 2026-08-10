#include "vela/physics/CarrierStatistics.h"
#include "vela/core/PhysicsCallCounters.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace vela {

namespace {

constexpr Real SqrtPi = 1.7724538509055160273;
constexpr Real BednarczykCoefficient = 0.75 * SqrtPi;

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
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
    if (config.model == "boltzmann")
        return CarrierStatisticsModel::Boltzmann;
    if (config.model == "fermi_dirac")
        return CarrierStatisticsModel::FermiDirac;
    throw std::invalid_argument(
        "CarrierStatisticsConfig.model must be 'boltzmann' or 'fermi_dirac'.");
}

bool usesFermiDirac(const CarrierStatisticsConfig& config)
{
    return usesFermiDirac(carrierStatisticsModel(config));
}

Real fermiDiracHalf(Real eta)
{
    ++physicsCallCounters.fermiDiracHalf;
    if (!std::isfinite(eta))
        return eta < 0.0 ? 0.0 : std::numeric_limits<Real>::max();
    if (eta < -40.0)
        return std::exp(eta);

    const Real exponential = std::exp(std::max(-eta, -700.0));
    const Real shifted = eta + 1.0;
    const Real gaussian = std::exp(-0.17 * shifted * shifted);
    const Real v = std::pow(eta, 4.0) + 50.0
        + 33.6 * eta * (1.0 - 0.68 * gaussian);
    const Real denominator = exponential
        + BednarczykCoefficient * std::pow(v, -0.375);
    return 1.0 / denominator;
}

Real fermiDiracHalfDerivative(Real eta)
{
    ++physicsCallCounters.fermiDiracHalfDerivative;
    if (!std::isfinite(eta))
        return 0.0;
    if (eta < -40.0)
        return std::exp(eta);

    const Real exponential = std::exp(std::max(-eta, -700.0));
    const Real shifted = eta + 1.0;
    const Real gaussian = std::exp(-0.17 * shifted * shifted);
    const Real bracket = 1.0 - 0.68 * gaussian;
    const Real bracketDerivative = 0.2312 * shifted * gaussian;
    const Real v = std::pow(eta, 4.0) + 50.0 + 33.6 * eta * bracket;
    const Real vDerivative = 4.0 * eta * eta * eta
        + 33.6 * (bracket + eta * bracketDerivative);
    const Real power = std::pow(v, -0.375);
    const Real denominator = exponential + BednarczykCoefficient * power;
    const Real denominatorDerivative = -exponential
        - 0.375 * BednarczykCoefficient * std::pow(v, -1.375) * vDerivative;
    return std::max<Real>(0.0, -denominatorDerivative / (denominator * denominator));
}

Real inverseFermiDiracHalf(Real value)
{
    ++physicsCallCounters.inverseFermiDiracHalf;
    if (!(value > 0.0) || !std::isfinite(value)) {
        if (value == std::numeric_limits<Real>::infinity())
            return std::numeric_limits<Real>::infinity();
        throw std::invalid_argument(
            "inverseFermiDiracHalf: argument must be positive and finite.");
    }

    Real lower = -500.0;
    Real upper = value < 1.0
        ? 2.0
        : std::max<Real>(2.0, std::pow(value / 0.752252778063675, 2.0 / 3.0) * 1.5);
    while (fermiDiracHalf(upper) < value && upper < 1.0e6)
        upper *= 2.0;

    Real eta = value < 0.5
        ? std::log(value)
        : std::pow(value / 0.752252778063675, 2.0 / 3.0);
    eta = std::clamp(eta, lower, upper);
    for (int iteration = 0; iteration < 80; ++iteration) {
        const Real function = fermiDiracHalf(eta) - value;
        if (function > 0.0)
            upper = eta;
        else
            lower = eta;
        if (std::abs(function) <= 2.0e-13 * std::max<Real>(1.0, value))
            break;

        const Real derivative = fermiDiracHalfDerivative(eta);
        const Real newton = derivative > 0.0 ? eta - function / derivative
                                             : 0.5 * (lower + upper);
        eta = (newton > lower && newton < upper && std::isfinite(newton))
            ? newton : 0.5 * (lower + upper);
    }
    return eta;
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
    return Nc * fermiDiracHalf(eta);
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
    return Nv * fermiDiracHalf(eta);
}

Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  CarrierStatisticsModel model)
{
    const Real n = electronDensity(ni, Nc, psi, phin, Vt, model);
    if (!usesFermiDirac(model) || n == 0.0)
        return n;
    const Real eta = (psi - phin) / Vt + std::log(ni / Nc);
    return Nc * fermiDiracHalfDerivative(eta);
}

Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              CarrierStatisticsModel model)
{
    const Real p = holeDensity(ni, Nv, psi, phip, Vt, model);
    if (!usesFermiDirac(model) || p == 0.0)
        return p;
    const Real eta = (phip - psi) / Vt + std::log(ni / Nv);
    return Nv * fermiDiracHalfDerivative(eta);
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
    return psi - Vt * (inverseFermiDiracHalf(n / Nc) - std::log(ni / Nc));
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
    return psi + Vt * (inverseFermiDiracHalf(p / Nv) - std::log(ni / Nv));
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
    // The remaining arithmetic and its ordering are unchanged, so every
    // intermediate value stays bit-identical to the previous implementation.
    const Real logNiNc = std::log(ni / Nc);
    const Real logNiNv = std::log(ni / Nv);
    const auto electronEta = [&](Real potential) {
        return (potential - 0.0) / Vt + logNiNc;
    };
    const auto holeEta = [&](Real potential) {
        return (0.0 - potential) / Vt + logNiNv;
    };
    const auto charge = [&](Real potential) {
        return Nc * fermiDiracHalf(electronEta(potential))
            - Nv * fermiDiracHalf(holeEta(potential))
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
        n = Nc * fermiDiracHalf(etaN);
        p = Nv * fermiDiracHalf(etaP);
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
        const Real dn = n == 0.0 ? n : Nc * fermiDiracHalfDerivative(etaN);
        const Real dp = p == 0.0 ? p : Nv * fermiDiracHalfDerivative(etaP);
        const Real derivative = (dn + dp) / Vt;
        const Real newton = derivative > 0.0 ? potential - residual / derivative
                                             : 0.5 * (lower + upper);
        potential = (newton > lower && newton < upper && std::isfinite(newton))
            ? newton : 0.5 * (lower + upper);
    }
    if (!haveFinalDensities) {
        n = Nc * fermiDiracHalf(electronEta(potential));
        p = Nv * fermiDiracHalf(holeEta(potential));
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


double intrinsicDensity(const Material& material, double temperature_K)
{
    return material.atTemperature(temperature_K).ni;
}

} // namespace vela
