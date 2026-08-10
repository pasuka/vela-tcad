#pragma once

#include "vela/material/Material.h"
#include <cstdint>
#include <string>

namespace vela {

struct CarrierStatisticsConfig {
    /// "boltzmann" or "fermi_dirac".  Boltzmann remains the compatibility default.
    std::string model = "boltzmann";
};

/**
 * @brief Resolved carrier-statistics model.
 *
 * `CarrierStatisticsConfig::model` is a string because it is part of the
 * configuration schema.  Hot assembly loops resolve it once and pass this
 * enumeration instead, so no `std::string` comparison happens per carrier
 * density evaluation.  The resolved overloads are numerically identical to the
 * configuration overloads; only the model lookup is hoisted.
 */
enum class CarrierStatisticsModel : std::uint8_t {
    Boltzmann,
    FermiDirac,
};

struct EquilibriumCarrierState {
    Real potential = 0.0; ///< Electrostatic potential relative to the intrinsic reference [V].
    Real n = 0.0;         ///< Electron density [m^-3].
    Real p = 0.0;         ///< Hole density [m^-3].
};

/// Validate `config.model` and return the resolved enumeration.
CarrierStatisticsModel carrierStatisticsModel(const CarrierStatisticsConfig& config);

bool usesFermiDirac(const CarrierStatisticsConfig& config);

constexpr bool usesFermiDirac(CarrierStatisticsModel model) noexcept
{
    return model == CarrierStatisticsModel::FermiDirac;
}

/**
 * @brief Normalized complete Fermi-Dirac integral F_{1/2}(eta).
 *
 * Uses the Bednarczyk analytic approximation (maximum relative error 0.4%).
 * The normalization is 2/sqrt(pi) times the defining integral, so the
 * non-degenerate limit is exp(eta).
 */
Real fermiDiracHalf(Real eta);
Real fermiDiracHalfDerivative(Real eta);
Real inverseFermiDiracHalf(Real value);

/**
 * @brief Boltzmann carrier statistics.
 *
 * Boltzmann (non-degenerate) approximation:
 *   n = ni * exp((psi - phin) / Vt)
 *   p = ni * exp((phip - psi) / Vt)
 *
 * Exponent arguments are clamped to [-500, 500] to prevent overflow.
 *
 * @param ni    Intrinsic carrier concentration [m^-3]
 * @param psi   Electrostatic potential [V]
 * @param phin  Electron quasi-Fermi potential [V]
 * @param phip  Hole quasi-Fermi potential [V]
 * @param Vt    Thermal voltage kT/q [V]
 */
double electronDensity(double ni, double psi, double phin, double Vt);
double holeDensity    (double ni, double psi, double phip, double Vt);

/// Carrier densities with an explicit statistics model and density of states.
Real electronDensity(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                     const CarrierStatisticsConfig& config);
Real holeDensity(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                 const CarrierStatisticsConfig& config);

/// Resolved-model overloads for hot assembly loops.
Real electronDensity(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                     CarrierStatisticsModel model);
Real holeDensity(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                 CarrierStatisticsModel model);

/// Derivatives with respect to the dimensionless reduced Fermi energy.
Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  const CarrierStatisticsConfig& config);
Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              const CarrierStatisticsConfig& config);

Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  CarrierStatisticsModel model);
Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              CarrierStatisticsModel model);

Real electronQuasiFermiPotential(Real ni, Real Nc, Real psi, Real n, Real Vt,
                                 const CarrierStatisticsConfig& config);
Real holeQuasiFermiPotential(Real ni, Real Nv, Real psi, Real p, Real Vt,
                             const CarrierStatisticsConfig& config);

/**
 * @brief Charge-neutral equilibrium state used by ideal Ohmic contacts.
 *
 * Solves n(psi)-p(psi)=netDoping numerically for Fermi-Dirac statistics.
 * The Boltzmann branch retains the closed-form legacy solution.
 */
EquilibriumCarrierState equilibriumCarrierState(
    Real netDoping, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config);

EquilibriumCarrierState equilibriumCarrierState(
    Real netDoping, Real ni, Real Nc, Real Nv, Real Vt,
    CarrierStatisticsModel model);

/// Equilibrium n0*p0 at the same local charge imbalance n-p.
Real equilibriumCarrierProduct(
    Real n, Real p, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config);

Real equilibriumCarrierProduct(
    Real n, Real p, Real ni, Real Nc, Real Nv, Real Vt,
    CarrierStatisticsModel model);

/// Temperature-adjusted intrinsic density for a material using temperature_K.
double intrinsicDensity(const Material& material, double temperature_K);

} // namespace vela
