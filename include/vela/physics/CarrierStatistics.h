#pragma once

#include "vela/material/Material.h"

namespace vela {

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

/// Temperature-adjusted intrinsic density for a material using temperature_K.
double intrinsicDensity(const Material& material, double temperature_K);

} // namespace vela
