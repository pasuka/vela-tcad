#pragma once

namespace vela::constants {

// Elementary charge [C]
constexpr double q    = 1.602176634e-19;

// Boltzmann constant [J/K]
constexpr double kb   = 1.380649e-23;

// Permittivity of free space [F/m]
constexpr double eps0 = 8.8541878128e-12;

// Planck constant [J·s]
constexpr double h    = 6.62607015e-34;

// Electron rest mass [kg]
constexpr double m0   = 9.1093837015e-31;

// Reference temperature [K]
constexpr double T0   = 300.0;

// Thermal voltage at 300 K [V]:  kT/q
constexpr double Vt_300 = kb * T0 / q;

} // namespace vela::constants
