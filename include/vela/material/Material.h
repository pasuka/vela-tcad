#pragma once

#include "vela/core/Types.h"
#include <optional>
#include <string>

namespace vela {

/**
 * @brief Physical and electrical properties of a semiconductor or insulator.
 */
struct Material {
    std::string name;

    Real eps_r = 1.0;   ///< Relative permittivity [-]
    Real ni    = 0.0;   ///< Intrinsic carrier concentration [m^-3]
    Real mun   = 0.0;   ///< Electron mobility [m^2/V/s]
    Real mup   = 0.0;   ///< Hole mobility [m^2/V/s]

    std::optional<Real> bandgap_eV;           ///< Band gap energy [eV]
    std::optional<Real> electron_affinity_eV; ///< Electron affinity [eV]
    std::optional<Real> Nc_m3;                ///< Effective conduction-band DOS [m^-3]
    std::optional<Real> Nv_m3;                ///< Effective valence-band DOS [m^-3]
    std::optional<Real> temperature_K;        ///< Material parameter temperature [K]
    /// Optional Sentaurus potential-based density-gradient parameters. When
    /// absent, solver-level controls provide the backward-compatible values.
    std::optional<Real> electron_quantum_gamma;
    std::optional<Real> electron_quantum_dos_mass_ratio;
    /// Independent mass for the Eq. 231 gradient coefficient.  If absent,
    /// the DOS mass preserves the legacy shared-mass behavior.
    std::optional<Real> electron_quantum_coefficient_mass_ratio;

    /// Optional thermal properties reserved for the lattice-temperature
    /// equation.  Keeping them in the material contract prevents phase-B
    /// implementations from silently sourcing unit-ambiguous defaults.
    std::optional<Real> thermal_conductivity_W_per_m_K;
    std::optional<Real> specific_heat_J_per_kg_K;
    std::optional<Real> mass_density_kg_per_m3;

    /// Return a copy with ni and low-field mobilities scaled to temperature_K.
    Material atTemperature(Real targetTemperature_K) const;
};

} // namespace vela
