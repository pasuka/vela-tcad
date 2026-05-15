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

    /// Return a copy with ni and low-field mobilities scaled to temperature_K.
    Material atTemperature(Real targetTemperature_K) const;
};

} // namespace vela
