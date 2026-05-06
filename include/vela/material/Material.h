#pragma once

#include "vela/core/Types.h"
#include <string>

namespace vela {

/**
 * @brief Physical and electrical properties of a semiconductor or insulator.
 */
struct Material {
    std::string name;

    Real eps_r = 1.0;   ///< Relative permittivity
    Real ni    = 0.0;   ///< Intrinsic carrier concentration [m^-3]
    Real mun   = 0.0;   ///< Electron mobility [m^2/V/s]
    Real mup   = 0.0;   ///< Hole mobility [m^2/V/s]
};

} // namespace vela
