#include "vela/core/ScalingSystem.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <stdexcept>

namespace vela {

using namespace constants;

ScalingSystem::ScalingSystem(Real temperature,
                             Real reference_conc,
                             Real relative_permittivity,
                             Real reference_mobility)
{
    if (temperature <= 0.0)
        throw std::invalid_argument("Temperature must be positive.");
    if (reference_conc <= 0.0)
        throw std::invalid_argument("Reference concentration must be positive.");
    if (relative_permittivity <= 0.0)
        throw std::invalid_argument("Relative permittivity must be positive.");
    if (reference_mobility <= 0.0)
        throw std::invalid_argument("Reference mobility must be positive.");

    // Thermal voltage [V]
    V0_ = kb * temperature / q;

    // Concentration scale [m^-3]
    C0_ = reference_conc;

    // Debye length [m]:  L0 = sqrt(eps0 * eps_r * kT / (q^2 * C0))
    L0_ = std::sqrt(eps0 * relative_permittivity * kb * temperature
                    / (q * q * C0_));

    // Mobility scale [m^2/V/s]
    mu0_ = reference_mobility;

    // Diffusivity scale [m^2/s]:  D0 = mu0 * V0
    D0_ = mu0_ * V0_;

    // Time scale [s]:  t0 = L0^2 / D0
    t0_ = (L0_ * L0_) / D0_;

    // Current density scale [A/m^2]:  J0 = q * D0 * C0 / L0
    J0_ = q * D0_ * C0_ / L0_;

    // Recombination rate scale [m^-3/s]:  R0 = D0 * C0 / L0^2
    R0_ = D0_ * C0_ / (L0_ * L0_);
}

} // namespace vela
