#pragma once

#include "vela/core/Types.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>

namespace vela {

/**
 * @brief Debye scaling system for the drift-diffusion equations.
 *
 * All internal solver variables are kept in dimensionless form.
 * This class computes the reference scales from physical parameters and
 * provides helpers to convert between physical and scaled quantities.
 */
class ScalingSystem {
public:
    /**
     * @param temperature          Lattice temperature [K]
     * @param reference_conc       Reference carrier concentration [m^-3]
     * @param relative_permittivity Relative permittivity of the semiconductor
     * @param reference_mobility   Reference carrier mobility [m^2/V/s]
     */
    ScalingSystem(Real temperature,
                  Real reference_conc,
                  Real relative_permittivity,
                  Real reference_mobility);

    // ------------------------------------------------------------------
    // Reference scales (SI units unless stated otherwise)
    // ------------------------------------------------------------------

    Real V0()  const { return V0_;  }  ///< Potential scale [V]
    Real C0()  const { return C0_;  }  ///< Concentration scale [m^-3]
    Real L0()  const { return L0_;  }  ///< Length scale (Debye length) [m]
    Real mu0() const { return mu0_; }  ///< Mobility scale [m^2/V/s]
    Real D0()  const { return D0_;  }  ///< Diffusivity scale [m^2/s]
    Real t0()  const { return t0_;  }  ///< Time scale [s]
    Real J0()  const { return J0_;  }  ///< Current density scale [A/m^2]
    Real R0()  const { return R0_;  }  ///< Recombination rate scale [m^-3/s]

    // ------------------------------------------------------------------
    // Scale (physical → dimensionless)
    // ------------------------------------------------------------------
    Real scalePotential     (Real phi)  const { return phi  / V0_; }
    Real scaleConcentration (Real n)    const { return n    / C0_; }
    Real scaleLength        (Real x)    const { return x    / L0_; }
    Real scaleMobility      (Real mu)   const { return mu   / mu0_; }
    Real scaleCurrentDensity(Real J)    const { return J    / J0_; }
    Real scaleRecombination (Real R)    const { return R    / R0_; }

    // ------------------------------------------------------------------
    // Unscale (dimensionless → physical)
    // ------------------------------------------------------------------
    Real unscalePotential     (Real phi_s) const { return phi_s  * V0_; }
    Real unscaleConcentration (Real n_s)   const { return n_s    * C0_; }
    Real unscaleLength        (Real x_s)   const { return x_s    * L0_; }
    Real unscaleMobility      (Real mu_s)  const { return mu_s   * mu0_; }
    Real unscaleCurrentDensity(Real J_s)   const { return J_s    * J0_; }
    Real unscaleRecombination (Real R_s)   const { return R_s    * R0_; }

private:
    Real V0_, C0_, L0_, mu0_, D0_, t0_, J0_, R0_;
};

} // namespace vela
