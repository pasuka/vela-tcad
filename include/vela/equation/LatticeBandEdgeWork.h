#pragma once
#include "vela/core/Types.h"
#include <array>
namespace vela {
struct LatticeBandEdgeWork {
    Real power_W_per_m;
    /// Order: Jn, Jp [A/m], Ec0, Ec1, Ev0, Ev1 [eV].
    std::array<Real,6> derivative;
};
/// Default lattice heat without RecGenHeat/Peltier/Thermodynamic.
/// Currents are conventional signed edge currents oriented node0 -> node1,
/// per unit out-of-plane width. Energies in eV are numerically volts after /q.
/// Use the same current and its Jacobian as in continuity, not abs(J) or a
/// separately reconstructed drift-only current. Negative local work is allowed.
LatticeBandEdgeWork latticeBandEdgeWork(Real electronCurrent_A_per_m,
    Real holeCurrent_A_per_m,Real Ec0_eV,Real Ec1_eV,Real Ev0_eV,Real Ev1_eV);
} // namespace vela
