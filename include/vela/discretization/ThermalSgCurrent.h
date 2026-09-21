#pragma once
#include "vela/physics/SiliconThermalPhysics.h"
namespace vela {
struct ThermalSgCurrent {
    Real current_A_per_m=0.;
    /// psi,fn,fp,T at node0, then at node1. Mobility held fixed here.
    std::array<Real,8> derivative{};
    Real mobilityDerivative=0.;
};
/// Experimental nonisothermal generalized SG conductance. Uses the arithmetic
/// edge temperature and DOS-normalized Fermi secant, with flat QF preserved.
/// At uniform T/Nc/Nv it recovers the existing generalized SG formula.
/// Positive current points node0 -> node1; geometricWeight is dual-length/edge-length.
/// This API alone does not qualify nonisothermal LDMOS transport.
ThermalSgCurrent thermalSgCurrent(const SiliconThermalState& a,const SiliconThermalResult& pa,
    const SiliconThermalState& b,const SiliconThermalResult& pb,
    Real mobility_m2_per_Vs,Real geometricWeight,bool electron);
} // namespace vela
