#pragma once
#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/core/PhysicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
// Evaluate the actual backward-Euler defect, not its linear prediction. Both
// vectors use the same fixed row scales. Acceptance is an intermediate-step
// decision only: a small defect does not qualify the steady residual.
struct ElectrothermalPseudoTrial {bool accepted=false;Real norm=0.,modelError=0.;};
inline ElectrothermalPseudoTrial electrothermalPseudoTrial(
    const VectorXd& initial, const VectorXd& defect, Real alpha) {
    if(initial.size()!=defect.size() || !initial.allFinite() || !defect.allFinite() ||
       !(alpha>0. && alpha<=1. && std::isfinite(alpha)))
        throw std::invalid_argument("Invalid pseudo transient defect trial");
    const Real before=initial.norm(),after=defect.norm();
    const Real error=(defect-(1.-alpha)*initial).norm()/std::max(alpha*before,Real(1e-300));
    if(!std::isfinite(before)||!std::isfinite(after)||!std::isfinite(error))
        throw std::invalid_argument("Nonfinite pseudo transient defect norm");
    return {before>0. && after<=(1.-1e-4*alpha)*before,after,error};
}
inline Real electrothermalStorageDifference(Real oldDensity,Real newDensity,Real area,bool hole) {
    if(!(oldDensity>=0. && newDensity>=0. && area>=0.) || !std::isfinite(oldDensity) ||
       !std::isfinite(newDensity) || !std::isfinite(area))
        throw std::invalid_argument("Invalid pseudo storage difference");
    return static_cast<Real>((hole?1.L:-1.L)*constants::q*area*
        (static_cast<long double>(newDensity)-oldDensity));
}
// Storage per unit device depth [C/m]. Its derivative divided by tau [s]
// has the same units and sign convention as the steady continuity Jacobian.
// Poisson and lattice temperature remain algebraic. A finite-exchange contact
// has a free hole continuity row, unlike an ideal Dirichlet contact.
inline std::array<ThermalQuantity,4> electrothermalCarrierStorage(
    const SiliconThermalResult& physics, Real continuityArea,
    const std::array<bool,2>& active) {
    if(!std::isfinite(continuityArea) || continuityArea<0.)
        throw std::invalid_argument("Invalid continuity storage area");
    std::array<ThermalQuantity,4> result{};
    for(int k=0;k<2;++k)if(active[k] && continuityArea>0.) {
        const auto& density=k==0?physics.electrons_m3:physics.holes_m3;
        const Real scale=(k==0?-1.:1.)*constants::q*continuityArea;
        auto& row=result[k+1];row.value=scale*density.value;
        for(int j=0;j<4;++j)row.derivative[j]=scale*density.derivative[j];
        if(!std::isfinite(row.value))throw std::invalid_argument("Nonfinite carrier storage");
        for(Real v:row.derivative)if(!std::isfinite(v))
            throw std::invalid_argument("Nonfinite carrier storage derivative");
    }
    return result;
}
}
