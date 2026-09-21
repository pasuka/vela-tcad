#pragma once
#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/core/PhysicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
// Local positivity projection: no node can reduce the step of another node.
// The representability floor is capped by oldDensity to preserve zero steps.
inline Real electrothermalProjectedDensity(Real oldDensity,Real relative,Real alpha) {
    if(!(oldDensity>0.) || !std::isfinite(oldDensity) || !std::isfinite(relative) ||
       !std::isfinite(alpha) || alpha<0. || alpha>1.)
        throw std::invalid_argument("Invalid local density projection");
    const Real target=oldDensity*(1.+alpha*relative);
    if(!std::isfinite(target))throw std::invalid_argument("Nonfinite density target");
    return std::max(target,std::max(.01*oldDensity,std::min(oldDensity,Real(1e-250))));
}
// Invert the SAME local statistics at the candidate psi/T. Keep the reference
// and return its QF increment, avoiding subtraction of two large physical QFs.
inline Real electrothermalDensityQf(const SiliconThermalResult& atCandidate,
                                   const SiliconThermalState& candidate,
                                   Real density, bool hole) {
    const Real dos=hole?atCandidate.Nv_m3.value:atCandidate.Nc_m3.value;
    const Real ratio=density/dos;
    if(!(ratio>0.) || !std::isfinite(ratio))
        throw std::invalid_argument("Density update must remain positive and finite");
    const Real eta=inverseFermiDiracHalf(ratio);
    const Real vt=constants::kb/constants::q*candidate.temperature_K;
    const Real qf=hole?candidate.holeQf_V+vt*(eta-atCandidate.holeEta.value):
        candidate.electronQf_V+vt*(atCandidate.electronEta.value-eta);
    if(!std::isfinite(qf))throw std::runtime_error("Nonfinite density-coordinate update");
    return qf;
}
}
