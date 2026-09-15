#pragma once
#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
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
