#pragma once

#include "vela/solver/ExtendedSparseLu.h"

namespace vela {

// Authoritative [psi, n, p] blocks in mesh order and in the solver's internal
// coordinates. Circuit state is in V and A/um, independent of those scales.
struct ExtendedPopulationState {
    ExtendedVector coordinates;
    long double potentialScale_V = 1;
    long double potentialOrigin_V = 0;
    long double concentrationScale_m3 = 1;
    long double deviceBias_V = 0;
    long double supplyBias_V = 0;
    long double drainCurrent_A_per_um = 0;
};

} // namespace vela
