#pragma once

#include "vela/core/Types.h"
#include <cmath>
#include <limits>

namespace vela::detail {

// Preserve ordinary evaluation order. A finite zero derivative has exactly zero
// action, even if forming flux/mobility first would overflow. Check all inputs
// before this shortcut so that an invalid flux cannot be silently discarded.
inline Real mobilityChainProduct(Real flux, Real mobility, Real derivative,
                                 Real voltageScale = 1.0, Real response = 1.0)
{
    if (!std::isfinite(flux) || !std::isfinite(mobility) || mobility <= 0.0 ||
        !std::isfinite(derivative) || !std::isfinite(voltageScale) || voltageScale <= 0.0 ||
        !std::isfinite(response))
        return std::numeric_limits<Real>::quiet_NaN();
    if (derivative == 0.0 || response == 0.0)
        return 0.0;
    return flux / mobility / voltageScale * derivative * response;
}

} // namespace vela::detail
