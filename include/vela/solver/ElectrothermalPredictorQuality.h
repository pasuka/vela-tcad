#pragma once
#include <array>
#include <cmath>
#include <algorithm>

namespace vela::experimental {
// Common scaled norm, carrier row ratio, then the three electrical block norms.
// Floors are the original acceptance limits, not new relaxed solver tolerances.
inline bool improvesElectrothermalPredictor(const std::array<double,5>& trial,
                                            const std::array<double,5>& baseline,
                                            const std::array<double,5>& floors) {
    if(!(trial[0]<baseline[0]))return false;
    for(unsigned k=0;k<5;++k)
        if(!std::isfinite(trial[k]) || trial[k]<0. || !std::isfinite(baseline[k]) ||
           !std::isfinite(floors[k]) || floors[k]<0. || trial[k]>std::max(baseline[k],floors[k]))return false;
    return true;
}
}
