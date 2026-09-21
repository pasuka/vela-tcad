#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
// Use the successful physical advance, not an unexecuted portion clipped by
// an exact output target. Zero-bias qualification keeps the original startup.
inline double electrothermalActualStep(double planned,double actual,int updates,int growth,
                                      double minimum,double maximum) {
    if(!(minimum>0. && maximum>=minimum && planned>=minimum && planned<=maximum && actual>=0.) ||
       !std::isfinite(planned) || !std::isfinite(actual) || !std::isfinite(maximum) || updates<0 || growth<1)
        throw std::invalid_argument("Invalid accepted electrothermal step");
    double step=actual>0.?std::min(planned,actual):planned;
    if(updates<=growth)step*=1.5;
    else if(updates>20)step*=.5;
    return std::clamp(step,minimum,maximum);
}
}
