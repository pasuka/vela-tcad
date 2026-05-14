#pragma once

#include "vela/core/Types.h"
#include <string>

namespace vela {

enum class CurveSweepMode {
    IV,
    CVQuasistatic,
    BVReverse,
};

struct CurveSweepStepDiagnostics {
    Real attemptedStep = 0.0;
    Real acceptedStep = 0.0;
    int retryCount = 0;
};

std::string toString(CurveSweepMode mode);
CurveSweepMode curveSweepModeFromString(const std::string& mode);

} // namespace vela
