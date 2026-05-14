#pragma once

#include "vela/core/Types.h"
#include <functional>

namespace vela::detail {

struct DCSweepStepControlConfig {
    Real start = 0.0;
    Real stop = 0.0;
    Real step = 0.0;
    Real minStep = 0.0;
    Real maxStep = 0.0;
    Real growthFactor = 1.0;
    Real shrinkFactor = 0.5;
    int maxRetries = 5;
    bool stopOnFailure = true;
};

struct DCSweepStepControlEvent {
    Real voltage = 0.0;
    bool converged = false;
    Real attemptedStep = 0.0;
    Real acceptedStep = 0.0;
    int retryCount = 0;
};

using DCSweepStepAttempt =
    std::function<bool(Real voltage, Real attemptedStep, int retryCount)>;
using DCSweepStepRecorder = std::function<void(const DCSweepStepControlEvent& event)>;

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record);

} // namespace vela::detail
