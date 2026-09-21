#pragma once

#include "vela/core/Types.h"
#include <functional>
#include <string>

namespace vela::detail {

struct DCSweepStepControlConfig {
    Real start = 0.0;
    Real stop = 0.0;
    Real step = 0.0;
    Real initialStep = 0.0;
    Real minStep = 0.0;
    Real maxStep = 0.0;
    Real growthFactor = 1.0;
    std::string growthMode = "fixed";
    int newtonIterationLimit = 0;
    Real shrinkFactor = 0.5;
    int maxRetries = 5;
    bool stopOnFailure = true;
    std::function<bool()> stopRequested;
};

struct DCSweepStepControlState {
    Real adaptiveStep = 0.0;
    bool initialized = false;
};

struct DCSweepStepControlEvent {
    Real voltage = 0.0;
    bool converged = false;
    Real attemptedStep = 0.0;
    Real acceptedStep = 0.0;
    int retryCount = 0;
    std::string failureReason;
    int newtonIterations = -1;
    Real growthFactor = 1.0;
    Real nextStepMagnitude = 0.0;
};

struct DCSweepStepAttemptResult {
    bool converged = false;
    int newtonIterations = -1;
    bool recovered = false;

    // Preserve existing bool callbacks in fixed-growth users of the controller.
    DCSweepStepAttemptResult(bool ok, int iterations = -1, bool usedRecovery = false)
        : converged(ok), newtonIterations(iterations), recovered(usedRecovery) {}
};

using DCSweepStepAttempt =
    std::function<DCSweepStepAttemptResult(Real voltage, Real attemptedStep, int retryCount)>;
using DCSweepStepRecorder = std::function<void(const DCSweepStepControlEvent& event)>;

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record,
                           DCSweepStepControlState* state = nullptr);

} // namespace vela::detail
