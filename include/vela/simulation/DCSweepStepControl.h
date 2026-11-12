#pragma once

#include "vela/core/Types.h"
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace vela::detail {

struct DCSweepStepControlEvent;
struct DCSweepStepControlConfig {
    /// Read-only observer of every attempt, including intermediate failed retries.
    std::function<void(const DCSweepStepControlEvent&)> attemptObserved;
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
    /// Quasistationary exits when division by Decrement crosses MinStep.
    bool stopBelowMinimum = false;
    /// Retain the proposed step after a first-attempt output-point landing
    /// that shortened the step.
    bool preserveStepAfterOutputClip = false;
    std::function<bool()> stopRequested;
    /// Optional state-dependent cap on the proposed step magnitude. Empty
    /// means no representable step respecting both state and minimum limits.
    std::function<std::optional<Real>(Real proposed, Real remaining)> stepLimit;
    /// Biases strictly between start and stop that the controller must land
    /// on exactly (merged into the nominal target sequence).  Points outside
    /// the swept interval are ignored. A nominal target at the same voltage
    /// also counts as an output point.
    std::vector<Real> mandatoryPoints;
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
    Real proposedStepMagnitude = 0.0;
    bool targetClipped = false;
    bool recovered = false;
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
