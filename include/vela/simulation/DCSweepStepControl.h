#pragma once

#include "vela/core/Types.h"
#include <functional>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace vela::detail {

template<class Scalar>
struct BasicDCSweepStepControlEvent;

template<class Scalar>
struct BasicDCSweepStepControlConfig {
    /// Read-only observer of every attempt, including intermediate failed retries.
    std::function<void(const BasicDCSweepStepControlEvent<Scalar>&)> attemptObserved;
    Scalar start = 0.0;
    Scalar stop = 0.0;
    Scalar step = 0.0;
    Scalar initialStep = 0.0;
    Scalar minStep = 0.0;
    Scalar maxStep = 0.0;
    Scalar growthFactor = 1.0;
    std::string growthMode = "fixed";
    int newtonIterationLimit = 0;
    Scalar shrinkFactor = 0.5;
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
    std::function<std::optional<Scalar>(Scalar proposed, Scalar remaining)> stepLimit;
    /// Biases strictly between start and stop that the controller must land
    /// on exactly (merged into the nominal target sequence).  Points outside
    /// the swept interval are ignored. A nominal target at the same voltage
    /// also counts as an output point.
    std::vector<Scalar> mandatoryPoints;
};

template<class Scalar>
struct BasicDCSweepStepControlState {
    Scalar adaptiveStep = 0.0;
    bool initialized = false;
};

template<class Scalar>
struct BasicDCSweepStepControlEvent {
    Scalar voltage = 0.0;
    bool converged = false;
    Scalar attemptedStep = 0.0;
    Scalar acceptedStep = 0.0;
    int retryCount = 0;
    std::string failureReason;
    int newtonIterations = -1;
    Scalar growthFactor = 1.0;
    Scalar nextStepMagnitude = 0.0;
    Scalar proposedStepMagnitude = 0.0;
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

template<class Scalar>
using BasicDCSweepStepAttempt =
    std::function<DCSweepStepAttemptResult(Scalar voltage, Scalar attemptedStep, int retryCount)>;
template<class Scalar>
using BasicDCSweepStepRecorder =
    std::function<void(const BasicDCSweepStepControlEvent<Scalar>& event)>;

using DCSweepStepControlConfig = BasicDCSweepStepControlConfig<Real>;
using DCSweepStepControlState = BasicDCSweepStepControlState<Real>;
using DCSweepStepControlEvent = BasicDCSweepStepControlEvent<Real>;
using DCSweepStepAttempt = BasicDCSweepStepAttempt<Real>;
using DCSweepStepRecorder = BasicDCSweepStepRecorder<Real>;

template<class Scalar>
void runTypedDCSweepStepControl(
    const BasicDCSweepStepControlConfig<Scalar>& cfg,
    const std::type_identity_t<BasicDCSweepStepAttempt<Scalar>>& attempt,
    const std::type_identity_t<BasicDCSweepStepRecorder<Scalar>>& record,
    std::type_identity_t<BasicDCSweepStepControlState<Scalar>*> state = nullptr);

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record,
                           DCSweepStepControlState* state = nullptr);

} // namespace vela::detail
