#include "vela/simulation/DCSweepStepControl.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela::detail {

namespace {

template<class Scalar>
void validateDCSweepStepControlConfig(const BasicDCSweepStepControlConfig<Scalar>& cfg)
{
    const auto requireFinite = [](Scalar value, const char* name) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string("DCSweep step control: ") + name +
                                        " must be finite.");
        }
    };

    requireFinite(cfg.start, "start");
    requireFinite(cfg.stop, "stop");
    requireFinite(cfg.step, "step");
    requireFinite(cfg.initialStep, "initialStep");
    requireFinite(cfg.minStep, "minStep");
    requireFinite(cfg.maxStep, "maxStep");
    requireFinite(cfg.growthFactor, "growthFactor");
    requireFinite(cfg.shrinkFactor, "shrinkFactor");
    if (cfg.growthMode != "fixed" && cfg.growthMode != "newton_iterations") {
        throw std::invalid_argument(
            "DCSweep step control: growthMode must be 'fixed' or 'newton_iterations'.");
    }
    if (cfg.growthMode == "newton_iterations" && cfg.newtonIterationLimit <= 0) {
        throw std::invalid_argument(
            "DCSweep step control: newtonIterationLimit must be positive for newton_iterations growth.");
    }

    if (cfg.step == 0.0)
        throw std::invalid_argument("DCSweep step control: step must be non-zero.");
    if (cfg.start != cfg.stop && cfg.start + cfg.step == cfg.start)
        throw std::invalid_argument("DCSweep step control: nominal step is not representable at start.");
    if ((cfg.stop - cfg.start) * cfg.step < 0.0) {
        throw std::invalid_argument(
            "DCSweep step control: step sign must move start toward stop.");
    }
    if (cfg.minStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: minStep must be positive.");
    if (cfg.maxStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: maxStep must be positive.");
    if (cfg.minStep > cfg.maxStep) {
        throw std::invalid_argument(
            "DCSweep step control: minStep must not exceed maxStep.");
    }
    if (cfg.initialStep < 0.0) {
        throw std::invalid_argument(
            "DCSweep step control: initialStep must be positive when specified.");
    }
    const Scalar resolvedInitialStep =
        cfg.initialStep > 0.0 ? cfg.initialStep : std::abs(cfg.step);
    if (resolvedInitialStep < cfg.minStep || resolvedInitialStep > cfg.maxStep) {
        throw std::invalid_argument(
            "DCSweep step control: initialStep must lie within [minStep, maxStep].");
    }
    if (cfg.growthFactor < 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: growthFactor must be at least 1.");
    }
    if (cfg.shrinkFactor <= 0.0 || cfg.shrinkFactor >= 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: shrinkFactor must be greater than 0 and less than 1.");
    }
    if (cfg.maxRetries < 0) {
        throw std::invalid_argument(
            "DCSweep step control: maxRetries must be non-negative.");
    }
}

} // namespace

template<class Scalar>
void runTypedDCSweepStepControl(
    const BasicDCSweepStepControlConfig<Scalar>& cfg,
    const std::type_identity_t<BasicDCSweepStepAttempt<Scalar>>& attempt,
    const std::type_identity_t<BasicDCSweepStepRecorder<Scalar>>& record,
    std::type_identity_t<BasicDCSweepStepControlState<Scalar>*> state)
{
    validateDCSweepStepControlConfig(cfg);

    Scalar previousVoltage = cfg.start;
    const Scalar direction = (cfg.step > 0.0) ? 1.0 : -1.0;
    // The wide controller consumes authoritative decimal targets. A fixed
    // voltage tolerance can swallow distinct extended80 checkpoints.
    constexpr Scalar tolerance = std::is_same_v<Scalar, Real> ? Scalar{1e-12} : Scalar{0};
    BasicDCSweepStepControlState<Scalar> localState;
    BasicDCSweepStepControlState<Scalar>& controlState = state != nullptr ? *state : localState;
    if (!controlState.initialized) {
        const Scalar resolvedInitialStep =
            cfg.initialStep > 0.0 ? cfg.initialStep : std::abs(cfg.step);
        controlState.adaptiveStep = std::min(resolvedInitialStep, cfg.maxStep);
        controlState.initialized = true;
    } else if (!std::isfinite(controlState.adaptiveStep) ||
               controlState.adaptiveStep <= 0.0) {
        throw std::invalid_argument(
            "DCSweep step control: persisted adaptiveStep must be finite and positive.");
    }
    Scalar& adaptiveStep = controlState.adaptiveStep;

    auto limitedTarget = [&](Scalar target, Scalar stepMagnitude) {
        const Scalar remaining = direction * (target - previousVoltage);
        if constexpr (!std::is_same_v<Scalar, Real>)
            if (stepMagnitude >= remaining)
                return target;
        return previousVoltage + direction * std::min(stepMagnitude, remaining);
    };

    auto advanceToward = [&](Scalar target, bool outputPoint) -> bool {
        int retryCount = 0;
        Scalar trialStep = std::min(adaptiveStep, cfg.maxStep);
        Scalar lastAttempted = 0.0;
        Scalar lastCandidate = previousVoltage;

        while (true) {
            const Scalar remaining = direction * (target - previousVoltage);
            if (remaining <= tolerance)
                return true;

            Scalar stepMagnitude = std::min(trialStep, remaining);
            if (cfg.stepLimit) {
                const auto limited = cfg.stepLimit(stepMagnitude, remaining);
                if (!limited || !std::isfinite(*limited) || *limited <= 0.0 ||
                    *limited > stepMagnitude * (1.0 + 1e-12) ||
                    *limited < std::min(cfg.minStep, remaining)) {
                    record({previousVoltage, false, 0.0, 0.0, retryCount,
                            "state_step_limit_below_minimum_or_invalid"});
                    return false;
                }
                stepMagnitude = std::min(stepMagnitude, *limited);
            }
            const Scalar candidate = limitedTarget(target, stepMagnitude);
            const Scalar attemptedStep = candidate - previousVoltage;
            if (candidate == previousVoltage) {
                record({candidate, false, 0.0, 0.0, retryCount, "step_not_representable"});
                return false;
            }
            const DCSweepStepAttemptResult outcome = attempt(candidate, attemptedStep, retryCount);
            lastAttempted = attemptedStep;
            lastCandidate = candidate;

            if (outcome.converged) {
                Scalar growth = cfg.growthFactor;
                Scalar growthBase = cfg.stepLimit ? std::abs(attemptedStep) : trialStep;
                if (cfg.growthMode == "newton_iterations") {
                    if (outcome.newtonIterations < 0) {
                        throw std::invalid_argument(
                            "DCSweep step control: accepted newton_iterations step requires "
                            "a non-negative Newton iteration count.");
                    }
                    const Scalar effort = static_cast<Scalar>(std::max(outcome.newtonIterations - 1, 0)) /
                        (0.75 * static_cast<Scalar>(cfg.newtonIterationLimit));
                    growth = outcome.recovered ? 1.0 :
                        1.0 + (cfg.growthFactor - 1.0) * std::max(Scalar(0), Scalar(1) - effort);
                    // Work was measured on the actual (possibly clipped) step.
                    // Keep the legacy unclipped-base behavior only in fixed mode.
                    growthBase = std::abs(attemptedStep);
                }
                const bool preserveProposal = cfg.preserveStepAfterOutputClip && outputPoint &&
                    retryCount == 0 && remaining < trialStep - tolerance && !cfg.stepLimit;
                const Scalar nextStep = preserveProposal ? trialStep : cfg.growthMode == "newton_iterations"
                    ? std::clamp(growthBase * growth, cfg.minStep, cfg.maxStep)
                    : std::min(cfg.maxStep, growthBase * growth);
                BasicDCSweepStepControlEvent<Scalar> event{candidate,true,attemptedStep,attemptedStep,retryCount,{},
                    outcome.newtonIterations,preserveProposal ? 1.0 : growth,
                    nextStep,trialStep,remaining<trialStep,outcome.recovered};
                if(cfg.attemptObserved)cfg.attemptObserved(event);
                record(event);
                previousVoltage = candidate;
                adaptiveStep = nextStep;
                return true;
            }

            if(cfg.attemptObserved) {
                const Scalar next=stepMagnitude*cfg.shrinkFactor;
                const bool blocked=retryCount>=cfg.maxRetries || stepMagnitude<=cfg.minStep ||
                    (cfg.stopBelowMinimum&&next<cfg.minStep);
                cfg.attemptObserved({candidate,false,attemptedStep,0.,retryCount,"non_convergence",
                    outcome.newtonIterations,1.,blocked?0.:std::max(next,cfg.minStep),
                    trialStep,remaining<trialStep,outcome.recovered});
            }
            std::string failureReason;
            if (retryCount >= cfg.maxRetries) {
                failureReason = "non_convergence";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            const Scalar shrunken = stepMagnitude * cfg.shrinkFactor;
            // A failed step above the floor still gets one attempt at the
            // floor, even when geometric shrinking would skip over it. A
            // clipped final/mandatory interval at or below the floor cannot
            // be retried at a larger step without crossing its target.
            if (stepMagnitude <= cfg.minStep ||
                (cfg.stopBelowMinimum && shrunken < cfg.minStep)) {
                failureReason = "min_step_exhausted";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            trialStep = std::max(shrunken, cfg.minStep);
            ++retryCount;
        }
    };

    bool blockedByFailedStep = false;
    // Merge requested output voltages into the nominal target sequence.
    std::vector<Scalar> pendingMandatory;
    for (const Scalar point : cfg.mandatoryPoints) {
        if (!std::isfinite(point))
            throw std::invalid_argument("DCSweep step control: mandatory points must be finite.");
        if (direction * (point - cfg.start) > tolerance && direction * (cfg.stop - point) > tolerance)
            pendingMandatory.push_back(point);
    }
    std::sort(pendingMandatory.begin(), pendingMandatory.end(),
              [&](Scalar a, Scalar b) { return direction * a < direction * b; });
    std::size_t nextMandatory = 0;
    auto reachOutput = [&](Scalar target) -> bool {
        while (direction * (previousVoltage - target) < -tolerance) {
            if (!advanceToward(target, true))
                return false;
            if (cfg.stopRequested && cfg.stopRequested())
                return false;
        }
        return true;
    };
    auto reachMandatoryBefore = [&](Scalar target) -> bool {
        while (nextMandatory < pendingMandatory.size() &&
               direction * (pendingMandatory[nextMandatory] - target) < -tolerance) {
            const Scalar point = pendingMandatory[nextMandatory++];
            if (direction * (point - previousVoltage) <= tolerance)
                continue;
            if (!reachOutput(point))
                return false;
        }
        return true;
    };
    Scalar nominalTarget = cfg.start + cfg.step;
    while (!blockedByFailedStep && direction * (nominalTarget - cfg.stop) <= tolerance) {
        if (!reachMandatoryBefore(nominalTarget)) {
            if (cfg.stopOnFailure || (cfg.stopRequested && cfg.stopRequested()))
                return;
            blockedByFailedStep = true;
            break;
        }
        const bool isOutput = nextMandatory < pendingMandatory.size() &&
            std::abs(pendingMandatory[nextMandatory] - nominalTarget) <= tolerance;
        while (direction * (previousVoltage - nominalTarget) < -tolerance) {
            if (!advanceToward(nominalTarget, isOutput)) {
                if (cfg.stopOnFailure)
                    return;
                blockedByFailedStep = true;
                break;
            }
            if (cfg.stopRequested && cfg.stopRequested())
                return;
        }
        nominalTarget += cfg.step;
    }

    if (!blockedByFailedStep && !reachMandatoryBefore(cfg.stop)) {
        if (cfg.stopOnFailure || (cfg.stopRequested && cfg.stopRequested()))
            return;
        blockedByFailedStep = true;
    }
    while (!blockedByFailedStep && direction * (previousVoltage - cfg.stop) < -tolerance) {
        if (!advanceToward(cfg.stop, false)) {
            if (cfg.stopOnFailure)
                return;
            break;
        }
        if (cfg.stopRequested && cfg.stopRequested())
            return;
    }
}

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record,
                           DCSweepStepControlState* state)
{
    runTypedDCSweepStepControl(cfg, attempt, record, state);
}

template void runTypedDCSweepStepControl(
    const BasicDCSweepStepControlConfig<double>&,
    const BasicDCSweepStepAttempt<double>&, const BasicDCSweepStepRecorder<double>&,
    BasicDCSweepStepControlState<double>*);
template void runTypedDCSweepStepControl(
    const BasicDCSweepStepControlConfig<long double>&,
    const BasicDCSweepStepAttempt<long double>&, const BasicDCSweepStepRecorder<long double>&,
    BasicDCSweepStepControlState<long double>*);

} // namespace vela::detail
