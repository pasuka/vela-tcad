#pragma once

#include "vela/core/Types.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace vela {

/// Configuration for pseudo-arclength continuation.
///
/// The defaults keep the feature disabled so existing decks that do not opt in
/// retain their voltage-parameterized stepping behavior unchanged.
struct PseudoArclengthConfig {
    bool enabled = false;
    /// Initial arclength step length (Delta s). Must be > 0 when enabled.
    Real initialStep = 0.0;
    /// Lower bound for the adaptive arclength step. Below this a step is a hard failure.
    Real minStep = 0.0;
    /// Upper bound for the adaptive arclength step.
    Real maxStep = 0.0;
    /// Multiplier applied to Delta s after a clean corrector convergence (>= 1).
    Real growthFactor = 1.2;
    /// Multiplier applied to Delta s after a corrector failure (0 < f < 1).
    Real shrinkFactor = 0.5;
    /// Maximum corrector (bordered-Newton) iterations per attempted step.
    int maxCorrectorIterations = 20;
    /// Corrector convergence tolerance on max(||F||_inf, |N|).
    Real correctorTolerance = 1.0e-8;
    /// Maximum number of shrink retries before the step is abandoned.
    int maxStepRetries = 8;
    /// Weight (theta) of the parameter component in the arclength norm.
    /// Must be > 0.
    Real parameterScale = 1.0;
    /// Weight applied to state-space inner products. A value of 0 selects
    /// the mesh-size-independent default 1 / x.size().
    Real stateWeight = 0.0;
    /// Initial damping factor for corrector backtracking. Must be in (0, 1].
    Real dampingFactor = 1.0;
    /// Maximum number of residual-monotone backtracking halvings per corrector step.
    int maxLineSearchSteps = 10;
    /// Relative residual increase tolerated by the corrector line search. Zero
    /// preserves strict monotonic decrease. This is useful with deliberately
    /// approximate (for example frozen-source) Jacobians.
    Real lineSearchRelativeIncreaseTolerance = 0.0;
    /// Optional absolute cap on a bordered-Newton continuation-parameter update.
    /// Zero disables the cap. This is a numerical trust region on lambda, not a
    /// clamp on the accepted branch or on the device state.
    Real maxParameterUpdate = 0.0;
};

/// A point on the solution branch: state vector x and continuation parameter lambda.
struct ArclengthState {
    VectorXd x;
    Real lambda = 0.0;
};

struct ArclengthTrial {
    ArclengthState state;
    Real effectiveAlpha = 1.;
    int protectedElectrons = 0, protectedHoles = 0;
};

/// Unit tangent direction in (x, lambda) space (already scaled by parameterScale).
struct ArclengthTangent {
    VectorXd xDot;
    Real lambdaDot = 0.0;
};

/// Value and first derivatives of a scalar functional evaluated on an
/// arclength state.  The load-line solver uses this for the terminal-current
/// row of a directly assembled bordered Newton system.
struct ArclengthScalarLinearization {
    Real value = 0.0;
    VectorXd stateDerivative;
    Real parameterDerivative = 0.0;
};

struct ArclengthScalarFunctional {
    std::function<Real(const VectorXd&, Real)> value;
    std::function<ArclengthScalarLinearization(const VectorXd&, Real)> linearize;
    std::function<ArclengthScalarLinearization(const VectorXd&,Real,const VectorXd&)> linearizeInColumns;
};

/// Outcome of a single attempted arclength step.
struct ArclengthStepResult {
    bool converged = false;
    ArclengthState state;
    Real arclengthStep = 0.0;   ///< Delta s actually accepted.
    int correctorIterations = 0;
    int retries = 0;
    Real residualNorm = 0.0;    ///< max(||F||_inf, |N|) at the accepted point.
    Real equationResidualNorm = 0.0; ///< ||F||_inf at the last corrector state.
    Real arclengthResidual = 0.0;    ///< Signed N at the last corrector state.
    Real bestTrialResidualNorm = std::numeric_limits<Real>::infinity();
    Real bestTrialEquationResidualNorm = std::numeric_limits<Real>::infinity();
    Real bestTrialArclengthResidual = std::numeric_limits<Real>::infinity();
    Real bestTrialAlpha = 0.0;
    Real smallestTrialAlpha = 0.0;
    int lineSearchAttempts = 0;
    Real rawStateUpdateNorm = 0.0;
    Real limitedStateUpdateNorm = 0.0;
    Real updateLimitChangeNorm = 0.0;
    Real rawParameterUpdate = 0.0;
    Real limitedParameterUpdate = 0.0;
    std::string failureReason;
};

/// Callbacks describing the nonlinear continuation system F(x, lambda) = 0.
///
/// All three callbacks must be supplied. The continuation core never assembles
/// the Jacobian directly; instead it delegates linear solves to solveJacobian so
/// the device-level integration can reuse the existing sparse Jacobian assembly
/// and factorization.
struct ArclengthSystem {
    /// Residual F(x, lambda), size n.
    std::function<VectorXd(const VectorXd&, Real)> residual;
    /// Parameter sensitivity dF/dlambda evaluated at (x, lambda), size n.
    std::function<VectorXd(const VectorXd&, Real)> parameterDerivative;
    /// Optional direct assembly of J=dF/dx.  Schur-complement continuation
    /// only needs solveJacobian, while a circuit border must be able to retain
    /// J inside one augmented sparse matrix even when J itself is singular.
    std::function<SparseMatrixd(const VectorXd& x, Real lambda)> jacobian;
    /// Assemble J diag(columns) directly, retaining the physical state x.
    /// Used only for depleted population columns in the bordered solve.
    std::function<SparseMatrixd(const VectorXd&,Real,const VectorXd&)> jacobianInColumns;
    std::function<SparseMatrixd(const VectorXd&,Real,const VectorXd&)> bootstrapJacobianInColumns;
    /// Solve J(x, lambda) * y = b where J = dF/dx. Writes the solution into y and
    /// returns false if the Jacobian is singular or the solve otherwise failed.
    std::function<bool(const VectorXd& x, Real lambda, const VectorXd& b, VectorXd& y)>
        solveJacobian;
    /// Optional batched solve for two right-hand sides at the same (x, lambda).
    /// Device integrations can use this to assemble and factor J only once.
    std::function<bool(const VectorXd& x,
                       Real lambda,
                       const VectorXd& b1,
                       const VectorXd& b2,
                       VectorXd& y1,
                       VectorXd& y2)>
        solveJacobianPair;
    /// Optional in-place limiter for the bordered-Newton update before line search.
    std::function<void(const VectorXd& x, VectorXd& deltaX, Real& deltaLambda)>
        limitUpdate;
    /// Coordinate-aware finite update, evaluated independently for every alpha.
    std::function<ArclengthTrial(const ArclengthState&, const VectorXd&, Real, Real)> makeTrial;
    /// Optional coordinate-aware secant predictor; invalid output rolls back.
    std::function<VectorXd(const VectorXd&, const VectorXd&, Real)> predictState;
    std::function<std::vector<int>(const VectorXd&, Real)> prescribedRows;
    std::function<bool(const VectorXd&, Real)> convergenceSatisfied;
    /// Positive physical column scales, including the final parameter column.
    std::function<VectorXd(const ArclengthState&)> columnScales;
    std::function<Real(const VectorXd&)> equationNorm;
    std::function<SparseMatrixd(const VectorXd&,Real)> bootstrapJacobian;
    bool automaticDerivativeRecovery=false;
};

/// Pseudo-arclength (Keller) continuation engine.
///
/// The engine implements a tangent predictor followed by a bordered-Newton
/// corrector that augments the system with one arclength constraint equation
///   N(x, lambda) = w*x_dot . (x - x0) + theta^2 * lambda_dot * (lambda - lambda0) - Delta s.
/// Because lambda is treated as an unknown, the augmented system remains
/// nonsingular at turning points where dlambda/d(arclength) -> 0, allowing the
/// branch to be traced across a fold that voltage-parameterized stepping cannot.
class PseudoArclengthContinuation {
public:
    PseudoArclengthContinuation(ArclengthSystem system, PseudoArclengthConfig config)
        : system_(std::move(system)), config_(std::move(config))
    {
        if(system_.makeTrial)
            throw std::invalid_argument("PseudoArclengthContinuation: coordinate-aware finite updates are currently qualified only for coupled load-line Newton");
        if (!system_.residual || !system_.parameterDerivative || !system_.solveJacobian) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: residual, parameterDerivative, and "
                "solveJacobian callbacks are all required.");
        }
        if (!(config_.parameterScale > 0.0) || !std::isfinite(config_.parameterScale)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: parameterScale must be finite and positive.");
        }
        if (config_.stateWeight < 0.0 || !std::isfinite(config_.stateWeight)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: stateWeight must be finite and non-negative.");
        }
        if (!(config_.dampingFactor > 0.0) || config_.dampingFactor > 1.0 ||
            !std::isfinite(config_.dampingFactor)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: dampingFactor must be finite and in (0, 1].");
        }
        if (config_.maxLineSearchSteps < 0) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: maxLineSearchSteps must be non-negative.");
        }
        if (config_.lineSearchRelativeIncreaseTolerance < 0.0 ||
            !std::isfinite(config_.lineSearchRelativeIncreaseTolerance)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: lineSearchRelativeIncreaseTolerance "
                "must be finite and non-negative.");
        }
        if (config_.maxParameterUpdate < 0.0 ||
            !std::isfinite(config_.maxParameterUpdate)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: maxParameterUpdate must be finite and "
                "non-negative.");
        }
    }

    const PseudoArclengthConfig& config() const { return config_; }

    /// Compute the unit tangent at a converged point.
    ///
    /// `directionSign` selects the initial travel direction (sign of lambda_dot)
    /// when no previous tangent is available. When `previous` is non-null the sign
    /// is instead chosen so the new tangent has a non-negative inner product with
    /// it, keeping the traversal moving consistently along the branch.
    ArclengthTangent computeTangent(const ArclengthState& point,
                                    Real directionSign,
                                    const ArclengthTangent* previous = nullptr) const
    {
        const VectorXd fLambda = system_.parameterDerivative(point.x, point.lambda);
        VectorXd z(point.x.size());
        if (!system_.solveJacobian(point.x, point.lambda, fLambda, z)) {
            throw std::runtime_error(
                "PseudoArclengthContinuation: Jacobian solve failed while computing tangent.");
        }
        const Real theta2 = config_.parameterScale * config_.parameterScale;
        const Real stateWeight = stateWeightForSize(point.x.size());
        const Real denom = std::sqrt(stateWeight * z.squaredNorm() + theta2);
        if (!(denom > 0.0) || !std::isfinite(denom)) {
            throw std::runtime_error(
                "PseudoArclengthContinuation: degenerate tangent normalization.");
        }
        Real sign = directionSign >= 0.0 ? 1.0 : -1.0;
        ArclengthTangent tangent;
        tangent.lambdaDot = sign / denom;
        tangent.xDot = -tangent.lambdaDot * z;
        if (previous != nullptr) {
            const Real inner =
                stateWeight * tangent.xDot.dot(previous->xDot) +
                theta2 * tangent.lambdaDot * previous->lambdaDot;
            if (inner < 0.0) {
                tangent.xDot = -tangent.xDot;
                tangent.lambdaDot = -tangent.lambdaDot;
            }
        }
        return tangent;
    }

    /// Form a normalized initial tangent from two already-converged branch
    /// points.  The secant is oriented toward directionSign in lambda and uses
    /// the same weighted norm as the analytic tangent and corrector constraint.
    ArclengthTangent computeSecantTangent(const ArclengthState& previous,
                                          const ArclengthState& current,
                                          Real directionSign) const
    {
        if (previous.x.size() != current.x.size() || current.x.size() <= 0) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: secant states must have equal, non-zero size.");
        }
        ArclengthTangent tangent;
        tangent.xDot = current.x - previous.x;
        tangent.lambdaDot = current.lambda - previous.lambda;
        if (!tangent.xDot.allFinite() || !std::isfinite(tangent.lambdaDot)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: secant difference must be finite.");
        }
        const Real theta2 = config_.parameterScale * config_.parameterScale;
        const Real stateWeight = stateWeightForSize(current.x.size());
        const Real norm = std::sqrt(
            stateWeight * tangent.xDot.squaredNorm() +
            theta2 * tangent.lambdaDot * tangent.lambdaDot);
        if (!(norm > 0.0) || !std::isfinite(norm)) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: secant points must be distinct.");
        }
        tangent.xDot /= norm;
        tangent.lambdaDot /= norm;
        const Real sign = directionSign >= 0.0 ? 1.0 : -1.0;
        if (sign * tangent.lambdaDot < 0.0) {
            tangent.xDot = -tangent.xDot;
            tangent.lambdaDot = -tangent.lambdaDot;
        }
        return tangent;
    }

    /// Attempt a single arclength step from `anchor` along `tangent`.
    ///
    /// The corrector keeps the predictor tangent fixed (standard pseudo-arclength).
    /// On corrector failure the step length is shrunk by `shrinkFactor` and retried
    /// up to `maxStepRetries` times or until it falls below `minStep`.
    ArclengthStepResult step(const ArclengthState& anchor,
                             const ArclengthTangent& tangent,
                             Real arclengthStep) const
    {
        ArclengthStepResult result;
        const Real theta2 = config_.parameterScale * config_.parameterScale;
        const Real stateWeight = stateWeightForSize(anchor.x.size());
        Real deltaS = arclengthStep;
        for (int retry = 0; retry <= config_.maxStepRetries; ++retry) {
            result.retries = retry;
            if (!(deltaS > 0.0) || deltaS < config_.minStep) {
                constexpr const char* minStepReason =
                    "arclength step shrank below min_step";
                if (result.failureReason.empty()) {
                    result.failureReason = minStepReason;
                } else {
                    result.failureReason += "; ";
                    result.failureReason += minStepReason;
                }
                result.arclengthStep = deltaS;
                return result;
            }
            // A retry is a fresh corrector attempt.  Retain the previous cause
            // only when the next shrink crosses minStep (handled above), not
            // when a later retry fails for a different reason.
            result.failureReason.clear();

            // Tangent predictor.
            ArclengthState current;
            current.x = anchor.x + deltaS * tangent.xDot;
            current.lambda = anchor.lambda + deltaS * tangent.lambdaDot;

            bool correctorOk = false;
            int iter = 0;
            Real residualNorm = std::numeric_limits<Real>::infinity();
            for (iter = 0; iter < config_.maxCorrectorIterations; ++iter) {
                VectorXd f;
                try {
                    f = system_.residual(current.x, current.lambda);
                } catch (const std::exception& error) {
                    result.failureReason =
                        std::string("residual evaluation failed: ") + error.what();
                    break;
                }
                const Real arclengthResidual = arclengthResidualNorm(
                    anchor, tangent, current, deltaS, stateWeight, theta2);
                const Real equationResidualNorm = infinityNorm(f);
                residualNorm = std::max(equationResidualNorm, std::abs(arclengthResidual));
                result.equationResidualNorm = equationResidualNorm;
                result.arclengthResidual = arclengthResidual;
                result.bestTrialResidualNorm = std::numeric_limits<Real>::infinity();
                result.bestTrialEquationResidualNorm =
                    std::numeric_limits<Real>::infinity();
                result.bestTrialArclengthResidual =
                    std::numeric_limits<Real>::infinity();
                result.bestTrialAlpha = 0.0;
                result.smallestTrialAlpha = 0.0;
                result.lineSearchAttempts = 0;
                result.rawStateUpdateNorm = 0.0;
                result.limitedStateUpdateNorm = 0.0;
                result.updateLimitChangeNorm = 0.0;
                result.rawParameterUpdate = 0.0;
                result.limitedParameterUpdate = 0.0;
                if (!std::isfinite(residualNorm)) {
                    result.failureReason = "non-finite corrector residual";
                    break;
                }
                if (residualNorm <= config_.correctorTolerance) {
                    correctorOk = true;
                    break;
                }

                VectorXd fLambda;
                try {
                    fLambda = system_.parameterDerivative(current.x, current.lambda);
                } catch (const std::exception& error) {
                    result.failureReason =
                        std::string("parameter derivative failed: ") + error.what();
                    break;
                }
                VectorXd a(current.x.size());
                VectorXd zStep(current.x.size());
                const VectorXd negF = -f;
                const bool pairSolved = system_.solveJacobianPair
                    ? system_.solveJacobianPair(
                          current.x,
                          current.lambda,
                          negF,
                          fLambda,
                          a,
                          zStep)
                    : (system_.solveJacobian(
                           current.x, current.lambda, negF, a) &&
                       system_.solveJacobian(
                           current.x, current.lambda, fLambda, zStep));
                if (!pairSolved) {
                    result.failureReason =
                        "Jacobian solve failed for bordered corrector RHS pair";
                    break;
                }

                const Real denom =
                    theta2 * tangent.lambdaDot - stateWeight * tangent.xDot.dot(zStep);
                if (!std::isfinite(denom) || std::abs(denom) <
                        std::numeric_limits<Real>::min()) {
                    result.failureReason = "singular bordered corrector denominator";
                    break;
                }
                Real deltaLambda =
                    (-arclengthResidual - stateWeight * tangent.xDot.dot(a)) / denom;
                if (!std::isfinite(deltaLambda)) {
                    result.failureReason = "non-finite bordered corrector update";
                    break;
                }
                if (config_.maxParameterUpdate > 0.0) {
                    deltaLambda = std::clamp(
                        deltaLambda,
                        -config_.maxParameterUpdate,
                        config_.maxParameterUpdate);
                }
                // Keep the state and parameter components on the same
                // bordered-Newton solve after applying the lambda trust
                // region.  Reusing deltaX from the uncapped deltaLambda would
                // violate J*deltaX + F_lambda*deltaLambda = -F.
                VectorXd deltaX = a - zStep * deltaLambda;
                if (!deltaX.allFinite()) {
                    result.failureReason = "non-finite bordered corrector update";
                    break;
                }
                const VectorXd rawDeltaX = deltaX;
                const Real rawDeltaLambda = deltaLambda;
                result.rawStateUpdateNorm = infinityNorm(rawDeltaX);
                result.rawParameterUpdate = rawDeltaLambda;
                if (system_.limitUpdate) {
                    system_.limitUpdate(current.x, deltaX, deltaLambda);
                    if (!deltaX.allFinite() || !std::isfinite(deltaLambda)) {
                        result.failureReason = "non-finite limited corrector update";
                        break;
                    }
                }
                result.limitedStateUpdateNorm = infinityNorm(deltaX);
                result.updateLimitChangeNorm = infinityNorm(deltaX - rawDeltaX);
                result.limitedParameterUpdate = deltaLambda;
                bool accepted = false;
                Real alpha = config_.dampingFactor;
                for (int ls = 0; ls <= config_.maxLineSearchSteps; ++ls) {
                    result.lineSearchAttempts = ls + 1;
                    result.smallestTrialAlpha = alpha;
                    ArclengthState trial;
                    trial.x = current.x + alpha * deltaX;
                    trial.lambda = current.lambda + alpha * deltaLambda;
                    VectorXd trialF;
                    try {
                        trialF = system_.residual(trial.x, trial.lambda);
                    } catch (const std::exception&) {
                        alpha *= 0.5;
                        continue;
                    }
                    const Real trialArclengthResidual = arclengthResidualNorm(
                        anchor, tangent, trial, deltaS, stateWeight, theta2);
                    const Real trialEquationResidualNorm = infinityNorm(trialF);
                    const Real trialNorm = std::max(
                        trialEquationResidualNorm, std::abs(trialArclengthResidual));
                    if (std::isfinite(trialNorm) &&
                        trialNorm < result.bestTrialResidualNorm) {
                        result.bestTrialResidualNorm = trialNorm;
                        result.bestTrialEquationResidualNorm =
                            trialEquationResidualNorm;
                        result.bestTrialArclengthResidual = trialArclengthResidual;
                        result.bestTrialAlpha = alpha;
                    }
                    const Real acceptedResidual = residualNorm *
                        (1.0 + config_.lineSearchRelativeIncreaseTolerance);
                    if (std::isfinite(trialNorm) && trialNorm <= acceptedResidual) {
                        current = std::move(trial);
                        accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if (!accepted) {
                    result.failureReason = "corrector line search rejected update";
                    break;
                }
            }

            if (correctorOk) {
                result.converged = true;
                result.state = current;
                result.arclengthStep = deltaS;
                result.correctorIterations = iter;
                result.residualNorm = residualNorm;
                return result;
            }
            result.correctorIterations = iter;
            result.residualNorm = residualNorm;
            if (result.failureReason.empty()) {
                result.failureReason =
                    "corrector failed to converge within iteration budget";
            }
            deltaS *= config_.shrinkFactor;
        }

        if (result.failureReason.empty())
            result.failureReason = "corrector failed to converge within retry budget";
        result.arclengthStep = deltaS;
        return result;
    }

    /// Suggest the next arclength step length after a step result, clamped to bounds.
    Real nextStep(const ArclengthStepResult& result) const
    {
        Real next = result.converged
            ? result.arclengthStep * config_.growthFactor
            : result.arclengthStep * config_.shrinkFactor;
        if (next > config_.maxStep)
            next = config_.maxStep;
        if (next < config_.minStep)
            next = config_.minStep;
        return next;
    }

private:
    Real stateWeightForSize(Eigen::Index size) const
    {
        if (size <= 0) {
            throw std::invalid_argument(
                "PseudoArclengthContinuation: state vector must be non-empty.");
        }
        if (config_.stateWeight > 0.0)
            return config_.stateWeight;
        return 1.0 / static_cast<Real>(size);
    }

    static Real arclengthResidualNorm(const ArclengthState& anchor,
                                      const ArclengthTangent& tangent,
                                      const ArclengthState& current,
                                      Real deltaS,
                                      Real stateWeight,
                                      Real theta2)
    {
        return stateWeight * tangent.xDot.dot(current.x - anchor.x) +
               theta2 * tangent.lambdaDot * (current.lambda - anchor.lambda) -
               deltaS;
    }

    static Real infinityNorm(const VectorXd& v)
    {
        Real norm = 0.0;
        for (Index i = 0; i < static_cast<Index>(v.size()); ++i) {
            const Real a = std::abs(v(static_cast<Eigen::Index>(i)));
            if (!std::isfinite(a))
                return std::numeric_limits<Real>::infinity();
            norm = std::max(norm, a);
        }
        return norm;
    }

    ArclengthSystem system_;
    PseudoArclengthConfig config_;
};

} // namespace vela
