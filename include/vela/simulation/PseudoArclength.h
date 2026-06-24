#pragma once

#include "vela/core/Types.h"

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

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
};

/// A point on the solution branch: state vector x and continuation parameter lambda.
struct ArclengthState {
    VectorXd x;
    Real lambda = 0.0;
};

/// Unit tangent direction in (x, lambda) space (already scaled by parameterScale).
struct ArclengthTangent {
    VectorXd xDot;
    Real lambdaDot = 0.0;
};

/// Outcome of a single attempted arclength step.
struct ArclengthStepResult {
    bool converged = false;
    ArclengthState state;
    Real arclengthStep = 0.0;   ///< Delta s actually accepted.
    int correctorIterations = 0;
    int retries = 0;
    Real residualNorm = 0.0;    ///< max(||F||_inf, |N|) at the accepted point.
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
    /// Solve J(x, lambda) * y = b where J = dF/dx. Writes the solution into y and
    /// returns false if the Jacobian is singular or the solve otherwise failed.
    std::function<bool(const VectorXd& x, Real lambda, const VectorXd& b, VectorXd& y)>
        solveJacobian;
    /// Optional in-place limiter for the bordered-Newton update before line search.
    std::function<void(const VectorXd& x, VectorXd& deltaX, Real& deltaLambda)>
        limitUpdate;
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
                result.failureReason = "arclength step shrank below min_step";
                result.arclengthStep = deltaS;
                return result;
            }

            // Tangent predictor.
            ArclengthState current;
            current.x = anchor.x + deltaS * tangent.xDot;
            current.lambda = anchor.lambda + deltaS * tangent.lambdaDot;

            bool correctorOk = false;
            int iter = 0;
            Real residualNorm = std::numeric_limits<Real>::infinity();
            for (iter = 0; iter < config_.maxCorrectorIterations; ++iter) {
                const VectorXd f = system_.residual(current.x, current.lambda);
                const Real arclengthResidual = arclengthResidualNorm(
                    anchor, tangent, current, deltaS, stateWeight, theta2);
                residualNorm = std::max(infinityNorm(f), std::abs(arclengthResidual));
                if (!std::isfinite(residualNorm))
                    break;
                if (residualNorm <= config_.correctorTolerance) {
                    correctorOk = true;
                    break;
                }

                const VectorXd fLambda =
                    system_.parameterDerivative(current.x, current.lambda);
                VectorXd a(current.x.size());
                VectorXd zStep(current.x.size());
                const VectorXd negF = -f;
                if (!system_.solveJacobian(current.x, current.lambda, negF, a))
                    break;
                if (!system_.solveJacobian(current.x, current.lambda, fLambda, zStep))
                    break;

                const Real denom =
                    theta2 * tangent.lambdaDot - stateWeight * tangent.xDot.dot(zStep);
                if (!std::isfinite(denom) || std::abs(denom) <
                        std::numeric_limits<Real>::min()) {
                    break;
                }
                Real deltaLambda =
                    (-arclengthResidual - stateWeight * tangent.xDot.dot(a)) / denom;
                VectorXd deltaX = a - zStep * deltaLambda;
                if (!deltaX.allFinite() || !std::isfinite(deltaLambda))
                    break;
                if (system_.limitUpdate) {
                    system_.limitUpdate(current.x, deltaX, deltaLambda);
                    if (!deltaX.allFinite() || !std::isfinite(deltaLambda))
                        break;
                }

                bool accepted = false;
                Real alpha = config_.dampingFactor;
                for (int ls = 0; ls <= config_.maxLineSearchSteps; ++ls) {
                    ArclengthState trial;
                    trial.x = current.x + alpha * deltaX;
                    trial.lambda = current.lambda + alpha * deltaLambda;
                    const VectorXd trialF = system_.residual(trial.x, trial.lambda);
                    const Real trialArclengthResidual = arclengthResidualNorm(
                        anchor, tangent, trial, deltaS, stateWeight, theta2);
                    const Real trialNorm = std::max(
                        infinityNorm(trialF), std::abs(trialArclengthResidual));
                    if (std::isfinite(trialNorm) && trialNorm < residualNorm) {
                        current = std::move(trial);
                        accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if (!accepted)
                    break;
            }

            if (correctorOk) {
                result.converged = true;
                result.state = current;
                result.arclengthStep = deltaS;
                result.correctorIterations = iter;
                result.residualNorm = residualNorm;
                return result;
            }
            deltaS *= config_.shrinkFactor;
        }

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
