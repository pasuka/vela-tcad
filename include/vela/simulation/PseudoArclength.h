#pragma once

#include "vela/core/Types.h"

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

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
    /// Weight (theta) of the parameter component in the arclength norm
    /// ||x_dot||^2 + theta^2 * lambda_dot^2 = 1. Must be > 0.
    Real parameterScale = 1.0;
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
};

/// Pseudo-arclength (Keller) continuation engine.
///
/// The engine implements a tangent predictor followed by a bordered-Newton
/// corrector that augments the system with one arclength constraint equation
///   N(x, lambda) = x_dot . (x - x0) + theta^2 * lambda_dot * (lambda - lambda0) - Delta s.
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
        const Real denom = std::sqrt(z.squaredNorm() + theta2);
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
                tangent.xDot.dot(previous->xDot) + theta2 * tangent.lambdaDot * previous->lambdaDot;
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
                const Real arclengthResidual =
                    tangent.xDot.dot(current.x - anchor.x) +
                    theta2 * tangent.lambdaDot * (current.lambda - anchor.lambda) - deltaS;
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

                const Real denom = theta2 * tangent.lambdaDot - tangent.xDot.dot(zStep);
                if (!std::isfinite(denom) || std::abs(denom) <
                        std::numeric_limits<Real>::min()) {
                    break;
                }
                const Real deltaLambda =
                    (-arclengthResidual - tangent.xDot.dot(a)) / denom;
                const VectorXd deltaX = a - zStep * deltaLambda;
                if (!deltaX.allFinite() || !std::isfinite(deltaLambda))
                    break;
                current.x += deltaX;
                current.lambda += deltaLambda;
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
