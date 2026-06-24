#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/simulation/PseudoArclength.h"

#include <cmath>

using namespace vela;
using Catch::Approx;

namespace {

// Unit circle F(x, lambda) = x^2 + lambda^2 - 1 has a turning point (fold) at
// (x, lambda) = (0, 1) where d(lambda)/d(arclength) -> 0. A voltage-parameterized
// continuation that advances lambda by a fixed step cannot pass this fold because
// the upper-left arc requires lambda to decrease while x becomes negative.
ArclengthSystem makeCircleSystem()
{
    ArclengthSystem system;
    system.residual = [](const VectorXd& x, Real lambda) {
        VectorXd f(1);
        f(0) = x(0) * x(0) + lambda * lambda - 1.0;
        return f;
    };
    system.parameterDerivative = [](const VectorXd&, Real lambda) {
        VectorXd d(1);
        d(0) = 2.0 * lambda;
        return d;
    };
    system.solveJacobian = [](const VectorXd& x, Real, const VectorXd& b, VectorXd& y) {
        const Real jac = 2.0 * x(0);
        if (std::abs(jac) < 1.0e-14)
            return false;
        y.resize(1);
        y(0) = b(0) / jac;
        return true;
    };
    return system;
}

PseudoArclengthConfig makeCircleConfig()
{
    PseudoArclengthConfig config;
    config.enabled = true;
    config.initialStep = 0.2;
    config.minStep = 1.0e-4;
    config.maxStep = 0.25;
    config.growthFactor = 1.1;
    config.shrinkFactor = 0.5;
    config.maxCorrectorIterations = 40;
    config.correctorTolerance = 1.0e-11;
    config.maxStepRetries = 12;
    config.parameterScale = 1.0;
    return config;
}

} // namespace

TEST_CASE("PseudoArclength: continuation crosses a turning point", "[arclength]")
{
    PseudoArclengthContinuation continuation(makeCircleSystem(), makeCircleConfig());

    ArclengthState point;
    point.x = VectorXd::Constant(1, 1.0);
    point.lambda = 0.0;

    ArclengthTangent tangent = continuation.computeTangent(point, +1.0);
    // Initial tangent at (1, 0): branch advances lambda with x stationary.
    REQUIRE(tangent.lambdaDot == Approx(1.0).margin(1.0e-9));
    REQUIRE(tangent.xDot(0) == Approx(0.0).margin(1.0e-9));

    Real step = continuation.config().initialStep;
    Real maxLambda = point.lambda;
    bool crossedFold = false;
    for (int i = 0; i < 200 && !crossedFold; ++i) {
        ArclengthStepResult result = continuation.step(point, tangent, step);
        REQUIRE(result.converged);
        // Each accepted point must lie on the unit circle.
        REQUIRE(result.state.x(0) * result.state.x(0) +
                    result.state.lambda * result.state.lambda ==
                Approx(1.0).margin(1.0e-8));
        point = result.state;
        maxLambda = std::max(maxLambda, point.lambda);
        ArclengthTangent previous = tangent;
        tangent = continuation.computeTangent(
            point, previous.lambdaDot >= 0.0 ? 1.0 : -1.0, &previous);
        step = continuation.nextStep(result);
        if (point.x(0) < -0.1 && point.lambda > 0.5)
            crossedFold = true;
    }

    // The fold at lambda = 1 was reached and crossed onto the x < 0 branch.
    REQUIRE(crossedFold);
    REQUIRE(maxLambda == Approx(1.0).margin(0.05));
    REQUIRE(maxLambda <= 1.0 + 1.0e-6);
}

TEST_CASE("PseudoArclength: step fails when length falls below min_step",
          "[arclength]")
{
    PseudoArclengthConfig config = makeCircleConfig();
    config.minStep = 1.0;       // larger than the requested step length
    config.maxStepRetries = 0;  // no retries
    PseudoArclengthContinuation continuation(makeCircleSystem(), config);

    ArclengthState point;
    point.x = VectorXd::Constant(1, 1.0);
    point.lambda = 0.0;
    ArclengthTangent tangent = continuation.computeTangent(point, +1.0);

    ArclengthStepResult result = continuation.step(point, tangent, 0.2);
    REQUIRE_FALSE(result.converged);
    REQUIRE(result.failureReason == "arclength step shrank below min_step");
}

TEST_CASE("PseudoArclength: missing callbacks are rejected", "[arclength]")
{
    ArclengthSystem incomplete;
    incomplete.residual = [](const VectorXd& x, Real) { return x; };
    // parameterDerivative and solveJacobian intentionally left empty.
    REQUIRE_THROWS_AS(
        PseudoArclengthContinuation(incomplete, makeCircleConfig()),
        std::invalid_argument);
}

TEST_CASE("PseudoArclength: parameterScale must be positive", "[arclength]")
{
    PseudoArclengthConfig config = makeCircleConfig();
    config.parameterScale = 0.0;
    REQUIRE_THROWS_AS(
        PseudoArclengthContinuation(makeCircleSystem(), config),
        std::invalid_argument);
}

TEST_CASE("PseudoArclength: singular tangent solve is reported", "[arclength]")
{
    ArclengthSystem singular = makeCircleSystem();
    singular.solveJacobian = [](const VectorXd&, Real, const VectorXd&, VectorXd&) {
        return false;
    };
    PseudoArclengthContinuation continuation(singular, makeCircleConfig());

    ArclengthState point;
    point.x = VectorXd::Constant(1, 1.0);
    point.lambda = 0.0;
    REQUIRE_THROWS_AS(continuation.computeTangent(point, +1.0), std::runtime_error);
}
