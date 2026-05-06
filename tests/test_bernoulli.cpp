#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/discretization/Bernoulli.h"
#include <cmath>

using namespace vela;
using Catch::Approx;

TEST_CASE("Bernoulli: B(0) ≈ 1", "[bernoulli]")
{
    REQUIRE(bernoulli(0.0) == Approx(1.0).epsilon(1.0e-12));
}

TEST_CASE("Bernoulli: tiny x uses Taylor", "[bernoulli]")
{
    // |x| < 1e-10: B(x) = 1 - x/2 + x²/12
    const double x = 1.0e-12;
    const double expected = 1.0 - x * 0.5 + x * x / 12.0;
    REQUIRE(bernoulli(x) == Approx(expected).epsilon(1.0e-12));

    const double xn = -1.0e-12;
    const double expectedn = 1.0 - xn * 0.5 + xn * xn / 12.0;
    REQUIRE(bernoulli(xn) == Approx(expectedn).epsilon(1.0e-12));
}

TEST_CASE("Bernoulli: B(x) + B(-x) = x * (exp(x)+1)/(exp(x)-1)", "[bernoulli]")
{
    // Identity: B(x) + B(-x) = x + x*exp(x)/(exp(x)-1) = x*(1+exp(x))/(exp(x)-1)
    // A simpler derived identity: B(x) - B(-x) = -x  (verify B(-x) = B(x) + x)
    //   B(x) = x/(exp(x)-1)
    //   B(-x) = -x/(exp(-x)-1) = x*exp(x)/(exp(x)-1)
    //   B(-x) - B(x) = x  →  B(x) - B(-x) = -x
    for (double x : {-5.0, -1.0, -0.1, 0.1, 1.0, 5.0}) {
        const double diff = bernoulli(x) - bernoulli(-x);
        REQUIRE(diff == Approx(-x).epsilon(1.0e-10));
    }
}

TEST_CASE("Bernoulli: non-negative and finite for a range of inputs", "[bernoulli]")
{
    // B(x) is mathematically positive for all finite x.
    // For very large x (>~700), x*exp(-x) underflows to 0.0 in double
    // precision – this is correct (no overflow, NaN, or negative result).
    for (double x : {-1000.0, -500.0, -100.0, -10.0, -1.0, 0.0,
                      1.0, 10.0, 100.0, 500.0, 1000.0}) {
        const double val = bernoulli(x);
        REQUIRE(std::isfinite(val));
        REQUIRE(val >= 0.0);
    }
    // Values in the well-behaved range must be strictly positive
    for (double x : {-100.0, -10.0, -1.0, 0.0, 1.0, 10.0, 100.0}) {
        REQUIRE(bernoulli(x) > 0.0);
    }
}

TEST_CASE("Bernoulli: large positive x does not overflow or underflow", "[bernoulli]")
{
    // B(x) ≈ x * exp(-x) for x > 500 → should be tiny but positive
    const double val = bernoulli(600.0);
    REQUIRE(val > 0.0);
    REQUIRE(std::isfinite(val));
}

TEST_CASE("Bernoulli: large negative x does not overflow", "[bernoulli]")
{
    // B(x) ≈ -x for x < -500
    const double val = bernoulli(-600.0);
    REQUIRE(val == Approx(600.0).epsilon(1.0e-6));
    REQUIRE(std::isfinite(val));
}
