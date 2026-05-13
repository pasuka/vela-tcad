#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/numerics/LineSearch.h"

#include <vector>

using namespace vela;

TEST_CASE("BacktrackingLineSearch: backtracks until caller accepts candidate", "[line_search]")
{
    LineSearchConfig cfg;
    cfg.initialDamping = 1.0;
    cfg.minDamping = 0.1;
    cfg.reduction = 0.5;
    cfg.sufficientDecrease = 0.0;
    cfg.maxBacktracks = 4;

    BacktrackingLineSearch search(cfg);
    VectorXd x(1);
    x << 0.0;
    VectorXd step(1);
    step << 1.0;
    VectorXd currentResidual(1);
    currentResidual << 1.0;

    std::vector<Real> candidates;
    const LineSearchResult result = search.search(
        x,
        step,
        currentResidual,
        [&](const VectorXd& candidate) {
            candidates.push_back(candidate(0));
            VectorXd residual(1);
            residual << 0.25;
            return residual;
        },
        [](const VectorXd& candidate, const VectorXd&) {
            return candidate(0) <= 0.25;
        });

    REQUIRE(result.accepted);
    REQUIRE(result.damping == Catch::Approx(0.25));
    REQUIRE(result.x(0) == Catch::Approx(0.25));
    REQUIRE(result.residualNorm == Catch::Approx(0.25));
    REQUIRE(candidates == std::vector<Real>{1.0, 0.5, 0.25});
}

TEST_CASE("BacktrackingLineSearch: rejection restores original state", "[line_search]")
{
    LineSearchConfig cfg;
    cfg.initialDamping = 1.0;
    cfg.minDamping = 0.25;
    cfg.reduction = 0.5;
    cfg.maxBacktracks = 10;

    BacktrackingLineSearch search(cfg);
    VectorXd x(1);
    x << -2.0;
    VectorXd step(1);
    step << 4.0;
    VectorXd currentResidual(1);
    currentResidual << 3.0;

    int residualCalls = 0;
    const LineSearchResult result = search.search(
        x,
        step,
        currentResidual,
        [&](const VectorXd&) {
            ++residualCalls;
            VectorXd residual(1);
            residual << 1.0;
            return residual;
        },
        [](const VectorXd&, const VectorXd&) {
            return false;
        });

    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.damping == Catch::Approx(0.0));
    REQUIRE(result.residualNorm == Catch::Approx(currentResidual.norm()));
    REQUIRE((result.x - x).norm() == Catch::Approx(0.0));
    REQUIRE((result.residual - currentResidual).norm() == Catch::Approx(0.0));
    REQUIRE(residualCalls == 3);
}

TEST_CASE("BacktrackingLineSearch: min damping bounds failed backtracking", "[line_search]")
{
    LineSearchConfig cfg;
    cfg.initialDamping = 1.0;
    cfg.minDamping = 0.25;
    cfg.reduction = 0.5;
    cfg.maxBacktracks = 8;

    BacktrackingLineSearch search(cfg);
    VectorXd x(1);
    x << 0.0;
    VectorXd step(1);
    step << 1.0;
    VectorXd currentResidual(1);
    currentResidual << 1.0;

    int residualCalls = 0;
    const LineSearchResult result = search.search(
        x,
        step,
        currentResidual,
        [&](const VectorXd&) {
            ++residualCalls;
            VectorXd residual(1);
            residual << 2.0;
            return residual;
        });

    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.damping == Catch::Approx(0.0));
    REQUIRE(result.x(0) == Catch::Approx(0.0));
    REQUIRE(residualCalls == 3);
}
