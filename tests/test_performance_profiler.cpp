#include <catch2/catch_test_macros.hpp>

#include "vela/core/PerformanceProfiler.h"
#include "vela/physics/CarrierStatistics.h"

#include <chrono>
#include <nlohmann/json.hpp>

using namespace vela;

TEST_CASE("PerformanceProfiler aggregates deterministic stages and counters",
          "[performance_profiler]")
{
    PerformanceProfiler profiler({true, "unused.json"});
    profiler.recordStage("newton.jacobian", std::chrono::nanoseconds(10));
    profiler.recordStage("newton.jacobian", std::chrono::nanoseconds(30));
    profiler.increment("newton.updates", 2);
    profiler.observe("linear.rows", 12.0);
    profiler.observe("linear.rows", 18.0);
    profiler.recordNewtonSolve(true, 2, 1.0, 1.0e-9, "abstol",
                               std::chrono::nanoseconds(50));

    const nlohmann::json json = profiler.toJson();
    REQUIRE(json.at("counters").at("newton.updates") == 4);
    REQUIRE(json.at("counters").at("newton.solve_calls") == 1);
    const auto& stage = json.at("stages").at(0);
    REQUIRE(stage.at("name") == "newton.jacobian");
    REQUIRE(stage.at("calls") == 2);
    REQUIRE(stage.at("total_ns") == 40);
    REQUIRE(stage.at("p50_ns") == 10);
    REQUIRE(stage.at("p95_ns") == 30);
    REQUIRE(json.at("observations").at("linear.rows").at("min") == 12.0);
    REQUIRE(json.at("observations").at("linear.rows").at("max") == 18.0);
    REQUIRE(json.at("newton_solves").at(0).at("iterations") == 2);
}

TEST_CASE("ScopedPerformanceTimer attributes physics calls to its stage",
          "[performance_profiler]")
{
    PerformanceProfiler profiler({true, "unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        {
            ScopedPerformanceTimer outer("stage.outer");
            (void)fermiDiracHalf(0.5);
            {
                ScopedPerformanceTimer inner("stage.inner");
                (void)fermiDiracHalf(0.25);
                (void)fermiDiracHalfDerivative(0.25);
            }
        }
        {
            ScopedPerformanceTimer other("stage.other");
            (void)inverseFermiDiracHalf(1.0);
        }
    }

    const nlohmann::json counters = profiler.toJson().at("counters");
    // Stage counters are inclusive: the outer stage also sees the inner calls.
    REQUIRE(counters.at("stage.inner.fermi_dirac_half_calls") == 1);
    REQUIRE(counters.at("stage.inner.fermi_dirac_half_derivative_calls") == 1);
    REQUIRE(counters.at("stage.outer.fermi_dirac_half_calls") == 2);
    REQUIRE(counters.at("stage.other.inverse_fermi_dirac_half_calls") == 1);
    REQUIRE_FALSE(counters.contains("stage.other.equilibrium_state_solves"));
}

TEST_CASE("ScopedPerformanceTimer reports equilibrium state solve effort",
          "[performance_profiler]")
{
    const CarrierStatisticsConfig fermiDirac{"fermi_dirac"};
    PerformanceProfiler profiler({true, "unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        ScopedPerformanceTimer timer("stage.node_sources");
        const auto state = equilibriumCarrierState(
            1.0e22, 1.0e16, 2.8e25, 1.04e25, 0.02585, fermiDirac);
        REQUIRE(state.n > 0.0);
        REQUIRE(state.p > 0.0);
    }

    const nlohmann::json counters = profiler.toJson().at("counters");
    REQUIRE(counters.at("stage.node_sources.equilibrium_state_solves") == 1);
    REQUIRE(counters.at("stage.node_sources.equilibrium_state_iterations") >= 1);
}

TEST_CASE("Performance profiling config is disabled by default and validates output",
          "[performance_profiler]")
{
    REQUIRE_FALSE(performanceProfilingConfigFromJson(nlohmann::json::object()).enabled);
    const auto config = performanceProfilingConfigFromJson({
        {"performance_profiling", {
            {"enabled", true},
            {"json_file", "profile.json"}
        }}
    });
    REQUIRE(config.enabled);
    REQUIRE(config.jsonFile == "profile.json");
    REQUIRE_THROWS(performanceProfilingConfigFromJson({
        {"performance_profiling", {{"enabled", true}, {"json_file", ""}}}
    }));
}

TEST_CASE("ScopedPerformanceTimer records only inside an active scope",
          "[performance_profiler]")
{
    PerformanceProfiler profiler({true, "unused.json"});
    {
        ScopedPerformanceTimer inactive("inactive");
    }
    {
        ActivePerformanceProfilerScope active(&profiler);
        ScopedPerformanceTimer timer("active");
    }
    const nlohmann::json json = profiler.toJson();
    REQUIRE(json.at("stages").size() == 1);
    REQUIRE(json.at("stages").at(0).at("name") == "active");
}
