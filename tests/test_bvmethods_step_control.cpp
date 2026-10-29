#include "../scripts/diagnostics/BvMethodsAdaptiveStep.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

using nlohmann::json;
using namespace vela::diagnostics;

namespace {
json resistorConfig()
{
    return {{"sweep", {{"start", 0.}, {"stop", 100000.},
        {"initial_step", .01}, {"min_step", 1e-5}, {"max_step", 2500.},
        {"external_circuit", {{"mode", "series_resistor"}, {"resistance_ohm_um", 1e7}}}}}};
}
}

TEST_CASE("BV predictor reads explicit device volts without borrowing parser defaults", "[bvmethods_step_control]")
{
    auto config = resistorConfig();
    const auto fallback = resolveInnerStepLimit(config, json::object());
    REQUIRE(fallback.volts == Catch::Approx(.02));
    REQUIRE(fallback.source == "fallback");
    REQUIRE(fallback.reason.find("external-voltage") != std::string::npos);

    config["sweep"]["external_circuit"]["max_inner_voltage_step_V"] = .075;
    const auto scripted = resolveInnerStepLimit(config, json::object());
    REQUIRE(scripted.volts == Catch::Approx(.075));
    REQUIRE(scripted.source == "sweep.external_circuit.max_inner_voltage_step_V");
    const auto override = resolveInnerStepLimit(config, {{"inner_step_control", {{"max_step_V", .004}}}});
    REQUIRE(override.volts == Catch::Approx(.004));
    REQUIRE(override.source == "manifest.inner_step_control.max_step_V");
}

TEST_CASE("BV normalized inner steps scale with the actual device span", "[bvmethods_step_control]")
{
    const auto config = resistorConfig();
    for (double span : {.5, 20., 1000.}) {
        for (double origin : {-100., 0., 100.}) {
            for (double sign : {-1., 1.}) {
                const json manifest = {{"inner_step_control", {{"max_step_fraction", .003},
                    {"initial_voltage_V", origin}, {"goal_voltage_V", origin + sign * span}}}};
                const auto result = resolveInnerStepLimit(config, manifest);
                REQUIRE(result.volts == Catch::Approx(.003 * span));
                REQUIRE(result.source.find("fraction") != std::string::npos);
            }
        }
    }
    const auto noInnerSpan = resolveInnerStepLimit(config,
        {{"inner_step_control", {{"max_step_fraction", .025}}}});
    REQUIRE(noInnerSpan.volts == Catch::Approx(.02));
    REQUIRE(noInnerSpan.source == "fallback");
    REQUIRE(noInnerSpan.reason.find("device-voltage span") != std::string::npos);
    const auto zeroSpan = resolveInnerStepLimit(config,
        {{"inner_step_control", {{"max_step_fraction", .025},
          {"initial_voltage_V", 4.}, {"goal_voltage_V", 4.}}}});
    REQUIRE(zeroSpan.source == "fallback");
}

TEST_CASE("BV device-voltage sweep steps retain their converted volt units", "[bvmethods_step_control]")
{
    // The frontend has already converted a normalized MaxStep to physical V.
    const json config = {{"sweep", {{"start", 0.}, {"stop", 10.}, {"max_step", .05}}}};
    const auto direct = resolveInnerStepLimit(config, json::object());
    REQUIRE(direct.volts == Catch::Approx(.05));
    REQUIRE(direct.source == "sweep.max_step");
    const auto fraction = resolveInnerStepLimit(config,
        {{"inner_step_control", {{"max_step_fraction", .003}}}});
    REQUIRE(fraction.volts == Catch::Approx(.03));
    auto resistor = config;
    resistor["sweep"]["external_circuit"] = {{"resistance_ohm_um", 1e7}};
    REQUIRE(resolveInnerStepLimit(resistor, json::object()).source == "fallback");
}

TEST_CASE("BV bad explicit inner controls fail rather than selecting a default", "[bvmethods_step_control]")
{
    for (double value : {0., -1., std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
        REQUIRE_THROWS_AS(resolveInnerStepLimit(resistorConfig(),
            {{"inner_step_control", {{"max_step_V", value}}}}), std::invalid_argument);
        auto config = resistorConfig();
        config["sweep"]["external_circuit"]["max_inner_voltage_step_V"] = value;
        REQUIRE_THROWS_AS(resolveInnerStepLimit(config, json::object()), std::invalid_argument);
    }
    REQUIRE_THROWS_AS(resolveInnerStepLimit(resistorConfig(),
        {{"inner_step_control", {{"max_step_fraction", 1.01}}}}), std::invalid_argument);
    REQUIRE_THROWS_AS(resolveInnerStepLimit(resistorConfig(),
        {{"inner_step_control", {{"max_step_V", true}}}}), std::invalid_argument);
    REQUIRE_THROWS_AS(resolveInnerStepLimit(resistorConfig(),
        {{"inner_step_control", {{"max_stpe_V", .1}}}}), std::invalid_argument);
}

TEST_CASE("BV secant guard permits large source steps and reversed device motion", "[bvmethods_step_control]")
{
    OuterStepControl outer;
    outer.minimum = 1e-12;
    for (double slope : {1., .1, 1e-4, -1e-4}) {
        const auto step = guardedOuterStep(1000., outer, .02, slope, 1.);
        REQUIRE(step.has_value());
        REQUIRE(std::abs(*step * slope) <= .02 * (1. + 1e-14));
        REQUIRE(*step <= 1000.);
    }
    REQUIRE(*guardedOuterStep(1000., outer, .02, 1e-4, 1.) == Catch::Approx(200.));
    REQUIRE(*guardedOuterStep(1000., outer, .02, 0., 1.) == Catch::Approx(1000.));
    // A minimum source step must not silently override a smaller inner limit.
    outer.minimum = .01;
    REQUIRE_FALSE(guardedOuterStep(1., outer, .001, 1., 1.).has_value());
    REQUIRE_THROWS_AS(guardedOuterStep(1., outer, .02, .1, 0.), std::invalid_argument);
}

TEST_CASE("BV outer controls use script values and preserve restart overrides", "[bvmethods_step_control]")
{
    auto config = resistorConfig();
    config["sweep"].update({{"growth_factor", 1.6}, {"shrink_factor", .4}, {"step_growth_mode", "fixed"}});
    config["sweep"]["external_circuit"].update({{"coupled_min_outer_step_V", 1e-4},
        {"coupled_max_outer_step_V", 10.}, {"coupled_outer_growth_factor", 1.2},
        {"coupled_max_step_retries", 7}});
    const auto control = resolveOuterStepControl(config, {{"initial_step_V", .03}});
    REQUIRE(control.initial == Catch::Approx(.03));
    REQUIRE(control.minimum == Catch::Approx(1e-4));
    REQUIRE(control.maximum == Catch::Approx(10.));
    REQUIRE(control.growth == Catch::Approx(1.2));
    REQUIRE(control.shrink == Catch::Approx(.4));
    REQUIRE(control.growthMode == "fixed");
    REQUIRE(control.maxRetries == 7);
    config["sweep"]["external_circuit"]["coupled_max_outer_step_V"] = 0.;
    REQUIRE(resolveOuterStepControl(config, json::object()).maximum == Catch::Approx(2500.));
    REQUIRE_THROWS_AS(resolveOuterStepControl(config, {{"initial_step_V", 1e-6}}), std::invalid_argument);
}

TEST_CASE("BV report-only predictor preserves outer advances even when its shadow guard rejects", "[bvmethods_step_control]")
{
    OuterStepControl control;control.minimum=.01;control.maximum=2500.;
    const auto enforced=predictorStepDecision(1.,control,.02,1.,1.,false);
    const auto report=predictorStepDecision(1.,control,.02,1.,1.,true);
    REQUIRE(*enforced.applied==Catch::Approx(.02));
    REQUIRE(report.guarded==enforced.guarded);
    REQUIRE(*report.applied==1.);
    const auto impossible=predictorStepDecision(1.,control,.001,1.,1.,false);
    const auto shadow=predictorStepDecision(1.,control,.001,1.,1.,true);
    REQUIRE_FALSE(impossible.applied.has_value());
    REQUIRE_FALSE(shadow.guarded.has_value());
    REQUIRE(*shadow.applied==1.);
    REQUIRE(*predictorStepDecision(4000.,control,.02,-1e-6,1.,true).applied==2500.);
}
