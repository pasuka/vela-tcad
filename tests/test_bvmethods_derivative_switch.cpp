#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../scripts/diagnostics/BvMethodsDerivativeSwitch.h"
using namespace vela::diagnostics;

TEST_CASE("Derivative switch retains feasible restoration through large update errors", "[derivative_switch]") {
    DerivativeSwitch control;
    REQUIRE(control.requestFull(true));
    // No raw-update comparison exists in this API. Both physical states pass
    // their original gates even if their update RMS differ substantially.
    REQUIRE(control.fullCandidate({true,true,1e-12,2e-12,0.}));
    REQUIRE(control.mode()==DerivativeMode::Trial);
    REQUIRE(control.requestFull(true));
    REQUIRE(control.fullCandidate({true,true,2e-12,1e-15,0.}));
    REQUIRE(control.mode()==DerivativeMode::Full);
}
TEST_CASE("Derivative switch requires persistent evidence before rollback", "[derivative_switch]") {
    DerivativeSwitch control({},true);
    REQUIRE(control.fullCandidate({true,false,1.,1.2,0.}));
    REQUIRE(control.mode()==DerivativeMode::Full);
    REQUIRE_FALSE(control.fullCandidate({true,false,1.2,1.4,0.}));
    REQUIRE(control.mode()==DerivativeMode::Bootstrap);
    REQUIRE(control.reason=="confirmed_model_failure");
}
TEST_CASE("Derivative switch cooldown prevents immediate retry of invalid candidates", "[derivative_switch]") {
    DerivativeSwitch control;
    REQUIRE(control.requestFull(true));
    REQUIRE_FALSE(control.fullCandidate({false,false,1.,0.,0.}));
    control.bootstrapApplied();
    REQUIRE_FALSE(control.requestFull(true));
    control.bootstrapApplied();
    REQUIRE(control.requestFull(true));
}
TEST_CASE("Derivative model ratio uses actual applied step and handles roundoff", "[derivative_switch]") {
    const double alpha=.25,phi=2.;
    const auto q=derivativeQuality({true,false,phi,.8*phi,(1-alpha)*(1-alpha)*phi},.1);
    REQUIRE(q.reliable);
    REQUIRE(q.predicted==Catch::Approx((2*alpha-alpha*alpha)*phi));
    REQUIRE(q.ratio==Catch::Approx(.2/(2*alpha-alpha*alpha)));
    const auto floor=derivativeQuality({true,true,1e-28,1e-28,1e-28},.1);
    REQUIRE_FALSE(floor.reliable);REQUIRE(floor.adequate);
    REQUIRE_FALSE(derivativeQuality({true,false,1.,.5,2.},.1).adequate);
}
TEST_CASE("Derivative switch bounds trial windows and probes stalled bootstrap", "[derivative_switch]") {
    DerivativeSwitchConfig config;config.maximumTrialSteps=2;
    DerivativeSwitch control(config);
    REQUIRE_FALSE(control.requestFull(false));control.bootstrapApplied();
    REQUIRE_FALSE(control.requestFull(false));control.bootstrapApplied();
    REQUIRE(control.requestFull(false));
    REQUIRE(control.fullCandidate({true,false,1.,.5,0.}));
    REQUIRE_FALSE(control.fullCandidate({true,false,.5,.6,0.}));
    REQUIRE(control.reason=="restore_window_exhausted");
    config.minimumTrialSteps=3;
    REQUIRE_THROWS_AS(DerivativeSwitch(config),std::invalid_argument);
}

TEST_CASE("Linear failure restoration is bounded until a state update", "[derivative_switch]") {
    DerivativeSwitch control;
    REQUIRE(control.probeAfterBootstrapLinearFailure());
    REQUIRE(control.mode()==DerivativeMode::Trial);
    REQUIRE(control.reason=="bootstrap_linear_failure_restore");
    REQUIRE_FALSE(control.probeAfterBootstrapLinearFailure());
    REQUIRE_FALSE(control.fullCandidate({}));
    REQUIRE(control.mode()==DerivativeMode::Bootstrap);
    REQUIRE_FALSE(control.probeAfterBootstrapLinearFailure());
    control.stateAdvanced();
    REQUIRE(control.probeAfterBootstrapLinearFailure());
    DerivativeSwitch full({},true);
    REQUIRE_FALSE(full.probeAfterBootstrapLinearFailure());
}

TEST_CASE("Model and linear failure probes share one unchanged-state allowance", "[derivative_switch]") {
    DerivativeSwitch c;
    REQUIRE(c.probeAfterBootstrapModelFailure());
    REQUIRE(c.reason=="bootstrap_model_non_decrease_restore");
    REQUIRE_FALSE(c.probeAfterBootstrapLinearFailure());
    REQUIRE_FALSE(c.fullCandidate({}));
    REQUIRE_FALSE(c.probeAfterBootstrapModelFailure());
    REQUIRE_FALSE(c.probeAfterBootstrapLinearFailure());
    c.stateAdvanced();
    REQUIRE(c.probeAfterBootstrapLinearFailure());
    REQUIRE_FALSE(c.fullCandidate({}));
    REQUIRE_FALSE(c.probeAfterBootstrapModelFailure());
}
