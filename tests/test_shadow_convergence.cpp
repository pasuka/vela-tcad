#include <catch2/catch_test_macros.hpp>
#include "vela/solver/ShadowConvergenceDiagnostics.h"
using namespace vela;
using namespace vela::shadow_convergence;

TEST_CASE("Shadow states do not turn absent or uncalibrated evidence into success", "[shadow_convergence]") {
    REQUIRE(combine({}) == "not_evaluated");
    REQUIRE(combine({"pass", "not_evaluated"}) == "not_evaluated");
    REQUIRE(combine({"not_calibrated", "not_evaluated"}) == "not_calibrated");
    REQUIRE(combine({"not_calibrated", "fail"}) == "fail");
    REQUIRE(combine({"pass", "pass"}) == "pass");
    REQUIRE(uncalibrated(0.)["status"] == "not_calibrated");
}
TEST_CASE("Compensated shadow defect retains a small term between cancelling terms", "[shadow_convergence]") {
    SparseMatrixd J(1, 3);
    J.insert(0, 0) = 1e30; J.insert(0, 1) = 1.; J.insert(0, 2) = -1e30;
    VectorXd d = VectorXd::Ones(3), r = VectorXd::Zero(1);
    const auto q = linearQuality(J, r, d);
    REQUIRE(q["beta_comp"]["value"].get<double>() > 4.99e-31);
    REQUIRE(q["beta_comp"]["value"].get<double>() < 5.01e-31);
    REQUIRE(q["worst_row"] == 0);
    REQUIRE_FALSE(q["forward_error_guaranteed"].get<bool>());
}
TEST_CASE("Zero system is not a calibrated quality certificate", "[shadow_convergence]") {
    SparseMatrixd J(2, 2);
    const auto q = linearQuality(J, VectorXd::Zero(2), VectorXd::Zero(2));
    REQUIRE(q["beta_norm"]["status"] == "not_evaluated");
    REQUIRE(q["beta_comp"]["status"] == "not_calibrated");
    REQUIRE(q["exact_zero_denominator_rows"] == 2);
}
