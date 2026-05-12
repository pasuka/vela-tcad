#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/core/ScalingSystem.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>

using namespace vela;
using namespace vela::constants;

// Helper: build a default ScalingSystem at 300 K with typical Si parameters
static ScalingSystem makeSi300()
{
    return ScalingSystem(
        300.0,       // temperature [K]
        1.0e22,      // reference concentration [m^-3]
        11.7,        // relative permittivity (Si)
        0.135        // reference mobility [m^2/V/s]
    );
}

TEST_CASE("ScalingSystem: thermal voltage at 300 K", "[scaling]")
{
    ScalingSystem sc = makeSi300();

    // kT/q at 300 K ~= 0.025852 V
    REQUIRE(sc.V0() == Catch::Approx(0.025852).epsilon(1e-3));
}

TEST_CASE("ScalingSystem: Debye length is positive", "[scaling]")
{
    ScalingSystem sc = makeSi300();
    REQUIRE(sc.L0() > 0.0);
}

TEST_CASE("ScalingSystem: scale/unscale potential round-trip", "[scaling]")
{
    ScalingSystem sc = makeSi300();

    const Real phi_physical = 0.7;  // e.g. built-in potential of a p-n junction [V]
    Real phi_scaled   = sc.scalePotential(phi_physical);
    Real phi_restored = sc.unscalePotential(phi_scaled);

    REQUIRE(phi_restored == Catch::Approx(phi_physical).epsilon(1e-12));
}

TEST_CASE("ScalingSystem: scale/unscale concentration round-trip", "[scaling]")
{
    ScalingSystem sc = makeSi300();

    const Real n_physical = 5.0e21;
    Real n_scaled   = sc.scaleConcentration(n_physical);
    Real n_restored = sc.unscaleConcentration(n_scaled);

    REQUIRE(n_restored == Catch::Approx(n_physical).epsilon(1e-12));
}

TEST_CASE("ScalingSystem: J0 and R0 are positive", "[scaling]")
{
    ScalingSystem sc = makeSi300();
    REQUIRE(sc.J0() > 0.0);
    REQUIRE(sc.R0() > 0.0);
}

TEST_CASE("ScalingSystem: invalid inputs throw", "[scaling]")
{
    REQUIRE_THROWS_AS(ScalingSystem(-1.0, 1e22, 11.7, 0.135),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(ScalingSystem(300.0, -1.0, 11.7, 0.135),
                      std::invalid_argument);
}
