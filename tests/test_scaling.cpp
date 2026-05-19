#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/core/ScalingSystem.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/core/PhysicalConstants.h"
#include <nlohmann/json.hpp>
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

TEST_CASE("UnitScalingSystem: positive reference scales on PN-like deck", "[scaling]")
{
    UnitScalingSystem::AutoInputs inputs;
    inputs.maxAbsNetDoping_m3 = 1.0e23;
    inputs.niFloor_m3 = 1.0e16;
    inputs.meshMaxLength_m = 2.0e-6;
    inputs.maxMobility_m2_V_s = 0.135;

    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        300.0,
        constants::eps0 * 11.7,
        inputs,
        UnitScalingReferenceConfig{});

    REQUIRE(sc.V0() > 0.0);
    REQUIRE(sc.C0() > 0.0);
    REQUIRE(sc.L0() > 0.0);
    REQUIRE(sc.mu0() > 0.0);
    REQUIRE(sc.lambda2() > 0.0);
    REQUIRE(sc.J0() > 0.0);
    REQUIRE(sc.R0() > 0.0);
}

TEST_CASE("UnitScalingSystem: auto C0 uses SI-equivalent concentration from unit_scaling input", "[scaling]")
{
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};

    UnitScalingSystem::AutoInputs inputs;
    inputs.maxAbsNetDoping_m3 = unitScaling.concentrationToSI(1.0e17);
    inputs.niFloor_m3 = unitScaling.concentrationToSI(1.0e10);
    inputs.meshMaxLength_m = unitScaling.lengthToSI(1.0);
    inputs.maxMobility_m2_V_s = unitScaling.mobilityToSI(1000.0);

    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        300.0,
        constants::eps0 * 11.7,
        inputs,
        UnitScalingReferenceConfig{});

    REQUIRE(sc.C0() == Catch::Approx(1.0e23).epsilon(1e-12));
}

TEST_CASE("UnitScalingSystem: scale/unscale round-trip for core quantities", "[scaling]")
{
    UnitScalingSystem::AutoInputs inputs;
    inputs.maxAbsNetDoping_m3 = 5.0e22;
    inputs.niFloor_m3 = 1.0e16;
    inputs.meshMaxLength_m = 1.5e-6;
    inputs.maxMobility_m2_V_s = 0.1;

    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        325.0,
        constants::eps0 * 11.7,
        inputs,
        UnitScalingReferenceConfig{});

    const Real potential = 0.82;
    const Real length = 8.0e-7;
    const Real concentration = 2.5e22;
    const Real field = 4.0e5;
    const Real currentDensity = 3.0e3;

    REQUIRE(sc.unscalePotential(sc.scalePotential(potential))
            == Catch::Approx(potential).epsilon(1e-12));
    REQUIRE(sc.unscaleLength(sc.scaleLength(length))
            == Catch::Approx(length).epsilon(1e-12));
    REQUIRE(sc.unscaleConcentration(sc.scaleConcentration(concentration))
            == Catch::Approx(concentration).epsilon(1e-12));
    REQUIRE(sc.unscaleElectricField(sc.scaleElectricField(field))
            == Catch::Approx(field).epsilon(1e-12));
    REQUIRE(sc.unscaleCurrentDensity(sc.scaleCurrentDensity(currentDensity))
            == Catch::Approx(currentDensity).epsilon(1e-12));
}

TEST_CASE("UnitScalingSystem reference config: supports auto and explicit values", "[scaling]")
{
    const nlohmann::json autoCfg = {
        {"scaling", {
            {"mode", "unit_scaling"},
            {"characteristic_length_um", "auto"},
            {"reference_concentration_cm3", "auto"},
            {"reference_mobility_cm2_V_s", "auto"},
        }}
    };
    const UnitScalingReferenceConfig autoRefs = parseUnitScalingReferenceConfig(autoCfg);
    REQUIRE(autoRefs.characteristicLength_m == std::nullopt);
    REQUIRE(autoRefs.referenceConcentration_m3 == std::nullopt);
    REQUIRE(autoRefs.referenceMobility_m2_V_s == std::nullopt);

    const nlohmann::json explicitCfg = {
        {"scaling", {
            {"mode", "unit_scaling"},
            {"characteristic_length_um", 2.5},
            {"reference_concentration_cm3", 5.0e16},
            {"reference_mobility_cm2_V_s", 900.0},
        }}
    };
    const UnitScalingReferenceConfig explicitRefs =
        parseUnitScalingReferenceConfig(explicitCfg);
    REQUIRE(explicitRefs.characteristicLength_m == Catch::Approx(2.5e-6));
    REQUIRE(explicitRefs.referenceConcentration_m3 == Catch::Approx(5.0e22));
    REQUIRE(explicitRefs.referenceMobility_m2_V_s == Catch::Approx(9.0e-2));
}
