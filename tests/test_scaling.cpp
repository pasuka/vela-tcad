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

TEST_CASE("UnitScalingConfig unit_scaling keeps TCAD values internal", "[scaling]")
{
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};

    REQUIRE(unitScaling.lengthToInternal(0.01) == Catch::Approx(0.01));
    REQUIRE(unitScaling.concentrationToInternal(1.0e17) == Catch::Approx(1.0e17));
    REQUIRE(unitScaling.sheetDensityToInternal(2.5e11) == Catch::Approx(2.5e11));
    REQUIRE(unitScaling.mobilityToInternal(1000.0) == Catch::Approx(1000.0));
    REQUIRE(unitScaling.velocityToInternal(1.07e7) == Catch::Approx(1.07e7));
    REQUIRE(unitScaling.electricFieldToInternal(3.0e5) == Catch::Approx(3.0e5));
    REQUIRE(unitScaling.inverseLengthToInternal(1.2e6) == Catch::Approx(1.2e6));
    REQUIRE(unitScaling.surfaceFieldCoefficientToInternal(0.7) == Catch::Approx(0.7));

    REQUIRE(unitScaling.unitSystem().internalLengthToMeters(0.01) ==
            Catch::Approx(1.0e-8));
    REQUIRE(unitScaling.unitSystem().internalConcentrationToM3(1.0e17) ==
            Catch::Approx(1.0e23));
    REQUIRE(unitScaling.velocityToSI(1.07e7) == Catch::Approx(1.07e5));
    REQUIRE(unitScaling.unitSystem().internalElectricFieldToVPerM(3.0e5) ==
            Catch::Approx(3.0e7));
}

TEST_CASE("PhysicalUnitSystem exposes TCAD composite factors", "[scaling]")
{
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};
    const PhysicalUnitSystem& units = unitScaling.unitSystem();

    REQUIRE(units.chargeAreaFactor() == Catch::Approx(1.0e-6));
    REQUIRE(units.chargeLineFactor() == Catch::Approx(1.0e-2));
    REQUIRE(units.chargeVolumeFactor() == Catch::Approx(1.0e-12));
    REQUIRE(units.chargeSheetFactor() == Catch::Approx(1.0e-8));
    REQUIRE(units.fieldFromCoordinateDeltaFactor() == Catch::Approx(1.0e4));
    REQUIRE(units.continuitySourceIntegralFactor() == Catch::Approx(1.0e-4));
    REQUIRE(units.currentPerInternalDepthFactor() == Catch::Approx(1.0e-6));
}
TEST_CASE("PhysicalUnitSystem 2D charge integrals match SI representation", "[scaling][poisson][charge]")
{
    const PhysicalUnitSystem units = PhysicalUnitSystem::tcadInternal();

    const Real legacyVolumeChargePerDepth = 1.0e23 * 1.0e-12;
    const Real tcadVolumeChargePerDepth =
        1.0e17 * 1.0 * units.chargeAreaFactor();
    REQUIRE(tcadVolumeChargePerDepth ==
            Catch::Approx(legacyVolumeChargePerDepth).epsilon(1e-12));

    const Real legacySheetChargePerDepth = 1.0e15 * 1.0e-6;
    const Real tcadSheetChargePerDepth =
        1.0e11 * 1.0 * units.chargeLineFactor();
    REQUIRE(tcadSheetChargePerDepth ==
            Catch::Approx(legacySheetChargePerDepth).epsilon(1e-12));
}

TEST_CASE("UnitScalingSystem: unit_scaling references stay in TCAD internal units", "[scaling]")
{
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};

    UnitScalingSystem::AutoInputs tcadInputs;
    tcadInputs.maxAbsNetDoping_m3 = 1.0e17;
    tcadInputs.niFloor_m3 = 1.0e10;
    tcadInputs.meshMaxLength_m = 1.0;
    tcadInputs.maxMobility_m2_V_s = 1000.0;

    const UnitScalingSystem tcad = UnitScalingSystem::fromInputs(
        300.0,
        constants::eps0 * 11.7,
        tcadInputs,
        UnitScalingReferenceConfig{},
        unitScaling.unitSystem());

    REQUIRE(tcad.C0() == Catch::Approx(1.0e17).epsilon(1e-12));
    REQUIRE(tcad.L0() == Catch::Approx(1.0).epsilon(1e-12));
    REQUIRE(tcad.mu0() == Catch::Approx(1000.0).epsilon(1e-12));

    UnitScalingSystem::AutoInputs legacyInputs;
    legacyInputs.maxAbsNetDoping_m3 = 1.0e23;
    legacyInputs.niFloor_m3 = 1.0e16;
    legacyInputs.meshMaxLength_m = 1.0e-6;
    legacyInputs.maxMobility_m2_V_s = 0.1;

    const UnitScalingSystem legacy = UnitScalingSystem::fromInputs(
        300.0,
        constants::eps0 * 11.7,
        legacyInputs,
        UnitScalingReferenceConfig{},
        UnitScalingConfig{}.unitSystem());

    REQUIRE(tcad.lambda2() == Catch::Approx(legacy.lambda2()).epsilon(1e-12));
    REQUIRE(tcad.E0() == Catch::Approx(legacy.E0() / 100.0).epsilon(1e-12));
    REQUIRE(tcad.J0() == Catch::Approx(legacy.J0() / 1.0e4).epsilon(1e-12));
    REQUIRE(tcad.R0() == Catch::Approx(legacy.R0() / 1.0e6).epsilon(1e-12));
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
    REQUIRE(explicitRefs.characteristicLength_m == Catch::Approx(2.5));
    REQUIRE(explicitRefs.referenceConcentration_m3 == Catch::Approx(5.0e16));
    REQUIRE(explicitRefs.referenceMobility_m2_V_s == Catch::Approx(900.0));
}

// ---------------------------------------------------------------------------
// Deck format_version 2 bridge
// ---------------------------------------------------------------------------

TEST_CASE("Deck format version defaults to 0 when absent", "[scaling][format_version]")
{
    const nlohmann::json cfg = nlohmann::json::object();
    REQUIRE(parseDeckFormatVersion(cfg) == 0);
}

TEST_CASE("Deck format version accepts 2", "[scaling][format_version]")
{
    const nlohmann::json cfg = {{"format_version", 2}};
    REQUIRE(parseDeckFormatVersion(cfg) == 2);
}

TEST_CASE("Deck format version rejects unsupported values", "[scaling][format_version]")
{
    REQUIRE_THROWS_AS(parseDeckFormatVersion(nlohmann::json{{"format_version", 1}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(parseDeckFormatVersion(nlohmann::json{{"format_version", 3}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(parseDeckFormatVersion(nlohmann::json{{"format_version", "2"}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(parseDeckFormatVersion(nlohmann::json{{"format_version", 2.5}}),
                      std::invalid_argument);
}

TEST_CASE("format_version 2 selects the TCAD internal unit system without scaling",
          "[scaling][format_version]")
{
    const nlohmann::json cfg = {{"format_version", 2}};
    const UnitScalingConfig scaling = parseUnitScalingConfig(cfg);
    REQUIRE(scaling.isUnitScaling());
    REQUIRE(scaling.unitSystem().concentrationM3PerInternal() == Catch::Approx(1.0e6));
    REQUIRE(scaling.unitSystem().lengthMPerInternal() == Catch::Approx(1.0e-6));
}

TEST_CASE("format_version 2 is physically equivalent to legacy scaling.mode unit_scaling",
          "[scaling][format_version]")
{
    const UnitScalingConfig viaVersion = parseUnitScalingConfig(
        nlohmann::json{{"format_version", 2}});
    const UnitScalingConfig viaScaling = parseUnitScalingConfig(
        nlohmann::json{{"scaling", {{"mode", "unit_scaling"}}}});
    REQUIRE(viaVersion.isUnitScaling() == viaScaling.isUnitScaling());
    REQUIRE(viaVersion.unitSystem().concentrationM3PerInternal() ==
            Catch::Approx(viaScaling.unitSystem().concentrationM3PerInternal()));
    REQUIRE(viaVersion.unitSystem().lengthMPerInternal() ==
            Catch::Approx(viaScaling.unitSystem().lengthMPerInternal()));
    REQUIRE(viaVersion.unitSystem().chargeAreaFactor() ==
            Catch::Approx(viaScaling.unitSystem().chargeAreaFactor()));
}

TEST_CASE("format_version 2 rejects a scaling block", "[scaling][format_version]")
{
    const nlohmann::json cfg = {
        {"format_version", 2},
        {"scaling", {{"mode", "unit_scaling"}}},
    };
    REQUIRE_THROWS_AS(parseUnitScalingConfig(cfg), std::invalid_argument);
    REQUIRE_THROWS_AS(parseUnitScalingReferenceConfig(cfg), std::invalid_argument);
}

TEST_CASE("format_version 2 reads normalization references from solver.normalization",
          "[scaling][format_version]")
{
    const nlohmann::json autoCfg = {
        {"format_version", 2},
        {"solver", {{"normalization", {
            {"characteristic_length_um", "auto"},
            {"reference_concentration_cm3", "auto"},
            {"reference_mobility_cm2_V_s", "auto"},
        }}}},
    };
    const UnitScalingReferenceConfig autoRefs = parseUnitScalingReferenceConfig(autoCfg);
    REQUIRE(autoRefs.characteristicLength_m == std::nullopt);
    REQUIRE(autoRefs.referenceConcentration_m3 == std::nullopt);
    REQUIRE(autoRefs.referenceMobility_m2_V_s == std::nullopt);

    const nlohmann::json explicitCfg = {
        {"format_version", 2},
        {"solver", {{"normalization", {
            {"characteristic_length_um", 2.5},
            {"reference_concentration_cm3", 5.0e16},
            {"reference_mobility_cm2_V_s", 900.0},
        }}}},
    };
    const UnitScalingReferenceConfig explicitRefs =
        parseUnitScalingReferenceConfig(explicitCfg);
    REQUIRE(explicitRefs.characteristicLength_m == Catch::Approx(2.5));
    REQUIRE(explicitRefs.referenceConcentration_m3 == Catch::Approx(5.0e16));
    REQUIRE(explicitRefs.referenceMobility_m2_V_s == Catch::Approx(900.0));
}

TEST_CASE("solver.normalization matches the legacy scaling reference block",
          "[scaling][format_version]")
{
    const nlohmann::json v2Cfg = {
        {"format_version", 2},
        {"solver", {{"normalization", {
            {"characteristic_length_um", 2.5},
            {"reference_concentration_cm3", 5.0e16},
            {"reference_mobility_cm2_V_s", 900.0},
        }}}},
    };
    const nlohmann::json legacyCfg = {
        {"scaling", {
            {"mode", "unit_scaling"},
            {"characteristic_length_um", 2.5},
            {"reference_concentration_cm3", 5.0e16},
            {"reference_mobility_cm2_V_s", 900.0},
        }},
    };
    const UnitScalingReferenceConfig v2Refs = parseUnitScalingReferenceConfig(v2Cfg);
    const UnitScalingReferenceConfig legacyRefs = parseUnitScalingReferenceConfig(legacyCfg);
    REQUIRE(v2Refs.characteristicLength_m == legacyRefs.characteristicLength_m);
    REQUIRE(v2Refs.referenceConcentration_m3 == legacyRefs.referenceConcentration_m3);
    REQUIRE(v2Refs.referenceMobility_m2_V_s == legacyRefs.referenceMobility_m2_V_s);
}

TEST_CASE("solver.normalization is rejected without format_version 2",
          "[scaling][format_version]")
{
    const nlohmann::json cfg = {
        {"solver", {{"normalization", {{"characteristic_length_um", 2.5}}}}},
    };
    REQUIRE_THROWS_AS(parseUnitScalingReferenceConfig(cfg), std::invalid_argument);
}

TEST_CASE("format_version 2 rejects invalid normalization values",
          "[scaling][format_version]")
{
    const auto make = [](const nlohmann::json& value) {
        return nlohmann::json{
            {"format_version", 2},
            {"solver", {{"normalization", {{"characteristic_length_um", value}}}}},
        };
    };
    REQUIRE_THROWS_AS(parseUnitScalingReferenceConfig(make(0.0)), std::invalid_argument);
    REQUIRE_THROWS_AS(parseUnitScalingReferenceConfig(make(-1.0)), std::invalid_argument);
    REQUIRE_THROWS_AS(parseUnitScalingReferenceConfig(make("default")), std::invalid_argument);
}
