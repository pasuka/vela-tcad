#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "vela/material/MaterialDatabase.h"
#include "vela/core/PhysicalConstants.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <fstream>
#include <filesystem>

using namespace vela;

TEST_CASE("Built-in materials match the emitted no-BGN 300 K snapshot", "[material][builtin_alignment]")
{
    std::ifstream input(VELA_DEFAULT_MATERIALS_FILE);
    REQUIRE(input.good());
    const auto snapshot = nlohmann::json::parse(input);
    REQUIRE(snapshot.at("schema") == "vela.templates_ldmos.materials.v1");
    REQUIRE(snapshot.at("materials").size() == 2);
    const auto mode = GENERATE(UnitScalingMode::LegacySI, UnitScalingMode::UnitScaling);
    const UnitScalingConfig scaling(mode);
    const auto& units = scaling.unitSystem();
    const MaterialDatabase builtin(scaling);
    const MaterialDatabase loaded(VELA_DEFAULT_MATERIALS_FILE, scaling);
    for (const auto& entry : snapshot.at("materials")) {
        const auto name = entry.at("name").get<std::string>();
        CAPTURE(name, scaling.isUnitScaling());
        const auto& actual = builtin.getMaterial(name);
        const auto& imported = loaded.getMaterial(name);
        const auto equal = [](Real value, Real expected) {
            REQUIRE(value == Catch::Approx(expected).epsilon(1e-12));
        };
        equal(actual.eps_r, entry.at("eps_r").get<Real>());
        REQUIRE(actual.bandgap_eV.has_value());
        REQUIRE(actual.electron_affinity_eV.has_value());
        REQUIRE(actual.temperature_K.has_value());
        equal(*actual.bandgap_eV, entry.at("bandgap_eV").get<Real>());
        equal(*actual.electron_affinity_eV, entry.at("electron_affinity_eV").get<Real>());
        equal(*actual.temperature_K, entry.at("temperature_K").get<Real>());
        equal(imported.eps_r, actual.eps_r);
        equal(*imported.bandgap_eV, *actual.bandgap_eV);
        equal(*imported.electron_affinity_eV, *actual.electron_affinity_eV);
        if (name == "Si") {
            REQUIRE(actual.Nc_m3.has_value());
            REQUIRE(actual.Nv_m3.has_value());
            equal(units.internalConcentrationToM3(actual.ni),
                  entry.at("intrinsic_carrier_density_cm3").get<Real>() * 1e6);
            equal(units.internalConcentrationToM3(*actual.Nc_m3),
                  entry.at("conduction_band_density_of_states_cm3").get<Real>() * 1e6);
            equal(units.internalConcentrationToM3(*actual.Nv_m3),
                  entry.at("valence_band_density_of_states_cm3").get<Real>() * 1e6);
            equal(units.internalMobilityToM2PerVS(actual.mun),
                  entry.at("electron_mobility_cm2_per_V_s").get<Real>() * 1e-4);
            equal(units.internalMobilityToM2PerVS(actual.mup),
                  entry.at("hole_mobility_cm2_per_V_s").get<Real>() * 1e-4);
            equal(imported.ni, actual.ni);
            equal(imported.mun, actual.mun);
            equal(imported.mup, actual.mup);
            equal(*imported.Nc_m3, *actual.Nc_m3);
            equal(*imported.Nv_m3, *actual.Nv_m3);
            REQUIRE_FALSE(actual.electron_quantum_gamma.has_value());
        }
    }
}

TEST_CASE("Built-in silicon satisfies Boltzmann mass action across its temperature law", "[material][builtin_alignment]")
{
    const auto mode = GENERATE(UnitScalingMode::LegacySI, UnitScalingMode::UnitScaling);
    const auto temperature = GENERATE(250.0, 300.0, 350.0);
    const MaterialDatabase db{UnitScalingConfig(mode)};
    const auto si = db.getMaterial("Si", temperature);
    REQUIRE(si.Nc_m3.has_value());
    REQUIRE(si.Nv_m3.has_value());
    REQUIRE(si.bandgap_eV.has_value());
    const Real thermal_eV = constants::kb * temperature / constants::q;
    // Logarithms avoid squaring large densities and test the physical relation.
    const Real inferred_gap = thermal_eV * (std::log(*si.Nc_m3 / si.ni)
                                         + std::log(*si.Nv_m3 / si.ni));
    REQUIRE(inferred_gap == Catch::Approx(*si.bandgap_eV).epsilon(1e-12));
}

TEST_CASE("Frozen partial material baselines preserve their pre-alignment effective silicon", "[material][builtin_alignment]")
{
    const auto relative = GENERATE(
        "reference_tcad/bvmethods_sentaurus2018/vela/materials_sentaurus2018.json",
        "reference_tcad/pn2d_sentaurus2018/source/pn2d_sentaurus2018_iv_materials.json",
        "reference_tcad/schottky_charon_sentaurus2018/vela/materials_mu_n_sentaurus2018.json",
        "reference_tcad/schottky_charon_sentaurus2018/vela/materials_mu_n_nc_eg_sentaurus2018.json",
        "reference_tcad/singledevice_sentaurus2018/vela/materials_sentaurus2018.json",
        "tests/fixtures/pn2d_bv_rejected_transition/materials.json");
    CAPTURE(relative);
    const std::filesystem::path root(VELA_SOURCE_DIR);
    const auto original = root / relative;
    const auto frozen = original.parent_path() / (original.stem().string() + "_legacy_300k.json");
    const UnitScalingConfig scaling(UnitScalingMode::UnitScaling);
    MaterialDatabase old(root.string() +
        "/reference_tcad/schottky_charon_sentaurus2018/vela/materials_legacy_builtin_300k.json", scaling);
    const auto& legacy = old.getMaterial("Si");
    REQUIRE(legacy.ni == Catch::Approx(1e10).epsilon(1e-12));
    REQUIRE(legacy.mun == Catch::Approx(1350.0).epsilon(1e-12));
    REQUIRE(legacy.mup == Catch::Approx(480.0).epsilon(1e-12));
    REQUIRE(*legacy.Nc_m3 == Catch::Approx(2.8e19).epsilon(1e-12));
    REQUIRE(*legacy.Nv_m3 == Catch::Approx(1.04e19).epsilon(1e-12));
    REQUIRE(*legacy.bandgap_eV == Catch::Approx(1.12).epsilon(1e-12));
    old.loadJson(original.string(), scaling);
    const MaterialDatabase completed(frozen.string(), scaling);
    const auto& expected = old.getMaterial("Si");
    const auto& actual = completed.getMaterial("Si");
    REQUIRE(actual.ni == Catch::Approx(expected.ni).epsilon(1e-12));
    REQUIRE(actual.mun == Catch::Approx(expected.mun).epsilon(1e-12));
    REQUIRE(actual.mup == Catch::Approx(expected.mup).epsilon(1e-12));
    REQUIRE(*actual.Nc_m3 == Catch::Approx(*expected.Nc_m3).epsilon(1e-12));
    REQUIRE(*actual.Nv_m3 == Catch::Approx(*expected.Nv_m3).epsilon(1e-12));
    REQUIRE(actual.bandgap_eV == expected.bandgap_eV);
    REQUIRE(actual.electron_affinity_eV == expected.electron_affinity_eV);
    REQUIRE(actual.electron_quantum_gamma == expected.electron_quantum_gamma);
    REQUIRE(actual.electron_quantum_dos_mass_ratio == expected.electron_quantum_dos_mass_ratio);
    REQUIRE(actual.electron_quantum_coefficient_mass_ratio == expected.electron_quantum_coefficient_mass_ratio);
}
