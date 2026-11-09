#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/material/MaterialDatabase.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <cmath>
using namespace vela;

TEST_CASE("Constants contract reproduces thermal and Poisson reference coefficients", "[constants_profile]")
{
    const bool native = constants::profileName == "sentaurus_u2022_12";
    const double vt = native ? 0.025851995266484913 : 0.025851999786435535;
    CHECK(constants::Vt_300 == Catch::Approx(vt).epsilon(2e-15));
    CHECK(constants::q / constants::eps0 == Catch::Approx(
        native ? 1.8095276817781392e-8 : 1.8095128179727827e-8).epsilon(2e-15));
    // Material mass action must use the same thermal voltage as the core.
    const MaterialDatabase db;
    const auto& si = db.getMaterial("Si");
    const double inferredGap = constants::Vt_300 *
        (std::log(*si.Nc_m3 / si.ni) + std::log(*si.Nv_m3 / si.ni));
    CHECK(inferredGap == Catch::Approx(*si.bandgap_eV).epsilon(2e-15));
    const double prefactor = 2 * std::pow(2 * constants::dosPi * constants::m0 *
        constants::kb * 300 / (constants::h * constants::h), 1.5) / 1e6;
    CHECK(prefactor == Catch::Approx(native ? 2.5093895119478055e+19 : 2.5094122252280967e19).epsilon(2e-15));
}

TEST_CASE("Mixed deck and material constants contracts fail before assembly", "[constants_profile]")
{
    const std::string selected(constants::profileName);
    const std::string other = selected == "codata2018" ? "sentaurus_u2022_12" : "codata2018";
    CHECK_NOTHROW(parseUnitScalingConfig(nlohmann::json{{"physical_constants_profile", selected}}));
    CHECK_THROWS_AS(parseUnitScalingConfig(nlohmann::json{{"physical_constants_profile", other}}), std::invalid_argument);
    CHECK_THROWS_AS(parseUnitScalingConfig(nlohmann::json{{"physical_constants_profile", "typo"}}), std::invalid_argument);
    const auto path = std::filesystem::current_path() / "constants_profile_material_test.json";
    struct Cleanup { std::filesystem::path p; ~Cleanup(){std::error_code e;std::filesystem::remove(p,e);} } cleanup{path};
    const auto write = [&](const std::string& profile) {
        std::ofstream out(path);
        out << nlohmann::json{{"physical_constants_profile",profile},
            {"materials", nlohmann::json::array({{{"name","Si"},{"eps_r",11.7}}})}};
    };
    MaterialDatabase db;
    write(selected); CHECK_NOTHROW(db.loadJson(path.string()));
    write(other); CHECK_THROWS_AS(db.loadJson(path.string()), std::runtime_error);
    { std::ofstream out(path);out << R"({"materials":[{"name":"Si","eps_r":11.7}]})"; }
    if (selected == "codata2018") CHECK_NOTHROW(db.loadJson(path.string()));
    else CHECK_THROWS_AS(db.loadJson(path.string()), std::runtime_error);
}
