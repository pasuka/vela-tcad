#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/io/SdeScriptReader.h"

using namespace Catch::Matchers;
using namespace vela;

TEST_CASE("sde parser reads rectangle and constant doping placement", "[sde][preprocess]")
{
    const std::string script = R"(
(sdegeo:create-rectangle (position 0 0 0) (position 1.0 2.0 0) "Silicon" "R.Si")
(sdedr:define-constant-profile "CP.N" "PhosphorusActiveConcentration" 1.0e17)
(sdedr:define-constant-profile-region "Place.N" "CP.N" "R.Si")
)";

    SdeScriptReader reader;
    const DeviceIr2D ir = reader.parseText(script, "unit.sde");

    REQUIRE(ir.geometry.size() == 1);
    REQUIRE(ir.regions.size() == 1);
    REQUIRE(ir.dopingProfiles.size() == 1);

    CHECK(ir.regions.front().name == "R.Si");
    CHECK(ir.regions.front().material == "Silicon");

    const auto& profile = ir.dopingProfiles.front();
    CHECK(profile.kind == DopingProfileKind::Constant);
    CHECK(profile.targetRegion == "R.Si");
    CHECK(profile.donors_cm3 == Catch::Approx(1.0e17));
    CHECK(profile.acceptors_cm3 == Catch::Approx(0.0));
}

TEST_CASE("sde parser rejects unsupported phase-1 command", "[sde][preprocess]")
{
    const std::string unsupported = R"(
(sdegeo:create-circle (position 0 0 0) 1.0 "Silicon" "R.Si")
)";

    SdeScriptReader reader;
    REQUIRE_THROWS_WITH(
        reader.parseText(unsupported, "unsupported.sde"),
        ContainsSubstring("unsupported command"));
}

TEST_CASE("sde parser reads gaussian and refinement subset", "[sde][preprocess]")
{
    const std::string script = R"(
(sdegeo:create-rectangle (position 0 0 0) (position 2.0 1.0 0) "Silicon" "R.Si")
(sdedr:define-gaussian-profile "GP.P" "BoronActiveConcentration" "PeakPos" 0.5 "PeakVal" 5.0e18 "ValueAtDepth" 1.0e16 "Depth" 0.2)
(sdedr:define-analytical-profile-region "Place.GP" "GP.P" "R.Si")
(sdedr:define-refinement-size "Ref.Si" 0.05 0.05 0.0 0.0)
(sdedr:define-refinement-region "Place.Ref.Si" "Ref.Si" "R.Si")
(sdedr:define-refinement-function "Ref.Si" "DopingConcentration" "MaxTransDiff" 2.0e16)
)";

    SdeScriptReader reader;
    const DeviceIr2D ir = reader.parseText(script, "gaussian_refine.sde");

    REQUIRE(ir.dopingProfiles.size() == 1);
    const auto& profile = ir.dopingProfiles.front();
    CHECK(profile.kind == DopingProfileKind::Gaussian);
    CHECK(profile.targetRegion == "R.Si");
    CHECK(profile.gaussianPeak_cm3 == Catch::Approx(5.0e18));
    CHECK(profile.gaussianValueAtDepth_cm3 == Catch::Approx(1.0e16));
    CHECK(profile.gaussianPeakPosUm.x_um == Catch::Approx(0.5));
    CHECK(profile.gaussianSigmaXUm > 0.0);
    CHECK(profile.gaussianActsOnDonors == false);

    REQUIRE(ir.meshControl.regionTargetSizeUm.contains("R.Si"));
    CHECK(ir.meshControl.regionTargetSizeUm.at("R.Si") == Catch::Approx(0.05));
    CHECK(ir.meshControl.refineByDopingGradient);
    CHECK(ir.meshControl.dopingGradientThresholdCm3PerUm == Catch::Approx(2.0e16));
}

TEST_CASE("sde parser rejects invalid gaussian subset parameters", "[sde][preprocess]")
{
    const std::string script = R"(
(sdegeo:create-rectangle (position 0 0 0) (position 2.0 1.0 0) "Silicon" "R.Si")
(sdedr:define-gaussian-profile "GP.Bad" "BoronActiveConcentration" "PeakPos" 0.5 "PeakVal" 1.0e16 "ValueAtDepth" 1.0e16 "Depth" 0.2)
(sdedr:define-analytical-profile-region "Place.Bad" "GP.Bad" "R.Si")
)";

    SdeScriptReader reader;
    REQUIRE_THROWS_WITH(
        reader.parseText(script, "invalid_gaussian.sde"),
        ContainsSubstring("ValueAtDepth must be smaller than PeakVal"));
}

TEST_CASE("sde parser allows negative gaussian peak position", "[sde][preprocess]")
{
    const std::string script = R"(
(sdegeo:create-rectangle (position -1 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdedr:define-gaussian-profile "GP.N" "PhosphorusActiveConcentration" "PeakPos" -0.25 "PeakVal" 1.0e19 "ValueAtDepth" 1.0e16 "Depth" 0.2)
(sdedr:define-analytical-profile-region "Place.GP" "GP.N" "R.Si")
)";

    SdeScriptReader reader;
    const DeviceIr2D ir = reader.parseText(script, "negative_peak.sde");
    REQUIRE(ir.dopingProfiles.size() == 1);
    CHECK(ir.dopingProfiles.front().gaussianPeakPosUm.x_um == Catch::Approx(-0.25));
}

TEST_CASE("sde parser reads polygon geometry", "[sde][preprocess]")
{
    const std::string script = R"(
(sdegeo:create-polygon (list (position 0 0 0) (position 1 0 0) (position 1 1 0) (position 0 1 0)) "Silicon" "R.Poly")
)";

    SdeScriptReader reader;
    const DeviceIr2D ir = reader.parseText(script, "polygon.sde");

    REQUIRE(ir.geometry.size() == 1);
    REQUIRE(ir.regions.size() == 1);
    CHECK(ir.geometry.front().kind == GeometryPrimitiveKind::Polygon);
    CHECK(ir.geometry.front().points.size() == 4);
    CHECK(ir.regions.front().name == "R.Poly");
}
