#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/RecombinationModel.h"

#include <cmath>
#include <stdexcept>
#include <limits>

using namespace vela;

TEST_CASE("SRH recombination is near zero when n*p equals ni squared", "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"srh"}, 1.0e-7, 2.0e-7));
    const Real ni = 1.0e16;
    const Real n = 2.0e21;
    const Real p = ni * ni / n;

    REQUIRE(model.srhRate(n, p, ni) == Catch::Approx(0.0).margin(1.0e6));
    REQUIRE(model.totalRate(n, p, ni) == Catch::Approx(0.0).margin(1.0e6));
}

TEST_CASE("Auger recombination increases at high carrier concentration", "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"auger"}));
    const Real ni = 1.0e16;

    const Real low = model.augerRate(1.0e21, 1.0e21, ni);
    const Real high = model.augerRate(1.0e24, 1.0e24, ni);

    REQUIRE(low > 0.0);
    REQUIRE(high > low);
}

TEST_CASE("Auger linearization remains finite for extreme initializer carriers",
          "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"auger"}));
    const Real ni = 1.0e16;
    const Real n = 1.0e234;
    const Real p = 1.0e234;

    const RecombinationLinearization electron = model.electronLinearization(n, p, ni);
    const RecombinationLinearization hole = model.holeLinearization(n, p, ni);

    REQUIRE(std::isfinite(electron.diagonal));
    REQUIRE(std::isfinite(electron.rhs));
    REQUIRE(std::isfinite(hole.diagonal));
    REQUIRE(std::isfinite(hole.rhs));
}

TEST_CASE("Total recombination is SRH plus Auger", "[recombination]")
{
    RecombinationModel total(recombinationModelConfig({"srh", "auger"}));
    const Real n = 1.0e22;
    const Real p = 2.0e22;
    const Real ni = 1.0e16;

    REQUIRE(total.totalRate(n, p, ni) ==
            Catch::Approx(total.srhRate(n, p, ni) + total.augerRate(n, p, ni)));
}

TEST_CASE("Default bandgap narrowing interface returns zero", "[bgn]")
{
    NoBandgapNarrowing bgn;
    REQUIRE(bgn.deltaEg(1.0e25, 1.0e24, 1.0e20) == Catch::Approx(0.0));
}

TEST_CASE("Slotboom bandgap narrowing grows effective intrinsic density", "[bgn]")
{
    BandgapNarrowingConfig cfg;
    cfg.model = "slotboom";
    SlotboomBandgapNarrowing bgn(cfg);

    const Real low = bgn.deltaEg(1.0e21, 0.0, 0.0);
    const Real high = bgn.deltaEg(1.0e25, 0.0, 0.0);

    REQUIRE(low >= 0.0);
    REQUIRE(high > low);
    REQUIRE(high == Catch::Approx(0.0833788).epsilon(1.0e-5));

    const Real ni = 1.0e16;
    const Real Vt = 0.025852;
    const Real niEff = effectiveIntrinsicDensity(ni, Vt, high);
    REQUIRE(niEff > ni);
    REQUIRE(niEff == Catch::Approx(ni * std::exp(high / (2.0 * Vt))));
}

TEST_CASE("Effective intrinsic density caps overflow", "[bgn]")
{
    const Real ni = 1.0e16;
    const Real Vt = 0.025852;
    const Real niEff = effectiveIntrinsicDensity(ni, Vt, 1.0e6);

    REQUIRE(std::isfinite(niEff));
    REQUIRE(niEff == std::numeric_limits<Real>::max());
}

TEST_CASE("Bandgap narrowing factory validates model names", "[bgn]")
{
    REQUIRE(makeBandgapNarrowingModel(bandgapNarrowingConfig("none"))->deltaEg(1.0e25, 0.0, 0.0)
            == Catch::Approx(0.0));
    REQUIRE(makeBandgapNarrowingModel(bandgapNarrowingConfig("slotboom"))->deltaEg(1.0e25, 0.0, 0.0)
            > 0.0);
    REQUIRE_THROWS_AS(makeBandgapNarrowingModel(bandgapNarrowingConfig("unknown")),
                      std::invalid_argument);
}

TEST_CASE("CarrierStatistics intrinsic density uses temperature_K material path", "[temperature]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    const Real ni300 = intrinsicDensity(si, 300.0);
    const Real ni450 = intrinsicDensity(si, 450.0);
    REQUIRE(ni300 == Catch::Approx(si.ni));
    REQUIRE(ni450 > ni300);
}
