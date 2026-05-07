#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/RecombinationModel.h"

#include <cmath>

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
