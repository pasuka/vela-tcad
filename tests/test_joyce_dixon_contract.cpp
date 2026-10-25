#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/experimental/JoyceDixonStatistics.h"
#include <cmath>
#include <limits>

using vela::experimental::JoyceDixonStatistics;
using Catch::Approx;

TEST_CASE("Experimental Joyce-Dixon contract is monotone and invertible",
          "[carrier_statistics][joyce_dixon]")
{
    double previous = 0.;
    for (double eta = -80.; eta <= 80.; eta += .25) {
        const auto jet = JoyceDixonStatistics::evaluate(eta);
        REQUIRE(jet.value > previous);
        REQUIRE(jet.first > 0.);
        REQUIRE(JoyceDixonStatistics::inverse(jet.value)
                == Approx(eta).margin(8e-13));
        previous = jet.value;
    }
    const auto dilute = JoyceDixonStatistics::evaluate(-30.);
    REQUIRE(dilute.value == Approx(std::exp(-30.)).epsilon(1e-12));
    const auto dense = JoyceDixonStatistics::evaluate(1e6);
    REQUIRE(dense.value == Approx(4./(3.*std::sqrt(std::acos(-1.)))
                                  * std::pow(1e6,1.5)).epsilon(1e-7));
    REQUIRE_THROWS_AS(JoyceDixonStatistics::inverse(0.), std::invalid_argument);
    REQUIRE_THROWS_AS(JoyceDixonStatistics::inverse(-1.), std::invalid_argument);
    REQUIRE(std::isinf(JoyceDixonStatistics::inverse(
        std::numeric_limits<double>::infinity())));
}

TEST_CASE("Experimental Joyce-Dixon eta derivatives follow its own population",
          "[carrier_statistics][joyce_dixon]")
{
    for (const double eta : {-30., -8., -1., 0., 2., 4., 5., 8., 30.}) {
        const double h = 1e-4;
        const auto mid = JoyceDixonStatistics::evaluate(eta);
        const auto lo = JoyceDixonStatistics::evaluate(eta-h);
        const auto hi = JoyceDixonStatistics::evaluate(eta+h);
        REQUIRE((hi.value-lo.value)/(2*h)
                == Approx(mid.first).epsilon(3e-8));
        REQUIRE((hi.first-lo.first)/(2*h)
                == Approx(mid.second).epsilon(3e-8));
    }
    const auto join = JoyceDixonStatistics::join();
    const auto lo = JoyceDixonStatistics::evaluate(join.eta-1e-8);
    const auto hi = JoyceDixonStatistics::evaluate(join.eta+1e-8);
    REQUIRE(lo.value == Approx(hi.value).epsilon(1e-8));
    REQUIRE(lo.first == Approx(hi.first).epsilon(1e-8));
}

TEST_CASE("Experimental statistics and FD-BGN recover measured reference states",
          "[carrier_statistics][joyce_dixon]")
{
    // Independent energy/population observations; source and scope are recorded
    // in the BVMethods statistics-contract validation report. Not a global fit.
    REQUIRE(JoyceDixonStatistics::inverse(4.189432330541544)
            == Approx(2.836420040347931).margin(5e-13));
    REQUIRE(JoyceDixonStatistics::inverse(9.061482309140672)
            == Approx(5.082184251016262).margin(5e-13));
    const double vt = .025851995266484913;
    const auto correction = [](double ratio) {
        return std::max(0., JoyceDixonStatistics::inverse(ratio)-std::log(ratio));
    };
    const double fd = vt*(correction(2.5999999999999938e20/2.8582520073866568e19)
                         +correction(9.999999999999986e17/3.104628932474839e19));
    REQUIRE(fd == Approx(.07496092103813459).margin(1e-13));
}

#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/BandgapNarrowing.h"
TEST_CASE("Joyce-Dixon public density and FD-BGN dispatch use one inverse", "[bvmethods_profile]")
{
    using namespace vela;
    const CarrierStatisticsConfig cfg{"fermi_dirac_joyce_dixon"};
    const auto model=carrierStatisticsModel(cfg);
    const double vt=.02585,ni=1e16,nc=2.8e25,nv=3.1e25;
    for(double eta:{-30.,-2.,1.,5.,8.,20.}) {
        const double psi=vt*(eta-std::log(ni/nc)),fn=.003;
        const double n=electronDensity(ni,nc,psi+fn,fn,vt,cfg);
        REQUIRE(inverseFermiDiracHalf(n/nc,model)==Approx(eta).margin(1e-12));
        REQUIRE(electronQuasiFermiPotential(ni,nc,psi+fn,n,vt,cfg)==Approx(fn).margin(1e-13));
        const double h=1e-7;
        const double fd=(electronDensity(ni,nc,psi+fn+h,fn,vt,cfg)-electronDensity(ni,nc,psi+fn-h,fn,vt,cfg))/(2*h);
        REQUIRE(electronDensityDerivativeEta(ni,nc,psi+fn,fn,vt,cfg)/vt==Approx(fd).epsilon(1e-8));
        const double p=holeDensity(ni,nv,-psi+fn,fn,vt,cfg);
        REQUIRE(holeQuasiFermiPotential(ni,nv,-psi+fn,p,vt,cfg)==Approx(fn).margin(1e-13));
    }
    const double nd=2e26,na=3e25;
    const double expected=vt*(inverseFermiDiracHalf(nd/nc,model)-std::log(nd/nc)+inverseFermiDiracHalf(na/nv,model)-std::log(na/nv));
    REQUIRE(fermiStatisticsBandgapCorrection(nd,na,nc,nv,vt,model)==Approx(expected).epsilon(1e-14));
    REQUIRE(fermiStatisticsBandgapCorrection(nd,na,nc,nv,vt)!=Approx(expected).epsilon(1e-6));
    REQUIRE(fermiDiracHalf(3.,CarrierStatisticsModel::FermiDirac)==fermiDiracHalf(3.));
    REQUIRE(fermiDiracHalfDerivative(3.,CarrierStatisticsModel::FermiDirac)==fermiDiracHalfDerivative(3.));
}
