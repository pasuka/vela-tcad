#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/experimental/JoyceDixonStatistics.h"
#include <cmath>
#include <limits>
#include <array>
#include <algorithm>

using vela::experimental::JoyceDixonStatistics;
using Catch::Approx;

TEST_CASE("Joyce-Dixon DOS-scaled dilute density keeps a smooth inverse and slope",
          "[carrier_statistics][joyce_dixon]")
{
    const double dos=3.104657033455896e19;
    for(double eta:{-700.,-731.,-750.}){
        const double value=JoyceDixonStatistics::density(eta,dos);
        REQUIRE(std::isnormal(value));
        REQUIRE(JoyceDixonStatistics::inverseDensity(value,dos)==Approx(eta).margin(2e-12));
        const double slope=JoyceDixonStatistics::densityDerivativeEta(eta,dos);
        const double h=1e-3;
        const double fd=(JoyceDixonStatistics::density(eta+h,dos)-JoyceDixonStatistics::density(eta-h,dos))/(2*h);
        REQUIRE(fd/slope==Approx(1.).epsilon(2e-7));
        REQUIRE(JoyceDixonStatistics::densitySlopeFromPopulation(value,dos)/value==Approx(1.));
    }
    REQUIRE(std::exp(-750.)==0.); // The dimensionless intermediate underflows.
    REQUIRE(JoyceDixonStatistics::density(-750.,dos)>0.); // The physical population does not.
}

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
#include "vela/core/PhysicalConstants.h"
#include "vela/physics/BandgapNarrowing.h"
TEST_CASE("Joyce-Dixon public QF mapping remains differentiable beyond normalized underflow", "[bvmethods_profile]")
{
    using namespace vela;
    const CarrierStatisticsConfig cfg{"fermi_dirac_joyce_dixon"};
    const double dos=3.104657033455896e19,ni=1.46e10,vt=constants::Vt_300,psi=.313;
    for(double eta:{-731.,-750.}){
        const double fp=psi+vt*(eta-std::log(ni/dos));
        const double p=holeDensity(ni,dos,psi,fp,vt,cfg);
        REQUIRE(std::isnormal(p));
        REQUIRE(holeQuasiFermiPotential(ni,dos,psi,p,vt,cfg)==Approx(fp).margin(2e-13));
        const double h=1e-6;
        const double fd=(holeDensity(ni,dos,psi,fp+h,vt,cfg)-holeDensity(ni,dos,psi,fp-h,vt,cfg))/(2*h);
        REQUIRE(fd/(holeDensityDerivativeEta(ni,dos,psi,fp,vt,cfg)/vt)==Approx(1.).epsilon(2e-7));
    }
}

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

TEST_CASE("Joyce-Dixon inverse preserves independently reverted Fermi series", "[carrier_statistics][joyce_dixon]") {
    // 80-digit reference: revert sum_k (-1)^(k+1) z^k/k^(3/2), then
    // expand log(z). Dense references use the same mathematical C1 join.
    // This checks coefficient accuracy independently of a double round trip.
    const std::array<std::array<double,2>,6> reference{{{0.5,-0.5175897360762486},{1,0.3487472610045682},{4,2.7296702158669968},{8,4.648909057646425},{9,5.057642331812939},{16,7.555030873204023}}};
    for(const auto& p:reference){
        const double eta=JoyceDixonStatistics::inverse(p[0]);
        REQUIRE(std::abs(eta-p[1])<=4*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(p[1])));
    }
}

TEST_CASE("JD primitive gamma retains the dilute value and population tangent", "[carrier_statistics][joyce_dixon]") {
    using S=JoyceDixonStatistics;
    const double a1=0.353553390593273762200422;
    for(double r:{0.,1e-300,1e-30,1e-18}){
        REQUIRE(std::isfinite(S::primitiveLogGamma(r)));
        REQUIRE(S::primitiveLogGammaSlope(r)==Approx(-a1).epsilon(2e-15));
        if(r>0.)REQUIRE(S::primitiveLogGamma(r)/r==Approx(-a1).epsilon(2e-14));
        else REQUIRE(S::primitiveLogGamma(r)==0.);
    }
    for(double r:{.1,1.,4.,8.,9.,16.}){
        const double h=r*1e-5;
        const double fd=(S::primitiveLogGamma(r+h)-S::primitiveLogGamma(r-h))/(2*h);
        REQUIRE(fd==Approx(S::primitiveLogGammaSlope(r)).epsilon(2e-7));
        REQUIRE(S::primitiveLogGamma(r)==Approx(std::log(r)-S::inverse(r)).epsilon(1e-12).margin(1e-15));
    }
    REQUIRE_THROWS(S::primitiveLogGamma(-1.));
}
