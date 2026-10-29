#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/SlotPopulationNumerics.h"

#include "vela/physics/BandgapNarrowing.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <span>

#include <cmath>
#include <limits>
#include <stdexcept>

using Catch::Approx;
using namespace vela;

namespace finite_step_reference {
#include "../src/physics/FermiHalfCoefficients.inc"
using Wide=boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200>>;
// Independent high-precision Clenshaw value/first derivative of the same
// production approximation. Tests do not compare against another FD model.
std::pair<Wide,Wide> polynomial(Wide eta) {
 std::size_t segment=0;while(eta>FermiBounds[segment+1])++segment;
 Wide scale=Wide(2)/(Wide(FermiBounds[segment+1])-FermiBounds[segment]);
 Wide t=(eta-FermiBounds[segment])*scale-1,b1=0,b2=0,d1=0,d2=0;
 auto c=FermiCoefficients[segment];
 for(std::size_t k=c.size()-1;k>0;--k){
  Wide d=2*b1+2*t*d1-d2,b=Wide(c[k])+2*t*b1-b2;
  d2=d1;d1=d;b2=b1;b1=b;
 }
 return {Wide(c[0])+t*b1-b2,(b1+t*d1-d2)*scale};
}
}

TEST_CASE("Small concentration steps retain nonlinear remainder and tiny QF directions", "[slot_population][finite_update]") {
 using namespace finite_step_reference;
 using namespace vela::slot_population;
 for(double eta:{-6.,-2.,2.,6.,12.,24.,48.})for(double t:{-1e-4,-1e-12,1e-14,1e-6}){
  const auto base=polynomial(Wide(eta));Wide d=t;
  for(int k=0;k<10;++k){auto f=polynomial(Wide(eta)+d);d-=(f.first-base.first-base.second*Wide(t))/f.second;}
  const Wide exact=d-Wide(t);const double got=inverseTangentRemainder(eta,t);
  CHECK(abs((Wide(got)-exact)/exact)<Wide("2e-11"));
 }
 // The dilute exact inverse is log1p; cancellation does not set the step to zero.
 for(double t:{-1e-4,-1e-14,1e-14,1e-4}){
  const Wide exact=log(1+Wide(t))-Wide(t);
  CHECK(abs((Wide(inverseTangentRemainder(-700.,t))-exact)/exact)<Wide("1e-14"));
 }
 const double vt=.025851999786435535,xp=-4.3236142232904707e-16/vt,xq=1.0846837444346201e-29/vt;
 for(double sign:{-1.,1.}) {
  const double step=tangentQuasiFermiStep(-2.,sign,vt,xp,xq,vt*xp);
  const Wide t=Wide(sign)*(Wide(xq)-xp);const auto base=polynomial(Wide(-2));Wide d=t;
  for(int k=0;k<10;++k){auto f=polynomial(Wide(-2)+d);d-=(f.first-base.first-base.second*t)/f.second;}
  const Wide expected=Wide(vt)*xq+Wide(sign)*vt*(d-t);
  CHECK(abs((Wide(step)-expected)/expected)<Wide("1e-11"));
  CHECK(step!=0.);
  const double written=tangentQuasiFermiStep(-2.,sign,vt,xp,xq,0.);
  CHECK(written==Catch::Approx(step-vt*xp).epsilon(1e-14));
 }
 CHECK(inverseTangentRemainder(0.,0.)==0.);
 CHECK_THROWS(inverseTangentRemainder(0.,.1));
}

TEST_CASE("Effective-band SG keeps Fermi populations and conservative thermal weights",
          "[carrier_statistics][sg][effective_band]")
{
    const Real vt = constants::Vt_300;
    const Real n0 = 2.8e25 * fermiDiracHalf(9.0);
    const Real n1 = 2.8e25 * fermiDiracHalf(10.0);
    const Real q0 = .02, q1 = .025, coef = .03;
    const Real a = std::log(n1/n0) + (q1-q0)/vt;
    // Independent two-population SG expression (no cancellation in this case).
    const auto B=[](Real v) { return v/std::expm1(v); };
    const Real expected = coef*(B(-a)*n0-B(a)*n1);
    const Real flux = sgElectronFermiDiracContinuityFlux(n0,n1,9.,10.,.123,q0,q1,vt,coef,true);
    REQUIRE(flux == Approx(expected).epsilon(2e-13));
    REQUIRE(sgElectronFermiDiracContinuityFlux(n1,n0,10.,9.,-.123,q1,q0,vt,coef,true)
        == Approx(-flux).epsilon(2e-13));
    REQUIRE(sgHoleFermiDiracContinuityFlux(n0,n1,9.,10.,-.123,-q0,-q1,vt,coef,true)
        == Approx(flux).epsilon(2e-13));
    REQUIRE(sgElectronFermiDiracContinuityFlux(n0,n1,9.,10.,.123,q0,q0,vt,coef,true)==0.);
    REQUIRE(sgHoleFermiDiracContinuityFlux(n0,n1,9.,10.,.123,q0,q0,vt,coef,true)==0.);
    REQUIRE(sgElectronFermiDiracQuantumContinuityFlux(n0,n1,n0,n1,9.,10.,.123,q0,q1,vt,coef,true)
        == Approx(flux).epsilon(2e-13));
    REQUIRE(static_cast<Real>(sgElectronFermiDiracQuantumContinuityFluxLongDouble(
        n1,n0,n1,9.,10.,.123,static_cast<long double>(q1)-q0,vt,coef,true))
        == Approx(flux).epsilon(2e-13));
    REQUIRE(static_cast<Real>(sgHoleFermiDiracContinuityFluxLongDouble(
        n0,n1,9.,10.,-.123,static_cast<long double>(q0)-q1,vt,coef,true))
        == Approx(flux).epsilon(2e-13));
    // Degenerate population law has not been replaced by Boltzmann statistics.
    REQUIRE(std::abs(flux/sgElectronFermiDiracContinuityFlux(n0,n1,9.,10.,.123,q0,q1,vt,coef)-1.)>1e-3);
    for (Real eta : {-30.,-20.}) {
        const Real c0=2.8e25*fermiDiracHalf(eta), c1=2.8e25*fermiDiracHalf(eta+.4);
        const Real drift=vt*.4+q1-q0;
        const Real original=sgElectronFermiDiracContinuityFlux(c0,c1,eta,eta+.4,drift,q0,q1,vt,coef);
        const Real effective=sgElectronFermiDiracContinuityFlux(c0,c1,eta,eta+.4,drift,q0,q1,vt,coef,true);
        REQUIRE(effective==Approx(original).epsilon(2e-8));
    }
    REQUIRE(std::isfinite(sgEffectiveBandArgument(1e-300,1e300,.1)));
    REQUIRE(std::isfinite(sgElectronFermiDiracContinuityFlux(1e-200,1e200,0.,0.,0.,0.,20.,vt,1e-200,true)));
}

TEST_CASE("Fermi transport choice is explicit and validated", "[carrier_statistics][effective_band]")
{
    REQUIRE(CarrierStatisticsConfig{}.transport=="generalized_einstein");
    REQUIRE(usesFermiDirac(CarrierStatisticsConfig{"fermi_dirac","effective_band_sg"}));
    REQUIRE_THROWS(carrierStatisticsModel(CarrierStatisticsConfig{"boltzmann","effective_band_sg"}));
    REQUIRE_THROWS(carrierStatisticsModel(CarrierStatisticsConfig{"fermi_dirac","misspelled"}));
}

TEST_CASE("Fermi-Dirac half integral has the documented normalization",
          "[carrier_statistics][fermi_dirac]")
{
    REQUIRE(fermiDiracHalf(0.0) == Approx(0.765147).epsilon(4.0e-3));
    REQUIRE(fermiDiracHalf(-12.0) == Approx(std::exp(-12.0)).epsilon(1.0e-3));
    REQUIRE(fermiDiracHalfDerivative(-12.0) ==
            Approx(std::exp(-12.0)).epsilon(2.0e-3));
    REQUIRE(fermiDiracHalfDerivative(4.0) > 0.0);
}

namespace {
long double referenceFermiHalf(long double eta, bool derivative = false)
{
    const auto integrand = [&](long double t) {
        const long double z = t*t - eta;
        const long double f = z > 0 ? std::exp(-z)/(1+std::exp(-z))
                                    : 1/(1+std::exp(z));
        const long double complement = z > 0 ? 1/(1+std::exp(-z))
                                             : std::exp(z)/(1+std::exp(z));
        return t*t*f*(derivative ? complement : 1.L);
    };
    const auto integral = boost::math::quadrature::gauss_kronrod<long double,61>::integrate(
        integrand, 0.L, std::sqrt(std::max(0.L,eta)+80.L), 12, 1e-16L);
    return 4/std::sqrt(std::acos(-1.L))*integral;
}
}

TEST_CASE("Fermi functions and derivatives agree with independent defining integrals",
          "[carrier_statistics][fermi_dirac]")
{
    for (Real eta=-40.;eta<=100.;eta+=.75) {
        INFO("eta=" << eta);
        CHECK(fermiDiracHalf(eta)==Approx(static_cast<Real>(referenceFermiHalf(eta))).epsilon(2e-12));
        CHECK(fermiDiracHalfDerivative(eta)==Approx(static_cast<Real>(referenceFermiHalf(eta,true))).epsilon(2e-11));
        CHECK(fermiDiracHalfDerivative(eta)>0.);
    }
    for (Real boundary : {-8.,-4.,0.,4.,8.,16.,32.,64.}) {
        for (Real eta : {std::nextafter(boundary,-INFINITY),boundary,std::nextafter(boundary,INFINITY)}) {
            INFO("boundary=" << boundary);
            CHECK(fermiDiracHalf(eta)==Approx(static_cast<Real>(referenceFermiHalf(eta))).epsilon(2e-12));
            CHECK(fermiDiracHalfDerivative(eta)==Approx(static_cast<Real>(referenceFermiHalf(eta,true))).epsilon(2e-11));
        }
    }
}

TEST_CASE("Fermi inverse remains accurate across extreme carrier density ratios",
          "[carrier_statistics][fermi_dirac]")
{
    for (Real value : {1e-300,1e-100,1e-20,1e-12,1e-8,.01,.5,1.,10.,1e3,1e20,1e200}) {
        const Real eta=inverseFermiDiracHalf(value);
        INFO("density ratio=" << value);
        REQUIRE(std::isfinite(eta));
        CHECK(std::abs(fermiDiracHalf(eta)/value-1.)<2e-12);
        CHECK(fermiDiracHalfDerivative(eta)>0.);
    }
    CHECK(std::isnan(fermiDiracHalf(std::numeric_limits<Real>::quiet_NaN())));
    CHECK(fermiDiracHalf(-INFINITY)==0.);
}

TEST_CASE("Fermi BGN correction uses the same accurately inverted statistics in SI and TCAD",
          "[carrier_statistics][fermi_dirac][bgn]")
{
    const Real vt=constants::kb*300./constants::q;
    for (Real ratio : {1e-6,.01,1.,10.}) {
        long double lo=std::log(ratio),hi=30.;
        for (int i=0;i<65;++i) {
            const auto mid=(lo+hi)/2;
            if (referenceFermiHalf(mid)>ratio) hi=mid; else lo=mid;
        }
        const Real expected=vt*static_cast<Real>((lo+hi)/2-std::log(static_cast<long double>(ratio)));
        const Real si=fermiStatisticsBandgapCorrection(ratio*2.8e25,0.,2.8e25,3.1e25,vt);
        const Real tcad=fermiStatisticsBandgapCorrection(ratio*2.8e19,0.,2.8e19,3.1e19,vt);
        CHECK(si==Approx(expected).margin(2e-14));
        CHECK(tcad==Approx(si).margin(2e-15));
    }
}

TEST_CASE("Fermi-Dirac half inverse round trips from nondegenerate to degenerate",
          "[carrier_statistics][fermi_dirac]")
{
    for (const Real eta : {-20.0, -5.0, 0.0, 2.0, 8.0, 25.0}) {
        const Real value = fermiDiracHalf(eta);
        REQUIRE(inverseFermiDiracHalf(value) == Approx(eta).margin(2.0e-9));
    }
}

TEST_CASE("Fermi-Dirac density and quasi-Fermi inversion are consistent",
          "[carrier_statistics][fermi_dirac]")
{
    CarrierStatisticsConfig statistics{"fermi_dirac"};
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.68e17;
    const Real Nc = 2.8e25;
    const Real Nv = 1.04e25;
    const Real psi = 0.61;
    const Real phin = 0.013;
    const Real phip = -0.007;

    const Real n = electronDensity(ni, Nc, psi, phin, Vt, statistics);
    const Real p = holeDensity(ni, Nv, psi, phip, Vt, statistics);
    REQUIRE(electronQuasiFermiPotential(ni, Nc, psi, n, Vt, statistics) ==
            Approx(phin).margin(2.0e-12));
    REQUIRE(holeQuasiFermiPotential(ni, Nv, psi, p, Vt, statistics) ==
            Approx(phip).margin(2.0e-12));
}

TEST_CASE("Fermi-Dirac Ohmic state solves charge neutrality numerically",
          "[carrier_statistics][fermi_dirac][ohmic]")
{
    CarrierStatisticsConfig statistics{"fermi_dirac"};
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.683405723e17;
    const Real Nc = 2.8e25;
    const Real Nv = 1.04e25;
    const Real netDoping = 3.256916861e26;

    const EquilibriumCarrierState state = equilibriumCarrierState(
        netDoping, ni, Nc, Nv, Vt, statistics);
    REQUIRE((state.n - state.p) / netDoping == Approx(1.0).epsilon(2.0e-12));
    REQUIRE(state.n > Nc);
    REQUIRE(state.p > 0.0);

    const EquilibriumCarrierState boltzmann = equilibriumCarrierState(
        netDoping, ni, Nc, Nv, Vt, CarrierStatisticsConfig{});
    REQUIRE(state.potential > boltzmann.potential);
}

TEST_CASE("Resolved carrier statistics model reproduces the configuration overloads bit for bit",
          "[carrier_statistics][fermi_dirac]")
{
    const CarrierStatisticsConfig boltzmann{"boltzmann"};
    const CarrierStatisticsConfig fermiDirac{"fermi_dirac"};
    REQUIRE(carrierStatisticsModel(boltzmann) == CarrierStatisticsModel::Boltzmann);
    REQUIRE(carrierStatisticsModel(fermiDirac) == CarrierStatisticsModel::FermiDirac);
    REQUIRE_FALSE(usesFermiDirac(CarrierStatisticsModel::Boltzmann));
    REQUIRE(usesFermiDirac(CarrierStatisticsModel::FermiDirac));
    REQUIRE_THROWS_AS(carrierStatisticsModel(CarrierStatisticsConfig{"maxwell"}),
                      std::invalid_argument);

    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.683405723e17;
    const Real Nc = 2.8e25;
    const Real Nv = 1.04e25;

    for (const CarrierStatisticsConfig& config : {boltzmann, fermiDirac}) {
        const CarrierStatisticsModel model = carrierStatisticsModel(config);
        for (const Real psi : {-0.8, -0.05, 0.0, 0.35, 1.1}) {
            for (const Real phi : {-0.4, 0.0, 0.6}) {
                REQUIRE(electronDensity(ni, Nc, psi, phi, Vt, model) ==
                        electronDensity(ni, Nc, psi, phi, Vt, config));
                REQUIRE(holeDensity(ni, Nv, psi, phi, Vt, model) ==
                        holeDensity(ni, Nv, psi, phi, Vt, config));
                REQUIRE(electronDensityDerivativeEta(ni, Nc, psi, phi, Vt, model) ==
                        electronDensityDerivativeEta(ni, Nc, psi, phi, Vt, config));
                REQUIRE(holeDensityDerivativeEta(ni, Nv, psi, phi, Vt, model) ==
                        holeDensityDerivativeEta(ni, Nv, psi, phi, Vt, config));
            }
        }
        for (const Real netDoping : {-3.2569e26, -1.0e20, 0.0, 1.0e20, 3.2569e26}) {
            const EquilibriumCarrierState fromModel =
                equilibriumCarrierState(netDoping, ni, Nc, Nv, Vt, model);
            const EquilibriumCarrierState fromConfig =
                equilibriumCarrierState(netDoping, ni, Nc, Nv, Vt, config);
            REQUIRE(fromModel.potential == fromConfig.potential);
            REQUIRE(fromModel.n == fromConfig.n);
            REQUIRE(fromModel.p == fromConfig.p);
            // For Fermi-Dirac the returned densities must equal a fresh
            // evaluation at the returned potential; the solver reuses the
            // densities it already evaluated for the final residual.
            if (usesFermiDirac(model)) {
                REQUIRE(fromModel.n ==
                        electronDensity(ni, Nc, fromModel.potential, 0.0, Vt, model));
                REQUIRE(fromModel.p ==
                        holeDensity(ni, Nv, fromModel.potential, 0.0, Vt, model));
            }
        }
    }
}

TEST_CASE("Generalized Fermi-Dirac SG preserves flat quasi-Fermi equilibrium",
          "[carrier_statistics][fermi_dirac][sg]")
{
    CarrierStatisticsConfig statistics{"fermi_dirac"};
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni0 = 1.0e16;
    const Real ni1 = 2.0e17;
    const Real Nc = 2.8e25;
    const Real Nv = 1.04e25;
    const Real psi0 = 0.42;
    const Real psi1 = 0.63;
    const Real qf = 0.017;
    const Real n0 = electronDensity(ni0, Nc, psi0, qf, Vt, statistics);
    const Real n1 = electronDensity(ni1, Nc, psi1, qf, Vt, statistics);
    const Real etaN0 = (psi0 - qf) / Vt + std::log(ni0 / Nc);
    const Real etaN1 = (psi1 - qf) / Vt + std::log(ni1 / Nc);
    const Real electronDrift = psi1 - psi0 + Vt * std::log(ni1 / ni0);
    REQUIRE(sgElectronFermiDiracContinuityFlux(
        n0, n1, etaN0, etaN1, electronDrift, qf, qf, Vt, 3.0) == 0.0);

    const Real p0 = holeDensity(ni0, Nv, psi0, qf, Vt, statistics);
    const Real p1 = holeDensity(ni1, Nv, psi1, qf, Vt, statistics);
    const Real etaP0 = (qf - psi0) / Vt + std::log(ni0 / Nv);
    const Real etaP1 = (qf - psi1) / Vt + std::log(ni1 / Nv);
    const Real holeDrift = psi1 - psi0 + Vt * std::log(ni0 / ni1);
    REQUIRE(sgHoleFermiDiracContinuityFlux(
        p0, p1, etaP0, etaP1, holeDrift, qf, qf, Vt, 3.0) == 0.0);

    const Real adjacentQf = std::nextafter(qf, std::numeric_limits<Real>::infinity());
    const Real adjacentN1 = electronDensity(
        ni1, Nc, psi1, adjacentQf, Vt, statistics);
    const Real adjacentEtaN1 = (psi1 - adjacentQf) / Vt + std::log(ni1 / Nc);
    const Real adjacentElectronFlux = sgElectronFermiDiracContinuityFlux(
        n0, adjacentN1, etaN0, adjacentEtaN1, electronDrift,
        qf, adjacentQf, Vt, 3.0);
    REQUIRE(std::isfinite(adjacentElectronFlux));
    REQUIRE(adjacentElectronFlux != 0.0);
    const Real quantumFactor = 0.25;
    const Real quantumElectronFlux =
        sgElectronFermiDiracQuantumContinuityFlux(
            quantumFactor * n0, quantumFactor * adjacentN1,
            n0, adjacentN1, etaN0, adjacentEtaN1, electronDrift,
            qf, adjacentQf, Vt, 3.0);
    REQUIRE(quantumElectronFlux / adjacentElectronFlux ==
            Catch::Approx(quantumFactor).epsilon(1.0e-12));

    const Real adjacentP1 = holeDensity(
        ni1, Nv, psi1, adjacentQf, Vt, statistics);
    const Real adjacentEtaP1 = (adjacentQf - psi1) / Vt + std::log(ni1 / Nv);
    const Real adjacentHoleFlux = sgHoleFermiDiracContinuityFlux(
        p0, adjacentP1, etaP0, adjacentEtaP1, holeDrift,
        qf, adjacentQf, Vt, 3.0);
    REQUIRE(std::isfinite(adjacentHoleFlux));
    REQUIRE(adjacentHoleFlux != 0.0);
}

TEST_CASE("Generalized SRH state retains a near-equilibrium quasi-Fermi split",
          "[carrier_statistics][srh][equilibrium]")
{
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.68e16;
    const Real split = 1.0e-18;
    const GeneralizedSrhCarrierState state = generalizedSrhCarrierState(
        1.0e21, 1.0e21, ni, 2.8e25, 1.04e25, split, Vt,
        CarrierStatisticsConfig{"boltzmann"});

    REQUIRE(state.electronDegeneracy == 1.0);
    REQUIRE(state.holeDegeneracy == 1.0);
    REQUIRE(state.excessProduct != 0.0);
    REQUIRE(state.excessProduct ==
            Approx(ni * ni * std::expm1(split / Vt)).epsilon(1.0e-14));
}

TEST_CASE("Generalized SRH state resolves deep-depletion generation",
          "[carrier_statistics][srh][deep_depletion]")
{
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.68e16;
    const Real split = -1.6;
    const GeneralizedSrhCarrierState state = generalizedSrhCarrierState(
        1.0, 1.0, ni, 2.8e25, 1.04e25, split, Vt,
        CarrierStatisticsConfig{"boltzmann"});

    REQUIRE(std::isfinite(state.excessProduct));
    REQUIRE(state.excessProduct < 0.0);
    REQUIRE(state.excessProduct == Approx(-ni * ni).epsilon(1.0e-14));
}

TEST_CASE("Generalized Fermi SRH factors reproduce the carrier product identity",
          "[carrier_statistics][fermi_dirac][srh]")
{
    const CarrierStatisticsConfig statistics{"fermi_dirac"};
    const Real Vt = constants::kb * 300.0 / constants::q;
    const Real ni = 1.68e16;
    const Real Nc = 2.8e25;
    const Real Nv = 1.04e25;
    const Real psi = 0.75;
    const Real phin = -0.05;
    const Real phip = 0.02;
    const Real n = electronDensity(ni, Nc, psi, phin, Vt, statistics);
    const Real p = holeDensity(ni, Nv, psi, phip, Vt, statistics);
    const GeneralizedSrhCarrierState state = generalizedSrhCarrierState(
        n, p, ni, Nc, Nv, phip - phin, Vt, statistics);

    REQUIRE(state.electronDegeneracy > 0.0);
    REQUIRE(state.holeDegeneracy > 0.0);
    REQUIRE(state.electronDegeneracy != Approx(1.0).margin(1.0e-3));
    REQUIRE(n * p ==
            Approx(state.equilibriumProduct * std::exp((phip - phin) / Vt))
                .epsilon(3.0e-10));
    REQUIRE(state.excessProduct ==
            Approx(n * p - state.equilibriumProduct).epsilon(3.0e-10));
}

TEST_CASE("Fermi half second derivative is positive and matches derivative slopes", "[statistics][fermi][thermal]") {
    for(Real eta:{-30.,-8.1,-8.,-7.9,-4.,0.,4.,8.,16.,32.,63.9,64.,64.1,100.,1e4}) {
        const Real second=fermiDiracHalfSecondDerivative(eta);
        CHECK(second>0.);
        for(Real fraction:{1.,.25}) {
            Real step=1e-4*std::max(1.,std::abs(eta))*fraction;
            const Real fd=(fermiDiracHalfDerivative(eta+step)-fermiDiracHalfDerivative(eta-step))/(2.*step);
            CHECK(second==Catch::Approx(fd).epsilon(2e-6));
        }
    }
}

TEST_CASE("SLOT relative Fermi update preserves identity and depleted recovery", "[slot_population]") {
    for(double eta : {-100.,-40.,-8.,0.,4.,16.}) {
        double z=vela::fermiDiracHalf(eta);
        REQUIRE(vela::slot_population::relativeInverse(z,0.)==0.);
        for(double change : {-0.999,-0.5,1e-12,0.5}) {
            double d=vela::slot_population::relativeInverse(z,change);
            REQUIRE(std::isfinite(d));
            REQUIRE(d*change>0.);
            REQUIRE(vela::fermiDiracHalf(eta+d)/z == Catch::Approx(1+change).epsilon(1e-9));
        }
    }
    REQUIRE_THROWS(vela::slot_population::relativeInverse(1.,-1.));
    REQUIRE_THROWS(vela::slot_population::fermiChange(0.,1e20));
    REQUIRE_THROWS(vela::slot_population::relativeInverse(0.,.1));
}

TEST_CASE("Dilute Fermi densities avoid subnormal intermediate loss", "[carrier_statistics][fermi_dirac][dilute_mapping]") {
    constexpr auto model=CarrierStatisticsModel::FermiDirac;
    // ni=DOS, Vt=1 and psi=0 expose eta exactly in both public carrier APIs.
    // The extended reference retains exp(eta) precision before multiplication.
    for(double dos:{3.1e19,2.8e25})for(double eta:{-699.999,-700.,-700.001,-729.,-750.}) {
        const double expected=static_cast<double>(static_cast<long double>(dos)*std::exp(static_cast<long double>(eta)));
        const double n=electronDensity(dos,dos,0.,-eta,1.,model);
        const double p=holeDensity(dos,dos,0.,eta,1.,model);
        REQUIRE(expected>0.);
        CHECK(std::abs(n/expected-1.)<3e-13);
        CHECK(std::abs(p/expected-1.)<3e-13);
        CHECK(std::abs(electronDensityDerivativeEta(dos,dos,0.,-eta,1.,model)/expected-1.)<3e-13);
        CHECK(std::abs(holeDensityDerivativeEta(dos,dos,0.,eta,1.,model)/expected-1.)<3e-13);
        const double h=1e-3;
        const double slope=(holeDensity(dos,dos,0.,eta+h,1.,model)-holeDensity(dos,dos,0.,eta-h,1.,model))/(2*h);
        CHECK(std::abs(slope/p-1.)<2e-7);
    }
    // Ordinary statistics remain bitwise unchanged, including the degenerate region.
    for(double eta:{-700.,-30.,-8.,0.,20.}) {
        CHECK(holeDensity(3.1e19,3.1e19,0.,eta,1.,model)==3.1e19*fermiDiracHalf(eta));
        CHECK(holeDensityDerivativeEta(3.1e19,3.1e19,0.,eta,1.,model)==3.1e19*fermiDiracHalfDerivative(eta));
    }
    // Recorded 6251 scale: construct the reference without first rounding p/DOS.
    const double target=5.678873e-298,dos=3.1e19;
    const double eta=static_cast<double>(std::log(static_cast<long double>(target)/dos));
    CHECK(std::abs(holeDensity(dos,dos,0.,eta,1.,model)/target-1.)<3e-13);
    CHECK(holeDensity(dos,dos,0.,-1000.,1.,model)==0.); // no artificial positive floor
}

TEST_CASE("Relative dilute Fermi inverse keeps ratios through subnormal endpoints", "[slot_population][dilute_mapping]") {
    using vela::slot_population::relativeInverse;
    for(double z:{std::numeric_limits<double>::denorm_min(),1.83e-317,1e-300,1e-20,1e-14}) {
        CHECK(relativeInverse(z,0.)==0.);
        for(double rel:{-.999,-.5,1e-12,.5,10.}) {
            const double d=relativeInverse(z,rel);
            // Compare a ratio, not an absolute tiny residual whose tolerance can vanish.
            CHECK(std::abs(std::expm1(d)-rel)<2e-12*std::abs(rel));
        }
    }
    // Starting dilute is insufficient: a non-dilute endpoint must retain Fermi statistics.
    const double z=1e-14,rel=1e13-1.;
    const double d=relativeInverse(z,rel);
    CHECK(std::abs(fermiDiracHalf(inverseFermiDiracHalf(z)+d)/(z*(1+rel))-1.)<1e-11);
    CHECK(std::abs(d-std::log1p(rel))>1e-3);
    // A projected finite density step must survive both inverse and forward mapping.
    constexpr double dos=3.1e19,p=5.678873e-295,relProjection=-.999;
    const double eta=static_cast<double>(std::log(static_cast<long double>(p)/dos));
    const double nextEta=eta+relativeInverse(p/dos,relProjection);
    const double next=holeDensity(dos,dos,0.,nextEta,1.,CarrierStatisticsModel::FermiDirac);
    CHECK(std::abs(next/(p*(1+relProjection))-1.)<3e-13);
}

TEST_CASE("SLOT primary-only acceptance separates physical reports from stopping", "[slot_population]") {
    using namespace vela::slot_population;
    for(bool closure:{false,true})for(bool boundary:{false,true}) {
        CHECK(accepts(1e-10,20.,closure,boundary,AcceptancePolicy::PrimaryOnly));
        CHECK(accepts(1e-7,.9,closure,boundary,AcceptancePolicy::PrimaryOnly));
        CHECK_FALSE(accepts(1.01e-10,1.,closure,boundary,AcceptancePolicy::PrimaryOnly));
        CHECK(accepts(1e-10,20.,closure,boundary,AcceptancePolicy::PrimaryWithGuards)==(closure&&boundary));
    }
    CHECK_FALSE(accepts(-1.,.1,true,true,AcceptancePolicy::PrimaryOnly));
    CHECK_FALSE(accepts(1e-12,std::numeric_limits<double>::infinity(),true,true,AcceptancePolicy::PrimaryOnly));
}
