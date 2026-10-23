#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/experimental/ReferencedPopulation.h"
#include "vela/physics/experimental/FermiSgDifferential.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/physics/BandgapNarrowing.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

using Catch::Approx;
using namespace vela;

TEST_CASE("Referenced population preserves finite coupled increments before charge formation",
          "[carrier_statistics][referenced_population]")
{
    using W = boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    const auto forbiddenFermi = [](const W&) -> W { throw std::logic_error("unexpected Fermi callback"); };
    PopulationReference r{1e10, 1e19, 0.025, 0.025, 0., 0., false, false};
    const W psi = -17.2, qf = W("8e-11"), dq = W("9e-16");
    const W base = referencedPopulation(r, psi, qf, forbiddenFermi);
    for (const W h : {W("0.001"), W("0.0003"), W("0.0001")}) {
        const W dp = W("0.37")*dq*h, df = dq*h;
        const W increment = referencedPopulationIncrement(r, psi, qf,
            W(psi+dp), W(qf+df), forbiddenFermi);
        const W expected = base*expm1(W((df-dp)*W(r.potentialScale)/W(r.thermalVoltage)));
        REQUIRE(abs((increment-expected)/expected) < W("1e-65"));
        REQUIRE(increment > 0);
        const W reverse = referencedPopulationIncrement(r, W(psi+dp), W(qf+df), psi, qf, forbiddenFermi);
        REQUIRE(reverse == -increment);
    }
    // Double collapse loses the same response; a zero difference is not a zero derivative.
    const auto collapsedBefore = holeDensity(r.intrinsic, double(psi)*r.potentialScale,
        double(qf)*r.potentialScale, r.thermalVoltage);
    const auto collapsedAfter = holeDensity(r.intrinsic, double(psi)*r.potentialScale,
        double(qf+dq)*r.potentialScale, r.thermalVoltage);
    REQUIRE(collapsedBefore == collapsedAfter);
    REQUIRE(referencedPopulationIncrement(r, psi, qf, psi, W(qf+dq), forbiddenFermi)>0);
    REQUIRE(referencedPopulationIncrement(r, psi, qf, W(psi+1), W(qf+1), forbiddenFermi)==0);
    // Finite loss saturates at the old density; this is not a tangent update.
    const W loss=referencedPopulationIncrement(r, psi, qf, psi, W(qf-100), forbiddenFermi);
    REQUIRE(loss<0);REQUIRE(-loss<=base);
    r.thermalVoltage=0;
    REQUIRE_THROWS(referencedPopulation(r, psi, qf, forbiddenFermi));
}

TEST_CASE("Referenced Fermi population uses a value increment and retains both potentials",
          "[carrier_statistics][referenced_population]")
{
    using W = boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    // Independent defining fugacity series, far from the degenerate boundary.
    auto half = [](const W& eta) -> W {
        const W z=exp(eta);W power=z,sum=0;
        for(int k=1;k<=32;++k){sum+=(k%2?1:-1)*power/(W(k)*sqrt(W(k)));power*=z;}
        return sum;
    };
    PopulationReference r{1e10, 1e19, .025, .025, .01, 0., true, true};
    const W psi=.1,qf=.2,delta("1e-20");
    const W increment=referencedPopulationIncrement(r,psi,qf,W(psi+delta),W(qf+2*delta),half);
    REQUIRE(increment<0);
    REQUIRE(referencedPopulationIncrement(r,psi,qf,W(psi+delta),W(qf+delta),half)==0);
    const W base=referencedPopulation(r,psi,qf,half);
    const double ordinary=electronDensity(r.intrinsic,r.densityOfStates,
        double(psi)*r.potentialScale-r.quasiFermiReference,double(qf)*r.potentialScale,
        r.thermalVoltage,CarrierStatisticsModel::FermiDirac);
    REQUIRE(double(base)==Approx(ordinary).epsilon(1e-13));
    r.densityOfStates=0;
    REQUIRE_THROWS(referencedPopulation(r,psi,qf,half));
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

TEST_CASE("Fermi SG differential retains equilibrium and coupled directions", "[carrier_statistics][sg_differential]") {
    using W=boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    auto statistics=[](const W& eta) { W f=exp(eta);return StatisticsJet<W>{f,f,f}; };
    auto evaluate=[&](const std::array<W,4>& x,bool electron) {
        return fermiSgDifferential(x[0],x[1],x[2],x[3],W(".025"),W(3),W(0),W(0),W(1),W(1),W(0),electron,statistics);
    };
    for(bool electron:{false,true})for(const W drop:{W(0),W("1e-20"),W(".01")}) {
        std::array<W,4>x{W(".11"),W(".13"),W(".03"),W(".03")+drop};
        auto a=evaluate(x,electron);
        REQUIRE(abs(a.partial[0]+a.partial[1]+a.partial[2]+a.partial[3])<W("1e-80"));
        if(drop==0){REQUIRE(a.value==0);REQUIRE(a.partial[2]!=0);}
        std::array<W,4> reverse{x[1],x[0],x[3],x[2]};
        REQUIRE(abs(evaluate(reverse,electron).value+a.value)<W("1e-80"));
        std::array<W,4> d{W(".2"),W("-.3"),W(".7"),W("-.1")};W jv=0;
        for(int k=0;k<4;++k)jv+=a.partial[k]*d[k];
        for(W h:{W("1e-6"),W("3e-7"),W("1e-7")}) {
            auto xp=x,xm=x;for(int k=0;k<4;++k){xp[k]+=h*d[k];xm[k]-=h*d[k];}
            W fd=(evaluate(xp,electron).value-evaluate(xm,electron).value)/(2*h);
            REQUIRE(abs((fd-jv)/jv)<W("1e-8"));
        }
    }
    std::array<W,4> flat{W(".1"),W(".1"),W(0),W(0)};
    auto limit=evaluate(flat,true);REQUIRE(limit.midpoint);REQUIRE(limit.value==0);
    REQUIRE(limit.partial[2]==-limit.partial[3]);
}

TEST_CASE("Fermi SG secant and midpoint derivatives match independent series values", "[carrier_statistics][sg_differential]") {
    using W=boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    auto statistics=[](const W& eta) {
        W z=exp(eta),p=z,f=0,d=0,dd=0;
        for(int k=1;k<=120;++k){W t=(k%2?1:-1)*p/(W(k)*sqrt(W(k)));f+=t;d+=k*t;dd+=k*k*t;p*=z;}
        return StatisticsJet<W>{f,d,dd};
    };
    for(bool midpoint:{false,true}) {
        std::array<W,4>x{W(".1"),W(".11"),W(".2"),midpoint?W(".21"):W(".207")};
        auto eval=[&](auto y){return fermiSgDifferential(y[0],y[1],y[2],y[3],W(".025"),W(2),W(0),W(0),W(1),W(1),W(0),true,statistics);};
        auto a=eval(x);REQUIRE(a.midpoint==midpoint);
        std::array<W,4>d{W(".2"),W("-.3"),W("-.1"),W(".4")};W tangent=0;
        for(int k=0;k<4;++k)tangent+=a.partial[k]*d[k];
        for(W h:{W("1e-6"),W("3e-7"),W("1e-7")}) {
            auto xp=x,xm=x;for(int k=0;k<4;++k){xp[k]+=h*d[k];xm[k]-=h*d[k];}
            W fd=(eval(xp).value-eval(xm).value)/(2*h);
            REQUIRE(abs((fd-tangent)/tangent)<W("1e-8"));
        }
    }
}


#include "vela/numerics/ReferencedDirection.h"

TEST_CASE("Referenced directions preserve sub ULP state and boundary changes", "[carrier_statistics][reference_increment]") {
    using W=boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    const double base=10., direction=1e-16, scale=.025;
    for(double h:{1e-3,3e-4,1e-4}) {
        REQUIRE(base+h*direction==base);
        W plus=referencedDirectionValue<W>(base,h,direction);
        W minus=referencedDirectionValue<W>(base,-h,direction);
        REQUIRE(abs((plus-minus)/(2*W(h))-W(direction))<W("1e-90"));
        W bp=referencedBoundaryValue<W>(base,h,direction,scale);
        W bm=referencedBoundaryValue<W>(base,-h,direction,scale);
        REQUIRE(abs((bp-bm)/(2*W(h))-W(direction)/W(scale))<W("1e-90"));
        // A matched state and contact direction leaves their constraint fixed.
        W xp=referencedDirectionValue<W>(base,h,1.);
        W target=referencedBoundaryValue<W>(base,h,scale,scale);
        REQUIRE(abs(xp-target)<W("1e-90"));
    }
    REQUIRE_THROWS_AS(referencedBoundaryValue<W>(0.,1.,1.,0.),std::invalid_argument);
    REQUIRE_THROWS_AS(referencedBoundaryValue<W>(std::numeric_limits<double>::infinity(),1.,1.,1.),std::invalid_argument);
    REQUIRE_THROWS_AS(referencedDirectionValue<W>(0.,1.,std::numeric_limits<double>::infinity()),std::invalid_argument);
}


#include "vela/numerics/ExpansionProductSum.h"

TEST_CASE("Expansion operator preserves coefficients and cancellation before rounding", "[carrier_statistics][expansion_operator]") {
    using W=boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    std::array<W,4> a{W("10000000000000000.25"),W("-10000000000000000"),W("1e-12"),W("-1e-12")};
    std::array<double,4> d{1.,1.,1.3,.7};
    W expected=0;for(int i=0;i<4;++i)expected+=a[i]*W(d[i]);
    for(bool reverse:{false,true}) {
        ExpansionProductSum accurate, rounded;
        for(int j=0;j<4;++j){int i=reverse?3-j:j;accurate.addProduct(expansionCoefficient(a[i]),d[i]);rounded.addProduct({double(a[i]),0.},d[i]);}
        REQUIRE(abs(W(accurate.value())-expected)<W("1e-16"));
        REQUIRE(abs(W(rounded.value())-expected)>W(".24"));
    }
    // Products whose low bits disappear before cancellation need FMA errors too.
    ExpansionProductSum product;
    double x=1.+std::ldexp(1.,-27),y=1.-std::ldexp(1.,-27);
    product.addProduct({x,0.},y);product.add(-1.);
    REQUIRE(product.value()==-std::ldexp(1.,-54));
    REQUIRE_THROWS_AS(product.add(std::numeric_limits<double>::infinity()),std::overflow_error);
}


TEST_CASE("Expansion operator retains sub ULP direction corrections", "[carrier_statistics][expansion_operator]") {
    using W=boost::multiprecision::cpp_dec_float_100;
    using namespace vela::experimental;
    W a("10000000000000000.25"),b("-10000000000000000");
    W x=W(1)+W("3.75e-17"),y=W(1)-W("3.75e-17");
    REQUIRE(double(x)==1.);REQUIRE(double(y)==1.);
    ExpansionProductSum paired,collapsed;
    paired.addProduct(expansionCoefficient(a),expansionCoefficient(x));
    paired.addProduct(expansionCoefficient(b),expansionCoefficient(y));
    collapsed.addProduct(expansionCoefficient(a),double(x));
    collapsed.addProduct(expansionCoefficient(b),double(y));
    W expected=a*x+b*y;
    REQUIRE(abs(W(paired.value())-expected)<W("1e-16"));
    REQUIRE(abs(W(collapsed.value())-expected)>W(".74"));
}
