#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

using Catch::Approx;
using namespace vela;

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
