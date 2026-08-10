#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/discretization/ScharfetterGummel.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using Catch::Approx;
using namespace vela;

TEST_CASE("Bednarczyk Fermi-Dirac half integral has the documented normalization",
          "[carrier_statistics][fermi_dirac]")
{
    REQUIRE(fermiDiracHalf(0.0) == Approx(0.765147).epsilon(4.0e-3));
    REQUIRE(fermiDiracHalf(-12.0) == Approx(std::exp(-12.0)).epsilon(1.0e-3));
    REQUIRE(fermiDiracHalfDerivative(-12.0) ==
            Approx(std::exp(-12.0)).epsilon(2.0e-3));
    REQUIRE(fermiDiracHalfDerivative(4.0) > 0.0);
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

    const Real adjacentP1 = holeDensity(
        ni1, Nv, psi1, adjacentQf, Vt, statistics);
    const Real adjacentEtaP1 = (adjacentQf - psi1) / Vt + std::log(ni1 / Nv);
    const Real adjacentHoleFlux = sgHoleFermiDiracContinuityFlux(
        p0, adjacentP1, etaP0, adjacentEtaP1, holeDrift,
        qf, adjacentQf, Vt, 3.0);
    REQUIRE(std::isfinite(adjacentHoleFlux));
    REQUIRE(adjacentHoleFlux != 0.0);
}
