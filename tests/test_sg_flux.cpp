#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace vela;
using Catch::Approx;

TEST_CASE("Homogeneous SG preserves sub-femtovolt conductive response", "[sg][pn2d]")
{
    // An ohmic quasi-neutral region carries current through increments far
    // smaller than a thermal voltage. Its linear response must not disappear
    // when two separately evaluated exponentials both round to one.
    const double vt = 0.02585;
    const double ni = 1.075e10;
    const double psi = 0.415;
    const double coef = 0.73;
    const double density = ni * std::exp(psi / vt);
    for (const double drop : {1e-12, 1e-16, 1e-20, -1e-20}) {
        const double expected = coef * density * drop / vt;
        const double electron = sgElectronContinuityFluxFromQuasiFermiStable(
            ni, psi, psi, 0.0, drop, vt, coef);
        const double hole = sgHoleContinuityFluxFromQuasiFermiStable(
            ni, -psi, -psi, 0.0, drop, vt, coef);
        REQUIRE(electron == Approx(expected).epsilon(1e-9));
        REQUIRE(hole == Approx(-expected).epsilon(1e-9));
        REQUIRE(sgElectronContinuityFluxFromQuasiFermiStable(
            ni, psi, psi, drop, 0.0, vt, coef) == Approx(-electron).epsilon(1e-9));
    }
    REQUIRE(sgElectronContinuityFluxFromQuasiFermiStable(
        ni, psi, psi, 0.0, 0.0, vt, coef) == 0.0);
}

// ---------------------------------------------------------------------------
// Pure diffusion (psi_j = psi_i -> dpsi = 0, u = 0, B(0) = 1)
// ---------------------------------------------------------------------------

TEST_CASE("SG electron flux: pure diffusion (dpsi=0)", "[sg]")
{
    const double Vt = 0.02585; // ~300 K thermal voltage
    const double mu = 0.135;   // Si electron mobility
    const double h  = 1.0e-6;  // 1 um edge

    const double n0 = 1.0e23;  // high concentration (n-side)
    const double n1 = 1.0e16;  // low concentration (p-side)

    // With dpsi = 0: Jn = mu*Vt/h * (n1 - n0)
    const double Jn = sgElectronFlux(n0, n1, 0.0, Vt, mu, h);
    const double expected = mu * Vt / h * (n1 - n0);
    REQUIRE(Jn == Approx(expected).epsilon(1.0e-10));

    // Direction: electrons diffuse from high n0 to low n1 (i->j direction
    // has conventional current flowing j->i) -> Jn < 0
    REQUIRE(Jn < 0.0);
}

TEST_CASE("SG hole flux: pure diffusion (dpsi=0)", "[sg]")
{
    const double Vt = 0.02585;
    const double mu = 0.048;
    const double h  = 1.0e-6;

    const double p0 = 1.0e16;  // low (n-side)
    const double p1 = 1.0e23;  // high (p-side)

    // With dpsi = 0: Jp = mup*Vt/h * (p1 - p0)  (same as diffusion)
    const double Jp = sgHoleFlux(p0, p1, 0.0, Vt, mu, h);
    const double expected = mu * Vt / h * (p1 - p0);
    REQUIRE(Jp == Approx(expected).epsilon(1.0e-10));
    REQUIRE(Jp > 0.0); // holes diffuse from p1->p0 -> conventional current i->j positive
}

// ---------------------------------------------------------------------------
// Equilibrium: with Boltzmann distribution and built-in field, Jn = 0
// ---------------------------------------------------------------------------

TEST_CASE("SG electron flux: zero at equilibrium", "[sg]")
{
    // At thermal equilibrium: n = ni * exp(psi/Vt)
    // So B(u)*n1 - B(-u)*n0 should be zero when psi_j = Vt*ln(n1/ni)
    // and psi_i = Vt*ln(n0/ni).
    const double Vt  = 0.02585;
    const double mu  = 0.135;
    const double h   = 1.0e-6;
    const double ni  = 1.0e16;
    const double n0  = 1.0e23;  // n-side
    const double n1  = 1.0e10;  // p-side (minority)

    const double psi0 = Vt * std::log(n0 / ni);
    const double psi1 = Vt * std::log(n1 / ni);
    const double dpsi = psi1 - psi0;
    const double u    = dpsi / Vt;

    // Jn = mu*Vt/h * (B(u)*n1 - B(-u)*n0)
    //    = mu*Vt/h * (B(u)*ni*exp(psi1/Vt) - B(-u)*ni*exp(psi0/Vt))
    // At equilibrium this equals zero.
    const double Jn = sgElectronFlux(n0, n1, dpsi, Vt, mu, h);
    REQUIRE(std::abs(Jn) < 1.0e8); // should be essentially zero (relative to large n)
    // Normalised check: |Jn| << |diffusion flux|
    const double diffFlux = mu * Vt / h * (n1 + n0); // magnitude scale
    REQUIRE(std::abs(Jn) / diffFlux < 1.0e-10);
}

// ---------------------------------------------------------------------------
// Antisymmetry: swapping i and j should negate the flux
// ---------------------------------------------------------------------------

TEST_CASE("SG fluxes: antisymmetry J_ji = -J_ij", "[sg]")
{
    const double Vt   = 0.02585;
    const double mu   = 0.135;
    const double h    = 1.0e-6;
    const double n0   = 2.0e20;
    const double n1   = 5.0e18;
    const double dpsi = 0.3;   // psi_j - psi_i

    const double Jn_ij = sgElectronFlux(n0, n1,  dpsi, Vt, mu, h);
    const double Jn_ji = sgElectronFlux(n1, n0, -dpsi, Vt, mu, h);

    REQUIRE(Jn_ij + Jn_ji == Approx(0.0).margin(std::abs(Jn_ij) * 1.0e-12));

    const double mup  = 0.048;
    const double p0   = 1.0e21;
    const double p1   = 1.0e18;
    const double Jp_ij = sgHoleFlux(p0, p1,  dpsi, Vt, mup, h);
    const double Jp_ji = sgHoleFlux(p1, p0, -dpsi, Vt, mup, h);

    REQUIRE(Jp_ij + Jp_ji == Approx(0.0).margin(std::abs(Jp_ij) * 1.0e-12));
}

// ---------------------------------------------------------------------------
// No NaN / Inf for extreme inputs
// ---------------------------------------------------------------------------

TEST_CASE("SG fluxes: finite for large dpsi", "[sg]")
{
    const double Vt = 0.02585;
    const double mu = 0.135;
    const double h  = 1.0e-6;

    for (double dpsi : {-50.0, -20.0, 20.0, 50.0}) {
        const double Jn = sgElectronFlux(1.0e16, 1.0e16, dpsi, Vt, mu, h);
        REQUIRE(std::isfinite(Jn));
        const double Jp = sgHoleFlux   (1.0e16, 1.0e16, dpsi, Vt, 0.048, h);
        REQUIRE(std::isfinite(Jp));
    }
}

TEST_CASE("SG quasi-Fermi fluxes stay finite for large potentials", "[sg][coupled]")
{
    const double Vt = constants::Vt_300;
    const double ni = 1.0e16;
    const double coef = 1.0;

    const double electronFlux = sgElectronContinuityFluxFromQuasiFermi(
        ni,
        50.0,
        0.0,
        0.0,
        0.0,
        Vt,
        coef);
    REQUIRE(std::isfinite(electronFlux));
    REQUIRE(electronFlux == Approx(0.0));

    const double holeFlux = sgHoleContinuityFluxFromQuasiFermi(
        ni,
        -50.0,
        0.0,
        0.0,
        0.0,
        Vt,
        coef);
    REQUIRE(std::isfinite(holeFlux));
    REQUIRE(holeFlux == Approx(0.0));
}

TEST_CASE("SG quasi-Fermi fluxes cancel flat QF at large absolute bias with electric field",
          "[sg][coupled]")
{
    const double Vt = constants::Vt_300;
    const double ni = 1.0e16;
    const double coef = 1.0;
    const double dpsi = 1.0;
    const double psi0 = -20.0;
    const double psi1 = psi0 + dpsi;
    const double phin = -12.8;
    const double phip = -12.8;

    const double electronFlux = sgElectronContinuityFluxFromQuasiFermi(
        ni,
        psi1,
        phin,
        phin,
        dpsi,
        Vt,
        coef);
    const double holeFlux = sgHoleContinuityFluxFromQuasiFermi(
        ni,
        psi0,
        phip,
        phip,
        dpsi,
        Vt,
        coef);

    REQUIRE(std::isfinite(electronFlux));
    REQUIRE(std::isfinite(holeFlux));
    REQUIRE(electronFlux == Approx(0.0).margin(1.0e-30));
    REQUIRE(holeFlux == Approx(0.0).margin(1.0e-30));
}

TEST_CASE("SG quasi-Fermi fluxes cancel flat QF with variable intrinsic density", "[sg][coupled][bgn]")
{
    const double Vt = constants::Vt_300;
    const double coef = 1.0;
    const double ni0 = 1.0e16;
    const double ni1 = 1.5e16;
    const double psi0 = 0.013;
    const double psi1 = -0.021;
    const double dpsi = psi1 - psi0;
    const double phin = 0.004;
    const double phip = -0.003;

    const double electronFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
        ni0,
        ni1,
        psi0,
        psi1,
        phin,
        phin,
        Vt,
        coef);
    const double holeFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
        ni0,
        ni1,
        psi0,
        psi1,
        phip,
        phip,
        Vt,
        coef);

    REQUIRE(electronFlux == Approx(0.0).margin(1.0e-30));
    REQUIRE(holeFlux == Approx(0.0).margin(1.0e-30));
}

TEST_CASE("SG variable-ni quasi-Fermi flux matches density form at large absolute bias",
          "[sg][coupled][bgn]")
{
    const double Vt = constants::Vt_300;
    const double coef = 1.0;
    const double ni = 1.65563e16;

    const double psi0 = -13.203650871693659;
    const double psi1 = -13.203650871693661;
    const double phin0 = -12.79890541782345;
    const double phin1 = -12.800000000000001;
    const double phip0 = -12.799999999999999;
    const double phip1 = -12.800000000000001;

    const double n0 = ni * std::exp((psi0 - phin0) / Vt);
    const double n1 = ni * std::exp((psi1 - phin1) / Vt);
    const double p0 = ni * std::exp((phip0 - psi0) / Vt);
    const double p1 = ni * std::exp((phip1 - psi1) / Vt);

    const double electronDensityFlux =
        sgElectronContinuityFlux(n0, n1, psi1 - psi0, Vt, coef);
    const double electronPlainQfFlux = sgElectronContinuityFluxFromQuasiFermi(
        ni,
        psi1,
        phin0,
        phin1,
        psi1 - psi0,
        Vt,
        coef);
    const double electronQfFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
        ni,
        ni,
        psi0,
        psi1,
        phin0,
        phin1,
        Vt,
        coef);
    const double holeDensityFlux =
        sgHoleContinuityFlux(p0, p1, psi1 - psi0, Vt, coef);
    const double holePlainQfFlux = sgHoleContinuityFluxFromQuasiFermi(
        ni,
        psi0,
        phip0,
        phip1,
        psi1 - psi0,
        Vt,
        coef);
    const double holeQfFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
        ni,
        ni,
        psi0,
        psi1,
        phip0,
        phip1,
        Vt,
        coef);

    REQUIRE(std::isfinite(electronQfFlux));
    REQUIRE(std::isfinite(holeQfFlux));
    REQUIRE(electronPlainQfFlux == Approx(electronDensityFlux).epsilon(1.0e-12));
    REQUIRE(electronQfFlux == Approx(electronDensityFlux).epsilon(1.0e-12));
    REQUIRE(holePlainQfFlux == Approx(holeDensityFlux).epsilon(1.0e-12));
    const long double holeEta =
        (static_cast<long double>(psi1) - static_cast<long double>(psi0))
        / static_cast<long double>(Vt);
    const long double holeRight =
        static_cast<long double>(bernoulli(static_cast<Real>(-holeEta)))
        * static_cast<long double>(ni)
        * std::exp(
            (static_cast<long double>(phip1) - static_cast<long double>(psi1))
            / static_cast<long double>(Vt));
    const long double holeExpected =
        static_cast<long double>(coef) * holeRight
        * std::expm1(
            (static_cast<long double>(phip0) - static_cast<long double>(phip1))
            / static_cast<long double>(Vt));
    REQUIRE(holeQfFlux ==
            Approx(static_cast<Real>(holeExpected)).epsilon(1.0e-12));
}

TEST_CASE("SG variable-ni electron flux decomposition reconstructs the production flux",
          "[sg][diagnostics][bgn]")
{
    const Real Vt = constants::Vt_300;
    const Real ni0 = 1.0e16;
    const Real ni1 = 1.7e16;
    const Real psi0 = -0.14;
    const Real psi1 = 0.23;
    const Real phin0 = -0.031;
    const Real phin1 = 0.047;
    const Real coef = 3.25e-4;

    for (const bool includeNiGradientDrift : {false, true}) {
        const SgElectronVariableNiFluxDecomposition decomposition =
            sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
                ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef,
                includeNiGradientDrift);
        const Real production = sgElectronContinuityFluxFromQuasiFermiVariableNi(
            ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef,
            includeNiGradientDrift);

        REQUIRE(decomposition.ni0 == ni0);
        REQUIRE(decomposition.ni1 == ni1);
        REQUIRE(decomposition.psi0 == psi0);
        REQUIRE(decomposition.psi1 == psi1);
        REQUIRE(decomposition.phin0 == phin0);
        REQUIRE(decomposition.phin1 == phin1);
        REQUIRE(decomposition.coef == coef);
        REQUIRE(decomposition.includeNiGradientDrift == includeNiGradientDrift);
        REQUIRE(decomposition.eta == Approx(
            (psi1 - psi0) / Vt
            + (includeNiGradientDrift ? std::log(ni1 / ni0) : 0.0)));
        REQUIRE(decomposition.leftTerm ==
                Approx(decomposition.bernoulliMinusEta * decomposition.n0));
        REQUIRE(decomposition.rightTerm ==
                Approx(decomposition.bernoulliEta * decomposition.n1));
        REQUIRE(decomposition.signedDifference ==
                Approx(decomposition.leftTerm - decomposition.rightTerm));
        REQUIRE(decomposition.stableFactorizedFlux == production);
        REQUIRE(std::isfinite(decomposition.stableFactorizedFlux));
        REQUIRE(std::isfinite(decomposition.highPrecisionReferenceFlux));
        REQUIRE(std::isfinite(decomposition.highPrecisionReferenceTermScale));
        const Real referenceScale = std::max({
            std::abs(production),
            decomposition.highPrecisionReferenceTermScale,
            Real{1.0e-300}});
        REQUIRE(std::abs(production - decomposition.highPrecisionReferenceFlux)
                / referenceScale <= 1.0e-6);
        REQUIRE(std::isfinite(decomposition.cancellationCondition));
        REQUIRE_FALSE(decomposition.node0ExponentClampedLow);
        REQUIRE_FALSE(decomposition.node0ExponentClampedHigh);
        REQUIRE_FALSE(decomposition.node1ExponentClampedLow);
        REQUIRE_FALSE(decomposition.node1ExponentClampedHigh);
    }
}

TEST_CASE("SG variable-ni electron flux decomposition is oriented and handles large eta",
          "[sg][diagnostics][bgn]")
{
    const Real Vt = constants::Vt_300;
    const Real ni0 = 8.0e15;
    const Real ni1 = 2.0e16;
    const Real psi0 = 0.0;
    const Real psi1 = 2.5;
    const Real phin0 = 0.10;
    const Real phin1 = 2.40;
    const Real coef = 1.0e-6;

    const auto forward = sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
        ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef, true);
    const auto reverse = sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
        ni1, ni0, psi1, psi0, phin1, phin0, Vt, coef, true);

    REQUIRE(forward.eta > 50.0);
    REQUIRE(reverse.eta < -50.0);
    REQUIRE(forward.stableFactorizedFlux ==
            sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef, true));
    REQUIRE(reverse.stableFactorizedFlux ==
            sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni1, ni0, psi1, psi0, phin1, phin0, Vt, coef, true));
    REQUIRE(reverse.reconstructedFlux ==
            Approx(-forward.reconstructedFlux).epsilon(1.0e-12));
    REQUIRE(reverse.stableFactorizedFlux ==
            Approx(-forward.stableFactorizedFlux).epsilon(1.0e-12));
}

TEST_CASE("SG variable-ni electron flux decomposition reports finite exact-cancellation policy",
          "[sg][diagnostics][bgn]")
{
    const Real Vt = constants::Vt_300;
    const Real phin = 0.004;
    const auto decomposition =
        sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
            1.0e16, 1.5e16, 0.013, -0.021, phin, phin,
            Vt, 1.0, true);

    REQUIRE(decomposition.flatQuasiFermiShortCircuit);
    REQUIRE(decomposition.reconstructedFlux == 0.0);
    REQUIRE(decomposition.stableFactorizedFlux == 0.0);
    REQUIRE(decomposition.highPrecisionReferenceFlux == 0.0);
    REQUIRE(std::isfinite(decomposition.highPrecisionReferenceTermScale));
    REQUIRE(std::isfinite(decomposition.cancellationCondition));
    REQUIRE(decomposition.cancellationCondition == std::numeric_limits<Real>::max());
}

TEST_CASE("SG variable-ni electron flux decomposition remains stable under severe cancellation",
          "[sg][diagnostics][bgn]")
{
    const Real Vt = constants::Vt_300;
    const Real ni0 = 1.0e16;
    const Real ni1 = 1.4e16;
    const Real psi0 = 9.99;
    const Real psi1 = 10.0;
    const Real phin0 = 0.0;
    const Real phin1 = 1.0e-12;
    const Real coef = 1.0e-170;
    const auto decomposition =
        sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
            ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef, true);

    const long double eta =
        (static_cast<long double>(psi1) - static_cast<long double>(psi0))
            / static_cast<long double>(Vt)
        + std::log(static_cast<long double>(ni1) / static_cast<long double>(ni0));
    const long double bernoulliEta = eta / std::expm1(eta);
    const long double expected = static_cast<long double>(coef) * bernoulliEta
        * static_cast<long double>(ni1)
        * std::exp(static_cast<long double>(psi1) / static_cast<long double>(Vt))
        * (std::exp(-static_cast<long double>(phin0) / static_cast<long double>(Vt))
           - std::exp(-static_cast<long double>(phin1) / static_cast<long double>(Vt)));

    REQUIRE(decomposition.cancellationCondition > 1.0e9);
    REQUIRE(std::isfinite(decomposition.stableFactorizedFlux));
    REQUIRE(std::isfinite(decomposition.highPrecisionReferenceFlux));
    REQUIRE(decomposition.highPrecisionReferenceFlux ==
            Approx(static_cast<Real>(expected)).epsilon(1.0e-8));
    const Real productionError = std::abs(
        decomposition.reconstructedFlux - decomposition.highPrecisionReferenceFlux);
    const Real stableError = std::abs(
        decomposition.stableFactorizedFlux - decomposition.highPrecisionReferenceFlux);
    REQUIRE(stableError <= productionError);
    REQUIRE(stableError <= 0.1 * productionError);
    REQUIRE(decomposition.stableFactorizedFlux ==
            sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni0, ni1, psi0, psi1, phin0, phin1, Vt, coef, true));
}

TEST_CASE("SG variable-ni electron flux decomposition reports endpoint exponent clamps",
          "[sg][diagnostics][bgn]")
{
    const Real Vt = constants::Vt_300;
    const auto decomposition =
        sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
            1.0e16, 1.0e16, 13.0, -13.0, 0.0, 0.1,
            Vt, 1.0, true);

    REQUIRE(decomposition.node0ExponentClampedHigh);
    REQUIRE_FALSE(decomposition.node0ExponentClampedLow);
    REQUIRE(decomposition.node1ExponentClampedLow);
    REQUIRE_FALSE(decomposition.node1ExponentClampedHigh);
    REQUIRE(std::isfinite(decomposition.reconstructedFlux));
    REQUIRE(std::isfinite(decomposition.stableFactorizedFlux));
    REQUIRE(std::isfinite(decomposition.highPrecisionReferenceFlux));
    REQUIRE(std::isfinite(decomposition.highPrecisionReferenceTermScale));
    REQUIRE(decomposition.stableFactorizedFlux ==
            sgElectronContinuityFluxFromQuasiFermiVariableNi(
                1.0e16, 1.0e16, 13.0, -13.0, 0.0, 0.1,
                Vt, 1.0, true));
}


static DeviceMesh makeSingleSiliconTriangleMesh()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6;  n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.25e-6; n2.y = 0.8e-6;  mesh.addNode(n2);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Region r0; r0.id = 0; r0.name = "silicon"; r0.material = "Si"; r0.cell_ids = {0};
    mesh.addRegion(r0);

    mesh.buildEdges();
    return mesh;
}

TEST_CASE("mobility doping bases distinguish net, total, and cell reconstruction",
          "[sg][mobility][doping-basis]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e23, 1.0e23);
    doping.setNodeDoping(1, 1.0e23, 0.0);
    doping.setNodeDoping(2, 1.0e23, 0.0);

    MobilityModelConfig netConfig;
    MobilityModelConfig totalConfig;
    totalConfig.dopingConcentrationBasis = "total_impurity";
    MobilityModelConfig cellConfig;
    cellConfig.dopingConcentrationBasis =
        "cell_reconstructed_total_impurity";

    REQUIRE(detail::nodeMobilityDopingConcentration(
                mesh, doping, 0, 0, &netConfig) == Approx(0.0));
    REQUIRE(detail::nodeMobilityDopingConcentration(
                mesh, doping, 0, 0, &totalConfig) == Approx(2.0e23));
    REQUIRE(detail::nodeMobilityDopingConcentration(
                mesh, doping, 0, 0, &cellConfig) == Approx(4.0e23 / 3.0));

    const auto edgeIt = std::find_if(
        mesh.edges().begin(), mesh.edges().end(), [](const Edge& edge) {
            return edge.n0 == 0 && edge.n1 == 1;
        });
    REQUIRE(edgeIt != mesh.edges().end());
    REQUIRE(detail::edgeMobilityDopingConcentration(
                mesh, doping, *edgeIt, 0, &netConfig) == Approx(0.5e23));
    REQUIRE(detail::edgeMobilityDopingConcentration(
                mesh, doping, *edgeIt, 0, &totalConfig) == Approx(1.5e23));
    REQUIRE(detail::edgeMobilityDopingConcentration(
                mesh, doping, *edgeIt, 0, &cellConfig) == Approx(4.0e23 / 3.0));
}
TEST_CASE("CoupledDDAssembler BGN continuity residuals vanish for flat quasi-Fermi levels",
          "[sg][coupled][bgn]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e24, 0.0);
    doping.setNodeDoping(1, 0.0, 1.0e23);
    doping.setNodeDoping(2, 1.0e24, 1.0e24);

    BandgapNarrowingConfig bgn;
    bgn.model = "slotboom";
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               recombinationModelConfig({"none"}),
                               bgn);

    REQUIRE(coupled.intrinsicDensity()[0] != Approx(coupled.intrinsicDensity()[1]));

    CoupledDDState state;
    state.psi.resize(3);
    state.phin.resize(3);
    state.phip.resize(3);
    state.psi << 0.020, -0.010, 0.030;
    state.phin << 0.0, 0.0, 0.0;
    state.phip << 0.0, 0.0, 0.0;

    const VectorXd residual = coupled.residual(coupled.pack(state), CoupledDDBoundaryConditions{});
    const int N = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < N; ++i) {
        REQUIRE(residual(N + i) == Approx(0.0).margin(1.0e-18));
        REQUIRE(residual(2 * N + i) == Approx(0.0).margin(1.0e-18));
    }
}

TEST_CASE("SG continuity residuals match DDAssembler and CoupledDDAssembler", "[sg][dd][coupled]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});

    DDAssembler dd(mesh,
                   matdb,
                   doping,
                   constants::Vt_300,
                   MobilityModelConfig{},
                   noRecombination);
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               noRecombination);

    CoupledDDState state;
    state.psi.resize(3);
    state.phin.resize(3);
    state.phip.resize(3);
    state.psi << 0.020, -0.010, 0.030;
    state.phin << 0.005, -0.002, 0.010;
    state.phip << -0.004, 0.006, -0.008;

    const VectorXd x = coupled.pack(state);
    const VectorXd n = coupled.electronDensity(x);
    const VectorXd p = coupled.holeDensity(x);
    const VectorXd coupledResidual = coupled.residual(x, CoupledDDBoundaryConditions{});

    dd.assembleElectronContinuity(state.psi, n, p);
    const VectorXd ddElectronResidual = dd.matrix() * n - dd.rhs();

    dd.assembleHoleContinuity(state.psi, n, p);
    const VectorXd ddHoleResidual = dd.matrix() * p - dd.rhs();

    const int N = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < N; ++i) {
        const double electronScale = std::max(1.0, std::abs(ddElectronResidual(i)));
        const double holeScale = std::max(1.0, std::abs(ddHoleResidual(i)));
        REQUIRE(coupledResidual(N + i) / electronScale ==
                Approx(ddElectronResidual(i) / electronScale).epsilon(1.0e-12).margin(1.0e-12));
        REQUIRE(coupledResidual(2 * N + i) / holeScale ==
                Approx(ddHoleResidual(i) / holeScale).epsilon(1.0e-12).margin(1.0e-12));
    }
}

TEST_CASE("CoupledDDAssembler validates doping size before BGN construction", "[sg][dd][coupled][bgn]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel shortDoping(mesh.numNodes() - 1);
    BandgapNarrowingConfig bgn;
    bgn.model = "slotboom";

    REQUIRE_THROWS_AS(DDAssembler(mesh,
                                  matdb,
                                  shortDoping,
                                  constants::Vt_300,
                                  MobilityModelConfig{},
                                  recombinationModelConfig({"none"}),
                                  bgn),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(CoupledDDAssembler(mesh,
                                         matdb,
                                         shortDoping,
                                         constants::Vt_300,
                                         MobilityModelConfig{},
                                         recombinationModelConfig({"none"}),
                                         bgn),
                      std::invalid_argument);
}

TEST_CASE("DDAssembler rejects incomplete unit scaling references", "[sg][dd][scaling]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0e23;
    scaling.mu0 = 0.135;
    scaling.D0 = scaling.mu0 * scaling.V0;
    scaling.L0 = 0.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;

    REQUIRE_THROWS_AS(DDAssembler(mesh,
                                  matdb,
                                  doping,
                                  constants::Vt_300,
                                  MobilityModelConfig{},
                                  recombinationModelConfig({"none"}),
                                  BandgapNarrowingConfig{},
                                  ImpactIonizationModelConfig{},
                                  {},
                                  {},
                                  scaling),
                      std::invalid_argument);

    scaling.L0 = std::numeric_limits<Real>::quiet_NaN();
    REQUIRE_THROWS_AS(DDAssembler(mesh,
                                  matdb,
                                  doping,
                                  constants::Vt_300,
                                  MobilityModelConfig{},
                                  recombinationModelConfig({"none"}),
                                  BandgapNarrowingConfig{},
                                  ImpactIonizationModelConfig{},
                                  {},
                                  {},
                                  scaling),
                      std::invalid_argument);
}

static DeviceMesh makeSingleSiliconTriangleMeshMicrometers()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0;  n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0;  n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.25; n2.y = 0.8; mesh.addNode(n2);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Region r0; r0.id = 0; r0.name = "silicon"; r0.material = "Si"; r0.cell_ids = {0};
    mesh.addRegion(r0);

    mesh.buildEdges();
    return mesh;
}

TEST_CASE("CoupledDDAssembler unit scaling Poisson residual matches DDAssembler",
          "[sg][dd][coupled][scaling]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMeshMicrometers();
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};
    MaterialDatabase matdb(unitScaling);
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e17, 0.0);
    doping.setNodeDoping(1, 0.0, 7.5e16);
    doping.setNodeDoping(2, 2.0e16, 0.0);

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0e17;
    scaling.mu0 = 1350.0;
    scaling.D0 = scaling.mu0 * scaling.V0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.unitSystem = unitScaling.unitSystem();
    scaling.chargeAreaFactor = unitScaling.unitSystem().chargeAreaFactor();
    scaling.chargeLineFactor = unitScaling.unitSystem().chargeLineFactor();
    scaling.fieldFromCoordinateDeltaFactor =
        unitScaling.unitSystem().fieldFromCoordinateDeltaFactor();

    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});
    DDAssembler dd(mesh,
                   matdb,
                   doping,
                   constants::Vt_300,
                   MobilityModelConfig{},
                   noRecombination,
                   BandgapNarrowingConfig{},
                   ImpactIonizationModelConfig{},
                   {},
                   {},
                   scaling);
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               noRecombination,
                               BandgapNarrowingConfig{},
                               ImpactIonizationModelConfig{},
                               {},
                               {},
                               scaling);

    CoupledDDState state;
    state.psi.resize(3);
    state.phin.resize(3);
    state.phip.resize(3);
    state.psi << 0.020, -0.010, 0.015;
    state.phin << 0.001, -0.002, 0.0005;
    state.phip << -0.0015, 0.0025, -0.00025;

    VectorXd x = coupled.pack(state);
    x.segment(0, 3) /= scaling.V0;
    x.segment(3, 3) /= scaling.V0;
    x.segment(6, 3) /= scaling.V0;
    const VectorXd n = coupled.electronDensity(x);
    const VectorXd p = coupled.holeDensity(x);
    const VectorXd nScaled = n / scaling.C0;
    const VectorXd pScaled = p / scaling.C0;
    const VectorXd psiScaled = state.psi / scaling.V0;

    dd.assemblePoissonWithCarriers(nScaled, pScaled, psiScaled);
    const VectorXd ddResidual = dd.matrix() * psiScaled - dd.rhs();
    const VectorXd coupledResidual = coupled.residual(x, CoupledDDBoundaryConditions{});

    for (int row = 0; row < 3; ++row)
        REQUIRE(coupledResidual(row) == Approx(ddResidual(row)).epsilon(1.0e-12).margin(1.0e-12));
}
TEST_CASE("CoupledDDAssembler unit scaling continuity residual matches DDAssembler",
          "[sg][dd][coupled][scaling]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMeshMicrometers();
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};
    MaterialDatabase matdb(unitScaling);
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e17, 0.0);
    doping.setNodeDoping(1, 0.0, 7.5e16);
    doping.setNodeDoping(2, 2.0e16, 0.0);

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0e17;
    scaling.mu0 = 1350.0;
    scaling.D0 = scaling.mu0 * scaling.V0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.unitSystem = unitScaling.unitSystem();
    scaling.chargeAreaFactor = unitScaling.unitSystem().chargeAreaFactor();
    scaling.chargeLineFactor = unitScaling.unitSystem().chargeLineFactor();
    scaling.fieldFromCoordinateDeltaFactor =
        unitScaling.unitSystem().fieldFromCoordinateDeltaFactor();
    scaling.currentDensityLineIntegralFactor =
        unitScaling.unitSystem().currentDensityAM2PerInternal() *
        unitScaling.unitSystem().lengthMPerInternal();

    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});
    DDAssembler dd(mesh,
                   matdb,
                   doping,
                   constants::Vt_300,
                   MobilityModelConfig{},
                   noRecombination,
                   BandgapNarrowingConfig{},
                   ImpactIonizationModelConfig{},
                   {},
                   {},
                   scaling);
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               noRecombination,
                               BandgapNarrowingConfig{},
                               ImpactIonizationModelConfig{},
                               {},
                               {},
                               scaling);

    CoupledDDState state;
    state.psi.resize(3);
    state.phin.resize(3);
    state.phip.resize(3);
    state.psi << 0.020, -0.010, 0.015;
    state.phin << 0.001, -0.002, 0.0005;
    state.phip << -0.0015, 0.0025, -0.00025;

    VectorXd x = coupled.pack(state);
    x.segment(0, 3) /= scaling.V0;
    x.segment(3, 3) /= scaling.V0;
    x.segment(6, 3) /= scaling.V0;
    const VectorXd n = coupled.electronDensity(x);
    const VectorXd p = coupled.holeDensity(x);
    const VectorXd nScaled = n / scaling.C0;
    const VectorXd pScaled = p / scaling.C0;
    const VectorXd psiScaled = state.psi / scaling.V0;

    dd.assembleElectronContinuity(psiScaled, nScaled, pScaled);
    const VectorXd ddElectronResidual = dd.matrix() * nScaled - dd.rhs();
    dd.assembleHoleContinuity(psiScaled, nScaled, pScaled);
    const VectorXd ddHoleResidual = dd.matrix() * pScaled - dd.rhs();
    const VectorXd coupledResidual = coupled.residual(x, CoupledDDBoundaryConditions{});

    for (int row = 0; row < 3; ++row) {
        REQUIRE(coupledResidual(3 + row) ==
                Approx(ddElectronResidual(row)).epsilon(1.0e-12).margin(1.0e-12));
        REQUIRE(coupledResidual(6 + row) ==
                Approx(ddHoleResidual(row)).epsilon(1.0e-12).margin(1.0e-12));
    }
}

static DeviceMesh makeContactedSiliconSquareMesh(Real sideLength)
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0;        n0.y = 0.0;        mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = sideLength; n1.y = 0.0;        mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = sideLength; n2.y = sideLength; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0;        n3.y = sideLength; mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region r0; r0.id = 0; r0.name = "silicon"; r0.material = "Si"; r0.cell_ids = {0, 1};
    mesh.addRegion(r0);

    Contact left; left.id = 0; left.name = "left"; left.region_id = 0; left.node_ids = {0, 3};
    mesh.addContact(left);
    Contact right; right.id = 1; right.name = "right"; right.region_id = 0; right.node_ids = {1, 2};
    mesh.addContact(right);

    mesh.buildEdges();
    return mesh;
}

static DDScalingSpec makeTcadCurrentScalingSpec(const UnitScalingConfig& unitScaling)
{
    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0e17;
    scaling.mu0 = 1350.0;
    scaling.D0 = scaling.mu0 * scaling.V0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.unitSystem = unitScaling.unitSystem();
    scaling.chargeAreaFactor = unitScaling.unitSystem().chargeAreaFactor();
    scaling.chargeLineFactor = unitScaling.unitSystem().chargeLineFactor();
    scaling.fieldFromCoordinateDeltaFactor =
        unitScaling.unitSystem().fieldFromCoordinateDeltaFactor();
    scaling.currentDensityLineIntegralFactor =
        unitScaling.unitSystem().currentDensityAM2PerInternal() *
        unitScaling.unitSystem().lengthMPerInternal();
    return scaling;
}

TEST_CASE("ContactCurrent unit scaling terminal current matches legacy SI",
          "[sg][contact_current][scaling]")
{
    DeviceMesh legacyMesh = makeContactedSiliconSquareMesh(1.0e-6);
    DeviceMesh scaledMesh = makeContactedSiliconSquareMesh(1.0);
    const UnitScalingConfig unitScaling{UnitScalingMode::UnitScaling};
    MaterialDatabase legacyMatdb;
    MaterialDatabase scaledMatdb(unitScaling);

    DopingModel legacyDoping(legacyMesh.numNodes());
    DopingModel scaledDoping(scaledMesh.numNodes());
    for (Index node = 0; node < legacyMesh.numNodes(); ++node) {
        legacyDoping.setNodeDoping(node, 1.0e23, 5.0e21);
        scaledDoping.setNodeDoping(node, 1.0e17, 5.0e15);
    }

    DDSolution legacy;
    DDSolution scaled;
    legacy.psi.resize(4); legacy.phin.resize(4); legacy.phip.resize(4);
    scaled.psi.resize(4); scaled.phin.resize(4); scaled.phip.resize(4);
    legacy.psi << 0.000, 0.045, 0.038, -0.006;
    legacy.phin << -0.010, 0.006, 0.004, -0.012;
    legacy.phip << 0.008, -0.004, -0.006, 0.010;
    scaled.psi = legacy.psi;
    scaled.phin = legacy.phin;
    scaled.phip = legacy.phip;

    const Real legacyNi = legacyMatdb.getMaterial("Si").ni;
    const Real scaledNi = scaledMatdb.getMaterial("Si").ni;
    legacy.n.resize(4); legacy.p.resize(4); scaled.n.resize(4); scaled.p.resize(4);
    for (int i = 0; i < 4; ++i) {
        legacy.n(i) = legacyNi * std::exp((legacy.psi(i) - legacy.phin(i)) / constants::Vt_300);
        legacy.p(i) = legacyNi * std::exp((legacy.phip(i) - legacy.psi(i)) / constants::Vt_300);
        scaled.n(i) = scaledNi * std::exp((scaled.psi(i) - scaled.phin(i)) / constants::Vt_300);
        scaled.p(i) = scaledNi * std::exp((scaled.phip(i) - scaled.psi(i)) / constants::Vt_300);
    }

    const MobilityModelConfig mobility = mobilityModelConfig("constant");
    ContactCurrent legacyCurrent(legacyMesh, legacyMatdb, legacyDoping, mobility, constants::T0);
    ContactCurrent scaledCurrent(scaledMesh,
                                 scaledMatdb,
                                 scaledDoping,
                                 mobility,
                                 constants::T0,
                                 makeTcadCurrentScalingSpec(unitScaling));

    const ContactCurrentDetailedResult scaledDetailed = scaledCurrent.computeDetailed(scaled, "left");
    REQUIRE(!scaledDetailed.edges.empty());
    for (const ContactCurrentEdgeDiagnostic& edge : scaledDetailed.edges) {
        REQUIRE(edge.edgeLength_m > 0.0);
        REQUIRE(edge.edgeLength_m < 2.0e-6);
        REQUIRE(edge.edgeCouple_m > 0.0);
        REQUIRE(edge.edgeCouple_m < 2.0e-6);
    }

    const ContactCurrentResult legacyLeft = legacyCurrent.compute(legacy, "left");
    const ContactCurrentResult scaledLeft = scaledDetailed.totals;

    const Real scale = std::max(1.0, std::abs(legacyLeft.totalCurrent));
    REQUIRE(scaledLeft.electronCurrent / scale ==
            Approx(legacyLeft.electronCurrent / scale).epsilon(1.0e-12).margin(1.0e-12));
    REQUIRE(scaledLeft.holeCurrent / scale ==
            Approx(legacyLeft.holeCurrent / scale).epsilon(1.0e-12).margin(1.0e-12));
    REQUIRE(scaledLeft.totalCurrent / scale ==
            Approx(legacyLeft.totalCurrent / scale).epsilon(1.0e-12).margin(1.0e-12));
}
TEST_CASE("Slotboom BGN uses total impurity density for compensated nodes", "[sg][coupled][bgn]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        doping.setNodeDoping(i, 1.0e24, 1.0e24);

    BandgapNarrowingConfig bgn;
    bgn.model = "slotboom";
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               recombinationModelConfig({"none"}),
                               bgn);

    const Material& si = matdb.getMaterial("Si");
    for (Real niEff : coupled.intrinsicDensity())
        REQUIRE(niEff > si.ni);
}

TEST_CASE("CoupledDDAssembler BGN residuals use variable-ni quasi-Fermi fluxes",
          "[sg][dd][coupled][bgn]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e24, 0.0);
    doping.setNodeDoping(1, 0.0, 1.0e23);
    doping.setNodeDoping(2, 1.0e24, 1.0e24);

    BandgapNarrowingConfig bgn;
    bgn.model = "slotboom";
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});

    DDAssembler dd(mesh,
                   matdb,
                   doping,
                   constants::Vt_300,
                   MobilityModelConfig{},
                   noRecombination,
                   bgn);
    CoupledDDAssembler coupled(mesh,
                               matdb,
                               doping,
                               constants::Vt_300,
                               MobilityModelConfig{},
                               noRecombination,
                               bgn);

    REQUIRE(coupled.intrinsicDensity()[0] != Approx(coupled.intrinsicDensity()[1]));
    REQUIRE(coupled.intrinsicDensity()[2] > coupled.intrinsicDensity()[0]);

    CoupledDDState state;
    state.psi.resize(3);
    state.phin.resize(3);
    state.phip.resize(3);
    state.psi << 0.020, -0.010, 0.030;
    state.phin << 0.005, -0.002, 0.010;
    state.phip << -0.004, 0.006, -0.008;

    const VectorXd x = coupled.pack(state);
    const VectorXd n = coupled.electronDensity(x);
    const VectorXd p = coupled.holeDensity(x);
    const VectorXd coupledResidual = coupled.residual(x, CoupledDDBoundaryConditions{});

    dd.assembleElectronContinuity(state.psi, n, p);
    const VectorXd ddElectronResidual = dd.matrix() * n - dd.rhs();

    dd.assembleHoleContinuity(state.psi, n, p);
    const VectorXd ddHoleResidual = dd.matrix() * p - dd.rhs();

    bool sawElectronDifference = false;
    bool sawHoleDifference = false;
    const int N = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < N; ++i) {
        REQUIRE(std::isfinite(coupledResidual(N + i)));
        REQUIRE(std::isfinite(coupledResidual(2 * N + i)));
        const double electronScale = std::max({1.0, std::abs(coupledResidual(N + i)),
                                               std::abs(ddElectronResidual(i))});
        const double holeScale = std::max({1.0, std::abs(coupledResidual(2 * N + i)),
                                           std::abs(ddHoleResidual(i))});
        sawElectronDifference = sawElectronDifference ||
            std::abs(coupledResidual(N + i) - ddElectronResidual(i)) / electronScale > 1.0e-6;
        sawHoleDifference = sawHoleDifference ||
            std::abs(coupledResidual(2 * N + i) - ddHoleResidual(i)) / holeScale > 1.0e-6;
    }
    REQUIRE(sawElectronDifference);
    REQUIRE(sawHoleDifference);
}

struct AssemblySystem {
    SparseMatrixd A;
    VectorXd b;
};

static AssemblySystem assembleReferencePoissonWithFreshGeometry(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    const VectorXd& n,
    const VectorXd& p,
    const VectorXd& psi)
{
    const Index N = mesh.numNodes();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto vol = detail::computeNodeVolumes(mesh);
    const auto couple = detail::computeEdgeCouplings(mesh);

    AssemblySystem system{SparseMatrixd(static_cast<int>(N), static_cast<int>(N)),
                          VectorXd::Zero(static_cast<int>(N))};

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.numEdges() * 4 + N);

    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const Real eps = detail::edgeEpsilon(edgeCells, mesh, matdb, e);
        const Real G = eps * couple[e] / h;

        const auto i = static_cast<int>(edge.n0);
        const auto j = static_cast<int>(edge.n1);
        triplets.emplace_back(i, i, G);
        triplets.emplace_back(j, j, G);
        triplets.emplace_back(i, j, -G);
        triplets.emplace_back(j, i, -G);
    }

    system.A.setFromTriplets(triplets.begin(), triplets.end());

    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        const Real diagCarrier = constants::q * (n(ii) + p(ii)) / Vt * vol[i];
        system.A.coeffRef(ii, ii) += diagCarrier;
        system.b(ii) = constants::q * (p(ii) - n(ii) + doping.netDoping(i)) * vol[i]
                       + diagCarrier * psi(ii);
    }

    return system;
}

static AssemblySystem assembleReferenceContinuityWithFreshGeometry(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    const MobilityModelConfig& mobilityConfig,
    const RecombinationModelConfig& recombinationConfig,
    CarrierType carrier,
    const VectorXd& psi,
    const VectorXd& nOld,
    const VectorXd& pOld)
{
    const Index N = mesh.numNodes();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto vol = detail::computeNodeVolumes(mesh);
    const auto couple = detail::computeEdgeCouplings(mesh);
    const Real temperature_K = Vt * constants::q / constants::kb;
    const auto ni = detail::buildNodeNi(mesh, matdb, temperature_K);
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);
    const auto mobility = makeMobilityModel(mobilityConfig);
    const RecombinationModel recombination(recombinationConfig);

    AssemblySystem system{SparseMatrixd(static_cast<int>(N), static_cast<int>(N)),
                          VectorXd::Zero(static_cast<int>(N))};

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.numEdges() * 4 + N);

    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const auto i = static_cast<int>(edge.n0);
        const auto j = static_cast<int>(edge.n1);
        const Real electricField = std::abs((psi(j) - psi(i)) / h);
        const Real mu = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, e, carrier, electricField);
        if (mu <= 0.0) continue;

        const Real coef = mu * Vt * couple[e] / h;
        const Real dpsi = psi(j) - psi(i);
        const SGEdgeWeights weights = sgEdgeWeights(dpsi, Vt);

        if (carrier == CarrierType::Electron) {
            triplets.emplace_back(i, i, coef * weights.b_minus);
            triplets.emplace_back(j, j, coef * weights.b_plus);
            triplets.emplace_back(i, j, -coef * weights.b_plus);
            triplets.emplace_back(j, i, -coef * weights.b_minus);
        } else {
            triplets.emplace_back(i, i, coef * weights.b_plus);
            triplets.emplace_back(j, j, coef * weights.b_minus);
            triplets.emplace_back(i, j, -coef * weights.b_minus);
            triplets.emplace_back(j, i, -coef * weights.b_plus);
        }
    }

    system.A.setFromTriplets(triplets.begin(), triplets.end());

    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        const RecombinationLinearization linearization =
            carrier == CarrierType::Electron
                ? recombination.electronLinearization(nOld(ii), pOld(ii), ni[i])
                : recombination.holeLinearization(nOld(ii), pOld(ii), ni[i]);
        system.A.coeffRef(ii, ii) += linearization.diagonal * vol[i];
        system.b(ii) += linearization.rhs * vol[i];
    }

    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        if (system.A.coeff(ii, ii) == 0.0) {
            system.A.coeffRef(ii, ii) = 1.0;
            system.b(ii) = 0.0;
        }
    }

    return system;
}

static void requireSystemsMatch(const SparseMatrixd& lhsA,
                                const VectorXd& lhsB,
                                const SparseMatrixd& rhsA,
                                const VectorXd& rhsB)
{
    REQUIRE(lhsA.rows() == rhsA.rows());
    REQUIRE(lhsA.cols() == rhsA.cols());
    REQUIRE(lhsB.size() == rhsB.size());

    for (int row = 0; row < lhsA.rows(); ++row) {
        for (int col = 0; col < lhsA.cols(); ++col) {
            const double lhs = lhsA.coeff(row, col);
            const double rhs = rhsA.coeff(row, col);
            const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
            REQUIRE(lhs / scale == Approx(rhs / scale).epsilon(1.0e-14).margin(1.0e-14));
        }
    }

    for (int row = 0; row < lhsB.size(); ++row) {
        const double lhs = lhsB(row);
        const double rhs = rhsB(row);
        const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
        REQUIRE(lhs / scale == Approx(rhs / scale).epsilon(1.0e-14).margin(1.0e-14));
    }
}

TEST_CASE("DDAssembler cached geometry matches fresh reference assembly", "[sg][dd][cache]")
{
    DeviceMesh mesh = makeSingleSiliconTriangleMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    const MobilityModelConfig mobilityConfig{};
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});

    DDAssembler cached(mesh,
                       matdb,
                       doping,
                       constants::Vt_300,
                       mobilityConfig,
                       noRecombination);

    VectorXd psi(3);
    VectorXd n(3);
    VectorXd p(3);
    psi << 0.020, -0.010, 0.030;
    n << 1.0e16, 2.0e16, 4.0e16;
    p << 3.0e15, 1.5e15, 2.5e15;

    cached.assemblePoissonWithCarriers(n, p, psi);
    const AssemblySystem referencePoisson = assembleReferencePoissonWithFreshGeometry(
        mesh, matdb, doping, constants::Vt_300, n, p, psi);
    requireSystemsMatch(cached.matrix(), cached.rhs(), referencePoisson.A, referencePoisson.b);

    cached.assembleElectronContinuity(psi, n, p);
    const AssemblySystem referenceElectrons = assembleReferenceContinuityWithFreshGeometry(
        mesh, matdb, doping, constants::Vt_300, mobilityConfig, noRecombination,
        CarrierType::Electron, psi, n, p);
    requireSystemsMatch(cached.matrix(), cached.rhs(), referenceElectrons.A, referenceElectrons.b);

    cached.assembleHoleContinuity(psi, n, p);
    const AssemblySystem referenceHoles = assembleReferenceContinuityWithFreshGeometry(
        mesh, matdb, doping, constants::Vt_300, mobilityConfig, noRecombination,
        CarrierType::Hole, psi, n, p);
    requireSystemsMatch(cached.matrix(), cached.rhs(), referenceHoles.A, referenceHoles.b);
}
