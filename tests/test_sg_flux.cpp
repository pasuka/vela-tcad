#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"
#include <cmath>

using namespace vela;
using Catch::Approx;

// ---------------------------------------------------------------------------
// Pure diffusion (ψ_j = ψ_i → dpsi = 0, u = 0, B(0) = 1)
// ---------------------------------------------------------------------------

TEST_CASE("SG electron flux: pure diffusion (dpsi=0)", "[sg]")
{
    const double Vt = 0.02585; // ~300 K thermal voltage
    const double mu = 0.135;   // Si electron mobility
    const double h  = 1.0e-6;  // 1 µm edge

    const double n0 = 1.0e23;  // high concentration (n-side)
    const double n1 = 1.0e16;  // low concentration (p-side)

    // With dpsi = 0: Jn = μ*Vt/h * (n1 - n0)
    const double Jn = sgElectronFlux(n0, n1, 0.0, Vt, mu, h);
    const double expected = mu * Vt / h * (n1 - n0);
    REQUIRE(Jn == Approx(expected).epsilon(1.0e-10));

    // Direction: electrons diffuse from high n0 to low n1 (i→j direction
    // has conventional current flowing j→i) → Jn < 0
    REQUIRE(Jn < 0.0);
}

TEST_CASE("SG hole flux: pure diffusion (dpsi=0)", "[sg]")
{
    const double Vt = 0.02585;
    const double mu = 0.048;
    const double h  = 1.0e-6;

    const double p0 = 1.0e16;  // low (n-side)
    const double p1 = 1.0e23;  // high (p-side)

    // With dpsi = 0: Jp = μp*Vt/h * (p1 - p0)  (same as diffusion)
    const double Jp = sgHoleFlux(p0, p1, 0.0, Vt, mu, h);
    const double expected = mu * Vt / h * (p1 - p0);
    REQUIRE(Jp == Approx(expected).epsilon(1.0e-10));
    REQUIRE(Jp > 0.0); // holes diffuse from p1→p0 → conventional current i→j positive
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

    // Jn = μ*Vt/h * (B(u)*n1 - B(-u)*n0)
    //    = μ*Vt/h * (B(u)*ni*exp(psi1/Vt) - B(-u)*ni*exp(psi0/Vt))
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
    const double dpsi = 0.3;   // ψ_j - ψ_i

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
