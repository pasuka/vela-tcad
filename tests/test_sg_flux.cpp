#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/RecombinationModel.h"
#include <algorithm>
#include <cmath>

using namespace vela;
using Catch::Approx;

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
