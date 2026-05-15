#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/RecombinationModel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

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

TEST_CASE("CoupledDDAssembler BGN residuals match DDAssembler with nonuniform ni", "[sg][dd][coupled][bgn]")
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
