#include "vela/equation/CoupledDDAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/physics/CarrierStatistics.h"
#include <Eigen/Sparse>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace vela {
namespace {

Real bernoulliDerivative(Real x)
{
    const Real ax = std::abs(x);
    if (ax < 1.0e-8)
        return -0.5 + x / 6.0 - x * x * x / 180.0;
    if (x > 500.0)
        return (1.0 - x) * std::exp(-x);
    if (x < -500.0)
        return -1.0;

    const Real em1 = std::expm1(x);
    const Real ex = em1 + 1.0;
    return (em1 - x * ex) / (em1 * em1);
}

} // namespace

CoupledDDAssembler::CoupledDDAssembler(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       double Vt,
                                       double taun,
                                       double taup)
    : CoupledDDAssembler(mesh,
                         matdb,
                         doping,
                         Vt,
                         MobilityModelConfig{},
                         recombinationModelConfig({"srh"}, taun, taup))
{}

CoupledDDAssembler::CoupledDDAssembler(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    const MobilityModelConfig& mobilityConfig,
    const RecombinationModelConfig& recombinationConfig)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , Vt_(Vt)
    , mobility_(makeMobilityModel(mobilityConfig))
    , recombination_(recombinationConfig)
    , ni_(detail::buildNodeNi(mesh, matdb))
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , vol_(detail::computeNodeVolumes(mesh))
    , couple_(detail::computeEdgeCouplings(mesh))
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            "CoupledDDAssembler: doping model size does not match mesh node count.");
}

VectorXd CoupledDDAssembler::pack(const CoupledDDState& state) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N)
        throw std::invalid_argument("CoupledDDAssembler::pack: state vector size mismatch.");

    VectorXd x(3 * N);
    x.segment(psiOffset(), N) = state.psi;
    x.segment(phinOffset(), N) = state.phin;
    x.segment(phipOffset(), N) = state.phip;
    return x;
}

CoupledDDState CoupledDDAssembler::unpack(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::unpack: vector size mismatch.");

    CoupledDDState state;
    state.psi = x.segment(psiOffset(), N);
    state.phin = x.segment(phinOffset(), N);
    state.phip = x.segment(phipOffset(), N);
    return state;
}

VectorXd CoupledDDAssembler::electronDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd n(N);
    for (int i = 0; i < N; ++i)
        n(i) = vela::electronDensity(ni_[static_cast<Index>(i)],
                                     x(psiOffset() + i),
                                     x(phinOffset() + i),
                                     Vt_);
    return n;
}

VectorXd CoupledDDAssembler::holeDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd p(N);
    for (int i = 0; i < N; ++i)
        p(i) = vela::holeDensity(ni_[static_cast<Index>(i)],
                                 x(psiOffset() + i),
                                 x(phipOffset() + i),
                                 Vt_);
    return p;
}

bool CoupledDDAssembler::hasPositiveFiniteCarriers(const VectorXd& x) const
{
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    for (int i = 0; i < n.size(); ++i) {
        if (!std::isfinite(n(i)) || !std::isfinite(p(i)) || n(i) <= 0.0 || p(i) <= 0.0)
            return false;
    }
    return true;
}

VectorXd CoupledDDAssembler::residual(const VectorXd& x,
                                      const CoupledDDBoundaryConditions& bcs) const
{
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::residual: vector size mismatch.");

    // n and p are needed for Poisson source and configured recombination.
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    VectorXd r = VectorXd::Zero(3 * N);
    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    // Pre-compute per-node exponentials used in the balanced SG flux formulas.
    // Using the Bernoulli identity B(-u) = B(u)*exp(u), the standard SG flux
    //   B(-u)*n_i - B(u)*n_j   (electrons, ni=n0)
    //   B(u)*p_i  - B(-u)*p_j  (holes,     ni=n0)
    // can be rewritten without catastrophic cancellation at equilibrium as:
    //   B(u) * ni_i * exp(psi_j/Vt) * [exp(-phin_i/Vt) - exp(-phin_j/Vt)]
    //   B(u) * ni_i * exp(-psi_i/Vt) * [exp(phip_i/Vt) - exp(phip_j/Vt)]
    // At equilibrium (phin=phip=0), the bracketed differences are exactly 0.
    std::vector<Real> expNegPhin(static_cast<std::size_t>(N));
    std::vector<Real> expPhip(static_cast<std::size_t>(N));
    std::vector<Real> expPsi(static_cast<std::size_t>(N));
    std::vector<Real> expNegPsi(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) {
        expNegPhin[static_cast<std::size_t>(k)] = std::exp(-x(phinOffset() + k) / Vt_);
        expPhip[static_cast<std::size_t>(k)]    = std::exp( x(phipOffset()  + k) / Vt_);
        expPsi[static_cast<std::size_t>(k)]     = std::exp( x(psiOffset()   + k) / Vt_);
        expNegPsi[static_cast<std::size_t>(k)]  = std::exp(-x(psiOffset()   + k) / Vt_);
    }

    const auto& edgeCells = edgeCells_;
    const auto& vol = vol_;
    const auto& couple = couple_;

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real dpsi = x(psiOffset() + j) - x(psiOffset() + i);
        const Real u = dpsi / Vt_;
        const Real Bu = bernoulli(u);

        const Real eps = detail::edgeEpsilon(edgeCells, mesh_, matdb_, e);
        const Real G = eps * couple[e] / h;
        const Real psiFlux = G * (x(psiOffset() + i) - x(psiOffset() + j));
        r(psiOffset() + i) += psiFlux;
        r(psiOffset() + j) -= psiFlux;

        const Real mun = detail::edgeMobility(
            edgeCells, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Electron);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            // Balanced electron SG flux (avoids catastrophic cancellation):
            //   B(-u)*n_i - B(u)*n_j = B(u)*ni_i*exp(psi_j/Vt)
            //                           * [exp(-phin_i/Vt) - exp(-phin_j/Vt)]
            const Real coef = mun * Vt_ * couple[e] / h;
            const Real ni_i = ni_[static_cast<Index>(i)];
            const Real nFlux = coef * Bu * ni_i
                               * expPsi[static_cast<std::size_t>(j)]
                               * (expNegPhin[static_cast<std::size_t>(i)]
                                  - expNegPhin[static_cast<std::size_t>(j)]);
            r(phinOffset() + i) += nFlux;
            r(phinOffset() + j) -= nFlux;
        }

        const Real mup = detail::edgeMobility(
            edgeCells, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Hole);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            // Balanced hole SG flux (avoids catastrophic cancellation):
            //   B(u)*p_i - B(-u)*p_j = B(u)*ni_i*exp(-psi_i/Vt)
            //                           * [exp(phip_i/Vt) - exp(phip_j/Vt)]
            const Real coef = mup * Vt_ * couple[e] / h;
            const Real ni_i = ni_[static_cast<Index>(i)];
            const Real pFlux = coef * Bu * ni_i
                               * expNegPsi[static_cast<std::size_t>(i)]
                               * (expPhip[static_cast<std::size_t>(i)]
                                  - expPhip[static_cast<std::size_t>(j)]);
            r(phipOffset() + i) += pFlux;
            r(phipOffset() + j) -= pFlux;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        r(psiOffset() + ii) -= constants::q *
            (p(ii) - n(ii) + doping_.netDoping(i)) * vol[i];

        const Real ni = ni_[i];
        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            // Compute n*p - ni^2 via the identity n*p = ni^2 * exp((phip-phin)/Vt),
            // i.e. n*p - ni^2 = ni^2 * expm1((phip-phin)/Vt). This avoids
            // catastrophic cancellation when phip ~= phin (near equilibrium), where
            // the naive form n*p - ni^2 is dominated by floating-point rounding
            // in exp(+u) * exp(-u) != 1.
            const Real dPhi = x(phipOffset() + ii) - x(phinOffset() + ii);
            const Real excessProduct = (ni > 0.0)
                ? ni * ni * std::expm1(dPhi / Vt_)
                : n(ii) * p(ii);
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);
            if (R != 0.0) {
                r(phinOffset() + ii) += R * vol[i];
                r(phipOffset() + ii) += R * vol[i];
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    // Insulating nodes such as SiO2 (mun = mup = ni = 0) can have no
    // transport or recombination contribution in the continuity equations.
    // Pin those otherwise unconstrained quasi-Fermi unknowns to avoid zero
    // residual/Jacobian rows. Explicit boundary conditions below take
    // precedence over this internal gauge constraint.
    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            r(phinOffset() + ii) = x(phinOffset() + ii);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            r(phipOffset() + ii) = x(phipOffset() + ii);
    }

    for (const auto& [node, value] : bcs.psi)
        r(psiOffset() + static_cast<int>(node)) = x(psiOffset() + static_cast<int>(node)) - value;
    for (const auto& [node, value] : bcs.phin)
        r(phinOffset() + static_cast<int>(node)) = x(phinOffset() + static_cast<int>(node)) - value;
    for (const auto& [node, value] : bcs.phip)
        r(phipOffset() + static_cast<int>(node)) = x(phipOffset() + static_cast<int>(node)) - value;

    return r;
}

SparseMatrixd CoupledDDAssembler::assembleJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::assembleJacobian: vector size mismatch.");

    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);

    std::vector<bool> constrainedRows(static_cast<std::size_t>(3 * N), false);
    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(psiOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phinOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phipOffset() + static_cast<int>(node))] = true;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(N) * 27);
    auto add = [&](int row, int col, Real value) {
        if (value != 0.0 && !constrainedRows[static_cast<std::size_t>(row)])
            triplets.emplace_back(row, col, value);
    };

    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);
    std::vector<Real> expNegPhin(static_cast<std::size_t>(N));
    std::vector<Real> expPhip(static_cast<std::size_t>(N));
    std::vector<Real> expPsi(static_cast<std::size_t>(N));
    std::vector<Real> expNegPsi(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) {
        expNegPhin[static_cast<std::size_t>(k)] = std::exp(-x(phinOffset() + k) / Vt_);
        expPhip[static_cast<std::size_t>(k)]    = std::exp( x(phipOffset()  + k) / Vt_);
        expPsi[static_cast<std::size_t>(k)]     = std::exp( x(psiOffset()   + k) / Vt_);
        expNegPsi[static_cast<std::size_t>(k)]  = std::exp(-x(psiOffset()   + k) / Vt_);
    }

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real dpsi = x(psiOffset() + j) - x(psiOffset() + i);
        const Real u = dpsi / Vt_;
        const Real Bu = bernoulli(u);
        const Real dBu = bernoulliDerivative(u);

        const Real eps = detail::edgeEpsilon(edgeCells_, mesh_, matdb_, e);
        const Real G = eps * couple_[e] / h;
        add(psiOffset() + i, psiOffset() + i,  G);
        add(psiOffset() + i, psiOffset() + j, -G);
        add(psiOffset() + j, psiOffset() + i, -G);
        add(psiOffset() + j, psiOffset() + j,  G);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Electron);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mun * Vt_ * couple_[e] / h;
            const Real factor = coef * ni_[static_cast<Index>(i)]
                              * expPsi[static_cast<std::size_t>(j)];
            const Real diff = expNegPhin[static_cast<std::size_t>(i)]
                            - expNegPhin[static_cast<std::size_t>(j)];
            const Real dF_dpsi_i = factor * (-dBu / Vt_) * diff;
            const Real dF_dpsi_j = factor * ((dBu + Bu) / Vt_) * diff;
            const Real dF_dphin_i = factor * Bu
                                  * (-expNegPhin[static_cast<std::size_t>(i)] / Vt_);
            const Real dF_dphin_j = factor * Bu
                                  * ( expNegPhin[static_cast<std::size_t>(j)] / Vt_);

            add(phinOffset() + i, psiOffset() + i, dF_dpsi_i);
            add(phinOffset() + i, psiOffset() + j, dF_dpsi_j);
            add(phinOffset() + i, phinOffset() + i, dF_dphin_i);
            add(phinOffset() + i, phinOffset() + j, dF_dphin_j);
            add(phinOffset() + j, psiOffset() + i, -dF_dpsi_i);
            add(phinOffset() + j, psiOffset() + j, -dF_dpsi_j);
            add(phinOffset() + j, phinOffset() + i, -dF_dphin_i);
            add(phinOffset() + j, phinOffset() + j, -dF_dphin_j);
        }

        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Hole);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mup * Vt_ * couple_[e] / h;
            const Real factor = coef * ni_[static_cast<Index>(i)]
                              * expNegPsi[static_cast<std::size_t>(i)];
            const Real diff = expPhip[static_cast<std::size_t>(i)]
                            - expPhip[static_cast<std::size_t>(j)];
            const Real dF_dpsi_i = factor * (-(dBu + Bu) / Vt_) * diff;
            const Real dF_dpsi_j = factor * (dBu / Vt_) * diff;
            const Real dF_dphip_i = factor * Bu
                                  * ( expPhip[static_cast<std::size_t>(i)] / Vt_);
            const Real dF_dphip_j = factor * Bu
                                  * (-expPhip[static_cast<std::size_t>(j)] / Vt_);

            add(phipOffset() + i, psiOffset() + i, dF_dpsi_i);
            add(phipOffset() + i, psiOffset() + j, dF_dpsi_j);
            add(phipOffset() + i, phipOffset() + i, dF_dphip_i);
            add(phipOffset() + i, phipOffset() + j, dF_dphip_j);
            add(phipOffset() + j, psiOffset() + i, -dF_dpsi_i);
            add(phipOffset() + j, psiOffset() + j, -dF_dpsi_j);
            add(phipOffset() + j, phipOffset() + i, -dF_dphip_i);
            add(phipOffset() + j, phipOffset() + j, -dF_dphip_j);
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        const Real dni_dpsi = n(ii) / Vt_;
        const Real dni_dphin = -n(ii) / Vt_;
        const Real dpi_dpsi = -p(ii) / Vt_;
        const Real dpi_dphip = p(ii) / Vt_;

        add(psiOffset() + ii, psiOffset() + ii,
            -constants::q * (dpi_dpsi - dni_dpsi) * vol_[i]);
        add(psiOffset() + ii, phinOffset() + ii,
            -constants::q * (-dni_dphin) * vol_[i]);
        add(psiOffset() + ii, phipOffset() + ii,
            -constants::q * dpi_dphip * vol_[i]);

        const Real ni = ni_[i];
        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            const Real dPhi = x(phipOffset() + ii) - x(phinOffset() + ii);
            const Real np = (ni > 0.0) ? ni * ni * std::exp(dPhi / Vt_) : n(ii) * p(ii);
            const Real excessProduct = (ni > 0.0)
                ? ni * ni * std::expm1(dPhi / Vt_)
                : n(ii) * p(ii);
            const auto deriv = recombination_.totalRateDerivativesFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);

            const Real dExcess_dphin = -np / Vt_;
            const Real dExcess_dphip =  np / Vt_;
            const Real dR_dpsi = deriv.dRateDn * dni_dpsi + deriv.dRateDp * dpi_dpsi;
            const Real dR_dphin = deriv.dRateDn * dni_dphin
                                 + deriv.dRateDExcess * dExcess_dphin;
            const Real dR_dphip = deriv.dRateDp * dpi_dphip
                                 + deriv.dRateDExcess * dExcess_dphip;

            add(phinOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i]);
            add(phinOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i]);
            add(phinOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i]);
            add(phipOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i]);
            add(phipOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i]);
            add(phipOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i]);

            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);
            if (R != 0.0) {
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            add(phinOffset() + ii, phinOffset() + ii, 1.0);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            add(phipOffset() + ii, phipOffset() + ii, 1.0);
    }

    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        triplets.emplace_back(psiOffset() + static_cast<int>(node),
                              psiOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        triplets.emplace_back(phinOffset() + static_cast<int>(node),
                              phinOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        triplets.emplace_back(phipOffset() + static_cast<int>(node),
                              phipOffset() + static_cast<int>(node), 1.0);
    }

    SparseMatrixd J(3 * N, 3 * N);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

SparseMatrixd CoupledDDAssembler::finiteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    if (relativeStep <= 0.0)
        throw std::invalid_argument(
            "CoupledDDAssembler::finiteDifferenceJacobian: relativeStep must be positive.");

    const int M = x.size();
    const VectorXd r0 = residual(x, bcs);
    std::vector<Eigen::Triplet<double>> triplets;
    // Heuristic: for a 2-D mesh with M = 3*N unknowns and sparse connectivity
    // (average ~6-7 neighbours per node), the Jacobian has O(N * 7 * 9) ~= 63*N
    // non-zeros.  Reserving M * 7 avoids the M^2 allocation that would OOM for
    // any realistically-sized mesh while still keeping reallocations rare.
    triplets.reserve(static_cast<std::size_t>(M) * 7);

    // Minimum absolute perturbation to prevent h == 0 when x(col) == 0.
    constexpr Real minAbsStep = 1.0e-15;

    for (int col = 0; col < M; ++col) {
        const Real h = std::max(relativeStep * std::max(1.0, std::abs(x(col))), minAbsStep);
        VectorXd xp = x;
        xp(col) += h;
        const VectorXd rp = residual(xp, bcs);
        const VectorXd dr = (rp - r0) / h;
        for (int row = 0; row < M; ++row) {
            if (dr(row) != 0.0)
                triplets.emplace_back(row, col, dr(row));
        }
    }

    SparseMatrixd J(M, M);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

} // namespace vela
