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

CoupledDDAssembler::CoupledDDAssembler(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       double Vt,
                                       double taun,
                                       double taup)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , Vt_(Vt)
    , taun_(taun)
    , taup_(taup)
    , ni_(detail::buildNodeNi(mesh, matdb))
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

    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    VectorXd r = VectorXd::Zero(3 * N);
    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    const auto edgeCells = detail::buildEdgeCellMap(mesh_);
    const auto vol = detail::computeNodeVolumes(mesh_);
    const auto couple = detail::computeEdgeCouplings(mesh_);

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real dpsi = x(psiOffset() + j) - x(psiOffset() + i);
        const Real u = dpsi / Vt_;
        const Real Bu = bernoulli(u);
        const Real Bmu = bernoulli(-u);

        const Real eps = detail::edgeEpsilon(edgeCells, mesh_, matdb_, e);
        const Real G = eps * couple[e] / h;
        const Real psiFlux = G * (x(psiOffset() + i) - x(psiOffset() + j));
        r(psiOffset() + i) += psiFlux;
        r(psiOffset() + j) -= psiFlux;

        const Real mun = detail::edgeAvgMaterialProp(
            edgeCells[e], mesh_, matdb_, &Material::mun, 0.0);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mun * Vt_ * couple[e] / h;
            r(phinOffset() + i) += coef * (Bmu * n(i) - Bu * n(j));
            r(phinOffset() + j) += coef * (Bu * n(j) - Bmu * n(i));
        }

        const Real mup = detail::edgeAvgMaterialProp(
            edgeCells[e], mesh_, matdb_, &Material::mup, 0.0);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mup * Vt_ * couple[e] / h;
            r(phipOffset() + i) += coef * (Bu * p(i) - Bmu * p(j));
            r(phipOffset() + j) += coef * (Bmu * p(j) - Bu * p(i));
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        r(psiOffset() + ii) -= constants::q *
            (p(ii) - n(ii) + doping_.netDoping(i)) * vol[i];

        const Real ni = ni_[i];
        const Real D = taup_ * (n(ii) + ni) + taun_ * (p(ii) + ni);
        if (D > 1.0e-100) {
            const Real R = (n(ii) * p(ii) - ni * ni) / D;
            r(phinOffset() + ii) += R * vol[i];
            r(phipOffset() + ii) += R * vol[i];
            hasElectronContribution[static_cast<std::size_t>(ii)] = true;
            hasHoleContribution[static_cast<std::size_t>(ii)] = true;
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

SparseMatrixd CoupledDDAssembler::finiteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    const int M = x.size();
    const VectorXd r0 = residual(x, bcs);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(M) * static_cast<std::size_t>(M));

    for (int col = 0; col < M; ++col) {
        const Real h = relativeStep * std::max(1.0, std::abs(x(col)));
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
