#include "vela/equation/PoissonAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include <Eigen/Sparse>
#include <stdexcept>

namespace vela {

PoissonAssembler::PoissonAssembler(const DeviceMesh&      mesh,
                                   const MaterialDatabase& matdb,
                                   const DopingModel&      doping)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(VectorXd::Zero(static_cast<int>(mesh.numNodes())))
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            "PoissonAssembler: doping model size does not match mesh node count.");
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

void PoissonAssembler::assemble()
{
    const Index N = mesh_.numNodes();

    const auto edgeCells = detail::buildEdgeCellMap(mesh_);
    const auto vol       = detail::computeNodeVolumes(mesh_);
    const auto couple    = detail::computeEdgeCouplings(mesh_);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // ---- Off-diagonal terms from edge fluxes ----
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1e-30) continue; // degenerate edge guard

        const Real eps = detail::edgeEpsilon(edgeCells, mesh_, matdb_, e);
        const Real G   = eps * couple[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        triplets.emplace_back(i, i,  G);
        triplets.emplace_back(j, j,  G);
        triplets.emplace_back(i, j, -G);
        triplets.emplace_back(j, i, -G);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // ---- Source term: rhs_i = +q * netDoping_i * vol_i ----
    for (Index i = 0; i < N; ++i)
        b_(static_cast<int>(i)) = constants::q * doping_.netDoping(i) * vol[i];
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

void PoissonAssembler::applyDirichlet(
    const std::unordered_map<Index, Real>& bcs)
{
    detail::applyDirichletBC(A_, b_, bcs);
}

} // namespace vela
