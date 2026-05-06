#include "vela/equation/PoissonAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

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
// Geometry helpers
// ---------------------------------------------------------------------------

std::vector<Real> PoissonAssembler::computeNodeVolumes() const
{
    const Index N = mesh_.numNodes();
    std::vector<Real> vol(N, 0.0);

    for (Index c = 0; c < mesh_.numCells(); ++c) {
        const auto& cell = mesh_.getCell(c);
        if (cell.node_ids.size() < 3) continue;

        const Node& n0 = mesh_.getNode(cell.node_ids[0]);
        const Node& n1 = mesh_.getNode(cell.node_ids[1]);
        const Node& n2 = mesh_.getNode(cell.node_ids[2]);

        // Signed area via cross product; take absolute value
        Real area = std::abs(
            (n1.x - n0.x) * (n2.y - n0.y) -
            (n2.x - n0.x) * (n1.y - n0.y)) * 0.5;

        // Distribute one-third of triangle area to each corner node
        vol[cell.node_ids[0]] += area / 3.0;
        vol[cell.node_ids[1]] += area / 3.0;
        vol[cell.node_ids[2]] += area / 3.0;
    }
    return vol;
}

std::vector<Real> PoissonAssembler::computeEdgeCouplings() const
{
    // Simplified approximation: couple = edge_length.
    // Replace with circumcenter-distance formula when Voronoi is computed.
    const Index E = mesh_.numEdges();
    std::vector<Real> couple(E);
    for (Index e = 0; e < E; ++e)
        couple[e] = mesh_.getEdge(e).length;
    return couple;
}

// ---------------------------------------------------------------------------
// Mapping helpers
// ---------------------------------------------------------------------------

void PoissonAssembler::buildEdgeCellMap()
{
    edgeCells_.assign(mesh_.numEdges(), {});

    // Build a sorted-pair → edge-id map for O(1) lookup
    std::unordered_map<Index, Index> pairToEdge; // key = n0*N + n1 (n0 < n1)
    const Index N = mesh_.numNodes();
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        pairToEdge[edge.n0 * N + edge.n1] = e;
    }

    for (Index c = 0; c < mesh_.numCells(); ++c) {
        const auto& cell = mesh_.getCell(c);
        if (cell.node_ids.size() < 3) continue;

        // Three edges per triangle
        const Index nids[3] = {
            cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};

        for (int k = 0; k < 3; ++k) {
            Index a = nids[k];
            Index b = nids[(k + 1) % 3];
            if (a > b) std::swap(a, b);

            auto it = pairToEdge.find(a * N + b);
            if (it != pairToEdge.end())
                edgeCells_[it->second].push_back(c);
        }
    }
}

Real PoissonAssembler::edgeEpsilon(Index edgeId) const
{
    // Average epsilon [F/m] over cells adjacent to this edge
    if (edgeCells_.empty())
        throw std::logic_error("PoissonAssembler: buildEdgeCellMap() not called.");

    const auto& cells = edgeCells_[edgeId];
    if (cells.empty())
        return constants::eps0; // fallback

    Real epsR_sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh_.getRegion(mesh_.getCell(c).region_id);
        Real epsR = 1.0;
        if (matdb_.hasMaterial(region.material))
            epsR = matdb_.getMaterial(region.material).eps_r;
        epsR_sum += epsR;
    }
    return (epsR_sum / cells.size()) * constants::eps0;
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

void PoissonAssembler::assemble()
{
    const Index N = mesh_.numNodes();
    buildEdgeCellMap();

    const std::vector<Real> vol    = computeNodeVolumes();
    const std::vector<Real> couple = computeEdgeCouplings();

    // Collect triplets for the sparse matrix
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // ---- Off-diagonal terms from edge fluxes ----
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1e-30) continue; // degenerate edge guard

        const Real eps = edgeEpsilon(e);
        const Real G   = eps * couple[e] / h; // conductance [F/m / m * m = F/m]

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        triplets.emplace_back(i, i,  G);
        triplets.emplace_back(j, j,  G);
        triplets.emplace_back(i, j, -G);
        triplets.emplace_back(j, i, -G);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // ---- Source term: rhs_i = +q * netDoping_i * vol_i ----
    for (Index i = 0; i < N; ++i) {
        b_(static_cast<int>(i)) =
            constants::q * doping_.netDoping(i) * vol[i];
    }
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

void PoissonAssembler::applyDirichlet(
    const std::unordered_map<Index, Real>& bcs)
{
    // For each constrained node i with prescribed value v:
    //   1. Subtract A(k,i)*v from rhs(k) for all free rows k
    //   2. Zero column i and row i
    //   3. Set A(i,i) = 1, rhs(i) = v
    //
    // A_ is stored in CSC format (Eigen default): InnerIterator(A_, j)
    // iterates over non-zeros in column j.

    A_.makeCompressed();

    // Step 1: propagate Dirichlet values into free-node RHS via column i
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        for (SparseMatrixd::InnerIterator it(A_, i); it; ++it) {
            const int k = static_cast<int>(it.row());
            if (k == i) continue;
            if (bcs.count(static_cast<Index>(k)) == 0)
                b_(k) -= it.value() * value;
        }
    }

    // Step 2 & 3: zero the Dirichlet rows and columns, fix diagonal and RHS.
    // Build a set of Dirichlet column indices for fast lookup.
    std::unordered_set<int> dirichletCols;
    for (const auto& [nodeId, _] : bcs)
        dirichletCols.insert(static_cast<int>(nodeId));

    // Iterate over every column; zero entries in Dirichlet rows, and zero
    // Dirichlet columns entirely (except the diagonal of Dirichlet nodes).
    for (int col = 0; col < A_.outerSize(); ++col) {
        const bool colIsDirichlet = dirichletCols.count(col) > 0;
        for (SparseMatrixd::InnerIterator it(A_, col); it; ++it) {
            const int row = static_cast<int>(it.row());
            const bool rowIsDirichlet = dirichletCols.count(row) > 0;

            if (rowIsDirichlet || colIsDirichlet) {
                // Keep the diagonal of Dirichlet nodes as 1
                it.valueRef() = (row == col) ? 1.0 : 0.0;
            }
        }
    }

    // Remove the explicit zeros from the sparse structure
    A_.prune(0.0);

    // Set prescribed RHS values
    for (const auto& [nodeId, value] : bcs)
        b_(static_cast<int>(nodeId)) = value;
}

} // namespace vela
