#include "vela/equation/DDAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/physics/SRHRecombination.h"
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace vela {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

DDAssembler::DDAssembler(const DeviceMesh&       mesh,
                         const MaterialDatabase& matdb,
                         const DopingModel&      doping,
                         double                  Vt,
                         double                  taun,
                         double                  taup)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , Vt_(Vt)
    , taun_(taun)
    , taup_(taup)
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(VectorXd::Zero(static_cast<int>(mesh.numNodes())))
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            "DDAssembler: doping model size does not match mesh node count.");
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

std::vector<Real> DDAssembler::computeNodeVolumes() const
{
    const Index N = mesh_.numNodes();
    std::vector<Real> vol(N, 0.0);
    for (Index c = 0; c < mesh_.numCells(); ++c) {
        const auto& cell = mesh_.getCell(c);
        if (cell.node_ids.size() < 3) continue;
        const Node& n0 = mesh_.getNode(cell.node_ids[0]);
        const Node& n1 = mesh_.getNode(cell.node_ids[1]);
        const Node& n2 = mesh_.getNode(cell.node_ids[2]);
        const Real area = std::abs(
            (n1.x - n0.x) * (n2.y - n0.y) -
            (n2.x - n0.x) * (n1.y - n0.y)) * 0.5;
        vol[cell.node_ids[0]] += area / 3.0;
        vol[cell.node_ids[1]] += area / 3.0;
        vol[cell.node_ids[2]] += area / 3.0;
    }
    return vol;
}

std::vector<Real> DDAssembler::computeEdgeCouplings() const
{
    const Index E = mesh_.numEdges();
    std::vector<Real> couple(E);
    for (Index e = 0; e < E; ++e)
        couple[e] = mesh_.getEdge(e).length;
    return couple;
}

void DDAssembler::buildEdgeCellMap()
{
    edgeCells_.assign(mesh_.numEdges(), {});
    std::unordered_map<Index, Index> pairToEdge;
    const Index N = mesh_.numNodes();
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        pairToEdge[edge.n0 * N + edge.n1] = e;
    }
    for (Index c = 0; c < mesh_.numCells(); ++c) {
        const auto& cell = mesh_.getCell(c);
        if (cell.node_ids.size() < 3) continue;
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

// ---------------------------------------------------------------------------
// Material helpers (average over edge-adjacent cells)
// ---------------------------------------------------------------------------

static Real avgMaterialProp(const std::vector<Index>& cells,
                            const DeviceMesh&          mesh,
                            const MaterialDatabase&    matdb,
                            Real Material::*           prop,
                            Real                       fallback)
{
    if (cells.empty()) return fallback;
    Real sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        Real val = fallback;
        if (matdb.hasMaterial(region.material))
            val = matdb.getMaterial(region.material).*prop;
        sum += val;
    }
    return sum / static_cast<Real>(cells.size());
}

Real DDAssembler::edgeEpsilon(Index edgeId) const
{
    if (edgeCells_.empty())
        throw std::logic_error("DDAssembler: buildEdgeCellMap() not called.");
    return avgMaterialProp(edgeCells_[edgeId], mesh_, matdb_,
                           &Material::eps_r, 1.0) * constants::eps0;
}

Real DDAssembler::edgeMun(Index edgeId) const
{
    if (edgeCells_.empty())
        throw std::logic_error("DDAssembler: buildEdgeCellMap() not called.");
    return avgMaterialProp(edgeCells_[edgeId], mesh_, matdb_,
                           &Material::mun, 0.0);
}

Real DDAssembler::edgeMup(Index edgeId) const
{
    if (edgeCells_.empty())
        throw std::logic_error("DDAssembler: buildEdgeCellMap() not called.");
    return avgMaterialProp(edgeCells_[edgeId], mesh_, matdb_,
                           &Material::mup, 0.0);
}

Real DDAssembler::nodeNi(Index nodeId) const
{
    // Return ni from any cell that contains this node.
    // Since all cells here use the same material (Si), a simple lookup
    // over all cells is sufficient for the prototype.
    for (Index c = 0; c < mesh_.numCells(); ++c) {
        const auto& cell = mesh_.getCell(c);
        for (Index nid : cell.node_ids) {
            if (nid == nodeId) {
                const auto& region = mesh_.getRegion(cell.region_id);
                if (matdb_.hasMaterial(region.material))
                    return matdb_.getMaterial(region.material).ni;
                break;
            }
        }
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

void DDAssembler::applyDirichlet(const std::unordered_map<Index, Real>& bcs)
{
    A_.makeCompressed();

    // Step 1: propagate prescribed values into free-node RHS via column i
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        for (SparseMatrixd::InnerIterator it(A_, i); it; ++it) {
            const int k = static_cast<int>(it.row());
            if (k == i) continue;
            if (bcs.count(static_cast<Index>(k)) == 0)
                b_(k) -= it.value() * value;
        }
    }

    // Step 2 & 3: zero Dirichlet rows/cols, set diagonal and RHS
    std::unordered_set<int> dirichletCols;
    for (const auto& [nodeId, _] : bcs)
        dirichletCols.insert(static_cast<int>(nodeId));

    for (int col = 0; col < A_.outerSize(); ++col) {
        const bool colIsDirichlet = dirichletCols.count(col) > 0;
        for (SparseMatrixd::InnerIterator it(A_, col); it; ++it) {
            const int  row          = static_cast<int>(it.row());
            const bool rowIsDirichlet = dirichletCols.count(row) > 0;
            if (rowIsDirichlet || colIsDirichlet)
                it.valueRef() = (row == col) ? 1.0 : 0.0;
        }
    }
    A_.prune(0.0);

    for (const auto& [nodeId, value] : bcs)
        b_(static_cast<int>(nodeId)) = value;
}

// ---------------------------------------------------------------------------
// Poisson with carriers
// ---------------------------------------------------------------------------

void DDAssembler::assemblePoissonWithCarriers(const VectorXd& n,
                                              const VectorXd& p,
                                              const VectorXd& psi)
{
    const Index N = mesh_.numNodes();
    buildEdgeCellMap();

    const std::vector<Real> vol    = computeNodeVolumes();
    const std::vector<Real> couple = computeEdgeCouplings();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // Edge flux terms (same structure as PoissonAssembler)
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real eps = edgeEpsilon(e);
        const Real G   = eps * couple[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        triplets.emplace_back(i, i,  G);
        triplets.emplace_back(j, j,  G);
        triplets.emplace_back(i, j, -G);
        triplets.emplace_back(j, i, -G);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // Linearised carrier charge:
    //   q*(n+p)/Vt added to diagonal; adjusted RHS keeps the equation
    //   equivalent to the nonlinear form around the current ψ.
    for (Index i = 0; i < N; ++i) {
        const int   ii   = static_cast<int>(i);
        const Real  ni_v = n(ii);
        const Real  pi_v = p(ii);
        const Real  vol_i = vol[i];

        // Linearisation diagonal contribution
        const Real diagCarrier = constants::q * (ni_v + pi_v) / Vt_ * vol_i;
        A_.coeffRef(ii, ii) += diagCarrier;

        // RHS: q*(p-n+Nd-Na)*vol + linearisation shift term
        b_(ii) = constants::q *
                 (pi_v - ni_v + doping_.netDoping(i)) * vol_i
                 + diagCarrier * psi(ii);
    }
}

// ---------------------------------------------------------------------------
// Electron continuity
// ---------------------------------------------------------------------------

void DDAssembler::assembleElectronContinuity(const VectorXd& psi,
                                             const VectorXd& n_old,
                                             const VectorXd& p_old)
{
    const Index N = mesh_.numNodes();
    buildEdgeCellMap();

    const std::vector<Real> vol    = computeNodeVolumes();
    const std::vector<Real> couple = computeEdgeCouplings();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // SG matrix entries from all edges
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real mun    = edgeMun(e);
        if (mun <= 0.0) continue;

        const Real coef  = mun * Vt_ * couple[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        const Real dpsi = psi(j) - psi(i);
        const Real u    = dpsi / Vt_;
        const Real Bu   = bernoulli( u);
        const Real Bmu  = bernoulli(-u);

        // Diagonal: coef * B(-u) for n_i
        triplets.emplace_back(i, i,  coef * Bmu);
        triplets.emplace_back(j, j,  coef * Bu);
        // Off-diagonal: -coef * B(+u) for n_j
        triplets.emplace_back(i, j, -coef * Bu);
        triplets.emplace_back(j, i, -coef * Bmu);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // SRH source term (linearised w.r.t. n to maintain diagonal dominance):
    //   R = n*p/D - ni²/D   (D = τp*(n+ni) + τn*(p+ni))
    // Move n-term to LHS diagonal; ni² term goes to RHS.
    for (Index i = 0; i < N; ++i) {
        const int   ii    = static_cast<int>(i);
        const Real  ni_v  = nodeNi(i);
        const Real  n_v   = n_old(ii);
        const Real  p_v   = p_old(ii);
        const Real  vol_i = vol[i];

        const Real D = taup_ * (n_v + ni_v) + taun_ * (p_v + ni_v);
        if (D < 1.0e-100) continue;

        // p/D * vol → positive diagonal addition
        A_.coeffRef(ii, ii) += (p_v / D) * vol_i;

        // ni²/D * vol → non-negative RHS
        b_(ii) += (ni_v * ni_v / D) * vol_i;
    }
}

// ---------------------------------------------------------------------------
// Hole continuity
// ---------------------------------------------------------------------------

void DDAssembler::assembleHoleContinuity(const VectorXd& psi,
                                         const VectorXd& n_old,
                                         const VectorXd& p_old)
{
    const Index N = mesh_.numNodes();
    buildEdgeCellMap();

    const std::vector<Real> vol    = computeNodeVolumes();
    const std::vector<Real> couple = computeEdgeCouplings();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // SG matrix entries for holes
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real mup   = edgeMup(e);
        if (mup <= 0.0) continue;

        const Real coef = mup * Vt_ * couple[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        const Real dpsi = psi(j) - psi(i);
        const Real u    = dpsi / Vt_;
        const Real Bu   = bernoulli( u);
        const Real Bmu  = bernoulli(-u);

        // Hole SG flux from i to j: Jp = μp*Vt/h*[B(-u)*p_j - B(+u)*p_i]
        // Diagonal: coef * B(+u) for p_i
        triplets.emplace_back(i, i,  coef * Bu);
        triplets.emplace_back(j, j,  coef * Bmu);
        // Off-diagonal: -coef * B(-u) for p_j
        triplets.emplace_back(i, j, -coef * Bmu);
        triplets.emplace_back(j, i, -coef * Bu);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // SRH source term (linearised w.r.t. p):
    //   R = n*p/D - ni²/D
    // Move p-term to LHS; ni² term goes to RHS.
    for (Index i = 0; i < N; ++i) {
        const int   ii    = static_cast<int>(i);
        const Real  ni_v  = nodeNi(i);
        const Real  n_v   = n_old(ii);
        const Real  p_v   = p_old(ii);
        const Real  vol_i = vol[i];

        const Real D = taup_ * (n_v + ni_v) + taun_ * (p_v + ni_v);
        if (D < 1.0e-100) continue;

        // n/D * vol → positive diagonal addition
        A_.coeffRef(ii, ii) += (n_v / D) * vol_i;

        // ni²/D * vol → non-negative RHS
        b_(ii) += (ni_v * ni_v / D) * vol_i;
    }
}

} // namespace vela
