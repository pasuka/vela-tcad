#include "vela/equation/DDAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include <Eigen/Sparse>
#include <stdexcept>

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
    , ni_(detail::buildNodeNi(mesh, matdb))
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(VectorXd::Zero(static_cast<int>(mesh.numNodes())))
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            "DDAssembler: doping model size does not match mesh node count.");
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

void DDAssembler::applyDirichlet(const std::unordered_map<Index, Real>& bcs)
{
    detail::applyDirichletBC(A_, b_, bcs);
}

// ---------------------------------------------------------------------------
// Poisson with carriers
// ---------------------------------------------------------------------------

void DDAssembler::assemblePoissonWithCarriers(const VectorXd& n,
                                              const VectorXd& p,
                                              const VectorXd& psi)
{
    const Index N = mesh_.numNodes();

    const auto edgeCells = detail::buildEdgeCellMap(mesh_);
    const auto vol       = detail::computeNodeVolumes(mesh_);
    const auto couple    = detail::computeEdgeCouplings(mesh_);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // Edge flux terms
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

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

    // Linearised carrier charge:
    //   q*(n+p)/Vt added to diagonal; adjusted RHS keeps the equation
    //   equivalent to the nonlinear form around the current ψ.
    for (Index i = 0; i < N; ++i) {
        const int  ii     = static_cast<int>(i);
        const Real ni_v   = n(ii);
        const Real pi_v   = p(ii);
        const Real vol_i  = vol[i];

        const Real diagCarrier = constants::q * (ni_v + pi_v) / Vt_ * vol_i;
        A_.coeffRef(ii, ii) += diagCarrier;

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

    const auto edgeCells = detail::buildEdgeCellMap(mesh_);
    const auto vol       = detail::computeNodeVolumes(mesh_);
    const auto couple    = detail::computeEdgeCouplings(mesh_);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // SG matrix entries from all edges
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real mun = detail::edgeAvgMaterialProp(
            edgeCells[e], mesh_, matdb_, &Material::mun, 0.0);
        if (mun <= 0.0) continue; // skip insulator edges

        const Real coef = mun * Vt_ * couple[e] / h;

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
        const int  ii    = static_cast<int>(i);
        const Real ni_i  = ni_[i];
        const Real n_v   = n_old(ii);
        const Real p_v   = p_old(ii);
        const Real vol_i = vol[i];

        const Real D = taup_ * (n_v + ni_i) + taun_ * (p_v + ni_i);
        if (D < 1.0e-100) continue;

        // p/D * vol → positive diagonal addition
        A_.coeffRef(ii, ii) += (p_v / D) * vol_i;

        // ni²/D * vol → non-negative RHS
        b_(ii) += (ni_i * ni_i / D) * vol_i;
    }

    // Guard: if any diagonal is still zero (insulator node with all edges
    // skipped due to zero mobility), pin that node to n = 0 to avoid
    // a singular system.
    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        if (A_.coeff(ii, ii) == 0.0) {
            A_.coeffRef(ii, ii) = 1.0;
            b_(ii) = 0.0;
        }
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

    const auto edgeCells = detail::buildEdgeCellMap(mesh_);
    const auto vol       = detail::computeNodeVolumes(mesh_);
    const auto couple    = detail::computeEdgeCouplings(mesh_);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // SG matrix entries for holes
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real mup = detail::edgeAvgMaterialProp(
            edgeCells[e], mesh_, matdb_, &Material::mup, 0.0);
        if (mup <= 0.0) continue; // skip insulator edges

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
        const int  ii    = static_cast<int>(i);
        const Real ni_i  = ni_[i];
        const Real n_v   = n_old(ii);
        const Real p_v   = p_old(ii);
        const Real vol_i = vol[i];

        const Real D = taup_ * (n_v + ni_i) + taun_ * (p_v + ni_i);
        if (D < 1.0e-100) continue;

        // n/D * vol → positive diagonal addition
        A_.coeffRef(ii, ii) += (n_v / D) * vol_i;

        // ni²/D * vol → non-negative RHS
        b_(ii) += (ni_i * ni_i / D) * vol_i;
    }

    // Guard: pin insulator nodes (zero-diagonal) to p = 0
    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        if (A_.coeff(ii, ii) == 0.0) {
            A_.coeffRef(ii, ii) = 1.0;
            b_(ii) = 0.0;
        }
    }
}

} // namespace vela
