#include "vela/equation/PoissonAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace vela {

PoissonAssembler::PoissonAssembler(
    const DeviceMesh&      mesh,
    const MaterialDatabase& matdb,
    const DopingModel&      doping,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    std::vector<PoissonNeumannBoundarySpec> neumannBoundaries)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , fixedCharges_(std::move(fixedCharges))
    , sheetCharges_(std::move(sheetCharges))
    , neumannBoundaries_(std::move(neumannBoundaries))
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

    // ---- Mobile-free source term: rhs_i = +q * netDoping_i * vol_i ----
    for (Index i = 0; i < N; ++i)
        b_(static_cast<int>(i)) = constants::q * doping_.netDoping(i) * vol[i];

    // ---- Region fixed charge: q * fixed_charge_m3 * cell_area / 3 ----
    detail::addFixedAndInterfaceChargeToRhs(
        mesh_, edgeCells, fixedCharges_, sheetCharges_, b_, "PoissonAssembler");

    // ---- Neumann boundary conditions ----
    // For each boundary segment defined by a polyline of node IDs, compute the
    // RHS contribution from the normal displacement D.n [C/m^2].
    // For each edge in the polyline: rhs += D_n * edge_length / 2 to each endpoint.
    for (const auto& neumannSpec : neumannBoundaries_) {
        if (neumannSpec.node_ids.size() < 2) continue;

        for (size_t i = 0; i + 1 < neumannSpec.node_ids.size(); ++i) {
            const Index n0 = neumannSpec.node_ids[i];
            const Index n1 = neumannSpec.node_ids[i + 1];

            if (n0 >= N || n1 >= N) {
                throw std::out_of_range(
                    "PoissonAssembler: Neumann boundary node ID out of range.");
            }

            const Node& node0 = mesh_.getNode(n0);
            const Node& node1 = mesh_.getNode(n1);
            const Real dx = node1.x - node0.x;
            const Real dy = node1.y - node0.y;
            const Real edgeLength = std::sqrt(dx * dx + dy * dy);

            if (edgeLength < 1e-30) continue; // Skip degenerate edges

            const Real endpointContribution = neumannSpec.normalDisplacement * edgeLength * 0.5;
            b_(static_cast<int>(n0)) += endpointContribution;
            b_(static_cast<int>(n1)) += endpointContribution;
        }
    }
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
