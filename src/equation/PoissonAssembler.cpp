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

namespace {

Real triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    if (cell.node_ids.size() < 3) return 0.0;

    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);

    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
}

struct RegionPairKey {
    std::string first;
    std::string second;

    bool operator==(const RegionPairKey& other) const
    {
        return first == other.first && second == other.second;
    }
};

struct RegionPairKeyHash {
    std::size_t operator()(const RegionPairKey& key) const
    {
        const std::hash<std::string> hash;
        return hash(key.first) ^ (hash(key.second) << 1U);
    }
};

RegionPairKey makeRegionPairKey(std::string a, std::string b)
{
    if (b < a)
        std::swap(a, b);
    return RegionPairKey{std::move(a), std::move(b)};
}

std::unordered_map<std::string, Real> fixedChargeByRegion(
    const std::vector<RegionFixedChargeSpec>& fixedCharges)
{
    std::unordered_map<std::string, Real> fixedByRegion;
    for (const auto& spec : fixedCharges) {
        const auto [_, inserted] = fixedByRegion.emplace(spec.region, spec.fixedCharge);
        if (!inserted)
            throw std::invalid_argument(
                "PoissonAssembler: duplicate fixed_charge_m3 for region '" +
                spec.region + "'.");
    }
    return fixedByRegion;
}

std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetChargeByRegionPair(
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges)
{
    std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetByRegionPair;
    for (const auto& spec : sheetCharges)
        sheetByRegionPair[makeRegionPairKey(spec.region0, spec.region1)] += spec.totalSheetCharge();
    return sheetByRegionPair;
}

} // namespace

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
    const auto fixedByRegion = fixedChargeByRegion(fixedCharges_);

    if (!fixedByRegion.empty()) {
        for (Index c = 0; c < mesh_.numCells(); ++c) {
            const Cell& cell = mesh_.getCell(c);
            const Region& region = mesh_.getRegion(cell.region_id);
            auto it = fixedByRegion.find(region.name);
            if (it == fixedByRegion.end()) continue;

            const Real nodeCharge = constants::q * it->second * triangleArea(mesh_, cell) / 3.0;
            for (Index nid : cell.node_ids)
                b_(static_cast<int>(nid)) += nodeCharge;
        }
    }

    // ---- Interface sheet charge on shared-node region-pair edges ----
    // Allocation rule: q * sheet_charge_m2 * edge_length / 2 to each endpoint.
    // Pre-index configured region pairs so assembly scans edges once. Multiple
    // sheet specs for the same unordered region pair are intentionally summed.
    const auto sheetByRegionPair = sheetChargeByRegionPair(sheetCharges_);
    if (!sheetByRegionPair.empty()) {
        for (Index e = 0; e < mesh_.numEdges(); ++e) {
            const auto& cells = edgeCells[e];
            if (cells.size() != 2) continue;

            const Region& r0 = mesh_.getRegion(mesh_.getCell(cells[0]).region_id);
            const Region& r1 = mesh_.getRegion(mesh_.getCell(cells[1]).region_id);
            const auto it = sheetByRegionPair.find(makeRegionPairKey(r0.name, r1.name));
            if (it == sheetByRegionPair.end()) continue;

            const Edge& edge = mesh_.getEdge(e);
            const Real endpointCharge = constants::q * it->second * edge.length * 0.5;
            b_(static_cast<int>(edge.n0)) += endpointCharge;
            b_(static_cast<int>(edge.n1)) += endpointCharge;
        }
    }

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
