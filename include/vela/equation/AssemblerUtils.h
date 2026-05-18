#pragma once

/**
 * @file AssemblerUtils.h
 * @brief Shared FVM assembly helpers used by PoissonAssembler and DDAssembler.
 *
 * All functions are free (non-member) to avoid duplication between assemblers.
 * Include this header from assembler .cpp files; do not expose it as part of
 * the public library API.
 */

#include "vela/core/Types.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/BandgapNarrowing.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>
#include <utility>

namespace vela::detail {


// ---------------------------------------------------------------------------
// Fixed and interface charge helpers
// ---------------------------------------------------------------------------

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
        std::size_t seed = hash(key.first);
        seed ^= hash(key.second) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

inline RegionPairKey makeRegionPairKey(std::string a, std::string b)
{
    if (b < a)
        std::swap(a, b);
    return RegionPairKey{std::move(a), std::move(b)};
}

inline Real triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    if (cell.node_ids.size() < 3) return 0.0;

    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);

    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
}

inline std::unordered_map<std::string, Real> fixedChargeByRegion(
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::string& context)
{
    std::unordered_map<std::string, Real> fixedByRegion;
    for (const auto& spec : fixedCharges) {
        const auto [_, inserted] = fixedByRegion.emplace(spec.region, spec.fixedCharge);
        if (!inserted)
            throw std::invalid_argument(
                context + ": duplicate fixed_charge_m3 for region '" + spec.region + "'.");
    }
    return fixedByRegion;
}

inline std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetChargeByRegionPair(
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges)
{
    std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetByRegionPair;
    for (const auto& spec : sheetCharges)
        sheetByRegionPair[makeRegionPairKey(spec.region0, spec.region1)] += spec.totalSheetCharge();
    return sheetByRegionPair;
}

inline VectorXd computeFixedAndInterfaceChargeRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    const std::string& context)
{
    VectorXd contribution = VectorXd::Zero(static_cast<int>(mesh.numNodes()));

    const auto fixedByRegion = fixedChargeByRegion(fixedCharges, context);
    if (!fixedByRegion.empty()) {
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            const Region& region = mesh.getRegion(cell.region_id);
            auto it = fixedByRegion.find(region.name);
            if (it == fixedByRegion.end()) continue;

            const Real nodeCharge = constants::q * it->second * triangleArea(mesh, cell) / 3.0;
            for (Index nid : cell.node_ids)
                contribution(static_cast<int>(nid)) += nodeCharge;
        }
    }

    const auto sheetByRegionPair = sheetChargeByRegionPair(sheetCharges);
    if (!sheetByRegionPair.empty()) {
        for (Index e = 0; e < mesh.numEdges(); ++e) {
            const auto& cells = edgeCells[e];
            if (cells.size() != 2) continue;

            const Region& r0 = mesh.getRegion(mesh.getCell(cells[0]).region_id);
            const Region& r1 = mesh.getRegion(mesh.getCell(cells[1]).region_id);
            const auto it = sheetByRegionPair.find(makeRegionPairKey(r0.name, r1.name));
            if (it == sheetByRegionPair.end()) continue;

            const Edge& edge = mesh.getEdge(e);
            const Real endpointCharge = constants::q * it->second * edge.length * 0.5;
            contribution(static_cast<int>(edge.n0)) += endpointCharge;
            contribution(static_cast<int>(edge.n1)) += endpointCharge;
        }
    }

    return contribution;
}

inline void addFixedAndInterfaceChargeToRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    VectorXd& rhs,
    const std::string& context)
{
    rhs += computeFixedAndInterfaceChargeRhs(mesh, edgeCells, fixedCharges, sheetCharges, context);
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Return precomputed per-node control-volume areas.
inline std::vector<Real> computeNodeVolumes(const DeviceMesh& mesh)
{
    const Index N = mesh.numNodes();
    std::vector<Real> vol(N, 0.0);
    for (Index i = 0; i < N; ++i)
        vol[i] = mesh.getNode(i).volume;
    return vol;
}

/// Return precomputed per-edge box coupling lengths.
inline std::vector<Real> computeEdgeCouplings(const DeviceMesh& mesh)
{
    const Index E = mesh.numEdges();
    std::vector<Real> couple(E, 0.0);
    for (Index e = 0; e < E; ++e)
        couple[e] = mesh.getEdge(e).couple;
    return couple;
}

/// Return a per-node max adjacent-edge electric-field magnitude [V/m].
inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh)
{
    std::vector<Real> maxField(mesh.numNodes(), 0.0);
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        if (edge.length <= 1.0e-30)
            continue;
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real edgeField = std::abs((psi(j) - psi(i)) / edge.length);
        maxField[edge.n0] = std::max(maxField[edge.n0], edgeField);
        maxField[edge.n1] = std::max(maxField[edge.n1], edgeField);
    }
    return maxField;
}

/// Build edge -> adjacent cell ids map.
inline std::vector<std::vector<Index>> buildEdgeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> edgeCells(mesh.numEdges());
    std::unordered_map<Index, Index> pairToEdge;
    const Index N = mesh.numNodes();
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        pairToEdge[edge.n0 * N + edge.n1] = e;
    }
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell = mesh.getCell(c);
        if (cell.node_ids.size() < 3) continue;
        const Index nids[3] = {
            cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        for (int k = 0; k < 3; ++k) {
            Index a = nids[k];
            Index b = nids[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            auto it = pairToEdge.find(a * N + b);
            if (it != pairToEdge.end())
                edgeCells[it->second].push_back(c);
        }
    }
    return edgeCells;
}

// ---------------------------------------------------------------------------
// Material helpers
// ---------------------------------------------------------------------------

/// Return the average value of a material property over edge-adjacent cells.
/// Falls back to @p fallback only when the edge has no adjacent cells.
/// Throws if any adjacent cell references an unknown material.
inline Real edgeAvgMaterialProp(
    const std::vector<Index>& cells,
    const DeviceMesh&          mesh,
    const MaterialDatabase&    matdb,
    Real Material::*           prop,
    Real                       fallback)
{
    if (cells.empty()) return fallback;
    Real sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        sum += matdb.getMaterial(region.material).*prop;
    }
    return sum / static_cast<Real>(cells.size());
}


/// Build one temperature-adjusted material per mesh cell for hot-path reuse.
inline std::vector<Material> buildCellMaterials(const DeviceMesh&       mesh,
                                                const MaterialDatabase& matdb,
                                                Real                    temperature_K)
{
    std::vector<Material> materials;
    materials.reserve(mesh.numCells());
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        materials.push_back(matdb.getMaterial(region.material, temperature_K));
    }
    return materials;
}


inline Real cellCentroidPotential(const DeviceMesh& mesh, const VectorXd& psi, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return std::numeric_limits<Real>::quiet_NaN();

    Real sum = 0.0;
    for (Index nodeId : cell.node_ids)
        sum += psi(static_cast<int>(nodeId));
    return sum / static_cast<Real>(cell.node_ids.size());
}

inline std::pair<Real, Real> cellCentroid(const DeviceMesh& mesh, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};

    Real x = 0.0;
    Real y = 0.0;
    for (Index nodeId : cell.node_ids) {
        const Node& node = mesh.getNode(nodeId);
        x += node.x;
        y += node.y;
    }
    const Real invCount = 1.0 / static_cast<Real>(cell.node_ids.size());
    return {x * invCount, y * invCount};
}

inline Real estimateSurfaceNormalField(const std::vector<Index>& cells,
                                       const DeviceMesh& mesh,
                                       const VectorXd& psi,
                                       Index edgeId,
                                       Index cellId)
{
    const Edge& edge = mesh.getEdge(edgeId);
    if (edge.length <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();

    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    const Real normalX = -(n1.y - n0.y) / edge.length;
    const Real normalY =  (n1.x - n0.x) / edge.length;

    const auto [cx, cy] = cellCentroid(mesh, cellId);
    const Real cellPhi = cellCentroidPotential(mesh, psi, cellId);
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cellPhi))
        return std::numeric_limits<Real>::quiet_NaN();

    Real maxField = std::numeric_limits<Real>::quiet_NaN();
    for (Index otherCellId : cells) {
        if (otherCellId == cellId)
            continue;
        const auto [ox, oy] = cellCentroid(mesh, otherCellId);
        const Real otherPhi = cellCentroidPotential(mesh, psi, otherCellId);
        if (!std::isfinite(ox) || !std::isfinite(oy) || !std::isfinite(otherPhi))
            continue;
        const Real normalDistance = std::abs((ox - cx) * normalX + (oy - cy) * normalY);
        if (normalDistance <= 1.0e-30)
            continue;
        const Real field = std::abs((otherPhi - cellPhi) / normalDistance);
        if (!std::isfinite(maxField) || field > maxField)
            maxField = field;
    }
    if (std::isfinite(maxField))
        return maxField;

    const Real edgePhi = 0.5 * (psi(static_cast<int>(edge.n0)) + psi(static_cast<int>(edge.n1)));
    const Real mx = 0.5 * (n0.x + n1.x);
    const Real my = 0.5 * (n0.y + n1.y);
    const Real normalDistance = std::abs((cx - mx) * normalX + (cy - my) * normalY);
    if (normalDistance <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();
    return std::abs((cellPhi - edgePhi) / normalDistance);
}

/// Return average model mobility [m^2/V/s] for edge @p edgeId.
inline Real edgeMobility(const std::vector<std::vector<Index>>& edgeCells,
                         const DeviceMesh&                       mesh,
                         const DopingModel&                      doping,
                         const MobilityModel&                    mobility,
                         const std::vector<Material>&            cellMaterials,
                         Index                                   edgeId,
                         CarrierType                             carrier,
                         Real                                    electricField,
                         const MobilityModelConfig*              mobilityConfig = nullptr,
                         const VectorXd*                         psi = nullptr)
{
    const auto& cells = edgeCells[edgeId];
    if (cells.empty()) return 0.0;

    const Edge& edge = mesh.getEdge(edgeId);
    const Real netDoping = 0.5 * (doping.netDoping(edge.n0) +
                                  doping.netDoping(edge.n1));

    const bool surfaceEnabled =
        mobilityConfig != nullptr && isSurfaceMobilityModel(*mobilityConfig);
    std::vector<std::string> adjacentRegionNames;
    if (surfaceEnabled) {
        adjacentRegionNames.reserve(cells.size());
        for (Index c : cells)
            adjacentRegionNames.push_back(mesh.getRegion(mesh.getCell(c).region_id).name);
    }

    Real sum = 0.0;
    Index contributingCells = 0;
    for (Index c : cells) {
        const Material& material = cellMaterials.at(static_cast<std::size_t>(c));
        const Real baseMobility = (carrier == CarrierType::Electron) ? material.mun : material.mup;
        if (baseMobility <= 0.0)
            continue;

        // Average only transport-capable cells.  This keeps oxide-only
        // edges pinned while preserving lateral semiconductor transport on
        // edges that lie along a semiconductor/oxide interface. Surface
        // mobility is enabled only on configured regions/interfaces; when no
        // normal-field estimate is available the NaN field disables the surface
        // factor while preserving any high-field velocity saturation.
        const Region& region = mesh.getRegion(mesh.getCell(c).region_id);
        const bool surfaceApplies = surfaceEnabled &&
            surfaceMobilityAppliesToRegionPair(*mobilityConfig, region.name, adjacentRegionNames);
        const Real surfaceNormalField = (surfaceApplies && psi != nullptr)
            ? estimateSurfaceNormalField(cells, mesh, *psi, edgeId, c)
            : std::numeric_limits<Real>::quiet_NaN();
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, netDoping, 0.0, 0.0, electricField, surfaceNormalField)
            : mobility.holeMobility(
                material, netDoping, 0.0, 0.0, electricField, surfaceNormalField);
        if (modelMobility <= 0.0)
            continue;

        sum += modelMobility;
        ++contributingCells;
    }
    if (contributingCells == 0)
        return 0.0;
    return sum / static_cast<Real>(contributingCells);
}

/// Return average dielectric constant [F/m] for edge @p edgeId.
inline Real edgeEpsilon(const std::vector<std::vector<Index>>& edgeCells,
                        const DeviceMesh&                       mesh,
                        const MaterialDatabase&                 matdb,
                        Index                                   edgeId)
{
    return edgeAvgMaterialProp(edgeCells[edgeId], mesh, matdb,
                               &Material::eps_r, 1.0) * constants::eps0;
}

// ---------------------------------------------------------------------------
// Per-node ni vector
// ---------------------------------------------------------------------------

/// Build per-node intrinsic concentration vector from the material database.
/// For interface nodes shared by multiple regions, uses the first-found value.
inline std::vector<Real> buildNodeNi(const DeviceMesh&       mesh,
                                     const MaterialDatabase& matdb,
                                     Real                    temperature_K = constants::T0)
{
    const Index N = mesh.numNodes();
    std::vector<Real> ni_v(N, 0.0);
    std::vector<bool> found(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Real ni_mat = matdb.getMaterial(region.material, temperature_K).ni;
        for (Index nid : cell.node_ids) {
            if (!found[nid]) {
                ni_v[nid]  = ni_mat;
                found[nid] = true;
            }
        }
    }
    return ni_v;
}

/// Validate that the doping model has one entry per mesh node.
inline void validateDopingMeshSize(const DeviceMesh& mesh,
                                   const DopingModel& doping,
                                   const std::string& context)
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            context + ": doping model size does not match mesh node count.");
}

/// Build per-node effective intrinsic concentration including bandgap narrowing.
inline std::vector<Real> buildEffectiveNodeNi(const DeviceMesh&       mesh,
                                              const MaterialDatabase& matdb,
                                              const DopingModel&      doping,
                                              const BandgapNarrowing& bgn,
                                              Real                    thermalVoltage)
{
    const Real temperature_K = thermalVoltage * constants::q / constants::kb;
    std::vector<Real> ni_v = buildNodeNi(mesh, matdb, temperature_K);
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        const Real totalImpurity = doping.donors(i) + doping.acceptors(i);
        const Real delta = bgn.deltaEg(totalImpurity, 0.0, 0.0);
        ni_v[i] = effectiveIntrinsicDensity(ni_v[i], thermalVoltage, delta);
    }
    return ni_v;
}

/// Validate inputs before building effective intrinsic concentrations.
inline std::vector<Real> buildValidatedEffectiveNodeNi(
    const std::string&             context,
    const DeviceMesh&              mesh,
    const MaterialDatabase&        matdb,
    const DopingModel&             doping,
    const BandgapNarrowingConfig&  bandgapNarrowingConfig,
    Real                           thermalVoltage)
{
    validateDopingMeshSize(mesh, doping, context);
    return buildEffectiveNodeNi(
        mesh,
        matdb,
        doping,
        *makeBandgapNarrowingModel(bandgapNarrowingConfig),
        thermalVoltage);
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

/**
 * @brief Apply strong Dirichlet BCs (row-replacement) to a sparse system.
 *
 * For each constrained node i with prescribed value v:
 *   1. Reduce rhs(k) by A(k,i)*v for all free rows k.
 *   2. Zero column i and row i.
 *   3. Set A(i,i) = 1, rhs(i) = v.
 *
 * The diagonal A(i,i) = 1 is always set explicitly via coeffRef so that
 * the constraint is enforced even if node i had no prior stiffness entries
 * (e.g. insulator nodes in continuity assemblies where all adjacent edges
 * were skipped due to zero mobility).
 */
inline void applyDirichletBC(SparseMatrixd&                         A,
                              VectorXd&                              b,
                              const std::unordered_map<Index, Real>& bcs)
{
    A.makeCompressed();

    // Step 1: propagate prescribed values into free-node RHS
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        for (SparseMatrixd::InnerIterator it(A, i); it; ++it) {
            const int k = static_cast<int>(it.row());
            if (k == i) continue;
            if (bcs.count(static_cast<Index>(k)) == 0)
                b(k) -= it.value() * value;
        }
    }

    // Step 2 & 3: zero Dirichlet rows/cols
    std::unordered_set<int> dirichletSet;
    for (const auto& [nodeId, _] : bcs)
        dirichletSet.insert(static_cast<int>(nodeId));

    for (int col = 0; col < A.outerSize(); ++col) {
        const bool colIsDirichlet = dirichletSet.count(col) > 0;
        for (SparseMatrixd::InnerIterator it(A, col); it; ++it) {
            const int  row          = static_cast<int>(it.row());
            const bool rowIsDirichlet = dirichletSet.count(row) > 0;
            if (rowIsDirichlet || colIsDirichlet)
                it.valueRef() = (row == col) ? 1.0 : 0.0;
        }
    }
    A.prune(0.0);

    // Explicitly ensure A(i,i) = 1 even if the node had no prior entries.
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        A.coeffRef(i, i) = 1.0;
        b(i) = value;
    }
}

} // namespace vela::detail
