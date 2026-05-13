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
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include <Eigen/Sparse>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <stdexcept>

namespace vela::detail {

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


/// Return average model mobility [m^2/V/s] for edge @p edgeId.
inline Real edgeMobility(const std::vector<std::vector<Index>>& edgeCells,
                         const DeviceMesh&                       mesh,
                         const MaterialDatabase&                 matdb,
                         const DopingModel&                      doping,
                         const MobilityModel&                    mobility,
                         Index                                   edgeId,
                         CarrierType                             carrier)
{
    const auto& cells = edgeCells[edgeId];
    if (cells.empty()) return 0.0;

    const Edge& edge = mesh.getEdge(edgeId);
    const Real netDoping = 0.5 * (doping.netDoping(edge.n0) +
                                  doping.netDoping(edge.n1));

    Real sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        const Material& material = matdb.getMaterial(region.material);
        if (carrier == CarrierType::Electron)
            sum += mobility.electronMobility(material, netDoping, 0.0, 0.0);
        else
            sum += mobility.holeMobility(material, netDoping, 0.0, 0.0);
    }
    return sum / static_cast<Real>(cells.size());
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
                                     const MaterialDatabase& matdb)
{
    const Index N = mesh.numNodes();
    std::vector<Real> ni_v(N, 0.0);
    std::vector<bool> found(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Real ni_mat = matdb.getMaterial(region.material).ni;
        for (Index nid : cell.node_ids) {
            if (!found[nid]) {
                ni_v[nid]  = ni_mat;
                found[nid] = true;
            }
        }
    }
    return ni_v;
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
