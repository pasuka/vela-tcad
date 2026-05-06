#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include <unordered_map>
#include <vector>

namespace vela {

/**
 * @brief FVM/Box assembler for the electrostatic Poisson equation.
 *
 * Discretises:
 *   -div(eps * grad(psi)) = q * (p - n + Nd - Na)
 *
 * with n = p = 0 in the current (linear) stage, reducing to:
 *   -div(eps * grad(psi)) = q * NetDoping
 *
 * Edge-flux formulation (Box method)
 * ------------------------------------
 * For each mesh edge connecting nodes i and j:
 *
 *   flux_ij = eps_ij * (couple_ij / h_ij) * (psi_j - psi_i)
 *
 * where h_ij is the edge length and couple_ij is the Voronoi coupling
 * length (initially approximated as h_ij; replace with real Voronoi
 * circumcenter distance once that is computed).
 *
 * Node volumes
 * ------------
 * Approximated as one-third of the sum of the areas of all triangles
 * that share the node.  Replace with real Voronoi dual-cell areas later.
 *
 * Usage
 * -----
 * @code
 *   PoissonAssembler asm(mesh, matdb, doping);
 *   asm.assemble();
 *   asm.applyDirichlet({{ nodeId, biasValue }, ...});
 *   auto psi = LinearSolver().solve(asm.matrix(), asm.rhs());
 * @endcode
 */
class PoissonAssembler {
public:
    PoissonAssembler(const DeviceMesh&      mesh,
                     const MaterialDatabase& matdb,
                     const DopingModel&      doping);

    /**
     * @brief Assemble the global stiffness matrix and RHS.
     *
     * May be called multiple times (e.g. after changing doping); each call
     * resets and rebuilds both matrix and RHS.
     */
    void assemble();

    /**
     * @brief Enforce Dirichlet boundary conditions (strong row-replacement).
     *
     * @param bcs  Map of nodeId → prescribed potential value [V].
     *
     * For each constrained node i with value v:
     *  - rhs of every free neighbour k is reduced by A(k,i)*v
     *  - column i and row i are zeroed
     *  - A(i,i) = 1, rhs(i) = v
     *
     * Must be called after assemble().
     */
    void applyDirichlet(const std::unordered_map<Index, Real>& bcs);

    const SparseMatrixd& matrix() const { return A_; }
    const VectorXd&      rhs()    const { return b_; }

private:
    // ---- geometry helpers ----

    /// Compute per-node control-volume areas (1/3 of adjacent triangle areas).
    std::vector<Real> computeNodeVolumes() const;

    /// Compute per-edge Voronoi coupling lengths.
    /// Current approximation: couple = edge_length.
    /// Replace body with circumcenter-distance formula when ready.
    std::vector<Real> computeEdgeCouplings() const;

    /// Return the average dielectric constant [F/m] for a given edge.
    Real edgeEpsilon(Index edgeId) const;

    // ---- data ----
    const DeviceMesh&       mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel&      doping_;

    SparseMatrixd A_;
    VectorXd      b_;

    /// edge id → cells that contain this edge
    std::vector<std::vector<Index>> edgeCells_;

    void buildEdgeCellMap();
};

} // namespace vela
