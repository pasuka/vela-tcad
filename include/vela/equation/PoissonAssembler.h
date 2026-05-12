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
 * where h_ij is the edge length and couple_ij is the precomputed box
 * coupling length from DeviceMesh::buildBoxGeometry().
 *
 * Node volumes
 * ------------
 * Computed as one-third of the sum of the areas of all Tri3 cells that
 * share the node (barycentric control volume).
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
    * @param bcs  Map of nodeId -> prescribed potential value [V].
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
    const DeviceMesh&       mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel&      doping_;

    SparseMatrixd A_;
    VectorXd      b_;
};

} // namespace vela
