#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include <string>
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
 *   -div(eps * grad(psi)) = q * (NetDoping + fixed_charge)
 *
 * Fixed region charge is supplied as an elementary-charge number density
 * [m^-3]. At most one fixed-charge spec may target a given region. Sheet
 * interface charge is supplied as an elementary-charge number density [m^-2];
 * multiple specs for the same unordered interface pair are summed. Both are
 * multiplied by q during RHS assembly.
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
struct RegionFixedChargeSpec {
    std::string region;          ///< Region name (matches Region::name)
    Real        fixedCharge = 0; ///< Fixed charge density [m^-3], in units of q
};

struct InterfaceSheetChargeSpec {
    std::string region0;          ///< First region name adjacent to the interface
    std::string region1;          ///< Second region name adjacent to the interface
    Real        sheetCharge = 0;  ///< Legacy total sheet charge density [m^-2], in units of q
    Real        fixedCharge = 0;  ///< Fixed interface charge density [m^-2], in units of q
    Real        trapDensity = 0;  ///< Interface trap density [m^-2], in units of q when occupied
    Real        trapOccupancy = 0; ///< Occupied trap fraction [-]

    Real totalSheetCharge() const { return sheetCharge + fixedCharge + trapDensity * trapOccupancy; }
};

/**
 * @brief Neumann boundary condition specification for Poisson equation.
 *
 * Specifies normal displacement D·n [C/m^2] on a boundary segment defined by
 * a polyline of node IDs. The RHS contribution for each edge in the polyline is:
 *
 *   rhs_contribution = normalDisplacement * edge_length / 2
 *
 * distributed equally to the two endpoint nodes.
 *
 * Sign convention:
 *   - Positive normalDisplacement: outward flux (field pointing out of domain)
 *   - Negative normalDisplacement: inward flux (field pointing into domain)
 *   - Zero normalDisplacement: insulating/symmetry boundary
 */
struct PoissonNeumannBoundarySpec {
    std::vector<Index> node_ids;           ///< Polyline defining the boundary segment
    Real               normalDisplacement; ///< Normal displacement D·n [C/m^2]
};

class PoissonAssembler {
public:
    PoissonAssembler(const DeviceMesh&      mesh,
                     const MaterialDatabase& matdb,
                     const DopingModel&      doping,
                     std::vector<RegionFixedChargeSpec> fixedCharges = {},
                     std::vector<InterfaceSheetChargeSpec> sheetCharges = {},
                     std::vector<PoissonNeumannBoundarySpec> neumannBoundaries = {});

    /**
     * @brief Assemble the global stiffness matrix and RHS.
     *
     * May be called multiple times (e.g. after changing doping); each call
     * resets and rebuilds both matrix and RHS.
     *
     * Fixed-charge specs with duplicate region names are rejected.
     *
     * Sheet-charge allocation rule for shared-node interfaces:
     * every mesh edge whose two adjacent cells belong to the requested
     * unordered region pair is treated as an interface segment. The segment charge
     * q * sheet_charge_m2 * edge_length is divided equally between the
     * two endpoint control volumes. This is equivalent to a line-length
     * dual-area split for a conforming shared-node 2-D mesh with unit depth.
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
    std::vector<RegionFixedChargeSpec> fixedCharges_;
    std::vector<InterfaceSheetChargeSpec> sheetCharges_;
    std::vector<PoissonNeumannBoundarySpec> neumannBoundaries_;

    SparseMatrixd A_;
    VectorXd      b_;
};

} // namespace vela
