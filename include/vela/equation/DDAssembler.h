#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include <vector>
#include <unordered_map>

namespace vela {

/**
 * @brief FVM assembler for the steady-state drift-diffusion equations.
 *
 * Provides three assembly routines used by the Gummel iteration:
 *
 *  1. assemblePoissonWithCarriers(n, p)
 *     Assembles the linearised electrostatic Poisson equation:
 *       -div(ε·grad(ψ)) + q*(n+p)/Vt·ψ = q*(p-n+Nd-Na) + q*(n+p)/Vt·ψ_old
 *     including free-carrier charge.
 *
 *  2. assembleElectronContinuity(psi, n_old, p_old)
 *     Assembles the Scharfetter-Gummel electron continuity equation
 *     (solving for n) with SRH recombination linearised w.r.t. n.
 *
 *  3. assembleHoleContinuity(psi, n_old, p_old)
 *     Assembles the SG hole continuity equation (solving for p) with
 *     SRH recombination linearised w.r.t. p.
 *
 * Dirichlet boundary conditions are applied via applyDirichlet().
 *
 * SRH parameters are set at construction time.
 */
class DDAssembler {
public:
    /**
     * @param mesh         The device mesh.
     * @param matdb        Material database (must contain all region materials).
     * @param doping       Per-node donor / acceptor concentrations [m^-3].
     * @param Vt           Thermal voltage kT/q [V].
     * @param taun         Electron SRH lifetime [s] (uniform).
     * @param taup         Hole SRH lifetime [s] (uniform).
     */
    DDAssembler(const DeviceMesh&       mesh,
                const MaterialDatabase& matdb,
                const DopingModel&      doping,
                double                  Vt,
                double                  taun,
                double                  taup);

    // ------------------------------------------------------------------
    // Assembly
    // ------------------------------------------------------------------

    /// Assemble the linearised Poisson equation including carrier charge.
    /// @param n    Current electron concentration per node [m^-3].
    /// @param p    Current hole concentration per node [m^-3].
    /// @param psi  Current potential per node [V] (used for linearisation).
    void assemblePoissonWithCarriers(const VectorXd& n,
                                     const VectorXd& p,
                                     const VectorXd& psi);

    /// Assemble the electron continuity matrix and RHS.
    /// Solves for n_new given ψ, p (from previous Gummel step).
    void assembleElectronContinuity(const VectorXd& psi,
                                    const VectorXd& n_old,
                                    const VectorXd& p_old);

    /// Assemble the hole continuity matrix and RHS.
    /// Solves for p_new given ψ, n (from previous Gummel step).
    void assembleHoleContinuity(const VectorXd& psi,
                                const VectorXd& n_old,
                                const VectorXd& p_old);

    // ------------------------------------------------------------------
    // Dirichlet boundary conditions (strong row-replacement)
    // ------------------------------------------------------------------
    void applyDirichlet(const std::unordered_map<Index, Real>& bcs);

    // ------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------
    const SparseMatrixd& matrix() const { return A_; }
    const VectorXd&      rhs()    const { return b_; }

private:
    const DeviceMesh&       mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel&      doping_;
    double                  Vt_;
    double                  taun_;
    double                  taup_;

    std::vector<Real> ni_; ///< Per-node intrinsic concentration [m^-3]

    SparseMatrixd A_;
    VectorXd      b_;
};

} // namespace vela
