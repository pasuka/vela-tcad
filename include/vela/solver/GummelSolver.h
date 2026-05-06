#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include <string>
#include <unordered_map>

namespace vela {

/**
 * @brief Result of a Gummel iteration.
 */
struct DDSolution {
    VectorXd psi;   ///< Electrostatic potential [V]
    VectorXd phin;  ///< Electron quasi-Fermi potential [V]
    VectorXd phip;  ///< Hole quasi-Fermi potential [V]
    VectorXd n;     ///< Electron concentration [m^-3]
    VectorXd p;     ///< Hole concentration [m^-3]
    int      iters; ///< Number of Gummel iterations performed
};

/**
 * @brief Configuration for the Gummel iteration driver.
 */
struct GummelConfig {
    int    maxIter     = 50;    ///< Maximum number of outer Gummel iterations
    double reltol      = 1.0e-6; ///< Relative convergence tolerance (||Δψ||/||ψ||)
    double dampingPsi  = 1.0;   ///< Damping factor for Poisson update (0 < α ≤ 1)
    double taun        = 1.0e-7; ///< Electron SRH lifetime [s]
    double taup        = 1.0e-7; ///< Hole SRH lifetime [s]
};

/**
 * @brief Gummel self-consistent iteration for steady-state drift-diffusion.
 *
 * Solves the coupled system:
 *   Poisson:            -div(ε·grad(ψ)) = q*(p - n + Nd - Na)
 *   Electron continuity: div(Fn) = R_SRH   (Fn = particle flux density [m⁻²s⁻¹])
 *   Hole continuity:    -div(Fp) = R_SRH   (Fp = particle flux density [m⁻²s⁻¹])
 *
 * Note: The Scharfetter–Gummel helpers return particle flux densities
 * (units [m⁻²s⁻¹]); the factor q is not included in those routines.
 *
 * using decoupled (Gummel) linearisation.
 *
 * Ohmic contact boundary conditions:
 *   ψ_contact = V_bias + ψ_bi   (ψ_bi computed from charge neutrality)
 *   n_contact = n_eq             (charge-neutral equilibrium concentration)
 *   p_contact = ni² / n_eq
 *   phin = phip = V_bias
 *
 * @param mesh           Device mesh.
 * @param matdb          Material database.
 * @param doping         Per-node doping concentrations.
 * @param contactBiases  Map of contact name → applied bias voltage [V].
 * @param cfg            Iteration settings.
 * @return               Converged (or last) solution.
 */
DDSolution runGummel(const DeviceMesh&                         mesh,
                     const MaterialDatabase&                    matdb,
                     const DopingModel&                         doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                         cfg = {});

/**
 * @brief Write a DDSolution to a VTK file.
 *
 * Fields written: Potential, ElectronQuasiFermi, HoleQuasiFermi,
 *                 Electrons, Holes, NetDoping.
 */
void writeDDSolutionVTK(const std::string&    filename,
                        const DeviceMesh&     mesh,
                        const DopingModel&    doping,
                        const DDSolution&     sol);

} // namespace vela
