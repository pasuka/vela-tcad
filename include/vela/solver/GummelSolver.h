#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <unordered_map>
#include <vector>

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
    int      iters = 0; ///< Number of Gummel iterations performed
    bool     converged = false; ///< True if the Gummel convergence criteria were met
};

/**
 * @brief Configuration for the Gummel iteration driver.
 */
struct GummelConfig {
    int    maxIter     = 50;    ///< Maximum number of outer Gummel iterations
    double reltol      = 1.0e-6; ///< Relative convergence tolerance (||dpsi||/||psi||)
    double abstol      = 0.0;   ///< Absolute update tolerance across psi, n, and p
    double temperature_K = constants::T0; ///< Lattice temperature [K]
    double dampingPsi  = 1.0;   ///< Damping factor for Poisson update (0 < alpha <= 1)
    double taun        = 1.0e-7; ///< Electron SRH lifetime [s]
    double taup        = 1.0e-7; ///< Hole SRH lifetime [s]
    std::string mobility = "constant"; ///< "constant" or "caughey_thomas"
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
};

/**
 * @brief Gummel self-consistent iteration for steady-state drift-diffusion.
 *
 * Solves the coupled system:
 *   Poisson:            -div(eps*grad(psi)) = q*(p - n + Nd - Na)
 *   Electron continuity: div(Fn) = R_SRH   (Fn = particle flux density [m^-2 s^-1])
 *   Hole continuity:    -div(Fp) = R_SRH   (Fp = particle flux density [m^-2 s^-1])
 *
 * Note: The Scharfetter-Gummel helpers return particle flux densities
 * (units [m^-2 s^-1]); the factor q is not included in those routines.
 *
 * using decoupled (Gummel) linearisation.
 *
 * Ohmic contact boundary conditions:
 *   psi_contact = V_bias + psi_bi   (psi_bi computed from charge neutrality)
 *   n_contact = n_eq             (charge-neutral equilibrium concentration)
 *   p_contact = ni^2 / n_eq
 *   phin = phip = V_bias
 *
 * @param mesh           Device mesh.
 * @param matdb          Material database.
 * @param doping         Per-node doping concentrations.
 * @param contactBiases  Map of contact name -> applied bias voltage [V].
 * @param cfg            Iteration settings.
 * @return               Converged (or last) solution.
 */
GummelConfig gummelConfigFromJson(const nlohmann::json& cfg);

DDSolution runGummel(const DeviceMesh&                         mesh,
                     const MaterialDatabase&                    matdb,
                     const DopingModel&                         doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                         cfg = {});

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess);

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
