#pragma once

#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

/// Per-contact physics metadata routed to the DD solvers.
///
/// Bias for each contact is still passed through ``contactBiases`` so the
/// existing DC sweep code can update it cheaply.  This auxiliary map carries
/// the parsed ``ContactBoundarySpec`` (Schottky barrier, work function, etc.)
/// for any contact whose type is not the default Ohmic.  Contacts missing
/// from the map fall back to the legacy Ohmic Dirichlet construction.
using ContactSpecsMap = std::unordered_map<std::string, ContactBoundarySpec>;


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
    double augerCn     = 2.8e-43; ///< Electron Auger coefficient [m^6/s]
    double augerCp     = 9.9e-44; ///< Hole Auger coefficient [m^6/s]
    double carrierFloor = 1.0; ///< Minimum solved carrier concentration [m^-3] for quasi-Fermi consistency.
    MobilityModelConfig mobility{}; ///< Mobility model configuration
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    ImpactIonizationModelConfig impactIonization; ///< Avalanche generation model.
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
    UnitScalingConfig inputScaling{}; ///< Input-unit mode from top-level config.
    UnitScalingReferenceConfig unitScalingRefs{}; ///< Optional reference overrides.
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
GummelConfig gummelConfigFromJson(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

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

/// Schottky-aware overloads.  ``contactSpecs`` selects the per-contact
/// boundary model; any contact missing from the map falls back to the
/// legacy Ohmic Dirichlet construction so existing decks keep working.
DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg);

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess);


DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     std::vector<RegionFixedChargeSpec>           fixedCharges,
                     std::vector<InterfaceSheetChargeSpec>        sheetCharges);

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess,
                     std::vector<RegionFixedChargeSpec>           fixedCharges,
                     std::vector<InterfaceSheetChargeSpec>        sheetCharges);


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
