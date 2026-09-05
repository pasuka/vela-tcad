#pragma once

#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/equation/PoissonChargeVolume.h"
#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/BandToBandTunnelingModel.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/RecombinationModel.h"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

namespace vela {

struct CoupledDDEdgeFluxDiagnostic;

/// Per-contact physics metadata routed to the DD solvers.
///
/// Bias for each contact is still passed through ``contactBiases`` so the
/// existing DC sweep code can update it cheaply.  This auxiliary map carries
/// the parsed ``ContactBoundarySpec`` (metal-gate role, Schottky barrier, etc.)
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
    VectorXd phinIncrement; ///< Optional cancellation-free electron QF increment [V].
    VectorXd phipIncrement; ///< Optional cancellation-free hole QF increment [V].
    VectorXd electronQuantumPotential; ///< Electron density-gradient correction [V].
    /// Continuous Phi/q state for all-material density-gradient restart [V].
    /// Lambda is reconstructed per material from this field.
    VectorXd electronQuantumPotentialLike;
    Real electronQfReference_V = 0.0; ///< Reference for phinIncrement.
    Real holeQfReference_V = 0.0; ///< Reference for phipIncrement.
    /// Optional per-node references for cancellation-free quasi-Fermi states.
    /// When absent, the scalar references above apply to every node.
    VectorXd electronQfReference;
    VectorXd holeQfReference;
    /// Optional original normalized coordinates for extended-Poisson restart.
    /// A hint is usable only while its physical fields and scale still match.
    VectorXd packedState;
    VectorXd packedLow; ///< Optional tails of all three normalized blocks.
    std::string packedMeshFingerprint;
    bool hasConsistentSplitPackedState() const;
    Real splitPotentialAt(int block,int node) const;
    Real splitQuasiFermiDifferenceAt(int node) const;
    Real packedPotentialScale_V = 0.0;
    VectorXd n;     ///< Electron concentration [m^-3]
    VectorXd p;     ///< Hole concentration [m^-3]
    int      iters = 0; ///< Number of Gummel iterations performed
    bool     converged = false; ///< True if the Gummel convergence criteria were met

    bool hasReferencedElectronQuasiFermi() const
    {
        return phinIncrement.size() == phin.size();
    }
    bool hasConsistentPackedState() const
    {
        if(packedLow.size()!=0)return hasConsistentSplitPackedState();
        const int count=static_cast<int>(psi.size());
        if(count==0 || packedState.size()!=3*count ||
           !std::isfinite(packedPotentialScale_V) || packedPotentialScale_V<=0. ||
           !packedState.allFinite() || !hasReferencedElectronQuasiFermi() ||
           !hasReferencedHoleQuasiFermi() || phin.size()!=count || phip.size()!=count)
            return false;
        // Physical CSV values historically canonicalize subnormals to zero.
        // Do not emit a competing exact-coordinate hint for that case.
        const auto normal=[](Real z){return z==0. || std::isnormal(z);};
        for(int i=0;i<count;++i) {
            const Real en=electronQuasiFermiAt(i),hp=holeQuasiFermiAt(i);
            const Real tol=32.*std::numeric_limits<Real>::epsilon();
            if(!std::isfinite(en) || !std::isfinite(hp) ||
               std::abs(phin(i)-en)>tol*std::max({Real{1.},std::abs(phin(i)),std::abs(en)}) ||
               std::abs(phip(i)-hp)>tol*std::max({Real{1.},std::abs(phip(i)),std::abs(hp)}))return false;
            if(!normal(psi(i)) || !normal(phinIncrement(i)) || !normal(phipIncrement(i)))return false;
            for(int b=0;b<3;++b)if(!normal(packedState(b*count+i)))return false;
            if(psi(i)!=packedState(i)*packedPotentialScale_V ||
               phinIncrement(i)!=packedState(count+i)*packedPotentialScale_V ||
               phipIncrement(i)!=packedState(2*count+i)*packedPotentialScale_V)return false;
        }
        return true;
    }
    bool hasReferencedHoleQuasiFermi() const
    {
        return phipIncrement.size() == phip.size();
    }
    Real electronQuasiFermiReferenceAt(int node) const
    {
        return electronQfReference.size() == phin.size()
            ? electronQfReference(node) : electronQfReference_V;
    }
    Real holeQuasiFermiReferenceAt(int node) const
    {
        return holeQfReference.size() == phip.size()
            ? holeQfReference(node) : holeQfReference_V;
    }
    Real electronQuasiFermiAt(int node) const
    {
        if(packedLow.size()!=0)return splitPotentialAt(1,node);
        return hasReferencedElectronQuasiFermi()
            ? electronQuasiFermiReferenceAt(node) + phinIncrement(node)
            : phin(node);
    }
    Real holeQuasiFermiAt(int node) const
    {
        if(packedLow.size()!=0)return splitPotentialAt(2,node);
        return hasReferencedHoleQuasiFermi()
            ? holeQuasiFermiReferenceAt(node) + phipIncrement(node)
            : phip(node);
    }
    Real holeMinusElectronQuasiFermiAt(int node) const
    {
        if(packedLow.size()!=0)return splitQuasiFermiDifferenceAt(node);
        if (hasReferencedElectronQuasiFermi() &&
            hasReferencedHoleQuasiFermi()) {
            return static_cast<Real>(
                (static_cast<long double>(holeQuasiFermiReferenceAt(node)) -
                 static_cast<long double>(electronQuasiFermiReferenceAt(node))) +
                (static_cast<long double>(phipIncrement(node)) -
                 static_cast<long double>(phinIncrement(node))));
        }
        return holeQuasiFermiAt(node) - electronQuasiFermiAt(node);
    }
};

/**
 * @brief Configuration for the Gummel iteration driver.
 */
struct GummelConfig {
    int    maxIter     = 50;    ///< Maximum number of outer Gummel iterations
    double reltol      = 1.0e-6; ///< Relative convergence tolerance (||dpsi||/||psi||)
    double abstol      = 0.0;   ///< Absolute update tolerance across psi, n, and p
    double temperature_K = constants::T0; ///< Lattice temperature [K]
    PoissonChargeVolumePolicy poissonChargeVolumePolicy =
        PoissonChargeVolumePolicy::Global;
    double dampingPsi  = 1.0;   ///< Damping factor for Poisson update (0 < alpha <= 1)
    double taun        = 1.0e-5; ///< Electron SRH lifetime [s]
    double taup        = 3.0e-6; ///< Hole SRH lifetime [s]
    double augerCn     = 2.90e-43; ///< Electron Auger coefficient [m^6/s]
    double augerCp     = 1.028e-43; ///< Hole Auger coefficient [m^6/s]
    std::string augerExcessProduct = "generalized_fermi"; ///< "generalized_fermi" or "classical_np".
    bool augerWithGeneration = true; ///< Compatibility default; false clips Auger generation.
    AugerDensityDependenceConfig augerDensityDependence{};
    double carrierFloor = 1.0; ///< Minimum solved carrier concentration [m^-3] for quasi-Fermi consistency.
    MobilityModelConfig mobility{}; ///< Mobility model configuration
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    SRHDopingDependenceConfig srhDopingDependence{}; ///< Sentaurus SRH(DopingDep).
    BandToBandTunnelingConfig bandToBand{}; ///< Local pair generation, including Sentaurus E2.
    ImpactIonizationModelConfig impactIonization; ///< Avalanche generation model.
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
    CarrierStatisticsConfig carrierStatistics; ///< Gummel currently supports Boltzmann only.
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

/// Contact-model-aware overloads.  ``contactSpecs`` selects the per-contact
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
                        const DDSolution&     sol,
                        UnitScalingConfig scaling = {});

void writeDDSolutionVTK(const std::string& filename,
                        const DeviceMesh& mesh,
                        const MaterialDatabase& matdb,
                        const DopingModel& doping,
                        const DDSolution& sol,
                        const MobilityModelConfig& mobilityConfig,
                        const RecombinationModelConfig& recombinationConfig,
                        const ImpactIonizationModelConfig& impactIonizationConfig,
                        const BandgapNarrowingConfig& bandgapNarrowingConfig,
                        Real temperature_K = constants::T0,
                        UnitScalingConfig scaling = {},
                        const CarrierStatisticsConfig& carrierStatistics = {},
                        const std::vector<CoupledDDEdgeFluxDiagnostic>*
                            sgEdgeFluxDiagnostics = nullptr,
                        bool writeCellFirstSgCurrentDiagnostics = false);

} // namespace vela
