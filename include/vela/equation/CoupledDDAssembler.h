#pragma once

#include "vela/core/Types.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/equation/DDAssembler.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/DensityGradientQuantumPotential.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/RecombinationModel.h"
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vela {

struct CoupledDDState {
    VectorXd psi;
    VectorXd phin;
    VectorXd phip;
};

struct CoupledDDThermionicBoundary {
    Real electronVelocity = 0.0;
    Real holeVelocity = 0.0;
    Real electronEquilibriumDensity = 0.0;
    Real holeEquilibriumDensity = 0.0;
    Real boundaryMeasure = 0.0;
};

struct CoupledDDBoundaryConditions {
    std::unordered_map<Index, Real> psi;
    std::unordered_map<Index, Real> phin;
    std::unordered_map<Index, Real> phip;
    std::unordered_map<Index, CoupledDDThermionicBoundary> thermionic;
};

// Diagnostic-only frozen inputs for one-stage feedback substitutions.  The
// production residual/Jacobian path never supplies this object.
struct CoupledDDFeedbackStateSubstitution {
    bool replaceElectronDensity = false;
    bool replaceHoleDensity = false;
    bool replaceElectronQuasiFermi = false;
    bool replaceHoleQuasiFermi = false;
    VectorXd electronDensity;
    VectorXd holeDensity;
    VectorXd electronQuasiFermi_V;
    VectorXd holeQuasiFermi_V;
    // JVP diagnostics: hold each edge's mobility at the unperturbed state.
    bool replaceTransportMobility = false;
    VectorXd electronEdgeMobility;
    VectorXd holeEdgeMobility;
};

struct CoupledDDCarrierTermDiagnostic {
    Index nodeId = 0;
    Real electronDensity_m3 = 0.0;
    Real holeDensity_m3 = 0.0;
    bool electronContinuityActive = true;
    bool holeContinuityActive = true;
    Real electronFlux = 0.0;
    Real holeFlux = 0.0;
    Real electronFluxAbsSum = 0.0;
    Real holeFluxAbsSum = 0.0;
    Real electronRecombination = 0.0;
    Real holeRecombination = 0.0;
    Real electronImpact = 0.0;
    Real holeImpact = 0.0;
    Real impactElectronSource = 0.0;
    Real impactHoleSource = 0.0;
    Real impactCombinedSource = 0.0;
    Real electronGauge = 0.0;
    Real holeGauge = 0.0;
    Real electronBoundary = 0.0;
    Real holeBoundary = 0.0;
    Real electronResidual = 0.0;
    Real holeResidual = 0.0;
};

// Diagnostic-only decomposition of the Poisson row assembled by residual().
// All contribution fields use the production residual's scaled coordinates;
// boundaryReplacement is the amount by which a Dirichlet contact row replaces
// the unconstrained physical equation.
struct CoupledDDPoissonTermDiagnostic {
    Index nodeId = 0;
    Real intrinsicDensity = 0.0;
    Real electronDensityOfStates = 0.0;
    Real holeDensityOfStates = 0.0;
    Real suppliedElectronDensity = 0.0;
    Real suppliedHoleDensity = 0.0;
    Real reconstructedElectronDensity = 0.0;
    Real reconstructedHoleDensity = 0.0;
    Real electronRequiredQuasiFermiPotential = 0.0;
    Real holeRequiredQuasiFermiPotential = 0.0;
    Real electronQuasiFermiMappingError = 0.0;
    Real holeQuasiFermiMappingError = 0.0;
    bool hasSuppliedCarrierState = false;
    Real dielectricFlux = 0.0;
    Real electronCharge = 0.0;
    Real holeCharge = 0.0;
    Real dopingCharge = 0.0;
    Real fixedInterfaceCharge = 0.0;
    Real unconstrainedResidual = 0.0;
    Real boundaryReplacement = 0.0;
    Real productionResidual = 0.0;
    Real closureError = 0.0;
    bool constrained = false;
};

struct CoupledDDEdgeFluxDiagnostic {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real x0 = 0.0;
    Real y0 = 0.0;
    Real x1 = 0.0;
    Real y1 = 0.0;
    Real length_m = 0.0;
    Real couple_m = 0.0;
    Real netDopingAvg_m3 = 0.0;
    Real ni0_m3 = 0.0;
    Real ni1_m3 = 0.0;
    Real electronDensity0_m3 = 0.0;
    Real electronDensity1_m3 = 0.0;
    Real holeDensity0_m3 = 0.0;
    Real holeDensity1_m3 = 0.0;
    Real psi0_V = 0.0;
    Real psi1_V = 0.0;
    Real phin0_V = 0.0;
    Real phin1_V = 0.0;
    Real phip0_V = 0.0;
    Real phip1_V = 0.0;
    Real electricField_V_m = 0.0;
    Real electronMobilityField_V_m = 0.0;
    Real electronMobility_m2_V_s = 0.0;
    Real holeMobility_m2_V_s = 0.0;
    // Exact dimensionless factors used by the production electron SG branch.
    // For Fermi statistics, eta0/eta1 are reduced chemical potentials and
    // generalizedEinsteinFactor modifies both drift and quasi-Fermi arguments.
    Real electronEta0 = 0.0;
    Real electronEta1 = 0.0;
    Real electronDriftPotential_V = 0.0;
    Real electronGeneralizedEinsteinFactor = 1.0;
    Real electronBernoulliArgument = 0.0;
    Real electronQuasiFermiArgument = 0.0;
    Real electronBernoulliPlus = 1.0;
    Real electronBernoulliMinus = 1.0;
    // Signed Scharfetter-Gummel continuity edge flux (added to node0's residual
    // and subtracted from node1's), identical to the residual edge loop.
    Real electronFlux = 0.0;
    Real holeFlux = 0.0;
    // Diagnostic response to one internal length unit of edge coupling.  This
    // remains defined when the production coupling is zero, allowing a direct
    // fixed-state replay of an external box-coefficient oracle.
    Real electronFluxPerInternalCouple = 0.0;
    Real holeFluxPerInternalCouple = 0.0;
    // The same integrated particle line flux before residual nondimensionalization
    // (particles per metre of out-of-plane depth per second).
    Real electronParticleLineFlux_per_m_s = 0.0;
    Real holeParticleLineFlux_per_m_s = 0.0;
    Real electronParticleLineFluxPerInternalCouple_per_m_s = 0.0;
    Real holeParticleLineFluxPerInternalCouple_per_m_s = 0.0;
};

// Diagnostic-only decomposition of one transport edge derivative with respect
// to one same-carrier quasi-Fermi endpoint.  The record is emitted once for
// each endpoint row so contact-row replacement can be audited independently
// from the physical edge operator.
struct CoupledDDTransportEdgeJacobianDiagnostic {
    Index edgeId = 0;
    std::string carrier;
    Index node0 = 0;
    Index node1 = 0;
    Index rowNode = 0;
    Index columnNode = 0;
    int rowEndpoint = 0;
    int columnEndpoint = 0;
    Real lengthInternal = 0.0;
    Real coupleInternal = 0.0;
    Real qfpDriveInternal = 0.0;
    Real dQfpDriveDColumnInternal = 0.0;
    Real mobilityInternal = 0.0;
    Real dMobilityDColumnInternal = 0.0;
    Real bernoulliNode0 = 0.0;
    Real bernoulliNode1 = 0.0;
    Real carrierDensityNode0Internal = 0.0;
    Real carrierDensityNode1Internal = 0.0;
    Real fluxPhysical = 0.0;
    Real fluxScaled = 0.0;
    Real rowSign = 0.0;
    Real productionFrozenMobilityDerivativePhysical = 0.0;
    Real frozenMobilityFiniteDifferenceDerivativePhysical = 0.0;
    Real liveMobilityFiniteDifferenceDerivativePhysical = 0.0;
    Real mobilityResponseDerivativePhysical = 0.0;
    Real liveMinusFrozenFiniteDifferenceDerivativePhysical = 0.0;
    Real bernoulliQfpDerivativePhysical = 0.0;
    Real carrierPopulationDerivativePhysical = 0.0;
    Real productionRowDerivativeScaled = 0.0;
    Real liveMobilityRowDerivativeScaled = 0.0;
    bool rowConstrained = false;
    bool columnConstrained = false;
    Real contactEliminatedProductionEdgeDerivative = 0.0;
    Real contactIdentityEntry = 0.0;
    Real continuityRowWeight = 1.0;
    Real solverProductionEdgeDerivative = 0.0;
    Real solverContactIdentityEntry = 0.0;
};

struct CarrierDiagonalFloorRegularizationConfig {
    bool enabled = false;
    Real scale = 1.0;
    Real minorityDensityRatio = 1.0;
};

class CoupledDDAssembler {
    struct FixedJacobianData;
public:
    /// Single-owner sequential cache. Stores no numeric Jacobian or state.
    struct StructureCache {
    private:
        friend class CoupledDDAssembler;
        std::vector<std::uint64_t> key;
        std::shared_ptr<const FixedJacobianData> data;
    };
    void setStructureCache(std::shared_ptr<StructureCache> cache) {
        structureCache_ = std::move(cache);
    }
    void setFermiNodeCache(bool enabled) { fermiNodeCacheEnabled_ = enabled; }
    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       double taun,
                       double taup);

    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       double taun,
                       double taup,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       const MobilityModelConfig& mobilityConfig,
                       const RecombinationModelConfig& recombinationConfig,
                       const BandgapNarrowingConfig& bandgapNarrowingConfig = {},
                       const ImpactIonizationModelConfig& impactIonizationConfig = {},
                       std::vector<RegionFixedChargeSpec> fixedCharges = {},
                       std::vector<InterfaceSheetChargeSpec> sheetCharges = {},
                       DDScalingSpec scaling = {},
                       CarrierDiagonalFloorRegularizationConfig carrierDiagonalFloor = {},
                       CarrierStatisticsConfig carrierStatistics = {},
                       DensityGradientQuantumPotentialConfig electronQuantumPotential = {});

    VectorXd pack(const CoupledDDState& state) const;
    CoupledDDState unpack(const VectorXd& x) const;

    /// Store quasi-Fermi unknowns as contact-referenced increments.  The
    /// public CoupledDDState and boundary-condition values remain scaled
    /// absolute potentials; only the internal packed coordinates are shifted.
    void setQuasiFermiReferences(Real electronReference_V,
                                 Real holeReference_V);
    /// Store independent quasi-Fermi reference coordinates at every node.
    /// References are constant coordinate shifts; physical fields and the
    /// assembled equations remain invariant under a consistent repartition.
    void setQuasiFermiReferenceFields(const VectorXd& electronReference_V,
                                     const VectorXd& holeReference_V);
    /// Repartition reference plus increment without discarding the low part.
    /// Callers must invalidate any external caches keyed only by coordinates.
    void recenterQuasiFermiState(VectorXd& state);
    Real electronQuasiFermiReference() const { return electronQfReference_V_; }
    Real holeQuasiFermiReference() const { return holeQfReference_V_; }
    Real electronQuasiFermiReferenceAt(Index node) const;
    Real holeQuasiFermiReferenceAt(Index node) const;
    VectorXd electronQuasiFermiReferenceField() const;
    VectorXd holeQuasiFermiReferenceField() const;

    /// Set a frozen electron quantum correction [V] for one outer
    /// density-gradient iteration.  The coupled unknown count remains 3*N.
    void setElectronQuantumPotential(const VectorXd& potential_V);
    const VectorXd& electronQuantumPotential() const
    {
        return electronQuantumPotential_V_;
    }
    bool electronQuantumPotentialEnabled() const
    {
        return electronQuantumPotentialConfig_.enabled;
    }

    VectorXd residual(const VectorXd& x,
                      const CoupledDDBoundaryConditions& bcs) const;

    VectorXd feedbackSubstitutionResidual(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        const CoupledDDFeedbackStateSubstitution& substitution) const;

    std::vector<CoupledDDPoissonTermDiagnostic> poissonTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    std::vector<CoupledDDPoissonTermDiagnostic> poissonTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        const VectorXd& suppliedElectronDensity,
        const VectorXd& suppliedHoleDensity) const;

    SparseMatrixd assembleJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;

    SparseMatrixd finiteDifferenceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep = 1.0e-6) const;

    // Assemble only the impact-ionization source contribution. These direct
    // paths avoid subtracting two full continuity residuals/Jacobians when the
    // avalanche block is many orders of magnitude smaller than transport.
    VectorXd impactIonizationSourceResidual(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    SparseMatrixd impactIonizationSourceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    SparseMatrixd impactIonizationSourceFiniteDifferenceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep) const;
    // Independent high-precision reference for branch-sensitive source audits.
    // Callers must choose a step below every nonzero active-branch margin.
    // Exact-zero abs/norm branches use the symmetric semismooth derivative.
    SparseMatrixd impactIonizationSourceBranchResolvedFiniteDifferenceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep) const;

    VectorXd electronDensity(const VectorXd& x) const;
    VectorXd holeDensity(const VectorXd& x) const;
    std::vector<CoupledDDCarrierTermDiagnostic> carrierContinuityTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    // Terms in the continuity equations actually solved by Newton.  In
    // postprocess_only mode the observable impact-ionization source is omitted
    // so convergence checks and row scaling remain identical to avalanche-off.
    std::vector<CoupledDDCarrierTermDiagnostic>
    carrierContinuityEquationTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    std::vector<CoupledDDCarrierTermDiagnostic>
    feedbackSubstitutionCarrierContinuityTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        const CoupledDDFeedbackStateSubstitution& substitution) const;

    std::vector<CoupledDDEdgeFluxDiagnostic> sgEdgeFluxDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;
    std::vector<CoupledDDTransportEdgeJacobianDiagnostic>
    transportEdgeJacobianDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real physicalFiniteDifferenceStep_V = 1.0e-7) const;

    bool hasPositiveFiniteCarriers(const VectorXd& x) const;
    Index numNodes() const { return mesh_.numNodes(); }
    const std::vector<Real>& intrinsicDensity() const { return ni_; }
    const std::vector<Real>& electronDensityOfStates() const { return Nc_; }
    const std::vector<Real>& holeDensityOfStates() const { return Nv_; }
    const CarrierStatisticsConfig& carrierStatistics() const { return carrierStatistics_; }
    bool usesScaledState() const { return scaling_.enabled; }
    Real potentialScale() const { return scaling_.V0; }
    Real concentrationScale() const { return scaling_.C0; }
    Real continuityResidualScale() const { return scaling_.enabled ? scaling_.C0 * scaling_.D0 : 1.0; }
    std::string impactIonizationConfigurationFingerprint() const;
    std::string impactIonizationActiveBranchFingerprint(const VectorXd& x) const;

private:
    VectorXd residualImpl(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        const CoupledDDFeedbackStateSubstitution* substitution) const;
    std::vector<CoupledDDCarrierTermDiagnostic>
    carrierContinuityTermDiagnosticsImpl(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        const CoupledDDFeedbackStateSubstitution* substitution,
        bool includeImpactIonization) const;

    template <typename Scalar>
    SparseMatrixd impactIonizationSourceFiniteDifferenceJacobianImpl(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep) const;

    int psiOffset() const { return 0; }
    int phinOffset() const { return static_cast<int>(mesh_.numNodes()); }
    int phipOffset() const { return 2 * static_cast<int>(mesh_.numNodes()); }

    using JacobianStorageIndex = SparseMatrixd::StorageIndex;
    static constexpr JacobianStorageIndex invalidJacobianOffset = -1;
    static constexpr std::size_t maxEdgeAvalancheStencilNodes = 8;
    static constexpr std::size_t maxEdgeAvalancheDerivativeColumns =
        2 + 2 * maxEdgeAvalancheStencilNodes;
    enum class MobilityDopingBasis : std::uint8_t {
        NetDoping,
        TotalImpurity,
        CellReconstructedTotalImpurity,
    };
    struct EdgeAssemblyKernel {
        Index n0 = 0;
        Index n1 = 0;
        Real length = 0.0;
        Real coupling = 0.0;
        Real poissonCoupling = 0.0;
        Point2 tangent = Point2::Zero();
        Real avalancheSourceArea = 0.0;
        Real electronLogNiNc0 = 0.0;
        Real electronLogNiNc1 = 0.0;
        Real electronDriftOffset = 0.0;
        Real holeLogNiNv0 = 0.0;
        Real holeLogNiNv1 = 0.0;
        Real holeDriftOffset = 0.0;
        bool activeTransport = false;
        std::array<Index, maxEdgeAvalancheStencilNodes> avalancheStencilNodes{};
        std::uint8_t avalancheStencilNodeCount = 0;
        bool touchesTransportContact = false;
        Real vectorGradientArea = 0.0;
        std::array<Point2, maxEdgeAvalancheStencilNodes> vectorGradientSensitivities{};
        std::vector<Real> electronLowFieldMobilities;
        std::vector<Real> holeLowFieldMobilities;
    };
    struct NodalCurrentReconstructionTerm {
        Index edgeId = 0;
        Point2 tangent = Point2::Zero();
        Real weight = 0.0;
    };
    struct NodalCurrentReconstructionKernel {
        std::vector<NodalCurrentReconstructionTerm> terms;
        Real a00 = 0.0;
        Real a01 = 0.0;
        Real a11 = 0.0;
        Real determinant = 0.0;
        Real fallbackWeight = 0.0;
        bool useLeastSquares = false;
    };

    void buildEdgeAssemblyKernels();
    void buildNodalCurrentReconstructionKernels();
    Real cachedEdgeMobility(
        Index edgeId,
        CarrierType carrier,
        Real drivingField,
        const VectorXd* psi) const;
    bool usesSentaurusExponentialQuantumCoupling() const;
    Real electronTransportPotential(Index node, Real psiRelative_V) const;
    Real electronDensityAt(Index node,
                           Real psiRelative_V,
                           Real phinIncrement_V) const;
    Real electronDensityDerivativeEtaAt(Index node,
                                        Real psiRelative_V,
                                        Real phinIncrement_V) const;
    bool usesSentaurusDefaultSrhDensityCoupling() const;
    Real electronSrhDensityAt(Index node,
                              Real psiRelative_V,
                              Real phinIncrement_V) const;
    Real electronSrhDensityDerivativeEtaAt(Index node,
                                           Real psiRelative_V,
                                           Real phinIncrement_V) const;
    Real nodeRecombinationRate(Index node,
                               Real electronTransportDensity,
                               Real electronSrhDensity,
                               Real holeDensity,
                               Real quasiFermiSplitting_V,
                               const GeneralizedSrhCarrierState* srhState = nullptr,
                               const GeneralizedSrhCarrierState* augerState = nullptr) const;
    void rebuildFixedJacobianPattern(
        const std::vector<bool>& constrainedRows,
        bool includeCellStencil,
        std::uint64_t boundarySignature) const;
    JacobianStorageIndex fixedJacobianOffset(int row, int col) const;
    std::vector<std::uint64_t> jacobianStructureKey(
        const std::vector<bool>& constrainedRows, bool includeCellStencil) const;

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    double Vt_;
    void refreshIalMobility(const VectorXd& x) const;
    mutable MobilityModelConfig mobilityConfig_;
    std::unique_ptr<MobilityModel> mobility_;
    RecombinationModelConfig recombinationConfig_;
    RecombinationModel recombination_;
    BandgapNarrowingConfig bandgapNarrowingConfig_;
    ImpactIonizationModelConfig impactIonizationConfig_;
    std::unique_ptr<ImpactIonizationModel> impactIonization_;
    bool impactIonizationEnabled_ = false;
    bool impactIonizationCoupled_ = false;
    bool bgnEnabled_ = false;
    CarrierStatisticsConfig carrierStatistics_;
    DensityGradientQuantumPotentialConfig electronQuantumPotentialConfig_;
    VectorXd electronQuantumPotential_V_;
    /// Resolved once in the constructor; hot loops use this instead of the
    /// string-valued configuration. Numerically identical: only the model lookup is hoisted.
    CarrierStatisticsModel carrierStatisticsModel_ = CarrierStatisticsModel::Boltzmann;
    bool usesFermiDirac_ = false;
    std::vector<Real> ni_;
    std::vector<Real> Nc_;
    std::vector<Real> Nv_;
    std::vector<Material> cellMaterials_;

    // Mesh-derived quantities cached at construction time.
    std::vector<std::vector<Index>> edgeCells_;
    std::vector<std::vector<Index>> nodeCells_;
    std::vector<std::vector<Index>> cellEdges_;
    std::vector<std::vector<Index>> nodeEdges_;
    std::vector<bool> contactNodes_;
    std::vector<EdgeAssemblyKernel> edgeAssemblyKernels_;
    std::vector<NodalCurrentReconstructionKernel>
        nodalCurrentReconstructionKernels_;
    std::vector<Real> vol_;
    std::vector<Real> poissonChargeVol_;
    std::vector<Real> poissonCouple_;
    std::vector<Real> couple_;
    VectorXd fixedInterfaceChargeRhs_;
    DDScalingSpec scaling_;
    CarrierDiagonalFloorRegularizationConfig carrierDiagonalFloor_;
    Real electronQfReference_V_ = 0.0;
    Real holeQfReference_V_ = 0.0;
    VectorXd electronQfReferenceField_V_;
    VectorXd holeQfReferenceField_V_;
    MobilityDopingBasis mobilityDopingBasis_ = MobilityDopingBasis::NetDoping;
    bool surfaceMobilityEnabled_ = false;
    bool highFieldMobilityEnabled_ = false;
    bool qfMobilityEnabled_ = false;
    bool vectorQfMobilityEnabled_ = false;
    bool transportMobilityDerivativeEnabled_ = false;
    mutable bool hasObservedJacobianPattern_ = false;
    mutable std::uint64_t lastObservedJacobianPattern_ = 0;
    mutable bool hasFixedJacobianPattern_ = false;
    mutable std::uint64_t fixedJacobianBoundarySignature_ = 0;
    mutable bool fixedJacobianIncludesCellStencil_ = false;
    struct FixedJacobianData {
        SparseMatrixd pattern;
        std::unordered_map<std::uint64_t, JacobianStorageIndex> offsets;
        std::vector<std::array<JacobianStorageIndex, 36>> edgeScatter;
        std::vector<std::array<JacobianStorageIndex, 4 * maxEdgeAvalancheDerivativeColumns>> avalancheScatter;
        std::vector<std::array<JacobianStorageIndex, 9>> nodeScatter;
        std::vector<std::array<JacobianStorageIndex, 81>> cellScatter;
        std::vector<bool> constrainedRows;
    };
    mutable std::shared_ptr<const FixedJacobianData> fixedJacobian_;
    std::shared_ptr<StructureCache> structureCache_;
    mutable std::vector<std::uint64_t> structureBaseKey_;
    bool fermiNodeCacheEnabled_ = false;
    mutable bool hasCachedActiveBranchFingerprint_ = false;
    mutable VectorXd cachedActiveBranchFingerprintState_;
    mutable std::string cachedActiveBranchFingerprint_;
};

} // namespace vela
