#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/numerics/LineSearch.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/BandToBandTunnelingModel.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/simulation/PseudoArclength.h"
#include "vela/solver/GummelSolver.h"
#include <functional>
#include <limits>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

struct NewtonCarrierRowConvergenceConfig {
    std::string mode = "off"; ///< "off", "report", or "enforce".
    Real epsRow = 1.0e-3;
    Real scaleFloor = 1.0e-300;
    Real minSourceScaleFraction = 1.0e-3;
    Real minSourceScale = 0.0;
    Real minSourceGlobalFraction = 0.0; ///< Optional fraction of the maximum carrier source scale; 0 disables.
    Real minCarrierDensity_m3 = 0.0; ///< Optional density qualification; 0 disables.
    Real minFluxScaleFraction = 0.0; ///< Optional fraction of the maximum carrier flux scale; 0 disables.
    Real minFluxScale = 0.0; ///< Optional absolute carrier flux qualification; 0 disables.
    int minEnforceMaxIter = 200;
    std::string diagnosticCsvFile;
    std::string traceCsvFile;
    std::vector<Index> traceNodes;
    int traceFirstIterations = 10;
    int traceEveryIterations = 10;
};

struct NewtonLocalUpdateDiagnosticsConfig {
    bool enabled = false; ///< Write opt-in per-node Newton update and linear-residual diagnostics.
    std::string csvFile;
    std::vector<Index> nodes;
    int firstIterations = 10;
    int everyIterations = 10;
};

struct NewtonContinuityRowScalingConfig {
    bool enabled = false;
    Real fluxFraction = 1.0e-3;
    Real scaleFloor = 1.0e-30;
    Real minSourceScale = 0.0;
    Real minWeight = 1.0e-12;
    Real maxWeight = 1.0e12;
};

struct NewtonLinearEquilibrationConfig {
    std::string mode = "off"; ///< "off" or one-pass "l2_row_column".
};

struct NewtonGlobalContinuityClosureConfig {
    std::string mode = "off"; ///< "off", "report", or "enforce".
    Real tolerance = 1.0e-2;
    Real sourceFloor = 1.0e-12;
};

struct NewtonGlobalContinuityCarrierClosure {
    bool qualified = false;
    Real contactFlux = 0.0;
    Real integratedSource = 0.0;
    Real mismatch = 0.0;
    Real ratio = 0.0;
};

struct NewtonGlobalContinuityClosureEvaluation {
    bool enabled = false;
    bool enforced = false;
    bool satisfied = true;
    NewtonGlobalContinuityCarrierClosure electron;
    NewtonGlobalContinuityCarrierClosure hole;
};

struct NewtonCarrierRowConvergenceViolation {
    Index nodeId = 0;
    std::string carrier;
    Real residual = 0.0;
    Real scale = 0.0;
    Real ratio = 0.0;
    Real carrierDensity_m3 = 0.0;
    Real flux = 0.0;
    Real recombination = 0.0;
    Real impact = 0.0;
    bool densityQualified = false;
    bool fluxQualified = false;
    bool sourceQualified = false;
};

struct NewtonCarrierRowConvergenceEvaluation {
    bool enabled = false;
    bool enforced = false;
    bool satisfied = true;
    Real epsRow = 0.0;
    Real maxRatio = 0.0;
    Index maxRatioNode = -1;
    std::string maxRatioCarrier;
    std::size_t qualifiedRowCount = 0;
    std::size_t ignoredRowCount = 0;
    std::vector<NewtonCarrierRowConvergenceViolation> violations;
};

struct NewtonCarrierRowRecoveryConfig {
    std::string mode = "off"; ///< "off" or "gummel_density".
    int maxAttempts = 1; ///< Maximum density passes per recovery cycle.
    int maxCycles = 1; ///< Maximum Newton/recovery cycles at one bias point.
    Real densityChangeReltol = 1.0e-8; ///< Stop density passes when max relative density change falls below this value.
};

struct NewtonCarrierRowRecoveryResult {
    DDSolution solution;
    bool attempted = false;
    std::string mode;
    int electronRowsUpdated = 0;
    int holeRowsUpdated = 0;
    int densityPasses = 0;
    int cyclesAttempted = 0;
    bool densityConverged = false;
    Real maxDensityRelativeChange = 0.0;
    Real maxPsiDelta_V = 0.0;
    Real maxCarrierDensityRatio = 0.0;
};

struct NewtonBlockAbsoluteConvergenceConfig {
    std::string mode = "off"; ///< "off" or "enforce".
    Real psiResidualCeiling = 0.0;
    Real electronResidualCeiling = 0.0;
    Real holeResidualCeiling = 0.0;
};

struct NewtonConfig {
    int maxIter = 20;
    Real reltol = 1.0e-8;
    Real abstol = 1.0e-18;
    Real temperature_K = constants::T0; ///< Lattice temperature [K]
    PoissonChargeVolumePolicy poissonChargeVolumePolicy =
        PoissonChargeVolumePolicy::Global;
    Real dampingFactor = 1.0;
    bool lineSearch = true;
    std::string lineSearchMode = "merit"; ///< "merit" or "block_filter".
    Real residualFilterGamma = 1.0e-4;
    Real residualFilterEnvelopeFactor = 2.0;
    bool verbose = true;
    bool warmStart = false; ///< Preserve supplied quasi-Fermi potentials instead of resetting interiors.
    bool diagnostics = false; ///< Store detailed line-search diagnostics in NewtonResult::history.
    Real maxUpdate = 0.0; ///< Optional infinity-norm cap on one Newton update in solver unknown units.
    Real quasiFermiUpdateLimit_V = 0.0; ///< Optional physical-voltage cap on phin/phip Newton updates.
    Real quasiFermiUpdateLimitMinority_V = 0.0; ///< Optional tighter physical-voltage cap on the minority-carrier quasi-Fermi update per node; 0 disables (uses the global cap for both carriers).
    std::string quasiFermiUpdateLimitMode = "componentwise_poisson_recorrect"; ///< "componentwise_poisson_recorrect" or direction-preserving "uniform_trust_region".
    Real quasiFermiTrustRegionGrowthFactor = 1.0; ///< Accepted-step radius growth; 1 disables adaptive expansion.
    Real quasiFermiTrustRegionMaxMultiplier = 1.0; ///< Maximum multiplier applied to configured quasi-Fermi limits.
    Real quasiFermiTrustRegionExpansionThreshold = 0.75; ///< Minimum actual/predicted residual decrease ratio for radius expansion.
    Real quasiFermiTrustRegionShrinkFactor = 1.0; ///< Radius contraction after a backtracked accepted step; 1 disables contraction.
    Real quasiFermiTrustRegionMinMultiplier = 1.0; ///< Positive minimum multiplier used by adaptive contraction.
    Real stallResidualFloor = 1.0e-9; ///< Residual floor for accepting line-search stalls as solved.
    Real poissonLineSearchStallResidualFloor = 1.0e-6; ///< Poisson-block floor for near-flat line-search stalls.
    Real poissonLineSearchStallRelativeIncrease = 1.0e-5; ///< Allowed best rejected residual increase at the Poisson floor.
    Real poissonLineSearchStallCarrierResidualFloor = 1.0e-6; ///< Carrier-block ceiling for Poisson-floor stall acceptance.
    Real poissonLineSearchStallContactMajorityQfDropLimit_V = 5.0e-11; ///< Maximum contact-edge majority-carrier quasi-Fermi drop allowed for Poisson-floor stall acceptance; 0 disables.
    NewtonBlockAbsoluteConvergenceConfig blockAbsoluteConvergence{}; ///< Optional authoritative per-block absolute convergence contract.
    Real contactMajorityQfBranchDropLimit_V = 0.0; ///< Maximum contact-edge majority-carrier quasi-Fermi drop allowed for every convergence path and best-iterate eligibility; 0 disables.
    std::vector<std::string> contactMajorityQfBranchGuardContacts; ///< Optional contact scope for the branch guard; empty means all transport contacts.
    bool carrierRowQualifiedStallAcceptance = false; ///< Accept a non-decreasing line-search stall within the configured block/contact floors only when enforced local carrier rows are all satisfied.
    Real carrierRegularizationScale = 0.0; ///< Optional carrier-row diagonal regularization scale.
    CarrierDiagonalFloorRegularizationConfig carrierDiagonalFloor{}; ///< Optional absolute floor for depleted minority carrier-row diagonals.
    NewtonCarrierRowConvergenceConfig carrierRowConvergence{}; ///< Optional per-carrier-row local residual convergence check.
    NewtonLocalUpdateDiagnosticsConfig localUpdateDiagnostics{}; ///< Opt-in raw/capped/applied Newton-step trace.
    NewtonCarrierRowRecoveryConfig carrierRowRecovery{}; ///< Optional recovery pass for locally unbalanced carrier rows.
    NewtonContinuityRowScalingConfig continuityRowScaling{}; ///< Optional source-aware left row equilibration.
    NewtonLinearEquilibrationConfig linearEquilibration{}; ///< Optional two-sided scaling of the linear system.
    NewtonGlobalContinuityClosureConfig globalContinuityClosure{}; ///< Optional contact-flux versus integrated-source convergence check.
    Real finiteDifferenceStep = 1.0e-6;
    std::string jacobian = "analytic"; ///< "analytic" or "finite_difference"
    std::string quasiFermiReference = "none"; ///< "none", "contact_majority", or "contact_basin"
    /// Recenter warm-start QF coordinates on the supplied physical state so
    /// sub-ULP corrections are represented as local increments.
    bool quasiFermiRecenterOnInitialState = false;
    std::string residualNorm = "block"; ///< "block" or "l2" convergence/line-search norm
    Real residualWeightPsi = 1.0;
    Real residualWeightPhin = 1.0;
    Real residualWeightPhip = 1.0;
    Real residualScalePsi = 0.0;  ///< <= 0 selects max(initial psi-block residual norm, 1)
    Real residualScalePhin = 0.0; ///< <= 0 selects max(initial electron-continuity residual norm, 1)
    Real residualScalePhip = 0.0; ///< <= 0 selects max(initial hole-continuity residual norm, 1)
    std::string contactBoundaryReconstruction = "dominant_signed_contact_mean";
    bool contactBoundaryMinorityElectronRelaxation = true;
    Real contactBoundaryMinorityElectronRelaxationBiasThreshold_V = 0.1;
    bool contactBoundaryMinorityElectronRelaxationTwoTerminalOnly = true;
    std::string contactBoundaryMinorityElectronRelaxationContactSide = "p_contact_only";
    Real contactBoundaryMinorityElectronRelaxationStrength = 1.0;
    UnitScalingConfig inputScaling{}; ///< Input-unit mode from top-level config.
    UnitScalingReferenceConfig unitScalingRefs{}; ///< Optional reference overrides.
    Real taun = 1.0e-5;
    Real taup = 3.0e-6;
    Real augerCn = 2.90e-43; ///< Electron Auger coefficient [m^6/s]
    Real augerCp = 1.028e-43; ///< Hole Auger coefficient [m^6/s]
    std::string augerExcessProduct = "generalized_fermi"; ///< "generalized_fermi" or "classical_np".
    MobilityModelConfig mobility{}; ///< Mobility model configuration
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    SRHDopingDependenceConfig srhDopingDependence{}; ///< Sentaurus SRH(DopingDep).
    BandToBandTunnelingConfig bandToBand{}; ///< Local pair generation, including Sentaurus E2.
    ImpactIonizationModelConfig impactIonization; ///< Avalanche generation model.
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
    CarrierStatisticsConfig carrierStatistics{}; ///< Boltzmann or Fermi-Dirac density/contact/transport statistics.
    DensityGradientQuantumPotentialConfig electronQuantumPotential{};
    ///< Optional Sentaurus eQuantumPotential outer iteration.
};

struct NewtonBlockResidualInfo {
    Real psi = 0.0;
    Real phin = 0.0;
    Real phip = 0.0;
    Real combined = 0.0;
};

struct NewtonResidualPeak {
    Index nodeId = 0;
    Real signedResidual = 0.0;
    Real absoluteResidual = 0.0;
};

struct NewtonIterationInfo {
    int iter = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
    Real relativeResidualNorm = 0.0;
    Real contactMajorityQfDrop = 0.0;
    Real rawStepNorm = 0.0;
    int lineSearchAttempts = 0;
    bool lineSearchAccepted = false;
    NewtonBlockResidualInfo blockResiduals;
    NewtonBlockResidualInfo rowScaledBlockResiduals;
    NewtonResidualPeak topPoissonResidual;
    NewtonResidualPeak topElectronResidual;
    NewtonResidualPeak topHoleResidual;
    std::string sourceJacobianActiveBranchFingerprint;
    std::string event;
    NewtonCarrierRowConvergenceEvaluation carrierRowConvergence;
    std::vector<LineSearchIterationInfo> lineSearchHistory;
};

struct NewtonCarrierDiagnostics {
    bool positiveFinite = true;
    Real minElectronDensity = 0.0;
    Real minHoleDensity = 0.0;
    int nonfiniteElectronCount = 0;
    int nonfiniteHoleCount = 0;
    int nonpositiveElectronCount = 0;
    int nonpositiveHoleCount = 0;
};

struct NewtonTopResidualNode {
    Index nodeId = 0;
    Real x = 0.0;
    Real y = 0.0;
    Real poissonResidual = 0.0;
    Real absPoissonResidual = 0.0;
    Real donors = 0.0;
    Real acceptors = 0.0;
    Real netDoping = 0.0;
    Real effectiveIntrinsicDensity = 0.0;
};

struct NewtonTopCarrierResidualNode {
    Index nodeId = 0;
    Real x = 0.0;
    Real y = 0.0;
    Real residual = 0.0;
    Real absResidual = 0.0;
};

struct NewtonFailureDiagnostics {
    std::string failureReason;
    int failedIteration = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
    int lineSearchAttempts = 0;
    std::string lineSearchFailureReason;
    NewtonBlockResidualInfo blockResiduals;
    NewtonCarrierDiagnostics carrierDiagnostics;
    Real maxContactMajorityQfDrop = 0.0;
    Real bestRejectedContactMajorityQfDrop = 0.0;
    std::vector<LineSearchIterationInfo> lineSearchHistory;
    std::vector<NewtonTopResidualNode> topPoissonResidualNodes;
    std::vector<NewtonTopCarrierResidualNode> topElectronResidualNodes;
    std::vector<NewtonTopCarrierResidualNode> topHoleResidualNodes;
};

struct DensityGradientOuterIterationInfo {
    int iteration = 0;
    int innerIterations = 0;
    bool innerConverged = false;
    Real innerResidualInfinityNorm = 0.0;
    Real innerLastUpdateInfinityNorm_V = 0.0;
    Index innerMaxUpdateNode = 0;
    Real innerMaxUpdateNodeValue_V = 0.0;
    Real rawChangeInfinityNorm_V = 0.0;
    Real appliedChangeInfinityNorm_V = 0.0;
    Real potentialInfinityNorm_V = 0.0;
    Real relaxation = 1.0;
};

struct NewtonResult {
    DDSolution solution;
    DDSolution bestSolution; ///< Lowest-residual accepted iterate that passes the configured branch guard.
    bool hasBestSolution = false;
    int bestIteration = 0;
    Real bestResidualNorm = std::numeric_limits<Real>::infinity();
    NewtonBlockResidualInfo bestBlockNorms;
    Real bestContactMajorityQfDrop = std::numeric_limits<Real>::infinity();
    bool converged = false;
    int iters = 0;
    Real initialResidualNorm = 0.0;
    Real finalResidualNorm = 0.0;
    std::string convergenceReason;
    NewtonBlockResidualInfo finalBlockNorms;
    NewtonCarrierRowConvergenceEvaluation finalCarrierRowConvergence;
    NewtonGlobalContinuityClosureEvaluation finalGlobalContinuityClosure;
    NewtonCarrierRowRecoveryResult carrierRowRecovery;
    std::vector<NewtonIterationInfo> history;
    std::vector<NewtonIterationInfo> trace;
    std::vector<DensityGradientOuterIterationInfo> electronQuantumOuterHistory;
    NewtonFailureDiagnostics failureDiagnostics;
};

struct NewtonResidualEvaluation {
    VectorXd raw;
    NewtonBlockResidualInfo blockNorms;
    std::vector<Real> intrinsicDensity;
    bool scaledState = false;
    Real potentialScale = 1.0;
};

struct NewtonPoissonTermEvaluation {
    std::vector<CoupledDDPoissonTermDiagnostic> rows;
    bool scaledState = false;
    Real potentialScale = 1.0;
};

struct NewtonStepEvaluation {
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonFeedbackSubstitutionEvaluation {
    std::string variant;
    bool replacesDensity = false;
    bool replacesQuasiFermi = false;
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation productionTrialResidual;
    std::vector<CoupledDDCarrierTermDiagnostic> carrierTerms;
    VectorXd desiredResidual;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    VectorXd carrierOnlyDeltaPhin;
    VectorXd carrierOnlyDeltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
    Real carrierOnlyRawStepNorm = 0.0;
    Real carrierOnlyStepNorm = 0.0;
};

struct NewtonMatrixConditionEstimate {
    Real largestSingularValue = 0.0;
    Real smallestResolvedSingularValue = 0.0;
    Real resolvedConditionNumber = 0.0;
    Index numericalRank = 0;
    Index rows = 0;
    Index columns = 0;
};

struct NewtonSchurLoopComponentEvaluation {
    std::string name;
    SparseMatrixd jacobianQfpPsi;
    SparseMatrixd effectiveLoop;
    VectorXd leaveOutDeltaPhin;
    VectorXd leaveOutDeltaPhip;
    VectorXd onlyDeltaPhin;
    VectorXd onlyDeltaPhip;
};

struct NewtonPoissonQfpCrossBlockEvaluation {
    NewtonResidualEvaluation residual;
    SparseMatrixd jacobianPsiPsi;
    SparseMatrixd jacobianPsiQfp;
    SparseMatrixd jacobianQfpPsi;
    SparseMatrixd jacobianQfpQfp;
    SparseMatrixd effectiveSchurLoop;
    std::vector<NewtonSchurLoopComponentEvaluation> loopComponents;
    VectorXd targetDeltaPhin;
    VectorXd targetDeltaPhip;
    VectorXd independentDeltaPsi;
    VectorXd independentDeltaPhin;
    VectorXd independentDeltaPhip;
    VectorXd noPsiQfpDeltaPsi;
    VectorXd noPsiQfpDeltaPhin;
    VectorXd noPsiQfpDeltaPhip;
    VectorXd noQfpPsiDeltaPsi;
    VectorXd noQfpPsiDeltaPhin;
    VectorXd noQfpPsiDeltaPhip;
    VectorXd schurDeltaPsi;
    VectorXd schurDeltaPhin;
    VectorXd schurDeltaPhip;
    VectorXd fullRawDeltaPsi;
    VectorXd fullRawDeltaPhin;
    VectorXd fullRawDeltaPhip;
    VectorXd fullCappedDeltaPsi;
    VectorXd fullCappedDeltaPhin;
    VectorXd fullCappedDeltaPhip;
    VectorXd psiQfpProduct;
    VectorXd qfpPsiProduct;
    VectorXd qfpFiniteDifferenceDirectionPhin;
    VectorXd qfpFiniteDifferenceDirectionPhip;
    VectorXd psiFiniteDifferenceDirection;
    VectorXd analyticPsiQfpDirectionalDerivative;
    VectorXd finiteDifferencePsiQfpDirectionalDerivative;
    VectorXd analyticQfpPsiDirectionalDerivative;
    VectorXd finiteDifferenceQfpPsiDirectionalDerivative;
    NewtonMatrixConditionEstimate jacobianPsiPsiCondition;
    NewtonMatrixConditionEstimate jacobianPsiPsiEquilibratedCondition;
    NewtonMatrixConditionEstimate jacobianQfpQfpCondition;
    NewtonMatrixConditionEstimate jacobianQfpQfpEquilibratedCondition;
    NewtonMatrixConditionEstimate schurCondition;
    NewtonMatrixConditionEstimate schurEquilibratedCondition;
    NewtonMatrixConditionEstimate effectiveSchurLoopCondition;
    Real jacobianPsiPsiNorm = 0.0;
    Real jacobianPsiQfpNorm = 0.0;
    Real jacobianQfpPsiNorm = 0.0;
    Real jacobianQfpQfpNorm = 0.0;
    Real fullLinearClosureNorm = 0.0;
    Real schurClosureNorm = 0.0;
    Real schurRelativeClosure = 0.0;
    Real loopComponentClosureNorm = 0.0;
    Real finiteDifferenceRelativeStep = 0.0;
    Real psiQfpDirectionalDerivativeRelativeError = 0.0;
    Real qfpPsiDirectionalDerivativeRelativeError = 0.0;
};

struct NewtonDirectionalDerivativeEvaluation {
    NewtonResidualEvaluation residual;
    VectorXd perturbationPsi;
    VectorXd perturbationPhin;
    VectorXd perturbationPhip;
    VectorXd analyticJv;
    VectorXd finiteDifferenceJv;
    VectorXd forwardResidual;
    VectorXd backwardResidual;
    Real perturbationNorm = 0.0;
    Real analyticNorm = 0.0;
    Real finiteDifferenceNorm = 0.0;
    Real absoluteError = 0.0;
    Real relativeError = 0.0;
};

struct NewtonBlockStepEvaluation {
    std::string mode;
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonPoissonBlockInitialization {
    DDSolution coldInitial;
    DDSolution poissonBlockInitial;
    NewtonBlockResidualInfo coldBlockResiduals;
    NewtonBlockResidualInfo poissonBlockResiduals;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonRegularizedCarrierStepEvaluation {
    Real regularizationScale = 0.0;
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
    Real regularizationDiagonalNorm = 0.0;
};

struct NewtonCarrierRowDiagnostic {
    Index nodeId = 0;
    Real electronResidual = 0.0;
    Real holeResidual = 0.0;
    Real electronDiagonal = 0.0;
    Real holeDiagonal = 0.0;
    Real electronRowAbsSum = 0.0;
    Real holeRowAbsSum = 0.0;
    Real electronOffdiagAbsSum = 0.0;
    Real holeOffdiagAbsSum = 0.0;
    Real electronRowL2Norm = 0.0;
    Real holeRowL2Norm = 0.0;
    Real rawDeltaPhin_V = 0.0;
    Real rawDeltaPhip_V = 0.0;
    Real cappedDeltaPhin_V = 0.0;
    Real cappedDeltaPhip_V = 0.0;
};

struct NewtonCarrierRowDiagnosticsEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<NewtonCarrierRowDiagnostic> rows;
    Real potentialScale = 1.0;
    Real rawCarrierStepNorm = 0.0;
    Real cappedCarrierStepNorm = 0.0;
};

struct NewtonPoissonLinearRowDiagnostic {
    Index nodeId = 0;
    Real residual = 0.0;
    Real diagonal = 0.0;
    Real poissonRowL2Norm = 0.0;
    Real poissonColumnL2Norm = 0.0;
    Real fullRowL2Norm = 0.0;
    Real fullColumnL2Norm = 0.0;
    Real rawDeltaPsi_V = 0.0;
    Real equilibratedDeltaPsi_V = 0.0;
};

struct NewtonPoissonLinearDiagnosticsEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<NewtonPoissonLinearRowDiagnostic> rows;
    std::vector<Index> focusPatchNodes;
    NewtonMatrixConditionEstimate focusPatchRawCondition;
    NewtonMatrixConditionEstimate focusPatchEquilibratedCondition;
    Index focusNode = 0;
    Real potentialScale = 1.0;
    Real poissonRowNormSpread = 0.0;
    Real poissonColumnNormSpread = 0.0;
    Real fullRowNormSpread = 0.0;
    Real fullColumnNormSpread = 0.0;
    Real rawStepNorm = 0.0;
    Real equilibratedStepNorm = 0.0;
    Real rawLinearClosureNorm = 0.0;
    Real equilibratedLinearClosureNorm = 0.0;
    Real rawRelativeLinearClosure = 0.0;
    Real equilibratedRelativeLinearClosure = 0.0;
    Real relativeStepDifference = 0.0;
};

struct NewtonCarrierBlockColumnDiagnostic {
    std::string carrier;
    Index nodeId = 0;
    Index reducedColumn = 0;
    Real diagonal = 0.0;
    Real columnL2Norm = 0.0;
    Real electronRowL2Norm = 0.0;
    Real holeRowL2Norm = 0.0;
    Real diagonalFraction = 0.0;
    Real crossCarrierRowFraction = 0.0;
    Real continuityRowWeight = 1.0;
    Real residual = 0.0;
    Real fullDeltaQfp_V = 0.0;
};

struct NewtonCarrierBlockSingularModeDiagnostic {
    Index modeIndex = 0;
    Real singularValue = 0.0;
    Real relativeSingularValue = 0.0;
    Real rhsProjection = 0.0;
    Real rhsEnergyFraction = 0.0;
    Real stepAmplitude = 0.0;
    Real stepEnergyFraction = 0.0;
    Real transportJacobianProjection = 0.0;
    Real recombinationJacobianProjection = 0.0;
    Real avalancheDiagonalJacobianProjection = 0.0;
    Real avalancheCrossJacobianProjection = 0.0;
    Real jacobianProjectionClosure = 0.0;
    Real transportRhsProjection = 0.0;
    Real recombinationRhsProjection = 0.0;
    Real avalancheRhsProjection = 0.0;
    Real rhsProjectionClosure = 0.0;
    Real noCrossCarrierStepAmplitude = 0.0;
    Real noRecombinationStepAmplitude = 0.0;
    Real noAvalancheStepAmplitude = 0.0;
    Real transportOnlyStepAmplitude = 0.0;
    Real rightElectronFraction = 0.0;
    Real leftElectronFraction = 0.0;
    std::string topRightCarrier;
    Index topRightNode = 0;
    Real topRightValue = 0.0;
    std::string topLeftCarrier;
    Index topLeftNode = 0;
    Real topLeftValue = 0.0;
};

struct NewtonCarrierBlockSolveVariantEvaluation {
    std::string name;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real scaledStepNorm = 0.0;
    Real physicalStepNorm_V = 0.0;
    Real relativeDifferenceFromFull = 0.0;
    Real cosineWithFull = 0.0;
    Real relativeLinearClosure = 0.0;
};

struct NewtonCarrierBlockDecompositionEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<NewtonCarrierBlockColumnDiagnostic> columns;
    std::vector<NewtonCarrierBlockSingularModeDiagnostic> singularModes;
    std::vector<NewtonCarrierBlockSolveVariantEvaluation> solveVariants;
    NewtonMatrixConditionEstimate rawCondition;
    NewtonMatrixConditionEstimate rowScaledCondition;
    NewtonMatrixConditionEstimate l2EquilibratedCondition;
    Index freeElectronUnknowns = 0;
    Index freeHoleUnknowns = 0;
    Real electronElectronNorm = 0.0;
    Real electronHoleNorm = 0.0;
    Real holeElectronNorm = 0.0;
    Real holeHoleNorm = 0.0;
    Real crossCarrierNormFraction = 0.0;
    Real recombinationCrossNorm = 0.0;
    Real avalancheCrossNorm = 0.0;
    Real transportCrossNorm = 0.0;
    Real freeColumnNormSpread = 0.0;
    Real freeRowNormSpread = 0.0;
    Real rowWeightSpread = 0.0;
};

struct NewtonCarrierTermDiagnosticsEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<CoupledDDCarrierTermDiagnostic> rows;
};

struct NewtonJacobianBlockAuditRow {
    std::string block;
    std::string configurationFingerprint;
    std::string activeBranchFingerprint;
    Real analyticNorm = 0.0;
    Real fdNorm = 0.0;
    Real diffNorm = 0.0;
    Real relDiff = 0.0;
    Real analyticPsiColumnNorm = 0.0;
    Real fdPsiColumnNorm = 0.0;
    Real diffPsiColumnNorm = 0.0;
    Real relPsiColumnDiff = 0.0;
    Real analyticPhinColumnNorm = 0.0;
    Real fdPhinColumnNorm = 0.0;
    Real diffPhinColumnNorm = 0.0;
    Real relPhinColumnDiff = 0.0;
    Real analyticPhipColumnNorm = 0.0;
    Real fdPhipColumnNorm = 0.0;
    Real diffPhipColumnNorm = 0.0;
    Real relPhipColumnDiff = 0.0;
    Real analyticElectronPhinNorm = 0.0;
    Real fdElectronPhinNorm = 0.0;
    Real diffElectronPhinNorm = 0.0;
    Real relElectronPhinDiff = 0.0;
    Real analyticElectronPhipNorm = 0.0;
    Real fdElectronPhipNorm = 0.0;
    Real diffElectronPhipNorm = 0.0;
    Real relElectronPhipDiff = 0.0;
    Real analyticHolePhinNorm = 0.0;
    Real fdHolePhinNorm = 0.0;
    Real diffHolePhinNorm = 0.0;
    Real relHolePhinDiff = 0.0;
    Real analyticHolePhipNorm = 0.0;
    Real fdHolePhipNorm = 0.0;
    Real diffHolePhipNorm = 0.0;
    Real relHolePhipDiff = 0.0;
};

class NewtonSolver {
public:
    NewtonSolver(const DeviceMesh& mesh,
                 const MaterialDatabase& matdb,
                 const DopingModel& doping,
                 const std::unordered_map<std::string, Real>& contactBiases,
                 NewtonConfig cfg = {},
                 std::vector<RegionFixedChargeSpec> fixedCharges = {},
                 std::vector<InterfaceSheetChargeSpec> sheetCharges = {},
                 ContactSpecsMap contactSpecs = {});

    NewtonResult solve() const;
    NewtonPoissonBlockInitialization buildPoissonBlockInitialization() const;
    NewtonResult solve(const DDSolution& initial) const;
    /// Solve only the nonlinear Poisson block. Interior quasi-Fermi fields are
    /// reconstructed from the nearest majority-carrier contact basin before
    /// each bias step. This is Vela's explicit approximation of the contact-
    /// potential extrapolation observed in Sentaurus ABA after equilibrium.
    NewtonResult solvePoissonOnly(const DDSolution& initial) const;
    NewtonResidualEvaluation evaluateResidual(const DDSolution& state) const;
    NewtonPoissonTermEvaluation evaluatePoissonTerms(
        const DDSolution& state) const;
    Real maxContactMajorityQuasiFermiDrop(const DDSolution& state) const;
    Real maxContactMajorityQuasiFermiDrop(
        const DDSolution& state,
        const std::vector<std::string>& contacts) const;
    NewtonStepEvaluation evaluateStep(const DDSolution& state) const;
    std::vector<NewtonFeedbackSubstitutionEvaluation>
    evaluateFeedbackSubstitutions(
        const DDSolution& state,
        const DDSolution& replacementState) const;
    NewtonPoissonQfpCrossBlockEvaluation
    evaluatePoissonQfpCrossBlockDecomposition(
        const DDSolution& state,
        const DDSolution& replacementState) const;
    NewtonDirectionalDerivativeEvaluation evaluateDirectionalDerivative(
        const DDSolution& state,
        const DDSolution& physicalPerturbation) const;
    NewtonBlockStepEvaluation evaluateBlockStep(
        const DDSolution& state,
        const std::string& mode) const;
    NewtonRegularizedCarrierStepEvaluation evaluateRegularizedCarrierStep(
        const DDSolution& state,
        Real regularizationScale) const;
    NewtonCarrierRowDiagnosticsEvaluation evaluateCarrierRowDiagnostics(
        const DDSolution& state) const;

    /// Diagnose Poisson rows and compare the production linear solve with a
    /// mathematically equivalent one-pass two-sided L2 equilibration.  This is
    /// read-only and does not alter solve().
    NewtonPoissonLinearDiagnosticsEvaluation evaluatePoissonLinearDiagnostics(
        const DDSolution& state,
        Index focusNode) const;
    NewtonCarrierBlockDecompositionEvaluation
    evaluateCarrierBlockDecomposition(const DDSolution& state) const;
    NewtonCarrierTermDiagnosticsEvaluation evaluateCarrierTermDiagnostics(
        const DDSolution& state,
        bool solvedEquationTerms = false) const;
    std::vector<NewtonJacobianBlockAuditRow> evaluateJacobianBlockAudit(
        const DDSolution& state,
        Real finiteDifferenceStep = 1.0e-7,
        std::vector<std::string> blocks = {},
        const std::string& finiteDifferenceMode = "double_symmetric",
        std::vector<Index> columnNodes = {}) const;
    std::vector<CoupledDDEdgeFluxDiagnostic> evaluateSgEdgeFluxDiagnostics(
        const DDSolution& state) const;
    std::vector<CoupledDDTransportEdgeJacobianDiagnostic>
    evaluateTransportEdgeJacobianDiagnostics(
        const DDSolution& state,
        Real physicalFiniteDifferenceStep_V = 1.0e-7) const;

    /// Build a pseudo-arclength continuation system over the coupled drift-diffusion
    /// residual, using the bias on `activeContact` (in volts) as the continuation
    /// parameter lambda. The returned callbacks operate on the scaled packed state
    /// vector [psi, phin, phip] produced by CoupledDDAssembler::pack, reusing the
    /// exact assembler, boundary-condition construction, and Jacobian assembly of
    /// the standard Newton solve. `biasFiniteDifferenceStep_V` sizes the central
    /// finite difference used to estimate dF/dV.
    ArclengthSystem makeArclengthSystem(const std::string& activeContact,
                                        Real biasFiniteDifferenceStep_V = 1.0e-4) const;

    /// Convert a physical-unit DDSolution into the scaled packed state vector and
    /// back, consistent with makeArclengthSystem. These help drive the continuation
    /// from a converged Newton solution and interpret the corrected state.
    VectorXd packArclengthState(const DDSolution& state) const;
    DDSolution unpackArclengthState(const VectorXd& x) const;
    /// Return a reusable unpacker backed by one cached assembler. This avoids
    /// rebuilding immutable mesh/material kernels for repeated terminal-current
    /// evaluations inside bordered circuit Newton.
    std::function<DDSolution(const VectorXd&)>
    makeArclengthStateUnpacker() const;
    /// Build a terminal-current functional on the same packed state and
    /// assembler as makeArclengthSystem. currentScale converts the assembler's
    /// internal current-per-depth units (excluding its line factor) to the
    /// caller's signed terminal-current units.
    ArclengthScalarFunctional makeArclengthContactCurrentFunctional(
        const std::string& contactName,
        Real currentScale) const;

private:
    void configureQuasiFermiReferences(CoupledDDAssembler& assembler) const;
    CoupledDDBoundaryConditions buildBoundaryConditions(
        const CoupledDDAssembler& assembler) const;
    CoupledDDBoundaryConditions buildBoundaryConditions(
        const CoupledDDAssembler& assembler,
        const std::unordered_map<std::string, Real>& contactBiases) const;
    std::shared_ptr<CoupledDDAssembler> makeArclengthAssembler() const;
    DDScalingSpec buildScalingSpec() const;
    DDSolution buildInitialGuess(const CoupledDDAssembler& assembler,
                                 const CoupledDDBoundaryConditions& bcs) const;
    DDSolution makeSolution(const CoupledDDAssembler& assembler,
                            const VectorXd& x,
                            int iters) const;
    NewtonResult solveClassicalWithFrozenElectronQuantumPotential(
        const DDSolution& initial,
        const VectorXd& electronQuantumPotential_V) const;

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    std::unordered_map<std::string, Real> contactBiases_;
    ContactSpecsMap contactSpecs_;
    NewtonConfig cfg_;
    std::vector<RegionFixedChargeSpec> fixedCharges_;
    std::vector<InterfaceSheetChargeSpec> sheetCharges_;
};

NewtonConfig newtonConfigFromJson(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

NewtonCarrierRowConvergenceEvaluation evaluateCarrierRowConvergence(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const NewtonCarrierRowConvergenceConfig& cfg);

NewtonGlobalContinuityClosureEvaluation evaluateGlobalContinuityClosure(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes,
    const NewtonGlobalContinuityClosureConfig& cfg);

NewtonCarrierRowRecoveryResult recoverCarrierRowsWithGummelDensity(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    const NewtonConfig& cfg,
    const DDSolution& state,
    const std::vector<NewtonCarrierRowConvergenceViolation>& violations,
    const NewtonCarrierRowRecoveryConfig& recovery,
    const ContactSpecsMap& contactSpecs = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg = {});


NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges,
                       ContactSpecsMap contactSpecs = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges,
                       ContactSpecsMap contactSpecs = {});

NewtonResult runNewtonPoissonOnly(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    const DDSolution& initial,
    const NewtonConfig& cfg,
    std::vector<RegionFixedChargeSpec> fixedCharges = {},
    std::vector<InterfaceSheetChargeSpec> sheetCharges = {},
    ContactSpecsMap contactSpecs = {});

} // namespace vela
