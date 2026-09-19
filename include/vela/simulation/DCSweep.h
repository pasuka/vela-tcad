#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/core/Types.h"
#include "vela/simulation/CurveSweep.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/post/ContactCurrent.h"
#include "vela/post/TerminalCharge.h"
#include "vela/post/StoredCharge.h"
#include "vela/simulation/PseudoArclength.h"
#include "vela/simulation/QfBoundsGuard.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace vela {

struct BVReverseCriteria {
    Real maxElectricField_V_per_m = 0.0;
    Real currentJumpRatio = 0.0;
    bool nonConvergenceBreakdown = true;
};

struct ContactEdgeDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
    std::string csvFile;
};

struct TerminalBalanceDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
    std::string csvFile;
};

struct SrhBalanceDiagnosticsConfig {
    bool enabled = false;
    std::string material = "Si";
    std::string drainContact = "drain";
    std::string substrateContact = "substrate";
    std::vector<std::string> kclContacts;
    std::string csvFile;
    /// A point is resolved only when |Id|/|KCL residual| reaches this margin.
    Real resolutionMarginRatio = 10.0;
};

struct TransportDiagnosticsConfig {
    bool enabled = false;
};

struct ContinuityBalanceDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
    std::string csvFile;
};

struct SgAvalancheEdgeDiagnosticsConfig {
    bool enabled = false;
    std::string csvFile;
};

struct PathIonizationDiagnosticsConfig {
    bool enabled = false;
    std::string csvFile;
    /// Optional ordered, one-row-per-segment path trace.  This preserves the
    /// geometry and local field/alpha support needed for external path audits.
    std::string segmentsCsvFile;
    std::size_t maxPaths = 1;
    std::size_t breakRank = 0;
    Real breakValue = 1.0;
    std::string drivingForce = "solver";
    /// Stop tracing once the electrostatic field falls below this SI value.
    /// Zero retains the complete monotone path to the mesh boundary.
    Real stopField_V_per_m = 0.0;
    /// Optional carrier-specific support cutoffs on the already traced shared
    /// geometric path. Zero inherits the full path selected by stopField.
    Real electronStopField_V_per_m = 0.0;
    Real holeStopField_V_per_m = 0.0;
    /// carrier_integral_arithmetic or carrier_alpha_length_arithmetic.
    std::string meanDefinition = "carrier_integral_arithmetic";
    /// path_mean (validated default) or carrier_integrals (audit alternative).
    std::string breakOrdering = "path_mean";
    std::string tracingMode = "edge_graph"; ///< edge_graph or continuous_cell.
    /// Retain one path per local maximum by default; narrower compatibility
    /// modes can merge nearby maxima.
    std::string pathRetention = "distinct_local_maxima";
    /// Reconstructed nodal maxima are the validated P1 path seeds; element
    /// maxima remain available as an audit alternative.
    std::string seedMode = "nodal_local_maxima";
    /// Path tangent: reconstructed field, true SG current, or an explicit
    /// quasi-Fermi-gradient diagnostic direction.
    std::string tracingVector = "electric_field";
    /// Carrier-current directions below this fraction of the global nodal
    /// current maximum are numerically unobservable and fall back to E.
    Real tracingCurrentRelativeFloor = 1.0e-8;
    /// QF gradients below this fraction of the corresponding global nodal
    /// maximum are directionally unobservable and fall back to electric field
    /// in sentaurus_eparallel_adaptive mode.
    Real tracingQfRelativeFloor = 5.1e-3;
    std::string tracingDirection = "bidirectional"; ///< bidirectional/along_vector/opposite_vector.
    Real seedField_V_per_m = 0.0; ///< Optional local-maximum qualification.
};

struct TriangleGssSourceDiagnosticsConfig {
    bool enabled = false;
    std::string csvFile;
};

struct BVProcessProbeDiagnosticsConfig {
    bool enabled = false;
    std::string csvFile;
};

struct AvalancheInternalSourceCurrentAuditConfig {
    bool enabled = false;
    std::string csvFile;
    std::string summaryFile;
};

struct ReleaseBVConfigAuditConfig {
    bool enabled = false;
    std::string csvFile;
    std::string summaryFile;
    Real diagnosticReferenceAScale = 2.0;
    Real diagnosticReferenceBScale = 1.05;
    std::string diagnosticReferenceSourceMappingMode = "edge_F_edge_alpha_edge_G_to_node";
    Real diagnosticReferenceQGFull_A_per_um = 0.0;
    Real diagnosticReferenceQGJunction_A_per_um = 0.0;
};

struct TerminalCurrentMethodCompareDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
    std::string csvFile;
};

struct NewtonHistoryDiagnosticsConfig {
    bool enabled = false;
    std::string csvFile;
    std::string attemptsCsvFile;
    std::string iterationsCsvFile;
    /// Optional directory for parent, initial, and final states of rejected
    /// nonlinear attempts. Empty keeps the diagnostic disabled.
    std::string rejectedStateDirectory;
};

struct ContactCurrentQfFloorDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
};

struct SweepDiagnosticsConfig {
    TerminalBalanceDiagnosticsConfig terminalBalance;
    SrhBalanceDiagnosticsConfig srhBalance;
    ContactEdgeDiagnosticsConfig contactEdge;
    TransportDiagnosticsConfig transport;
    ContinuityBalanceDiagnosticsConfig continuityBalance;
    SgAvalancheEdgeDiagnosticsConfig sgAvalancheEdges;
    PathIonizationDiagnosticsConfig pathIonizationIntegrals;
    TriangleGssSourceDiagnosticsConfig triangleGssSources;
    BVProcessProbeDiagnosticsConfig bvProcessProbe;
    AvalancheInternalSourceCurrentAuditConfig avalancheInternalSourceCurrentAudit;
    ReleaseBVConfigAuditConfig releaseBVConfigAudit;
    TerminalCurrentMethodCompareDiagnosticsConfig terminalCurrentMethodCompare;
    NewtonHistoryDiagnosticsConfig newtonHistory;
    QfBoundsDiagnosticsConfig qfBounds;
    ContactCurrentQfFloorDiagnosticsConfig contactCurrentQfFloor;
};

struct SweepPredictorConfig {
    std::string mode = "none";
    std::vector<std::string> fields;
    Real maxExtrapolationRatio = 2.0;
};

struct SweepBranchAcceptanceConfig {
    bool terminalCurrentConsistency = false;
    Real minTerminalCurrentRatio = 0.0;
    bool terminalKcl = false;
    std::vector<std::string> terminalKclContacts;
    Real minCurrentToKclRatio = 0.0;
    bool psiPhinJump = false;
    Real maxPsiPhinJump_V = 0.0;
    bool carrierDensityJump = false;
    Real maxElectronDensityJumpDex = 0.0;
    Real maxElectronDensityJumpP95AbsDex =
        std::numeric_limits<Real>::infinity();
};

struct SweepArclengthConfig {
    bool enabled = false;
    /// Arclength predictor type. Only "tangent" is currently supported.
    std::string predictor = "tangent";
    /// Numerical parameters forwarded to PseudoArclengthContinuation. The bias
    /// voltage acts as the continuation parameter lambda.
    PseudoArclengthConfig core;
    /// Finite-difference step (in volts) used to estimate dF/dV at the active contact.
    Real biasFiniteDifferenceStep_V = 1.0e-4;
    /// Optional impact-ionization source Jacobian used only by the arclength
    /// tangent and corrector. Empty preserves the main Newton configuration.
    std::string sourceJacobianMode;
    /// Optional earlier converged state used with the sweep start state to form
    /// the first branch tangent by a secant.  This avoids relying on a single,
    /// potentially ill-scaled Jacobian tangent at the restart point.
    std::string initialSecantStateFile;
    Real initialSecantBias_V = std::numeric_limits<Real>::quiet_NaN();
};

struct SweepContinuationConfig {
    SweepPredictorConfig predictor;
    SweepBranchAcceptanceConfig branchAcceptance;
    SweepArclengthConfig arclength;
};

struct SweepInitializationConfig {
    std::string mode = "none";
    std::string diagnosticCsv;
    std::string writeStateFile;
};

struct ExternalResistorControlConfig {
    bool enabled = false;
    std::string solver = "nested_scalar";
    Real resistance_ohm_um = 0.0;
    Real currentDirection = 1.0;
    Real initialInnerVoltage_V = 0.0;
    Real residualTolerance_V = 1.0e-6;
    Real coupledVoltageCoefficient = 1.0;
    Real voltageTolerance_V = 1.0e-8;
    Real maxInnerVoltageStep_V = 0.1;
    int maxBracketSteps = 200;
    int maxIterations = 40;
    Real coupledEquationTolerance = 1.0e-8;
    Real currentDirectionalStep = 1.0e-5;
    Real coupledDampingFactor = 1.0;
    int coupledMaxLineSearchSteps = 12;
    Real coupledInitialOuterStep_V = 0.0;
    Real coupledMinOuterStep_V = 1.0e-6;
    Real coupledMaxOuterStep_V = 0.0;
    Real coupledOuterGrowthFactor = 1.5;
    Real coupledOuterShrinkFactor = 0.5;
    int coupledMaxStepRetries = 12;
    bool coupledApplyDeviceUpdateLimit = true;
    std::string coupledLineSearchMode = "merit";
    std::string coupledLinearSolver = "direct_bordered";
    Real coupledFilterGamma = 1.0e-4;
    Real coupledFilterEnvelopeFactor = 2.0;
    bool coupledInexactDeviceForcingEnabled = false;
    Real coupledInexactDeviceToleranceMax = 1.0e-8;
    Real coupledInexactLoadActivationRatio = 10.0;
    bool coupledInexactZeroConvergedResidual = false;
    bool coupledLinearizationAudit = false;
};

struct VoltageToCurrentControlConfig {
    bool enabled = false;
    Real switchVoltage_V = 0.0;
    Real currentDirection = 1.0;
    std::vector<Real> currentPoints_A_per_um;
    Real currentTolerance_A_per_um = 1.0e-10;
    Real voltageTolerance_V = 1.0e-8;
    Real maxInnerVoltageStep_V = 0.05;
    int maxBracketSteps = 200;
    int maxIterations = 40;
};

struct BoundaryControlPersistenceConfig {
    std::string evaluationCsv;
    std::string checkpointDirectory;
    bool resume = false;
    Real predictorMaxStepFactor = 4.0;
    int preferredMaxEvaluations = 3;
    bool adaptiveDeviceContinuation = false;
};

struct DCSweepConfig {
    CurveSweepMode mode = CurveSweepMode::IV;
    std::string contact;
    Real start = 0.0;
    Real stop = 0.0;
    Real step = 0.0;
    Real initialStep = 0.0;
    std::vector<Real> biasPoints;
    Real minStep = 0.0;
    Real maxStep = 0.0;
    Real growthFactor = 1.0;
    /// fixed (legacy) or newton_iterations for ordinary Newton voltage sweeps.
    std::string stepGrowthMode = "fixed";
    Real shrinkFactor = 0.5;
    int maxRetries = 5;
    bool stopOnFailure = true;
    std::string currentContact;
    bool writeVtk = false;
    std::string vtkPrefix;
    std::string csvFile = "dc_sweep.csv";
    std::string initialStateFile;
    std::string writeStateFile;
    bool frozenStateComputeCurrent = false;
    SweepInitializationConfig initialization;
    ExternalResistorControlConfig externalResistor;
    VoltageToCurrentControlConfig voltageToCurrent;
    BoundaryControlPersistenceConfig boundaryControl;
    std::string writeStateEveryPointPrefix;
    /// Optional diagnostic history for accepted adaptive steps between
    /// requested bias points. writeStateFile remains the rolling checkpoint.
    std::string writeStateEveryAcceptedStepPrefix;
    std::string chargeContact;
    std::vector<std::string> chargeRegions;
    Real chargeContactRadius = 0.0;
    bool chargePerMeter = true;
    Real chargeDepth_m = 1.0;
    std::vector<TerminalChargeConfig> terminalCharges;
    bool storedChargeEnabled = false;
    StoredChargeConfig storedCharge;
    BVReverseCriteria breakdown;
    SweepDiagnosticsConfig diagnostics;
    SweepContinuationConfig continuation;
    UnitScalingConfig scaling;
};

struct DCSweepPoint {
    Real voltage = 0.0;
    Real bias = 0.0;
    Real innerVoltage_V = 0.0;
    Real outerVoltage_V = 0.0;
    Real seriesResistance_ohm_um = 0.0;
    Real loadLineResidual_V = 0.0;
    Real targetCurrent_A_per_um = 0.0;
    Real currentBoundaryResidual_A_per_um = 0.0;
    int boundaryControlEvaluations = 0;
    std::string boundaryControlMode;
    Real electronCurrent = 0.0;
    Real electronDriftCurrent = 0.0;
    Real electronDiffusionCurrent = 0.0;
    Real holeCurrent = 0.0;
    Real holeDriftCurrent = 0.0;
    Real holeDiffusionCurrent = 0.0;
    Real totalCurrent = 0.0;
    bool converged = false;
    int iterations = 0;
    std::string solverMethod;
    int gummelIterations = 0;
    int newtonIterations = 0;
    std::string handoffStage;
    std::string newtonConvergenceReason;
    Real finalPsiResidualNorm = 0.0;
    Real finalElectronContinuityResidualNorm = 0.0;
    Real finalHoleContinuityResidualNorm = 0.0;
    int carrierRowViolations = 0;
    Real carrierRowMaxRatio = 0.0;
    bool globalContinuityClosureSatisfied = true;
    Real globalElectronContinuityClosureRatio = 0.0;
    Real globalHoleContinuityClosureRatio = 0.0;
    Real globalElectronContactFlux = 0.0;
    Real globalHoleContactFlux = 0.0;
    Real globalElectronIntegratedSource = 0.0;
    Real globalHoleIntegratedSource = 0.0;
    bool carrierRowRecoveryAttempted = false;
    int carrierRowRecoveryElectronRows = 0;
    int carrierRowRecoveryHoleRows = 0;
    int carrierRowRecoveryDensityPasses = 0;
    int carrierRowRecoveryCycles = 0;
    Real carrierRowRecoveryMaxDensityRelativeChange = 0.0;
    Real carrierRowRecoveryMaxPsiDelta_V = 0.0;
    Real carrierRowRecoveryMaxDensityRatio = 0.0;
    Real attemptedStep = 0.0;
    Real acceptedStep = 0.0;
    int retryCount = 0;
    Real terminalCharge = 0.0;
    Real capacitance = 0.0;
    std::vector<std::pair<std::string, Real>> terminalChargeValues;
    std::vector<std::pair<std::string, Real>> terminalCapacitanceValues;
    std::vector<std::pair<std::string, Real>> extraFields;
    Real maxElectricField = 0.0;
    Real currentJumpRatio = 0.0;
    bool breakdownDetected = false;
    Real breakdownVoltage = 0.0;
    std::string breakdownCriterion;
    bool failed = false;
    Real lastStableBias = 0.0;
    Real failedBias = 0.0;
    std::string failureReason;
    std::string newtonFailureClass;
    std::string failureDiagnosticsJson;
    NewtonFailureDiagnostics newtonFailureDiagnostics;
    std::string validationDiagnostics;
    int qfBoundsViolations = 0;
    bool qfBoundsRecovered = false;
    std::string predictorMode;
    bool predictedInitialState = false;
    std::string branchAcceptanceStatus;
    std::string branchAcceptanceReason;
    Real terminalCurrentConsistencyRatio = 1.0;
    Real terminalKclResidual = 0.0;
    Real currentToTerminalKclRatio = 0.0;
    Real psiPhinMaxJump_V = 0.0;
    Real electronDensityJumpMedianDex = 0.0;
    Real electronDensityJumpP95AbsDex = 0.0;
    Real electronDensityJumpMaxAbsDex = 0.0;
    Index electronDensityJumpMaxNode = -1;
    std::string outputCsv;
    std::string outputVtk;
};

struct TerminalKclAcceptanceEvaluation {
    Real residual = 0.0;
    Real currentMagnitude = 0.0;
    Real currentToResidualRatio = 0.0;
    bool satisfied = true;
};

namespace detail {

TerminalKclAcceptanceEvaluation evaluateTerminalKclAcceptance(
    Real monitoredCurrent,
    const std::vector<Real>& terminalCurrents,
    Real minCurrentToResidualRatio);

} // namespace detail

struct ReleaseBVConfigAuditMetadata {
    bool enabled = false;
    std::string model;
    std::string couplingMode;
    std::string drivingForce;
    std::string parameterSet;
    Real aScale = 1.0;
    Real bScale = 1.0;
    Real switchField_V_per_cm = 0.0;
    Real minimumField_V_per_cm = 0.0;
    std::string smoothing;
    Real electronRefDens_cm3 = 0.0;
    Real holeRefDens_cm3 = 0.0;
    std::string sourceMappingMode;
    std::string currentMagnitudeMode;
    std::string lambdaAva;
    Real depth2D_um = 1.0;
    std::string currentNormalization;
    std::string qGNormalization;
    std::string auditCsvFile;
    std::string auditSummaryFile;
};

struct DCSweepResult {
    DeviceMesh mesh;
    std::vector<DCSweepPoint> points;
    std::optional<ReleaseBVConfigAuditMetadata> releaseBVConfigAudit;
};

class DCSweep {
public:
    /// Opt-in cache for sequential requests. Each solve receives fresh copies;
    /// no accepted/rejected nonlinear state is carried between requests.
    explicit DCSweep(bool reusePreparedInputs = false, bool reuseLinearAnalysis = false)
        : reusePreparedInputs_(reusePreparedInputs || reuseLinearAnalysis),
          reuseLinearAnalysis_(reuseLinearAnalysis) {}
    /// Clear on rejected/invalid worker requests; never retain rejected states.
    void clearLinearContext() const { sequentialLinearSolver_.reset(); }
    std::vector<DCSweepPoint> run(const std::string& configFile) const;
    DCSweepResult runWithResult(const std::string& configFile) const;
private:
    struct PreparedInputs;
    bool reusePreparedInputs_;
    bool reuseLinearAnalysis_;
    mutable std::shared_ptr<LinearSolver> sequentialLinearSolver_;
    mutable std::shared_ptr<PreparedInputs> preparedInputs_;
};

} // namespace vela
