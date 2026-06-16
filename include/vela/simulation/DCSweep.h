#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/core/Types.h"
#include "vela/simulation/CurveSweep.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/post/ContactCurrent.h"
#include "vela/post/TerminalCharge.h"
#include "vela/post/StoredCharge.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"
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

struct TransportDiagnosticsConfig {
    bool enabled = false;
};

struct ContinuityBalanceDiagnosticsConfig {
    bool enabled = false;
    std::vector<std::string> contacts;
    std::string csvFile;
};

struct SweepDiagnosticsConfig {
    TerminalBalanceDiagnosticsConfig terminalBalance;
    ContactEdgeDiagnosticsConfig contactEdge;
    TransportDiagnosticsConfig transport;
    ContinuityBalanceDiagnosticsConfig continuityBalance;
};

struct DCSweepConfig {
    CurveSweepMode mode = CurveSweepMode::IV;
    std::string contact;
    Real start = 0.0;
    Real stop = 0.0;
    Real step = 0.0;
    std::vector<Real> biasPoints;
    Real minStep = 0.0;
    Real maxStep = 0.0;
    Real growthFactor = 1.0;
    Real shrinkFactor = 0.5;
    int maxRetries = 5;
    bool stopOnFailure = true;
    std::string currentContact;
    bool writeVtk = false;
    std::string vtkPrefix;
    std::string csvFile = "dc_sweep.csv";
    std::string initialStateFile;
    std::string writeStateFile;
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
    UnitScalingConfig scaling;
};

struct DCSweepPoint {
    Real voltage = 0.0;
    Real bias = 0.0;
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
    std::string outputCsv;
    std::string outputVtk;
};

struct DCSweepResult {
    DeviceMesh mesh;
    std::vector<DCSweepPoint> points;
};

class DCSweep {
public:
    std::vector<DCSweepPoint> run(const std::string& configFile) const;
    DCSweepResult runWithResult(const std::string& configFile) const;
};

} // namespace vela
