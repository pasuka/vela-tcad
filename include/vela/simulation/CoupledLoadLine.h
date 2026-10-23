#pragma once

#include "vela/simulation/PseudoArclength.h"
#include "vela/solver/LinearSolver.h"
#include "vela/solver/LinearRefinement.h"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace vela {

struct CoupledLoadLineConfig {
    int maxIterations = 40;
    Real equationTolerance = 1.0e-8;
    Real loadLineTolerance_V = 1.0e-6;
    /// Coefficient of the continuation voltage in the scalar boundary row:
    /// voltageCoefficient * lambda + R * I - outer = 0.  The default 1
    /// recovers a physical series resistor; 0 gives a directly augmented
    /// current-control row without changing the device equations.
    Real voltageCoefficient = 1.0;
    Real currentDirectionalStep = 1.0e-5;
    Real dampingFactor = 1.0;
    int maxLineSearchSteps = 12;
    Real maxVoltageUpdate_V = 1.0;
    bool applyDeviceUpdateLimit = true;
    std::string lineSearchMode = "merit";
    Real filterGamma = 1.0e-4;
    /// Maximum growth of each residual block relative to its solve-entry
    /// scale. A large circuit residual must never relax the device-PDE cap,
    /// and repeated accepted steps must not ratchet either envelope upward.
    Real filterEnvelopeFactor = 2.0;
    /// Optional block-inexact forcing.  Far from scalar-row closure the device
    /// block keeps equationTolerance.  As |load residual| enters
    /// inexactLoadActivationRatio * loadLineTolerance_V, the effective device
    /// tolerance grows geometrically toward inexactDeviceToleranceMax.  The
    /// same scale is used by convergence, merit, and residual-filter tests.
    bool inexactDeviceForcingEnabled = false;
    Real inexactDeviceToleranceMax = 1.0e-8;
    Real inexactLoadActivationRatio = 10.0;
    /// Linearization of the coupled device/circuit Newton step.  The direct
    /// bordered mode factors [J F_V; R I_x 1+R I_V] as one matrix and therefore
    /// remains applicable when the voltage-controlled device block J is singular.
    std::string linearSolver = "schur";
    /// Inexact-Newton block policy: once the PDE block is already within its
    /// requested tolerance, do not amplify its small residual through a
    /// near-singular device Jacobian while correcting the circuit equation.
    bool zeroConvergedDeviceResidual = false;
    /// Opt-in JVP audit for the raw bordered direction.  This evaluates an
    /// extra device Jacobian and two terminal-current residuals per Newton
    /// iteration, so production runs keep it disabled.
    bool linearizationAuditEnabled = false;
    /// Opt-in correction of the direct bordered solve using a multiprecision
    /// defect of the original unscaled system. Zero preserves the existing
    /// solve path. This is a linear accuracy control, not a nonlinear criterion.
    int directRefinementSteps = 0;
};

inline Real coupledEffectiveEquationTolerance(
    const CoupledLoadLineConfig& config,
    Real loadResidual)
{
    if (!config.inexactDeviceForcingEnabled)
        return config.equationTolerance;
    const Real loadRatio = std::abs(loadResidual) /
        config.loadLineTolerance_V;
    if (!std::isfinite(loadRatio) ||
        loadRatio >= config.inexactLoadActivationRatio) {
        return config.equationTolerance;
    }
    const Real closure = std::clamp(
        1.0 - loadRatio / config.inexactLoadActivationRatio,
        Real{0.0}, Real{1.0});
    const Real logStrict = std::log(config.equationTolerance);
    const Real logMaximum = std::log(config.inexactDeviceToleranceMax);
    return std::exp(logStrict + closure * (logMaximum - logStrict));
}

inline Real coupledLinearizedDeviceTolerance(
    const CoupledLoadLineConfig& config,
    Real loadResidual)
{
    // Tangential zero-RHS following is safe only while the PDE block is inside
    // the *current* forcing tolerance.  Using the eventual maximum cap here
    // lets a trial leave the active tolerance and then repeatedly suppresses
    // the corrective PDE right-hand side.
    return coupledEffectiveEquationTolerance(config, loadResidual);
}

inline bool coupledResidualFilterAccepts(
    const CoupledLoadLineConfig& config,
    Real equationResidual,
    Real loadResidual,
    Real trialEquationResidual,
    Real trialLoadResidual,
    Real alpha,
    Real anchoredEquationEnvelopeScale =
        std::numeric_limits<Real>::quiet_NaN(),
    Real anchoredLoadEnvelopeScale =
        std::numeric_limits<Real>::quiet_NaN())
{
    const Real equationScaled = equationResidual /
        coupledEffectiveEquationTolerance(config, loadResidual);
    const Real loadScaled =
        std::abs(loadResidual) / config.loadLineTolerance_V;
    const Real trialEquationScaled = trialEquationResidual /
        coupledEffectiveEquationTolerance(config, trialLoadResidual);
    const Real trialLoadScaled =
        std::abs(trialLoadResidual) / config.loadLineTolerance_V;
    if (!std::isfinite(trialEquationScaled) ||
        !std::isfinite(trialLoadScaled)) {
        return false;
    }

    // Independent block envelopes are essential for the very different
    // scales in a high-resistance load line.  The former shared max envelope
    // allowed an O(1 V) circuit residual to authorize an O(1) PDE residual
    // even when the device block had already converged to about 1e-6.
    const Real equationEnvelope = config.filterEnvelopeFactor *
        (std::isfinite(anchoredEquationEnvelopeScale)
            ? anchoredEquationEnvelopeScale
            : std::max(Real{1.0}, equationScaled));
    const Real loadEnvelope = config.filterEnvelopeFactor *
        (std::isfinite(anchoredLoadEnvelopeScale)
            ? anchoredLoadEnvelopeScale
            : std::max(Real{1.0}, loadScaled));
    if (trialEquationScaled > equationEnvelope ||
        trialLoadScaled > loadEnvelope) {
        return false;
    }

    const Real sufficientDecrease = config.filterGamma * alpha;
    const bool equationImproves = trialEquationScaled <=
        (1.0 - sufficientDecrease) * equationScaled;
    const bool loadImproves = trialLoadScaled <=
        (1.0 - sufficientDecrease) * loadScaled;
    return equationImproves || loadImproves;
}

struct CoupledLoadLineResult {
    bool converged = false;
    ArclengthState state;
    Real directedCurrent_A_per_um = 0.0;
    Real equationResidualNorm = std::numeric_limits<Real>::infinity();
    Real loadLineResidual_V = std::numeric_limits<Real>::infinity();
    Real effectiveEquationTolerance = 0.0;
    bool inexactConvergenceUsed = false;
    Real merit = std::numeric_limits<Real>::infinity();
    int iterations = 0;
    int lineSearchAttempts = 0;
    std::string failureReason;
    struct IterationTrace {
        struct TrialAudit {
            Real alpha = 0.0;
            bool accepted = false;
            VectorXd appliedStateUpdate;
            Real appliedVoltageUpdate_V = 0.0;
            VectorXd predictedDeviceResidual;
            VectorXd actualDeviceResidual;
            Real predictedLoadResidual_V = 0.0;
            Real actualLoadResidual_V = 0.0;
        };
        int iteration = 0;
        Real innerVoltage_V = 0.0;
        Real equationResidualNorm = std::numeric_limits<Real>::infinity();
        Eigen::Index equationResidualIndex = -1;
        Real effectiveEquationTolerance = 0.0;
        Real loadLineResidual_V = std::numeric_limits<Real>::infinity();
        Real merit = std::numeric_limits<Real>::infinity();
        Real proposedVoltageUpdate_V = 0.0;
        Real rawStateUpdateNorm = 0.0;
        Eigen::Index rawStateUpdateIndex = -1;
        Real rawVoltageUpdate_V = 0.0;
        Real limitedStateUpdateNorm = 0.0;
        Real acceptedAlpha = 0.0;
        Real bestTrialMerit = std::numeric_limits<Real>::infinity();
        Real bestTrialEquationResidualNorm =
            std::numeric_limits<Real>::infinity();
        Eigen::Index bestTrialEquationResidualIndex = -1;
        Real bestTrialLoadLineResidual_V =
            std::numeric_limits<Real>::infinity();
        bool deviceResidualZeroed = false;
        Real deviceLinearizationResidualNorm =
            std::numeric_limits<Real>::quiet_NaN();
        Real circuitLinearizationResidual =
            std::numeric_limits<Real>::quiet_NaN();
        Real analyticCurrentDirectionalDerivative =
            std::numeric_limits<Real>::quiet_NaN();
        Real finiteDifferenceCurrentDirectionalDerivative =
            std::numeric_limits<Real>::quiet_NaN();
        Real currentDirectionalRelativeError =
            std::numeric_limits<Real>::quiet_NaN();
        Real deviceJvpAnalyticNorm =
            std::numeric_limits<Real>::quiet_NaN();
        Real deviceJvpFiniteDifferenceNorm =
            std::numeric_limits<Real>::quiet_NaN();
        Real deviceJvpErrorNorm =
            std::numeric_limits<Real>::quiet_NaN();
        Real deviceJvpRelativeError =
            std::numeric_limits<Real>::quiet_NaN();
        Eigen::Index deviceJvpErrorIndex = -1;
        Real deviceJvpAnalyticAtError =
            std::numeric_limits<Real>::quiet_NaN();
        Real deviceJvpFiniteDifferenceAtError =
            std::numeric_limits<Real>::quiet_NaN();
        // Populated only by the opt-in audit. Retain actual representable
        // differences, not alpha*delta, to expose writeback rounding losses.
        VectorXd rawStateUpdate;
        VectorXd rawPredictedDeviceResidual;
        Real componentwiseBackwardError = std::numeric_limits<Real>::quiet_NaN();
        std::vector<TrialAudit> trials;
    };
    std::vector<IterationTrace> trace;
};

/// Newton solver for the fully coupled stationary series-resistor system
///
///   F(x, Vin) = 0
///   Vin + R * I(x, Vin) - Vout = 0.
///
/// The device Jacobian is reused through ArclengthSystem. The circuit row is
/// differentiated by directional finite differences of the terminal-current
/// functional, so no dense current-gradient vector is formed.
class CoupledLoadLineNewton {
public:
    using CurrentFunctional = std::function<Real(const VectorXd&, Real)>;

    CoupledLoadLineNewton(ArclengthSystem system,
                          CurrentFunctional current,
                          CoupledLoadLineConfig config,
                          ArclengthScalarFunctional currentLinearization = {})
        : system_(std::move(system))
        , current_(std::move(current))
        , config_(std::move(config))
        , currentLinearization_(std::move(currentLinearization))
    {
        if (!system_.residual || !system_.parameterDerivative ||
            !system_.solveJacobian || !current_) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: residual, parameter derivative, Jacobian "
                "solve, and current callbacks are required.");
        }
        if (config_.maxIterations <= 0 || config_.maxLineSearchSteps < 0) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: iteration limits are invalid.");
        }
        if (config_.directRefinementSteps < 0 || config_.directRefinementSteps > 8 ||
            (config_.directRefinementSteps > 0 &&
             config_.linearSolver != "direct_bordered" &&
             config_.linearSolver != "direct_bordered_qr")) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: directRefinementSteps requires an explicit "
                "direct solver and a count in [0,8].");
        }
        requirePositive(config_.equationTolerance, "equationTolerance");
        requirePositive(config_.loadLineTolerance_V, "loadLineTolerance_V");
        requirePositive(config_.currentDirectionalStep, "currentDirectionalStep");
        requirePositive(config_.dampingFactor, "dampingFactor");
        requirePositive(config_.maxVoltageUpdate_V, "maxVoltageUpdate_V");
        if (!std::isfinite(config_.voltageCoefficient)) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: voltageCoefficient must be finite.");
        }
        if (config_.dampingFactor > 1.0) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: dampingFactor must not exceed 1.");
        }
        if (config_.lineSearchMode != "merit" &&
            config_.lineSearchMode != "residual_filter") {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: lineSearchMode must be 'merit' or "
                "'residual_filter'.");
        }
        if (config_.linearSolver != "schur" &&
            config_.linearSolver != "direct_bordered" &&
            config_.linearSolver != "direct_bordered_qr" &&
            config_.linearSolver != "auto") {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: linearSolver must be 'schur', "
                "'direct_bordered', 'direct_bordered_qr', or 'auto'.");
        }
        if ((config_.linearSolver == "direct_bordered" ||
             config_.linearSolver == "direct_bordered_qr") &&
            (!system_.jacobian || !currentLinearization_.linearize)) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: direct_bordered requires direct Jacobian "
                "assembly and current-linearization callbacks.");
        }
        requirePositive(config_.filterGamma, "filterGamma");
        requirePositive(config_.filterEnvelopeFactor, "filterEnvelopeFactor");
        if (config_.filterGamma >= 1.0 || config_.filterEnvelopeFactor < 1.0) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: residual-filter controls are invalid.");
        }
        if (config_.inexactDeviceForcingEnabled &&
            (!(config_.inexactDeviceToleranceMax >= config_.equationTolerance) ||
             !std::isfinite(config_.inexactDeviceToleranceMax) ||
             !(config_.inexactLoadActivationRatio > 0.0) ||
             !std::isfinite(config_.inexactLoadActivationRatio))) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: inexact device-forcing controls are invalid.");
        }
    }

    CoupledLoadLineResult solve(const ArclengthState& initial,
                                Real outerVoltage_V,
                                Real resistance_ohm_um) const
    {
        if (initial.x.size() <= 0 || !initial.x.allFinite() ||
            !std::isfinite(initial.lambda)) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: initial state must be finite and non-empty.");
        }
        if (!(resistance_ohm_um > 0.0) || !std::isfinite(resistance_ohm_um) ||
            !std::isfinite(outerVoltage_V)) {
            throw std::invalid_argument(
                "CoupledLoadLineNewton: outer voltage and resistance are invalid.");
        }

        CoupledLoadLineResult result;
        result.state = initial;
        Real anchoredEquationEnvelopeScale =
            std::numeric_limits<Real>::quiet_NaN();
        Real anchoredLoadEnvelopeScale =
            std::numeric_limits<Real>::quiet_NaN();
        for (int iter = 0; iter <= config_.maxIterations; ++iter) {
            VectorXd f;
            Real current = std::numeric_limits<Real>::quiet_NaN();
            try {
                f = system_.residual(result.state.x, result.state.lambda);
                current = current_(result.state.x, result.state.lambda);
            } catch (const std::exception& error) {
                result.failureReason =
                    std::string("residual evaluation failed: ") + error.what();
                return result;
            }
            result.equationResidualNorm = infinityNorm(f);
            Eigen::Index maxResidualIndex = -1;
            if (f.size() > 0)
                f.cwiseAbs().maxCoeff(&maxResidualIndex);
            result.directedCurrent_A_per_um = current;
            result.loadLineResidual_V = config_.voltageCoefficient *
                result.state.lambda +
                resistance_ohm_um * current - outerVoltage_V;
            result.effectiveEquationTolerance =
                coupledEffectiveEquationTolerance(
                    config_, result.loadLineResidual_V);
            if (iter == 0) {
                anchoredEquationEnvelopeScale = std::max(
                    Real{1.0}, result.equationResidualNorm /
                        result.effectiveEquationTolerance);
                anchoredLoadEnvelopeScale = std::max(
                    Real{1.0}, std::abs(result.loadLineResidual_V) /
                        config_.loadLineTolerance_V);
            }
            result.merit = normalizedMerit(
                result.equationResidualNorm, result.loadLineResidual_V);
            result.iterations = iter;
            result.trace.push_back({
                iter,
                result.state.lambda,
                result.equationResidualNorm,
                maxResidualIndex,
                result.effectiveEquationTolerance,
                result.loadLineResidual_V,
                result.merit});
            if (!std::isfinite(result.merit)) {
                result.failureReason = "non-finite coupled residual";
                return result;
            }
            if (result.equationResidualNorm <=
                    result.effectiveEquationTolerance &&
                std::abs(result.loadLineResidual_V) <=
                    config_.loadLineTolerance_V) {
                result.inexactConvergenceUsed =
                    result.equationResidualNorm > config_.equationTolerance;
                result.converged = true;
                return result;
            }
            if (iter == config_.maxIterations)
                break;

            VectorXd fLambda;
            try {
                fLambda = system_.parameterDerivative(
                    result.state.x, result.state.lambda);
            } catch (const std::exception& error) {
                result.failureReason =
                    std::string("parameter derivative failed: ") + error.what();
                return result;
            }
            VectorXd deltaState;
            Real deltaVoltage = 0.0;
            bool updateSolved = false;
            const bool zeroDeviceResidual =
                config_.zeroConvergedDeviceResidual &&
                result.equationResidualNorm <=
                    coupledLinearizedDeviceTolerance(
                        config_, result.loadLineResidual_V);
            result.trace.back().deviceResidualZeroed = zeroDeviceResidual;
            const VectorXd linearizedDeviceResidual =
                zeroDeviceResidual
                ? VectorXd::Zero(f.size())
                : f;
            const bool canSolveDirect = system_.jacobian &&
                currentLinearization_.linearize;
            if (config_.linearSolver == "direct_bordered" ||
                config_.linearSolver == "direct_bordered_qr" ||
                (config_.linearSolver == "auto" && canSolveDirect)) {
                updateSolved = solveDirectBorderedUpdate(
                    result.state,
                    linearizedDeviceResidual,
                    fLambda, result.loadLineResidual_V,
                    resistance_ohm_um, deltaState, deltaVoltage,
                    result.failureReason);
                if (!updateSolved &&
                    (config_.linearSolver == "direct_bordered" ||
                     config_.linearSolver == "direct_bordered_qr")) {
                    if (result.failureReason.empty()) {
                        result.failureReason =
                            "direct bordered coupled Jacobian solve failed";
                    }
                    return result;
                }
            }
            if (!updateSolved) {
                updateSolved = solveSchurUpdate(
                    result.state, linearizedDeviceResidual, fLambda,
                    result.loadLineResidual_V,
                    resistance_ohm_um, deltaState, deltaVoltage,
                    result.failureReason);
            }
            if (!updateSolved)
                return result;
            if (!deltaState.allFinite() || !std::isfinite(deltaVoltage)) {
                result.failureReason = "non-finite coupled Newton update";
                return result;
            }
            SparseMatrixd auditJacobian;
            ArclengthScalarLinearization auditCurrent;
            bool hasAuditLinearization = false;
            if (config_.linearizationAuditEnabled && system_.jacobian &&
                currentLinearization_.linearize) {
                try {
                    const SparseMatrixd jacobian =
                        system_.jacobian(result.state.x, result.state.lambda);
                    const ArclengthScalarLinearization currentLinearization =
                        currentLinearization_.linearize(
                            result.state.x, result.state.lambda);
                    auditJacobian = jacobian;
                    auditCurrent = currentLinearization;
                    hasAuditLinearization = true;
                    result.trace.back().rawStateUpdate = deltaState;
                    result.trace.back().rawPredictedDeviceResidual =
                        f + jacobian * deltaState + fLambda * deltaVoltage;
                    // Evaluate the original, unscaled system in long double.
                    // This is a diagnostic, never an acceptance override.
                    std::vector<long double> sums(f.size()), denominators(f.size());
                    for (Eigen::Index i = 0; i < f.size(); ++i) {
                        const long double term = static_cast<long double>(fLambda(i)) * deltaVoltage;
                        sums[i] = static_cast<long double>(linearizedDeviceResidual(i)) + term;
                        denominators[i] = std::abs(static_cast<long double>(linearizedDeviceResidual(i))) + std::abs(term);
                    }
                    for (int k = 0; k < jacobian.outerSize(); ++k)
                        for (SparseMatrixd::InnerIterator e(jacobian, k); e; ++e) {
                            const long double term = static_cast<long double>(e.value()) * deltaState(e.col());
                            sums[e.row()] += term;
                            denominators[e.row()] += std::abs(term);
                        }
                    long double error = 0;
                    for (Eigen::Index i = 0; i < f.size(); ++i)
                        if (denominators[i] > 0)
                            error = std::max(error, std::abs(sums[i]) / denominators[i]);
                    const long double voltageTerm =
                        (static_cast<long double>(config_.voltageCoefficient) +
                         static_cast<long double>(resistance_ohm_um) * currentLinearization.parameterDerivative) * deltaVoltage;
                    long double loadSum = result.loadLineResidual_V + voltageTerm;
                    long double loadDenom = std::abs(static_cast<long double>(result.loadLineResidual_V)) + std::abs(voltageTerm);
                    for (Eigen::Index i = 0; i < deltaState.size(); ++i) {
                        const long double term = static_cast<long double>(resistance_ohm_um) *
                            currentLinearization.stateDerivative(i) * deltaState(i);
                        loadSum += term;
                        loadDenom += std::abs(term);
                    }
                    if (loadDenom > 0) error = std::max(error, std::abs(loadSum) / loadDenom);
                    result.trace.back().componentwiseBackwardError = static_cast<Real>(error);
                    const VectorXd deviceLinearizationResidual =
                        jacobian * deltaState + fLambda * deltaVoltage +
                        linearizedDeviceResidual;
                    const VectorXd analyticDeviceJvp =
                        jacobian * deltaState + fLambda * deltaVoltage;
                    const Real directionMagnitude = std::max(
                        infinityNorm(deltaState), std::abs(deltaVoltage));
                    VectorXd finiteDifferenceDeviceJvp =
                        VectorXd::Zero(analyticDeviceJvp.size());
                    if (directionMagnitude > 0.0) {
                        const Real epsilon =
                            config_.currentDirectionalStep / directionMagnitude;
                        const VectorXd fPlus = system_.residual(
                            result.state.x + epsilon * deltaState,
                            result.state.lambda + epsilon * deltaVoltage);
                        const VectorXd fMinus = system_.residual(
                            result.state.x - epsilon * deltaState,
                            result.state.lambda - epsilon * deltaVoltage);
                        finiteDifferenceDeviceJvp =
                            (fPlus - fMinus) / (2.0 * epsilon);
                    }
                    const Real analyticCurrentDirectional =
                        currentLinearization.stateDerivative.dot(deltaState) +
                        currentLinearization.parameterDerivative * deltaVoltage;
                    const Real finiteDifferenceCurrentDirectional =
                        directionalCurrentDerivative(
                            result.state, deltaState, deltaVoltage);
                    result.trace.back().deviceLinearizationResidualNorm =
                        infinityNorm(deviceLinearizationResidual);
                    result.trace.back().circuitLinearizationResidual =
                        result.loadLineResidual_V +
                        config_.voltageCoefficient * deltaVoltage +
                        resistance_ohm_um * analyticCurrentDirectional;
                    result.trace.back().analyticCurrentDirectionalDerivative =
                        analyticCurrentDirectional;
                    result.trace.back().finiteDifferenceCurrentDirectionalDerivative =
                        finiteDifferenceCurrentDirectional;
                    result.trace.back().currentDirectionalRelativeError =
                        std::abs(
                            analyticCurrentDirectional -
                            finiteDifferenceCurrentDirectional) /
                        std::max(
                            {std::abs(analyticCurrentDirectional),
                             std::abs(finiteDifferenceCurrentDirectional),
                             Real{1.0e-30}});
                    result.trace.back().deviceJvpAnalyticNorm =
                        infinityNorm(analyticDeviceJvp);
                    result.trace.back().deviceJvpFiniteDifferenceNorm =
                        infinityNorm(finiteDifferenceDeviceJvp);
                    result.trace.back().deviceJvpErrorNorm = infinityNorm(
                        analyticDeviceJvp - finiteDifferenceDeviceJvp);
                    const VectorXd deviceJvpError =
                        analyticDeviceJvp - finiteDifferenceDeviceJvp;
                    Eigen::Index deviceJvpErrorIndex = -1;
                    if (deviceJvpError.size() > 0)
                        deviceJvpError.cwiseAbs().maxCoeff(
                            &deviceJvpErrorIndex);
                    result.trace.back().deviceJvpErrorIndex =
                        deviceJvpErrorIndex;
                    if (deviceJvpErrorIndex >= 0) {
                        result.trace.back().deviceJvpAnalyticAtError =
                            analyticDeviceJvp(deviceJvpErrorIndex);
                        result.trace.back().deviceJvpFiniteDifferenceAtError =
                            finiteDifferenceDeviceJvp(deviceJvpErrorIndex);
                    }
                    result.trace.back().deviceJvpRelativeError =
                        result.trace.back().deviceJvpErrorNorm /
                        std::max(
                            {result.trace.back().deviceJvpAnalyticNorm,
                             result.trace.back().deviceJvpFiniteDifferenceNorm,
                             Real{1.0e-30}});
                } catch (const std::exception&) {
                    // The audit is diagnostic-only and must not change the
                    // nonlinear solver's accept/reject behavior.
                }
            }
            Eigen::Index rawStateUpdateIndex = -1;
            result.trace.back().rawStateUpdateNorm =
                deltaState.cwiseAbs().maxCoeff(&rawStateUpdateIndex);
            result.trace.back().rawStateUpdateIndex = rawStateUpdateIndex;
            result.trace.back().rawVoltageUpdate_V = deltaVoltage;
            if (std::abs(deltaVoltage) > config_.maxVoltageUpdate_V) {
                const Real scale =
                    config_.maxVoltageUpdate_V / std::abs(deltaVoltage);
                deltaState *= scale;
                deltaVoltage *= scale;
            }
            if (config_.applyDeviceUpdateLimit && system_.limitUpdate)
                system_.limitUpdate(result.state.x, deltaState, deltaVoltage);
            result.trace.back().limitedStateUpdateNorm = infinityNorm(deltaState);
            result.trace.back().proposedVoltageUpdate_V = deltaVoltage;

            bool accepted = false;
            Real alpha = config_.dampingFactor;
            for (int lineSearch = 0;
                 lineSearch <= config_.maxLineSearchSteps;
                 ++lineSearch) {
                ++result.lineSearchAttempts;
                ArclengthState trial;
                trial.x = result.state.x + alpha * deltaState;
                trial.lambda = result.state.lambda + alpha * deltaVoltage;
                try {
                    const VectorXd trialF =
                        system_.residual(trial.x, trial.lambda);
                    const Real trialCurrent = current_(trial.x, trial.lambda);
                    const Real trialLoadResidual = config_.voltageCoefficient *
                        trial.lambda +
                        resistance_ohm_um * trialCurrent - outerVoltage_V;
                    const Real trialMerit = normalizedMerit(
                        infinityNorm(trialF), trialLoadResidual);
                    const Real trialEquationResidual = infinityNorm(trialF);
                    if (trialMerit < result.trace.back().bestTrialMerit) {
                        Eigen::Index trialResidualIndex = -1;
                        if (trialF.size() > 0)
                            trialF.cwiseAbs().maxCoeff(&trialResidualIndex);
                        result.trace.back().bestTrialMerit = trialMerit;
                        result.trace.back().bestTrialEquationResidualNorm =
                            infinityNorm(trialF);
                        result.trace.back().bestTrialEquationResidualIndex =
                            trialResidualIndex;
                        result.trace.back().bestTrialLoadLineResidual_V =
                            trialLoadResidual;
                    }
                    const bool filterAccepted =
                        config_.lineSearchMode == "residual_filter" &&
                        residualFilterAccepts(
                            result.equationResidualNorm,
                            result.loadLineResidual_V,
                            trialEquationResidual,
                            trialLoadResidual,
                            alpha,
                            anchoredEquationEnvelopeScale,
                            anchoredLoadEnvelopeScale);
                    const bool acceptedByPolicy =
                        config_.lineSearchMode == "residual_filter"
                        ? filterAccepted
                        : trialMerit < result.merit;
                    if (hasAuditLinearization) {
                        try {
                            CoupledLoadLineResult::IterationTrace::TrialAudit audit;
                            audit.alpha = alpha;
                            audit.accepted = std::isfinite(trialMerit) && acceptedByPolicy;
                            audit.appliedStateUpdate = trial.x - result.state.x;
                            audit.appliedVoltageUpdate_V = trial.lambda - result.state.lambda;
                            audit.predictedDeviceResidual = f + auditJacobian * audit.appliedStateUpdate +
                                fLambda * audit.appliedVoltageUpdate_V;
                            audit.actualDeviceResidual = trialF;
                            audit.predictedLoadResidual_V = result.loadLineResidual_V +
                                config_.voltageCoefficient * audit.appliedVoltageUpdate_V + resistance_ohm_um *
                                (auditCurrent.stateDerivative.dot(audit.appliedStateUpdate) +
                                 auditCurrent.parameterDerivative * audit.appliedVoltageUpdate_V);
                            audit.actualLoadResidual_V = trialLoadResidual;
                            result.trace.back().trials.push_back(std::move(audit));
                        } catch (const std::exception&) {
                            // Audit allocation must not reject a valid trial.
                        }
                    }
                    if (std::isfinite(trialMerit) && acceptedByPolicy) {
                        result.state = std::move(trial);
                        result.trace.back().acceptedAlpha = alpha;
                        accepted = true;
                        break;
                    }
                } catch (const std::exception&) {
                    // Backtracking will try a smaller finite state.
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                const auto& trace = result.trace.back();
                result.failureReason =
                    "coupled line search rejected update; proposed_delta_inner_V=" +
                    diagnosticReal(trace.proposedVoltageUpdate_V) +
                    "; best_trial_equation_residual=" +
                    diagnosticReal(trace.bestTrialEquationResidualNorm) +
                    "; best_trial_equation_residual_index=" +
                    std::to_string(trace.bestTrialEquationResidualIndex) +
                    "; best_trial_load_residual_V=" +
                    diagnosticReal(trace.bestTrialLoadLineResidual_V) +
                    "; best_trial_merit=" +
                    diagnosticReal(trace.bestTrialMerit) +
                    "; raw_state_update_inf=" +
                    diagnosticReal(trace.rawStateUpdateNorm) +
                    "; raw_state_update_index=" +
                    std::to_string(trace.rawStateUpdateIndex) +
                    "; raw_delta_inner_V=" +
                    diagnosticReal(trace.rawVoltageUpdate_V) +
                    "; limited_state_update_inf=" +
                    diagnosticReal(trace.limitedStateUpdateNorm);
                return result;
            }
        }
        result.failureReason = "coupled Newton iteration budget exhausted";
        return result;
    }

    /// Extrapolate two already-converged load-line points to a requested outer
    /// voltage.  Both terminal currents are reevaluated with the same
    /// functional used by Newton, so no separately supplied outer voltage is
    /// needed for the previous point.
    ArclengthState secantPredict(const ArclengthState& previous,
                                 const ArclengthState& current,
                                 Real targetOuterVoltage_V,
                                 Real resistance_ohm_um) const
    {
        ArclengthState predicted = current;
        if (previous.x.size() != current.x.size() || current.x.size() <= 0)
            return predicted;
        const Real previousOuter = config_.voltageCoefficient * previous.lambda +
            resistance_ohm_um *
            current_(previous.x, previous.lambda);
        const Real currentOuter = config_.voltageCoefficient * current.lambda +
            resistance_ohm_um *
            current_(current.x, current.lambda);
        const Real outerDifference = currentOuter - previousOuter;
        if (!std::isfinite(previousOuter) || !std::isfinite(currentOuter) ||
            !std::isfinite(outerDifference) ||
            std::abs(outerDifference) <=
                100.0 * std::numeric_limits<Real>::epsilon()) {
            return predicted;
        }
        const Real ratio =
            (targetOuterVoltage_V - currentOuter) / outerDifference;
        if (!std::isfinite(ratio) || std::abs(ratio) > 4.0)
            return predicted;
        VectorXd stateIncrement = ratio * (current.x - previous.x);
        Real voltageIncrement =
            ratio * (current.lambda - previous.lambda);
        // Reuse the device's quasi-Fermi trust region for the predictor.  A
        // perfectly finite secant can otherwise extrapolate a low-density node
        // into a regime where both continuity derivatives underflow.  Scale
        // the complete bordered increment together so the branch direction is
        // retained rather than clipping individual state components.
        if (system_.limitUpdate)
            system_.limitUpdate(current.x, stateIncrement, voltageIncrement);
        predicted.x = current.x + stateIncrement;
        predicted.lambda = current.lambda + voltageIncrement;
        if (!predicted.x.allFinite() || !std::isfinite(predicted.lambda))
            return current;
        // Only police predictors extrapolated from an already converged
        // device state.  A secant that moves such a state far outside the
        // device residual envelope is not a useful globalization seed; fall
        // back to the last accepted point and let bordered Newton construct a
        // local tangent instead.
        try {
            const Real referenceTolerance =
                config_.inexactDeviceForcingEnabled
                ? config_.inexactDeviceToleranceMax
                : config_.equationTolerance;
            const Real currentResidual = infinityNorm(
                system_.residual(current.x, current.lambda));
            if (currentResidual <= referenceTolerance) {
                const Real predictedResidual = infinityNorm(
                    system_.residual(predicted.x, predicted.lambda));
                const Real allowedResidual = config_.filterEnvelopeFactor *
                    std::max(referenceTolerance, currentResidual);
                if (!std::isfinite(predictedResidual) ||
                    predictedResidual > allowedResidual) {
                    return current;
                }
            }
        } catch (const std::exception&) {
            return current;
        }
        return predicted;
    }

private:
    static void requirePositive(Real value, const char* name)
    {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument(
                std::string("CoupledLoadLineNewton: ") + name +
                " must be finite and positive.");
        }
    }

    static std::string diagnosticReal(Real value)
    {
        std::ostringstream stream;
        stream << std::scientific << std::setprecision(9) << value;
        return stream.str();
    }

    Real directionalCurrentDerivative(const ArclengthState& state,
                                       const VectorXd& deltaState,
                                       Real deltaVoltage) const
    {
        const Real magnitude = std::max(
            infinityNorm(deltaState), std::abs(deltaVoltage));
        if (!(magnitude > 0.0))
            return 0.0;
        const Real epsilon = config_.currentDirectionalStep / magnitude;
        const VectorXd plusState = state.x + epsilon * deltaState;
        const VectorXd minusState = state.x - epsilon * deltaState;
        const Real plus = current_(
            plusState, state.lambda + epsilon * deltaVoltage);
        const Real minus = current_(
            minusState, state.lambda - epsilon * deltaVoltage);
        return (plus - minus) / (2.0 * epsilon);
    }

    bool solveSchurUpdate(const ArclengthState& state,
                          const VectorXd& f,
                          const VectorXd& fLambda,
                          Real loadResidual,
                          Real resistance,
                          VectorXd& deltaState,
                          Real& deltaVoltage,
                          std::string& failureReason) const
    {
        VectorXd a(state.x.size());
        VectorXd z(state.x.size());
        const bool pairSolved = system_.solveJacobianPair
            ? system_.solveJacobianPair(
                  state.x, state.lambda, -f, fLambda, a, z)
            : (system_.solveJacobian(state.x, state.lambda, -f, a) &&
               system_.solveJacobian(state.x, state.lambda, fLambda, z));
        if (!pairSolved) {
            failureReason = "device Jacobian solve failed";
            return false;
        }
        const Real currentAlongA = directionalCurrentDerivative(state, a, 0.0);
        const Real currentAlongZ = directionalCurrentDerivative(state, z, 0.0);
        const VectorXd zero = VectorXd::Zero(state.x.size());
        const Real currentAlongVoltage =
            directionalCurrentDerivative(state, zero, 1.0);
        const Real denominator = config_.voltageCoefficient +
            resistance * currentAlongVoltage -
            resistance * currentAlongZ;
        if (!std::isfinite(denominator) ||
            std::abs(denominator) <=
                100.0 * std::numeric_limits<Real>::epsilon()) {
            failureReason = "singular coupled load-line denominator";
            return false;
        }
        deltaVoltage =
            (-loadResidual - resistance * currentAlongA) / denominator;
        deltaState = a - z * deltaVoltage;
        return deltaState.allFinite() && std::isfinite(deltaVoltage);
    }

    bool solveDirectBorderedUpdate(const ArclengthState& state,
                                   const VectorXd& f,
                                   const VectorXd& fLambda,
                                   Real loadResidual,
                                   Real resistance,
                                   VectorXd& deltaState,
                                   Real& deltaVoltage,
                                   std::string& failureReason) const
    {
        if (!system_.jacobian || !currentLinearization_.linearize) {
            failureReason = "direct bordered callbacks are unavailable";
            return false;
        }
        const SparseMatrixd deviceJacobian =
            system_.jacobian(state.x, state.lambda);
        const ArclengthScalarLinearization currentLinearization =
            currentLinearization_.linearize(state.x, state.lambda);
        const Eigen::Index n = state.x.size();
        if (deviceJacobian.rows() != n || deviceJacobian.cols() != n ||
            f.size() != n || fLambda.size() != n ||
            currentLinearization.stateDerivative.size() != n ||
            !currentLinearization.stateDerivative.allFinite() ||
            !std::isfinite(currentLinearization.parameterDerivative)) {
            failureReason = "direct bordered dimensions or current derivative are invalid";
            return false;
        }

        const VectorXd circuitStateDerivative =
            resistance * currentLinearization.stateDerivative;
        const Real circuitVoltageDerivative = config_.voltageCoefficient +
            resistance * currentLinearization.parameterDerivative;
        if (!circuitStateDerivative.allFinite() ||
            !std::isfinite(circuitVoltageDerivative)) {
            failureReason = "direct bordered circuit row is non-finite";
            return false;
        }

        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(static_cast<std::size_t>(
            deviceJacobian.nonZeros() + 2 * n + 1));
        for (int outer = 0; outer < deviceJacobian.outerSize(); ++outer) {
            for (SparseMatrixd::InnerIterator entry(deviceJacobian, outer);
                 entry; ++entry) {
                triplets.emplace_back(entry.row(), entry.col(), entry.value());
            }
        }
        for (Eigen::Index row = 0; row < n; ++row) {
            if (fLambda(row) != 0.0)
                triplets.emplace_back(row, n, fLambda(row));
            const Real value = circuitStateDerivative(row);
            if (value != 0.0)
                triplets.emplace_back(n, row, value);
        }
        triplets.emplace_back(n, n, circuitVoltageDerivative);

        SparseMatrixd rawAugmented(n + 1, n + 1);
        rawAugmented.setFromTriplets(triplets.begin(), triplets.end());
        rawAugmented.makeCompressed();
        VectorXd rawRhs(n + 1);
        rawRhs.head(n) = -f;
        rawRhs(n) = -loadResidual;

        // Two-sided max-norm equilibration is important for a TCAD device
        // bordered by a 1e12-ohm-um circuit row.  Row scaling alone leaves the
        // voltage column and weak contact-current columns many orders apart.
        // Include the right-hand side in each row norm.  A predictor can place
        // a carrier row where every derivative is almost underflow-small while
        // its nonlinear residual is still finite; Jacobian-only scaling would
        // then amplify the RHS toward DBL_MAX without changing the equation.
        VectorXd rowMatrixMagnitude = VectorXd::Zero(n + 1);
        for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
            for (SparseMatrixd::InnerIterator entry(rawAugmented, outer);
                 entry; ++entry) {
                if (!std::isfinite(entry.value())) {
                    failureReason = "direct bordered matrix contains a non-finite entry";
                    return false;
                }
                rowMatrixMagnitude(entry.row()) = std::max(
                    rowMatrixMagnitude(entry.row()), std::abs(entry.value()));
            }
        }
        VectorXd rowScale(n + 1);
        for (Eigen::Index row = 0; row <= n; ++row) {
            if (!(rowMatrixMagnitude(row) > 0.0) ||
                !std::isfinite(rowMatrixMagnitude(row))) {
                failureReason = "direct bordered matrix has structural zero row " +
                    std::to_string(row);
                return false;
            }
            const Real rowMagnitude = std::max(
                rowMatrixMagnitude(row), std::abs(rawRhs(row)));
            rowScale(row) = 1.0 / rowMagnitude;
        }
        VectorXd columnScale = VectorXd::Zero(n + 1);
        for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
            for (SparseMatrixd::InnerIterator entry(rawAugmented, outer);
                 entry; ++entry) {
                columnScale(entry.col()) = std::max(
                    columnScale(entry.col()),
                    std::abs(rowScale(entry.row()) * entry.value()));
            }
        }
        for (Eigen::Index column = 0; column <= n; ++column) {
            if (!(columnScale(column) > 0.0) ||
                !std::isfinite(columnScale(column))) {
                failureReason = "direct bordered matrix has structural zero column " +
                    std::to_string(column);
                return false;
            }
            columnScale(column) = 1.0 / columnScale(column);
        }

        std::vector<Eigen::Triplet<Real>> rowScaledTriplets;
        std::vector<Eigen::Triplet<Real>> scaledTriplets;
        rowScaledTriplets.reserve(
            static_cast<std::size_t>(rawAugmented.nonZeros()));
        scaledTriplets.reserve(static_cast<std::size_t>(rawAugmented.nonZeros()));
        for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
            for (SparseMatrixd::InnerIterator entry(rawAugmented, outer);
                 entry; ++entry) {
                const Real rowScaledValue =
                    rowScale(entry.row()) * entry.value();
                rowScaledTriplets.emplace_back(
                    entry.row(), entry.col(), rowScaledValue);
                scaledTriplets.emplace_back(
                    entry.row(), entry.col(),
                    rowScaledValue * columnScale(entry.col()));
            }
        }
        SparseMatrixd rowScaledAugmented(n + 1, n + 1);
        rowScaledAugmented.setFromTriplets(
            rowScaledTriplets.begin(), rowScaledTriplets.end());
        rowScaledAugmented.makeCompressed();
        SparseMatrixd augmented(n + 1, n + 1);
        augmented.setFromTriplets(
            scaledTriplets.begin(), scaledTriplets.end());
        augmented.makeCompressed();
        VectorXd rhs = VectorXd(rawRhs.array() * rowScale.array());
        VectorXd scaledUpdate;
        bool updateUsesColumnScaling = true;
        std::vector<std::string> backendDiagnostics;
        const auto recordCandidate = [&](const std::string& backend,
                                         const SparseMatrixd& matrix,
                                         const VectorXd& candidate) {
            const Real residual = candidate.size() == matrix.cols() &&
                    candidate.allFinite()
                ? infinityNorm(matrix * candidate - rhs)
                : std::numeric_limits<Real>::infinity();
            const Real updateNorm = candidate.size() > 0 && candidate.allFinite()
                ? infinityNorm(candidate)
                : std::numeric_limits<Real>::infinity();
            std::ostringstream detail;
            detail << backend << "_residual=" << std::scientific
                   << std::setprecision(9) << residual
                   << "," << backend << "_update_inf=" << updateNorm;
            backendDiagnostics.push_back(detail.str());
            return residual;
        };
        if (config_.linearSolver == "direct_bordered_qr") {
            if (solveSpqrSystem(rowScaledAugmented, rhs, scaledUpdate)) {
                const Real spqrResidual = recordCandidate(
                    "spqr", rowScaledAugmented, scaledUpdate);
                const Real spqrTolerance = 1.0e-8 *
                    std::max(Real{1.0}, infinityNorm(rhs));
                if (!std::isfinite(spqrResidual) ||
                    spqrResidual > spqrTolerance) {
                    scaledUpdate.resize(0);
                } else {
                    updateUsesColumnScaling = false;
                }
            } else {
                backendDiagnostics.emplace_back("spqr_failed");
            }
        }
        Eigen::SparseLU<SparseMatrixd> lu;
        if (scaledUpdate.size() != n + 1)
            lu.compute(augmented);
        if (scaledUpdate.size() != n + 1 && lu.info() == Eigen::Success) {
            scaledUpdate = lu.solve(rhs);
            recordCandidate("eigen_lu", augmented, scaledUpdate);
        } else if (scaledUpdate.size() != n + 1) {
            if (solveUmfPackSystem(augmented, rhs, scaledUpdate)) {
                const Real umfpackResidual = recordCandidate(
                    "umfpack", augmented, scaledUpdate);
                const Real umfpackTolerance = 1.0e-8 *
                    std::max(Real{1.0}, infinityNorm(rhs));
                if (!std::isfinite(umfpackResidual) ||
                    umfpackResidual > umfpackTolerance) {
                    scaledUpdate.resize(0);
                }
            } else {
                backendDiagnostics.emplace_back("umfpack_failed");
                scaledUpdate.resize(0);
            }
            if (scaledUpdate.size() != n + 1 &&
                solveSpqrSystem(rowScaledAugmented, rhs, scaledUpdate)) {
                const Real spqrResidual = recordCandidate(
                    "spqr", rowScaledAugmented, scaledUpdate);
                const Real spqrTolerance = 1.0e-8 *
                    std::max(Real{1.0}, infinityNorm(rhs));
                if (!std::isfinite(spqrResidual) ||
                    spqrResidual > spqrTolerance) {
                    scaledUpdate.resize(0);
                } else {
                    updateUsesColumnScaling = false;
                }
            } else if (scaledUpdate.size() != n + 1) {
                backendDiagnostics.emplace_back("spqr_failed");
            }
            if (scaledUpdate.size() != n + 1) {
                Eigen::SparseLU<SparseMatrixd> relaxedPivotLu;
                relaxedPivotLu.setPivotThreshold(1.0e-2);
                relaxedPivotLu.compute(augmented);
                if (relaxedPivotLu.info() == Eigen::Success) {
                    scaledUpdate = relaxedPivotLu.solve(rhs);
                    recordCandidate("relaxed_eigen_lu", augmented, scaledUpdate);
                    if (relaxedPivotLu.info() != Eigen::Success)
                        scaledUpdate.resize(0);
                } else {
                    backendDiagnostics.emplace_back("relaxed_eigen_lu_failed");
                }
            }
            if (scaledUpdate.size() != n + 1) {
                // Static-pivot regularization is applied only to obtain a
                // factorization.  The update is subsequently audited against
                // the original, unregularized bordered matrix.
                for (const Real shift : {1.0e-14, 1.0e-12, 1.0e-10}) {
                    SparseMatrixd regularized = augmented;
                    for (Eigen::Index diagonal = 0; diagonal <= n; ++diagonal)
                        regularized.coeffRef(diagonal, diagonal) += shift;
                    regularized.makeCompressed();
                    Eigen::SparseLU<SparseMatrixd> regularizedLu;
                    regularizedLu.compute(regularized);
                    if (regularizedLu.info() != Eigen::Success)
                        continue;
                    VectorXd candidate = regularizedLu.solve(rhs);
                    if (regularizedLu.info() != Eigen::Success ||
                        !candidate.allFinite()) {
                        continue;
                    }
                    std::ostringstream shiftName;
                    shiftName << "shift_" << std::scientific
                              << std::setprecision(1) << shift;
                    const Real candidateResidual = recordCandidate(
                        shiftName.str(), augmented, candidate);
                    const Real candidateTolerance = 1.0e-8 *
                        std::max(Real{1.0}, infinityNorm(rhs));
                    if (std::isfinite(candidateResidual) &&
                        candidateResidual <= candidateTolerance) {
                        scaledUpdate = std::move(candidate);
                        break;
                    }
                }
            }
            if (scaledUpdate.size() != n + 1) {
                Eigen::Index currentDerivativeNonzeros = 0;
                for (Eigen::Index i = 0; i < n; ++i) {
                    if (currentLinearization.stateDerivative(i) != 0.0)
                        ++currentDerivativeNonzeros;
                }
                std::ostringstream detail;
                detail << "direct bordered factorization failed; n=" << n + 1
                       << "; nnz=" << augmented.nonZeros()
                       << "; rhs_inf=" << std::scientific
                       << std::setprecision(9) << infinityNorm(rhs)
                       << "; current_derivative_nnz="
                       << currentDerivativeNonzeros
                       << "; current_derivative_inf="
                       << infinityNorm(currentLinearization.stateDerivative);
                for (const std::string& backend : backendDiagnostics)
                    detail << "; " << backend;
                failureReason = detail.str();
                return false;
            }
        }
        if (scaledUpdate.size() != n + 1 || !scaledUpdate.allFinite()) {
            failureReason =
                "direct bordered linear solve returned a non-finite update";
            return false;
        }
        const Real linearResidual =
            updateUsesColumnScaling
            ? infinityNorm(augmented * scaledUpdate - rhs)
            : infinityNorm(rowScaledAugmented * scaledUpdate - rhs);
        const Real linearTolerance = 1.0e-8 *
            std::max(Real{1.0}, infinityNorm(rhs));
        if (!std::isfinite(linearResidual) ||
            linearResidual > linearTolerance) {
            failureReason = "direct bordered linear solve residual too large: " +
                std::to_string(linearResidual);
            return false;
        }
        VectorXd update = updateUsesColumnScaling
            ? VectorXd(columnScale.array() * scaledUpdate.array())
            : scaledUpdate;
        if (config_.directRefinementSteps > 0) {
            // Always measure the original system, including a candidate obtained
            // from a regularized factorization or the row-scaled QR path.
            auto defect = linear_refinement::highPrecisionDefect(
                rawAugmented, -rawRhs, update);
            constexpr Real backwardTolerance = 1.0e-10;
            if (defect.backwardMax > backwardTolerance) {
                Eigen::SparseLU<SparseMatrixd> correctionLu;
                correctionLu.compute(augmented);
                if (correctionLu.info() != Eigen::Success) {
                    failureReason = "direct bordered refinement factorization failed";
                    return false;
                }
                for (int k = 0; k < config_.directRefinementSteps &&
                     defect.backwardMax > backwardTolerance; ++k) {
                    const VectorXd correction = correctionLu.solve(
                        (-defect.values.cwiseProduct(rowScale)).eval());
                    if (correctionLu.info() != Eigen::Success || !correction.allFinite()) {
                        failureReason = "direct bordered refinement correction failed";
                        return false;
                    }
                    VectorXd candidate = update + correction.cwiseProduct(columnScale);
                    if (!candidate.allFinite()) {
                        failureReason = "direct bordered refinement update is non-finite";
                        return false;
                    }
                    auto next = linear_refinement::highPrecisionDefect(
                        rawAugmented, -rawRhs, candidate);
                    if (next.backwardMax >= defect.backwardMax)
                        break;
                    update = std::move(candidate);
                    defect = std::move(next);
                }
            }
            if (defect.backwardMax > backwardTolerance) {
                failureReason = "direct bordered refinement did not meet componentwise accuracy";
                return false;
            }
            const Real refinedScaledResidual =
                infinityNorm(defect.values.cwiseProduct(rowScale));
            if (!std::isfinite(refinedScaledResidual) ||
                refinedScaledResidual > linearTolerance) {
                failureReason = "direct bordered refined residual exceeds the scaled tolerance";
                return false;
            }
        }
        deltaState = update.head(n);
        deltaVoltage = update(n);
        return true;
    }

    Real normalizedMerit(Real equationResidual, Real loadResidual) const
    {
        return std::max(
            equationResidual /
                coupledEffectiveEquationTolerance(config_, loadResidual),
            std::abs(loadResidual) / config_.loadLineTolerance_V);
    }

    bool residualFilterAccepts(Real equationResidual,
                               Real loadResidual,
                               Real trialEquationResidual,
                               Real trialLoadResidual,
                               Real alpha,
                               Real anchoredEquationEnvelopeScale,
                               Real anchoredLoadEnvelopeScale) const
    {
        return coupledResidualFilterAccepts(
            config_, equationResidual, loadResidual,
            trialEquationResidual, trialLoadResidual, alpha,
            anchoredEquationEnvelopeScale, anchoredLoadEnvelopeScale);
    }

    static Real infinityNorm(const VectorXd& values)
    {
        Real norm = 0.0;
        for (Eigen::Index i = 0; i < values.size(); ++i) {
            const Real magnitude = std::abs(values(i));
            if (!std::isfinite(magnitude))
                return std::numeric_limits<Real>::infinity();
            norm = std::max(norm, magnitude);
        }
        return norm;
    }

    ArclengthSystem system_;
    CurrentFunctional current_;
    CoupledLoadLineConfig config_;
    ArclengthScalarFunctional currentLinearization_;
};

} // namespace vela
