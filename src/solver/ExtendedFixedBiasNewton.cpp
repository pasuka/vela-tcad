#include "vela/solver/ExtendedFixedBiasNewton.h"

#include "vela/simulation/BorderedLoadLineSystem.h"
#include "vela/simulation/ExtendedBorderedLoadLine.h"
#include "vela/simulation/LoadLineResidual.h"
#include <cmath>
#include <stdexcept>

namespace vela {
namespace {

BasicDensityUpdateError<long double> physicalUpdateError(
    const ExtendedPopulationState& state, const ExtendedVector& direction,
    const ExtendedPopulationDeviceSystem& frame, const ExtendedFixedBiasUpdateReferences& config,
    const ExtendedSeriesResistorNewtonConfig* circuit, long double voltageDirection)
{
    const auto nodes = static_cast<Eigen::Index>(frame.electronFree.size());
    DensityUpdateAccumulator<long double> metric(config.digits);
    const auto densityToCm3 = state.concentrationScale_m3 * 1e-6L;
    for (Eigen::Index node = 0; node < nodes; ++node) {
        metric.add(0, state.coordinates[node] * state.potentialScale_V,
                   direction[node] * state.potentialScale_V, config.potential_V);
        if (frame.electronFree[node])
            metric.add(1, state.coordinates[nodes + node] * densityToCm3,
                       direction[nodes + node] * densityToCm3, config.electron_cm3);
        if (frame.holeFree[node])
            metric.add(2, state.coordinates[2 * nodes + node] * densityToCm3,
                       direction[2 * nodes + node] * densityToCm3, config.hole_cm3);
    }
    if (circuit)
        for (std::size_t i = 0; i < circuit->contactAliases; ++i)
            metric.add(3, state.deviceBias_V, voltageDirection, circuit->contactReference_V);
    return metric.finish();
}

ExtendedFixedBiasNewtonResult solveExtendedPopulationNewton(
    const ExtendedPopulationState& seed, const ExtendedFixedBiasNewtonConfig& config,
    const ExtendedSeriesResistorNewtonConfig* circuit,
    const ExtendedSeriesResistorAssembly& assemble,
    const ExtendedSeriesResistorAssembly& assembleAfterLocalProtection = {})
{
    const auto count = seed.coordinates.size();
    const auto nodes = count / 3;
    const auto validPositive = [](long double value) { return std::isfinite(value) && value > 0; };
    if (!assemble || count <= 0 || count % 3 || !seed.coordinates.allFinite() ||
        !validPositive(seed.potentialScale_V) || !validPositive(seed.concentrationScale_m3) ||
        !std::isfinite(seed.deviceBias_V) || !std::isfinite(seed.supplyBias_V) ||
        config.maxUpdates < 1 || config.maxBacktracks < 0 || config.refinementCorrections < 0 ||
        !validPositive(config.damping) || config.damping > 1 ||
        !validPositive(config.maxLinearBackwardError) || !validPositive(config.maxScaledLinearDefect) ||
        config.updateError.digits < 0 || config.updateError.digits > 15 ||
        !validPositive(config.updateError.potential_V) || !validPositive(config.updateError.electron_cm3) ||
        !validPositive(config.updateError.hole_cm3) ||
        (circuit && config.afterLocalProtectionJacobian) ||
        (config.rhsGrowth && (!circuit || !validPositive(config.rhsGrowth->factor) || config.rhsGrowth->factor < 1)))
        throw std::invalid_argument("Invalid extended fixed-bias Newton settings");
    validatePopulationUpdate(config.protection);

    ExtendedFixedBiasNewtonResult result;
    result.state = seed;
    auto frame = assemble(result.state);
    const auto validFrame = [&](const ExtendedSeriesResistorFrame& evaluation) {
        const auto& current = evaluation.device;
        if (circuit && (current.voltageDerivative.size() != count || !current.voltageDerivative.allFinite() ||
            evaluation.portDerivative.size() != count || !evaluation.portDerivative.allFinite() ||
            !std::isfinite(evaluation.current_A_per_um)))
            return false;
        return current.residual.size() == count && current.residual.allFinite() &&
            current.jacobian.rows() == count && current.jacobian.cols() == count &&
            std::ssize(current.electronFree) == nodes && std::ssize(current.holeFree) == nodes;
    };
    if (!validFrame(frame))
        throw std::invalid_argument("Invalid extended fixed-bias assembly");

    const ExtendedVector zero = ExtendedVector::Zero(count);
    const auto loadResidual = [&](const ExtendedPopulationState& state, long double current) {
        return circuit ? ArithmeticLoadLineResidual::evaluate(1.L, state.deviceBias_V,
            circuit->direction * circuit->resistance_ohm_um, current, state.supplyBias_V).voltage_V : 0.L;
    };
    ExtendedVector rhsGrowthScale;
    bool useRecoveryJacobian = false;
    for (int update = 0; update < config.maxUpdates; ++update) {
        if (useRecoveryJacobian && config.afterLocalProtectionJacobian) {
            auto replacement = config.afterLocalProtectionJacobian(result.state);
            if (replacement.rows() != count || replacement.cols() != count)
                throw std::invalid_argument("Invalid fixed-bias recovery Jacobian dimensions");
            frame.device.jacobian = std::move(replacement);
        }
        ExtendedVector columns = ExtendedVector::Ones(count + 1);
        for (Eigen::Index i = nodes; i < count; ++i)
            if (result.state.coordinates[i] > 0)
                columns[i] = result.state.coordinates[i];
        columns[count] = seed.potentialScale_V;
        auto eliminate = frame.device.prescribedPopulationRows;
        const auto load = loadResidual(result.state, frame.current_A_per_um);
        ExtendedVector circuitRow = circuit
            ? ExtendedVector((circuit->direction * circuit->resistance_ohm_um) * frame.portDerivative)
            : zero;
        auto system = detail::buildBorderedLoadLineSystem(frame.device.jacobian, frame.device.residual,
            circuit ? frame.device.voltageDerivative : zero, circuitRow, 1.L, load,
            columns, eliminate, result.failure);
        if (!system)
            return result;
        if (config.rhsGrowth && !result.rhsGrowth) {
            // Keep the initial augmented row scaling fixed throughout this supply trial.
            rhsGrowthScale = system->rowScale;
            const auto initial = system->rhs.cwiseAbs().maxCoeff();
            result.rhsGrowth = ExtendedRhsGrowthReport{initial, initial, std::nullopt};
        }
        const auto goal = (circuit || config.refineToLinearGoal) ? std::optional{ExtendedLinearAccuracyGoal{
            config.maxLinearBackwardError, config.maxScaledLinearDefect}} : std::nullopt;
        const auto linear = solveExtendedBorderedSystem(*system, eliminate,
            frame.device.retainedIndependentPotentialRows, config.refinementCorrections, goal,
            config.linearDefectMethod);
        result.factoredDimension = linear.factoredDimension;
        result.linearBackwardError = linear.backwardError;
        result.scaledLinearDefect = linear.scaledResidualInfinityNorm;
        if (!std::isfinite(result.linearBackwardError) ||
            result.linearBackwardError > config.maxLinearBackwardError ||
            !std::isfinite(result.scaledLinearDefect) ||
            result.scaledLinearDefect > config.maxScaledLinearDefect) {
            result.failure = "extended fixed-bias linear precision check failed";
            return result;
        }
        const auto scaledNorm = [&](const ExtendedVector& residual) {
            return residual.cwiseProduct(system->rowScale.head(count)).cwiseAbs().maxCoeff();
        };
        result.scaledResidual = scaledNorm(frame.device.residual);
        const ExtendedVector direction = linear.update.head(count);
        const auto residualBefore = result.scaledResidual;
        bool accepted = false;
        std::size_t protectedElectrons = 0, protectedHoles = 0;
        const int maxBacktracks = config.lineSearch ? config.maxBacktracks : 0;
        for (int attempt = 0; attempt <= maxBacktracks; ++attempt) {
            auto trial = result.state;
            BasicPopulationTrial<long double> protectedTrial;
            try {
                protectedTrial = makePopulationTrial(result.state.coordinates, direction,
                    std::ldexp(config.damping, -attempt), frame.device.electronFree, frame.device.holeFree,
                    config.protection, config.exactProtection ? &*config.exactProtection : nullptr);
                trial.coordinates = std::move(protectedTrial.state);
                if (circuit)
                    trial.deviceBias_V += protectedTrial.alpha * linear.update[count];
                const bool protectedWrite = !protectedTrial.protectedElectrons.empty()
                    || !protectedTrial.protectedHoles.empty();
                const auto& assembleNext = protectedWrite && assembleAfterLocalProtection
                    ? assembleAfterLocalProtection : assemble;
                auto next = assembleNext(trial);
                if (!validFrame(next))
                    continue;
                const auto merit = scaledNorm(next.device.residual);
                const auto nextLoad = loadResidual(trial, next.current_A_per_um);
                auto error = physicalUpdateError(trial, direction, next.device, config.updateError,
                                                 circuit, linear.update[count]);
                if (!std::isfinite(nextLoad))
                    continue;
                if (circuit)
                    trial.drainCurrent_A_per_um = next.current_A_per_um;
                const bool stationary = result.scaledResidual == 0 && merit == 0 && error.valid && error.rms == 0;
                if (!std::isfinite(merit) || (config.lineSearch && merit >= result.scaledResidual && !stationary))
                    continue;
                result.state = std::move(trial);
                frame = std::move(next);
                result.updateError = error;
                result.scaledResidual = merit;
                result.acceptedAlpha = protectedTrial.alpha;
                ++result.updates;
                protectedElectrons = protectedTrial.protectedElectrons.size();
                protectedHoles = protectedTrial.protectedHoles.size();
                if (auto& growth = result.rhsGrowth) {
                    const auto norm = std::max(
                        frame.device.residual.cwiseProduct(rhsGrowthScale.head(count)).cwiseAbs().maxCoeff(),
                        std::abs(nextLoad * rhsGrowthScale[count]));
                    growth->peakNorm = std::max(growth->peakNorm, norm);
                    if (!growth->firstExceededUpdate && growth->initialNorm > 0 &&
                        norm / growth->initialNorm > config.rhsGrowth->factor) {
                        growth->firstExceededUpdate = result.updates;
                        if (config.rhsGrowth->enforce) {
                            result.failure = "rhs_growth_limit";
                            return result;
                        }
                    }
                }
                accepted = true;
                break;
            } catch (const std::domain_error&) {
                // Re-evaluate the original tangent at the next trial damping.
            }
        }
        if (!accepted) {
            result.failure = "extended fixed-bias line search rejected update";
            return result;
        }
        if (config.observeAcceptedUpdate)
            config.observeAcceptedUpdate(result.state, {result.updates, useRecoveryJacobian,
                protectedElectrons, protectedHoles, result.acceptedAlpha, result.updateError.rms,
                residualBefore, result.scaledResidual, frame.device.residual});
        useRecoveryJacobian = (config.afterLocalProtectionJacobian || assembleAfterLocalProtection)
            && (protectedElectrons || protectedHoles);
        if (result.updateError.valid && result.updateError.rms < 1 &&
            (!circuit || !circuit->enforceLoadLine ||
             std::abs(loadResidual(result.state, frame.current_A_per_um)) <= circuit->tolerance_V)) {
            result.converged = true;
            return result;
        }
    }
    result.failure = "extended fixed-bias Newton iteration budget exhausted";
    return result;
}

} // namespace

ExtendedFixedBiasNewtonResult solveExtendedFixedBiasNewton(
    const ExtendedPopulationState& seed, const ExtendedFixedBiasNewtonConfig& config,
    const ExtendedFixedBiasAssembly& assemble)
{
    if (!assemble)
        throw std::invalid_argument("Missing extended fixed-bias assembly");
    return solveExtendedPopulationNewton(seed, config, nullptr, [&](const auto& state) {
        return ExtendedSeriesResistorFrame{assemble(state), {}, 0.L};
    });
}

ExtendedFixedBiasNewtonResult solveExtendedSeriesResistorNewton(
    const ExtendedPopulationState& seed, const ExtendedFixedBiasNewtonConfig& config,
    const ExtendedSeriesResistorNewtonConfig& circuit,
    const ExtendedSeriesResistorAssembly& assemble,
    const ExtendedSeriesResistorAssembly& assembleAfterLocalProtection)
{
    if (config.lineSearch || !std::isfinite(circuit.resistance_ohm_um) || circuit.resistance_ohm_um <= 0 ||
        !std::isfinite(circuit.tolerance_V) || circuit.tolerance_V <= 0 ||
        !std::isfinite(circuit.contactReference_V) || circuit.contactReference_V <= 0 ||
        !circuit.contactAliases || (circuit.direction != 1 && circuit.direction != -1) ||
        (assembleAfterLocalProtection && config.protection.protection != "local_nonpositive"))
        throw std::invalid_argument("Unsupported extended series-resistor Newton settings");
    return solveExtendedPopulationNewton(seed, config, &circuit, assemble, assembleAfterLocalProtection);
}

} // namespace vela
