#pragma once

#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

// Shared numerical policy for explicitly enabled coupled continuation.
// A resistor's external voltage span does not define a device-voltage span.
namespace vela::detail {

inline constexpr double defaultPredictedInnerStep_V = 0.02;

struct InnerStepLimit {
    double volts = defaultPredictedInnerStep_V;
    std::string source = "fallback";
    std::string reason;
};

inline double stepNumber(const nlohmann::json& object, const char* key)
{
    const auto& value = object.at(key);
    if (!value.is_number() || value.is_boolean())
        throw std::invalid_argument(std::string(key) + " must be a finite number");
    const double result = value.get<double>();
    if (!std::isfinite(result))
        throw std::invalid_argument(std::string(key) + " must be finite");
    return result;
}

inline double positiveStepNumber(const nlohmann::json& object, const char* key)
{
    const double result = stepNumber(object, key);
    if (!(result > 0.))
        throw std::invalid_argument(std::string(key) + " must be positive");
    return result;
}

inline const nlohmann::json& stepObject(const nlohmann::json& parent, const char* key)
{
    static const nlohmann::json empty = nlohmann::json::object();
    if (!parent.contains(key)) return empty;
    const auto& result = parent.at(key);
    if (!result.is_object())
        throw std::invalid_argument(std::string(key) + " must be an object");
    return result;
}

inline InnerStepLimit resolveInnerStepLimit(const nlohmann::json& config,
                                            const nlohmann::json& manifest)
{
    const auto& sweep = stepObject(config, "sweep");
    const auto& circuit = stepObject(sweep, "external_circuit");
    const bool external = !circuit.empty() && circuit.value("enabled", true);
    if (circuit.contains("max_inner_voltage_step_V"))
        positiveStepNumber(circuit, "max_inner_voltage_step_V");

    std::string reason = "no explicit device-voltage step limit";
    if (manifest.contains("inner_step_control")) {
        const auto& control = stepObject(manifest, "inner_step_control");
        for (const auto& item : control.items()) {
            if (item.key() != "max_step_V" && item.key() != "max_step_fraction" &&
                item.key() != "initial_voltage_V" && item.key() != "goal_voltage_V")
                throw std::invalid_argument("unknown inner_step_control field: " + item.key());
        }
        if (control.contains("max_step_V"))
            return {positiveStepNumber(control, "max_step_V"),
                    "manifest.inner_step_control.max_step_V", "explicit inner volts"};
        for (const char* key : {"initial_voltage_V", "goal_voltage_V"})
            if (control.contains(key)) stepNumber(control, key);
        if (control.contains("max_step_fraction")) {
            const double fraction = positiveStepNumber(control, "max_step_fraction");
            if (fraction > 1.)
                throw std::invalid_argument("max_step_fraction must not exceed one");
            std::optional<double> start, goal;
            if (control.contains("initial_voltage_V")) start = stepNumber(control, "initial_voltage_V");
            if (control.contains("goal_voltage_V")) goal = stepNumber(control, "goal_voltage_V");
            // A direct voltage sweep has the same device and source coordinates.
            // Never borrow these endpoints from an external-resistor sweep.
            if (!external) {
                if (!start && sweep.contains("start")) start = stepNumber(sweep, "start");
                if (!goal && sweep.contains("stop")) goal = stepNumber(sweep, "stop");
            }
            if (start && goal) {
                const double span = std::abs(*goal - *start);
                const double volts = span * fraction;
                if (!std::isfinite(span) || !std::isfinite(volts))
                    throw std::invalid_argument("derived inner voltage step overflows");
                if (volts > 0.)
                    return {volts, "manifest.inner_step_control.max_step_fraction * inner_span_V",
                            "fraction times absolute device-voltage span"};
                reason = "device-voltage span is zero or not representable";
            } else {
                reason = "normalized step lacks a complete device-voltage span";
            }
        } else {
            reason = "inner_step_control does not specify a step size";
        }
    }
    // Read presence from the original JSON: the ordinary boundary solver's
    // implicit 0.1 V default must not become an explicit predictor setting.
    if (circuit.contains("max_inner_voltage_step_V"))
        return {positiveStepNumber(circuit, "max_inner_voltage_step_V"),
                "sweep.external_circuit.max_inner_voltage_step_V", "explicit device-voltage scale"};
    if (!external && sweep.contains("max_step"))
        return {positiveStepNumber(sweep, "max_step"), "sweep.max_step",
                "direct-voltage sweep; already converted to volts"};
    if (external && !manifest.contains("inner_step_control"))
        reason = "external-voltage step controls do not determine an inner-voltage limit";
    return {defaultPredictedInnerStep_V, "fallback", reason};
}

struct OuterStepControl {
    double initial = 0.01, minimum = 1e-5, maximum = 2500.;
    double growth = 1.41, shrink = 0.5;
    int maxRetries = 20;
    std::string growthMode = "newton_iterations";
};

inline OuterStepControl resolveOuterStepControl(const nlohmann::json& config,
                                               const nlohmann::json& manifest)
{
    const auto& sweep = stepObject(config, "sweep");
    const auto& circuit = stepObject(sweep, "external_circuit");
    OuterStepControl result;
    const auto choose = [&](const char* specific, const char* general, double fallback, bool allowZero) {
        if (circuit.contains(specific)) {
            const double value = stepNumber(circuit, specific);
            if (value > 0.) return value;
            if (!allowZero || value < 0.)
                throw std::invalid_argument(std::string(specific) + " must be positive");
        }
        return sweep.contains(general) ? positiveStepNumber(sweep, general) : fallback;
    };
    result.minimum = choose("coupled_min_outer_step_V", "min_step", result.minimum, false);
    result.maximum = choose("coupled_max_outer_step_V", "max_step", result.maximum, true);
    result.initial = manifest.contains("initial_step_V")
        ? positiveStepNumber(manifest, "initial_step_V")
        : choose("coupled_initial_outer_step_V", "initial_step", result.initial, true);
    result.growth = choose("coupled_outer_growth_factor", "growth_factor", result.growth, false);
    result.shrink = choose("coupled_outer_shrink_factor", "shrink_factor", result.shrink, false);
    result.growthMode = sweep.value("step_growth_mode", result.growthMode);
    if (circuit.contains("coupled_max_step_retries")) {
        if (!circuit.at("coupled_max_step_retries").is_number_integer())
            throw std::invalid_argument("coupled_max_step_retries must be an integer");
        const auto value = circuit.at("coupled_max_step_retries").get<long long>();
        if (value < 0 || value > std::numeric_limits<int>::max())
            throw std::invalid_argument("coupled_max_step_retries is out of range");
        result.maxRetries = static_cast<int>(value);
    }
    if (result.minimum > result.maximum || result.initial < result.minimum ||
        result.initial > result.maximum || result.growth < 1. || result.shrink >= 1. ||
        (result.growthMode != "fixed" && result.growthMode != "newton_iterations"))
        throw std::invalid_argument("inconsistent outer adaptive-step controls");
    return result;
}

inline std::optional<double> guardedOuterStep(double proposed, const OuterStepControl& control,
                                               double innerLimit, double deltaInner, double deltaOuter)
{
    if (!std::isfinite(proposed) || proposed <= 0. || !std::isfinite(innerLimit) || innerLimit <= 0. ||
        !std::isfinite(deltaInner) || !std::isfinite(deltaOuter) || deltaOuter == 0.)
        throw std::invalid_argument("invalid inner-step predictor inputs");
    double allowed = std::min(proposed, control.maximum);
    if (deltaInner != 0.) {
        // Long double avoids intermediate overflow for a very small secant slope.
        const long double bound = static_cast<long double>(innerLimit) *
            std::abs(static_cast<long double>(deltaOuter)) / std::abs(static_cast<long double>(deltaInner));
        if (bound < allowed) allowed = static_cast<double>(bound);
    }
    // Do not silently exceed the requested inner limit by clamping upward.
    if (allowed < control.minimum) return std::nullopt;
    return allowed;
}

} // namespace vela::detail
