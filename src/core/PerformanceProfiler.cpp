#include "vela/core/PerformanceProfiler.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <utility>

namespace vela {

thread_local PhysicsCallCounters physicsCallCounters{};

namespace {

thread_local PerformanceProfiler* activeProfiler = nullptr;

std::int64_t percentile(std::vector<std::int64_t> values, double fraction)
{
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const std::size_t index = std::min(
        values.size() - 1,
        static_cast<std::size_t>(
            std::ceil(fraction * static_cast<double>(values.size()))) - 1);
    return values[index];
}

/**
 * @brief Emit the carrier-statistics call counts issued inside one stage.
 *
 * The emitted counters are stage-inclusive: a call made inside
 * `jacobian.edge_physics` is also counted for the enclosing `dd.jacobian`
 * scope, exactly like the stage timings.  They must therefore not be summed
 * across parent and child stages.
 */
void recordPhysicsCallDeltas(PerformanceProfiler& profiler,
                             std::string_view stage,
                             const PhysicsCallCounters& started)
{
    if (!profiler.enabled())
        return;
    const PhysicsCallCounters& current = physicsCallCounters;
    const auto emit = [&](std::string_view suffix,
                          std::uint64_t end) {
        if (end == begin)
            return;
        std::string name;
        name.reserve(stage.size() + suffix.size() + 1);
        name.append(stage);
        name.push_back('.');
        name.append(suffix);
        profiler.increment(name, end - begin);
    };
    emit("fermi_dirac_half_calls",
         started.fermiDiracHalf, current.fermiDiracHalf);
    emit("fermi_dirac_half_derivative_calls",
         started.fermiDiracHalfDerivative, current.fermiDiracHalfDerivative);
    emit("inverse_fermi_dirac_half_calls",
         started.inverseFermiDiracHalf, current.inverseFermiDiracHalf);
    emit("equilibrium_state_solves",
         started.equilibriumStateSolves, current.equilibriumStateSolves);
    emit("equilibrium_state_iterations",
         started.equilibriumStateIterations,
         current.equilibriumStateIterations);
}

} // namespace

PerformanceProfiler::PerformanceProfiler(PerformanceProfilingConfig config)
    : config_(std::move(config))
    , started_(std::chrono::steady_clock::now())
{}

void PerformanceProfiler::recordStage(std::string_view name,
                                      std::chrono::nanoseconds elapsed)
{
    if (!enabled())
        return;
    StageData& stage = stages_[std::string(name)];
    const std::int64_t value = elapsed.count();
    ++stage.calls;
    stage.totalNanoseconds += value;
    if (stage.calls == 1) {
        stage.minimumNanoseconds = value;
        stage.maximumNanoseconds = value;
    } else {
        stage.minimumNanoseconds = std::min(stage.minimumNanoseconds, value);
        stage.maximumNanoseconds = std::max(stage.maximumNanoseconds, value);
    }
    stage.samples.push_back(value);
}

void PerformanceProfiler::increment(std::string_view name, std::uint64_t amount)
{
    if (enabled())
        counters_[std::string(name)] += amount;
}

void PerformanceProfiler::observe(std::string_view name, double value)
{
    if (!enabled())
        return;
    ObservationData& observation = observations_[std::string(name)];
    ++observation.count;
    observation.total += value;
    observation.last = value;
    if (observation.count == 1) {
        observation.minimum = value;
        observation.maximum = value;
    } else {
        observation.minimum = std::min(observation.minimum, value);
        observation.maximum = std::max(observation.maximum, value);
    }
}

void PerformanceProfiler::recordNewtonSolve(bool converged,
                                             int iterations,
                                             double initialResidual,
                                             double finalResidual,
                                             std::string_view reason,
                                             std::chrono::nanoseconds elapsed)
{
    if (!enabled())
        return;
    newtonSolves_.push_back({
        converged,
        iterations,
        initialResidual,
        finalResidual,
        std::string(reason),
        elapsed.count()});
    increment("newton.solve_calls");
    increment("newton.updates", static_cast<std::uint64_t>(std::max(iterations, 0)));
}

nlohmann::json PerformanceProfiler::toJson() const
{
    nlohmann::json stages = nlohmann::json::array();
    for (const auto& [name, stage] : stages_) {
        stages.push_back({
            {"name", name},
            {"calls", stage.calls},
            {"total_ns", stage.totalNanoseconds},
            {"average_ns", stage.calls == 0
                ? 0.0
                : static_cast<double>(stage.totalNanoseconds) /
                    static_cast<double>(stage.calls)},
            {"p50_ns", percentile(stage.samples, 0.50)},
            {"p95_ns", percentile(stage.samples, 0.95)},
            {"min_ns", stage.minimumNanoseconds},
            {"max_ns", stage.maximumNanoseconds}
        });
    }

    nlohmann::json counters = nlohmann::json::object();
    for (const auto& [name, value] : counters_)
        counters[name] = value;

    nlohmann::json observations = nlohmann::json::object();
    for (const auto& [name, value] : observations_) {
        observations[name] = {
            {"count", value.count},
            {"last", value.last},
            {"min", value.minimum},
            {"max", value.maximum},
            {"average", value.count == 0
                ? 0.0
                : value.total / static_cast<double>(value.count)}
        };
    }

    nlohmann::json newtonSolves = nlohmann::json::array();
    for (std::size_t index = 0; index < newtonSolves_.size(); ++index) {
        const NewtonSolveData& solve = newtonSolves_[index];
        newtonSolves.push_back({
            {"index", index + 1},
            {"converged", solve.converged},
            {"iterations", solve.iterations},
            {"initial_residual", solve.initialResidual},
            {"final_residual", solve.finalResidual},
            {"reason", solve.reason},
            {"elapsed_ns", solve.elapsedNanoseconds}
        });
    }

    return {
        {"schema_version", 1},
        {"clock", "steady_clock"},
        {"enabled", enabled()},
        {"elapsed_ns", std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - started_).count()},
        {"stages", std::move(stages)},
        {"counters", std::move(counters)},
        {"observations", std::move(observations)},
        {"newton_solves", std::move(newtonSolves)}
    };
}

void PerformanceProfiler::writeJson() const
{
    if (!enabled())
        return;
    const std::filesystem::path path(config_.jsonFile);
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::trunc);
    if (!output.is_open())
        throw std::runtime_error(
            "PerformanceProfiler: cannot open output file: " + path.string());
    output << toJson().dump(2) << '\n';
}

ActivePerformanceProfilerScope::ActivePerformanceProfilerScope(
    PerformanceProfiler* profiler) noexcept
    : previous_(activeProfiler)
{
    activeProfiler = profiler;
}

ActivePerformanceProfilerScope::~ActivePerformanceProfilerScope()
{
    activeProfiler = previous_;
}

ScopedPerformanceTimer::ScopedPerformanceTimer(std::string_view stage) noexcept
    : profiler_(activePerformanceProfiler())
    , stage_(stage)
{
    if (profiler_ != nullptr) {
        startedCounters_ = physicsCallCounters;
        started_ = std::chrono::steady_clock::now();
    }
}

ScopedPerformanceTimer::~ScopedPerformanceTimer()
{
    if (profiler_ != nullptr) {
        profiler_->recordStage(
            stage_,
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started_));
        recordPhysicsCallDeltas(*profiler_, stage_, startedCounters_);
    }
}

PerformanceProfiler* activePerformanceProfiler() noexcept
{
    return activeProfiler;
}

void incrementPerformanceCounter(std::string_view name,
                                 std::uint64_t amount) noexcept
{
    if (activeProfiler != nullptr)
        activeProfiler->increment(name, amount);
}

void observePerformanceValue(std::string_view name, double value) noexcept
{
    if (activeProfiler != nullptr)
        activeProfiler->observe(name, value);
}

PerformanceProfilingConfig performanceProfilingConfigFromJson(
    const nlohmann::json& solverJson)
{
    PerformanceProfilingConfig config;
    if (!solverJson.contains("performance_profiling"))
        return config;
    const nlohmann::json& value = solverJson.at("performance_profiling");
    if (value.is_boolean()) {
        config.enabled = value.get<bool>();
        return config;
    }
    if (!value.is_object()) {
        throw std::invalid_argument(
            "solver.performance_profiling must be a boolean or object.");
    }
    config.enabled = value.value("enabled", true);
    config.jsonFile = value.value("json_file", config.jsonFile);
    if (config.enabled && config.jsonFile.empty()) {
        throw std::invalid_argument(
            "solver.performance_profiling.json_file must not be empty when enabled.");
    }
    return config;
}

} // namespace vela
