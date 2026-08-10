#pragma once

#include "vela/core/PhysicsCallCounters.h"

#include <chrono>
#include <cstdint>
#include <map>
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <string_view>
#include <vector>

namespace vela {

struct PerformanceProfilingConfig {
    bool enabled = false;
    std::string jsonFile = "performance_profile.json";
};

class PerformanceProfiler {
public:
    explicit PerformanceProfiler(PerformanceProfilingConfig config = {});

    bool enabled() const noexcept { return config_.enabled; }
    const std::string& jsonFile() const noexcept { return config_.jsonFile; }

    void recordStage(std::string_view name, std::chrono::nanoseconds elapsed);
    void increment(std::string_view name, std::uint64_t amount = 1);
    void observe(std::string_view name, double value);
    void recordNewtonSolve(bool converged,
                           int iterations,
                           double initialResidual,
                           double finalResidual,
                           std::string_view reason,
                           std::chrono::nanoseconds elapsed);

    nlohmann::json toJson() const;
    void writeJson() const;

private:
    struct StageData {
        std::uint64_t calls = 0;
        std::int64_t totalNanoseconds = 0;
        std::int64_t minimumNanoseconds = 0;
        std::int64_t maximumNanoseconds = 0;
        std::vector<std::int64_t> samples;
    };

    struct ObservationData {
        std::uint64_t count = 0;
        double total = 0.0;
        double minimum = 0.0;
        double maximum = 0.0;
        double last = 0.0;
    };

    struct NewtonSolveData {
        bool converged = false;
        int iterations = 0;
        double initialResidual = 0.0;
        double finalResidual = 0.0;
        std::string reason;
        std::int64_t elapsedNanoseconds = 0;
    };

    PerformanceProfilingConfig config_;
    std::chrono::steady_clock::time_point started_;
    std::map<std::string, StageData> stages_;
    std::map<std::string, std::uint64_t> counters_;
    std::map<std::string, ObservationData> observations_;
    std::vector<NewtonSolveData> newtonSolves_;
};

class ActivePerformanceProfilerScope {
public:
    explicit ActivePerformanceProfilerScope(PerformanceProfiler* profiler) noexcept;
    ~ActivePerformanceProfilerScope();

    ActivePerformanceProfilerScope(const ActivePerformanceProfilerScope&) = delete;
    ActivePerformanceProfilerScope& operator=(const ActivePerformanceProfilerScope&) = delete;

private:
    PerformanceProfiler* previous_ = nullptr;
};

class ScopedPerformanceTimer {
public:
    explicit ScopedPerformanceTimer(std::string_view stage) noexcept;
    ~ScopedPerformanceTimer();

    ScopedPerformanceTimer(const ScopedPerformanceTimer&) = delete;
    ScopedPerformanceTimer& operator=(const ScopedPerformanceTimer&) = delete;

private:
    PerformanceProfiler* profiler_ = nullptr;
    std::string_view stage_;
    std::chrono::steady_clock::time_point started_{};
    /// Counter snapshot taken on entry; the exit delta is the stage-inclusive
    /// carrier-statistics call count for this scope.
    PhysicsCallCounters startedCounters_{};
};

PerformanceProfiler* activePerformanceProfiler() noexcept;
void incrementPerformanceCounter(std::string_view name,
                                 std::uint64_t amount = 1) noexcept;
void observePerformanceValue(std::string_view name, double value) noexcept;
PerformanceProfilingConfig performanceProfilingConfigFromJson(
    const nlohmann::json& solverJson);

} // namespace vela
