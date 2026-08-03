#include "vela/core/RuntimeLog.h"

#include <nlohmann/json.hpp>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace vela {

namespace {

struct RuntimeLogState {
    std::mutex mutex;
    std::shared_ptr<spdlog::logger> logger;
    RuntimeLogProfile profile = RuntimeLogProfile::Default;
    RuntimeLogCliOverrides cliOverrides;
};

RuntimeLogState& runtimeLogState()
{
    static RuntimeLogState state;
    return state;
}

long long monotonicNowNs()
{
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();
}

std::string hostName()
{
    const char* value = std::getenv("COMPUTERNAME");
    if (value && *value)
        return value;
    value = std::getenv("HOSTNAME");
    if (value && *value)
        return value;
    return "unknown";
}

int processId()
{
#if defined(_WIN32)
    return _getpid();
#else
    return static_cast<int>(getpid());
#endif
}

spdlog::level::level_enum spdlogLevelFromString(const std::string& value)
{
    if (value == "trace") return spdlog::level::trace;
    if (value == "debug") return spdlog::level::debug;
    if (value == "info") return spdlog::level::info;
    if (value == "warn") return spdlog::level::warn;
    if (value == "error") return spdlog::level::err;
    if (value == "critical") return spdlog::level::critical;
    if (value == "off") return spdlog::level::off;
    throw std::invalid_argument("RuntimeLog: unknown log level '" + value + "'.");
}

std::filesystem::path defaultLogPath(const std::string& configFile)
{
    const std::filesystem::path cfgPath(configFile);
    const std::filesystem::path parent = cfgPath.parent_path().empty()
        ? std::filesystem::current_path()
        : cfgPath.parent_path();
    return parent / (cfgPath.stem().string() + ".log");
}

RuntimeLogConfig runtimeLogConfigFromJson(const nlohmann::json& cfg,
                                          const std::string& configFile)
{
    RuntimeLogConfig runtime;
    if (cfg.contains("runtime_log")) {
        const auto& node = cfg.at("runtime_log");
        if (!node.is_object())
            throw std::invalid_argument("RuntimeLog: runtime_log must be an object.");
        runtime.enabled = node.value("enabled", runtime.enabled);
        runtime.file = node.value("file", runtime.file);
        runtime.level = node.value("level", runtime.level);
        runtime.flushLevel = node.value("flush_level", runtime.flushLevel);
        runtime.append = node.value("append", runtime.append);
        runtime.profile = runtimeLogProfileFromString(
            node.value("profile", toString(runtime.profile)));
    }

    if (runtime.file.empty())
        runtime.file = defaultLogPath(configFile).string();
    else {
        std::filesystem::path filePath(runtime.file);
        if (filePath.is_relative()) {
            const std::filesystem::path cfgPath(configFile);
            const std::filesystem::path baseDir = cfgPath.parent_path().empty()
                ? std::filesystem::current_path()
                : cfgPath.parent_path();
            filePath = baseDir / filePath;
            runtime.file = filePath.string();
        }
    }

    {
        RuntimeLogState& state = runtimeLogState();
        std::lock_guard<std::mutex> lock(state.mutex);
        if (state.cliOverrides.enabled.has_value())
            runtime.enabled = *state.cliOverrides.enabled;
        if (state.cliOverrides.file.has_value())
            runtime.file = *state.cliOverrides.file;
        if (state.cliOverrides.profile.has_value())
            runtime.profile = *state.cliOverrides.profile;
    }

    // Validate levels even if disabled so user gets immediate config feedback.
    (void)spdlogLevelFromString(runtime.level);
    (void)spdlogLevelFromString(runtime.flushLevel);
    return runtime;
}

void runtimeLog(spdlog::level::level_enum level, const std::string& message)
{
    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    if (state.logger)
        state.logger->log(level, message);
}

} // namespace

RuntimeLogProfile runtimeLogProfileFromString(const std::string& value)
{
    if (value == "minimal")
        return RuntimeLogProfile::Minimal;
    if (value == "default")
        return RuntimeLogProfile::Default;
    if (value == "debug")
        return RuntimeLogProfile::Debug;
    throw std::invalid_argument("RuntimeLog: unknown profile '" + value + "'.");
}

std::string toString(RuntimeLogProfile profile)
{
    switch (profile) {
        case RuntimeLogProfile::Minimal:
            return "minimal";
        case RuntimeLogProfile::Default:
            return "default";
        case RuntimeLogProfile::Debug:
            return "debug";
    }
    return "default";
}

void setRuntimeLogCliOverrides(const RuntimeLogCliOverrides& overrides)
{
    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    state.cliOverrides = overrides;
}

void clearRuntimeLogCliOverrides()
{
    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    state.cliOverrides = RuntimeLogCliOverrides{};
}

RuntimeLogSession::RuntimeLogSession(bool active, std::string simulationType)
    : active_(active),
      simulationType_(std::move(simulationType)),
      startTicksNs_(monotonicNowNs())
{
}

RuntimeLogSession RuntimeLogSession::fromConfig(const nlohmann::json& cfg,
                                                const std::string& configFile,
                                                const std::string& simulationType)
{
    const RuntimeLogConfig runtimeCfg = runtimeLogConfigFromJson(cfg, configFile);
    if (!runtimeCfg.enabled)
        return RuntimeLogSession(false, simulationType);

    const std::filesystem::path filePath(runtimeCfg.file);
    if (!filePath.parent_path().empty())
        std::filesystem::create_directories(filePath.parent_path());

    auto sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
        filePath.string(), runtimeCfg.append);
    auto logger = std::make_shared<spdlog::logger>("vela.runtime", sink);
    logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
    logger->set_level(spdlogLevelFromString(runtimeCfg.level));
    logger->flush_on(spdlogLevelFromString(runtimeCfg.flushLevel));

    RuntimeLogState& state = runtimeLogState();
    {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.logger = std::move(logger);
        state.profile = runtimeCfg.profile;
    }

    RuntimeLogSession session(true, simulationType);
    runtimeLogInfo("=== Vela Runtime Log ===");
    runtimeLogInfo("version: " VELA_VERSION);
    runtimeLogInfo("simulation_type: " + simulationType);
    runtimeLogInfo("config_file: " + std::filesystem::absolute(configFile).string());
    runtimeLogInfo("log_file: " + std::filesystem::absolute(filePath).string());
    runtimeLogInfo("host: " + hostName());
    runtimeLogInfo("pid: " + std::to_string(processId()));
    runtimeLogInfo("profile: " + toString(runtimeCfg.profile));
    return session;
}

RuntimeLogSession::RuntimeLogSession(RuntimeLogSession&& other) noexcept
    : active_(other.active_),
      finished_(other.finished_),
      simulationType_(std::move(other.simulationType_)),
      startTicksNs_(other.startTicksNs_)
{
    other.active_ = false;
    other.finished_ = true;
}

RuntimeLogSession& RuntimeLogSession::operator=(RuntimeLogSession&& other) noexcept
{
    if (this == &other)
        return *this;
    if (active_ && !finished_)
        finish(false);
    active_ = other.active_;
    finished_ = other.finished_;
    simulationType_ = std::move(other.simulationType_);
    startTicksNs_ = other.startTicksNs_;
    other.active_ = false;
    other.finished_ = true;
    return *this;
}

RuntimeLogSession::~RuntimeLogSession()
{
    if (active_ && !finished_)
        finish(false);
}

void RuntimeLogSession::finish(bool success,
                               std::size_t successCount,
                               std::size_t failureCount)
{
    if (!active_ || finished_)
        return;

    const long long elapsedNs = monotonicNowNs() - startTicksNs_;
    const double elapsedS = static_cast<double>(elapsedNs) * 1.0e-9;
    runtimeLogInfo("status: " + std::string(success ? "success" : "failed"));
    runtimeLogInfo("elapsed_seconds: " + std::to_string(elapsedS));
    if (successCount > 0 || failureCount > 0) {
        runtimeLogInfo("success_count: " + std::to_string(successCount));
        runtimeLogInfo("failure_count: " + std::to_string(failureCount));
    }
    runtimeLogInfo("=== End Runtime Log ===");

    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    if (state.logger) {
        state.logger->flush();
        state.logger.reset();
    }
    finished_ = true;
}

bool runtimeLogEnabled()
{
    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    return static_cast<bool>(state.logger);
}

RuntimeLogProfile runtimeLogProfile()
{
    RuntimeLogState& state = runtimeLogState();
    std::lock_guard<std::mutex> lock(state.mutex);
    return state.profile;
}

bool runtimeLogAllows(RuntimeLogProfile required)
{
    return static_cast<int>(runtimeLogProfile()) >= static_cast<int>(required);
}

void runtimeLogInfo(const std::string& message)
{
    runtimeLog(spdlog::level::info, message);
}

void runtimeLogWarn(const std::string& message)
{
    runtimeLog(spdlog::level::warn, message);
}

void runtimeLogError(const std::string& message)
{
    runtimeLog(spdlog::level::err, message);
}

void runtimeLogDebug(const std::string& message)
{
    runtimeLog(spdlog::level::debug, message);
}

} // namespace vela

