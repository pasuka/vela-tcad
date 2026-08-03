#pragma once

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <optional>
#include <string>

namespace vela {

enum class RuntimeLogProfile {
    Minimal = 0,
    Default = 1,
    Debug = 2,
};

struct RuntimeLogConfig {
    bool enabled = true;
    std::string file;
    std::string level = "info";
    std::string flushLevel = "error";
    bool append = false;
    RuntimeLogProfile profile = RuntimeLogProfile::Default;
};

struct RuntimeLogCliOverrides {
    std::optional<bool> enabled;
    std::optional<std::string> file;
    std::optional<RuntimeLogProfile> profile;
};

class RuntimeLogSession {
public:
    RuntimeLogSession() = default;
    RuntimeLogSession(const RuntimeLogSession&) = delete;
    RuntimeLogSession& operator=(const RuntimeLogSession&) = delete;
    RuntimeLogSession(RuntimeLogSession&& other) noexcept;
    RuntimeLogSession& operator=(RuntimeLogSession&& other) noexcept;
    ~RuntimeLogSession();

    static RuntimeLogSession fromConfig(const nlohmann::json& cfg,
                                        const std::string& configFile,
                                        const std::string& simulationType);

    bool active() const { return active_; }
    void finish(bool success,
                std::size_t successCount = 0,
                std::size_t failureCount = 0);

private:
    explicit RuntimeLogSession(bool active, std::string simulationType);

    bool active_ = false;
    bool finished_ = false;
    std::string simulationType_;
    long long startTicksNs_ = 0;
};

RuntimeLogProfile runtimeLogProfileFromString(const std::string& value);
std::string toString(RuntimeLogProfile profile);

void setRuntimeLogCliOverrides(const RuntimeLogCliOverrides& overrides);
void clearRuntimeLogCliOverrides();

bool runtimeLogEnabled();
RuntimeLogProfile runtimeLogProfile();
bool runtimeLogAllows(RuntimeLogProfile required);

void runtimeLogInfo(const std::string& message);
void runtimeLogWarn(const std::string& message);
void runtimeLogError(const std::string& message);
void runtimeLogDebug(const std::string& message);

} // namespace vela

