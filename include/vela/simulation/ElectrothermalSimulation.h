#pragma once
#include <nlohmann/json_fwd.hpp>
#include <iosfwd>
#include <filesystem>
#include <memory>
namespace vela {
/// Sweep-local immutable preparation cache. Each point still owns its state,
/// temperature preparation, neutral roots and nonlinear solver. Not thread-safe.
class ElectrothermalPreparationContext {
    struct Impl;
    std::shared_ptr<Impl> prepared_;
    friend nlohmann::json solveElectrothermalPoint(const nlohmann::json&,
        std::ostream&, ElectrothermalPreparationContext*);
};
/// Solve the explicit SI four-equation silicon configuration in process.
/// Same residual, analytic Jacobian and stopping rules as the qualified probe.
nlohmann::json solveElectrothermalPoint(const nlohmann::json& config, std::ostream& progress,
    ElectrothermalPreparationContext* preparation=nullptr);
/// Explicit electrothermal DC deck; paths resolve relative to the deck file.
nlohmann::json runElectrothermalSweep(const nlohmann::json& config,
    const std::filesystem::path& configFile);
}
