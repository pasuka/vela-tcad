#pragma once
#include <nlohmann/json_fwd.hpp>
#include <iosfwd>
#include <filesystem>
#include <memory>
namespace vela {
namespace experimental { class ElectrothermalDirectSolver; }
/// Sweep-local immutable preparation cache. Each point still owns its state,
/// temperature preparation, neutral roots and nonlinear state. Linear analysis may
/// persist, but each Newton matrix is numerically refactorized. Not thread-safe.
class ElectrothermalPreparationContext {
public:
    void clear() { prepared_.reset(); clearLinearContext(); }
    void clearLinearContext() { linear_.reset(); }
private:
    struct Impl;
    std::shared_ptr<Impl> prepared_;
    std::shared_ptr<experimental::ElectrothermalDirectSolver> linear_;
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
