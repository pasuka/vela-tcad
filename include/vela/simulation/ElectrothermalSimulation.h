#pragma once
#include <nlohmann/json_fwd.hpp>
#include <iosfwd>
#include <filesystem>
namespace vela {
/// Solve the explicit SI four-equation silicon configuration in process.
/// Same residual, analytic Jacobian and stopping rules as the qualified probe.
nlohmann::json solveElectrothermalPoint(const nlohmann::json& config, std::ostream& progress);
/// Explicit electrothermal DC deck; paths resolve relative to the deck file.
nlohmann::json runElectrothermalSweep(const nlohmann::json& config,
    const std::filesystem::path& configFile);
}
