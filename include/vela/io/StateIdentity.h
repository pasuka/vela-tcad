#pragma once
#include "vela/mesh/DeviceMesh.h"
#include "vela/core/UnitScaling.h"
#include <string>
#include <string_view>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace vela {
std::string stateSha256(std::string_view bytes);
// Canonical little-endian encoding of node order, coordinates, topology,
// regions, contacts and length units. Excludes derived solver caches.
std::string stateMeshIdentity(const DeviceMesh& mesh, UnitScalingConfig scaling = UnitScalingConfig{});
std::string stateMeshIdentity(const DeviceMesh& mesh, double lengthMPerInternal);
// Input-only source hashes, prepared once per request, never at each state write.
nlohmann::json stateInputProvenance(const nlohmann::json& config,
                                  const std::filesystem::path& directory);
} // namespace vela
