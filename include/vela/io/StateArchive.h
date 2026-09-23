#pragma once

#include <cstdint>
#include <filesystem>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace vela {

// Physical SI fields, independent of solver coordinates and legacy codecs.
// Each field has nodeCount entries in the mesh's original node order.
struct StateArchive {
    std::uint64_t nodeCount = 0;
    std::map<std::string, std::vector<double>> fields;
    nlohmann::json metadata;
};

// Required metadata: mode (dd/electrothermal), mesh_sha256, potential_origin_V.
// Validates without altering values (including independent QF low bits).
void validateStateArchive(const StateArchive& state);

// The caller supplies an already prepared mesh identity; no mesh rehash per IO.
StateArchive readStateArchive(const std::filesystem::path& path,
                              std::uint64_t expectedNodes,
                              const std::string& expectedMeshSha256);

// Same-directory temporary, close, then checked replacement. A failed write
// leaves the previous destination intact. Does not claim power-loss durability.
void writeStateArchive(const std::filesystem::path& path, const StateArchive& state);

} // namespace vela
