#pragma once
#include "vela/io/StateArchive.h"

namespace vela {
// Replace only state vectors with HDF5 references. The JSON control record is
// committed by the caller after these immutable state files have been closed.
nlohmann::json packElectrothermalRecord(const std::filesystem::path& recordPath,
    const nlohmann::json& record, const nlohmann::json& metadata);
// Requires a mesh identity computed independently from the current input mesh.
nlohmann::json unpackElectrothermalRecord(const std::filesystem::path& recordPath,
    const nlohmann::json& record, std::uint64_t nodes, const std::string& meshSha256,
    double potentialOrigin);
} // namespace vela
