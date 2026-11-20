#pragma once

#include "vela/solver/ExtendedPopulationState.h"
#include <cstdint>
#include <filesystem>
#include <string>

namespace vela {

enum class ExtendedPopulationArchiveRole { Accepted, Diagnostic };

struct ExtendedPopulationArchiveIdentity {
    std::string meshSha256;
    std::string inputSha256;
    ExtendedPopulationArchiveRole role;
};

// Requires HDF5 StateArchive support and a binary 64-significand-bit long double.
// The schema is distinct from legacy StateArchive; a binary64 archive cannot
// silently become an authoritative extended continuation state.
// An accepted archive is reopened and checked exactly before atomic replacement.
void writeExtendedPopulationArchive(const std::filesystem::path& path,
                                    const ExtendedPopulationState& state,
                                    const ExtendedPopulationArchiveIdentity& identity);
ExtendedPopulationState readExtendedPopulationArchive(const std::filesystem::path& path,
                                                      std::uint64_t expectedNodes,
                                                      const ExtendedPopulationArchiveIdentity& expectedIdentity);

} // namespace vela
