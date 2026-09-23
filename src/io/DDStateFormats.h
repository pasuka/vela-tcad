#pragma once
#include <filesystem>
#include <string>
#include "vela/core/Types.h"

namespace vela {
// Internal canonical bridge; files contain actual typed arrays, not VDS blobs.
std::string readPublicDDState(const std::filesystem::path&, Index expected);
void writePublicDDState(const std::filesystem::path&, const std::string& canonical);
}
