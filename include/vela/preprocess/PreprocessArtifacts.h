#pragma once

#include "vela/preprocess/DeviceIr2D.h"
#include <filesystem>

namespace vela {

struct PreprocessArtifacts {
    std::filesystem::path sourceScript;
    std::filesystem::path outputDirectory;
    DeviceIr2D ir;
};

} // namespace vela

