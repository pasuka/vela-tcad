#pragma once

#include "vela/preprocess/PreprocessArtifacts.h"
#include <filesystem>
#include <string>

namespace vela {

class VelaPreprocessDriver {
public:
    PreprocessArtifacts buildIrOnly(const std::string& sdeScriptPath,
                                    const std::filesystem::path& outputDirectory) const;
};

} // namespace vela

