#include "vela/preprocess/VelaPreprocessDriver.h"

#include "vela/io/SdeScriptReader.h"
#include <filesystem>

namespace vela {

PreprocessArtifacts VelaPreprocessDriver::buildIrOnly(const std::string& sdeScriptPath,
                                                      const std::filesystem::path& outputDirectory) const
{
    SdeScriptReader reader;
    PreprocessArtifacts artifacts;
    artifacts.sourceScript = sdeScriptPath;
    artifacts.outputDirectory = outputDirectory;
    artifacts.ir = reader.parseFile(sdeScriptPath);
    std::filesystem::create_directories(outputDirectory);
    return artifacts;
}

} // namespace vela

