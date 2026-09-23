#include "vela/io/DDSolutionState.h"
#include "vela/io/DDSolutionCsv.h"
#include <stdexcept>

namespace vela {
namespace {
const nlohmann::json& context() {
    const auto* value = activeDDStateArchiveMetadata();
    if (!value) throw std::logic_error("DD state I/O requires an explicit mesh and reference context");
    return *value;
}
}
DDSolution readDDSolutionState(const std::filesystem::path& path, Index expected,
                              UnitScalingConfig scaling) {
    const auto& metadata = context();
    const auto state = readStateArchive(path, expected, metadata.at("mesh_sha256").get<std::string>());
    if (state.metadata.at("potential_origin_V") != metadata.at("potential_origin_V"))
        throw std::invalid_argument("HDF5 state reference origin differs; translate the seed explicitly");
    return restoreDDSolution(state, scaling);
}
void writeDDSolutionState(const std::filesystem::path& path, const DDSolution& solution,
                          UnitScalingConfig scaling) {
    writeStateArchive(path, archiveDDSolution(solution, context(), scaling));
}
} // namespace vela
