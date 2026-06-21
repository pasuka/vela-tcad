#pragma once

#include "vela/core/Types.h"
#include "vela/solver/GummelSolver.h"

#include <filesystem>

namespace vela {

DDSolution readDDSolutionStateCsv(const std::filesystem::path& path,
                                  Index expectedNodeCount);

void writeDDSolutionStateCsv(const std::filesystem::path& path,
                             const DDSolution& solution);

} // namespace vela
