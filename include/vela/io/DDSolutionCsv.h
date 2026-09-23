#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/solver/GummelSolver.h"

#include <filesystem>
#include <string>

namespace vela {

// Production restart I/O uses vela.state/2 HDF5 with DDStateArchiveScope.
// Explicit CSV helpers below are migration/import utilities only.
DDSolution readDDSolutionState(const std::filesystem::path& path,
                              Index expectedNodeCount,
                              UnitScalingConfig scaling = UnitScalingConfig{});
void writeDDSolutionState(const std::filesystem::path& path,
                          const DDSolution& solution,
                          UnitScalingConfig scaling = UnitScalingConfig{});

DDSolution readDDSolutionStateCsv(
    const std::filesystem::path& path,
    Index expectedNodeCount,
    UnitScalingConfig scaling = UnitScalingConfig{});

void writeDDSolutionStateCsv(
    const std::filesystem::path& path,
    const DDSolution& solution,
    UnitScalingConfig scaling = UnitScalingConfig{});

} // namespace vela
