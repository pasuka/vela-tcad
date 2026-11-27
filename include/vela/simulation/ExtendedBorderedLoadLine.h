#pragma once

#include "vela/solver/ExtendedSparseLu.h"
#include <span>
#include <optional>
#include <vector>

namespace vela {

namespace detail {
template<class Scalar> struct BorderedLoadLineSystem;
}

struct ExtendedBorderedDirection {
    ExtendedVector update;
    Eigen::Index factoredDimension;
    long double backwardErrorBeforeRestoration;
    long double backwardError;
    long double scaledResidualInfinityNorm;
};

enum class ExtendedLinearDefectMethod { mp100, eft };

struct ExtendedLinearAccuracyGoal {
    long double backwardError, scaledDefect;
};

/// Factor the scaled system, refine against the original matrix and recover
/// the physical coordinates without narrowing. Independent retained rows
/// remain in the factorization and are restored before auditing the result.
ExtendedBorderedDirection solveExtendedBorderedSystem(
    const detail::BorderedLoadLineSystem<long double>& system,
    const std::vector<int>& eliminatedRows,
    std::span<const int> retainedIndependentRows,
    int refinementCorrections,
    std::optional<ExtendedLinearAccuracyGoal> goal = {},
    ExtendedLinearDefectMethod defectMethod = ExtendedLinearDefectMethod::mp100);

} // namespace vela
