#include "vela/simulation/ExtendedBorderedLoadLine.h"

#include "vela/simulation/BorderedLoadLineSystem.h"
#include "vela/solver/ConstraintElimination.h"
#include "vela/solver/LinearRefinement.h"
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

namespace vela {

ExtendedBorderedDirection solveExtendedBorderedSystem(
    const detail::BorderedLoadLineSystem<long double>& system,
    const std::vector<int>& eliminatedRows,
    std::span<const int> retainedIndependentRows,
    int refinementCorrections,
    std::optional<ExtendedLinearAccuracyGoal> goal)
{
    const auto size = system.augmented.rows();
    if (size <= 0 || system.augmented.cols() != size ||
        system.rawAugmented.rows() != size || system.rawAugmented.cols() != size ||
        system.rawRhs.size() != size || system.rhs.size() != size ||
        system.rowScale.size() != size || system.columnScale.size() != size ||
        !system.rowScale.allFinite() || !system.columnScale.allFinite() ||
        (system.rowScale.array() <= 0).any() || (system.columnScale.array() <= 0).any() ||
        refinementCorrections < 0)
        throw std::invalid_argument("Invalid extended bordered system");

    const auto reduced = reduceLinearConstraints(system.augmented, system.rhs, eliminatedRows);
    const auto dimension = reduced.matrix.rows();
    std::optional<ExtendedSparseLu> factors;
    if (dimension > 0)
        factors.emplace(reduced.matrix);
    const auto recoverDirection = [&](const BasicConstraintReduction<long double>& reduction) {
        const auto reducedDirection = factors ? factors->solve(reduction.rhs) : ExtendedVector{};
        ExtendedVector unscaled = reduction.recover(reducedDirection).cwiseProduct(system.columnScale);
        return unscaled;
    };

    auto direction = recoverDirection(reduced);
    const ExtendedVector residual = -system.rawRhs;
    std::optional<linear_refinement::PreparedMpDefect<long double>> preparedDefect;
    if (!retainedIndependentRows.empty() || (!goal && refinementCorrections > 0))
        preparedDefect.emplace(system.rawAugmented);
    const auto evaluateDefect = [&] {
        if (preparedDefect)
            return preparedDefect->evaluate(residual, direction);
        return linear_refinement::computeHighPrecisionDefect(system.rawAugmented, residual, direction);
    };
    auto defect = evaluateDefect();
    auto before = defect.backwardMax;
    const auto restore = [&] {
        before = defect.backwardMax;
        if (!retainedIndependentRows.empty()) {
            restoreIndependentConstraintRows(system.rawAugmented, system.rawRhs, retainedIndependentRows, direction);
            defect = evaluateDefect();
        }
    };
    for (int iteration = 0; ; ++iteration) {
        if (goal || iteration == refinementCorrections)
            restore();
        if (iteration == refinementCorrections ||
            (goal && defect.backwardMax <= goal->backwardError &&
             defect.values.cwiseProduct(system.rowScale).cwiseAbs().maxCoeff() <= goal->scaledDefect))
            break;
        const ExtendedVector correctionRhs = -defect.values.cwiseProduct(system.rowScale);
        const auto correction = reduceLinearConstraints(system.augmented, correctionRhs, eliminatedRows);
        direction += recoverDirection(correction);
        defect = evaluateDefect();
    }
    const auto maximum = defect.values.cwiseProduct(system.rowScale).cwiseAbs().maxCoeff();
    if (!std::isfinite(maximum))
        throw std::runtime_error("Nonfinite extended bordered scaled residual");
    return {std::move(direction), dimension, before, defect.backwardMax, maximum};
}

} // namespace vela
