#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <unordered_set>
#include <vector>

namespace vela::experimental {

// Arithmetic-only report primitive. Not a Newton acceptance predicate. Callers
// must establish variable identity, units, active masks and update provenance.
enum class UpdateErrorReportMode { Off, Report };
enum class UpdateErrorEquation : std::size_t { Poisson, Electron, Hole, Contact };
enum class UpdateErrorValidity { Valid, CandidateUnrepresentable, UnderflowDegeneracy, NoAcceptedUpdate };
enum class UpdateErrorDeltaKind { Unspecified, RawLinearTangent, RawFiniteDifference, AppliedDifference };
enum class UpdateErrorDenominatorKind { Unspecified, Before, RawTangentAfter, RawFiniteAfter, AppliedAfter };

// The numerator and denominator are independent choices. In particular, a raw
// linear direction with an applied-state denominator is not an applied delta.
struct UpdateErrorProvenance {
    UpdateErrorDeltaKind delta = UpdateErrorDeltaKind::Unspecified;
    UpdateErrorDenominatorKind denominator = UpdateErrorDenominatorKind::Unspecified;
};

struct NewtonUpdateErrorReportConfig {
    UpdateErrorReportMode mode = UpdateErrorReportMode::Off;
    int digits = 5;
    // Explicit metric units: V, cm^-3, cm^-3, V for a voltage contact variable.
    // The caller must supply the actual reference for its contact formulation.
    std::array<double, 4> errRef{0.025852, 1e10, 1e10, 0.025852};
};

struct UpdateErrorRow {
    UpdateErrorEquation equation = UpdateErrorEquation::Poisson;
    std::size_t id = 0;
    bool active = true;
    double denominatorState = 0.;
    double delta = 0.;
    UpdateErrorValidity validity = UpdateErrorValidity::Valid;
    // Identity within this equation's independent solve coordinates, if known.
    // Several distinct metric entries may alias one unknown. Never deduplicate
    // their RMS contributions. Missing identity leaves the solve count unknown.
    std::optional<std::size_t> independentUnknownId;
};

struct UpdateErrorRowContribution {
    UpdateErrorRow row;
    std::optional<double> normalizedMagnitude;
};

struct UpdateErrorEquationReport {
    std::size_t count = 0;
    std::size_t invalidCount = 0;
    std::optional<double> rms;
    std::optional<std::size_t> largestContributor;
    std::optional<std::size_t> independentUnknownCount;
};

struct NewtonUpdateErrorReport {
    bool enabled = false;
    bool invalidConfig = false;
    bool completeEquationSet = false;
    UpdateErrorProvenance provenance;
    std::size_t metricEntryCount = 0;
    std::optional<std::size_t> independentUnknownCount;
    std::array<UpdateErrorEquationReport, 4> equations;
    std::optional<double> rms;
    std::vector<UpdateErrorRowContribution> rows;
};

namespace update_error_detail {
// Scaled sum of squares avoids squaring extreme finite contributions.
struct RmsAccumulator {
    double scale = 0., ssq = 0.;
    void add(double value) {
        if (value == 0.) return;
        if (scale < value) {
            const double ratio = scale / value;
            ssq = 1. + ssq * ratio * ratio;
            scale = value;
        } else {
            const double ratio = value / scale;
            ssq += ratio * ratio;
        }
    }
    double rms(std::size_t count) const {
        return scale * std::sqrt(ssq / static_cast<double>(count));
    }
};
}

// No callbacks or mutable solver objects: this cannot assemble or accept a step.
// Evaluate one (variable identity, update category, denominator time) combination.
inline NewtonUpdateErrorReport evaluateUpdateErrorReport(
    const NewtonUpdateErrorReportConfig& config,
    std::span<const UpdateErrorRow> rows,
    bool completeEquationSet = false,
    UpdateErrorProvenance provenance = {})
{
    NewtonUpdateErrorReport result;
    if (config.mode == UpdateErrorReportMode::Off) return result;
    result.enabled = true;
    result.completeEquationSet = completeEquationSet;
    result.provenance = provenance;
    const double amplification = std::pow(10., config.digits);
    if (!std::isfinite(amplification) || amplification == 0.) result.invalidConfig = true;
    for (double reference : config.errRef)
        if (!std::isfinite(reference) || reference <= 0.) result.invalidConfig = true;
    if (result.invalidConfig) return result;

    std::array<update_error_detail::RmsAccumulator, 4> sums;
    std::array<double, 4> largest{-1., -1., -1., -1.};
    std::array<std::unordered_set<std::size_t>, 4> unknowns;
    std::array<bool, 4> allIdentified{true, true, true, true};
    update_error_detail::RmsAccumulator total;
    std::size_t count = 0, invalid = 0;
    for (const auto& row : rows) {
        const auto k = static_cast<std::size_t>(row.equation);
        if (k >= result.equations.size()) {
            result.invalidConfig = true; result.rms.reset(); return result;
        }
        UpdateErrorRowContribution item{row, std::nullopt};
        if (row.active) {
            auto& equation = result.equations[k];
            ++equation.count; ++count;
            if (row.independentUnknownId) unknowns[k].insert(*row.independentUnknownId);
            else allIdentified[k] = false;
            bool valid = row.validity == UpdateErrorValidity::Valid &&
                std::isfinite(row.denominatorState) && std::isfinite(row.delta);
            double value = 0.;
            if (valid) {
                const double a = std::abs(row.denominatorState), b = config.errRef[k];
                // Scale all three operands to avoid overflowing abs(x)+reference.
                const double scale = std::max({a, b, std::abs(row.delta)});
                const double ratio = (std::abs(row.delta) / scale) / (a / scale + b / scale);
                value = ratio * amplification;
                valid = std::isfinite(value);
                if (valid && value == 0. && row.delta != 0.) {
                    valid = false;
                    item.row.validity = UpdateErrorValidity::UnderflowDegeneracy;
                }
            }
            if (!valid) {
                if (item.row.validity == UpdateErrorValidity::Valid)
                    item.row.validity = UpdateErrorValidity::CandidateUnrepresentable;
                ++equation.invalidCount; ++invalid;
            }
            else {
                item.normalizedMagnitude = value;
                sums[k].add(value); total.add(value);
                if (value > largest[k]) { largest[k] = value; equation.largestContributor = row.id; }
            }
        }
        result.rows.push_back(item);
    }
    result.metricEntryCount = count;
    std::size_t independentCount = 0;
    bool completeIdentity = true;
    for (std::size_t k = 0; k < result.equations.size(); ++k) {
        auto& equation = result.equations[k];
        if (allIdentified[k]) {
            equation.independentUnknownCount = unknowns[k].size();
            independentCount += unknowns[k].size();
        } else completeIdentity = false;
        if (equation.count && !equation.invalidCount) equation.rms = sums[k].rms(equation.count);
    }
    if (count && completeIdentity) result.independentUnknownCount = independentCount;
    if (count && !invalid) result.rms = total.rms(count);
    return result;
}

} // namespace vela::experimental
