#pragma once

#include "vela/core/Types.h"
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <stdexcept>

namespace vela {

// Opt-in stopping experiment. Concentration references are always physical
// cm^-3, independently of the device's internal unit/scaling choices.
struct DensityUpdateConvergenceConfig {
    std::string mode = "off"; // off, report, enforce (OR), enforce_joint (AND)
    int digits = 5;
    Real potentialReference_V = 0.025852;
    Real electronReference_cm3 = 1e10;
    Real holeReference_cm3 = 1e10;
    Real contactReference_V = 0.025852;
    std::string diagnosticCsv;
};

inline void validateDensityUpdateConvergence(const DensityUpdateConvergenceConfig& c)
{
    if (c.mode != "off" && c.mode != "report" && c.mode != "enforce" && c.mode != "enforce_joint")
        throw std::invalid_argument("density update convergence mode must be off, report, enforce or enforce_joint");
    if (c.digits < 0 || c.digits > 15)
        throw std::invalid_argument("density update convergence digits must be in [0,15]");
    for (Real reference : {c.potentialReference_V, c.electronReference_cm3,
                           c.holeReference_cm3, c.contactReference_V})
        if (!(reference > 0.) || !std::isfinite(reference))
            throw std::invalid_argument("density update convergence references must be positive and finite");
}

template<class Scalar>
struct BasicDensityUpdateError {
    bool valid = false;
    Scalar rms = std::numeric_limits<Scalar>::infinity();
    // Poisson, electron, hole, active-contact aliases (not independent DOFs).
    std::array<std::size_t, 4> counts{};
    std::array<Scalar, 4> blockRms{};
};

using DensityUpdateError = BasicDensityUpdateError<Real>;

template<class Scalar, class ReportScalar = Scalar>
class DensityUpdateAccumulator {
public:
    explicit DensityUpdateAccumulator(int digits) : amplification_(std::pow(Scalar{10},digits)) {}

    void add(int block, Scalar state, Scalar direction, Scalar reference) {
        const Scalar normalized = amplification_ * direction / (std::abs(state) + reference);
        valid_ = valid_ && std::isfinite(state) && std::isfinite(direction) &&
                std::isfinite(normalized);
        squares_[block] += normalized * normalized;
        ++counts_[block];
    }

    BasicDensityUpdateError<ReportScalar> finish(bool validState = true) const {
        BasicDensityUpdateError<ReportScalar> out;
        out.counts = counts_;
        Scalar total = 0;
        std::size_t count = 0;
        for (std::size_t block = 0; block < counts_.size(); ++block) {
            total += squares_[block];
            count += counts_[block];
            if (counts_[block])
                out.blockRms[block] = static_cast<ReportScalar>(std::sqrt(squares_[block]/counts_[block]));
        }
        if (count)
            out.rms = static_cast<ReportScalar>(std::sqrt(total/count));
        out.valid = validState && valid_ && std::isfinite(out.rms);
        return out;
    }
private:
    Scalar amplification_;
    std::array<Scalar,4> squares_{};
    std::array<std::size_t,4> counts_{};
    bool valid_ = true;
};

// Raw, uncapped linear tangent in the numerator; the actual accepted state
// in the denominator. Backtracking must never shrink the numerator into an
// artificial convergence signal. No accepted update => no update criterion.
using DensityUpdateErrorEvaluator = std::function<DensityUpdateError(
    const VectorXd& before, const VectorXd& rawDirection,
    const VectorXd& acceptedAfter, Real voltageAfter, Real rawVoltageDirection)>;

} // namespace vela
