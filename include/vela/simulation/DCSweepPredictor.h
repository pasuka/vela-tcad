#pragma once

#include "vela/simulation/DCSweep.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela::detail {

inline void requireMatchingFieldSize(const VectorXd& previous,
                                     const VectorXd& current,
                                     const char*     fieldName)
{
    if (previous.size() != current.size()) {
        throw std::invalid_argument(
            std::string("DCSweep predictor: field size mismatch for ") + fieldName);
    }
}

inline bool predictorUsesField(const std::vector<std::string>& fields,
                               const std::string&              field)
{
    return std::find(fields.begin(), fields.end(), field) != fields.end();
}

inline std::vector<std::string> effectivePredictorFields(const SweepPredictorConfig& config)
{
    if (!config.fields.empty())
        return config.fields;
    if (config.mode == "none" || config.mode == "constant")
        return {};
    return {"psi", "phin", "phip"};
}

inline Real maxPsiPhinJump(const DDSolution& previous, const DDSolution& current)
{
    if (previous.psi.size() != current.psi.size() ||
        previous.phin.size() != current.phin.size() ||
        previous.psi.size() != previous.phin.size() ||
        current.psi.size() != current.phin.size()) {
        throw std::invalid_argument(
            "DCSweep: psi-phin branch acceptance requires matching psi/phin field sizes.");
    }

    Real maxJump = 0.0;
    for (int i = 0; i < current.psi.size(); ++i) {
        const Real previousExponent = previous.psi(i) - previous.phin(i);
        const Real currentExponent = current.psi(i) - current.phin(i);
        const Real jump = std::abs(currentExponent - previousExponent);
        if (std::isfinite(jump))
            maxJump = std::max(maxJump, jump);
        else
            return std::numeric_limits<Real>::infinity();
    }
    return maxJump;
}

struct ElectronDensityJumpStats {
    Real medianDex = 0.0;
    Real p95AbsDex = 0.0;
    Real maxAbsDex = 0.0;
    Real maxSignedDex = 0.0;
    Index maxNode = -1;
};

inline ElectronDensityJumpStats electronDensityJumpStats(const DDSolution& previous,
                                                         const DDSolution& current)
{
    requireMatchingFieldSize(previous.n, current.n, "electron density");
    std::vector<Real> signedJumps;
    signedJumps.reserve(static_cast<std::size_t>(current.n.size()));
    ElectronDensityJumpStats stats;
    const Real eps = std::numeric_limits<Real>::min();
    for (int i = 0; i < current.n.size(); ++i) {
        const Real before = std::max(std::abs(previous.n(i)), eps);
        const Real after = std::max(std::abs(current.n(i)), eps);
        const Real jump = std::log10(after) - std::log10(before);
        if (!std::isfinite(jump)) {
            stats.medianDex = std::numeric_limits<Real>::infinity();
            stats.p95AbsDex = std::numeric_limits<Real>::infinity();
            stats.maxAbsDex = std::numeric_limits<Real>::infinity();
            stats.maxSignedDex = jump;
            stats.maxNode = i;
            return stats;
        }
        signedJumps.push_back(jump);
        const Real absJump = std::abs(jump);
        if (stats.maxNode == std::numeric_limits<Index>::max() ||
            absJump > stats.maxAbsDex) {
            stats.maxAbsDex = absJump;
            stats.maxSignedDex = jump;
            stats.maxNode = i;
        }
    }
    if (signedJumps.empty())
        return stats;

    auto medianValues = signedJumps;
    std::sort(medianValues.begin(), medianValues.end());
    const std::size_t n = medianValues.size();
    stats.medianDex = (n % 2 == 0)
        ? 0.5 * (medianValues[n / 2 - 1] + medianValues[n / 2])
        : medianValues[n / 2];

    std::vector<Real> absValues;
    absValues.reserve(signedJumps.size());
    for (Real value : signedJumps)
        absValues.push_back(std::abs(value));
    std::sort(absValues.begin(), absValues.end());
    const std::size_t p95Index = std::min(
        absValues.size() - 1,
        static_cast<std::size_t>(std::ceil(0.95 * static_cast<Real>(absValues.size()))) - 1);
    stats.p95AbsDex = absValues[p95Index];
    return stats;
}

inline std::string electronDensityJumpAcceptanceFailure(
    const SweepBranchAcceptanceConfig& config,
    const ElectronDensityJumpStats&    stats)
{
    if (!std::isfinite(stats.p95AbsDex) ||
        stats.p95AbsDex > config.maxElectronDensityJumpP95AbsDex) {
        return "electron_density_p95_jump_exceeded";
    }
    if (!std::isfinite(stats.maxAbsDex) ||
        stats.maxAbsDex > config.maxElectronDensityJumpDex) {
        return "electron_density_jump_exceeded";
    }
    return {};
}

inline DDSolution predictDCSweepInitialState(const SweepPredictorConfig& config,
                                             const DDSolution*           previous,
                                             const DDSolution&           current,
                                             Real                       previousBias,
                                             Real                       currentBias,
                                             Real                       targetBias,
                                             int                        retryCount = 0)
{
    DDSolution predicted = current;
    if (config.mode == "none" || config.mode == "constant" ||
        previous == nullptr || retryCount > 0) {
        return predicted;
    }

    const Real denominator = currentBias - previousBias;
    if (denominator == 0.0)
        return predicted;

    const Real unclampedRatio = (targetBias - currentBias) / denominator;
    const Real maxRatio = std::max(config.maxExtrapolationRatio, 1.0);
    const Real ratio = std::clamp(unclampedRatio, -maxRatio, maxRatio);
    const std::vector<std::string> fields = effectivePredictorFields(config);

    if (predictorUsesField(fields, "psi")) {
        requireMatchingFieldSize(previous->psi, current.psi, "psi");
        predicted.psi = current.psi + ratio * (current.psi - previous->psi);
    }
    if (predictorUsesField(fields, "phin")) {
        requireMatchingFieldSize(previous->phin, current.phin, "phin");
        predicted.phin = current.phin + ratio * (current.phin - previous->phin);
    }
    if (predictorUsesField(fields, "phip")) {
        requireMatchingFieldSize(previous->phip, current.phip, "phip");
        predicted.phip = current.phip + ratio * (current.phip - previous->phip);
    }

    return predicted;
}

} // namespace vela::detail
