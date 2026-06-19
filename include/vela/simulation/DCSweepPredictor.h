#pragma once

#include "vela/simulation/DCSweep.h"
#include <algorithm>
#include <cmath>
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

inline DDSolution predictDCSweepInitialState(const SweepPredictorConfig& config,
                                             const DDSolution*           previous,
                                             const DDSolution&           current,
                                             Real                       previousBias,
                                             Real                       currentBias,
                                             Real                       targetBias)
{
    DDSolution predicted = current;
    if (config.mode == "none" || config.mode == "constant" || previous == nullptr)
        return predicted;

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
