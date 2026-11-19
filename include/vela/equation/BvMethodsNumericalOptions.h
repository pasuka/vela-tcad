#pragma once

namespace vela {

enum class BvMethodsAvalancheWeights { TruncatedBox, CoefficientBox };
enum class BvMethodsAvalancheCurrent { RepresentativeVertex, GssLaux };
enum class BvMethodsAvalancheSwitch { Smooth, Sharp };

struct BvMethodsNumericalOptions {
    BvMethodsAvalancheWeights avalancheWeights = BvMethodsAvalancheWeights::TruncatedBox;
    BvMethodsAvalancheCurrent avalancheCurrent = BvMethodsAvalancheCurrent::RepresentativeVertex;
    BvMethodsAvalancheSwitch avalancheSwitch = BvMethodsAvalancheSwitch::Smooth;

    bool operator==(const BvMethodsNumericalOptions&) const = default;
};

} // namespace vela
