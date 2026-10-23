#pragma once

#include <cmath>
#include <stdexcept>

namespace vela::experimental {

// Convert BEFORE multiplying or adding. A caller-supplied wide Scalar retains
// the original binary64 direction even when the update is below the reference
// value's ULP. This does not redefine the direction from rounded end states.
template<class Scalar>
Scalar referencedDirectionValue(double reference, double step, double direction) {
    if (!std::isfinite(reference) || !std::isfinite(step) || !std::isfinite(direction))
        throw std::invalid_argument("Nonfinite referenced direction input");
    return Scalar(reference) + Scalar(step) * Scalar(direction);
}

// An affine physical boundary voltage converted to packed coordinates. The
// frozen target is already packed; only the physical voltage increment is scaled.
template<class Scalar>
Scalar referencedBoundaryValue(double packedReference, double step,
                               double voltageDirection, double potentialScale) {
    if (!std::isfinite(packedReference) || !std::isfinite(potentialScale) || potentialScale <= 0.0)
        throw std::invalid_argument("Invalid boundary potential scale");
    return Scalar(packedReference) +
        referencedDirectionValue<Scalar>(0.0, step, voltageDirection) / Scalar(potentialScale);
}

} // namespace vela::experimental
