#pragma once

#include <cmath>

namespace vela::detail {

// B(-a)c0 - B(a)c1 for a carrier-oriented Bernoulli argument.
template<class Scalar>
Scalar sgWeightedPopulationDifference(const Scalar& coefficient,
    const Scalar& from, const Scalar& to,
    const Scalar& backwardWeight, const Scalar& forwardWeight)
{
    return coefficient * (from * backwardWeight - to * forwardWeight);
}

// Evaluate the left/right difference through expm1 when the caller has a
// resolved log ratio. Keeping the weighted right term explicit lets the
// assembler choose how to form it in its working units and precision.
template<class Scalar, class ExpMinusOne>
Scalar sgFactorizedPopulationDifference(const Scalar& weightedRight,
    const Scalar& logLeftOverRight, ExpMinusOne expMinusOne)
{
    return weightedRight * expMinusOne(logLeftOverRight);
}

// For a = log(c1/c0) + delta the two forms are algebraically identical.
// With a large drive they cannot nearly cancel and subtracting the products
// avoids expm1 overflow and preserves the smaller endpoint's tangent.
template<class Scalar, class Bernoulli, class ExpMinusOne>
Scalar sgBalancedEdgeFlux(const Scalar& coefficient,
    const Scalar& from, const Scalar& to, const Scalar& delta,
    double deltaValue, const Scalar& argument, Bernoulli bernoulli,
    ExpMinusOne expMinusOne)
{
    if (std::abs(deltaValue) > 50.0)
        return sgWeightedPopulationDifference(
            coefficient, from, to, bernoulli(-argument), bernoulli(argument));
    return sgFactorizedPopulationDifference(
        coefficient * to * bernoulli(argument), delta, expMinusOne);
}

} // namespace vela::detail
