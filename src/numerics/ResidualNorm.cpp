#include "vela/numerics/ResidualNorm.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela {

ResidualNormValue ResidualNorm::compute(const VectorXd& r)
{
    ResidualNormValue value;
    value.l2 = r.norm();
    for (int i = 0; i < r.size(); ++i)
        value.linf = std::max(value.linf, std::abs(r(i)));
    return value;
}

ResidualBlockNormValue ResidualNorm::computeBlocks(const VectorXd& r, Index numNodes)
{
    const int N = static_cast<int>(numNodes);
    if (N < 0 || r.size() != 3 * N)
        throw std::invalid_argument(
            "ResidualNorm::computeBlocks: residual size must equal 3*numNodes.");

    ResidualBlockNormValue value;
    value.psi = r.segment(0, N).norm();
    value.phin = r.segment(N, N).norm();
    value.phip = r.segment(2 * N, N).norm();
    value.combined = std::sqrt(value.psi * value.psi
        + value.phin * value.phin
        + value.phip * value.phip);
    return value;
}

Real ResidualNorm::normalizedBlockL2(const ResidualBlockNormValue& current,
                                     const ResidualBlockNormValue& scale,
                                     const ResidualBlockWeights& weights)
{
    const auto normalizedSquare = [](Real value, Real denom, Real weight) {
        if (weight <= 0.0)
            return 0.0;
        const Real safeDenom = std::max(std::abs(denom), 1.0e-300);
        const Real normalized = value / safeDenom;
        return weight * normalized * normalized;
    };

    return std::sqrt(
        normalizedSquare(current.psi, scale.psi, weights.psi)
        + normalizedSquare(current.phin, scale.phin, weights.phin)
        + normalizedSquare(current.phip, scale.phip, weights.phip));
}

Real ResidualNorm::relative(Real current, Real initial)
{
    return current / std::max(initial, 1.0e-300);
}

} // namespace vela
