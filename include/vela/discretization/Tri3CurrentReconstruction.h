#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace vela::detail {

template<class Scalar>
std::array<Scalar, 2> tri3CurrentFromEdgePair(const std::array<double, 2>& a,
    const Scalar& ja, const std::array<double, 2>& b, const Scalar& jb)
{
    const double determinant = a[0] * b[1] - a[1] * b[0];
    if (std::abs(determinant) <= 1e-14)
        throw std::invalid_argument("Degenerate triangle current projection");
    return {(ja * Scalar(b[1]) - Scalar(a[1]) * jb) / Scalar(determinant),
            (Scalar(a[0]) * jb - ja * Scalar(b[0])) / Scalar(determinant)};
}

// Edges are oriented 0->1, 1->2, 2->0. Each edge vector averages the
// two pair reconstructions using the other edges' dual-face measures.
template<class Scalar>
std::array<Scalar, 2> gssLauxCurrent(const std::array<std::array<double, 2>, 3>& tangents,
    const std::array<double, 3>& lengths, const std::array<double, 3>& partialVolumes,
    const std::array<Scalar, 3>& current)
{
    const std::array<std::array<Scalar, 2>, 3> pairs{
        tri3CurrentFromEdgePair(tangents[0], current[0], tangents[1], current[1]),
        tri3CurrentFromEdgePair(tangents[0], current[0], tangents[2], current[2]),
        tri3CurrentFromEdgePair(tangents[1], current[1], tangents[2], current[2])};
    const auto pairIndex = [](int a, int b) {
        return std::min(a, b) == 0 ? std::max(a, b) - 1 : 2;
    };
    const double total = partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
    if (total <= 1e-300)
        throw std::invalid_argument("Triangle has no positive box partial volume");
    std::array<std::array<Scalar, 2>, 3> edgeVectors;
    for (int edge = 0; edge < 3; ++edge) {
        const int a = edge == 0 ? 1 : 0;
        const int b = edge == 2 ? 1 : 2;
        const auto& first = pairs[pairIndex(edge, a)];
        const auto& second = pairs[pairIndex(edge, b)];
        const double wa = 2.0 * partialVolumes[a] / lengths[a];
        const double wb = 2.0 * partialVolumes[b] / lengths[b];
        const double sum = wa + wb;
        for (int c = 0; c < 2; ++c)
            edgeVectors[edge][c] = sum > 0.0
                ? (Scalar(wa) * first[c] + Scalar(wb) * second[c]) / Scalar(sum)
                : Scalar(0.5) * (first[c] + second[c]);
    }
    return {(Scalar(partialVolumes[0]) * edgeVectors[0][0] + Scalar(partialVolumes[1]) * edgeVectors[1][0] +
             Scalar(partialVolumes[2]) * edgeVectors[2][0]) / Scalar(total),
            (Scalar(partialVolumes[0]) * edgeVectors[0][1] + Scalar(partialVolumes[1]) * edgeVectors[1][1] +
             Scalar(partialVolumes[2]) * edgeVectors[2][1]) / Scalar(total)};
}

} // namespace vela::detail
