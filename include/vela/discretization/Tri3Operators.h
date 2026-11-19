#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>

namespace vela::detail {

// The stencil interleaves x and y coefficients for vertices 1 and 2. Forming
// differences relative to vertex 0 also gives a zero gradient for constants.
template<class Scalar>
std::array<Scalar, 2> tri3Gradient(const std::array<Scalar, 3>& values,
    const std::array<double, 4>& stencil)
{
    const auto first = values[1] - values[0];
    const auto second = values[2] - values[0];
    return {first * Scalar(stencil[0]) + second * Scalar(stencil[2]),
            first * Scalar(stencil[1]) + second * Scalar(stencil[3])};
}

// The input measures can represent either geometric vertex boxes or
// coefficient-derived boxes. The caller selects the mesh convention.
template<class Scalar, class Measure = double, class VertexValue>
Scalar tri3VertexAverage(const std::array<double, 3>& measures,
    VertexValue vertexValue)
{
    Scalar total{};
    Measure measure{};
    for (std::size_t vertex = 0; vertex < measures.size(); ++vertex) {
        total = total + Scalar(measures[vertex]) * Scalar(vertexValue(vertex));
        measure += measures[vertex];
    }
    if (!(measure > 0.0))
        throw std::invalid_argument("Tri3 vertex average requires positive total measure");
    return total / Scalar(measure);
}

} // namespace vela::detail
