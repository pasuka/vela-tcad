#pragma once

#include <algorithm>
#include <cmath>

namespace vela::detail {

template<class Scalar>
Scalar slotboomNarrowing(Scalar concentration, Scalar reference, Scalar coefficient,
                       Scalar smoothing, Scalar offset)
{
    if (concentration <= 0 || coefficient == 0)
        return 0;
    const auto x = std::log(concentration / reference);
    return std::max(offset + coefficient * (x + std::sqrt(x * x + smoothing)), Scalar{0});
}

} // namespace vela::detail
