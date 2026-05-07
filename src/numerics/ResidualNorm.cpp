#include "vela/numerics/ResidualNorm.h"
#include <algorithm>
#include <cmath>

namespace vela {

ResidualNormValue ResidualNorm::compute(const VectorXd& r)
{
    ResidualNormValue value;
    value.l2 = r.norm();
    for (int i = 0; i < r.size(); ++i)
        value.linf = std::max(value.linf, std::abs(r(i)));
    return value;
}

Real ResidualNorm::relative(Real current, Real initial)
{
    return current / std::max(initial, 1.0e-300);
}

} // namespace vela
