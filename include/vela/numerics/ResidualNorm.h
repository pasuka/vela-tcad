#pragma once

#include "vela/core/Types.h"

namespace vela {

struct ResidualNormValue {
    Real l2 = 0.0;
    Real linf = 0.0;
};

class ResidualNorm {
public:
    static ResidualNormValue compute(const VectorXd& r);
    static Real relative(Real current, Real initial);
};

} // namespace vela
