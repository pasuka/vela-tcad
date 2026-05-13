#pragma once

#include "vela/core/Types.h"

namespace vela {

struct ResidualNormValue {
    Real l2 = 0.0;
    Real linf = 0.0;
};

struct ResidualBlockNormValue {
    Real psi = 0.0;
    Real phin = 0.0;
    Real phip = 0.0;
    Real combined = 0.0;
};

struct ResidualBlockWeights {
    Real psi = 1.0;
    Real phin = 1.0;
    Real phip = 1.0;
};

class ResidualNorm {
public:
    static ResidualNormValue compute(const VectorXd& r);
    static ResidualBlockNormValue computeBlocks(const VectorXd& r, Index numNodes);
    static Real normalizedBlockL2(const ResidualBlockNormValue& current,
                                  const ResidualBlockNormValue& scale,
                                  const ResidualBlockWeights& weights = {});
    static Real relative(Real current, Real initial);
};

} // namespace vela
