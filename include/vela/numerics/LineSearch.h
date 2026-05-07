#pragma once

#include "vela/core/Types.h"
#include <functional>

namespace vela {

struct LineSearchConfig {
    bool enabled = true;
    Real initialDamping = 1.0;
    Real minDamping = 1.0e-4;
    Real reduction = 0.5;
    Real sufficientDecrease = 1.0e-4;
    int maxBacktracks = 12;
};

struct LineSearchResult {
    VectorXd x;
    VectorXd residual;
    Real damping = 0.0;
    Real residualNorm = 0.0;
    bool accepted = false;
};

class BacktrackingLineSearch {
public:
    using ResidualFunction = std::function<VectorXd(const VectorXd&)>;
    using AcceptFunction = std::function<bool(const VectorXd&, const VectorXd&)>;

    explicit BacktrackingLineSearch(LineSearchConfig cfg = {});

    LineSearchResult search(const VectorXd& x,
                            const VectorXd& step,
                            const VectorXd& currentResidual,
                            const ResidualFunction& residualFunction,
                            const AcceptFunction& acceptFunction = {}) const;

private:
    LineSearchConfig cfg_;
};

} // namespace vela
