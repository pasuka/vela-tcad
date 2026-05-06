#include "vela/numerics/LineSearch.h"
#include <algorithm>
#include <cmath>

namespace vela {

BacktrackingLineSearch::BacktrackingLineSearch(LineSearchConfig cfg)
    : cfg_(cfg)
{}

LineSearchResult BacktrackingLineSearch::search(
    const VectorXd& x,
    const VectorXd& step,
    const VectorXd& currentResidual,
    const ResidualFunction& residualFunction,
    const AcceptFunction& acceptFunction) const
{
    const Real currentNorm = currentResidual.norm();
    // When line search is disabled the caller's initialDamping is applied as a
    // fixed damping factor; minDamping applies only during backtracking so it
    // must not constrain a disabled line search from below.
    Real alpha = cfg_.enabled
        ? std::clamp(cfg_.initialDamping, cfg_.minDamping, 1.0)
        : std::clamp(cfg_.initialDamping, 0.0, 1.0);

    const int attempts = cfg_.enabled ? std::max(1, cfg_.maxBacktracks + 1) : 1;
    for (int k = 0; k < attempts; ++k) {
        VectorXd candidate = x + alpha * step;
        VectorXd residual = residualFunction(candidate);
        const Real norm = residual.norm();
        const bool finite = candidate.allFinite() && residual.allFinite() && std::isfinite(norm);
        const bool acceptedByCaller = !acceptFunction || acceptFunction(candidate, residual);
        const Real target = (1.0 - cfg_.sufficientDecrease * alpha) * currentNorm;

        if (finite && acceptedByCaller && (!cfg_.enabled || norm <= target || norm < currentNorm)) {
            return {candidate, residual, alpha, norm, true};
        }

        alpha *= cfg_.reduction;
        if (alpha < cfg_.minDamping)
            break;
    }

    return {x, currentResidual, 0.0, currentNorm, false};
}

} // namespace vela
