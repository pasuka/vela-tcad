#include "vela/numerics/LineSearch.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace vela {

namespace {

std::string lineSearchRejectionReason(bool finite,
                                      bool acceptedByCaller,
                                      bool sufficientDecrease)
{
    if (!finite)
        return "nonfinite_residual";
    if (!acceptedByCaller)
        return "carrier_invalid";
    if (!sufficientDecrease)
        return "line_search_non_decrease";
    return "unknown_rejection";
}

} // namespace

BacktrackingLineSearch::BacktrackingLineSearch(LineSearchConfig cfg)
    : cfg_(cfg)
{}

LineSearchResult BacktrackingLineSearch::search(
    const VectorXd& x,
    const VectorXd& step,
    const VectorXd& currentResidual,
    const ResidualFunction& residualFunction,
    const AcceptFunction& acceptFunction,
    const NormFunction& normFunction) const
{
    const auto normOf = [&](const VectorXd& residual) {
        return normFunction ? normFunction(residual) : residual.norm();
    };
    const Real currentNorm = normOf(currentResidual);
    // When line search is disabled the caller's initialDamping is applied as a
    // fixed damping factor; minDamping applies only during backtracking so it
    // must not constrain a disabled line search from below.
    Real alpha = cfg_.enabled
        ? std::clamp(cfg_.initialDamping, cfg_.minDamping, 1.0)
        : std::clamp(cfg_.initialDamping, 0.0, 1.0);

    std::vector<LineSearchIterationInfo> history;
    if (cfg_.recordHistory) {
        const int reserve = cfg_.enabled ? std::max(1, cfg_.maxBacktracks + 1) : 1;
        history.reserve(static_cast<std::size_t>(reserve));
    }

    int attemptCount = 0;
    const int attempts = cfg_.enabled ? std::max(1, cfg_.maxBacktracks + 1) : 1;
    std::string lastRejectionReason;
    for (int k = 0; k < attempts; ++k) {
        VectorXd candidate = x + alpha * step;
        VectorXd residual = residualFunction(candidate);
        const Real norm = normOf(residual);
        const bool finite = candidate.allFinite() && residual.allFinite() && std::isfinite(norm);
        const bool acceptedByCaller = !acceptFunction || acceptFunction(candidate, residual);
        const Real target = (1.0 - cfg_.sufficientDecrease * alpha) * currentNorm;
        const bool sufficientDecrease = !cfg_.enabled || norm <= target || norm < currentNorm;
        const bool accepted = finite && acceptedByCaller && sufficientDecrease;
        const std::string rejectionReason = accepted
            ? std::string{}
            : lineSearchRejectionReason(finite, acceptedByCaller, sufficientDecrease);
        ++attemptCount;

        if (cfg_.recordHistory) {
            history.push_back({
                k,
                alpha,
                norm,
                target,
                finite,
                acceptedByCaller,
                sufficientDecrease,
                accepted,
                rejectionReason});
        }

        if (accepted) {
            return {candidate, residual, alpha, norm, true, attemptCount, {}, std::move(history)};
        }

        lastRejectionReason = rejectionReason;

        alpha *= cfg_.reduction;
        if (alpha < cfg_.minDamping)
            break;
    }

    if (lastRejectionReason.empty())
        lastRejectionReason = "line_search_rejected";
    return {x, currentResidual, 0.0, currentNorm, false, attemptCount,
            lastRejectionReason, std::move(history)};
}

} // namespace vela
