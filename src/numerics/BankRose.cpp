#include "vela/numerics/BankRose.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/PopulationDomainError.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace vela {

void validateBankRoseConfig(const BankRoseConfig& c)
{
    if (c.notDamped < 0 || c.maxTrials < 1 || c.maxTrials > 256 ||
        !std::isfinite(c.sufficientDecrease) || c.sufficientDecrease <= 0 ||
        c.sufficientDecrease >= 1 || !std::isfinite(c.initialK) || c.initialK < 0 ||
        !std::isfinite(c.seedK) || c.seedK <= 0 ||
        !std::isfinite(c.growth) || c.growth <= 1 ||
        !std::isfinite(c.decay) || c.decay <= 0 || c.decay >= 1 ||
        !std::isfinite(c.minDamping) || c.minDamping <= 0 || c.minDamping > 1) {
        throw std::invalid_argument("Invalid bank_rose parameters: require not_damped>=0, "
            "max_trials in [1,256], omega/decay in (0,1), finite K>=0, seed_K>0, "
            "growth>1 and min_damping in (0,1].");
    }
}

BankRoseDamping::BankRoseDamping(BankRoseConfig config)
    : cfg_(config), nextK_(config.initialK)
{
    validateBankRoseConfig(cfg_);
}

LineSearchResult BankRoseDamping::search(
    int iteration, const VectorXd& x, const VectorXd& step,
    const VectorXd& currentResidual,
    const BacktrackingLineSearch::ResidualFunction& residualFunction,
    const BacktrackingLineSearch::AcceptFunction& acceptFunction,
    const BacktrackingLineSearch::NormFunction& normFunction,
    const BacktrackingLineSearch::CandidateFunction& candidateFunction,
    bool recordHistory)
{
    if (iteration < 1 || x.size() != step.size() || !residualFunction)
        throw std::invalid_argument("BankRoseDamping requires a positive iteration, compatible direction and residual callback.");
    ScopedPerformanceTimer timer("newton.bank_rose");
    incrementPerformanceCounter("newton.bank_rose_calls");
    const auto normOf = [&](const VectorXd& r) { return normFunction ? normFunction(r) : r.stableNorm(); };
    const Real currentNorm = normOf(currentResidual);
    LineSearchResult out{x, currentResidual, 0.0, currentNorm, false, 0, {}, {}};
    if (!x.allFinite() || !step.allFinite() || !currentResidual.allFinite() ||
        !std::isfinite(currentNorm) || currentNorm < 0) {
        out.failureReason = "bank_rose_nonfinite_input";
        return out;
    }
    // Zero residual is no permission to move the state or declare device qualification.
    if (currentNorm == 0) {
        out.failureReason = "bank_rose_zero_residual";
        return out;
    }
    const bool active = iteration > cfg_.notDamped;
    Real k = active ? nextK_ : 0.0;
    Real bestNorm = std::numeric_limits<Real>::infinity();
    for (int trial = 0; trial < (active ? cfg_.maxTrials : 1); ++trial) {
        // Extended intermediate avoids avoidable overflow in K*norm; storage remains Real.
        const Real alpha = active ? static_cast<Real>(1.0L /
            (1.0L + static_cast<long double>(k) * currentNorm)) : 1.0;
        if (!std::isfinite(alpha) || alpha < cfg_.minDamping || alpha <= 0) {
            out.failureReason = "bank_rose_min_damping";
            break;
        }
        VectorXd candidate;
        VectorXd r;
        Real norm = std::numeric_limits<Real>::infinity();
        bool populationDomainError = false;
        try {
            candidate = candidateFunction ? candidateFunction(x, step, alpha) : VectorXd(x + alpha * step);
            if (candidate.allFinite()) r = residualFunction(candidate);
        } catch (const PopulationDomainError&) {
            // No residual exists for this trial. Use the same bounded K retry
            // as other invalid trials; never alter the parent or clamp density.
            populationDomainError = true;
            incrementPerformanceCounter("newton.bank_rose_population_domain_rejections");
        }
        if (!populationDomainError && r.size() == currentResidual.size()) norm = normOf(r);
        const bool finite = !populationDomainError && candidate.allFinite() && r.size() == currentResidual.size() &&
            r.allFinite() && std::isfinite(norm) && norm >= 0;
        const bool valid = finite && (!acceptFunction || acceptFunction(candidate, r));
        const long double target = (1.0L - static_cast<long double>(cfg_.sufficientDecrease) * alpha) * currentNorm;
        // Strict norm reduction also guards a rounded target equal to currentNorm.
        const bool decrease = !active || (norm < currentNorm && static_cast<long double>(norm) <= target);
        const bool accepted = valid && decrease;
        std::string reason = accepted ? "" : populationDomainError ? "bank_rose_population_domain" : !finite ? "bank_rose_nonfinite_trial" :
            !valid ? "bank_rose_invalid_trial" : "bank_rose_insufficient_decrease";
        ++out.attempts;
        incrementPerformanceCounter("newton.bank_rose_trials");
        if (recordHistory) {
            LineSearchIterationInfo info{trial, alpha, norm,
                active ? static_cast<Real>(target) : currentNorm, finite, valid, decrease, accepted, reason};
            info.bankRoseActive = active;
            info.bankRoseK = k;
            out.history.push_back(std::move(info));
        }
        if (accepted) {
            if (active) nextK_ = k * cfg_.decay;
            out.x = std::move(candidate);
            out.residual = std::move(r);
            out.damping = alpha;
            out.residualNorm = norm;
            out.accepted = true;
            out.failureReason.clear();
            incrementPerformanceCounter("newton.bank_rose_accepted");
            return out;
        }
        out.failureReason = reason;
        if (valid && norm < bestNorm) {
            bestNorm = norm;
            out.bestRejectedCandidate = true;
            out.bestRejectedX = candidate;
            out.bestRejectedResidual = r;
            out.bestRejectedDamping = alpha;
            out.bestRejectedResidualNorm = norm;
        }
        if (!active || trial + 1 == cfg_.maxTrials) break;
        if (k == 0) k = cfg_.seedK;
        else if (k > std::numeric_limits<Real>::max() / cfg_.growth) {
            out.failureReason = "bank_rose_k_overflow";
            break;
        } else k *= cfg_.growth;
    }
    incrementPerformanceCounter("newton.bank_rose_rejected");
    return out;
}
} // namespace vela
