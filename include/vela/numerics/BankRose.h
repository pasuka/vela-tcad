#pragma once

#include "vela/numerics/LineSearch.h"

namespace vela {

/// Bank--Rose Algorithm Global (1981), with bounded failure handling.
struct BankRoseConfig {
    bool enabled = false;
    int notDamped = 25; ///< One-based Newton iterations <= this take a full step.
    Real sufficientDecrease = 1.0e-4;
    Real initialK = 0.0;
    Real seedK = 1.0; ///< First positive K after rejecting K=0; inverse norm units.
    Real growth = 10.0;
    Real decay = 0.1;
    Real minDamping = 1.0e-8;
    int maxTrials = 32;
};

void validateBankRoseConfig(const BankRoseConfig& config);

/// State is local to one nonlinear solve. A failed search does not commit K.
/// Caller supplies the actual constrained direction and a fixed residual norm.
class BankRoseDamping {
public:
    explicit BankRoseDamping(BankRoseConfig config = {});
    Real nextK() const { return nextK_; }
    LineSearchResult search(
        int iteration, const VectorXd& x, const VectorXd& step,
        const VectorXd& currentResidual,
        const BacktrackingLineSearch::ResidualFunction& residualFunction,
        const BacktrackingLineSearch::AcceptFunction& acceptFunction = {},
        const BacktrackingLineSearch::NormFunction& normFunction = {},
        const BacktrackingLineSearch::CandidateFunction& candidateFunction = {},
        bool recordHistory = false);

private:
    BankRoseConfig cfg_;
    Real nextK_;
};

} // namespace vela
