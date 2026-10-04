#pragma once

#include "vela/core/Types.h"
#include <functional>
#include <string>
#include <vector>

namespace vela {

struct LineSearchConfig {
    bool enabled = true;
    Real initialDamping = 1.0;
    Real minDamping = 1.0e-4;
    Real reduction = 0.5;
    Real sufficientDecrease = 1.0e-4;
    int maxBacktracks = 12;
    bool recordHistory = false; ///< Store per-backtrack diagnostics in LineSearchResult::history.
};

struct LineSearchIterationInfo {
    int attempt = 0;
    Real damping = 0.0;
    Real residualNorm = 0.0;
    Real targetResidualNorm = 0.0;
    bool finite = false;
    bool acceptedByCaller = false;
    bool sufficientDecrease = false;
    bool accepted = false;
    std::string rejectionReason;
};

struct LineSearchResult {
    VectorXd x;
    VectorXd residual;
    Real damping = 0.0;
    Real residualNorm = 0.0;
    bool accepted = false;
    int attempts = 0;
    std::string failureReason;
    std::vector<LineSearchIterationInfo> history;
    bool bestRejectedCandidate = false;
    VectorXd bestRejectedX;
    VectorXd bestRejectedResidual;
    Real bestRejectedDamping = 0.0;
    Real bestRejectedResidualNorm = 0.0;
};

class BacktrackingLineSearch {
public:
    using CandidateFunction = std::function<VectorXd(const VectorXd&,const VectorXd&,Real)>;
    using ResidualFunction = std::function<VectorXd(const VectorXd&)>;
    using AcceptFunction = std::function<bool(const VectorXd&, const VectorXd&)>;
    using NormFunction = std::function<Real(const VectorXd&)>;
    using DecreaseAcceptFunction = std::function<bool(const VectorXd&, Real)>;

    explicit BacktrackingLineSearch(LineSearchConfig cfg = {});

    LineSearchResult search(const VectorXd& x,
                            const VectorXd& step,
                            const VectorXd& currentResidual,
                            const ResidualFunction& residualFunction,
                            const AcceptFunction& acceptFunction = {},
                            const NormFunction& normFunction = {},
                            const DecreaseAcceptFunction& decreaseAcceptFunction = {},
                            const CandidateFunction& candidateFunction = {}) const;

private:
    LineSearchConfig cfg_;
};

} // namespace vela
