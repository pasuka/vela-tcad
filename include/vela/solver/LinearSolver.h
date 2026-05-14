#pragma once

#include "vela/core/Types.h"

#include <Eigen/SparseLU>
#include <cstddef>
#include <vector>

namespace vela {

/**
 * @brief Sparse direct linear solver based on Eigen SparseLU.
 *
 * Wraps Eigen::SparseLU for the convenience of the rest of the solver
 * pipeline.  Swap this class for an iterative solver (e.g. BiCGSTAB)
 * without changing any call sites.
 */
class LinearSolver {
public:
    /**
     * @brief Solve the linear system A * x = b.
     *
     * Reuses Eigen's symbolic analysis when consecutive solves have the same
     * sparse structure. Numerical factorisation is still performed every call
     * so changed coefficient values are reflected in the solution.
     *
     * @throws std::invalid_argument if dimensions are inconsistent.
     * @throws std::runtime_error if the analysis, factorisation, or solve fails.
     */
    VectorXd solve(const SparseMatrixd& A, const VectorXd& b);

    /**
     * @brief Clear the cached sparse pattern and force re-analysis next solve.
     */
    void clearPatternCache();

    /**
     * @brief Number of symbolic pattern analyses performed by this instance.
     *
     * This is primarily useful for tests and lightweight performance
     * diagnostics; it does not count numerical factorisations.
     */
    std::size_t patternAnalysisCount() const noexcept;

private:
    using StorageIndex = SparseMatrixd::StorageIndex;

    bool patternMatches(const SparseMatrixd& A) const;
    void cachePattern(const SparseMatrixd& A);
    void analyzePatternIfNeeded(const SparseMatrixd& A);

    Eigen::SparseLU<SparseMatrixd> solver_;
    bool hasAnalyzedPattern_ = false;
    int cachedRows_ = 0;
    int cachedCols_ = 0;
    std::size_t cachedNonZeros_ = 0;
    std::vector<StorageIndex> cachedOuterStarts_;
    std::vector<StorageIndex> cachedInnerIndices_;
    std::size_t patternAnalysisCount_ = 0;
};

} // namespace vela
