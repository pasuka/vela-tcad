#pragma once

#include "vela/core/Types.h"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseLU>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace vela {
namespace detail { class DirectBackend; }

/// Rank-revealing sparse direct fallback for bordered systems whose principal
/// device block may be singular. Returns false instead of throwing when Eigen
/// cannot factor or solve the matrix.
bool solveSparseQrSystem(const SparseMatrixd& A,
                         const VectorXd& b,
                         VectorXd& x) noexcept;

/// Solve with SuiteSparse UMFPACK when Vela was configured with that optional
/// backend. Returns false when the backend is unavailable or the solve fails.
bool solveUmfPackSystem(const SparseMatrixd& A,
                        const VectorXd& b,
                        VectorXd& x) noexcept;

/// Rank-revealing SuiteSparseQR solve for large sparse bordered systems.
bool solveSpqrSystem(const SparseMatrixd& A,
                     const VectorXd& b,
                     VectorXd& x) noexcept;

/**
 * @brief Sparse direct linear solver based on Eigen SparseLU.
 *
 * Wraps Eigen::SparseLU for the convenience of the rest of the solver
 * pipeline.  Swap this class for an iterative solver (e.g. BiCGSTAB)
 * without changing any call sites.
 *
 * Experimental backend override: set the environment variable
 * VELA_LINEAR_SOLVER to one of "sparselu" (default), "sparseqr",
 * "umfpack", "sparselu_metis", "umfpack_metis", "mumps", "mumps_metis",
 * "superlu_mt", "superlu_mt_metis", "strumpack", "bicgstab_ilut",
 * "gmres_ilut", or "simplicial_ldlt". Optional direct backends retain symbolic
 * and identical-value numeric factors; the other alternate backends do not.
 */
class LinearSolver {
public:
    /// Empty selects VELA_LINEAR_SOLVER at construction (default: sparselu).
    explicit LinearSolver(const std::string& backend = "");
    ~LinearSolver();
    /**
     * @brief Solve the linear system A * x = b.
     *
     * Reuses Eigen's symbolic analysis when consecutive solves have the same
     * sparse structure. Numerical factorisation is repeated when coefficient
     * values change; an identical matrix reuses its factors for a new RHS.
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
    struct UmfPackState;
    std::string backend_;
    std::string captureDirectory_;
    /// VELA_LINEAR_FACTOR_STATISTICS=0 disables optional factor diagnostics.
    bool factorStatistics_ = true;
    std::unique_ptr<UmfPackState> umfpack_;
    std::unique_ptr<detail::DirectBackend> direct_;
    void createDirectBackend();
    VectorXd solveDirect(const SparseMatrixd& A, const VectorXd& b);
    VectorXd solveUmfPack(const SparseMatrixd& A, const VectorXd& b);
    using StorageIndex = SparseMatrixd::StorageIndex;
#if defined(VELA_SPARSELU_ORDERING_AMD)
    using SparseLUOrdering = Eigen::AMDOrdering<StorageIndex>;
#else
    using SparseLUOrdering = Eigen::COLAMDOrdering<StorageIndex>;
#endif

    class SparseLUSolver
        : public Eigen::SparseLU<SparseMatrixd, SparseLUOrdering> {
    public:
        bool analysisIsOk() const noexcept { return this->m_analysisIsOk; }
        /// Conventional scalar LU work implied by stored factor structure;
        /// excludes blocked-kernel padding, pivot searches and memory work.
        double structuralFlops() const {
            std::vector<double> lower(this->cols(), 0.), upper(this->cols(), 0.);
            for (Eigen::Index col = 0; col < this->cols(); ++col) {
                for (typename SparseLUSolver::SCMatrix::InnerIterator it(this->m_Lstore, col); it; ++it) {
                    if (it.row() > col) ++lower[col];
                    else if (it.row() < col) ++upper[it.row()];
                }
                for (Eigen::Map<SparseMatrixd>::InnerIterator it(this->m_Ustore, col); it; ++it)
                    if (it.row() < col) ++upper[it.row()];
            }
            double flops = 0.;
            for (Eigen::Index i = 0; i < this->cols(); ++i)
                flops += lower[i] * (1. + 2. * upper[i]);
            return flops;
        }
    };

    bool patternMatches(const SparseMatrixd& A) const;
    void cachePattern(const SparseMatrixd& A);
    void analyzePatternIfNeeded(const SparseMatrixd& A);
    /// True when `A` is bit-for-bit the matrix of the retained factorisation.
    bool valuesMatchFactorization(const SparseMatrixd& A) const;
    void cacheFactorizationValues(const SparseMatrixd& A);

    SparseLUSolver solver_;
    bool hasAnalyzedPattern_ = false;
    Eigen::Index cachedRows_ = 0;
    Eigen::Index cachedCols_ = 0;
    std::size_t cachedNonZeros_ = 0;
    std::vector<StorageIndex> cachedOuterStarts_;
    std::vector<StorageIndex> cachedInnerIndices_;
    /// Values of the matrix whose numeric factorisation `solver_` still holds.
    std::vector<double> cachedValues_;
    bool hasFactorization_ = false;
    std::size_t patternAnalysisCount_ = 0;
};

} // namespace vela
