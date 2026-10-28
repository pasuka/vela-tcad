#pragma once

#include "vela/core/Types.h"
#include "vela/solver/ConstraintElimination.h"

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
 * @brief Sparse linear solver with reusable direct-backend analysis.
 *
 * Experimental backend override: set the environment variable
 * VELA_LINEAR_SOLVER to one of "sparselu", "sparseqr",
 * "umfpack", "sparselu_metis", "umfpack_metis", "mumps", "mumps_metis",
 * "superlu_mt", "superlu_mt_metis", "strumpack", "bicgstab_ilut",
 * "gmres_ilut", "simplicial_ldlt", or "pardiso_mkl". The latter requires an
 * explicit absolute VELA_MKL_RT_LIBRARY (or VELA_MKL_RT_DLL) runtime path on
 * Windows/Linux and is never selected automatically. It uses binary64 LP64.
 * Optional direct backends retain symbolic
 * and identical-value numeric factors; the other alternate backends do not.
 */
class LinearSolver {
public:
    /// Available direct backends in preference order: UMFPACK, SparseLU,
    /// STRUMPACK, MUMPS, SuperLU_MT. SparseLU is always available.
    static std::vector<std::string> availableDirectBackends();
    /// Explicit environment override, otherwise the first available backend.
    static std::string selectedBackend();
    /// Empty selects selectedBackend(); a failed solve never changes backends.
    explicit LinearSolver(const std::string& backend = "");
    ~LinearSolver();
    const std::string& backend() const noexcept { return backend_; }
    /**
     * @brief Solve the linear system A * x = b.
     *
     * Prescribed diagonal rows are eliminated by default, with the full RHS
     * correction and back-substitution. The returned vector retains A.rows().
     * Reuses Eigen's symbolic analysis when consecutive reduced solves have the same
     * sparse structure. Numerical factorisation is repeated when coefficient
     * values change; an identical matrix reuses its factors for a new RHS.
     *
     * @throws std::invalid_argument if dimensions are inconsistent.
     * @throws std::runtime_error if the analysis, factorisation, or solve fails.
     */
    VectorXd solve(const SparseMatrixd& A, const VectorXd& b);
    /// Prescribed rows identified before a coordinate transformation or border.
    /// Preserves their affine coupling, RHS correction and full back-substitution.
    VectorXd solve(const SparseMatrixd& A, const VectorXd& b,
                   const std::vector<int>& constraintRows);
    /// Enabled by default; VELA_CONSTRAINT_ELIMINATION=0 selects the full solve.
    bool constraintEliminationEnabled() const noexcept { return eliminateConstraints_; }
    void setConstraintElimination(bool enabled) noexcept { eliminateConstraints_ = enabled; }
    Eigen::Index lastReducedDimension() const noexcept { return lastReducedDimension_; }

    /**
     * @brief Clear the cached sparse pattern and force re-analysis next solve.
     */
    void clearPatternCache();
    /// Force fresh numerical factorization while retaining symbolic analysis.
    void clearNumericCache();

    /**
     * @brief Number of symbolic pattern analyses performed by this instance.
     *
     * Counts backend analyses of the reduced matrix. An all-prescribed system
     * needs no factorization. This does not count numerical factorisations.
     */
    std::size_t patternAnalysisCount() const noexcept;

private:
    struct UmfPackState;
    bool eliminateConstraints_ = true;
    Eigen::Index lastReducedDimension_ = 0;
    VectorXd solveUnreduced(const SparseMatrixd& A, const VectorXd& b);
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
