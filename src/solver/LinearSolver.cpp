#include "vela/solver/LinearSolver.h"
#include "vela/core/PerformanceProfiler.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela {

namespace {

std::string backendFromEnvironment()
{
    const char* env = std::getenv("VELA_LINEAR_SOLVER");
    if (env == nullptr || *env == '\0')
        return "sparselu";
    return std::string(env);
}

std::string sparseMatrixDiagnostics(const SparseMatrixd& A, const VectorXd& b)
{
    std::vector<double> rowAbs(static_cast<std::size_t>(A.rows()), 0.0);
    std::vector<double> colAbs(static_cast<std::size_t>(A.cols()), 0.0);
    std::vector<int> rowNonzeroCount(static_cast<std::size_t>(A.rows()), 0);
    std::vector<int> rowOffdiagNonzeroCount(static_cast<std::size_t>(A.rows()), 0);
    int nonfiniteEntries = 0;
    std::vector<std::string> nonfiniteEntryPositions;
    for (int outer = 0; outer < A.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator it(A, outer); it; ++it) {
            const double value = it.value();
            if (!std::isfinite(value)) {
                ++nonfiniteEntries;
                if (nonfiniteEntryPositions.size() < 8) {
                    nonfiniteEntryPositions.push_back(
                        std::to_string(it.row()) + ":" + std::to_string(it.col()));
                }
                continue;
            }
            if (value != 0.0) {
                ++rowNonzeroCount[static_cast<std::size_t>(it.row())];
                if (it.row() != it.col())
                    ++rowOffdiagNonzeroCount[static_cast<std::size_t>(it.row())];
            }
            const double absValue = std::abs(value);
            rowAbs[static_cast<std::size_t>(it.row())] += absValue;
            colAbs[static_cast<std::size_t>(it.col())] += absValue;
        }
    }

    int zeroRows = 0;
    int zeroCols = 0;
    std::vector<int> zeroRowIndices;
    std::vector<int> zeroColIndices;
    for (std::size_t i = 0; i < rowAbs.size(); ++i) {
        if (rowAbs[i] == 0.0) {
            ++zeroRows;
            if (zeroRowIndices.size() < 8)
                zeroRowIndices.push_back(static_cast<int>(i));
        }
    }
    for (std::size_t i = 0; i < colAbs.size(); ++i) {
        if (colAbs[i] == 0.0) {
            ++zeroCols;
            if (zeroColIndices.size() < 8)
                zeroColIndices.push_back(static_cast<int>(i));
        }
    }

    double diagMinAbs = std::numeric_limits<double>::infinity();
    double diagMaxAbs = 0.0;
    int zeroDiagonal = 0;
    int nonfiniteDiagonal = 0;
    int singletonDiagonalRows = 0;
    for (int i = 0; i < A.rows(); ++i) {
        const double diag = A.coeff(i, i);
        if (!std::isfinite(diag)) {
            ++nonfiniteDiagonal;
            continue;
        }
        const double absDiag = std::abs(diag);
        diagMinAbs = std::min(diagMinAbs, absDiag);
        diagMaxAbs = std::max(diagMaxAbs, absDiag);
        if (absDiag == 0.0)
            ++zeroDiagonal;

        const bool onlyUnitDiagonal = std::abs(diag - 1.0) <= 1.0e-14 &&
            rowNonzeroCount[static_cast<std::size_t>(i)] == 1 &&
            rowOffdiagNonzeroCount[static_cast<std::size_t>(i)] == 0;
        if (onlyUnitDiagonal)
            ++singletonDiagonalRows;
    }
    if (!std::isfinite(diagMinAbs))
        diagMinAbs = 0.0;

    int rhsNonfinite = 0;
    std::vector<int> rhsNonfiniteIndices;
    for (int i = 0; i < b.size(); ++i) {
        if (!std::isfinite(b(i))) {
            ++rhsNonfinite;
            if (rhsNonfiniteIndices.size() < 8)
                rhsNonfiniteIndices.push_back(i);
        }
    }

    auto appendIndices = [](std::ostringstream& stream, const auto& values) {
        stream << '[';
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i != 0)
                stream << ',';
            stream << values[i];
        }
        stream << ']';
    };

    std::ostringstream out;
    out << " matrix_diagnostics={rows=" << A.rows()
        << ", cols=" << A.cols()
        << ", nnz=" << A.nonZeros()
        << ", nonfinite_entries=" << nonfiniteEntries
        << ", rhs_nonfinite=" << rhsNonfinite
        << ", zero_rows=" << zeroRows
        << ", zero_cols=" << zeroCols
        << ", nonfinite_entry_positions=";
    appendIndices(out, nonfiniteEntryPositions);
    out << ", rhs_nonfinite_indices=";
    appendIndices(out, rhsNonfiniteIndices);
    out << ", zero_row_indices=";
    appendIndices(out, zeroRowIndices);
    out << ", zero_col_indices=";
    appendIndices(out, zeroColIndices);
    out << ", diag_min_abs=" << diagMinAbs
        << ", diag_max_abs=" << diagMaxAbs
        << ", zero_diagonal=" << zeroDiagonal
        << ", nonfinite_diagonal=" << nonfiniteDiagonal
        << ", singleton_unit_diagonal_rows=" << singletonDiagonalRows
        << '}';
    return out.str();
}
VectorXd solveWithAlternateBackend(const std::string& backend,
                                   const SparseMatrixd& A,
                                   const VectorXd& b)
{
    auto checkSolved = [&backend](Eigen::ComputationInfo info, const char* phase) {
        if (info != Eigen::Success)
            throw std::runtime_error(
                "LinearSolver[" + backend + "]: " + phase + " failed.");
    };

    if (backend == "sparseqr") {
        Eigen::SparseQR<SparseMatrixd, Eigen::COLAMDOrdering<SparseMatrixd::StorageIndex>> qr;
        qr.compute(A);
        checkSolved(qr.info(), "factorisation");
        VectorXd x = qr.solve(b);
        checkSolved(qr.info(), "solve");
        return x;
    }
    if (backend == "bicgstab_ilut") {
        Eigen::BiCGSTAB<SparseMatrixd, Eigen::IncompleteLUT<double>> solver;
        solver.setTolerance(1e-14);
        solver.setMaxIterations(4000);
        solver.compute(A);
        checkSolved(solver.info(), "preconditioner setup");
        VectorXd x = solver.solve(b);
        checkSolved(solver.info(), "solve");
        return x;
    }
    if (backend == "gmres_ilut") {
        Eigen::GMRES<SparseMatrixd, Eigen::IncompleteLUT<double>> solver;
        solver.setTolerance(1e-14);
        solver.setMaxIterations(4000);
        solver.set_restart(200);
        solver.compute(A);
        checkSolved(solver.info(), "preconditioner setup");
        VectorXd x = solver.solve(b);
        checkSolved(solver.info(), "solve");
        return x;
    }
    if (backend == "simplicial_ldlt") {
        Eigen::SimplicialLDLT<SparseMatrixd> ldlt;
        ldlt.compute(A);
        checkSolved(ldlt.info(), "factorisation");
        VectorXd x = ldlt.solve(b);
        checkSolved(ldlt.info(), "solve");
        return x;
    }
    throw std::runtime_error(
        "LinearSolver: unknown VELA_LINEAR_SOLVER backend '" + backend + "'.");
}

} // namespace

VectorXd LinearSolver::solve(const SparseMatrixd& A, const VectorXd& b)
{
    ScopedPerformanceTimer totalTimer("linear.total");
    incrementPerformanceCounter("linear.solve_calls");
    observePerformanceValue("linear.rows", static_cast<double>(A.rows()));
    observePerformanceValue("linear.nonzeros", static_cast<double>(A.nonZeros()));
    if (A.rows() != A.cols())
        throw std::invalid_argument("LinearSolver: matrix must be square.");
    if (A.rows() != b.size())
        throw std::invalid_argument("LinearSolver: RHS size must match matrix rows.");

    SparseMatrixd compressed;
    const SparseMatrixd* matrix = &A;
    if (!A.isCompressed()) {
        compressed = A;
        compressed.makeCompressed();
        matrix = &compressed;
    }

    static const std::string backend = backendFromEnvironment();
    if (backend != "sparselu")
        return solveWithAlternateBackend(backend, *matrix, b);

    // A repeated solve with the same matrix (Newton line searches, feedback
    // state substitution) currently refactorises identical values.  Reusing the
    // retained factorisation is numerically identical: the same L and U are
    // used for the back-substitution.
    const bool reuseFactorization = hasFactorization_ && patternMatches(*matrix)
        && valuesMatchFactorization(*matrix);

    {
        ScopedPerformanceTimer timer("linear.analyze");
        analyzePatternIfNeeded(*matrix);
    }
    {
        ScopedPerformanceTimer timer("linear.factorize");
        if (reuseFactorization) {
            incrementPerformanceCounter("linear.factorize_cache_hits");
        } else {
            incrementPerformanceCounter("linear.factorize_calls");
            hasFactorization_ = false;
            cachedValues_.clear();
            solver_.factorize(*matrix);
        }
    }

    if (!reuseFactorization) {
        if (solver_.info() != Eigen::Success)
            throw std::runtime_error(
                "LinearSolver: SparseLU factorisation failed. "
                "Matrix may be singular or ill-conditioned." +
                sparseMatrixDiagnostics(*matrix, b));
        cacheFactorizationValues(*matrix);
    }

    const double factorNonzerosL = static_cast<double>(solver_.nnzL());
    const double factorNonzerosU = static_cast<double>(solver_.nnzU());
    observePerformanceValue("linear.factor_nonzeros_l", factorNonzerosL);
    observePerformanceValue("linear.factor_nonzeros_u", factorNonzerosU);
    observePerformanceValue(
        "linear.factor_fill_ratio",
        matrix->nonZeros() > 0
            ? (factorNonzerosL + factorNonzerosU) /
                static_cast<double>(matrix->nonZeros())
            : 0.0);

    VectorXd x;
    {
        ScopedPerformanceTimer timer("linear.solve");
        x = solver_.solve(b);
    }

    if (solver_.info() != Eigen::Success)
        throw std::runtime_error(
            "LinearSolver: SparseLU back-substitution failed." +
            sparseMatrixDiagnostics(*matrix, b));

    return x;
}

void LinearSolver::clearPatternCache()
{
    hasAnalyzedPattern_ = false;
    cachedRows_ = 0;
    cachedCols_ = 0;
    cachedNonZeros_ = 0;
    cachedOuterStarts_.clear();
    cachedInnerIndices_.clear();
    cachedValues_.clear();
    hasFactorization_ = false;
}

bool LinearSolver::valuesMatchFactorization(const SparseMatrixd& A) const
{
    const auto valueCount = static_cast<std::size_t>(A.nonZeros());
    if (!hasFactorization_ || cachedValues_.size() != valueCount)
        return false;
    if (valueCount == 0)
        return true;
    // Bitwise comparison: values that differ only in the sign of a zero, or any
    // NaN payload, are treated as a miss so the factorisation is redone.
    return std::memcmp(A.valuePtr(), cachedValues_.data(),
                       valueCount * sizeof(double)) == 0;
}

void LinearSolver::cacheFactorizationValues(const SparseMatrixd& A)
{
    const auto valueCount = static_cast<std::size_t>(A.nonZeros());
    cachedValues_.assign(A.valuePtr(), A.valuePtr() + valueCount);
    hasFactorization_ = true;
}

std::size_t LinearSolver::patternAnalysisCount() const noexcept
{
    return patternAnalysisCount_;
}

bool LinearSolver::patternMatches(const SparseMatrixd& A) const
{
    if (!hasAnalyzedPattern_)
        return false;
    if (A.rows() != cachedRows_ || A.cols() != cachedCols_)
        return false;
    if (static_cast<std::size_t>(A.nonZeros()) != cachedNonZeros_)
        return false;

    const auto outerCount = static_cast<std::size_t>(A.outerSize() + 1);
    if (cachedOuterStarts_.size() != outerCount)
        return false;
    if (!std::equal(A.outerIndexPtr(), A.outerIndexPtr() + outerCount,
                    cachedOuterStarts_.begin()))
        return false;

    const auto innerCount = static_cast<std::size_t>(A.nonZeros());
    if (cachedInnerIndices_.size() != innerCount)
        return false;
    return std::equal(A.innerIndexPtr(), A.innerIndexPtr() + innerCount,
                      cachedInnerIndices_.begin());
}

void LinearSolver::cachePattern(const SparseMatrixd& A)
{
    cachedRows_ = A.rows();
    cachedCols_ = A.cols();
    cachedNonZeros_ = static_cast<std::size_t>(A.nonZeros());

    const auto outerCount = static_cast<std::size_t>(A.outerSize() + 1);
    cachedOuterStarts_.assign(A.outerIndexPtr(), A.outerIndexPtr() + outerCount);

    const auto innerCount = static_cast<std::size_t>(A.nonZeros());
    cachedInnerIndices_.assign(A.innerIndexPtr(), A.innerIndexPtr() + innerCount);
    hasAnalyzedPattern_ = true;
}

void LinearSolver::analyzePatternIfNeeded(const SparseMatrixd& A)
{
    if (patternMatches(A)) {
        incrementPerformanceCounter("linear.analyze_cache_hits");
        return;
    }

    incrementPerformanceCounter("linear.analyze_calls");
    solver_.analyzePattern(A);
    if (!solver_.analysisIsOk()) {
        clearPatternCache();
        throw std::runtime_error(
            "LinearSolver: SparseLU symbolic analysis failed." +
            sparseMatrixDiagnostics(A, VectorXd::Zero(A.rows())));
    }

    cachePattern(A);
    ++patternAnalysisCount_;
}

} // namespace vela
