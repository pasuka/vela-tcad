#include "vela/solver/LinearSolver.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace vela {

namespace {

std::string backendFromEnvironment()
{
    const char* env = std::getenv("VELA_LINEAR_SOLVER");
    if (env == nullptr || *env == '\0')
        return "sparselu";
    return std::string(env);
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

    analyzePatternIfNeeded(*matrix);
    solver_.factorize(*matrix);

    if (solver_.info() != Eigen::Success)
        throw std::runtime_error(
            "LinearSolver: SparseLU factorisation failed. "
            "Matrix may be singular or ill-conditioned.");

    VectorXd x = solver_.solve(b);

    if (solver_.info() != Eigen::Success)
        throw std::runtime_error(
            "LinearSolver: SparseLU back-substitution failed.");

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
    if (patternMatches(A))
        return;

    solver_.analyzePattern(A);
    if (!solver_.analysisIsOk()) {
        clearPatternCache();
        throw std::runtime_error("LinearSolver: SparseLU symbolic analysis failed.");
    }

    cachePattern(A);
    ++patternAnalysisCount_;
}

} // namespace vela
