#include "vela/solver/LinearSolver.h"
#include "vela/solver/ExtendedSparseLu.h"
#include "DirectBackend.h"
#include "LinearCapture.h"
#if defined(VELA_HAS_UMFPACK)
#include <Eigen/UmfPackSupport>
#endif
#if defined(VELA_HAS_SPQR)
#include <Eigen/SPQRSupport>
#endif
#include "vela/core/PerformanceProfiler.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
extern "C" void openblas_set_num_threads(int);
extern "C" int openblas_get_num_threads();
#endif

namespace vela {

namespace {

// Use the same exact structure/value cache as other persistent direct solvers.
// Iterative refinement changes only the RHS; rebuilding these wide factors
// on every correction is unnecessary and can dominate a TCAD target.
class LongDoubleSparseLU final : public detail::DirectBackend {
    using Matrix=Eigen::SparseMatrix<long double>;
    using Vector=Eigen::Matrix<long double,Eigen::Dynamic,1>;
    Eigen::SparseLU<Matrix,Eigen::COLAMDOrdering<int>> solver_;
public:
    void analyze(const SparseMatrixd& a) override {
        if(std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits)
            throw std::runtime_error("sparselu_long_double requires long double wider than double.");
        Matrix wide=a.cast<long double>();solver_.analyzePattern(wide);
    }
    void factor(const SparseMatrixd& a) override {
        Matrix wide=a.cast<long double>();solver_.factorize(wide);
        if(solver_.info()!=Eigen::Success)throw std::runtime_error("long double factorisation failed");
    }
    VectorXd solve(const VectorXd& b) override {
        Vector wide=solver_.solve(b.cast<long double>());
        VectorXd result=wide.cast<double>();
        if(solver_.info()!=Eigen::Success || !wide.allFinite() || !result.allFinite())
            throw std::runtime_error("long double solve failed or cannot be represented");
        return result;
    }
};

bool sparseAllFinite(const SparseMatrixd& matrix) {
    // Uncompressed Eigen matrices can contain unused capacity between columns.
    // Only logical entries belong to the input system.
    for (int col=0; col<matrix.outerSize(); ++col)
        for (SparseMatrixd::InnerIterator entry(matrix,col); entry; ++entry)
            if (!std::isfinite(entry.value())) return false;
    return true;
}

int requestedBlasThreads() {
    const char* value=std::getenv("VELA_BLAS_THREADS");
    if(!value || !*value) return 0;
    if(std::strcmp(value,"1")==0) return 1;
    if(std::strcmp(value,"2")==0) return 2;
    if(std::strcmp(value,"4")==0) return 4;
    throw std::invalid_argument("VELA_BLAS_THREADS must be 1, 2 or 4");
}
void verifyBlasThreads([[maybe_unused]] bool recordObservation = true) {
    const int requested=requestedBlasThreads();
    if(!requested) return;
#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
    const int actual=openblas_get_num_threads();
    if(actual!=requested) throw std::runtime_error("OpenBLAS thread configuration changed during solve");
    if (recordObservation) observePerformanceValue("linear.openblas_threads",actual);
#else
    throw std::invalid_argument("OpenBLAS thread control unavailable in this build");
#endif
}
void configureBlasThreads() {
    const int requested=requestedBlasThreads();
    if(!requested) return;
#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
    openblas_set_num_threads(requested);
    verifyBlasThreads(false);
    static std::atomic<bool> logged{false};
    if(!logged.exchange(true)) std::fprintf(stderr,"VELA_BLAS_THREADS_VERIFIED requested=%d actual=%d\n",requested,requested);
#else
    verifyBlasThreads(false);
#endif
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
                                   const VectorXd& b,LinearSolver* retainedReduced=nullptr)
{
    auto checkSolved = [&backend](Eigen::ComputationInfo info, const char* phase) {
        if (info != Eigen::Success)
            throw std::runtime_error(
                "LinearSolver[" + backend + "]: " + phase + " failed.");
    };

    const bool keepPsi = backend == "pardiso_mkl_eliminated_keep_psi" ||
                         backend == "sparselu_eliminated_keep_psi";
    if (backend == "sparselu_eliminated" || backend == "sparselu_eliminated_long_double" ||
        backend == "umfpack_eliminated" || backend == "pardiso_mkl_eliminated" || keepPsi) {
        // Eliminate only exact numerical diagonal-singleton equations. Keep
        // their actual (possibly regularized) diagonal and affine increment.
        const int n = static_cast<int>(A.rows());
        // Explicitly opted-in three equal nodal blocks: [psi, phin, phip].
        // Keep the complete psi block, including its affine Dirichlet rows.
        // This is not valid for interleaved or additional-equation layouts.
        if (keepPsi && (n == 0 || n % 3 != 0))
            throw std::invalid_argument("keep_psi elimination requires three equal nodal blocks");
        const int retainedPrefix = keepPsi ? n / 3 : 0;
        std::vector<int> count(n, 0), map(n, -1);
        VectorXd diagonal = A.diagonal();
        for (int j=0;j<A.outerSize();++j)
            for (SparseMatrixd::InnerIterator it(A,j);it;++it) {
                if (!std::isfinite(it.value())) throw std::runtime_error("elimination: nonfinite matrix");
                if (it.value()!=0.) ++count[it.row()];
            }
        if (!b.allFinite()) throw std::runtime_error("elimination: nonfinite RHS");
        std::vector<long double> fixed(n, 0.L), corrected(n, 0.L);
        int nf=0;
        for (int i=0;i<n;++i) {
            if (i>=retainedPrefix && count[i]==1 && diagonal[i]!=0.) fixed[i]=static_cast<long double>(b[i])/diagonal[i];
            else map[i]=nf++;
            corrected[i]=b[i];
        }
        std::vector<Eigen::Triplet<double>> entries;
        entries.reserve(A.nonZeros());
        for (int j=0;j<A.outerSize();++j)
            for (SparseMatrixd::InnerIterator it(A,j);it;++it) {
                const int i=static_cast<int>(it.row());
                if (map[i]<0 || it.value()==0.) continue;
                if (map[j]<0) corrected[i]-=static_cast<long double>(it.value())*fixed[j];
                else entries.emplace_back(map[i],map[j],it.value());
            }
        VectorXd x(n), rhs(nf);
        for (int i=0;i<n;++i) {
            if(map[i]<0) x[i]=static_cast<double>(fixed[i]);
            else rhs[map[i]]=static_cast<double>(corrected[i]);
        }
        if (!rhs.allFinite()) throw std::runtime_error("elimination: nonfinite corrected RHS");
        observePerformanceValue("linear.eliminated_dofs",n-nf);
        observePerformanceValue("linear.reduced_dimension",nf);
        observePerformanceValue("linear.retained_prefix",retainedPrefix);
        if(nf) {
            SparseMatrixd reduced(nf,nf); reduced.setFromTriplets(entries.begin(),entries.end());
            reduced.makeCompressed();
            VectorXd freeX;
            if(backend=="sparselu_eliminated_long_double")
                freeX=solveWithAlternateBackend("sparselu_long_double",reduced,rhs);
            else if(backend=="umfpack_eliminated" || backend=="pardiso_mkl_eliminated" ||
                    backend=="pardiso_mkl_eliminated_keep_psi") {
                // Explicit backend selection avoids changing the caller's environment.
                // Reuse the same affine reduction and lift for both vendor solvers.
                if(retainedReduced)freeX=retainedReduced->solve(reduced,rhs);
                else {LinearSolver reducedSolver(backend=="umfpack_eliminated" ? "umfpack" : "pardiso_mkl");freeX=reducedSolver.solve(reduced,rhs);}
            }
            else {
                Eigen::SparseLU<SparseMatrixd,Eigen::COLAMDOrdering<int>> lu;
                lu.compute(reduced);checkSolved(lu.info(),"reduced factorisation");
                freeX=lu.solve(rhs);checkSolved(lu.info(),"reduced solve");
            }
            for(int i=0;i<n;++i) if(map[i]>=0) x[i]=freeX[map[i]];
        }
        if(!x.allFinite()) throw std::runtime_error("elimination: nonfinite solution");
        return x;
    }
    if (backend == "sparselu_long_double") {
        if (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
            throw std::runtime_error("sparselu_long_double requires long double wider than double.");
        const VectorXd x = solveExtendedSparseLu(A.cast<long double>(), b.cast<long double>()).cast<double>();
        if (!x.allFinite())
            throw std::runtime_error("sparselu_long_double produced an unrepresentable solution.");
        // Assembly and the public state remain binary64. No claim of a full
        // extended-precision nonlinear solve is made by this explicit backend.
        return x;
    }
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

struct LinearSolver::UmfPackState {
#if defined(VELA_HAS_UMFPACK)
    // Eigen retains matrix storage for iterative refinement. Own it for the
    // entire numeric-factor lifetime, including repeated RHS-only solves.
    SparseMatrixd matrix;
    class Solver : public Eigen::UmfPackLU<SparseMatrixd> {
    public:
        double statistic(int index) const { return this->m_umfpackInfo[index]; }
    } solver;
#endif
};

std::vector<std::string> LinearSolver::availableDirectBackends()
{
    std::vector<std::string> names;
#if defined(VELA_HAS_UMFPACK)
    names.emplace_back("umfpack");
#endif
    names.emplace_back("sparselu");
#if defined(VELA_HAS_STRUMPACK)
    names.emplace_back("strumpack");
#endif
#if defined(VELA_HAS_MUMPS)
    names.emplace_back("mumps");
#endif
#if defined(VELA_HAS_SUPERLU_MT)
    names.emplace_back("superlu_mt");
#endif
    return names;
}

std::string LinearSolver::selectedBackend()
{
    const char* env = std::getenv("VELA_LINEAR_SOLVER");
    return env && *env ? std::string(env) : availableDirectBackends().front();
}

LinearSolver::LinearSolver(const std::string& backend)
    : backend_(backend.empty() ? selectedBackend() : backend)
{
    configureBlasThreads();
    if (const char* flag = std::getenv("VELA_CONSTRAINT_ELIMINATION")) {
        if (std::strcmp(flag, "0") == 0) eliminateConstraints_ = false;
        else if (*flag && std::strcmp(flag, "1") != 0)
            throw std::invalid_argument("VELA_CONSTRAINT_ELIMINATION must be 0 or 1");
    }
    if(const char* directory=std::getenv("VELA_LINEAR_CAPTURE_DIR")) captureDirectory_=directory;
    if (const char* flag = std::getenv("VELA_LINEAR_FACTOR_STATISTICS")) {
        if (std::strcmp(flag, "0") == 0) factorStatistics_ = false;
        else if (*flag && std::strcmp(flag, "1") != 0)
            throw std::invalid_argument("VELA_LINEAR_FACTOR_STATISTICS must be 0 or 1");
    }
    if (backend_ == "sparselu_long_double" || backend_ == "pardiso_mkl" || backend_ == "sparselu_metis" || backend_ == "mumps" || backend_ == "mumps_metis" || backend_ == "superlu_mt" || backend_ == "superlu_mt_metis" || backend_ == "strumpack") createDirectBackend();
    if (backend_ == "umfpack" || backend_ == "umfpack_metis") {
#if defined(VELA_HAS_UMFPACK)
        umfpack_ = std::make_unique<UmfPackState>();
        if (backend_ == "umfpack_metis")
            umfpack_->solver.umfpackControl()[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
#else
        throw std::invalid_argument("LinearSolver: UMFPACK was not enabled in this build");
#endif
    }
}

LinearSolver::~LinearSolver() = default;

void LinearSolver::createDirectBackend() {
    if (backend_ == "sparselu_long_double") direct_=std::make_unique<LongDoubleSparseLU>();
    else if (backend_ == "pardiso_mkl") {
#ifdef VELA_HAS_MKL_PARDISO
        direct_ = detail::makeMklPardiso();
#else
        throw std::invalid_argument(
            "LinearSolver: pardiso_mkl was not enabled in this build; "
            "configure with VELA_ENABLE_MKL_PARDISO=ON");
#endif
    }
    else if (backend_ == "sparselu_metis") direct_ = detail::makeMetisSparseLU();
    else if (backend_ == "mumps" || backend_ == "mumps_metis") direct_ = detail::makeMumps(backend_ == "mumps_metis");
    else if (backend_ == "superlu_mt" || backend_ == "superlu_mt_metis") direct_ = detail::makeSuperLuMt(backend_ == "superlu_mt_metis");
    else if (backend_ == "strumpack") direct_ = detail::makeStrumpack();
}

VectorXd LinearSolver::solveDirect(const SparseMatrixd& a, const VectorXd& b) {
    try {
        if(a.rows()==0) throw std::invalid_argument("empty system");
        if (!b.allFinite() || !Eigen::Map<const VectorXd>(a.valuePtr(),a.nonZeros()).allFinite())
            throw std::runtime_error("nonfinite input");
        // Exact constraint substitution can leave one scalar. Avoid invoking
        // vendor graph analysis for that degenerate graph (MKL's Linux
        // analysis path can reject it); this is an exact division, not a
        // change of backend or a regularized fallback.
        if(a.rows()==1) {
            const Real pivot=a.coeff(0,0);
            if(pivot==0.)throw std::runtime_error("singular scalar system");
            VectorXd x(1);x[0]=b[0]/pivot;
            if(!x.allFinite())throw std::runtime_error("nonfinite scalar solution");
            incrementPerformanceCounter("linear.scalar_solve_calls");
            return x;
        }
        if (!direct_) createDirectBackend();
        const bool same=patternMatches(a), reuse=same && valuesMatchFactorization(a);
        if(!reuse) {
            // Some vendor symbolic paths assume every row and column is present.
            // Reject these provably singular inputs before calling that code.
            std::vector<bool> rows(static_cast<std::size_t>(a.rows()),false);
            for(int j=0;j<a.outerSize();++j) {
                bool nonzero=false;
                for(SparseMatrixd::InnerIterator it(a,j);it;++it) if(it.value()!=0.) {
                    rows[it.row()]=true;nonzero=true;
                }
                if(!nonzero) throw std::runtime_error("singular system: empty numerical column");
            }
            if(std::find(rows.begin(),rows.end(),false)!=rows.end())
                throw std::runtime_error("singular system: empty numerical row");
        }
        {
            ScopedPerformanceTimer timer("linear.analyze");
            if(same) incrementPerformanceCounter("linear.analyze_cache_hits");
            else {
                hasFactorization_=false;
                // A changed structure requires a new vendor instance.
                createDirectBackend();
                incrementPerformanceCounter("linear.analyze_calls");direct_->analyze(a);
                cachePattern(a);++patternAnalysisCount_;
            }
        }
        {
            ScopedPerformanceTimer timer("linear.factorize");
            if(reuse) incrementPerformanceCounter("linear.factorize_cache_hits");
            else {
                hasFactorization_=false;cachedValues_.clear();
                incrementPerformanceCounter("linear.factorize_calls");direct_->factor(a);
                cacheFactorizationValues(a);
            }
        }
        if(!reuse && factorStatistics_ && activePerformanceProfiler()) {
            ScopedPerformanceTimer timer("linear.factor_statistics");direct_->statistics();
        }
        ScopedPerformanceTimer timer("linear.solve");
        VectorXd x=direct_->solve(b);
        if(x.size()!=b.size() || !x.allFinite()) throw std::runtime_error("invalid solution");
        return x;
    } catch(const std::exception& e) {
        const std::string message=e.what();clearPatternCache();
        throw std::runtime_error("LinearSolver["+backend_+"]: "+message+sparseMatrixDiagnostics(a,b));
    }
}

VectorXd LinearSolver::solveUmfPack(const SparseMatrixd& A, const VectorXd& b)
{
#if defined(VELA_HAS_UMFPACK)
    const auto fail = [&](const char* phase) {
        clearPatternCache();
        throw std::runtime_error(std::string("LinearSolver: UMFPACK ") + phase +
            " failed." + sparseMatrixDiagnostics(A, b));
    };
    if (!b.allFinite() || !Eigen::Map<const VectorXd>(A.valuePtr(), A.nonZeros()).allFinite())
        fail("nonfinite input check");
    const bool samePattern = patternMatches(A);
    const bool reuse = samePattern && valuesMatchFactorization(A);
    if (!reuse) umfpack_->matrix = A;
    auto& lu = umfpack_->solver;
    {
        ScopedPerformanceTimer timer("linear.analyze");
        if (samePattern) incrementPerformanceCounter("linear.analyze_cache_hits");
        else {
            hasFactorization_ = false;
            incrementPerformanceCounter("linear.analyze_calls");
            lu.analyzePattern(umfpack_->matrix);
            if (lu.info() != Eigen::Success) fail("symbolic analysis");
            if (backend_ == "umfpack_metis" && lu.statistic(UMFPACK_ORDERING_USED) != UMFPACK_ORDERING_METIS)
                fail("requested METIS ordering unavailable");
            cachePattern(A);
            ++patternAnalysisCount_;
        }
    }
    {
        ScopedPerformanceTimer timer("linear.factorize");
        if (reuse) incrementPerformanceCounter("linear.factorize_cache_hits");
        else {
            hasFactorization_ = false;
            cachedValues_.clear();
            incrementPerformanceCounter("linear.factorize_calls");
            lu.factorize(umfpack_->matrix);
            if (lu.info() != Eigen::Success) fail("factorisation");
            cacheFactorizationValues(A);
        }
    }
    if (!reuse && factorStatistics_ && activePerformanceProfiler()) {
        ScopedPerformanceTimer diagnostics("linear.factor_statistics");
        const auto emit = [&](const char* name, int index) {
            const double value = lu.statistic(index);
            if (std::isfinite(value) && value >= 0.) observePerformanceValue(name, value);
        };
        emit("linear.numeric_factor_nonzeros_l", UMFPACK_LNZ);
        emit("linear.numeric_factor_nonzeros_u", UMFPACK_UNZ);
        emit("linear.umfpack_reported_flops", UMFPACK_FLOPS);
        emit("linear.umfpack_strategy_used", UMFPACK_STRATEGY_USED);
        emit("linear.umfpack_ordering_used", UMFPACK_ORDERING_USED);
        emit("linear.umfpack_max_front_entries", UMFPACK_MAX_FRONT_SIZE);
        observePerformanceValue("linear.numeric_factor_fill_ratio",
            (lu.statistic(UMFPACK_LNZ) + lu.statistic(UMFPACK_UNZ)) / A.nonZeros());
        observePerformanceValue("linear.umfpack_internal_peak_bytes",
            lu.statistic(UMFPACK_PEAK_MEMORY) * lu.statistic(UMFPACK_SIZE_OF_UNIT));
    }
    VectorXd x(A.cols());
    {
        ScopedPerformanceTimer timer("linear.solve");
        // Eigen's solve expression discards this adapter's boolean status.
        // Check it explicitly, then reject nonfinite solutions as well.
        if (!lu._solve_impl(b, x) || !x.allFinite()) fail("back-substitution");
    }
    if (factorStatistics_)
        observePerformanceValue("linear.umfpack_refinement_steps", lu.statistic(UMFPACK_IR_TAKEN));
    return x;
#else
    (void)A; (void)b;
    throw std::invalid_argument("LinearSolver: UMFPACK was not enabled in this build");
#endif
}

bool solveSparseQrSystem(const SparseMatrixd& A,
                         const VectorXd& b,
                         VectorXd& x) noexcept
{
    try {
        x = solveWithAlternateBackend("sparseqr", A, b);
        return x.size() == A.cols() && x.allFinite();
    } catch (const std::exception&) {
        return false;
    }
}

bool solveUmfPackSystem(const SparseMatrixd& A,
                        const VectorXd& b,
                        VectorXd& x) noexcept
{
#if defined(VELA_HAS_UMFPACK)
    try {
        Eigen::UmfPackLU<SparseMatrixd> solver;
        solver.umfpackControl()[UMFPACK_STRATEGY] =
            UMFPACK_STRATEGY_UNSYMMETRIC;
        solver.umfpackControl()[UMFPACK_PIVOT_TOLERANCE] = 1.0e-12;
        solver.compute(A);
        if (solver.info() == Eigen::InvalidInput)
            return false;
        x = solver.solve(b);
        return x.size() == A.cols() && x.allFinite();
    } catch (const std::exception&) {
        return false;
    }
#else
    (void)A;
    (void)b;
    (void)x;
    return false;
#endif
}

bool solveSpqrSystem(const SparseMatrixd& A,
                     const VectorXd& b,
                     VectorXd& x) noexcept
{
#if defined(VELA_HAS_SPQR)
    try {
        using LongSparseMatrix = Eigen::SparseMatrix<
            Real, Eigen::ColMajor, SuiteSparse_long>;
        LongSparseMatrix matrix = A;
        matrix.makeCompressed();
        VectorXd rhs = b;
        cholmod_sparse cholmodMatrix = Eigen::viewAsCholmod(matrix);
        cholmod_dense cholmodRhs = Eigen::viewAsCholmod(rhs);
        cholmod_common common;
        if (!cholmod_l_start(&common))
            return false;
        cholmod_dense* cholmodSolution =
            SuiteSparseQR_min2norm<Real, SuiteSparse_long>(
                SPQR_ORDERING_DEFAULT,
                1.0e-10,
                &cholmodMatrix,
                &cholmodRhs,
                &common);
        if (cholmodSolution == nullptr) {
            cholmod_l_finish(&common);
            return false;
        }
        x.resize(static_cast<Eigen::Index>(cholmodSolution->nrow));
        const Real* values = static_cast<const Real*>(cholmodSolution->x);
        for (Eigen::Index i = 0; i < x.size(); ++i)
            x(i) = values[i];
        cholmod_l_free_dense(&cholmodSolution, &common);
        cholmod_l_finish(&common);
        return x.size() == A.cols() && x.allFinite();
    } catch (const std::exception&) {
        return false;
    }
#else
    (void)A;
    (void)b;
    (void)x;
    return false;
#endif
}

VectorXd LinearSolver::solve(const SparseMatrixd& A, const VectorXd& b)
{
    if (!b.allFinite() || !sparseAllFinite(A))
        return solveUnreduced(A,b);
    // Named legacy elimination backends own their affine reduction, including
    // the three-block keep-psi contract. Do not eliminate psi a first time in
    // the new generic front end or reduce a second time in the vendor solver.
    if (backend_.find("_eliminated") != std::string::npos) {
        lastReducedDimension_ = A.rows();
        const bool keepPsi = backend_.find("_keep_psi") != std::string::npos;
        for (int row : diagonalConstraintRows(A))
            if (!keepPsi || row >= A.rows()/3) --lastReducedDimension_;
        return solveUnreduced(A,b);
    }
    return solve(A,b,eliminateConstraints_ ? diagonalConstraintRows(A) : std::vector<int>{});
}

VectorXd LinearSolver::solve(const SparseMatrixd& A, const VectorXd& b,
                             const std::vector<int>& constraints)
{
    if (A.rows() != A.cols() || A.rows() != b.size())
        throw std::invalid_argument("LinearSolver: invalid system dimensions");
    if (!b.allFinite() || !sparseAllFinite(A))
        return solveUnreduced(A,b);
    lastReducedDimension_ = A.rows();
    if (!eliminateConstraints_ || constraints.empty()) return solveUnreduced(A,b);
    ScopedPerformanceTimer timer("linear.constraint_elimination");
    const auto reduced = eliminateLinearConstraints(A,b,constraints);
    lastReducedDimension_ = reduced.matrix.rows();
    observePerformanceValue("linear.original_rows",static_cast<double>(A.rows()));
    observePerformanceValue("linear.eliminated_rows",static_cast<double>(constraints.size()));
    observePerformanceValue("linear.reduced_rows",static_cast<double>(reduced.matrix.rows()));
    VectorXd direction;
    try {
        if (reduced.matrix.rows()) {
            if (backend_.find("_eliminated") != std::string::npos) {
                const std::string base = backend_.substr(0,backend_.find("_eliminated")) +
                    (backend_.find("long_double") != std::string::npos ? "_long_double" : "");
                if (!eliminatedDirect_) eliminatedDirect_=std::make_unique<LinearSolver>(base);
                eliminatedDirect_->setConstraintElimination(false);
                direction=eliminatedDirect_->solve(reduced.matrix,reduced.rhs);
            } else direction=solveUnreduced(reduced.matrix,reduced.rhs);
        }
    } catch (const std::exception& error) {
        throw std::runtime_error(std::string(error.what())+"; original system"+sparseMatrixDiagnostics(A,b));
    }
    // The backend already recorded this solve; fully prescribed systems never
    // enter it. Keep the outer check without double-counting reduced solves.
    verifyBlasThreads(reduced.matrix.rows() == 0);
    return reduced.recover(direction);
}

VectorXd LinearSolver::solveUnreduced(const SparseMatrixd& A, const VectorXd& b)
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

    const auto finish=[&](VectorXd x) {
        verifyBlasThreads();
        if(!captureDirectory_.empty()) detail::captureLinearInput(captureDirectory_,*matrix,b,x);
        return x;
    };
    if (backend_ == "umfpack" || backend_ == "umfpack_metis")
        return finish(solveUmfPack(*matrix, b));
    if (backend_ == "sparselu_long_double" || backend_ == "pardiso_mkl" || backend_ == "sparselu_metis" || backend_ == "mumps" || backend_ == "mumps_metis" || backend_ == "superlu_mt" || backend_ == "superlu_mt_metis" || backend_ == "strumpack") return finish(solveDirect(*matrix,b));
    if (backend_ != "sparselu") {
        if(backend_=="umfpack_eliminated" || backend_=="pardiso_mkl_eliminated" || backend_=="pardiso_mkl_eliminated_keep_psi"){
            if(!eliminatedDirect_)eliminatedDirect_=std::make_unique<LinearSolver>(backend_=="umfpack_eliminated"?"umfpack":"pardiso_mkl");
            eliminatedDirect_->setConstraintElimination(false);
        }
        return finish(solveWithAlternateBackend(backend_, *matrix, b,eliminatedDirect_.get()));
    }

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
        if (factorStatistics_ && activePerformanceProfiler()) {
            ScopedPerformanceTimer diagnostics("linear.factor_statistics");
            observePerformanceValue("linear.numeric_factor_nonzeros_l", solver_.nnzL());
            observePerformanceValue("linear.numeric_factor_nonzeros_u", solver_.nnzU());
            observePerformanceValue("linear.numeric_factor_fill_ratio",
                double(solver_.nnzL() + solver_.nnzU()) / matrix->nonZeros());
            observePerformanceValue("linear.sparselu_structural_flops_estimate", solver_.structuralFlops());
#if defined(VELA_SPARSELU_ORDERING_AMD)
            incrementPerformanceCounter("linear.numeric_ordering_amd");
#else
            incrementPerformanceCounter("linear.numeric_ordering_colamd");
#endif
        }
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

    return finish(std::move(x));
}

void LinearSolver::clearNumericCache()
{
    if(eliminatedDirect_)eliminatedDirect_->clearNumericCache();
    cachedValues_.clear();
    hasFactorization_ = false;
}

void LinearSolver::clearPatternCache()
{
    eliminatedDirect_.reset();
    direct_.reset();
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
