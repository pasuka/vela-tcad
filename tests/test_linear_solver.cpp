#include "vela/solver/SlotPopulationNumerics.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PerformanceProfiler.h"
#include "vela/solver/ConstraintElimination.h"
#include "vela/solver/LinearSolver.h"
#include "vela/solver/ExtendedSparseLu.h"
#include "vela/solver/LinearRefinement.h"
#include "vela/solver/PopulationState.h"

TEST_CASE("Selected backend handles a scalar left by exact constraint substitution", "[linear_solver][constraints]") {
    using namespace vela;
    SparseMatrixd a(2,2);a.insert(0,0)=2.;a.insert(0,1)=1.;a.insert(1,0)=3.;a.insert(1,1)=4.;
    VectorXd exact(2);exact<<-3.,5.;
    LinearSolver solver;
    REQUIRE((solver.solve(a,a*exact,{0})-exact).norm()<1e-13);
    REQUIRE(solver.lastReducedDimension()==1);
    SparseMatrixd scalar(1,1);scalar.insert(0,0)=0.;
    REQUIRE_THROWS(solver.solve(scalar,VectorXd::Ones(1),{}));
}
#include "vela/solver/ElectrothermalPredictorQuality.h"
#include "vela/solver/ElectrothermalStepControl.h"
#include "vela/solver/ElectrothermalLocalPrediction.h"
#include "vela/solver/ElectrothermalResidualMixing.h"
#include "vela/solver/ElectrothermalDensityUpdate.h"
#include "vela/solver/ElectrothermalNaturalDamping.h"
#include "../src/tools/ElectrothermalIterationControl.h"

#include <Eigen/Sparse>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <limits>
#include <vector>
#include <cstdlib>

using namespace vela;

TEST_CASE("Constraint elimination preserves affine rows and nonzero boundary RHS", "[linear_solver][constraints]") {
    SparseMatrixd a(5,5);
    std::vector<Eigen::Triplet<double>> t={{0,0,2.},{0,4,-3.},
        {1,0,4.},{1,1,-5.},{1,4,6.},
        {2,0,2.},{2,1,3.},{2,2,8.},{2,3,-1.},{2,4,.2},
        {3,0,-1.},{3,2,2.},{3,3,6.},{3,4,.3},
        {4,0,.2},{4,1,-.5},{4,2,.7},{4,3,1.},{4,4,4.}};
    a.setFromTriplets(t.begin(),t.end());
    VectorXd exact(5); exact<<.31,-1.7,2.3,-.6,1.9;
    const VectorXd rhs=a*exact;
    const auto r=eliminateLinearConstraints(a,rhs,{1,0});
    CHECK(r.kept==std::vector<int>{2,3,4}); CHECK(r.matrix.rows()==3);
    LinearSolver solver("sparselu"); REQUIRE(solver.constraintEliminationEnabled());
    const auto x=solver.solve(a,rhs,{0,1});
    CHECK(solver.lastReducedDimension()==3);
    CHECK((x-exact).lpNorm<Eigen::Infinity>()<2e-15);
    CHECK((a*x-rhs).lpNorm<Eigen::Infinity>()<2e-14);
    CHECK((r.recover(solver.solve(r.matrix,r.rhs))-exact).norm()<3e-15);
    solver.setConstraintElimination(false);
    CHECK((solver.solve(a,rhs)-x).norm()<3e-15);
    CHECK(solver.lastReducedDimension()==5);
    solver.setConstraintElimination(true);
    CHECK((solver.solve(a,2.*rhs,{0,1})-2.*x).norm()<4e-15);
}

TEST_CASE("Default constraint elimination ignores stored zeros but preserves tiny couplings", "[linear_solver][constraints]") {
    SparseMatrixd a(3,3); a.insert(0,0)=2.;a.insert(0,1)=0.;
    a.insert(1,0)=1.;a.insert(1,1)=3.;a.insert(1,2)=1e-250;
    a.insert(2,1)=.5;a.insert(2,2)=4.;a.makeCompressed();
    CHECK(diagonalConstraintRows(a)==std::vector<int>{0});
    const VectorXd exact=VectorXd::LinSpaced(3,1.,3.);
    LinearSolver solver("sparselu");
    CHECK((solver.solve(a,a*exact)-exact).norm()<1e-14);
    CHECK(solver.lastReducedDimension()==2);
    SparseMatrixd diagonal(2,2);diagonal.insert(0,0)=2.;diagonal.insert(1,1)=-3.;
    VectorXd rhs(2);rhs<<6.,12.;
    const auto all=eliminateLinearConstraints(diagonal,rhs,{0,1});
    CHECK(all.matrix.rows()==0);CHECK((all.recover(VectorXd{})-solver.solve(diagonal,rhs)).norm()==0.);
    CHECK(solver.lastReducedDimension()==0);
}

TEST_CASE("Constraint substitution rejects cycles singular pivots and invalid metadata", "[linear_solver][constraints]") {
    SparseMatrixd a(2,2);a.insert(0,0)=2.;a.insert(0,1)=1.;a.insert(1,0)=1.;a.insert(1,1)=3.;
    const auto b=VectorXd::Ones(2).eval();
    CHECK_THROWS_AS(eliminateLinearConstraints(a,b,{0,1}),std::invalid_argument);
    CHECK_THROWS_AS(eliminateLinearConstraints(a,b,{0,0}),std::invalid_argument);
    CHECK_THROWS_AS(eliminateLinearConstraints(a,b,{2}),std::invalid_argument);
    a.coeffRef(0,0)=0.;CHECK_THROWS_AS(eliminateLinearConstraints(a,b,{0}),std::invalid_argument);
    a.coeffRef(0,0)=std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(diagonalConstraintRows(a),std::invalid_argument);
}

TEST_CASE("Diagonal substitution preserves ordered RHS corrections and sparse slots",
          "[linear_solver][constraints]") {
    SparseMatrixd matrix(5, 5);
    const std::vector<Eigen::Triplet<double>> entries{
        {0,0,1.}, {0,3,0.}, {1,0,1.e16}, {1,1,2.}, {1,2,-1.e16},
        {1,3,0.}, {1,4,-5.}, {2,2,1.}, {3,1,3.}, {3,3,5.},
        {3,4,0.}, {4,4,1.}};
    matrix.setFromTriplets(entries.begin(), entries.end());
    VectorXd rhs(5);
    rhs << 1., 1., 1., 3., 1.;
    const auto rows = diagonalConstraintRows(matrix);
    REQUIRE(rows == std::vector<int>{0, 2, 4});
    const auto reduced = detail::eliminateDiagonalLinearConstraints(matrix, rhs, rows);
    REQUIRE(reduced.matrix.nonZeros() == 4);
    CHECK(reduced.rhs[0] == 5.);
    CHECK(reduced.rhs[1] == 3.);
    LinearSolver solver("sparselu");
    const auto expected = solver.solve(matrix, rhs, rows);
    CHECK((solver.solve(matrix, rhs) - expected).norm() == 0.);
}

TEST_CASE("Constraint validation ignores unused uncompressed sparse capacity", "[linear_solver][constraints]") {
    SparseMatrixd a(3,3);a.reserve(Eigen::VectorXi::Constant(3,5));
    a.insert(0,0)=2.;a.insert(1,0)=.1;a.insert(1,1)=3.;a.insert(2,1)=.2;a.insert(2,2)=4.;
    REQUIRE_FALSE(a.isCompressed());
    for(int col=0;col<a.outerSize();++col)
        for(int k=a.outerIndexPtr()[col]+a.innerNonZeroPtr()[col];k<a.outerIndexPtr()[col+1];++k)
            a.valuePtr()[k]=std::numeric_limits<double>::quiet_NaN();
    const VectorXd exact=VectorXd::LinSpaced(3,1.,3.);
    LinearSolver solver("sparselu");
    CHECK((solver.solve(a,a*exact)-exact).norm()<1e-14);
}

#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
extern "C" int openblas_get_num_threads();
extern "C" void openblas_set_num_threads(int);
#endif

TEST_CASE("Explicit BLAS thread control is opt-in and observable", "[linear_solver][blas_threads]") {
    struct Guard {
        const bool existed=std::getenv("VELA_BLAS_THREADS")!=nullptr;
        const std::string old=existed?std::getenv("VELA_BLAS_THREADS"):"";
#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
        const int oldThreads=openblas_get_num_threads();
#endif
        void set(const char* value) {
#ifdef _WIN32
            _putenv_s("VELA_BLAS_THREADS",value?value:"");
#else
            if(value) setenv("VELA_BLAS_THREADS",value,1);else unsetenv("VELA_BLAS_THREADS");
#endif
        }
        ~Guard() {
            set(existed?old.c_str():nullptr);
#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
            openblas_set_num_threads(oldThreads);
#endif
        }
    } guard;
    guard.set("3");REQUIRE_THROWS_AS(LinearSolver("sparselu"),std::invalid_argument);
    guard.set("invalid");REQUIRE_THROWS_AS(LinearSolver("sparselu"),std::invalid_argument);
    guard.set("1");
#if defined(VELA_HAS_OPENBLAS_THREAD_CONTROL)
    SparseMatrixd a(3,3);a.insert(0,0)=2.;a.insert(0,1)=.3;a.insert(1,1)=4.;a.insert(2,0)=-.2;a.insert(2,2)=5.;a.makeCompressed();
    SparseMatrixd diagonal(3,3);diagonal.insert(0,0)=2.;diagonal.insert(1,1)=4.;diagonal.insert(2,2)=5.;
    VectorXd exact(3);exact<<1.,-2.,3.;
    const auto exercise = [&](const std::string& backend, int count,
                              bool eliminate, const SparseMatrixd& matrix) {
        INFO("backend=" << backend << ", threads=" << count << ", eliminate=" << eliminate
             << ", nonzeros=" << matrix.nonZeros());
        const auto requested=std::to_string(count);
        guard.set(requested.c_str());
        const VectorXd b=matrix*exact;
        const double tolerance=backend=="sparselu"?1e-14:1e-13;
        PerformanceProfiler measured({true,"unused.json"});
        {
            ActivePerformanceProfilerScope active(&measured);
            LinearSolver solver(backend);
            solver.setConstraintElimination(eliminate);
            REQUIRE(openblas_get_num_threads()==count);
            REQUIRE_FALSE(measured.toJson().at("observations").contains("linear.openblas_threads"));
            REQUIRE((solver.solve(matrix,b)-exact).norm()<tolerance);
            REQUIRE((solver.solve(matrix,2.*b)-2.*exact).norm()<tolerance);
            const auto seen=measured.toJson().at("observations").at("linear.openblas_threads");
            REQUIRE(seen.at("min")==count);REQUIRE(seen.at("max")==count);
            REQUIRE(seen.at("count")==2);
            const int changed=count==1?2:1;
            openblas_set_num_threads(changed);
            REQUIRE_THROWS(solver.solve(matrix,b));
            REQUIRE(measured.toJson().at("observations").at("linear.openblas_threads").at("count")==2);
            guard.set(nullptr);
            LinearSolver noControl("sparselu");
            REQUIRE(openblas_get_num_threads()==changed);
        }
    };
    std::vector<std::string> backends{"sparselu"};
#if defined(VELA_HAS_UMFPACK)
    backends.emplace_back("umfpack");
#endif
    for(const auto& backend : backends) {
        for(const int count : {1,2,4}) {
            exercise(backend,count,false,a);
            exercise(backend,count,true,a);
            exercise(backend,count,true,diagonal);
        }
    }
#else
    REQUIRE_THROWS_AS(LinearSolver("sparselu"),std::invalid_argument);
#endif
}

namespace {

SparseMatrixd makeSparseMatrix(
    int rows,
    int cols,
    const std::vector<Eigen::Triplet<double>>& triplets)
{
    SparseMatrixd A(rows, cols);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
}

} // namespace

TEST_CASE("Direct backend preference honors explicit overrides without failover", "[linear_solver][default]") {
    struct EnvironmentGuard {
        const std::string old = std::getenv("VELA_LINEAR_SOLVER") ? std::getenv("VELA_LINEAR_SOLVER") : "";
        void set(const char* value) {
#ifdef _WIN32
            _putenv_s("VELA_LINEAR_SOLVER", value);
#else
            setenv("VELA_LINEAR_SOLVER", value, 1);
#endif
        }
        ~EnvironmentGuard() { set(old.c_str()); }
    } env;
    std::vector<std::string> expected;
#if defined(VELA_HAS_UMFPACK)
    expected.emplace_back("umfpack");
#endif
    expected.emplace_back("sparselu");
#if defined(VELA_HAS_STRUMPACK)
    expected.emplace_back("strumpack");
#endif
#if defined(VELA_HAS_MUMPS)
    expected.emplace_back("mumps");
#endif
#if defined(VELA_HAS_SUPERLU_MT)
    expected.emplace_back("superlu_mt");
#endif
    REQUIRE(LinearSolver::availableDirectBackends() == expected);
    env.set("");
    LinearSolver solver;
    REQUIRE(solver.backend() == expected.front());
    auto a = makeSparseMatrix(2,2,{{0,0,3.},{0,1,-1.},{1,0,2.},{1,1,5.}});
    VectorXd exact(2); exact << 2., -1.;
    REQUIRE((solver.solve(a,a*exact)-exact).norm() < 1e-13);
    a.coeffRef(0,0)=4.;
    REQUIRE((solver.solve(a,a*exact)-exact).norm() < 1e-13);
    REQUIRE(solver.patternAnalysisCount()==1);
    env.set("sparselu");
    REQUIRE(LinearSolver().backend()=="sparselu");
    REQUIRE(solver.backend()==expected.front());
    env.set("unknown_backend");
    LinearSolver invalid;
    REQUIRE_THROWS(invalid.solve(a,a*exact));
    REQUIRE(invalid.backend()=="unknown_backend");
}

TEST_CASE("Optional direct backends preserve equations caches and failure recovery", "[linear_solver][backend_contract]") {
    std::vector<std::string> backends;
    if(std::numeric_limits<long double>::digits>std::numeric_limits<double>::digits)
        backends.push_back("sparselu_long_double");
#if defined(VELA_HAS_STRUMPACK)
    backends.push_back("strumpack");
#endif
#if defined(VELA_HAS_SUPERLU_MT)
    backends.push_back("superlu_mt");
#if defined(VELA_HAS_METIS)
    backends.push_back("superlu_mt_metis");
#endif
#endif
#if defined(VELA_HAS_MUMPS)
    backends.push_back("mumps");backends.push_back("mumps_metis");
#endif
#if defined(VELA_HAS_METIS)
    backends.push_back("sparselu_metis");
#endif
#if defined(VELA_HAS_UMFPACK)
    backends.push_back("umfpack_metis");
#endif
    for(const auto& backend:backends) {
        DYNAMIC_SECTION(backend) {
            LinearSolver solver(backend);
            REQUIRE_THROWS(solver.solve(SparseMatrixd(0,0),VectorXd(0)));
            auto a=makeSparseMatrix(4,4,{{0,0,0.},{0,1,2.},{1,0,3.},{1,1,4.},
                {1,2,-1.},{2,1,.3},{2,2,6.},{3,3,7.}});
            VectorXd exact(4);exact<<1.,-2.,.5,3.;
            const auto check=[&](const SparseMatrixd& m){
                VectorXd b=m*exact,x=solver.solve(m,b);
                REQUIRE((x-exact).norm()<1e-11);REQUIRE((m*x-b).norm()<1e-11);
            };
            PerformanceProfiler profiler({true,"unused.json"});
            {
                ActivePerformanceProfilerScope active(&profiler);
                check(a);auto saved=a;a.setZero();a.resize(1,1);
                REQUIRE((solver.solve(saved,2.*(saved*exact))-2.*exact).norm()<1e-11);
                a=saved;a.coeffRef(1,0)=-5.;check(a);
            }
            REQUIRE(solver.patternAnalysisCount()==1);
            const auto counts=profiler.toJson().at("counters");
            REQUIRE(counts.at("linear.factorize_calls")==2);
            REQUIRE(counts.at("linear.factorize_cache_hits")==1);
            // A new DC request must preserve the graph but refresh numeric
            // factors, even for a bitwise identical coefficient matrix.
            solver.clearNumericCache();
            PerformanceProfiler requestProfiler({true,"unused.json"});
            {
                ActivePerformanceProfilerScope active(&requestProfiler);
                check(a);
            }
            REQUIRE(solver.patternAnalysisCount()==1);
            const auto requestCounts=requestProfiler.toJson().at("counters");
            REQUIRE(requestCounts.at("linear.analyze_cache_hits")==1);
            REQUIRE(requestCounts.at("linear.factorize_calls")==1);
            REQUIRE_FALSE(requestCounts.contains("linear.factorize_cache_hits"));
            a.coeffRef(3,0)=.7;check(a);REQUIRE(solver.patternAnalysisCount()==2);
            solver.clearPatternCache();check(a);REQUIRE(solver.patternAnalysisCount()==3);
            auto singular=makeSparseMatrix(4,4,{{0,0,1.},{1,1,1.},{2,2,1.}});
            REQUIRE_THROWS(solver.solve(singular,exact));check(a);
            auto dependent=makeSparseMatrix(4,4,{{0,0,1.},{0,1,2.},{1,0,1.},{1,1,2.},{2,2,1.},{3,3,1.}});
            REQUIRE_THROWS(solver.solve(dependent,exact));check(a);
            auto bad=exact;bad[0]=std::numeric_limits<double>::quiet_NaN();
            REQUIRE_THROWS(solver.solve(a,bad));check(a);
            auto badA=a;badA.coeffRef(0,1)=std::numeric_limits<double>::infinity();
            REQUIRE_THROWS(solver.solve(badA,exact));check(a);
            // Severe row scaling must not turn the solve into a symmetric problem.
            auto scaled=a;
            const double scale[]={1e-12,1e8,1e-4,1e3};
            for(int j=0;j<scaled.outerSize();++j)
                for(SparseMatrixd::InnerIterator it(scaled,j);it;++it) it.valueRef()*=scale[it.row()];
            VectorXd rhs=scaled*exact,x=solver.solve(scaled,rhs);
            REQUIRE((x-exact).norm()<1e-10);
        }
    }
}

TEST_CASE("Disabling factor diagnostics preserves states and factor reuse", "[linear_solver][diagnostics]") {
    struct Environment {
        std::string old;
        bool existed;
        Environment() : existed(std::getenv("VELA_LINEAR_FACTOR_STATISTICS") != nullptr) {
            if (existed) old=std::getenv("VELA_LINEAR_FACTOR_STATISTICS");
        }
        void set(const char* value) {
#ifdef _WIN32
            _putenv_s("VELA_LINEAR_FACTOR_STATISTICS",value ? value : "");
#else
            if(value) setenv("VELA_LINEAR_FACTOR_STATISTICS",value,1);
            else unsetenv("VELA_LINEAR_FACTOR_STATISTICS");
#endif
        }
        ~Environment() { set(existed ? old.c_str() : nullptr); }
    } env;
    std::string backend="sparselu";
    SECTION("SparseLU") {}
#if defined(VELA_HAS_UMFPACK)
    SECTION("UMFPACK") { backend="umfpack"; }
#endif
    const auto a=makeSparseMatrix(3,3,{{0,0,4.},{0,1,.2},{1,1,5.},{1,2,-.3},{2,0,.1},{2,2,6.}});
    const VectorXd rhs=VectorXd::Ones(3);
    env.set("1"); LinearSolver on(backend);
    const VectorXd expected=on.solve(a,rhs);
    env.set("0"); LinearSolver off(backend);
    PerformanceProfiler profiler({true,"unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        REQUIRE((off.solve(a,rhs)-expected).norm()==0.);
        REQUIRE((off.solve(a,2.*rhs)-2.*expected).norm()<1e-14);
    }
    const auto json=profiler.toJson();
    REQUIRE(json.at("counters").at("linear.factorize_calls")==1);
    REQUIRE(json.at("counters").at("linear.factorize_cache_hits")==1);
    REQUIRE_FALSE(json.at("observations").contains("linear.numeric_factor_nonzeros_l"));
    REQUIRE_FALSE(json.at("observations").contains("linear.sparselu_structural_flops_estimate"));
    for(const auto& stage:json.at("stages")) REQUIRE(stage.at("name")!="linear.factor_statistics");
    env.set("invalid"); REQUIRE_THROWS_AS(LinearSolver(backend),std::invalid_argument);
}

TEST_CASE("Factor statistics preserve solves and count only fresh factors", "[linear_solver][diagnostics]") {
    std::string backend = "sparselu";
#if defined(VELA_HAS_UMFPACK)
    SECTION("UMFPACK") { backend = "umfpack"; }
#endif
    SECTION("SparseLU") {}
    LinearSolver measured(backend), plain(backend);
    const auto a = makeSparseMatrix(3,3,{{0,0,10.},{0,1,1.},{0,2,2.},
        {1,0,3.},{1,1,20.},{1,2,4.},{2,0,5.},{2,1,6.},{2,2,30.}});
    const VectorXd rhs = VectorXd::Ones(3);
    const VectorXd reference = plain.solve(a,rhs);
    PerformanceProfiler profiler({true,"unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        REQUIRE((measured.solve(a,rhs)-reference).norm()==0.);
        REQUIRE((measured.solve(a,2.*rhs)-2.*reference).norm()<1e-14);
    }
    const auto json=profiler.toJson();
    const auto obs=json.at("observations");
    REQUIRE(obs.at("linear.numeric_factor_nonzeros_l").at("count")==1);
    REQUIRE(obs.at("linear.numeric_factor_nonzeros_l").at("last")==6.);
    REQUIRE(obs.at("linear.numeric_factor_nonzeros_u").at("last")==6.);
    REQUIRE(obs.at("linear.numeric_factor_fill_ratio").at("last").get<double>()==Catch::Approx(12./9.));
    if(backend=="sparselu")
        REQUIRE(obs.at("linear.sparselu_structural_flops_estimate").at("last")==13.);
    else {
        REQUIRE(obs.at("linear.umfpack_reported_flops").at("last").get<double>()>=13.);
        REQUIRE(obs.at("linear.umfpack_internal_peak_bytes").at("last").get<double>()>0.);
    }
#if defined(_WIN32)
    REQUIRE(json.at("resources").at("process_peak_working_set_bytes").get<double>()>0.);
#endif
}

TEST_CASE("UMFPACK retains factors safely and recovers after rejected systems", "[linear_solver][umfpack]")
{
#if defined(VELA_HAS_UMFPACK)
    LinearSolver solver("umfpack");
    // This case audits the backend's full-pattern cache; reduction has its own cases.
    solver.setConstraintElimination(false);
    const VectorXd exact = (VectorXd(3) << 1., -2., 3.).finished();
    auto a = makeSparseMatrix(3,3,{{0,0,4.},{1,1,5.},{2,2,6.},{0,1,.7},{2,0,-.3}});
    const auto check = [&](const SparseMatrixd& m) {
        const VectorXd b = m * exact;
        const VectorXd x = solver.solve(m,b);
        REQUIRE((m*x-b).norm()<1e-12);
        REQUIRE((x-exact).norm()<1e-12);
    };
    PerformanceProfiler profiler({true,"unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        check(a);
        // Destroy caller storage, then reuse identical factors with a new RHS.
        const auto saved = a;
        a.setZero(); a.resize(100,100);
        const VectorXd b = saved * (2.*exact);
        REQUIRE((solver.solve(saved,b)-2.*exact).norm()<1e-12);
        a=saved; a.coeffRef(0,0)=8.;check(a);
    }
    REQUIRE(solver.patternAnalysisCount()==1);
    const auto counts=profiler.toJson().at("counters");
    REQUIRE(counts.at("linear.factorize_calls")==2);
    REQUIRE(counts.at("linear.factorize_cache_hits")==1);
    // Same nnz, different indices must invalidate symbolic analysis.
    a=makeSparseMatrix(3,3,{{0,0,4.},{1,1,5.},{2,2,6.},{1,0,.7},{0,2,-.3}});
    check(a);REQUIRE(solver.patternAnalysisCount()==2);
    a.coeffRef(2,1)=.2;REQUIRE_FALSE(a.isCompressed());check(a);
    REQUIRE(solver.patternAnalysisCount()==3);
    auto singular=makeSparseMatrix(3,3,{{0,0,1.},{2,2,1.}});
    REQUIRE_THROWS_WITH(solver.solve(singular,exact),Catch::Matchers::ContainsSubstring("zero_rows=1"));
    check(a);REQUIRE(solver.patternAnalysisCount()==5);
    auto bad=exact;bad[0]=std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_WITH(solver.solve(a,bad),Catch::Matchers::ContainsSubstring("nonfinite input"));
    check(a);
    solver.clearPatternCache();check(a);
    auto small=makeSparseMatrix(2,2,{{0,0,2.},{1,1,3.}});
    REQUIRE((solver.solve(small,VectorXd::Ones(2))-(VectorXd(2)<<.5,1./3.).finished()).norm()<1e-12);
#else
    REQUIRE_THROWS_AS(LinearSolver("umfpack"),std::invalid_argument);
#endif
}

TEST_CASE("LinearSolver reuses symbolic analysis for identical sparse pattern", "[linear_solver]")
{
    const SparseMatrixd A1 = makeSparseMatrix(3, 3, {
        {0, 0, 4.0}, {1, 0, -1.0},
        {0, 1, -1.0}, {1, 1, 4.0}, {2, 1, -1.0},
        {1, 2, -1.0}, {2, 2, 4.0},
    });
    const VectorXd xExpected1 = (VectorXd(3) << 1.0, 2.0, 3.0).finished();
    const VectorXd b1 = A1 * xExpected1;

    LinearSolver solver;
    const VectorXd x1 = solver.solve(A1, b1);
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE((A1 * x1 - b1).norm() == Catch::Approx(0.0).margin(1e-12));

    // Same row/column indices, different numerical values. This should reuse
    // the previous analyzePattern result but still perform a fresh factorize.
    const SparseMatrixd A2 = makeSparseMatrix(3, 3, {
        {0, 0, 6.0}, {1, 0, -2.0},
        {0, 1, -2.0}, {1, 1, 7.0}, {2, 1, -1.5},
        {1, 2, -1.5}, {2, 2, 5.0},
    });
    const VectorXd xExpected2 = (VectorXd(3) << -2.0, 0.5, 4.0).finished();
    const VectorXd b2 = A2 * xExpected2;

    const VectorXd x2 = solver.solve(A2, b2);
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE((A2 * x2 - b2).norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("LinearSolver re-analyzes when sparse pattern changes", "[linear_solver]")
{
    const SparseMatrixd diagonal = makeSparseMatrix(3, 3, {
        {0, 0, 2.0}, {1, 1, 3.0}, {2, 2, 4.0},
    });
    const SparseMatrixd coupled = makeSparseMatrix(3, 3, {
        {0, 0, 2.0}, {1, 0, 0.25},
        {0, 1, 0.5}, {1, 1, 3.0},
        {2, 2, 4.0},
    });

    LinearSolver solver;
    const VectorXd b = (VectorXd(3) << 2.0, 6.0, 12.0).finished();

    const VectorXd x1 = solver.solve(diagonal, b);
    REQUIRE(solver.patternAnalysisCount() == 0);
    REQUIRE(solver.lastReducedDimension() == 0);
    REQUIRE((diagonal * x1 - b).norm() == Catch::Approx(0.0).margin(1e-12));

    const VectorXd x2 = solver.solve(coupled, b);
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE(solver.lastReducedDimension() == 2);
    REQUIRE((coupled * x2 - b).norm() == Catch::Approx(0.0).margin(1e-12));

    solver.clearPatternCache();
    const VectorXd x3 = solver.solve(coupled, b);
    REQUIRE(solver.patternAnalysisCount() == 2);
    REQUIRE((coupled * x3 - b).norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("Eigen SparseLU profiles numeric factor fill", "[linear_solver][performance]")
{
    const SparseMatrixd A = makeSparseMatrix(3, 3, {
        {0, 0, 4.0}, {1, 0, -1.0},
        {0, 1, -1.0}, {1, 1, 4.0}, {2, 1, -1.0},
        {1, 2, -1.0}, {2, 2, 4.0},
    });
    const VectorXd expected = (VectorXd(3) << 1.0, 2.0, 3.0).finished();
    PerformanceProfiler profiler({true, "unused.json"});
    LinearSolver solver("sparselu");
    {
        ActivePerformanceProfilerScope active(&profiler);
        const VectorXd solution = solver.solve(A, A * expected);
        REQUIRE((solution - expected).norm() ==
                Catch::Approx(0.0).margin(1e-12));
    }

    const nlohmann::json observations = profiler.toJson().at("observations");
    REQUIRE(observations.at("linear.factor_nonzeros_l").at("last") > 0.0);
    REQUIRE(observations.at("linear.factor_nonzeros_u").at("last") > 0.0);
    REQUIRE(observations.at("linear.factor_fill_ratio").at("last") > 0.0);
}

TEST_CASE("LinearSolver reuses the numeric factorisation for an identical matrix",
          "[linear_solver][performance]")
{
    const SparseMatrixd A = makeSparseMatrix(3, 3, {
        {0, 0, 4.0}, {1, 0, -1.0},
        {0, 1, -1.0}, {1, 1, 4.0}, {2, 1, -1.0},
        {1, 2, -1.0}, {2, 2, 4.0},
    });
    const VectorXd x1 = (VectorXd(3) << 1.0, 2.0, 3.0).finished();
    const VectorXd x2 = (VectorXd(3) << -4.0, 0.25, 7.0).finished();

    SparseMatrixd perturbed = A;
    perturbed.coeffRef(1, 1) = 4.5;
    perturbed.makeCompressed();

    PerformanceProfiler profiler({true, "unused.json"});
    LinearSolver solver;
    {
        ActivePerformanceProfilerScope active(&profiler);
        // Same matrix, different right-hand sides: only one factorisation.
        REQUIRE((solver.solve(A, A * x1) - x1).norm() ==
                Catch::Approx(0.0).margin(1e-12));
        REQUIRE((solver.solve(A, A * x2) - x2).norm() ==
                Catch::Approx(0.0).margin(1e-12));
        // Changed values with the same pattern must refactorise.
        REQUIRE((solver.solve(perturbed, perturbed * x1) - x1).norm() ==
                Catch::Approx(0.0).margin(1e-12));
    }

    const nlohmann::json counters = profiler.toJson().at("counters");
    REQUIRE(counters.at("linear.solve_calls") == 3);
    REQUIRE(counters.at("linear.factorize_calls") == 2);
    REQUIRE(counters.at("linear.factorize_cache_hits") == 1);
    REQUIRE(counters.at("linear.analyze_calls") == 1);
}

TEST_CASE("LinearSolver rejects invalid dimensions", "[linear_solver]")
{
    const SparseMatrixd rectangular = makeSparseMatrix(2, 3, {
        {0, 0, 1.0}, {1, 1, 1.0},
    });
    const VectorXd rhs2 = VectorXd::Ones(2);

    LinearSolver solver;
    REQUIRE_THROWS_AS(solver.solve(rectangular, rhs2), std::invalid_argument);

    const SparseMatrixd square = makeSparseMatrix(2, 2, {
        {0, 0, 2.0}, {1, 1, 3.0},
    });
    const VectorXd rhsWrongSize = VectorXd::Ones(3);

    REQUIRE_THROWS_AS(solver.solve(square, rhsWrongSize), std::invalid_argument);
}

TEST_CASE("LinearSolver accepts uncompressed input and reuses its compressed pattern", "[linear_solver]")
{
    SparseMatrixd A1(3, 3);
    A1.reserve(Eigen::VectorXi::Constant(3, 3));
    A1.insert(0, 0) = 4.0;
    A1.insert(1, 0) = -1.0;
    A1.insert(0, 1) = -1.0;
    A1.insert(1, 1) = 4.0;
    A1.insert(2, 1) = -1.0;
    A1.insert(1, 2) = -1.0;
    A1.insert(2, 2) = 4.0;
    REQUIRE_FALSE(A1.isCompressed());

    const VectorXd xExpected1 = (VectorXd(3) << 1.0, 2.0, 3.0).finished();
    const VectorXd b1 = A1 * xExpected1;

    LinearSolver solver;
    const VectorXd x1 = solver.solve(A1, b1);
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE((A1 * x1 - b1).norm() == Catch::Approx(0.0).margin(1e-12));

    SparseMatrixd A2(3, 3);
    A2.reserve(Eigen::VectorXi::Constant(3, 3));
    A2.insert(0, 0) = 6.0;
    A2.insert(1, 0) = -2.0;
    A2.insert(0, 1) = -2.0;
    A2.insert(1, 1) = 7.0;
    A2.insert(2, 1) = -1.5;
    A2.insert(1, 2) = -1.5;
    A2.insert(2, 2) = 5.0;
    REQUIRE_FALSE(A2.isCompressed());

    const VectorXd xExpected2 = (VectorXd(3) << -2.0, 0.5, 4.0).finished();
    const VectorXd b2 = A2 * xExpected2;

    const VectorXd x2 = solver.solve(A2, b2);
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE((A2 * x2 - b2).norm() == Catch::Approx(0.0).margin(1e-12));
}
TEST_CASE("LinearSolver factorisation failures report sparse matrix diagnostics",
          "[linear_solver][diagnostics]")
{
    const SparseMatrixd singular = makeSparseMatrix(3, 3, {
        {0, 0, 2.0},
        {2, 2, -4.0},
    });
    const VectorXd rhs = VectorXd::Ones(3);

    LinearSolver solver;
    REQUIRE_THROWS_WITH(
        solver.solve(singular, rhs),
        Catch::Matchers::ContainsSubstring("zero_rows=1") &&
            Catch::Matchers::ContainsSubstring("zero_cols=1") &&
            Catch::Matchers::ContainsSubstring("zero_row_indices=[1]") &&
            Catch::Matchers::ContainsSubstring("zero_col_indices=[1]") &&
            Catch::Matchers::ContainsSubstring("nonfinite_entries=0") &&
            Catch::Matchers::ContainsSubstring("diag_min_abs=0"));
}
TEST_CASE("Experimental electrothermal symbolic reuse checks exact indices and always refactorizes", "[linear][electrothermal]") {
    std::string backend="sparselu_colamd";
    SECTION("SparseLU COLAMD"){}
    SECTION("SparseLU AMD"){backend="sparselu_amd";}
#if defined(VELA_HAS_UMFPACK)
    SECTION("UMFPACK"){backend="umfpack";}
#endif
    experimental::ElectrothermalDirectSolver cached(backend),plain(backend);
    VectorXd exact(3);exact<<1.,-2.,3.;
    auto a=makeSparseMatrix(3,3,{{0,0,4.},{1,1,5.},{2,2,6.},{0,1,.1}});
    const auto check=[&](SparseMatrixd matrix){
        const VectorXd rhs=matrix*exact;
        cached.compute(matrix,true);plain.compute(matrix,false);
        REQUIRE(cached.info()==Eigen::Success);REQUIRE(plain.info()==Eigen::Success);
        const VectorXd x=cached.solve(rhs),y=plain.solve(rhs);
        REQUIRE((matrix*x-rhs).norm()<1e-12);
        REQUIRE((x-y).norm()==0.);
    };
    check(a);REQUIRE(cached.analyses()==1);
    a.coeffRef(0,0)=8.;check(a);REQUIRE(cached.analyses()==1);
    // Same dimensions and nonzero count, different edge location.
    a=makeSparseMatrix(3,3,{{0,0,4.},{1,1,5.},{2,2,6.},{1,2,.1}});
    check(a);REQUIRE(cached.analyses()==2);
    // Also exercise uncompressed input and a changed nonzero count.
    a.coeffRef(2,0)=.2;check(a);REQUIRE(cached.analyses()==3);
    check(a);REQUIRE(cached.analyses()==3);
    REQUIRE(cached.factorizations()==5);REQUIRE(plain.analyses()==5);
    SparseMatrixd b=makeSparseMatrix(2,2,{{0,0,2.},{1,1,3.}});
    cached.compute(b,true);REQUIRE(cached.analyses()==4);
    VectorXd rhs(2);rhs<<2.,6.;REQUIRE((b*cached.solve(rhs)-rhs).norm()<1e-12);
    // Destroy/replace caller storage before solving: refinement must use the
    // solver's owned original coefficients, not a borrowed temporary matrix.
    b.setZero();b.resize(200,200);
    const auto kept=cached.solve(rhs);
    REQUIRE(kept[0]==Catch::Approx(1.));REQUIRE(kept[1]==Catch::Approx(2.));
    // Failed numerical factors must not qualify later symbolic reuse.
    b=makeSparseMatrix(2,2,{{0,0,2.},{1,1,3.}});
    b.coeffRef(1,1)=0.;cached.compute(b,true);
    REQUIRE(cached.info()!=Eigen::Success);
    const auto failedAnalyses=cached.analyses();
    b.coeffRef(1,1)=3.;cached.compute(b,true);
    REQUIRE(cached.info()==Eigen::Success);
    REQUIRE(cached.analyses()==failedAnalyses+1);
    REQUIRE((b*cached.solve(rhs)-rhs).norm()<1e-12);

}

TEST_CASE("Electrothermal direct solver rejects unavailable or unknown backends", "[linear][electrothermal]") {
    experimental::ElectrothermalDirectSolver defaultSolver;
#if defined(VELA_HAS_UMFPACK)
    CHECK(defaultSolver.backend()=="umfpack");
#else
    CHECK(defaultSolver.backend()=="sparselu_colamd");
#endif
    REQUIRE_THROWS_AS(experimental::ElectrothermalDirectSolver("unknown"),std::invalid_argument);
#if !defined(VELA_HAS_UMFPACK)
    REQUIRE_THROWS_AS(experimental::ElectrothermalDirectSolver("umfpack"),std::invalid_argument);
#endif
}

TEST_CASE("Experimental electrothermal stall watch rejects only sustained far-from-converged tiny progress", "[newton][electrothermal]") {
    experimental::ElectrothermalStagnationWatch disabled,watch(3);
    for(int i=0;i<10;++i)REQUIRE_FALSE(disabled.update(7.,6.99999,1e-6,true));
    REQUIRE_FALSE(watch.update(7.,6.99999,1e-6,true));
    REQUIRE_FALSE(watch.update(7.,6.99999,1e-6,true));
    // Useful progress resets the window even at a small step.
    REQUIRE_FALSE(watch.update(7.,6.,1e-6,true));
    REQUIRE_FALSE(watch.update(7.,6.99999,1e-6,true));
    REQUIRE_FALSE(watch.update(7.,6.99999,1e-6,true));
    REQUIRE(watch.update(7.,6.99999,1e-6,true));
    for(int i=0;i<10;++i){
        REQUIRE_FALSE(watch.update(1e-10,9.9999e-11,1e-6,true));
        REQUIRE_FALSE(watch.update(7.,6.99999,.1,true));
        REQUIRE_FALSE(watch.update(7.,7.,1e-6,false));
    }
}
TEST_CASE("Electrothermal predictor cannot hide worse carrier closure behind a better global norm", "[electrothermal][predictor]") {
    const std::array<double,5> baseline{10.,4.,3.,2.,1.},floors{1e-9,1.,.01,.01,.01};
    auto trial=baseline;trial[0]=5.;
    REQUIRE(experimental::improvesElectrothermalPredictor(trial,baseline,floors));
    for(int k=1;k<5;++k){auto worse=trial;worse[k]=2.*baseline[k];
        REQUIRE_FALSE(experimental::improvesElectrothermalPredictor(worse,baseline,floors));}
    REQUIRE_FALSE(experimental::improvesElectrothermalPredictor(baseline,baseline,floors));
    auto converged=baseline;converged[1]=.2;trial=converged;trial[0]=5.;trial[1]=.8;
    REQUIRE(experimental::improvesElectrothermalPredictor(trial,converged,floors));
    trial[1]=1.01;REQUIRE_FALSE(experimental::improvesElectrothermalPredictor(trial,converged,floors));
    trial=baseline;trial[0]=5.;trial[4]=std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(experimental::improvesElectrothermalPredictor(trial,baseline,floors));
}
TEST_CASE("Electrothermal growth follows accepted voltage advance and preserves startup and bounds", "[electrothermal][predictor]") {
    const auto next=[](double planned,double actual,int updates){return experimental::electrothermalActualStep(planned,actual,updates,12,.0001,4./3.);};
    REQUIRE(next(.1,0.,0)==Catch::Approx(.15));
    REQUIRE(next(.5,.1,5)==Catch::Approx(.15));
    REQUIRE(next(.5,.1,15)==Catch::Approx(.1));
    REQUIRE(next(.5,.1,21)==Catch::Approx(.05));
    REQUIRE(next(1.,1.,5)==Catch::Approx(4./3.));
    REQUIRE(next(.0001,.00001,21)==Catch::Approx(.0001));
    REQUIRE_THROWS(next(.1,-.1,5));
    REQUIRE_THROWS(next(.1,std::numeric_limits<double>::quiet_NaN(),5));
}
TEST_CASE("Electrothermal local prediction reproduces low degree fields and bounds amplification", "[electrothermal][predictor]") {
    const std::vector<double> bias{0.,.3,.8,1.};const double target=1.05;
    const auto fit=vela::experimental::electrothermalLocalWeights(bias,target);
    CHECK(fit.degree==2);
    double predicted=0.;for(std::size_t i=0;i<bias.size();++i)predicted+=fit.weights[i]*(2.+3.*bias[i]+4.*bias[i]*bias[i]);
    CHECK(predicted==Catch::Approx(2.+3.*target+4.*target*target).epsilon(1e-12));
    CHECK(fit.weights.sum()==Catch::Approx(1.).epsilon(1e-12));
    const auto guarded=vela::experimental::electrothermalLocalWeights({.375,.7125,1.21875,4./3.},2.4723958333333336);
    CHECK(guarded.degree==1);CHECK(guarded.amplification<=8.);
    CHECK_THROWS_AS(vela::experimental::electrothermalLocalWeights({0.,0.,1.},2.),std::invalid_argument);
    CHECK_THROWS_AS(vela::experimental::electrothermalLocalWeights({0.,.1,.2},20.),std::invalid_argument);
}

TEST_CASE("Electrothermal residual mixing solves a spanned linear problem and rejects empty changes", "[electrothermal][ngmres]") {
    VectorXd x(2),scale(2);x<<1.,1.;scale<<.01,10.;
    Eigen::Matrix2d matrix;matrix<<2.,1.,-1.,3.;
    std::vector<VectorXd> states{(VectorXd(2)<<2.,1.).finished(),(VectorXd(2)<<1.,2.).finished()};
    std::vector<VectorXd> residuals;for(const auto& state:states)residuals.push_back(matrix*state);
    auto mixed=experimental::electrothermalResidualMix(x,matrix*x,states,residuals,scale);
    CHECK(mixed.rank==2);CHECK((matrix*(x+mixed.direction)).norm()<1e-12);
    states.push_back(states.front());residuals.push_back(residuals.front());
    mixed=experimental::electrothermalResidualMix(x,matrix*x,states,residuals,scale);
    CHECK(mixed.rank==2);CHECK((matrix*(x+mixed.direction)).norm()<1e-12);
    CHECK_THROWS_AS(experimental::electrothermalResidualMix(x,matrix*x,{x},{matrix*x},scale),std::invalid_argument);
    scale[0]=0.;CHECK_THROWS_AS(experimental::electrothermalResidualMix(x,matrix*x,states,residuals,scale),std::invalid_argument);
    scale.setConstant(1e300);
    CHECK_THROWS_AS(experimental::electrothermalResidualMix(x,VectorXd::Constant(2,1e300),states,residuals,scale),std::invalid_argument);
}

TEST_CASE("Electrothermal local density projection preserves positive states without global damping", "[electrothermal]") {
    using experimental::electrothermalProjectedDensity;
    CHECK(electrothermalProjectedDensity(100.,-2.,1.)==1.);
    CHECK(electrothermalProjectedDensity(100.,.5,1.)==150.);
    CHECK(electrothermalProjectedDensity(100.,-2.,.1)==80.);
    for(Real n:{Real(1e-280),Real(1e-100),Real(1.),Real(1e26)}) {
        CHECK(electrothermalProjectedDensity(n,0.,1.)==n);
        CHECK(electrothermalProjectedDensity(n,-2.,0.)==n);
        CHECK(electrothermalProjectedDensity(n,-2.,1.)>0.);
    }
    CHECK_THROWS(electrothermalProjectedDensity(1.,std::numeric_limits<Real>::infinity(),1.));
    CHECK_THROWS(electrothermalProjectedDensity(0.,1.,1.));
}

TEST_CASE("Electrothermal natural corrector recognizes linear contraction and rejects growth", "[electrothermal]") {
    VectorXd d(2);d<<2.,-1.;
    for(Real alpha:{.1,.5,1.}) {
        const auto trial=experimental::electrothermalNaturalTrial(d,(1.-alpha)*d,alpha);
        CHECK(trial.decreasing);
        CHECK(trial.theta==Catch::Approx(1.-alpha));
    }
    const auto growth=experimental::electrothermalNaturalTrial(d,2.*d,1.);
    CHECK_FALSE(growth.decreasing);CHECK(growth.nextAlpha>0.);CHECK(growth.nextAlpha<=.5);
    const auto projected=experimental::electrothermalNaturalTrial(d,20.*d,.4,true);
    CHECK(projected.nextAlpha==.2);
    CHECK_THROWS(experimental::electrothermalNaturalTrial(VectorXd::Zero(2),d,1.));
}

TEST_CASE("Opt-in MKL PARDISO preserves equations and invalidates cached factors", "[linear_solver][mkl]") {
#if defined(_WIN32) || defined(__linux__)
    if(!std::getenv("VELA_MKL_RT_DLL") && !std::getenv("VELA_MKL_RT_LIBRARY")) {
        REQUIRE_THROWS_AS(LinearSolver("pardiso_mkl"),std::invalid_argument);
        return;
    }
    LinearSolver solver("pardiso_mkl");
    SparseMatrixd a(4,4);
    a.insert(0,0)=1.00000001; // Non-unit affine constraint.
    a.insert(1,0)=2.;a.insert(1,1)=4.;a.insert(1,2)=-1.;
    a.insert(2,0)=-3.;a.insert(2,1)=.5;a.insert(2,2)=6.;a.insert(2,3)=2.;
    a.insert(3,2)=1.;a.insert(3,3)=5.;a.makeCompressed();
    VectorXd truth(4);truth<<.25,-2.,3.,1.;
    auto verify=[&] {
        VectorXd b=a*truth;auto x=solver.solve(a,b);
        REQUIRE((x-truth).lpNorm<Eigen::Infinity>()<1e-12);
        REQUIRE((a*x-b).lpNorm<Eigen::Infinity>()<1e-12);
    };
    verify();truth<<1.,4.,-2.,.3;verify(); // Same factors, different RHS.
    a.coeffRef(1,1)=7.;verify(); // Same pattern, changed numerical values.
    a.coeffRef(0,2)=.2;a.makeCompressed();verify(); // New pattern.
    SparseMatrixd singular(4,4);singular.insert(0,0)=1.;singular.makeCompressed();
    REQUIRE_THROWS(solver.solve(singular,truth));
    verify(); // Recovery must not use stale factors after an error.
#else
    REQUIRE_THROWS_AS(LinearSolver("pardiso_mkl"),std::invalid_argument);
#endif
}

TEST_CASE("Strict elimination preserves affine constraints and changing systems", "[linear_solver][elimination]") {
    std::vector<const char*> backends{"sparselu_eliminated", "sparselu_eliminated_long_double"};
#ifdef VELA_HAS_UMFPACK
    backends.push_back("umfpack_eliminated");
#endif
#if defined(_WIN32) || defined(__linux__)
    if(std::getenv("VELA_MKL_RT_DLL") || std::getenv("VELA_MKL_RT_LIBRARY")) backends.push_back("pardiso_mkl_eliminated");
#endif
    for(const auto* backend : backends) {
        INFO(backend);
        if(std::string(backend)=="sparselu_eliminated_long_double" &&
           std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits) continue;
        LinearSolver solver(backend);
        SparseMatrixd a(4,4);
        a.insert(0,0)=1.00000001; a.insert(0,2)=0.;
        a.insert(1,0)=7.;a.insert(1,1)=4.;a.insert(1,2)=-2.;a.insert(1,3)=3.;
        a.insert(2,0)=-5.;a.insert(2,1)=1.;a.insert(2,2)=3.;
        a.insert(3,3)=-2.;a.makeCompressed();
        VectorXd exact(4);exact<<.25,2.,-3.,-.5;
        auto verify=[&]() {
            VectorXd b=a*exact, x=solver.solve(a,b);
            REQUIRE((x-exact).lpNorm<Eigen::Infinity>()<1e-12);
            REQUIRE((a*x-b).lpNorm<Eigen::Infinity>()<1e-12);
        };
        verify();exact<<0.,-1.,4.,2.;verify();
        a.coeffRef(0,1)=.1;a.makeCompressed();verify(); // changed constraint set
        SparseMatrixd diagonal(4,4);diagonal.setIdentity();
        REQUIRE((solver.solve(diagonal,exact)-exact).norm()==0.);
        SparseMatrixd zero(4,4);zero.makeCompressed();
        REQUIRE_THROWS(solver.solve(zero,exact));
    }
}

TEST_CASE("Carrier-only elimination retains affine potential constraints", "[linear_solver][elimination]") {
    std::vector<const char*> backends{"sparselu_eliminated_keep_psi"};
#if defined(_WIN32) || defined(__linux__)
    if(std::getenv("VELA_MKL_RT_DLL") || std::getenv("VELA_MKL_RT_LIBRARY")) backends.push_back("pardiso_mkl_eliminated_keep_psi");
#endif
    for(const auto* backend:backends) {
        INFO(backend);
        LinearSolver solver(backend);
        // Two nodes per block. Nonzero boundary increments must affect the
        // coupled free RHS; psi row zero must remain in the factorization.
        SparseMatrixd a(6,6);
        a.insert(0,0)=2.;
        a.insert(1,0)=3.;a.insert(1,1)=8.;a.insert(1,2)=-2.;a.insert(1,4)=1.;
        a.insert(2,2)=1.00000001;
        a.insert(3,0)=-1.;a.insert(3,1)=.5;a.insert(3,2)=2.;a.insert(3,3)=5.;
        a.insert(4,4)=-2.;
        a.insert(5,1)=2.;a.insert(5,4)=3.;a.insert(5,5)=7.;a.makeCompressed();
        VectorXd exact(6);exact<<.25,-2.,.7,3.,-.4,2.;
        auto check=[&](int eliminated) {
            PerformanceProfiler profiler({true,"unused.json"});
            VectorXd b=a*exact,x;
            {ActivePerformanceProfilerScope active(&profiler);x=solver.solve(a,b);}
            REQUIRE((x-exact).lpNorm<Eigen::Infinity>()<1e-12);
            REQUIRE((a*x-b).lpNorm<Eigen::Infinity>()<1e-12);
            const auto obs=profiler.toJson().at("observations");
            REQUIRE(obs.at("linear.eliminated_dofs").at("last")==eliminated);
            REQUIRE(obs.at("linear.reduced_dimension").at("last")==6-eliminated);
            REQUIRE(obs.at("linear.retained_prefix").at("last")==2);
        };
        check(2);exact<<-1.,.1,2.,-3.,.5,4.;check(2);
        a.coeffRef(2,1)=.25;a.makeCompressed();check(1);
        SparseMatrixd incompatible(4,4);incompatible.setIdentity();
        REQUIRE_THROWS_AS(solver.solve(incompatible,VectorXd::Ones(4)),std::invalid_argument);
    }
}

TEST_CASE("Vendor elimination reuses only identical reduced factors", "[linear_solver][elimination]") {
    std::vector<const char*> backends;
#ifdef VELA_HAS_UMFPACK
    backends.push_back("umfpack_eliminated");
#endif
    if(std::getenv("VELA_MKL_RT_DLL") || std::getenv("VELA_MKL_RT_LIBRARY"))
        backends.push_back("pardiso_mkl_eliminated_keep_psi");
    for(const auto* backend:backends) {
        INFO(backend);LinearSolver solver(backend);
        SparseMatrixd a(6,6);a.setIdentity();a.coeffRef(1,0)=2.;a.coeffRef(1,1)=4.;
        a.coeffRef(1,2)=-1.;a.coeffRef(3,1)=.5;a.coeffRef(3,3)=3.;a.makeCompressed();
        VectorXd exact(6);exact<<.2,1.,-.3,2.,.7,-1.;
        PerformanceProfiler profiler({true,"unused.json"});ActivePerformanceProfilerScope scope(&profiler);
        auto check=[&]{const VectorXd b=a*exact;const auto x=solver.solve(a,b);
            REQUIRE((x-exact).lpNorm<Eigen::Infinity>()<1e-12);};
        check();exact*=2.;check(); // Changed affine constraints, same reduced A.
        CHECK(profiler.toJson().at("counters").at("linear.factorize_calls")==1);
        CHECK(profiler.toJson().at("counters").at("linear.factorize_cache_hits")==1);
        a.coeffRef(1,1)=5.;check();
        CHECK(profiler.toJson().at("counters").at("linear.factorize_calls")==2);
        solver.clearNumericCache();check();
        CHECK(profiler.toJson().at("counters").at("linear.factorize_calls")==3);
        solver.clearPatternCache();check();
        CHECK(profiler.toJson().at("counters").at("linear.factorize_calls")==4);
    }
}

TEST_CASE("Extended SparseLU solves nonsymmetric systems without stale factors", "[linear_solver][long_double]") {
    SparseMatrixd a(3,3);
    a.insert(0,0)=1e-8; a.insert(0,1)=2e-8;
    a.insert(1,0)=-3.; a.insert(1,1)=7.; a.insert(1,2)=1.;
    a.insert(2,1)=2e8; a.insert(2,2)=9e8;
    a.makeCompressed();
    VectorXd exact(3); exact << 1.25,-2.5,3.75;
    LinearSolver solver("sparselu_long_double");
    if constexpr (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits) {
        REQUIRE_THROWS_WITH(solver.solve(a,a*exact), Catch::Matchers::ContainsSubstring("wider"));
    } else {
        REQUIRE((solver.solve(a,a*exact)-exact).norm()<1e-12);
        exact << -4.,2.,1.;
        REQUIRE((solver.solve(a,a*exact)-exact).norm()<1e-12);
        a.coeffRef(1,1)=11.;
        REQUIRE((solver.solve(a,a*exact)-exact).norm()<1e-12);
        SparseMatrixd singular(3,3); singular.insert(0,0)=1.; singular.makeCompressed();
        REQUIRE_THROWS(solver.solve(singular,exact));
    }
}

TEST_CASE("SLOT memory mixed linear guard rejects constrained and material errors", "[slot_population]") {
    vela::SparseMatrixd a(3,3);
    a.insert(0,0)=1.;a.insert(1,1)=1e-20;a.insert(1,2)=2e-20;a.insert(2,2)=2.;
    vela::VectorXd x(3);x<<0.,1e-15,0.;
    vela::VectorXd b=vela::VectorXd::Zero(3);
    REQUIRE(vela::slot_population::linearQuality(a,b,x,.02585,1).accepted==false);
    // Normwise guards remain independent of the mixed near-zero allowance.
    b=a*x;
    REQUIRE(vela::slot_population::linearQuality(a,b,x,.02585,1).accepted);
    x[0]=1e-10;
    REQUIRE_FALSE(vela::slot_population::linearQuality(a,b,x,.02585,1).accepted);
    x[0]=0.;x[1]=1.;
    REQUIRE_FALSE(vela::slot_population::linearQuality(a,b,x,.02585,1).accepted);
}

TEST_CASE("SLOT bounded row refinement repairs weak rows without relaxing guards", "[slot_population]") {
    using namespace vela::slot_population;
    SparseMatrixd a(9,9);std::vector<Eigen::Triplet<double>> t{
        {0,0,1.},{1,1,1.},{2,2,1.},{3,3,1.},{4,4,1.},{5,5,1.},
        {6,6,1e-17},{6,7,-.99e-17},{7,6,-1e-19},{7,7,1e-15},{8,8,1.}};
    a.setFromTriplets(t.begin(),t.end());
    VectorXd b=VectorXd::Zero(9),x=VectorXd::Zero(9);b[8]=x[8]=1.;x[6]=1e-11;
    REQUIRE_FALSE(linearQuality(a,b,x,.02585,3).accepted);
    auto r=refineMixedRows(a,b,x,.02585,3,4);
    REQUIRE(r.after.accepted);REQUIRE(linearQuality(a,b,x,.02585,3).accepted);
    CHECK(x.head(6).isZero());CHECK((a*x-b).norm()<1e-27);
    const auto good=x;auto identity=refineMixedRows(a,b,x,.02585,3,4);
    CHECK(identity.sweeps==0);CHECK((x.array()==good.array()).all());
    // Boundary errors cannot be hidden by the weak-row repair.
    x=VectorXd::Zero(9);x[8]=1.;x[0]=1e-3;x[6]=1e-11;auto original=x;
    auto rejected=refineMixedRows(a,b,x,.02585,3,4);
    REQUIRE_FALSE(rejected.after.accepted);CHECK((x.array()==original.array()).all());
    CHECK_THROWS(refineMixedRows(a,b,x,.02585,3,5));
    x=VectorXd::Zero(9);x[8]=1.;x[6]=1e-11;
    LinearSolver factorSolver("sparselu");factorSolver.solve(a,b);
    int calls=0;
    auto factors=refineWithFactors(a,b,x,.02585,3,3,[&](const VectorXd& rhs){++calls;return factorSolver.solve(a,rhs);});
    REQUIRE(factors.after.accepted);CHECK(calls==1);
    CHECK(factorSolver.patternAnalysisCount()==1);
    x[6]=1e-11;original=x;
    auto refused=refineWithFactors(a,b,x,.02585,3,1,[](const VectorXd& rhs){return VectorXd::Zero(rhs.size());});
    CHECK_FALSE(refused.after.accepted);CHECK((x.array()==original.array()).all());
}

TEST_CASE("SLOT diagnostic hold preserves population while potential and neighbors respond", "[slot_population]") {
    using namespace vela::slot_population;
    SparseMatrixd a(6,6);std::vector<Eigen::Triplet<Real>> t{
        {0,0,2.},{0,4,1.},{1,1,3.},{2,2,1.},{3,3,1.},
        {4,4,-2.},{4,5,1.},{5,5,4.},{5,4,-1.}};
    a.setFromTriplets(t.begin(),t.end());
    VectorXd b(6);b<<3.,6.,0.,0.,-20.,7.;
    auto constrained=populationConstraintMatrix(a,{4},2);auto rhs=b;rhs[4]=0.;
    const auto x=LinearSolver("sparselu").solve(constrained,rhs);
    CHECK(x[4]==Catch::Approx(x[0]).margin(1e-14));
    CHECK(std::abs(x[0])>.1);CHECK(std::abs(x[5])>.1);
    const auto originalResidual=(a*x-b).eval();
    for(int i:{0,1,2,3,5})CHECK(std::abs(originalResidual[i])<1e-12);
    CHECK(std::abs(originalResidual[4])>1.); // Held row is not claimed solved.
    CHECK((populationConstraintMatrix(a,{},2)-a).norm()==0.);
    CHECK(retainPopulationHold(true,-2.));CHECK(retainPopulationHold(true,0.));
    CHECK_FALSE(retainPopulationHold(true,.01));CHECK_FALSE(retainPopulationHold(false,-2.));
    CHECK_THROWS(populationConstraintMatrix(a,{0},2));
    CHECK_THROWS(populationConstraintMatrix(a,{4,4},2));
    // Relative Fermi recovery with dq=applied psi exactly preserves eta.
    const double dpsi=.0259*std::log1p(2./.0259);
    const double dq=dpsi+relativeInverse(1e-20,0.);
    CHECK(dq-dpsi==0.);
}

TEST_CASE("Extended SparseLU retains assembled coefficients and direction during population writeback",
          "[linear_solver][extended80]") {
    using namespace vela;
    if (std::numeric_limits<long double>::digits < 64) {
        SKIP("80-bit long double is unavailable");
    }
    const long double low = std::ldexp(1.0L, -60);
    ExtendedSparseMatrix matrix(3, 3);
    matrix.insert(0, 0) = 1 + low;
    matrix.insert(0, 1) = -1;
    matrix.insert(1, 1) = 1;
    matrix.insert(2, 2) = 1;
    ExtendedVector rhs(3);
    rhs << 0, low, -low;
    const auto direction = solveExtendedSparseLu(matrix, rhs);
    matrix.makeCompressed();
    CHECK(solveExtendedSparseLu(matrix, rhs) == direction);
    CHECK(direction[0] < low);
    CHECK(direction[1] == low);
    CHECK(direction[2] == -low);
    const ExtendedVector negativeRhs = -rhs;
    const auto defect = linear_refinement::computeHighPrecisionDefect(matrix, negativeRhs, direction);
    const ExtendedVector narrowed = direction.cast<double>().cast<long double>();
    const auto rounded = linear_refinement::computeHighPrecisionDefect(matrix, negativeRhs, narrowed);
    CHECK(defect.backwardMax < rounded.backwardMax);
    const ExtendedSparseLu factors(matrix);
    const ExtendedVector refined = linear_refinement::refineWithSolve(
        matrix, negativeRhs, ExtendedVector::Ones(3), ExtendedVector::Zero(3),
        [&](const ExtendedSparseMatrix&, const ExtendedVector& correction) {
            return factors.solve(correction);
        }, 1, [](int, const ExtendedVector&, const auto&) {});
    CHECK(refined == direction);
    const ExtendedVector parent = ExtendedVector::Ones(3);
    const auto trial = makePopulationTrial(parent, direction, 1.0L, {true}, {true}, {});
    CHECK(trial.state[0] - parent[0] == low);
    CHECK(trial.state[1] - parent[1] == low);
    CHECK(trial.state[2] - parent[2] == -low);
    CHECK_THROWS(solveExtendedSparseLu(matrix, ExtendedVector::Ones(2)));
    matrix.coeffRef(0, 0) = std::numeric_limits<long double>::infinity();
    CHECK_THROWS(solveExtendedSparseLu(matrix, rhs));
}

TEST_CASE("Wide affine elimination preserves contact and circuit corrections through recovery",
          "[linear_solver][constraints][extended80]") {
    if (std::numeric_limits<long double>::digits < 64) {
        SKIP("80-bit long double is unavailable");
    }
    const long double low = std::ldexp(1.0L, -60);
    ExtendedSparseMatrix matrix(3, 3);
    matrix.insert(0, 0) = 1;
    matrix.insert(0, 2) = -1;
    matrix.insert(1, 0) = -2;
    matrix.insert(1, 1) = 1;
    matrix.insert(2, 0) = 1;
    matrix.insert(2, 2) = 1;
    ExtendedVector rhs(3);
    rhs << 1 + low, -2, 1;
    const auto reduced = reduceLinearConstraints(matrix, rhs, {0, 1});
    REQUIRE(reduced.kept == std::vector<int>{2});
    CHECK(reduced.rhs[0] == -low);
    const auto direction = reduced.recover(solveExtendedSparseLu(reduced.matrix, reduced.rhs));
    CHECK(direction[2] == -low / 2);
    CHECK(direction[1] == low);
    CHECK(direction[0] - 1 == low / 2);
    CHECK(((matrix * direction).eval() == rhs));
}

TEST_CASE("Retained independent contact rows restore the exact correction",
          "[linear_solver][constraints][extended80]") {
    if (std::numeric_limits<long double>::digits < 64)
        SKIP("80-bit long double is unavailable");
    constexpr auto low = 0x1p-60L;
    ExtendedSparseMatrix matrix(3, 3);
    matrix.insert(0, 0) = 1;
    matrix.insert(0, 1) = 0;
    matrix.insert(1, 1) = 2;
    matrix.insert(2, 0) = 1;
    matrix.insert(2, 1) = 1;
    matrix.insert(2, 2) = 1;
    ExtendedVector rhs(3), direction(3);
    rhs << 0, 2 + 2 * low, 2 + low;
    direction << 0x1p-149L, 1 + low, 1;
    const ExtendedVector residual = -rhs;
    const auto raw = linear_refinement::computeHighPrecisionDefect(matrix, residual, direction);
    CHECK(raw.backwardMax == 1.L);
    CHECK(std::abs(raw.values[0]) < 1e-40L);
    const auto original = direction;
    const std::array<int, 2> coupled{0, 2};
    CHECK_THROWS(restoreIndependentConstraintRows(matrix, rhs, coupled, direction));
    CHECK(direction == original);
    const std::array<int, 2> fixed{0, 1};
    restoreIndependentConstraintRows(matrix, rhs, fixed, direction);
    CHECK(direction[0] == 0);
    CHECK(direction[1] == 1 + low);
    CHECK(direction[2] == original[2]);
    CHECK(matrix.rows() == 3);
    CHECK(linear_refinement::computeHighPrecisionDefect(matrix, residual, direction).backwardMax == 0.L);
}
