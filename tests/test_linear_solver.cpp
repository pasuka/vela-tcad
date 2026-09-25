#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PerformanceProfiler.h"
#include "vela/solver/LinearSolver.h"
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
    VectorXd exact(3);exact<<1.,-2.,3.;const VectorXd b=a*exact;
    LinearSolver solver("sparselu");REQUIRE(openblas_get_num_threads()==1);
    PerformanceProfiler profiler({true,"unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        REQUIRE((solver.solve(a,b)-exact).norm()<1e-14);
        REQUIRE((solver.solve(a,2.*b)-2.*exact).norm()<1e-14);
    }
    const auto observation=profiler.toJson().at("observations").at("linear.openblas_threads");
    REQUIRE(observation.at("min")==1);REQUIRE(observation.at("max")==1);
    REQUIRE(observation.at("count")==2);
    openblas_set_num_threads(2);REQUIRE_THROWS(solver.solve(a,b));
    guard.set(nullptr);LinearSolver noControl("sparselu");REQUIRE(openblas_get_num_threads()==2);
    for(const char* count : {"2", "4"}) {
        guard.set(count);
#if defined(VELA_HAS_UMFPACK)
        LinearSolver threaded("umfpack");
#else
        LinearSolver threaded("sparselu");
#endif
        const int expected=std::atoi(count);
        REQUIRE(openblas_get_num_threads()==expected);
        PerformanceProfiler measured({true,"unused.json"});
        {
            ActivePerformanceProfilerScope active(&measured);
            REQUIRE((threaded.solve(a,b)-exact).norm()<1e-13);
            REQUIRE((threaded.solve(a,2.*b)-2.*exact).norm()<1e-13);
        }
        const auto seen=measured.toJson().at("observations").at("linear.openblas_threads");
        REQUIRE(seen.at("min")==expected);REQUIRE(seen.at("max")==expected);
        REQUIRE(seen.at("count")==2);
        openblas_set_num_threads(1);REQUIRE_THROWS(threaded.solve(a,b));
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
    REQUIRE(solver.patternAnalysisCount() == 1);
    REQUIRE((diagonal * x1 - b).norm() == Catch::Approx(0.0).margin(1e-12));

    const VectorXd x2 = solver.solve(coupled, b);
    REQUIRE(solver.patternAnalysisCount() == 2);
    REQUIRE((coupled * x2 - b).norm() == Catch::Approx(0.0).margin(1e-12));

    solver.clearPatternCache();
    const VectorXd x3 = solver.solve(coupled, b);
    REQUIRE(solver.patternAnalysisCount() == 3);
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
