#include "vela/solver/LinearSolver.h"
#include "vela/core/PerformanceProfiler.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <limits>

using namespace vela;
using Catch::Matchers::ContainsSubstring;
namespace {
struct Environment {
    std::string name, old;
    bool existed;
    explicit Environment(const char* key) : name(key), existed(std::getenv(key) != nullptr) {
        if (existed) old = std::getenv(key);
    }
    void set(const char* value) const {
#ifdef _WIN32
        _putenv_s(name.c_str(), value ? value : "");
#else
        if (value) setenv(name.c_str(), value, 1); else unsetenv(name.c_str());
#endif
    }
    ~Environment() { set(existed ? old.c_str() : nullptr); }
};
}

#ifndef VELA_HAS_MKL_PARDISO
TEST_CASE("PARDISO disabled builds reject explicit selection", "[linear_solver][pardiso]") {
    const auto defaults = LinearSolver::availableDirectBackends();
    REQUIRE(std::find(defaults.begin(), defaults.end(), "pardiso_mkl") == defaults.end());
    Environment library("VELA_MKL_RT_LIBRARY"), selection("VELA_LINEAR_SOLVER");
    library.set("disabled-backends-must-not-load-this-library");
    REQUIRE_THROWS_WITH(LinearSolver("pardiso_mkl"), ContainsSubstring("VELA_ENABLE_MKL_PARDISO=ON"));
    selection.set("pardiso_mkl");
    REQUIRE_THROWS_WITH(LinearSolver(), ContainsSubstring("VELA_ENABLE_MKL_PARDISO=ON"));
    LinearSolver solver("sparselu");
    SparseMatrixd a(2,2); a.insert(0,0)=2.; a.insert(1,0)=1.; a.insert(1,1)=3.; a.makeCompressed();
    VectorXd b(2); b<<4.,-1.;
    REQUIRE((a*solver.solve(a,b)-b).norm()<1e-13);
}
#else
TEST_CASE("PARDISO requires explicit runtime selection", "[linear_solver][pardiso]") {
    Environment library("VELA_MKL_RT_LIBRARY"), legacy("VELA_MKL_RT_DLL");
    library.set(nullptr); legacy.set(nullptr);
    const auto defaults = LinearSolver::availableDirectBackends();
    REQUIRE(std::find(defaults.begin(), defaults.end(), "pardiso_mkl") == defaults.end());
    REQUIRE_THROWS_AS(LinearSolver("pardiso_mkl"), std::invalid_argument);
#ifdef VELA_TEST_MKL_LIBRARY
    library.set("relative_mkl_library");
    REQUIRE_THROWS_AS(LinearSolver("pardiso_mkl"), std::invalid_argument);
    const auto missing = std::filesystem::path(VELA_TEST_MKL_LIBRARY).parent_path()
        / "missing-vela-mkl-runtime" / "not-a-library";
    library.set(missing.string().c_str());
    REQUIRE_THROWS_WITH(LinearSolver("pardiso_mkl"), ContainsSubstring("runtime load failed"));
    // An absolute, loadable library lacking MKL symbols must also be rejected.
    // The second fixture contains no MKL entry points.
    library.set(VELA_TEST_MKL_NO_SYMBOLS);
    REQUIRE_THROWS_WITH(LinearSolver("pardiso_mkl"), ContainsSubstring("missing symbol"));
    legacy.set(VELA_TEST_MKL_LIBRARY); library.set("");
    REQUIRE_NOTHROW(LinearSolver("pardiso_mkl"));
    library.set("relative_mkl_library"); // Primary variable takes precedence.
    REQUIRE_THROWS_AS(LinearSolver("pardiso_mkl"), std::invalid_argument);
#endif
}

TEST_CASE("PARDISO propagates vendor errors and recovers without stale factors", "[linear_solver][pardiso]") {
#ifdef VELA_TEST_MKL_LIBRARY
    Environment library("VELA_MKL_RT_LIBRARY"), phase("VELA_TEST_MKL_FAIL_PHASE");
    library.set(VELA_TEST_MKL_LIBRARY); phase.set(nullptr);
    LinearSolver solver("pardiso_mkl");
    SparseMatrixd a(2,2); a.insert(0,0)=2.; a.insert(1,1)=4.; a.makeCompressed();
    VectorXd b(2); b<<6.,-8.;
    for (const char* failedPhase : {"11", "22", "33"}) {
        solver.clearPatternCache(); phase.set(failedPhase);
        REQUIRE_THROWS_WITH(solver.solve(a,b), ContainsSubstring(std::string("phase ")+failedPhase+" error -4"));
        REQUIRE(solver.backend()=="pardiso_mkl");
        phase.set(nullptr);
        REQUIRE((a*solver.solve(a,b)-b).norm()<1e-14);
    }
#else
    SKIP("Runtime adapter supported on Windows and Linux only");
#endif
}

TEST_CASE("PARDISO real runtime preserves nonsymmetric equations and cache boundaries", "[linear_solver][pardiso][mkl_runtime]") {
#ifdef VELA_TEST_MKL_LIBRARY
    const char* primary = std::getenv("VELA_MKL_RT_LIBRARY");
    const char* legacy = std::getenv("VELA_MKL_RT_DLL");
    if ((!primary || !*primary) && (!legacy || !*legacy)) SKIP("Set VELA_MKL_RT_LIBRARY to test installed oneMKL");
    LinearSolver solver("pardiso_mkl");
    SparseMatrixd a(4,4);
    a.insert(0,0)=1.00000001;
    a.insert(1,0)=2.; a.insert(1,1)=4.; a.insert(1,2)=-1.;
    a.insert(2,0)=-3.; a.insert(2,1)=.5; a.insert(2,2)=6.; a.insert(2,3)=2.;
    a.insert(3,2)=1.; a.insert(3,3)=5.; a.makeCompressed();
    VectorXd exact(4); exact<<.25,-2.,3.,1.;
    auto check = [&] {
        VectorXd b=a*exact, x=solver.solve(a,b);
        REQUIRE((x-exact).lpNorm<Eigen::Infinity>()<1e-12);
        REQUIRE((a*x-b).lpNorm<Eigen::Infinity>()<1e-12);
    };
    PerformanceProfiler profiler({true,"unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        check(); exact<<1.,4.,-2.,.3; check();
        a.coeffRef(1,1)=7.; check();
        solver.clearNumericCache(); check();
    }
    REQUIRE(solver.patternAnalysisCount()==1);
    const auto counters=profiler.toJson().at("counters");
    REQUIRE(counters.at("linear.factorize_calls")==3);
    REQUIRE(counters.at("linear.factorize_cache_hits")==1);
    a.coeffRef(0,2)=.2; a.makeCompressed(); check();
    REQUIRE(solver.patternAnalysisCount()==2);
    solver.clearPatternCache(); check(); REQUIRE(solver.patternAnalysisCount()==3);
    auto bad=a; bad.coeffRef(1,1)=std::numeric_limits<double>::infinity();
    REQUIRE_THROWS(solver.solve(bad,exact)); check();
    auto badRhs=exact; badRhs[0]=std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS(solver.solve(a,badRhs)); check();
    SparseMatrixd singular(4,4); singular.insert(0,0)=1.; singular.makeCompressed();
    REQUIRE_THROWS(solver.solve(singular,exact)); check();
    REQUIRE_THROWS(solver.solve(SparseMatrixd(0,0),VectorXd(0))); check();
    REQUIRE(solver.backend()=="pardiso_mkl");
#else
    SKIP("Runtime adapter supported on Windows and Linux only");
#endif
}
#endif
