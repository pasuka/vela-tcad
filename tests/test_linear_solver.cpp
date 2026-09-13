#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PerformanceProfiler.h"
#include "vela/solver/LinearSolver.h"
#include "../src/tools/ElectrothermalIterationControl.h"

#include <Eigen/Sparse>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <vector>

using namespace vela;

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

TEST_CASE("LinearSolver profiles numeric factor fill", "[linear_solver][performance]")
{
    const SparseMatrixd A = makeSparseMatrix(3, 3, {
        {0, 0, 4.0}, {1, 0, -1.0},
        {0, 1, -1.0}, {1, 1, 4.0}, {2, 1, -1.0},
        {1, 2, -1.0}, {2, 2, 4.0},
    });
    const VectorXd expected = (VectorXd(3) << 1.0, 2.0, 3.0).finished();
    PerformanceProfiler profiler({true, "unused.json"});
    LinearSolver solver;
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
    experimental::ElectrothermalSparseLU cached,plain;
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
