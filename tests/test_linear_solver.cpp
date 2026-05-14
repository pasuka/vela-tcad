#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/solver/LinearSolver.h"

#include <Eigen/Sparse>
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
