#include "vela/solver/LinearSolver.h"
#include <Eigen/SparseLU>
#include <stdexcept>

namespace vela {

VectorXd LinearSolver::solve(const SparseMatrixd& A, const VectorXd& b)
{
    Eigen::SparseLU<SparseMatrixd> solver;

    solver.analyzePattern(A);
    solver.factorize(A);

    if (solver.info() != Eigen::Success)
        throw std::runtime_error(
            "LinearSolver: SparseLU factorisation failed. "
            "Matrix may be singular or ill-conditioned.");

    VectorXd x = solver.solve(b);

    if (solver.info() != Eigen::Success)
        throw std::runtime_error(
            "LinearSolver: SparseLU back-substitution failed.");

    return x;
}

} // namespace vela
