#pragma once

#include "vela/core/Types.h"

namespace vela {

/**
 * @brief Sparse direct linear solver based on Eigen SparseLU.
 *
 * Wraps Eigen::SparseLU for the convenience of the rest of the solver
 * pipeline.  Swap this class for an iterative solver (e.g. BiCGSTAB)
 * without changing any call sites.
 */
class LinearSolver {
public:
    /**
     * @brief Solve the linear system A * x = b.
     *
     * @throws std::runtime_error if the factorisation or solve fails.
     */
    VectorXd solve(const SparseMatrixd& A, const VectorXd& b);
};

} // namespace vela
