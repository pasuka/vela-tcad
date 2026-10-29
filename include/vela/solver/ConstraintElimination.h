#pragma once

#include "vela/core/Types.h"
#include <vector>

namespace vela {

/// Rows with a finite nonzero diagonal and no numerical off-diagonal entries.
/// Explicitly stored zeros are ignored; no magnitude threshold is used.
std::vector<int> diagonalConstraintRows(const SparseMatrixd& matrix);

/// Exact affine substitution for prescribed Newton rows A dx = b.
/// The selected constraint block must have an acyclic dependency graph and
/// nonzero diagonal. No physical equation is inferred from a small coefficient.
struct ConstraintReduction {
    SparseMatrixd matrix;
    VectorXd rhs;
    std::vector<int> kept, eliminated;
    SparseMatrixd recoveryMatrix;
    VectorXd recoveryOffset;

    VectorXd recover(const VectorXd& reducedDirection) const;
};

ConstraintReduction eliminateLinearConstraints(
    const SparseMatrixd& matrix, const VectorXd& rhs,
    const std::vector<int>& constraintRows);

} // namespace vela
