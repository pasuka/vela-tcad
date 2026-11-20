#pragma once

#include "vela/core/Types.h"
#include <span>
#include <vector>

namespace vela {

/// Rows with a finite nonzero diagonal and no numerical off-diagonal entries.
/// Explicitly stored zeros are ignored; no magnitude threshold is used.
std::vector<int> diagonalConstraintRows(const SparseMatrixd& matrix);

/// Exact affine substitution for prescribed Newton rows A dx = b.
/// The selected constraint block must have an acyclic dependency graph and
/// nonzero diagonal. No physical equation is inferred from a small coefficient.
template<class Scalar>
struct BasicConstraintReduction {
    using Matrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    Matrix matrix;
    Vector rhs;
    std::vector<int> kept, eliminated;
    Matrix recoveryMatrix;
    Vector recoveryOffset;

    Vector recover(const Vector& reducedDirection) const;
};
using ConstraintReduction = BasicConstraintReduction<Real>;

template<class Scalar>
BasicConstraintReduction<Scalar> reduceLinearConstraints(
    const Eigen::SparseMatrix<Scalar>& matrix,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& rhs,
    const std::vector<int>& constraintRows);

/// Restore explicitly prescribed independent rows after unscaling the
/// direction. The factored system keeps these rows and its original size.
/// Every selected row must have a finite nonzero pivot and no numerical
/// off-diagonals. Re-evaluate the defect on the original system afterward.
template<class Scalar>
void restoreIndependentConstraintRows(
    const Eigen::SparseMatrix<Scalar>& matrix,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& rhs,
    std::span<const int> constraintRows,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& direction);

ConstraintReduction eliminateLinearConstraints(
    const SparseMatrixd& matrix, const VectorXd& rhs,
    const std::vector<int>& constraintRows);

} // namespace vela
