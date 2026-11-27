#pragma once

#include "vela/numerics/ScaledProduct.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace vela::detail {
template<class Scalar>
struct BorderedLoadLineSystem {
    using Matrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Matrix rawAugmented, rowScaledAugmented, augmented;
    Vector rawRhs, rhs, rowScale, columnScale;
};

/// Build both scalings in the same arithmetic as the original bordered system.
/// A caller may use the row-only system for a rank-revealing factorization.
template<class Scalar>
std::optional<BorderedLoadLineSystem<Scalar>> buildBorderedLoadLineSystem(
    const Eigen::SparseMatrix<Scalar>& deviceJacobian,
    const std::type_identity_t<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>& f,
    const std::type_identity_t<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>& fLambda,
    const std::type_identity_t<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>& circuitStateDerivative,
    std::type_identity_t<Scalar> circuitVoltageDerivative, std::type_identity_t<Scalar> loadResidual,
    const std::type_identity_t<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>& physicalColumns,
    std::vector<int>& prescribed, std::string& failureReason)
{
    using Matrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    const auto n = deviceJacobian.rows();
    if (n <= 0 || deviceJacobian.cols() != n || f.size() != n || fLambda.size() != n ||
        circuitStateDerivative.size() != n || physicalColumns.size() != n + 1) {
        failureReason = "direct bordered dimensions are invalid";
        return std::nullopt;
    }
    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>(
        deviceJacobian.nonZeros() + 2 * n + 1));
    for (int outer = 0; outer < deviceJacobian.outerSize(); ++outer) {
        for (typename Matrix::InnerIterator entry(deviceJacobian, outer);
             entry; ++entry) {
            triplets.emplace_back(entry.row(), entry.col(), entry.value());
        }
    }
    bool circuitDependsOnState = false;
    for (Eigen::Index row = 0; row < n; ++row) {
        if (fLambda(row) != 0.0)
            triplets.emplace_back(row, n, fLambda(row));
        const Scalar value = circuitStateDerivative(row);
        if (value != 0.0) {
            triplets.emplace_back(n, row, value);
            circuitDependsOnState = true;
        }
    }
    triplets.emplace_back(n, n, circuitVoltageDerivative);
    // A scalar circuit row prescribes the voltage increment exactly. A
    // rounded nonzero increment on a homogeneous row would have unit
    // componentwise backward error regardless of its magnitude.
    const int circuitRow = static_cast<int>(n);
    if (!circuitDependsOnState && circuitVoltageDerivative != 0.0 &&
        std::find(prescribed.begin(), prescribed.end(), circuitRow) == prescribed.end())
        prescribed.push_back(circuitRow);

    Matrix rawAugmented(n + 1, n + 1);
    rawAugmented.setFromTriplets(triplets.begin(), triplets.end());
    rawAugmented.makeCompressed();
    Vector rawRhs(n + 1);
    rawRhs.head(n) = -f;
    rawRhs(n) = -loadResidual;

    // Two-sided max-norm equilibration is important for a TCAD device
    // bordered by a 1e12-ohm-um circuit row.  Row scaling alone leaves the
    // voltage column and weak contact-current columns many orders apart.
    // Include the right-hand side in each row norm.  A predictor can place
    // a carrier row where every derivative is almost underflow-small while
    // its nonlinear residual is still finite; Jacobian-only scaling would
    // then amplify the RHS toward the scalar limit without changing the equation.
    Vector rowMatrixMagnitude = Vector::Zero(n + 1);
    std::vector<bool> rawRowSupport(n+1,false),rawColumnSupport(n+1,false);
    for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
        for (typename Matrix::InnerIterator entry(rawAugmented, outer);
             entry; ++entry) {
            if (!std::isfinite(entry.value())) {
                failureReason = "direct bordered matrix contains a non-finite entry: row="+std::to_string(entry.row())+" col="+std::to_string(entry.col());
                return std::nullopt;
            }
            if(entry.value()!=0.){rawRowSupport[entry.row()]=true;rawColumnSupport[entry.col()]=true;}
            rowMatrixMagnitude(entry.row()) = std::max(
                rowMatrixMagnitude(entry.row()), std::abs(entry.value()*physicalColumns[entry.col()]));
        }
    }
    Vector rowScale(n + 1);
    for (Eigen::Index row = 0; row <= n; ++row) {
        if (!(rowMatrixMagnitude(row) > 0.0) ||
            !std::isfinite(rowMatrixMagnitude(row))) {
            failureReason = std::string(std::isfinite(rowMatrixMagnitude(row))
                ? (rawRowSupport[row]?"direct bordered matrix has weighted row magnitude underflow row ":"direct bordered matrix has structural zero row ")
                : "direct bordered matrix has non-finite weighted row magnitude ") +
                std::to_string(row);
            return std::nullopt;
        }
        const Scalar rowMagnitude = std::max(
            rowMatrixMagnitude(row), std::abs(rawRhs(row)));
        rowScale(row) = 1.0 / rowMagnitude;
    }
    Vector columnScale = Vector::Zero(n + 1);
    for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
        for (typename Matrix::InnerIterator entry(rawAugmented, outer);
             entry; ++entry) {
            columnScale(entry.col()) = std::max(
                columnScale(entry.col()),
                std::abs(scaledProduct(rowScale(entry.row()),entry.value(),physicalColumns[entry.col()])));
        }
    }
    for (Eigen::Index column = 0; column <= n; ++column) {
        if (!(columnScale(column) > 0.0) ||
            !std::isfinite(columnScale(column))) {
            failureReason = std::string(std::isfinite(columnScale(column))
                ? (rawColumnSupport[column]?"direct bordered matrix has weighted column magnitude underflow column ":"direct bordered matrix has structural zero column ")
                : "direct bordered matrix has non-finite column scale ") +
                std::to_string(column);
            return std::nullopt;
        }
        columnScale(column) = physicalColumns[column] / columnScale(column);
    }

    std::vector<Eigen::Triplet<Scalar>> rowScaledTriplets;
    std::vector<Eigen::Triplet<Scalar>> scaledTriplets;
    rowScaledTriplets.reserve(
        static_cast<std::size_t>(rawAugmented.nonZeros()));
    scaledTriplets.reserve(static_cast<std::size_t>(rawAugmented.nonZeros()));
    for (int outer = 0; outer < rawAugmented.outerSize(); ++outer) {
        for (typename Matrix::InnerIterator entry(rawAugmented, outer);
             entry; ++entry) {
            const Scalar rowScaledValue =
                rowScale(entry.row()) * entry.value();
            rowScaledTriplets.emplace_back(
                entry.row(), entry.col(), rowScaledValue);
            scaledTriplets.emplace_back(
                entry.row(), entry.col(),
                scaledProduct(rowScale(entry.row()),entry.value(),columnScale(entry.col())));
        }
    }
    Matrix rowScaledAugmented(n + 1, n + 1);
    rowScaledAugmented.setFromTriplets(
        rowScaledTriplets.begin(), rowScaledTriplets.end());
    rowScaledAugmented.makeCompressed();
    Matrix augmented(n + 1, n + 1);
    augmented.setFromTriplets(
        scaledTriplets.begin(), scaledTriplets.end());
    augmented.makeCompressed();
    Vector rhs = Vector(rawRhs.array() * rowScale.array());
    return BorderedLoadLineSystem<Scalar>{std::move(rawAugmented), std::move(rowScaledAugmented),
        std::move(augmented), std::move(rawRhs), std::move(rhs), std::move(rowScale), std::move(columnScale)};
}
} // namespace vela::detail
