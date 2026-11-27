#include "vela/solver/ConstraintElimination.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <stdexcept>
#include <utility>

namespace vela {

std::vector<int> diagonalConstraintRows(const SparseMatrixd& a) {
    if (a.rows() != a.cols()) throw std::invalid_argument("constraints require a square matrix");
    std::vector<bool> off(a.rows(), false);
    VectorXd diagonal = VectorXd::Zero(a.rows());
    for (int k = 0; k < a.outerSize(); ++k)
        for (SparseMatrixd::InnerIterator v(a,k); v; ++v) {
            if (!std::isfinite(v.value())) throw std::invalid_argument("nonfinite constraint matrix");
            if (v.row() == v.col()) diagonal[v.row()] = v.value();
            else if (v.value() != 0.) off[v.row()] = true;
        }
    std::vector<int> rows;
    for (int i = 0; i < a.rows(); ++i)
        if (!off[i] && diagonal[i] != 0.) rows.push_back(i);
    return rows;
}

template<class Scalar>
typename BasicConstraintReduction<Scalar>::Vector
BasicConstraintReduction<Scalar>::recover(const Vector& direction) const {
    if (direction.size() != matrix.cols() || !direction.allFinite())
        throw std::invalid_argument("invalid reduced direction");
    Vector full = recoveryOffset + recoveryMatrix * direction;
    if (!full.allFinite()) throw std::runtime_error("nonfinite constraint recovery");
    return full;
}

ConstraintReduction detail::eliminateDiagonalLinearConstraints(
    const SparseMatrixd& matrix, const VectorXd& rhs, std::vector<int> rows)
{
    const int size = static_cast<int>(matrix.rows());
    const int reducedSize = size - static_cast<int>(rows.size());
    ConstraintReduction out;
    out.eliminated = std::move(rows);
    out.kept.reserve(reducedSize);
    out.recoveryOffset = VectorXd::Zero(size);
    std::vector<int> reducedIndex(size, -1);
    for (int row : out.eliminated)
        out.recoveryOffset[row] = rhs[row] / matrix.coeff(row, row);
    out.rhs.resize(reducedSize);
    auto prescribed = out.eliminated.begin();
    for (int row = 0; row < size; ++row) {
        if (prescribed != out.eliminated.end() && row == *prescribed) {
            ++prescribed;
            continue;
        }
        const int index = static_cast<int>(out.kept.size());
        reducedIndex[row] = index;
        out.kept.push_back(row);
        out.rhs[index] = rhs[row];
    }

    // The retained CSC entries are already sorted. Process prescribed columns
    // in ascending order so the RHS corrections keep their original order.
    out.matrix.resize(reducedSize, reducedSize);
    out.matrix.reserve(matrix.nonZeros());
    for (int col = 0; col < size; ++col) {
        const int targetColumn = reducedIndex[col];
        if (targetColumn >= 0)
            out.matrix.startVec(targetColumn);
        for (SparseMatrixd::InnerIterator entry(matrix, col); entry; ++entry) {
            const int targetRow = reducedIndex[entry.row()];
            if (targetRow < 0)
                continue;
            if (targetColumn >= 0)
                out.matrix.insertBack(targetRow, targetColumn) = entry.value();
            else if (entry.value() != 0.)
                out.rhs[targetRow] -= entry.value() * out.recoveryOffset[col];
        }
    }
    out.matrix.finalize();
    out.recoveryMatrix.resize(size, reducedSize);
    out.recoveryMatrix.reserve(reducedSize);
    for (int col = 0; col < reducedSize; ++col) {
        out.recoveryMatrix.startVec(col);
        out.recoveryMatrix.insertBack(out.kept[col], col) = 1.;
    }
    out.recoveryMatrix.finalize();
    if (!out.recoveryOffset.allFinite() || !out.rhs.allFinite())
        throw std::runtime_error("nonfinite constraint substitution");
    return out;
}

template<class Scalar>
BasicConstraintReduction<Scalar> reduceLinearConstraints(
    const Eigen::SparseMatrix<Scalar>& a,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b,
    const std::vector<int>& selected) {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    if (a.rows() != a.cols() || a.rows() != b.size() || !b.allFinite())
        throw std::invalid_argument("square finite constraint system required");
    const int n = static_cast<int>(a.rows());
    BasicConstraintReduction<Scalar> out;
    out.eliminated = selected;
    std::sort(out.eliminated.begin(), out.eliminated.end());
    if (std::adjacent_find(out.eliminated.begin(), out.eliminated.end()) != out.eliminated.end())
        throw std::invalid_argument("duplicate constraint row");
    std::vector<int> slot(n,-1), freeSlot(n,-1);
    for (int i : out.eliminated) {
        if (i < 0 || i >= n) throw std::invalid_argument("constraint row outside system");
        slot[i] = i;
    }
    for (int i = 0; i < n; ++i) if (slot[i] < 0) {
        freeSlot[i] = static_cast<int>(out.kept.size()); out.kept.push_back(i);
    }
    Eigen::SparseMatrix<Scalar,Eigen::RowMajor> rows = a;
    std::vector<int> dependencies(n,0);
    std::vector<std::vector<int>> dependents(n);
    std::queue<int> ready;
    for (int i : out.eliminated) {
        if (rows.coeff(i,i) == 0. || !std::isfinite(rows.coeff(i,i)))
            throw std::invalid_argument("constraint requires a finite nonzero pivot");
        for (typename decltype(rows)::InnerIterator v(rows,i); v; ++v) {
            if (!std::isfinite(v.value())) throw std::invalid_argument("nonfinite constraint coefficient");
            if (v.col() != i && slot[v.col()] >= 0 && v.value() != 0.) {
                ++dependencies[i]; dependents[v.col()].push_back(i);
            }
        }
        if (!dependencies[i]) ready.push(i);
    }
    // x_E = offset + T x_K. Dependency order also handles a contact population
    // constraint coupled to a contact potential and to the retained circuit.
    out.recoveryOffset = Vector::Zero(n);
    std::vector<std::map<int,Scalar>> expressions(n);
    std::size_t visited = 0;
    while (!ready.empty()) {
        const int i = ready.front(); ready.pop(); ++visited;
        const Scalar pivot = rows.coeff(i,i);
        Scalar offset = b[i] / pivot;
        for (typename decltype(rows)::InnerIterator v(rows,i); v; ++v) {
            if (v.col() == i || v.value() == 0.) continue;
            const Scalar scale = -v.value() / pivot;
            if (slot[v.col()] < 0) expressions[i][freeSlot[v.col()]] += scale;
            else {
                offset += scale * out.recoveryOffset[v.col()];
                for (const auto& [k,value] : expressions[v.col()]) expressions[i][k] += scale * value;
            }
        }
        if (!std::isfinite(offset)) throw std::runtime_error("nonfinite constraint offset");
        out.recoveryOffset[i] = offset;
        for (int j : dependents[i]) if (--dependencies[j] == 0) ready.push(j);
    }
    if (visited != out.eliminated.size())
        throw std::invalid_argument("cyclic prescribed constraints are unsupported");
    const int m = static_cast<int>(out.kept.size());
    out.rhs.resize(m);
    std::vector<Eigen::Triplet<Scalar>> entries, recovery;
    entries.reserve(a.nonZeros()); recovery.reserve(n + selected.size());
    for (int i : out.kept) {
        const int r = freeSlot[i]; out.rhs[r] = b[i];
        recovery.emplace_back(i,r,1.);
        for (typename decltype(rows)::InnerIterator v(rows,i); v; ++v) {
            if (!std::isfinite(v.value())) throw std::invalid_argument("nonfinite free coefficient");
            if (slot[v.col()] < 0) entries.emplace_back(r,freeSlot[v.col()],v.value());
            else if (v.value() != 0.) {
                out.rhs[r] -= v.value() * out.recoveryOffset[v.col()];
                for (const auto& [k,value] : expressions[v.col()])
                    if (value != 0.) entries.emplace_back(r,k,v.value()*value);
            }
        }
    }
    for (int i : out.eliminated) for (const auto& [k,value] : expressions[i])
        if (value != 0.) recovery.emplace_back(i,k,value);
    out.matrix.resize(m,m); out.matrix.setFromTriplets(entries.begin(),entries.end());
    out.recoveryMatrix.resize(n,m); out.recoveryMatrix.setFromTriplets(recovery.begin(),recovery.end());
    if (!out.rhs.allFinite() || !Eigen::Map<const Vector>(out.matrix.valuePtr(),out.matrix.nonZeros()).allFinite() ||
        !Eigen::Map<const Vector>(out.recoveryMatrix.valuePtr(),out.recoveryMatrix.nonZeros()).allFinite())
        throw std::runtime_error("nonfinite constraint substitution");
    return out;
}

ConstraintReduction eliminateLinearConstraints(
    const SparseMatrixd& matrix, const VectorXd& rhs, const std::vector<int>& constraintRows) {
    return reduceLinearConstraints(matrix, rhs, constraintRows);
}

template<class Scalar>
void restoreIndependentConstraintRows(
    const Eigen::SparseMatrix<Scalar>& matrix,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& rhs,
    std::span<const int> constraintRows,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& direction)
{
    if (matrix.rows() != matrix.cols() || rhs.size() != matrix.rows() ||
        direction.size() != rhs.size() || !rhs.allFinite() || !direction.allFinite())
        throw std::invalid_argument("Invalid system for independent constraint restoration");
    if (constraintRows.empty())
        return;
    std::vector<bool> selected(matrix.rows(), false);
    for (const int row : constraintRows) {
        if (row < 0 || row >= matrix.rows() || selected[row])
            throw std::invalid_argument("Invalid independent constraint row");
        selected[row] = true;
    }
    for (int col = 0; col < matrix.outerSize(); ++col)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator entry(matrix, col); entry; ++entry)
            if (selected[entry.row()] && entry.row() != entry.col() && entry.value() != 0)
                throw std::invalid_argument("Independent constraint row contains an off-diagonal entry");
    std::vector<Scalar> restored;
    restored.reserve(constraintRows.size());
    for (int row : constraintRows) {
        const auto diagonal = matrix.coeff(row, row);
        if (diagonal == 0 || !std::isfinite(diagonal))
            throw std::invalid_argument("Independent constraint row has an invalid pivot");
        const auto value = rhs[row] / diagonal;
        if (!std::isfinite(value))
            throw std::runtime_error("Nonfinite independent constraint direction");
        restored.push_back(value);
    }
    for (std::size_t i = 0; i < constraintRows.size(); ++i)
        direction[constraintRows[i]] = restored[i];
}

template struct BasicConstraintReduction<double>;
template struct BasicConstraintReduction<long double>;
template BasicConstraintReduction<double> reduceLinearConstraints(
    const SparseMatrixd&, const VectorXd&, const std::vector<int>&);
template BasicConstraintReduction<long double> reduceLinearConstraints(
    const Eigen::SparseMatrix<long double>&,
    const Eigen::Matrix<long double, Eigen::Dynamic, 1>&, const std::vector<int>&);
template void restoreIndependentConstraintRows<double>(
    const SparseMatrixd&, const VectorXd&, std::span<const int>, VectorXd&);
template void restoreIndependentConstraintRows<long double>(
    const Eigen::SparseMatrix<long double>&, const Eigen::Matrix<long double, Eigen::Dynamic, 1>&,
    std::span<const int>, Eigen::Matrix<long double, Eigen::Dynamic, 1>&);

} // namespace vela
