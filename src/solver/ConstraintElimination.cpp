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

VectorXd ConstraintReduction::recover(const VectorXd& direction) const {
    if (direction.size() != matrix.cols() || !direction.allFinite())
        throw std::invalid_argument("invalid reduced direction");
    VectorXd full = recoveryOffset + recoveryMatrix * direction;
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

ConstraintReduction eliminateLinearConstraints(
    const SparseMatrixd& a, const VectorXd& b, const std::vector<int>& selected) {
    if (a.rows() != a.cols() || a.rows() != b.size() || !b.allFinite())
        throw std::invalid_argument("square finite constraint system required");
    const int n = static_cast<int>(a.rows());
    ConstraintReduction out;
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
    Eigen::SparseMatrix<double,Eigen::RowMajor> rows = a;
    std::vector<int> dependencies(n,0);
    std::vector<std::vector<int>> dependents(n);
    std::queue<int> ready;
    for (int i : out.eliminated) {
        if (rows.coeff(i,i) == 0. || !std::isfinite(rows.coeff(i,i)))
            throw std::invalid_argument("constraint requires a finite nonzero pivot");
        for (decltype(rows)::InnerIterator v(rows,i); v; ++v) {
            if (!std::isfinite(v.value())) throw std::invalid_argument("nonfinite constraint coefficient");
            if (v.col() != i && slot[v.col()] >= 0 && v.value() != 0.) {
                ++dependencies[i]; dependents[v.col()].push_back(i);
            }
        }
        if (!dependencies[i]) ready.push(i);
    }
    // x_E = offset + T x_K. Dependency order also handles a contact population
    // constraint coupled to a contact potential and to the retained circuit.
    out.recoveryOffset = VectorXd::Zero(n);
    std::vector<std::map<int,double>> expressions(n);
    std::size_t visited = 0;
    while (!ready.empty()) {
        const int i = ready.front(); ready.pop(); ++visited;
        const double pivot = rows.coeff(i,i);
        double offset = b[i] / pivot;
        for (decltype(rows)::InnerIterator v(rows,i); v; ++v) {
            if (v.col() == i || v.value() == 0.) continue;
            const double scale = -v.value() / pivot;
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
    std::vector<Eigen::Triplet<double>> entries, recovery;
    entries.reserve(a.nonZeros()); recovery.reserve(n + selected.size());
    for (int i : out.kept) {
        const int r = freeSlot[i]; out.rhs[r] = b[i];
        recovery.emplace_back(i,r,1.);
        for (decltype(rows)::InnerIterator v(rows,i); v; ++v) {
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
    if (!out.rhs.allFinite() || !Eigen::Map<const VectorXd>(out.matrix.valuePtr(),out.matrix.nonZeros()).allFinite() ||
        !Eigen::Map<const VectorXd>(out.recoveryMatrix.valuePtr(),out.recoveryMatrix.nonZeros()).allFinite())
        throw std::runtime_error("nonfinite constraint substitution");
    return out;
}
} // namespace vela
