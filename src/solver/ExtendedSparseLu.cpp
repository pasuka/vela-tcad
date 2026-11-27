#include "vela/solver/ExtendedSparseLu.h"

#include "vela/solver/SparseLuFactors.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {

struct ExtendedSparseLu::Impl {
    detail::SparseLuFactors<long double> solver;
};

ExtendedSparseLu::ExtendedSparseLu()
    : impl_(std::make_unique<Impl>())
{
    if (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
        throw std::runtime_error("Extended SparseLU requires long double wider than double");
}

ExtendedSparseLu::ExtendedSparseLu(const ExtendedSparseMatrix& matrix)
    : ExtendedSparseLu()
{
    if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
        throw std::invalid_argument("Invalid extended SparseLU matrix dimensions");
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col)
        for (ExtendedSparseMatrix::InnerIterator entry(matrix, col); entry; ++entry)
            if (!std::isfinite(entry.value()))
                throw std::invalid_argument("Nonfinite extended SparseLU coefficient");

    if (matrix.isCompressed()) {
        impl_->solver.compute(matrix);
    } else {
        ExtendedSparseMatrix compressed = matrix;
        compressed.makeCompressed();
        impl_->solver.compute(compressed);
    }
    if (impl_->solver.info() != Eigen::Success)
        throw std::runtime_error("Extended SparseLU factorization failed: " +
                                 impl_->solver.lastErrorMessage());
}

ExtendedSparseLu::~ExtendedSparseLu() = default;

void ExtendedSparseLu::analyzePattern(const ExtendedSparseMatrix& matrix)
{
    impl_->solver.analyzePattern(matrix);
}

void ExtendedSparseLu::factorize(const ExtendedSparseMatrix& matrix)
{
    impl_->solver.factorize(matrix);
    if (impl_->solver.info() != Eigen::Success)
        throw std::runtime_error("Extended SparseLU factorization failed: " +
                                 impl_->solver.lastErrorMessage());
}

ExtendedVector ExtendedSparseLu::solve(const ExtendedVector& rhs) const
{
    if (rhs.size() != impl_->solver.rows() || !rhs.allFinite())
        throw std::invalid_argument("Invalid extended SparseLU RHS");
    ExtendedVector direction = impl_->solver.solve(rhs);
    if (impl_->solver.info() != Eigen::Success || !direction.allFinite())
        throw std::runtime_error("Extended SparseLU solve failed or produced a nonfinite solution");
    return direction;
}

ExtendedVector solveExtendedSparseLu(const ExtendedSparseMatrix& matrix,
                                    const ExtendedVector& rhs)
{
    return ExtendedSparseLu(matrix).solve(rhs);
}

} // namespace vela
