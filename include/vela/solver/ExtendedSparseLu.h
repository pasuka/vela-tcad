#pragma once

#include "vela/core/Types.h"
#include <memory>

namespace vela {

using ExtendedVector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using ExtendedSparseMatrix = Eigen::SparseMatrix<long double>;

/// Reuse factors when corrections use the same prepared matrix.
class ExtendedSparseLu {
public:
    explicit ExtendedSparseLu(const ExtendedSparseMatrix& matrix);
    ~ExtendedSparseLu();
    ExtendedVector solve(const ExtendedVector& rhs) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// Factor and solve an already prepared system without converting its values or
/// the returned direction to binary64. The caller owns scaling and constraints.
ExtendedVector solveExtendedSparseLu(const ExtendedSparseMatrix& matrix,
                                    const ExtendedVector& rhs);

} // namespace vela
