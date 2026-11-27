#pragma once

#include <Eigen/SparseLU>
#include <vector>

namespace vela::detail {

template<class Scalar, class Ordering = Eigen::COLAMDOrdering<int>>
class SparseLuFactors : public Eigen::SparseLU<Eigen::SparseMatrix<Scalar>, Ordering> {
public:
    bool analysisIsOk() const noexcept { return this->m_analysisIsOk; }

    // Conventional scalar LU work implied by the stored factors.
    double structuralFlops() const {
        using Supernodes = typename SparseLuFactors::SCMatrix;
        using Upper = Eigen::Map<Eigen::SparseMatrix<Scalar>>;
        std::vector<double> lower(this->cols(), 0.), upper(this->cols(), 0.);
        for (Eigen::Index col = 0; col < this->cols(); ++col) {
            for (typename Supernodes::InnerIterator it(this->m_Lstore, col); it; ++it) {
                if (it.row() > col) ++lower[col];
                else if (it.row() < col) ++upper[it.row()];
            }
            for (typename Upper::InnerIterator it(this->m_Ustore, col); it; ++it)
                if (it.row() < col) ++upper[it.row()];
        }
        double flops = 0.;
        for (Eigen::Index i = 0; i < this->cols(); ++i)
            flops += lower[i] * (1. + 2. * upper[i]);
        return flops;
    }
};

} // namespace vela::detail
