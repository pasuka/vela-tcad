#pragma once
#include "vela/core/Types.h"
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace vela::experimental {
// Diagnostic tool only. Numeric factorization always runs, including when
// matrix values are identical. Only exact compressed index arrays are reused.
class ElectrothermalSparseLU {
public:
    void compute(SparseMatrixd matrix,bool reuse) {
        matrix.makeCompressed();
        const bool same=reuse && lu_ && rows_==matrix.rows() && cols_==matrix.cols()
            && inner_.size()==static_cast<std::size_t>(matrix.nonZeros())
            && std::equal(outer_.begin(),outer_.end(),matrix.outerIndexPtr())
            && std::equal(inner_.begin(),inner_.end(),matrix.innerIndexPtr());
        if(!same){
            lu_=std::make_unique<Eigen::SparseLU<SparseMatrixd>>();
            lu_->analyzePattern(matrix);++analyses_;
            rows_=matrix.rows();cols_=matrix.cols();
            outer_.assign(matrix.outerIndexPtr(),matrix.outerIndexPtr()+matrix.outerSize()+1);
            inner_.assign(matrix.innerIndexPtr(),matrix.innerIndexPtr()+matrix.nonZeros());
        }
        lu_->factorize(matrix);++factorizations_;
    }
    Eigen::ComputationInfo info() const {return lu_->info();}
    VectorXd solve(const VectorXd& rhs) {return lu_->solve(rhs);}
    unsigned analyses() const {return analyses_;}
    unsigned factorizations() const {return factorizations_;}
private:
    std::unique_ptr<Eigen::SparseLU<SparseMatrixd>> lu_;
    Eigen::Index rows_=0,cols_=0;
    std::vector<SparseMatrixd::StorageIndex> outer_,inner_;
    unsigned analyses_=0,factorizations_=0;
};

// An opt-in early rejection, never a convergence criterion. Compare before
// and after under the SAME iteration row scaling, not across recenterings.
class ElectrothermalStagnationWatch {
public:
    explicit ElectrothermalStagnationWatch(unsigned window=0):window_(window){}
    bool update(double before,double after,double alpha,bool accepted) {
        const bool slow=accepted && std::isfinite(before) && std::isfinite(after)
            && std::isfinite(alpha) && before>1e-6 && after>=0. && after<before
            && alpha>0. && alpha<1e-4 && (before-after)/before<1e-3;
        consecutive_=slow?consecutive_+1:0;
        return window_>0 && consecutive_>=window_;
    }
private:
    unsigned window_=0,consecutive_=0;
};
} // namespace vela::experimental
