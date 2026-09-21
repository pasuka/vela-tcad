#pragma once
#include "vela/core/Types.h"
#include <Eigen/SparseLU>
#if defined(VELA_HAS_UMFPACK)
#include <Eigen/UmfPackSupport>
#endif
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela::experimental {
// Numeric factorization always runs, including when
// matrix values are identical. Only exact compressed index arrays are reused.
class ElectrothermalDirectSolver {
public:
    explicit ElectrothermalDirectSolver(std::string backend="sparselu_colamd"):backend_(std::move(backend)) {
        if(backend_!="sparselu_colamd" && backend_!="sparselu_amd" && backend_!="umfpack")
            throw std::invalid_argument("Unknown electrothermal linear solver: "+backend_);
#if !defined(VELA_HAS_UMFPACK)
        if(backend_=="umfpack")throw std::invalid_argument("Electrothermal UMFPACK was not enabled in this build");
#endif
    }
    const std::string& backend() const {return backend_;}
    void compute(SparseMatrixd matrix,bool reuse) {
        matrix.makeCompressed();
        // Eigen's UMFPACK adapter retains a reference to the input matrix for
        // iterative refinement during solve(), beyond this compute() call.
        matrix_=std::move(matrix);
        const bool same=reuse && lu_ && rows_==matrix_.rows() && cols_==matrix_.cols()
            && inner_.size()==static_cast<std::size_t>(matrix_.nonZeros())
            && std::equal(outer_.begin(),outer_.end(),matrix_.outerIndexPtr())
            && std::equal(inner_.begin(),inner_.end(),matrix_.innerIndexPtr());
        if(!same){
            lu_=makeSolver();
            lu_->analyzePattern(matrix_);++analyses_;
            rows_=matrix_.rows();cols_=matrix_.cols();
            outer_.assign(matrix_.outerIndexPtr(),matrix_.outerIndexPtr()+matrix_.outerSize()+1);
            inner_.assign(matrix_.innerIndexPtr(),matrix_.innerIndexPtr()+matrix_.nonZeros());
        }
        lu_->factorize(matrix_);++factorizations_;
    }
    Eigen::ComputationInfo info() const {return lu_->info();}
    VectorXd solve(const VectorXd& rhs) {return lu_->solve(rhs);}
    unsigned analyses() const {return analyses_;}
    unsigned factorizations() const {return factorizations_;}
private:
    struct Solver {
        virtual ~Solver()=default;
        virtual void analyzePattern(const SparseMatrixd&)=0;
        virtual void factorize(const SparseMatrixd&)=0;
        virtual Eigen::ComputationInfo info() const=0;
        virtual VectorXd solve(const VectorXd&)=0;
    };
    template<class LU> struct Adapter final:Solver {
        LU lu;
        void analyzePattern(const SparseMatrixd& a) override {lu.analyzePattern(a);}
        void factorize(const SparseMatrixd& a) override {lu.factorize(a);}
        Eigen::ComputationInfo info() const override {return lu.info();}
        VectorXd solve(const VectorXd& rhs) override {return lu.solve(rhs);}
    };
    std::unique_ptr<Solver> makeSolver() const {
        if(backend_=="sparselu_amd")return std::make_unique<Adapter<Eigen::SparseLU<SparseMatrixd,Eigen::AMDOrdering<SparseMatrixd::StorageIndex>>>>();
#if defined(VELA_HAS_UMFPACK)
        if(backend_=="umfpack")return std::make_unique<Adapter<Eigen::UmfPackLU<SparseMatrixd>>>();
#endif
        return std::make_unique<Adapter<Eigen::SparseLU<SparseMatrixd>>>();
    }
    std::string backend_;
    SparseMatrixd matrix_;
    std::unique_ptr<Solver> lu_;
    Eigen::Index rows_=0,cols_=0;
    std::vector<SparseMatrixd::StorageIndex> outer_,inner_;
    unsigned analyses_=0,factorizations_=0;
};
using ElectrothermalSparseLU=ElectrothermalDirectSolver;

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
