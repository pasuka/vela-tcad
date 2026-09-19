#include "DirectBackend.h"
#include "vela/core/PerformanceProfiler.h"
#include <stdexcept>
#if defined(VELA_HAS_STRUMPACK)
#include <StrumpackSparseSolver.hpp>
#include <omp.h>
#endif

namespace vela::detail {
#if defined(VELA_HAS_STRUMPACK)
class StrumpackBackend final : public DirectBackend {
    strumpack::SparseSolver<double,int> solver_{false};
    std::vector<int> ptr_,cols_,source_;
    std::vector<double> values_;
    int n_=0;
    bool factored_=false;
    static void check(strumpack::ReturnCode rc,const char* phase) {
        if(rc!=strumpack::ReturnCode::SUCCESS)
            throw std::runtime_error(std::string("STRUMPACK ")+phase+" code="+std::to_string(static_cast<int>(rc)));
    }
public:
    StrumpackBackend() {
        omp_set_dynamic(0);omp_set_num_threads(directBackendThreads());
        auto& opts=solver_.options();
        opts.set_compression(strumpack::CompressionType::NONE);
        opts.set_Krylov_solver(strumpack::KrylovSolver::DIRECT);
        opts.set_reordering_method(strumpack::ReorderingStrategy::METIS);
        opts.set_matching(strumpack::MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING);
        opts.disable_replace_tiny_pivots();opts.disable_gpu();
    }
    void analyze(const SparseMatrixd& a) override {
        n_=static_cast<int>(a.rows());ptr_.assign(n_+1,0);
        for(int k=0;k<a.nonZeros();++k)++ptr_[a.innerIndexPtr()[k]+1];
        for(int i=0;i<n_;++i)ptr_[i+1]+=ptr_[i];
        auto next=ptr_;cols_.resize(a.nonZeros());source_.resize(a.nonZeros());values_.resize(a.nonZeros());
        for(int j=0;j<a.cols();++j)for(int k=a.outerIndexPtr()[j];k<a.outerIndexPtr()[j+1];++k) {
            const int target=next[a.innerIndexPtr()[k]]++;cols_[target]=j;source_[target]=k;values_[target]=a.valuePtr()[k];
        }
        solver_.set_csr_matrix(n_,ptr_.data(),cols_.data(),values_.data(),false);
        check(solver_.reorder(),"reorder");
    }
    void factor(const SparseMatrixd& a) override {
        if(factored_) {
            for(std::size_t k=0;k<source_.size();++k)values_[k]=a.valuePtr()[source_[k]];
            solver_.update_matrix_values(n_,ptr_.data(),cols_.data(),values_.data(),false);
        }
        check(solver_.factor(),"factor");factored_=true;
    }
    VectorXd solve(const VectorXd& b) override {
        VectorXd x(b.size());check(solver_.solve(b.data(),x.data()),"solve");return x;
    }
    void statistics() const override {
        observePerformanceValue("linear.strumpack_factor_entries",solver_.factor_nonzeros());
        observePerformanceValue("linear.strumpack_factor_bytes",solver_.factor_memory());
    }
};
#endif
std::unique_ptr<DirectBackend> makeStrumpack() {
#if defined(VELA_HAS_STRUMPACK)
    return std::make_unique<StrumpackBackend>();
#else
    throw std::invalid_argument("STRUMPACK was not enabled in this build");
#endif
}
}
