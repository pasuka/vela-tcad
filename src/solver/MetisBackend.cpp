#include "DirectBackend.h"
#include "vela/core/PerformanceProfiler.h"
#include <Eigen/SparseLU>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#if defined(VELA_HAS_METIS)
#include <metis.h>
#endif

namespace vela::detail {
int directBackendThreads() {
    const char* s=std::getenv("VELA_LINEAR_THREADS");
    if(!s || !*s) return 1;
    if(std::string(s)=="1") return 1;
    if(std::string(s)=="2") return 2;
    if(std::string(s)=="4") return 4;
    throw std::invalid_argument("VELA_LINEAR_THREADS must be 1, 2 or 4");
}
std::vector<int> metisPermutation(const SparseMatrixd& a) {
#if defined(VELA_HAS_METIS)
    if (a.rows()!=a.cols() || a.rows()<1 || a.rows()>std::numeric_limits<int>::max())
        throw std::invalid_argument("METIS requires a nonempty square 32-bit system");
    idx_t n=static_cast<idx_t>(a.rows());
    std::vector<std::vector<idx_t>> graph(n);
    for(int col=0;col<a.outerSize();++col)
        for(SparseMatrixd::InnerIterator it(a,col);it;++it) if(it.row()!=it.col()) {
            graph[it.row()].push_back(static_cast<idx_t>(it.col()));
            graph[it.col()].push_back(static_cast<idx_t>(it.row()));
        }
    std::vector<idx_t> ptr(n+1),edges,perm(n),inverse(n);
    for(idx_t i=0;i<n;++i) {
        auto& row=graph[i];std::sort(row.begin(),row.end());
        row.erase(std::unique(row.begin(),row.end()),row.end());
        if(edges.size()+row.size()>static_cast<std::size_t>(std::numeric_limits<idx_t>::max()))
            throw std::overflow_error("METIS graph index overflow");
        edges.insert(edges.end(),row.begin(),row.end());ptr[i+1]=static_cast<idx_t>(edges.size());
    }
    idx_t opts[METIS_NOPTIONS];METIS_SetDefaultOptions(opts);
    opts[METIS_OPTION_SEED]=0;opts[METIS_OPTION_NUMBERING]=0;
    idx_t empty=0;
    const int rc=METIS_NodeND(&n,ptr.data(),edges.empty()?&empty:edges.data(),nullptr,opts,perm.data(),inverse.data());
    if(rc!=METIS_OK) throw std::runtime_error("METIS_NodeND failed: "+std::to_string(rc));
    std::vector<int> result(n);std::vector<bool> seen(n);
    for(idx_t i=0;i<n;++i) {
        if(perm[i]<0 || perm[i]>=n || seen[perm[i]] || inverse[perm[i]]!=i)
            throw std::runtime_error("METIS returned invalid permutation");
        seen[perm[i]]=true;result[i]=static_cast<int>(perm[i]);
    }
    return result;
#else
    (void)a;throw std::invalid_argument("METIS was not enabled in this build");
#endif
}

#if defined(VELA_HAS_METIS)
class MetisSparseLU final : public DirectBackend {
    SparseMatrixd matrix_;
    std::vector<int> perm_,source_;
    Eigen::SparseLU<SparseMatrixd,Eigen::NaturalOrdering<int>> solver_;
public:
    void analyze(const SparseMatrixd& a) override {
        perm_=metisPermutation(a);
        matrix_.resize(a.rows(),a.cols());matrix_.reserve(a.nonZeros());source_.clear();
        for(int j=0;j<a.cols();++j) {
            matrix_.startVec(j);
            for(int k=a.outerIndexPtr()[perm_[j]];k<a.outerIndexPtr()[perm_[j]+1];++k) {
                matrix_.insertBack(a.innerIndexPtr()[k],j)=a.valuePtr()[k];source_.push_back(k);
            }
        }
        matrix_.finalize();solver_.analyzePattern(matrix_);
    }
    void factor(const SparseMatrixd& a) override {
        for(std::size_t k=0;k<source_.size();++k) matrix_.valuePtr()[k]=a.valuePtr()[source_[k]];
        solver_.factorize(matrix_);
        if(solver_.info()!=Eigen::Success) throw std::runtime_error("SparseLU METIS factorization failed");
    }
    VectorXd solve(const VectorXd& b) override {
        VectorXd y=solver_.solve(b),x(b.size());
        if(solver_.info()!=Eigen::Success) throw std::runtime_error("SparseLU METIS solve failed");
        for(int j=0;j<b.size();++j) x[perm_[j]]=y[j];
        return x;
    }
    void statistics() const override {
        observePerformanceValue("linear.numeric_factor_nonzeros_l",solver_.nnzL());
        observePerformanceValue("linear.numeric_factor_nonzeros_u",solver_.nnzU());
    }
};
#endif
std::unique_ptr<DirectBackend> makeMetisSparseLU() {
#if defined(VELA_HAS_METIS)
    return std::make_unique<MetisSparseLU>();
#else
    throw std::invalid_argument("METIS was not enabled in this build");
#endif
}
}
