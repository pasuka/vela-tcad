#include "DirectBackend.h"
#include "vela/core/PerformanceProfiler.h"
#include <stdexcept>
#if defined(VELA_HAS_SUPERLU_MT)
#include <superlu_mt/slu_mt_ddefs.h>
#endif

namespace vela::detail {
#if defined(VELA_HAS_SUPERLU_MT)
class SuperLuMtBackend final : public DirectBackend {
    SparseMatrixd matrix_;
    SuperMatrix a_{},ac_{},l_{},u_{};
    superlumt_options_t options_{};
    Gstat_t stats_{};
    std::vector<int_t> pc_,pr_,tree_,counts_,partition_;
    int threads_=directBackendThreads();
    bool metis_,haveStats_=false,factored_=false;
    void prepare(yes_no_t repeat) {
        if(ac_.Store) {Destroy_CompCol_Permuted(&ac_);ac_.Store=nullptr;}
        pdgstrf_init(threads_,DOFACT,NOTRANS,repeat,sp_ienv(1),sp_ienv(2),1.,NO,0.,
            pc_.data(),pr_.data(),nullptr,0,&a_,&ac_,&options_,&stats_);
    }
public:
    explicit SuperLuMtBackend(bool metis):metis_(metis) {omp_set_dynamic(0);omp_set_num_threads(threads_);}
    ~SuperLuMtBackend() override {
        // options' structural arrays are owned by vectors, not pxgstrf_finalize.
        if(ac_.Store) Destroy_CompCol_Permuted(&ac_);
        if(a_.Store) Destroy_SuperMatrix_Store(&a_);
        if(l_.Store) Destroy_SuperNode_SCP(&l_);
        if(u_.Store) Destroy_CompCol_NCP(&u_);
        if(haveStats_) StatFree(&stats_);
    }
    void analyze(const SparseMatrixd& a) override {
        static_assert(sizeof(int_t)==sizeof(SparseMatrixd::StorageIndex));
        matrix_=a;const int n=static_cast<int>(a.rows());
        pc_.resize(n);pr_.resize(n);tree_.resize(n);counts_.resize(n);partition_.resize(n);
        options_.etree=tree_.data();options_.colcnt_h=counts_.data();options_.part_super_h=partition_.data();
        options_.SymmetricMode=NO;options_.PrintStat=NO;
        dCreate_CompCol_Matrix(&a_,n,n,static_cast<int>(a.nonZeros()),matrix_.valuePtr(),matrix_.innerIndexPtr(),matrix_.outerIndexPtr(),SLU_NC,SLU_D,SLU_GE);
        StatAlloc(n,threads_,sp_ienv(1),sp_ienv(2),&stats_);haveStats_=true;StatInit(n,threads_,&stats_);
        if(metis_) {const auto q=metisPermutation(a);for(int j=0;j<n;++j)pc_[q[j]]=j;}
        else get_perm_c(3,&a_,pc_.data()); // COLAMD; no numerical normal equations.
        prepare(NO);
    }
    void factor(const SparseMatrixd& a) override {
        std::copy(a.valuePtr(),a.valuePtr()+a.nonZeros(),matrix_.valuePtr());
        if(factored_) {StatInit(static_cast<int>(a.rows()),threads_,&stats_);prepare(YES);}
        int_t info=0;pdgstrf(&options_,&ac_,pr_.data(),&l_,&u_,&stats_,&info);
        if(info!=0) throw std::runtime_error("SuperLU_MT factorization info="+std::to_string(info));
        factored_=true;
    }
    VectorXd solve(const VectorXd& b) override {
        VectorXd x=b;SuperMatrix rhs{};int_t info=0;
        dCreate_Dense_Matrix(&rhs,static_cast<int>(b.size()),1,x.data(),static_cast<int>(b.size()),SLU_DN,SLU_D,SLU_GE);
        dgstrs(NOTRANS,&l_,&u_,pr_.data(),pc_.data(),&rhs,&stats_,&info);
        Destroy_SuperMatrix_Store(&rhs);
        if(info!=0) throw std::runtime_error("SuperLU_MT solve info="+std::to_string(info));return x;
    }
    void statistics() const override {
        observePerformanceValue("linear.numeric_factor_nonzeros_l",static_cast<SCPformat*>(l_.Store)->nnz);
        observePerformanceValue("linear.numeric_factor_nonzeros_u",static_cast<NCPformat*>(u_.Store)->nnz);
    }
};
#endif
std::unique_ptr<DirectBackend> makeSuperLuMt(bool metis) {
#if defined(VELA_HAS_SUPERLU_MT)
    return std::make_unique<SuperLuMtBackend>(metis);
#else
    (void)metis;throw std::invalid_argument("SuperLU_MT was not enabled in this build");
#endif
}
}
