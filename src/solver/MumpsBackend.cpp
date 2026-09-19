#include "DirectBackend.h"
#include "vela/core/PerformanceProfiler.h"
#include <stdexcept>
#if defined(VELA_HAS_MUMPS)
#include <dmumps_c.h>
#include <omp.h>
#endif

namespace vela::detail {
#if defined(VELA_HAS_MUMPS)
class MumpsBackend final : public DirectBackend {
    DMUMPS_STRUC_C id_{};
    std::vector<MUMPS_INT> rows_,cols_;
    std::vector<double> values_;
    bool initialized_=false;
    void run(int job) {
        id_.job=job;dmumps_c(&id_);
        if(id_.infog[0]<0)
            throw std::runtime_error("MUMPS job "+std::to_string(job)+" INFOG(1)="+
                std::to_string(id_.infog[0])+" INFOG(2)="+std::to_string(id_.infog[1]));
    }
public:
    explicit MumpsBackend(bool metis) {
        omp_set_dynamic(0);omp_set_num_threads(directBackendThreads());
        id_.sym=0;id_.par=1;id_.comm_fortran=-987654;run(-1);initialized_=true;
        id_.icntl[0]=-1;id_.icntl[1]=-1;id_.icntl[2]=-1;id_.icntl[3]=0;
        if(metis) id_.icntl[6]=5;
    }
    ~MumpsBackend() override {if(initialized_) {id_.job=-2;dmumps_c(&id_);}}
    void analyze(const SparseMatrixd& a) override {
        id_.n=static_cast<MUMPS_INT>(a.rows());id_.nz=static_cast<MUMPS_INT>(a.nonZeros());
        rows_.resize(a.nonZeros());cols_.resize(a.nonZeros());values_.assign(a.valuePtr(),a.valuePtr()+a.nonZeros());
        for(int j=0;j<a.outerSize();++j) for(int k=a.outerIndexPtr()[j];k<a.outerIndexPtr()[j+1];++k) {
            rows_[k]=a.innerIndexPtr()[k]+1;cols_[k]=j+1;
        }
        id_.irn=rows_.data();id_.jcn=cols_.data();id_.a=values_.data();run(1);
    }
    void factor(const SparseMatrixd& a) override {
        std::copy(a.valuePtr(),a.valuePtr()+a.nonZeros(),values_.begin());run(2);
    }
    VectorXd solve(const VectorXd& b) override {
        VectorXd x=b;id_.rhs=x.data();id_.nrhs=1;id_.lrhs=static_cast<MUMPS_INT>(b.size());
        try {run(3);} catch(...) {id_.rhs=nullptr;throw;}
        id_.rhs=nullptr;return x;
    }
    void statistics() const override {
        observePerformanceValue("linear.mumps_ordering_used",id_.infog[6]);
        if(id_.infog[19]>=0) observePerformanceValue("linear.mumps_factor_entries",id_.infog[19]);
        if(id_.rinfog[1]>=0) observePerformanceValue("linear.mumps_factor_flops",id_.rinfog[1]);
    }
};
#endif
std::unique_ptr<DirectBackend> makeMumps(bool metis) {
#if defined(VELA_HAS_MUMPS)
    return std::make_unique<MumpsBackend>(metis);
#else
    (void)metis;throw std::invalid_argument("MUMPS was not enabled in this build");
#endif
}
}
