#include "DirectBackend.h"
#include "vela/core/PerformanceProfiler.h"
#include <array>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#define VELA_MKL_CALL __cdecl
#elif defined(__linux__)
#include <dlfcn.h>
#define VELA_MKL_CALL
#endif

namespace vela::detail {
#if defined(_WIN32) || defined(__linux__)
static_assert(sizeof(int) == 4 && std::is_same_v<Real, double>,
              "MKL LP64 adapter requires 32-bit int and binary64 Real");
// Explicit runtime-only integration: no MKL dependency in ordinary builds.
// LP64 C ABI, real nonsymmetric matrices; all input/output scalars are double.
class MklPardiso final : public DirectBackend {
    using Pardiso = void(VELA_MKL_CALL*)(void**,const int*,const int*,const int*,const int*,
        const int*,const void*,const int*,const int*,int*,const int*,int*,const int*,void*,void*,int*);
    using Init = void(VELA_MKL_CALL*)(void**,const int*,int*);
    using SetInt = int(VELA_MKL_CALL*)(int);
#ifdef _WIN32
    HMODULE module_=nullptr;
#else
    void* module_=nullptr;
#endif
    Pardiso pardiso_=nullptr;
    SetInt threads_=nullptr;
    std::array<void*,64> pt_{};
    std::array<int,64> ip_{};
    std::vector<int> rows_,cols_,source_,perm_;
    std::vector<double> values_;
    int n_=0,threadCount_=1;
    bool started_=false;
    template<class T> T symbol(const char* name) {
#ifdef _WIN32
        auto p=GetProcAddress(module_,name);
#else
        auto p=dlsym(module_,name);
#endif
        if(!p) throw std::runtime_error(std::string("MKL missing symbol ")+name);
        return reinterpret_cast<T>(p);
    }
    void unload() noexcept {
        if(!module_) return;
#ifdef _WIN32
        FreeLibrary(module_);
#else
        dlclose(module_);
#endif
        module_=nullptr;
    }
    int call(int phase,double* b,double* x) noexcept {
        int one=1,type=11,msg=0,error=0;
        int previous=threads_(threadCount_);
        pardiso_(pt_.data(),&one,&one,&type,&phase,&n_,values_.data(),rows_.data(),
                 cols_.data(),perm_.data(),&one,ip_.data(),&msg,b,x,&error);
        threads_(previous);
        return error;
    }
    void checked(int phase,double* b=nullptr,double* x=nullptr) {
        int error=call(phase,b,x);
        if(error) throw std::runtime_error("MKL PARDISO phase "+std::to_string(phase)+" error "+std::to_string(error));
    }
public:
    MklPardiso() {
        const char* path=std::getenv("VELA_MKL_RT_LIBRARY");
        if(!path || !*path) path=std::getenv("VELA_MKL_RT_DLL"); // Existing Windows configurations.
        if(!path || !*path || !std::filesystem::path(path).is_absolute())
            throw std::invalid_argument("pardiso_mkl requires an absolute VELA_MKL_RT_LIBRARY (or VELA_MKL_RT_DLL) path");
        threadCount_=directBackendThreads();
#ifdef _WIN32
        module_=LoadLibraryExW(std::filesystem::path(path).c_str(),nullptr,LOAD_WITH_ALTERED_SEARCH_PATH);
        if(!module_) throw std::runtime_error("MKL runtime load failed, Windows error "+std::to_string(GetLastError()));
#else
        module_=dlopen(path,RTLD_NOW|RTLD_LOCAL);
        if(!module_) throw std::runtime_error(std::string("MKL runtime load failed: ")+dlerror());
#endif
        try {
            auto layer=symbol<SetInt>("MKL_Set_Interface_Layer");
            if(layer(0)!=0) throw std::runtime_error("MKL LP64 interface selection failed");
            threads_=symbol<SetInt>("MKL_Set_Num_Threads_Local");
            pardiso_=symbol<Pardiso>("pardiso");
            auto init=symbol<Init>("pardisoinit");int type=11;init(pt_.data(),&type,ip_.data());
            ip_[0]=1;ip_[1]=2;ip_[7]=2;ip_[9]=13;ip_[10]=1;ip_[12]=1;ip_[26]=1;ip_[27]=0;ip_[34]=0;
        } catch(...) {unload();throw;}
    }
    ~MklPardiso() override {
        if(started_) call(-1,nullptr,nullptr);
        unload();
    }
    void analyze(const SparseMatrixd& a) override {
        if(a.rows()>=std::numeric_limits<int>::max() || a.nonZeros()>=std::numeric_limits<int>::max())
            throw std::invalid_argument("MKL LP64 matrix too large");
        n_=static_cast<int>(a.rows());rows_.assign(n_+1,0);
        for(int k=0;k<a.nonZeros();++k)++rows_[a.innerIndexPtr()[k]+1];
        for(int i=0;i<n_;++i)rows_[i+1]+=rows_[i];
        auto next=rows_;cols_.resize(a.nonZeros());source_.resize(a.nonZeros());values_.resize(a.nonZeros());
        for(int j=0;j<a.cols();++j)for(int k=a.outerIndexPtr()[j];k<a.outerIndexPtr()[j+1];++k){
            int dst=next[a.innerIndexPtr()[k]]++;cols_[dst]=j+1;source_[dst]=k;values_[dst]=a.valuePtr()[k];
        }
        for(auto& r:rows_)++r;
        perm_.assign(n_,0);started_=true;checked(11);
    }
    void factor(const SparseMatrixd& a) override {
        for(std::size_t k=0;k<source_.size();++k)values_[k]=a.valuePtr()[source_[k]];
        checked(22);
    }
    VectorXd solve(const VectorXd& b) override {
        VectorXd rhs=b,x=VectorXd::Zero(n_);checked(33,rhs.data(),x.data());
        observePerformanceValue("linear.mkl_refinement_steps",ip_[6]);
        observePerformanceValue("linear.mkl_perturbed_pivots",ip_[13]);
        observePerformanceValue("linear.mkl_threads",threadCount_);
        return x;
    }
};
#endif
std::unique_ptr<DirectBackend> makeMklPardiso() {
#if defined(_WIN32) || defined(__linux__)
    return std::make_unique<MklPardiso>();
#else
    throw std::invalid_argument("pardiso_mkl runtime backend requires Windows or Linux");
#endif
}
} // namespace vela::detail
