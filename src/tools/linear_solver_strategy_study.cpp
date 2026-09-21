// Independent frozen-matrix experiment. No controls are added to production.
#include "LinearReplaySystem.h"
#include <numeric>
#if defined(VELA_HAS_UMFPACK)
#include <umfpack.h>
#endif

static std::string permutationHash(const std::vector<int>& p) {
    std::uint64_t h=14695981039346656037ULL;
    for(auto i:p) for(int k=0;k<4;++k) { h^=(static_cast<unsigned>(i)>>(8*k))&255;h*=1099511628211ULL; }
    return std::to_string(h);
}
static void gate(const System& s,const VectorXd& x,Json& record) {
    record["quality"]=quality(s,x);
    require(record["quality"]["raw"]["normwise_backward_error"].get<double>()<1e-12,"Raw error gate failed");
    require(record["quality"]["scaled"]["normwise_backward_error"].get<double>()<1e-12,"Scaled error gate failed");
}
template<class Ordering>
static Json eigenRun(const std::vector<System>& systems,std::vector<int>& q) {
    Eigen::SparseLU<SparseMatrixd,Ordering> lu;
    auto start=Clock::now();lu.analyzePattern(systems.front().matrix);
    const double analysis=seconds(start);
    q.resize(systems.front().matrix.cols());
    // Eigen stores old->new; UMFPACK Q[k] is the old index of new column k.
    for(int old=0;old<static_cast<int>(q.size());++old) q[lu.colsPermutation().indices()[old]]=old;
    Json result={{"analysis_seconds",analysis},{"actual_q_fnv1a64",permutationHash(q)},{"systems",Json::array()}};
    for(const auto& s:systems) {
        start=Clock::now();lu.factorize(s.matrix);const auto factor=seconds(start);
        require(lu.info()==Eigen::Success,"SparseLU factor failed");
        start=Clock::now();const VectorXd x=lu.solve(s.rhs);const VectorXd y=lu.solve(2.*s.rhs);
        const auto solve=seconds(start);
        require(lu.info()==Eigen::Success && (y-2.*x).norm()/std::max(x.norm(),1e-300)<1e-10,"SparseLU RHS gate failed");
        Json row={{"factor_seconds",factor},{"two_rhs_seconds",solve},{"lnz",lu.nnzL()},{"unz",lu.nnzU()},
            {"fill",double(lu.nnzL()+lu.nnzU())/s.matrix.nonZeros()}};
        gate(s,x,row);result["systems"].push_back(row);
    }
    return result;
}
#if defined(VELA_HAS_UMFPACK)
struct Handles {
    void* symbolic=nullptr;void* numeric=nullptr;
    ~Handles(){umfpack_di_free_numeric(&numeric);umfpack_di_free_symbolic(&symbolic);}
};
static Json umfRun(const std::vector<System>& systems,const std::string& mode,const std::vector<int>& q) {
    Handles h;double control[UMFPACK_CONTROL],info[UMFPACK_INFO];umfpack_di_defaults(control);
    const bool given=mode.find("given_")==0;
    if(given) {
        control[UMFPACK_FIXQ]=1;control[UMFPACK_SINGLETONS]=0;
        control[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_SYMMETRIC;
    }
    if(mode.find("unsym")!=std::string::npos)control[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_UNSYMMETRIC;
    if(mode.find("noscale")!=std::string::npos)control[UMFPACK_SCALE]=UMFPACK_SCALE_NONE;
    if(mode.find("pivot1")!=std::string::npos) {
        control[UMFPACK_PIVOT_TOLERANCE]=1.;control[UMFPACK_SYM_PIVOT_TOLERANCE]=1.;
    }
    const auto& a=systems.front().matrix;auto start=Clock::now();
    const int status=given?umfpack_di_qsymbolic(a.rows(),a.cols(),a.outerIndexPtr(),a.innerIndexPtr(),a.valuePtr(),q.data(),&h.symbolic,control,info)
        :umfpack_di_symbolic(a.rows(),a.cols(),a.outerIndexPtr(),a.innerIndexPtr(),a.valuePtr(),&h.symbolic,control,info);
    const auto analysis=seconds(start);require(status==UMFPACK_OK,"UMFPACK analysis failed");
    Json result={{"analysis_seconds",analysis},{"control",std::vector<double>(control,control+UMFPACK_CONTROL)},
        {"strategy_used",info[UMFPACK_STRATEGY_USED]},{"ordering_used",info[UMFPACK_ORDERING_USED]},
        {"input_q_fnv1a64",given?Json(permutationHash(q)):Json(nullptr)},{"systems",Json::array()}};
    for(const auto& s:systems) {
        const auto& m=s.matrix;umfpack_di_free_numeric(&h.numeric);start=Clock::now();
        const auto code=umfpack_di_numeric(m.outerIndexPtr(),m.innerIndexPtr(),m.valuePtr(),h.symbolic,&h.numeric,control,info);
        const auto factor=seconds(start);require(code==UMFPACK_OK,"UMFPACK factor failed");
        Json row={{"factor_seconds",factor},{"lnz",info[UMFPACK_LNZ]},{"unz",info[UMFPACK_UNZ]},
            {"fill",(info[UMFPACK_LNZ]+info[UMFPACK_UNZ])/m.nonZeros()},
            {"flops",info[UMFPACK_FLOPS]},{"internal_peak_bytes",info[UMFPACK_PEAK_MEMORY]*info[UMFPACK_SIZE_OF_UNIT]},
            {"off_diagonal_pivots",info[UMFPACK_NOFF_DIAG]}};
        std::vector<int> actual(m.cols());
        require(umfpack_di_get_numeric(nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,actual.data(),nullptr,nullptr,nullptr,h.numeric)==UMFPACK_OK,"Read Q failed");
        row["actual_q_fnv1a64"]=permutationHash(actual);row["matches_input_q"]=given?Json(q==actual):Json(nullptr);
        VectorXd x(m.cols()),y(m.cols()),rhs2=2.*s.rhs;start=Clock::now();
        const int sx=umfpack_di_solve(UMFPACK_A,m.outerIndexPtr(),m.innerIndexPtr(),m.valuePtr(),x.data(),s.rhs.data(),h.numeric,control,info);
        const int sy=umfpack_di_solve(UMFPACK_A,m.outerIndexPtr(),m.innerIndexPtr(),m.valuePtr(),y.data(),rhs2.data(),h.numeric,control,info);
        row["two_rhs_seconds"]=seconds(start);
        require(sx==UMFPACK_OK && sy==UMFPACK_OK && (y-2.*x).norm()/std::max(x.norm(),1e-300)<1e-10,"UMFPACK RHS gate failed");
        gate(s,x,row);result["systems"].push_back(row);
    }
    return result;
}
#endif
int main(int argc,char** argv) {
    try {
        require(argc>=4,"Usage: linear_solver_strategy_study MODE output.json capture.bin...");
        require(!std::filesystem::exists(argv[2]),"Refusing to overwrite results");
        std::vector<System> systems;for(int i=3;i<argc;++i)systems.emplace_back(argv[i]);
        const auto& a=systems.front().matrix;
        for(const auto& s:systems)require(s.matrix.rows()==a.rows() && s.matrix.nonZeros()==a.nonZeros() &&
            std::memcmp(a.outerIndexPtr(),s.matrix.outerIndexPtr(),4*(a.cols()+1))==0 &&
            std::memcmp(a.innerIndexPtr(),s.matrix.innerIndexPtr(),4*a.nonZeros())==0,"Pattern changed within study");
        const std::string mode=argv[1];Json result;std::vector<int> q;
        if(mode=="eigen_colamd")result=eigenRun<Eigen::COLAMDOrdering<int>>(systems,q);
        else if(mode=="eigen_amd")result=eigenRun<Eigen::AMDOrdering<int>>(systems,q);
        else {
#if defined(VELA_HAS_UMFPACK)
            const std::vector<std::string> allowed{"default","unsym","noscale","pivot1","given_colamd","given_amd","given_colamd_noscale","given_colamd_unsym","given_colamd_pivot1"};
            require(std::find(allowed.begin(),allowed.end(),mode)!=allowed.end(),"Unknown study mode");
            // Derive the same final (postordered) Eigen column order outside timed UMFPACK service.
            if(mode.find("given_")==0) {
                const auto prepare=[&]<class O>() {
                    Eigen::SparseLU<SparseMatrixd,O> lu;lu.analyzePattern(a);q.resize(a.cols());
                    for(int old=0;old<a.cols();++old)q[lu.colsPermutation().indices()[old]]=old;
                };
                if(mode=="given_amd")prepare.template operator()<Eigen::AMDOrdering<int>>();
                else prepare.template operator()<Eigen::COLAMDOrdering<int>>();
            }
            result=umfRun(systems,mode,q);
#else
            throw std::runtime_error("UMFPACK unavailable");
#endif
        }
        result["mode"]=mode;result["status"]="pass";
        std::ofstream stream(argv[2]);stream.exceptions(std::ios::badbit|std::ios::failbit);stream<<result.dump(2)<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
