// Isolated experiment: never link vela_core compiled with different Eigen macros.
#include <omp.h> // GCC's C++ OpenMP overloads must precede Eigen's extern-C UMFPACK include.
#include "LinearReplaySystem.h"
#include <Eigen/UmfPackSupport>
#include <StrumpackSparseSolver.hpp>
#include <cstdlib>
#include <map>

extern "C" {
char* openblas_get_config();
int openblas_get_parallel();
int openblas_get_num_threads();
void openblas_set_num_threads(int);
}

static int threadOption(const char* name) {
    const char* value=std::getenv(name);
    const std::string s=value?value:"1";
    require(s=="1" || s=="2" || s=="4","Thread count must be 1, 2 or 4");
    return std::stoi(s);
}
static bool samePattern(const SparseMatrixd& a,const SparseMatrixd& b) {
    return a.rows()==b.rows() && a.nonZeros()==b.nonZeros() &&
        std::memcmp(a.outerIndexPtr(),b.outerIndexPtr(),sizeof(int)*(a.cols()+1))==0 &&
        std::memcmp(a.innerIndexPtr(),b.innerIndexPtr(),sizeof(int)*a.nonZeros())==0;
}
static void check(strumpack::ReturnCode rc) {
    require(rc==strumpack::ReturnCode::SUCCESS,"STRUMPACK operation failed");
}
struct StudyUmfPack : Eigen::UmfPackLU<SparseMatrixd> {
    double statistic(int index) const {return m_umfpackInfo[index];}
};

// Observe vendor phase boundaries without replacing any numerical operation.
class TimedStrumpack : public strumpack::SparseSolver<double,int> {
    using Base = strumpack::SparseSolver<double,int>;
    Clock::time_point phaseStart_;
protected:
    void perf_counters_start() override {
        Base::perf_counters_start();phaseStart_=Clock::now();
    }
    void perf_counters_stop(const std::string& name) override {
        phases[name]+=seconds(phaseStart_);Base::perf_counters_stop(name);
    }
public:
    TimedStrumpack():Base(false) {}
    std::map<std::string,double> phases;
};

int main(int argc,char** argv) {
    Json result={{"status","running"},{"systems",Json::array()}};
    if(argc<4) {std::cerr<<"Usage: study mode output.json capture.bin...\n";return 2;}
    if(std::filesystem::exists(argv[2])) {std::cerr<<"Refusing overwrite\n";return 2;}
    int exitCode=0;
    try {
        const std::string mode=argv[1];
        const bool sparse=mode=="sparselu",umf=mode=="umfpack";
        const bool blr=mode=="blr6" || mode=="blr8";
        const bool hss=mode=="hss6" || mode=="hss8";
        const bool orderingMode=mode=="strumpack_amd" || mode=="strumpack_mmd" || mode=="strumpack_and";
        require(sparse||umf||blr||hss||mode=="strumpack"||mode=="amalg"||orderingMode,"Unknown mode");
        const char* coldEnv=std::getenv("VELA_STUDY_COLD_ANALYSIS");
        require(!coldEnv || std::string(coldEnv)=="0" || std::string(coldEnv)=="1","Invalid cold-analysis switch");
        const bool cold=coldEnv && std::string(coldEnv)=="1";
        const char* resetEnv=std::getenv("VELA_STUDY_RESET_PER_DIRECTORY");
        require(!resetEnv || std::string(resetEnv)=="0" || std::string(resetEnv)=="1","Invalid directory-reset switch");
        const bool resetDirectory=resetEnv && std::string(resetEnv)=="1";
        require(!(cold && resetDirectory),"Conflicting analysis policies");
        const auto ordering=mode=="strumpack_amd"?strumpack::ReorderingStrategy::AMD:
            mode=="strumpack_mmd"?strumpack::ReorderingStrategy::MMD:
            mode=="strumpack_and"?strumpack::ReorderingStrategy::AND:strumpack::ReorderingStrategy::METIS;
        const int solverThreads=threadOption("VELA_LINEAR_THREADS");
        const int blasThreads=threadOption("VELA_BLAS_THREADS");
        // This OpenBLAS build schedules its BLAS worker queue on the calling
        // OpenMP team. For serial sparse algorithms, size that team for BLAS.
        const int ompThreads=sparse||umf?blasThreads:solverThreads;
        omp_set_dynamic(0);omp_set_max_active_levels(1);omp_set_num_threads(ompThreads);
        openblas_set_num_threads(blasThreads);
        require(openblas_get_num_threads()==blasThreads && omp_get_max_threads()==ompThreads,
                "Runtime thread settings do not match requested settings");
        result.update({{"mode",mode},{"openblas_config",openblas_get_config()},
            {"openblas_parallel",openblas_get_parallel()},{"blas_threads",openblas_get_num_threads()},
            {"omp_threads",omp_get_max_threads()},{"omp_max_active_levels",omp_get_max_active_levels()}});
        result["cold_analysis"]=cold;
        result["reset_per_directory"]=resetDirectory;
        if(!sparse && !umf) result["ordering"]=strumpack::get_name(ordering);
        if(blr||hss) result["compression"]={{"relative_tolerance",mode.back()=='6'?1e-6:1e-8},
            {"absolute_tolerance",0.},{"min_separator",128},{"min_front",256},{"leaf_size",64},
            {"outer_solver","PREC_GMRES"},{"outer_relative_tolerance",1e-13},{"max_iterations",100}};
#ifdef EIGEN_USE_BLAS
        result["eigen_use_blas"]=true;
#else
        result["eigen_use_blas"]=false;
#endif
        Eigen::SparseLU<SparseMatrixd,Eigen::COLAMDOrdering<int>> eigen;
        StudyUmfPack umfpack;
        std::unique_ptr<TimedStrumpack> strum;
        Eigen::SparseMatrix<double,Eigen::RowMajor,int> csr;
        SparseMatrixd previous;
        std::filesystem::path previousDirectory;
        double total=0.;
        for(int index=3;index<argc;++index) {
            result["active_input"]=argv[index];
            System s(argv[index]);
            const auto directory=std::filesystem::path(argv[index]).parent_path();
            const bool analyze=cold || index==3 || !samePattern(previous,s.matrix) ||
                (resetDirectory && directory!=previousDirectory);
            const auto start=Clock::now();
            double preparation=0.,analysis=0.,factor=0.,reorderSeconds=0.;
            VectorXd x,y,rhs2=2.*s.rhs;
            std::int64_t entries=0;
            int iterations=0;
            if(sparse || umf) {
                auto run=[&](auto& solver) {
                    auto t=Clock::now();
                    if(analyze) solver.analyzePattern(s.matrix);
                    analysis=seconds(t);t=Clock::now();
                    solver.factorize(s.matrix);factor=seconds(t);
                    require(solver.info()==Eigen::Success,"Factorization failed");
                    x=solver.solve(s.rhs);require(solver.info()==Eigen::Success,"First RHS failed");
                    y=solver.solve(rhs2);require(solver.info()==Eigen::Success,"Second RHS failed");
                };
                if(sparse) {run(eigen);entries=eigen.nnzL()+eigen.nnzU();}
                else {run(umfpack);entries=static_cast<std::int64_t>(umfpack.statistic(UMFPACK_LNZ)+umfpack.statistic(UMFPACK_UNZ));}
            } else {
                auto t=Clock::now();csr=s.matrix;preparation=seconds(t);t=Clock::now();
                if(analyze) {
                    strum=std::make_unique<TimedStrumpack>();
                    auto& opts=strum->options();
                    opts.set_reordering_method(ordering);
                    opts.set_matching(strumpack::MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING);
                    opts.disable_replace_tiny_pivots();opts.disable_gpu();
                    opts.set_compression(blr?strumpack::CompressionType::BLR:hss?strumpack::CompressionType::HSS:strumpack::CompressionType::NONE);
                    opts.set_Krylov_solver(blr||hss?strumpack::KrylovSolver::PREC_GMRES:strumpack::KrylovSolver::DIRECT);
                    opts.set_rel_tol(1e-13);opts.set_abs_tol(0.);opts.set_maxit(100);opts.set_gmres_restart(30);
                    if(blr||hss) {
                        opts.set_compression_rel_tol(mode.back()=='6'?1e-6:1e-8);
                        opts.set_compression_abs_tol(0.);
                        opts.set_compression_min_sep_size(128);opts.set_compression_min_front_size(256);
                        opts.set_compression_leaf_size(64);
                    }
                    if(mode=="amalg") {opts.enable_MUMPS_SYMQAMD();opts.enable_agg_amalg();}
                    strum->set_csr_matrix(csr.rows(),csr.outerIndexPtr(),csr.innerIndexPtr(),csr.valuePtr(),false);
                    const auto reorderStart=Clock::now();check(strum->reorder());reorderSeconds=seconds(reorderStart);
                } else {
                    strum->phases.clear();
                    strum->update_matrix_values(csr.rows(),csr.outerIndexPtr(),csr.innerIndexPtr(),csr.valuePtr(),false);
                }
                analysis=seconds(t);t=Clock::now();check(strum->factor());factor=seconds(t);
                x.resize(s.rhs.size());y.resize(s.rhs.size());
                check(strum->solve(s.rhs.data(),x.data()));iterations=strum->Krylov_iterations();
                check(strum->solve(rhs2.data(),y.data()));iterations+=strum->Krylov_iterations();
                entries=strum->factor_nonzeros();
            }
            const double elapsed=seconds(start);total+=elapsed;
            require(openblas_get_num_threads()==blasThreads && omp_get_max_threads()==ompThreads &&
                    omp_get_max_active_levels()==1,"Runtime threads changed during matrix solve");
            auto q=quality(s,x);
            if(!s.rawAvailable) {
                q.erase("raw");q["relative_solver_step_difference_by_block"]=q["relative_raw_step_difference_by_block"];
                q.erase("relative_raw_step_difference_by_block");
            }
            q["raw_available"]=s.rawAvailable;
            const double repeat=(y-2.*x).norm()/std::max(x.norm(),1e-300);
            const bool pass=q["scaled"]["normwise_backward_error"].get<double>()<1e-12 &&
                (!s.rawAvailable || q["raw"]["normwise_backward_error"].get<double>()<1e-12) &&
                repeat<1e-10 && y.allFinite();
            result["systems"].push_back({{"input",argv[index]},{"analysis_reused",!analyze},
                {"preparation_seconds",preparation},{"analysis_seconds",analysis},{"factor_seconds",factor},
                {"total_seconds",elapsed},{"factor_entries",entries},{"two_rhs_krylov_iterations",iterations},{"quality",q},
                {"two_rhs_relative_difference",repeat},{"pass",pass}});
            auto& row=result["systems"].back();
            row["matrix_rows"]=s.matrix.rows();row["matrix_nonzeros"]=s.matrix.nonZeros();
            row["factor_fill_ratio"]=double(entries)/s.matrix.nonZeros();
            row["runtime_threads_verified"]=true;
            if(strum) {
                if(analyze) require(strum->phases.contains("nested dissection") &&
                    strum->phases.contains("symbolic factorization"),"Missing vendor phase timing");
                row["vendor_phase_seconds"]=strum->phases;
                row["reorder_seconds"]=reorderSeconds;
                // This remainder includes matching/equilibration/symmetrization
                // and uninstrumented overhead; it is not a matching-only timer.
                row["reorder_other_seconds"]=reorderSeconds-strum->phases["nested dissection"]-
                    strum->phases["symbolic factorization"];
            }
            result["linear_seconds"]=total;
            require(pass,"Original matrix/RHS accuracy gate failed");
            previous=s.matrix;
            previousDirectory=directory;
        }
        result.erase("active_input");result["status"]="pass";
    } catch(const std::exception& e) {result["status"]="failed";result["error"]=e.what();exitCode=1;}
    std::ofstream output(argv[2]);output<<result.dump(2)<<'\n';
    return output.good()?exitCode:2;
}
