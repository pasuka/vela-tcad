#include "LinearReplaySystem.h"
int main(int argc, char** argv) {
    std::string currentFile;
    try {
        require(argc >= 4, "Usage: linear_solver_replay backend output.json capture.bin...");
        require(!std::filesystem::exists(argv[2]), "Refusing to overwrite results");
        LinearSolver solver(argv[1]);
        Json output={{"backend",argv[1]},{"systems",Json::array()}};
        for(int i=3;i<argc;++i) {
            currentFile=argv[i];
            System system(argv[i]);
            PerformanceProfiler profiler({true,"unused.json"});
            VectorXd x, repeated;
            {
                ActivePerformanceProfilerScope active(&profiler);
                x=solver.solve(system.matrix,system.rhs);
                // A changed RHS must use retained factors, including refinement.
                repeated=solver.solve(system.matrix,2.*system.rhs);
            }
            auto q=quality(system,x);
            q["raw_available"]=system.rawAvailable;
            if(!(q["scaled"]["normwise_backward_error"].get<double>()<1e-12))
                throw std::runtime_error("Scaled backward error exceeds replay gate: "+q["scaled"].dump());
            require(!system.rawAvailable || q["raw"]["normwise_backward_error"].get<double>()<1e-12,
                    "Original backward error exceeds replay gate");
            if(!system.rawAvailable) {
                q.erase("raw");q["coordinate_scope"]="solver_input_only";
                q["relative_solver_step_difference_by_block"]=q.at("relative_raw_step_difference_by_block");
                q.erase("relative_raw_step_difference_by_block");
            }
            require((repeated-2.*x).norm()/std::max(x.norm(),1e-300)<1e-10,
                    "Retained factor RHS solution differs");
            output["systems"].push_back({{"file",argv[i]},{"rows",system.matrix.rows()},
                {"nnz",system.matrix.nonZeros()},{"quality",q},
                {"cumulative_analyses",solver.patternAnalysisCount()},
                {"profiling",profiler.toJson()}});
        }
        output["status"]="pass";
        std::ofstream stream(argv[2]);stream.exceptions(std::ios::badbit|std::ios::failbit);
        stream<<output.dump(2)<<'\n';
        std::cout<<"REPLAY_PASS "<<output["systems"].size()<<" systems\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<currentFile<<": "<<e.what()<<'\n';return 1;}
}
