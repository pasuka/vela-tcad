#pragma once
// Explicit diagnostics for the independent BVMethods experiment executable.
// No JSON simulation schema, process environment, or standard-entry default uses this scope.
#include "vela/core/Types.h"
#include <bit>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
namespace vela::experimental {
struct DCSweepExperiment {
    bool bootstrapModelGuardEnabled=false;
    bool reportModelQuality=false;
    std::string directory;
    std::string initialPackedStateFile;
    Real initialInnerVoltage_V=0.;
};
inline const DCSweepExperiment*& activeDCSweepExperiment() {
    static thread_local const DCSweepExperiment* active=nullptr;
    return active;
}
class DCSweepExperimentScope {
    const DCSweepExperiment* previous_;
public:
    explicit DCSweepExperimentScope(const DCSweepExperiment& options)
        :previous_(activeDCSweepExperiment()){activeDCSweepExperiment()=&options;}
    ~DCSweepExperimentScope(){activeDCSweepExperiment()=previous_;}
    DCSweepExperimentScope(const DCSweepExperimentScope&)=delete;
    DCSweepExperimentScope& operator=(const DCSweepExperimentScope&)=delete;
};
inline VectorXd readPackedExperiment(const std::string& file) {
    if constexpr(std::endian::native!=std::endian::little)
        throw std::runtime_error("experiment packed vectors require little-endian doubles");
    std::ifstream input(file,std::ios::binary|std::ios::ate);
    const auto bytes=input.tellg();
    if(!input||bytes<=0||bytes%sizeof(Real))throw std::runtime_error("invalid packed experiment input");
    VectorXd x(bytes/sizeof(Real));input.seekg(0);
    input.read(reinterpret_cast<char*>(x.data()),bytes);
    if(!input||!x.allFinite())throw std::runtime_error("invalid packed experiment values");
    return x;
}
inline void writePackedExperiment(const std::filesystem::path& file,const VectorXd& x) {
    std::ofstream output(file,std::ios::binary);
    output.write(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(Real));
    if(!output)throw std::runtime_error("cannot write packed experiment evidence");
}
} // namespace vela::experimental
