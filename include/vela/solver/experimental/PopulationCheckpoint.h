#pragma once
#include "vela/io/DDSolutionState.h"
#include "vela/io/StateIdentity.h"
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela::experimental {
namespace population_checkpoint_detail {
inline std::string bits(double value) {
    const auto u = std::bit_cast<std::uint64_t>(value);
    std::string s(16, '0');
    constexpr char hex[] = "0123456789abcdef";
    for (int i=0; i<16; ++i) s[15-i] = hex[(u >> (4*i)) & 15];
    return s;
}
inline double value(const nlohmann::json& encoded) {
    const auto s = encoded.get<std::string>();
    if (s.size()!=16) throw std::invalid_argument("Invalid population bit encoding");
    std::uint64_t u=0;
    for (const char c:s) {
        if (!((c>='0'&&c<='9')||(c>='a'&&c<='f')))
            throw std::invalid_argument("Invalid population bit encoding");
        u=(u<<4) | static_cast<unsigned>(c<='9'?c-'0':c-'a'+10);
    }
    return std::bit_cast<double>(u);
}
inline std::string archiveIdentity(const StateArchive& archive) {
    // Canonical hex words retain signed zero; metadata binds mesh, bias and inputs.
    nlohmann::json j={{"nodes",archive.nodeCount},{"metadata",archive.metadata}};
    j["fields"]=nlohmann::json::object();
    for (const auto& [name,values]:archive.fields) {
        auto& a=j["fields"][name]; a=nlohmann::json::array();
        for (double v:values) a.push_back(bits(v));
    }
    return stateSha256(j.dump());
}
inline double physical(double value, const UnitScalingConfig& scaling) {
    double si=scaling.unitSystem().internalConcentrationToM3(value);
    if(si!=0. && std::abs(si)<std::numeric_limits<double>::min()) si=0.;
    return si;
}
} // namespace population_checkpoint_detail

// Private sidecar, deliberately separate from the public SI archive schema.
// This preserves population storage only, not solver caches/history or all DDSolution fields.
inline nlohmann::json populationCheckpoint(const DDSolution& solution,
    const StateArchive& archive, UnitScalingConfig scaling = {}) {
    using namespace population_checkpoint_detail;
    validateStateArchive(archive);
    if(archive.metadata.at("mode")!="dd" ||
       solution.n.size()!=static_cast<Eigen::Index>(archive.nodeCount) ||
       solution.p.size()!=static_cast<Eigen::Index>(archive.nodeCount))
        throw std::invalid_argument("Population checkpoint shape/mode mismatch");
    // Bind the entire supplied state, not just its density arrays.
    if(archiveIdentity(archiveDDSolution(solution,archive.metadata,scaling))!=archiveIdentity(archive))
        throw std::invalid_argument("Population checkpoint state mismatch");
    nlohmann::json j={{"schema","vela.experimental.population-checkpoint.v1"},
        {"archive_sha256",archiveIdentity(archive)},
        {"concentration_scale_bits",bits(scaling.unitSystem().concentrationM3PerInternal())},
        {"length_scale_bits",bits(scaling.unitSystem().lengthMPerInternal())}};
    for (int block=0;block<2;++block) {
        const auto& v=block==0?solution.n:solution.p;
        auto& a=j[block==0?"electrons_bits":"holes_bits"]; a=nlohmann::json::array();
        for(double x:v) {
            if(!std::isfinite(x)||x<0.) throw std::invalid_argument("Invalid internal population");
            a.push_back(bits(x));
        }
    }
    j["payload_sha256"]=stateSha256(j.dump());
    return j;
}

inline DDSolution restorePopulationCheckpoint(const StateArchive& archive,
    nlohmann::json sidecar, UnitScalingConfig scaling = {}) {
    using namespace population_checkpoint_detail;
    const auto digest=sidecar.at("payload_sha256").get<std::string>();
    sidecar.erase("payload_sha256");
    if(stateSha256(sidecar.dump())!=digest ||
       sidecar.at("schema")!="vela.experimental.population-checkpoint.v1" ||
       sidecar.at("archive_sha256")!=archiveIdentity(archive) ||
       sidecar.at("concentration_scale_bits")!=bits(scaling.unitSystem().concentrationM3PerInternal()) ||
       sidecar.at("length_scale_bits")!=bits(scaling.unitSystem().lengthMPerInternal()))
        throw std::invalid_argument("Population checkpoint identity mismatch");
    auto result=restoreDDSolution(archive,scaling);
    for(int block=0;block<2;++block) {
        auto& v=block==0?result.n:result.p;
        const auto& a=sidecar.at(block==0?"electrons_bits":"holes_bits");
        const auto& si=archive.fields.at(block==0?"electrons_m3":"holes_m3");
        if(!a.is_array() || a.size()!=static_cast<std::size_t>(v.size()))
            throw std::invalid_argument("Population checkpoint count mismatch");
        for(Eigen::Index i=0;i<v.size();++i) {
            const double x=value(a[i]);
            if(!std::isfinite(x)||x<0. || bits(physical(x,scaling))!=bits(si[i]))
                throw std::invalid_argument("Population checkpoint SI mismatch");
            v[i]=x;
        }
    }
    return result;
}
} // namespace vela::experimental
