#include "vela/io/ElectrothermalState.h"
#include "vela/io/StateIdentity.h"
#include <array>
#include <fstream>
#include <set>
#include <stdexcept>

namespace vela {
namespace {
using json = nlohmann::json;
namespace fs = std::filesystem;
constexpr std::array<const char*, 5> keys{"state_interleaved", "referenced_state_interleaved",
    "electron_qf_reference_V", "hole_qf_reference_V", "temperature_K"};
std::string digest(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot read state archive " + p.string());
    return stateSha256(std::string(std::istreambuf_iterator<char>(f), {}));
}
StateArchive fromRecord(const json& value, json metadata) {
    const bool referenced = value.contains("referenced_state_interleaved");
    const auto x = value.at(referenced ? "referenced_state_interleaved" : "state_interleaved").get<std::vector<double>>();
    if (x.empty() || x.size()%4) throw std::invalid_argument("Invalid four-equation layout");
    StateArchive s; s.nodeCount=x.size()/4; s.metadata=std::move(metadata);
    s.metadata["mode"]="electrothermal";
    s.metadata["potential_origin_V"]=value.value("potential_origin_V",s.metadata.value("potential_origin_V",0.));
    s.metadata["record_fields"]=json::array();
    for (const auto* k: keys) if (value.contains(k)) s.metadata["record_fields"].push_back(k);
    const auto physical=value.contains("state_interleaved") ? value.at("state_interleaved").get<std::vector<double>>() : std::vector<double>{};
    if (!physical.empty() && physical.size()!=x.size()) throw std::invalid_argument("Incomplete physical state");
    for (const auto* k:{"psi","phin","phip","temperature_K"}) s.fields[k].resize(s.nodeCount);
    for (const auto* carrier:{"electron","hole"}) {
        const std::string ref=std::string(carrier)+"_qf_reference_V";
        const std::string inc=std::string(carrier)+"_qf_increment_V";
        if (!value.contains(ref)) throw std::invalid_argument("Explicit thermal QF references required");
        s.fields[ref]=value.at(ref).get<std::vector<double>>();
        if (s.fields[ref].size()!=s.nodeCount) throw std::invalid_argument("Incomplete thermal QF reference");
        s.fields[inc].resize(s.nodeCount);
    }
    for (std::size_t i=0;i<s.nodeCount;++i) {
        s.fields["psi"][i]=x[4*i]; s.fields["temperature_K"][i]=x[4*i+3];
        if (!physical.empty() && (physical[4*i]!=x[4*i] || physical[4*i+3]!=x[4*i+3]))
            throw std::invalid_argument("Inconsistent thermal psi or temperature representations");
        for (int k:{1,2}) {
            const std::string carrier=k==1?"electron":"hole", field=k==1?"phin":"phip";
            const double ref=s.fields.at(carrier+"_qf_reference_V")[i];
            s.fields[carrier+"_qf_increment_V"][i]=referenced?x[4*i+k]:x[4*i+k]-ref;
            s.fields[field][i]=physical.empty()?x[4*i+k]+ref:physical[4*i+k];
        }
    }
    if (value.contains("temperature_K") && value.at("temperature_K").get<std::vector<double>>()!=s.fields.at("temperature_K"))
        throw std::invalid_argument("Inconsistent temperature field");
    validateStateArchive(s); return s;
}
json pack(const fs::path& p, const json& record, const json& metadata, unsigned& sequence) {
    json out=record;
    if (out.contains("state_archive")) throw std::invalid_argument("State record already packed");
    if (out.contains("state_interleaved") || out.contains("referenced_state_interleaved")) {
        auto state=fromRecord(out,metadata);
        const auto path=p.parent_path()/(p.stem().string()+"_state_"+std::to_string(sequence++)+".h5");
        if (fs::exists(path)) throw std::runtime_error("Refusing to overwrite referenced state " + path.string());
        writeStateArchive(path,state);
        for (const auto* k:keys) out.erase(k);
        out["state_archive"]={{"file",path.filename().string()},{"sha256",digest(path)}};
    }
    if (out.contains("diagnostic_tangent_predictor"))
        out["diagnostic_tangent_predictor"]=pack(p,out.at("diagnostic_tangent_predictor"),metadata,sequence);
    if (out.contains("diagnostic_predictor_candidates"))
        for (auto& item:out.at("diagnostic_predictor_candidates")) item=pack(p,item,metadata,sequence);
    return out;
}
json unpack(const fs::path& p,const json& record,std::uint64_t nodes,const std::string& identity,double origin) {
    json out=record;
    if (out.contains("state_archive")) {
        for (const auto* k:keys) if (out.contains(k)) throw std::invalid_argument("Mixed inline and HDF5 thermal state");
        const auto ref=out.at("state_archive");
        const fs::path relative=ref.at("file").get<std::string>();
        const auto path=relative.is_absolute()?relative:p.parent_path()/relative;
        if (digest(path)!=ref.at("sha256").get<std::string>()) throw std::invalid_argument("Thermal state file digest mismatch");
        const auto s=readStateArchive(path,nodes,identity);
        if(s.metadata.at("mode")!="electrothermal" || s.metadata.at("potential_origin_V").get<double>()!=origin)
            throw std::invalid_argument("Thermal state mode or reference origin mismatch");
        const auto saved=s.metadata.at("record_fields").get<std::set<std::string>>();
        std::vector<double> physical(4*nodes),increment(4*nodes);
        for(std::size_t i=0;i<nodes;++i) {
            for(int k:{0,3})physical[4*i+k]=increment[4*i+k]=s.fields.at(k==0?"psi":"temperature_K")[i];
            for(int k:{1,2}) {
                physical[4*i+k]=s.fields.at(k==1?"phin":"phip")[i];
                increment[4*i+k]=s.fields.at(k==1?"electron_qf_increment_V":"hole_qf_increment_V")[i];
            }
        }
        for(const auto& k:saved) {
            if(k=="state_interleaved")out[k]=physical;
            else if(k=="referenced_state_interleaved")out[k]=increment;
            else if(k=="temperature_K"||k=="electron_qf_reference_V"||k=="hole_qf_reference_V")out[k]=s.fields.at(k);
            else throw std::invalid_argument("Unknown saved thermal field");
        }
        out.erase("state_archive");
    } else if(out.contains("state_interleaved")||out.contains("referenced_state_interleaved"))
        throw std::invalid_argument("Production thermal state requires HDF5; convert the input explicitly");
    if(out.contains("diagnostic_tangent_predictor"))
        out["diagnostic_tangent_predictor"]=unpack(p,out.at("diagnostic_tangent_predictor"),nodes,identity,origin);
    if(out.contains("diagnostic_predictor_candidates"))
        for(auto& item:out.at("diagnostic_predictor_candidates")) item=unpack(p,item,nodes,identity,origin);
    return out;
}
} // namespace
nlohmann::json packElectrothermalRecord(const fs::path& p,const json& value,const json& metadata) {
    unsigned sequence=0;return pack(p,value,metadata,sequence);
}
nlohmann::json unpackElectrothermalRecord(const fs::path& p,const json& value,std::uint64_t nodes,const std::string& sha,double origin) {
    return unpack(p,value,nodes,sha,origin);
}
} // namespace vela
