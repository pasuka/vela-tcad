#include "vela/io/StateIdentity.h"
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <fstream>

namespace vela {
std::string stateSha256(std::string_view input) {
    constexpr std::array<std::uint32_t,64> k{
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::array<std::uint32_t,8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                                0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    if (input.size() > std::numeric_limits<std::uint64_t>::max()/8)
        throw std::length_error("State identity input too large");
    std::vector<unsigned char> bytes(input.begin(),input.end());
    const auto bits=static_cast<std::uint64_t>(bytes.size())*8;
    bytes.push_back(0x80);
    while(bytes.size()%64!=56) bytes.push_back(0);
    for(int i=7;i>=0;--i) bytes.push_back(static_cast<unsigned char>(bits>>(8*i)));
    for(std::size_t offset=0;offset<bytes.size();offset+=64) {
        std::array<std::uint32_t,64> w{};
        for(int i=0;i<16;++i) for(int j=0;j<4;++j) w[i]=(w[i]<<8)|bytes[offset+4*i+j];
        for(int i=16;i<64;++i) {
            const auto x=w[i-15],y=w[i-2];
            w[i]=w[i-16]+(std::rotr(x,7)^std::rotr(x,18)^(x>>3))+w[i-7]+(std::rotr(y,17)^std::rotr(y,19)^(y>>10));
        }
        auto [a,b,c,d,e,f,g,z]=h;
        for(int i=0;i<64;++i) {
            const auto t1=z+(std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25))+((e&f)^((~e)&g))+k[i]+w[i];
            const auto t2=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
            z=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=z;
    }
    constexpr char hex[]="0123456789abcdef";
    std::string out;
    for(auto v:h) for(int shift=28;shift>=0;shift-=4) out+=hex[(v>>shift)&15];
    return out;
}

std::string stateMeshIdentity(const DeviceMesh& mesh, UnitScalingConfig scaling) {
    return stateMeshIdentity(mesh, scaling.unitSystem().lengthMPerInternal());
}
std::string stateMeshIdentity(const DeviceMesh& mesh, double lengthMPerInternal) {
    if (!(std::isfinite(lengthMPerInternal) && lengthMPerInternal > 0.))
        throw std::invalid_argument("State mesh length unit must be positive and finite");
    std::string bytes="vela.mesh/1";
    auto integer=[&](std::uint64_t v) {for(int i=0;i<8;++i)bytes.push_back(static_cast<char>(v>>(i*8)));};
    auto real=[&](double v) {integer(std::bit_cast<std::uint64_t>(v));};
    auto text=[&](const std::string& value) {integer(value.size());bytes+=value;};
    auto indices=[&](const auto& values) {integer(values.size());for(auto v:values)integer(v);};
    real(lengthMPerInternal);
    integer(mesh.numNodes());
    for(const auto& n:mesh.nodes()){integer(n.id);real(n.x);real(n.y);}
    integer(mesh.numCells());
    for(const auto& c:mesh.cells()){integer(c.id);integer(c.region_id);indices(c.node_ids);}
    integer(mesh.numRegions());
    for(const auto& r:mesh.regions()){integer(r.id);text(r.name);text(r.material);indices(r.cell_ids);}
    integer(mesh.numContacts());
    for(const auto& c:mesh.contacts()){
        integer(c.id);text(c.name);integer(c.region_id);indices(c.node_ids);
        integer(c.edge_node_ids.size());for(const auto& edge:c.edge_node_ids){integer(edge[0]);integer(edge[1]);}
    }
    return stateSha256(bytes);
}
nlohmann::json stateInputProvenance(const nlohmann::json& config,
                                  const std::filesystem::path& directory) {
    using json = nlohmann::json;
    json files = json::object();
    const auto visit = [&](auto&& self, const json& value, const std::string& key) -> void {
        if (value.is_object()) {
            for (const auto& [name, child] : value.items()) self(self, child, key + "/" + name);
        } else if (value.is_array()) {
            for (std::size_t i=0; i<value.size(); ++i) self(self, value[i], key + "/" + std::to_string(i));
        } else if (key.ends_with("_file") && value.is_string()) {
            std::filesystem::path path = value.get<std::string>();
            if (path.is_relative()) path = directory/path;
            std::ifstream file(path, std::ios::binary);
            if (!file) throw std::runtime_error("Cannot fingerprint state input: " + path.string());
            files[key] = stateSha256(std::string(std::istreambuf_iterator<char>(file), {}));
        }
    };
    // Deliberately excludes output paths, diagnostics and restart/history files.
    for (const char* key : {"mesh_file", "materials_file", "node_doping_file", "doping",
                            "mesh_geometry", "mobility", "mobility_SI"})
        if (config.contains(key)) visit(visit, config.at(key), key);
    if (config.contains("solver") && config.at("solver").contains("mobility"))
        visit(visit, config.at("solver").at("mobility"), "solver/mobility");
    return files;
}
} // namespace vela
