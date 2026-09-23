#include "DDStateFormats.h"
#include "dd_state_generated.h"
#include <highfive/highfive.hpp>
#include <bit>
#include <cstdint>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace vela {
namespace {
using Names=std::vector<std::string>;
constexpr auto units="V;m^-3";
[[noreturn]] void bad(const char* message) {throw std::runtime_error(std::string("DD state format: ")+message);}
std::uint64_t get(const std::string& b,std::size_t pos,unsigned width) {
    if(pos>b.size() || width>b.size()-pos) bad("truncated data");
    std::uint64_t v=0;for(unsigned i=0;i<width;++i) v|=std::uint64_t(static_cast<unsigned char>(b[pos+i]))<<(8*i);return v;
}
void put(std::string& b,std::uint64_t v,unsigned width) {for(unsigned i=0;i<width;++i)b.push_back(static_cast<char>((v>>(8*i))&255));}
Names names(unsigned flags) {
    if((flags&~7u) || ((flags&2u)&&!(flags&1u))) bad("unsupported fields");
    Names n{"psi","phin","phip","electrons_m3","holes_m3"};
    if(flags&1)n.emplace_back("electron_quantum_potential_V");
    if(flags&2)n.emplace_back("electron_quantum_potential_like_V");
    if(flags&4) for(const auto* s:{"electron_qf_increment_V","hole_qf_increment_V","electron_qf_reference_V","hole_qf_reference_V"}) n.emplace_back(s);
    return n;
}
unsigned flagsFor(const Names& n) {
    for(unsigned f:{0u,1u,3u,4u,5u,7u})if(names(f)==n)return f;
    bad("unsupported field names or ordering");
}
struct Table {std::uint64_t count=0;Names fields;std::vector<double> data;};
Table fromCanonical(const std::string& b) {
    Table t;t.count=get(b,16,8);t.fields=names(static_cast<unsigned>(get(b,12,4)));
    const auto size=t.count*t.fields.size();if(b.size()!=24+8*size)bad("invalid canonical data");
    t.data.reserve(size);for(std::size_t i=0;i<size;++i)t.data.push_back(std::bit_cast<double>(get(b,24+8*i,8)));return t;
}
std::string canonical(const Table& t,Index expected) {
    if(expected<0 || t.count!=static_cast<std::uint64_t>(expected) || t.data.size()!=t.count*t.fields.size())bad("node count or dimensions mismatch");
    const unsigned f=flagsFor(t.fields);std::string b="VELADS01";b.reserve(24+8*t.data.size());
    put(b,1,4);put(b,f,4);put(b,t.count,8);for(double x:t.data)put(b,std::bit_cast<std::uint64_t>(x),8);return b;
}
std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p,std::ios::binary|std::ios::ate);if(!in)bad("cannot open file");
    const auto size=in.tellg();if(size<0 || size>1024LL*1024*1024)bad("file size limit");
    std::string b(static_cast<std::size_t>(size),'\0');in.seekg(0);in.read(b.data(),size);if(!in)bad("read failed");return b;
}
void writeFile(const std::filesystem::path& p,const std::string& b) {
    std::ofstream out(p,std::ios::binary|std::ios::trunc);out.write(b.data(),b.size());out.close();if(!out)bad("write failed");
}
std::string npyHeader(const Table& t) {
    std::string h="{'descr': [";
    for(const auto& n:t.fields)h+="('"+n+"', '<f8'), ";
    h+="], 'fortran_order': False, 'shape': ("+std::to_string(t.count)+",), }";
    h.append((64-(10+h.size()+1)%64)%64,' ');h+='\n';return h;
}
Table readNpy(const std::string& b,Index expected) {
    if(b.size()<10 || b.substr(0,6)!=std::string("\x93NUMPY",6))bad("invalid NPY magic");
    const auto major=get(b,6,1);if((major!=1&&major!=2)||get(b,7,1)!=0)bad("unsupported NPY version");
    const auto offset=major==1?10u:12u;const auto length=get(b,8,major==1?2:4);
    if(length>65536 || offset+length>b.size())bad("invalid NPY header length");
    const auto h=b.substr(offset,length);std::smatch m;Table t;
    if(!std::regex_search(h,m,std::regex("'shape'\\s*:\\s*\\(\\s*([0-9]+)\\s*,\\s*\\)")))bad("NPY needs one-dimensional records");
    t.count=std::stoull(m[1]);if(expected<0||t.count!=static_cast<std::uint64_t>(expected))bad("NPY node count mismatch");
    if(!std::regex_search(h,std::regex("'fortran_order'\\s*:\\s*False")))bad("NPY record layout must be C order");
    if(!std::regex_search(h,m,std::regex("'descr'\\s*:\\s*\\[([^\\]]*)\\]")))bad("NPY needs named float64 fields");
    const std::string descr=m[1];const std::regex item("\\('([A-Za-z0-9_]+)'\\s*,\\s*'<f8'\\s*\\)");
    for(auto it=std::sregex_iterator(descr.begin(),descr.end(),item);it!=std::sregex_iterator();++it)t.fields.push_back((*it)[1]);
    if(!std::regex_match(std::regex_replace(descr,item,""),std::regex("[\\s,]*")))bad("unsupported NPY dtype");
    flagsFor(t.fields);const auto size=t.count*t.fields.size();
    if(b.size()!=offset+length+size*8)bad("NPY payload size mismatch");
    t.data.resize(size);
    for(std::size_t i=0;i<t.count;++i)for(std::size_t c=0;c<t.fields.size();++c)
        t.data[c*t.count+i]=std::bit_cast<double>(get(b,offset+length+8*(i*t.fields.size()+c),8));
    return t;
}
} // namespace

void writePublicDDState(const std::filesystem::path& p,const std::string& bytes) {
    const Table t=fromCanonical(bytes);
    if(p.extension()==".h5") {
        HighFive::File file(p.string(),HighFive::File::Truncate);
        file.createAttribute("schema",std::string("vela.ddstate/1"));file.createAttribute("units",std::string(units));
        file.createAttribute("node_count",t.count);file.createDataSet("field_names",t.fields);
        for(std::size_t c=0;c<t.fields.size();++c) {
            auto ds=file.createDataSet<double>(t.fields[c],HighFive::DataSpace(std::vector<std::size_t>{static_cast<std::size_t>(t.count)}));
            ds.write_raw(t.data.data()+c*t.count);
        }
        file.flush();return;
    }
    if(p.extension()==".vfb") {
        flatbuffers::FlatBufferBuilder b;std::vector<flatbuffers::Offset<flatbuffers::String>> n;
        for(const auto& name:t.fields)n.push_back(b.CreateString(name));
        const auto labels=b.CreateVector(n);const auto values=b.CreateVector(t.data);const auto u=b.CreateString(units);
        auto state=vela_state::CreateDDState(b,1,t.count,labels,values,u);vela_state::FinishDDStateBuffer(b,state);
        writeFile(p,std::string(reinterpret_cast<const char*>(b.GetBufferPointer()),b.GetSize()));return;
    }
    if(p.extension()==".npy") {
        const auto h=npyHeader(t);std::string b("\x93NUMPY\x01\x00",8);put(b,h.size(),2);b+=h;
        b.reserve(b.size()+8*t.data.size());
        for(std::size_t i=0;i<t.count;++i)for(std::size_t c=0;c<t.fields.size();++c)put(b,std::bit_cast<std::uint64_t>(t.data[c*t.count+i]),8);
        writeFile(p,b);return;
    }
    bad("unsupported extension");
}

std::string readPublicDDState(const std::filesystem::path& p,Index expected) {
    Table t;
    if(p.extension()==".h5") {
        HighFive::File f(p.string(),HighFive::File::ReadOnly);
        if(f.getAttribute("schema").read<std::string>()!="vela.ddstate/1" || f.getAttribute("units").read<std::string>()!=units)bad("HDF5 schema or units mismatch");
        t.count=f.getAttribute("node_count").read<std::uint64_t>();
        if(expected<0||t.count!=static_cast<std::uint64_t>(expected))bad("HDF5 node count mismatch");
        auto labels=f.getDataSet("field_names");if(labels.getElementCount()>11)bad("too many HDF5 fields");
        t.fields=labels.read<Names>();flagsFor(t.fields);t.data.resize(t.count*t.fields.size());
        for(std::size_t c=0;c<t.fields.size();++c) {
            auto ds=f.getDataSet(t.fields[c]);const auto type=ds.getDataType();
            if(ds.getDimensions()!=std::vector<std::size_t>{static_cast<std::size_t>(t.count)} ||
               H5Tget_class(type.getId())!=H5T_FLOAT || H5Tget_size(type.getId())!=8)bad("HDF5 field type or dimension mismatch");
            ds.read_raw(t.data.data()+c*t.count);
        }
    } else if(p.extension()==".vfb") {
        const auto b=readFile(p);flatbuffers::Verifier verifier(reinterpret_cast<const std::uint8_t*>(b.data()),b.size());
        if(!vela_state::VerifyDDStateBuffer(verifier))bad("invalid FlatBuffer");
        const auto* s=vela_state::GetDDState(b.data());
        if(s->version()!=1 || s->units()->str()!=units)bad("FlatBuffer schema or units mismatch");
        t.count=s->node_count();if(expected<0||t.count!=static_cast<std::uint64_t>(expected))bad("FlatBuffer node count mismatch");
        if(s->field_names()->size()>11)bad("too many FlatBuffer fields");
        for(const auto* n:*s->field_names())t.fields.push_back(n->str());flagsFor(t.fields);
        if(s->values()->size()!=t.count*t.fields.size())bad("FlatBuffer dimension mismatch");
        t.data.assign(s->values()->begin(),s->values()->end());
    } else if(p.extension()==".npy")t=readNpy(readFile(p),expected);
    else bad("unsupported extension");
    return canonical(t,expected);
}
} // namespace vela
