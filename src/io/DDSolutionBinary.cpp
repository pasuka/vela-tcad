#include "vela/io/DDSolutionCsv.h"
#ifdef VELA_ENABLE_HDF5_STATE
#include "vela/io/DDSolutionState.h"
#endif
#ifdef VELA_ENABLE_STATE_FORMATS
#include "DDStateFormats.h"
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace vela {
namespace {
constexpr std::array<char, 8> magic{'V','E','L','A','D','S','0','1'};
thread_local const ScopedDDStateInput* stateInput = nullptr;
static_assert(sizeof(Real) == 8 && std::numeric_limits<Real>::is_iec559);

[[noreturn]] void invalid(const char* why) {
    throw std::runtime_error(std::string("VDS1 restart state: ") + why);
}
void put(std::ostream& out, std::uint64_t value, int bytes) {
    for (int i=0;i<bytes;++i) out.put(static_cast<char>((value >> (8*i)) & 255));
}
std::uint64_t get(std::istream& in, int bytes) {
    std::uint64_t value=0;
    for (int i=0;i<bytes;++i) {
        const int c=in.get(); if (c==std::char_traits<char>::eof()) invalid("truncated file");
        value |= static_cast<std::uint64_t>(static_cast<unsigned char>(c)) << (8*i);
    }
    return value;
}
void checkQf(const DDSolution& s) {
    for (int i=0;i<s.psi.size();++i) {
        const Real n=s.electronQfReference(i)+s.phinIncrement(i);
        const Real p=s.holeQfReference(i)+s.phipIncrement(i);
        auto valid=[](Real a,Real b) {
            return std::isfinite(b) && std::abs(a-b)<=32*std::numeric_limits<Real>::epsilon()*
                std::max({Real{1},std::abs(a),std::abs(b)});
        };
        if (!valid(s.phin(i),n)||!valid(s.phip(i),p)) invalid("inconsistent referenced quasi-Fermi coordinates");
    }
}
DDSolution readBinary(std::istream& in, Index expected, UnitScalingConfig scaling) {
    if (get(in,4)!=1) invalid("unsupported version");
    const auto flags=get(in,4),count=get(in,8);
    if ((flags & ~std::uint64_t{7}) || ((flags&2)&&!(flags&1))) invalid("unsupported fields");
    if (count!=static_cast<std::uint64_t>(expected) || count>static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        invalid("node count mismatch");
    const auto fields=5+((flags&1)?1:0)+((flags&2)?1:0)+((flags&4)?4:0);
    const auto start=in.tellg();in.seekg(0,std::ios::end);const auto end=in.tellg();in.seekg(start);
    if (start<0||end<start||static_cast<std::uint64_t>(end-start)!=count*8*fields)
        invalid("payload length mismatch");
    DDSolution s;
    std::vector<VectorXd*> columns{&s.psi,&s.phin,&s.phip,&s.n,&s.p};
    s.electronQuantumPotential=VectorXd::Zero(static_cast<int>(count));
    if(flags&1) columns.push_back(&s.electronQuantumPotential);
    if(flags&2) columns.push_back(&s.electronQuantumPotentialLike);
    if(flags&4) {
        columns.insert(columns.end(),{&s.phinIncrement,&s.phipIncrement,&s.electronQfReference,&s.holeQfReference});
    }
    for(auto* column:columns) {
        column->resize(static_cast<int>(count));
        if constexpr (std::endian::native==std::endian::little) {
            in.read(reinterpret_cast<char*>(column->data()),column->size()*sizeof(Real));
            if(!in) invalid("truncated field");
        }
        for(int i=0;i<column->size();++i) {
            const Real v=std::endian::native==std::endian::little ? (*column)(i) : std::bit_cast<Real>(get(in,8));
            if(!std::isfinite(v)) invalid("non-finite field");
            (*column)(i)=v;
        }
    }
    const auto& units=scaling.unitSystem();
    for(int i=0;i<s.psi.size();++i) {
        s.n(i)=units.m3ToInternalConcentration(s.n(i));
        s.p(i)=units.m3ToInternalConcentration(s.p(i));
        if(!std::isfinite(s.n(i))||!std::isfinite(s.p(i))) invalid("density unit conversion overflow");
    }
    if(flags&4) {
        checkQf(s);
        if(count) {s.electronQfReference_V=s.electronQfReference(0);s.holeQfReference_V=s.holeQfReference(0);}
    }
    s.iters=0;s.converged=true;
    return s;
}
void writeBinary(std::ostream& out,const DDSolution& s,UnitScalingConfig scaling) {
    const auto count=s.psi.size();
    std::vector<VectorXd> columns{s.psi,s.phin,s.phip,s.n,s.p};
    for(const auto& c:columns) if(c.size()!=count) invalid("inconsistent field sizes");
    std::uint32_t flags=0;
    auto optional=[&](const VectorXd& c,std::uint32_t flag) {
        if(c.size()!=0 && c.size()!=count) invalid("partial optional field");
        if(c.size()==count && count>0) {flags|=flag;columns.push_back(c);}
    };
    optional(s.electronQuantumPotential,1);
    if(s.electronQuantumPotentialLike.size() && !(flags&1)) invalid("quantum-like field without quantum field");
    optional(s.electronQuantumPotentialLike,2);
    const bool qf=s.phinIncrement.size()!=0||s.phipIncrement.size()!=0;
    if(qf) {
        if(s.phinIncrement.size()!=count||s.phipIncrement.size()!=count) invalid("partial referenced coordinates");
        auto ref=[&](const VectorXd& c,Real scalar)->VectorXd {
            if(c.size()==count) return c;
            if(c.size()!=0) invalid("partial reference field");
            return VectorXd::Constant(count,scalar);
        };
        flags|=4;
        columns.push_back(s.phinIncrement);columns.push_back(s.phipIncrement);
        columns.push_back(ref(s.electronQfReference,s.electronQfReference_V));
        columns.push_back(ref(s.holeQfReference,s.holeQfReference_V));
    } else if(s.electronQfReference.size()||s.holeQfReference.size()) invalid("reference without increment");
    const auto& units=scaling.unitSystem();
    for(int i=0;i<count;++i) {
        columns[3](i)=units.internalConcentrationToM3(s.n(i));
        columns[4](i)=units.internalConcentrationToM3(s.p(i));
    }
    // Match the existing CSV writer's explicit subnormal-to-zero restart policy.
    for(auto& c:columns) for(auto& v:c) {
        if(!std::isfinite(v)) invalid("non-finite field");
        if(v!=0 && std::abs(v)<std::numeric_limits<Real>::min()) v=0;
    }
    if(qf) {
        DDSolution check=s;
        check.phin=columns[1];check.phip=columns[2];
        check.phinIncrement=columns[columns.size()-4];check.phipIncrement=columns[columns.size()-3];
        check.electronQfReference=columns[columns.size()-2];check.holeQfReference=columns.back();checkQf(check);
    }
    out.write(magic.data(),magic.size());put(out,1,4);put(out,flags,4);put(out,count,8);
    for(const auto& c:columns) {
        if constexpr (std::endian::native==std::endian::little)
            out.write(reinterpret_cast<const char*>(c.data()),c.size()*sizeof(Real));
        else for(Real v:c) put(out,std::bit_cast<std::uint64_t>(v),8);
    }
    if(!out) invalid("write failed");
}
} // namespace

ScopedDDStateInput::ScopedDDStateInput(std::filesystem::path name, std::string data)
    : path(name.lexically_normal()), bytes(std::move(data)) {
    if (stateInput) invalid("nested memory seed scope");
    if (!path.is_absolute() || path.extension()!=".vds") invalid("memory seed needs absolute .vds path");
    if (bytes.size()<24 || bytes.size()>16*1024*1024 || !std::equal(magic.begin(),magic.end(),bytes.begin()))
        invalid("invalid or oversized memory seed");
    stateInput=this;
}
ScopedDDStateInput::~ScopedDDStateInput() { stateInput=nullptr; }

DDSolution readDDSolutionState(const std::filesystem::path& path,Index expected,UnitScalingConfig scaling) {
#ifdef VELA_ENABLE_HDF5_STATE
    if (const auto* metadata = activeDDStateArchiveMetadata()) {
        const auto state = readStateArchive(path, expected, metadata->at("mesh_sha256").get<std::string>());
        if (state.metadata.at("potential_origin_V") != metadata->at("potential_origin_V"))
            throw std::invalid_argument("HDF5 state reference origin differs; translate the seed explicitly");
        return restoreDDSolution(state, scaling);
    }
#endif
    if(path.extension()==".h5" || path.extension()==".vfb" || path.extension()==".npy") {
#ifdef VELA_ENABLE_STATE_FORMATS
        std::istringstream in(readPublicDDState(path,expected),std::ios::in|std::ios::binary);
        in.seekg(magic.size());return readBinary(in,expected,scaling);
#else
        invalid("public state formats were not enabled at build time");
#endif
    }
    if (stateInput && std::filesystem::absolute(path).lexically_normal()==stateInput->path) {
        std::istringstream in(stateInput->bytes,std::ios::in|std::ios::binary);
        in.seekg(magic.size());
        return readBinary(in,expected,scaling);
    }
    std::ifstream in(path,std::ios::binary);std::array<char,8> header{};
    if(!in) invalid("cannot open input");
    in.read(header.data(),header.size());
    if(header==magic) return readBinary(in,expected,scaling);
    if(path.extension()==".vds") invalid("invalid magic");
    return readDDSolutionStateCsv(path,expected,scaling);
}
void writeDDSolutionState(const std::filesystem::path& path,const DDSolution& s,UnitScalingConfig scaling) {
#ifdef VELA_ENABLE_HDF5_STATE
    if (const auto* metadata = activeDDStateArchiveMetadata()) {
        writeStateArchive(path, archiveDDSolution(s, *metadata, scaling));
        return;
    }
#endif
    if(path.extension()==".h5" || path.extension()==".vfb" || path.extension()==".npy") {
#ifdef VELA_ENABLE_STATE_FORMATS
        std::ostringstream out(std::ios::out|std::ios::binary);writeBinary(out,s,scaling);
        if(!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
        writePublicDDState(path,out.str());return;
#else
        invalid("public state formats were not enabled at build time");
#endif
    }
    if(path.extension()!=".vds") {writeDDSolutionStateCsv(path,s,scaling);return;}
    if(!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path,std::ios::binary|std::ios::trunc);
    if(!out) invalid("cannot open output");
    writeBinary(out,s,scaling);out.close();if(!out) invalid("close failed");
}
} // namespace vela
