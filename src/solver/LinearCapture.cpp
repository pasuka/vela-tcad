#include "LinearCapture.h"
#include <atomic>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
namespace vela::detail {
void captureLinearInput(const std::string& directory,const SparseMatrixd& a,const VectorXd& b,const VectorXd& x) {
    if(directory.empty() || a.rows()<30000 || a.rows()%3!=0) return;
    static std::atomic<unsigned> count=0;
    const auto index=count.fetch_add(1);if(index>=4)return;
    static_assert(sizeof(double)==8 && sizeof(SparseMatrixd::StorageIndex)==4 && std::endian::native==std::endian::little);
    const auto path=std::filesystem::path(directory)/("input_"+std::to_string(index)+".bin");
    if(std::filesystem::exists(path))throw std::runtime_error("Refusing to overwrite linear capture");
    std::ofstream out(path,std::ios::binary);out.exceptions(std::ios::badbit|std::ios::failbit);
    const auto bytes=[&](const void* p,std::size_t n){out.write(static_cast<const char*>(p),static_cast<std::streamsize>(n));};
    const auto size=[&](std::uint64_t n){bytes(&n,8);};
    const auto vector=[&](const VectorXd& v){size(v.size());bytes(v.data(),8*v.size());};
    bytes("VELALU02",8);size(index);size(index);
    size(a.rows());size(a.cols());size(a.nonZeros());
    bytes(a.outerIndexPtr(),4*(a.cols()+1));bytes(a.innerIndexPtr(),4*a.nonZeros());bytes(a.valuePtr(),8*a.nonZeros());
    vector(b);vector(x);
}
}
