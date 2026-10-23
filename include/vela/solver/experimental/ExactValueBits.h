#pragma once
#include "vela/core/Types.h"
#include <boost/pfr/core.hpp>
#include <boost/pfr/core_name.hpp>
#include <bit>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

namespace vela::experimental {
namespace exact_bits_detail {
template<class> inline constexpr bool unsupported = false;
template<class> struct IsVector : std::false_type {};
template<class T, class A> struct IsVector<std::vector<T,A>> : std::true_type {};
inline void word(std::string& out, std::uint64_t value) {
    for (int i=0; i<8; ++i) out.push_back(static_cast<char>((value >> (8*i)) & 255));
}
inline void text(std::string& out, const std::string& value) {
    word(out,value.size()); out.append(value);
}
template<class T> void append(std::string& out, const T& value) {
    using U=std::remove_cvref_t<T>;
    if constexpr(std::is_same_v<U,double>) {
        out.push_back('d'); word(out,std::bit_cast<std::uint64_t>(value));
    } else if constexpr(std::is_same_v<U,bool>) {
        out.push_back('b'); out.push_back(value?'1':'0');
    } else if constexpr(std::is_integral_v<U>) {
        static_assert(sizeof(U)<=8);
        out.push_back(std::is_signed_v<U>?'i':'u'); out.push_back(sizeof(U));
        word(out,static_cast<std::uint64_t>(value));
    } else if constexpr(std::is_enum_v<U>) {
        out.push_back('e'); append(out,static_cast<std::underlying_type_t<U>>(value));
    } else if constexpr(std::is_same_v<U,std::string>) {
        out.push_back('s'); text(out,value);
    } else if constexpr(std::is_same_v<U,VectorXd>) {
        out.push_back('v'); word(out,value.size());
        for(Eigen::Index i=0;i<value.size();++i) word(out,std::bit_cast<std::uint64_t>(value[i]));
    } else if constexpr(std::is_same_v<U,SparseMatrixd>) {
        // Preserve shape, traversal order, stored zeros and coefficient bits.
        // Allocator capacity, pointers and padding are not numerical state.
        out.push_back('m'); word(out,value.rows()); word(out,value.cols());
        word(out,value.nonZeros()); out.push_back(value.isCompressed()?'1':'0');
        for(int k=0;k<value.outerSize();++k) {
            word(out,k);
            for(typename U::InnerIterator it(value,k);it;++it) {
                out.push_back('c'); word(out,it.row()); word(out,it.col());
                word(out,std::bit_cast<std::uint64_t>(it.value()));
            }
            out.push_back(';');
        }
    } else if constexpr(IsVector<U>::value) {
        out.push_back('a'); word(out,value.size());
        for(const auto& item:value) append(out,item);
    } else if constexpr(std::is_aggregate_v<U>) {
        out.push_back('r'); word(out,boost::pfr::tuple_size_v<U>);
        // Every declared field is visited, including fields added in the future.
        // Unsupported new field types fail compilation rather than being omitted.
        boost::pfr::for_each_field(value,[&](const auto& field,auto index) {
            text(out,std::string(boost::pfr::get_name<decltype(index)::value,U>()));
            append(out,field);
        });
    } else {
        static_assert(unsupported<U>,"Unsupported exact diagnostic field type");
    }
}
} // namespace exact_bits_detail

// Post-solve comparison encoding, not a persistence/restart format. Does not
// evaluate floating point expressions: signed zero, NaN payloads, infinities and
// subnormals retain their bit patterns. Never serialize object padding/pointers.
template<class T> std::string exactValueBits(const T& value) {
    std::string out="vela.experimental.logical-value-bits.v1\n";
    exact_bits_detail::append(out,value);
    return out;
}
} // namespace vela::experimental
