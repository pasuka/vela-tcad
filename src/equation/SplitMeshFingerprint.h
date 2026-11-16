#pragma once

#include <boost/version.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <span>
#include <sstream>
#include <string>
#include <string_view>

namespace vela::detail {
namespace mesh_sha1 {

inline constexpr std::size_t blockSize = 64;
using State = std::array<std::uint32_t, 5>;

// SHA-1 compression as specified in RFC 3174, section 6.1.
inline void compress(State& state, std::span<const std::byte, blockSize> block)
{
    std::array<std::uint32_t, 80> words;
    for (std::size_t i = 0; i < 16; ++i) {
        words[i] = 0;
        for (const auto byte : block.subspan(i * 4, 4))
            words[i] = (words[i] << 8) | std::to_integer<std::uint32_t>(byte);
    }
    for (std::size_t i = 16; i < words.size(); ++i)
        words[i] = std::rotl(words[i - 3] ^ words[i - 8] ^ words[i - 14] ^ words[i - 16], 1);

    auto [a, b, c, d, e] = state;
    for (std::size_t i = 0; i < words.size(); ++i) {
        std::uint32_t mix, constant;
        if (i < 20) {
            mix = (b & c) | (~b & d);
            constant = 0x5a827999;
        } else if (i < 40) {
            mix = b ^ c ^ d;
            constant = 0x6ed9eba1;
        } else if (i < 60) {
            mix = (b & c) | (b & d) | (c & d);
            constant = 0x8f1bbcdc;
        } else {
            mix = b ^ c ^ d;
            constant = 0xca62c1d6;
        }
        const auto next = std::rotl(a, 5) + mix + e + constant + words[i];
        e = d;
        d = c;
        c = std::rotl(b, 30);
        b = a;
        a = next;
    }
    const State working{a, b, c, d, e};
    for (std::size_t i = 0; i < state.size(); ++i)
        state[i] += working[i];
}

} // namespace mesh_sha1

inline std::string splitMeshFingerprint(std::string_view geometry)
{
    using mesh_sha1::blockSize;
    mesh_sha1::State state{0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0};
    auto remaining = std::as_bytes(std::span{geometry.data(), geometry.size()});
    while (remaining.size() >= blockSize) {
        mesh_sha1::compress(state, remaining.first<blockSize>());
        remaining = remaining.subspan<blockSize>();
    }

    // The tail is shorter than one block. RFC 3174 puts its marker and the
    // original 64-bit bit count in one or two zero-initialized blocks.
    constexpr std::size_t lengthBytes = 8;
    std::array<std::byte, 2 * blockSize> tail{};
    std::ranges::copy(remaining, tail.begin());
    tail[remaining.size()] = std::byte{0x80};
    const bool needsSecondBlock = remaining.size() >= blockSize - lengthBytes;
    const auto padding = std::span{tail};
    const auto finalBlock = needsSecondBlock ? padding.last<blockSize>() : padding.first<blockSize>();
    const auto length = finalBlock.last<lengthBytes>();
    auto bitCount = static_cast<std::uint64_t>(geometry.size()) * 8;
    for (std::size_t i = length.size(); i != 0; --i) {
        length[i - 1] = static_cast<std::byte>(bitCount & 0xff);
        bitCount >>= 8;
    }
    if (needsSecondBlock)
        mesh_sha1::compress(state, padding.first<blockSize>());
    mesh_sha1::compress(state, finalBlock);

    // Preserve checkpoint identifiers on Boost 1.84/1.85. The former UUID
    // helper returned five words there and 20 bytes starting in Boost 1.86.
    constexpr int hexWidth = BOOST_VERSION < 108600 ? 2 : 8;
    std::ostringstream fingerprint;
    fingerprint << "sha1:" << std::hex << std::setfill('0');
    for (const auto word : state)
        fingerprint << std::setw(hexWidth) << word;
    return fingerprint.str();
}

} // namespace vela::detail
