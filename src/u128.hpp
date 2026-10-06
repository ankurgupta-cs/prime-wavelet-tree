#pragma once
// u128.hpp -- minimal 128-bit unsigned helpers for the Carmichael project.
// Requires GCC or Clang on a 64-bit target (unsigned __int128). Not MSVC.

#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace cn {

using u32  = std::uint32_t;
using u64  = std::uint64_t;
using u128 = unsigned __int128;

inline constexpr u128 U128_MAX = ~static_cast<u128>(0);

// Decimal conversion. Values at the 10^24 scale of the CN table round-trip
// exactly; anything up to 2^128 - 1 (39 digits) is supported.
inline std::string to_string(u128 x) {
    if (x == 0) return "0";
    char buf[40];
    int i = 40;
    while (x > 0) {
        buf[--i] = static_cast<char>('0' + static_cast<unsigned>(x % 10));
        x /= 10;
    }
    return std::string(buf + i, buf + 40);
}

// Strict decimal parse: digits only; throws on empty input, non-digit, or
// overflow past 2^128 - 1.
inline u128 parse_u128(std::string_view s) {
    if (s.empty()) throw std::invalid_argument("parse_u128: empty string");
    u128 v = 0;
    for (char c : s) {
        if (c < '0' || c > '9') throw std::invalid_argument("parse_u128: non-digit");
        const unsigned d = static_cast<unsigned>(c - '0');
        if (v > (U128_MAX - d) / 10) throw std::out_of_range("parse_u128: overflow");
        v = v * 10 + d;
    }
    return v;
}

inline std::ostream& operator<<(std::ostream& os, u128 x) { return os << to_string(x); }

} // namespace cn
