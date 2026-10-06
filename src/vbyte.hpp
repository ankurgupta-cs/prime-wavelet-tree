#pragma once
// vbyte.hpp -- byte-aligned variable-length integer coding for ORC1.
//
// Little-endian base-128: each byte carries 7 payload bits (low bits first);
// the high bit is a continuation flag. Values up to u128 are supported so the
// same codec covers factors (< 2^40), deltas, and any future n-typed fields.
// Layout follows Plaisance-Kurz-Lemire (arXiv 1503.07387); scalar codec now,
// SIMD later if query latency ever asks for it.

#include <cstdint>
#include <vector>

#include "u128.hpp"

namespace cn {

inline constexpr int kVByteMax128 = 19;   // ceil(128 / 7) bytes worst case

inline void vbyte_append(std::vector<std::uint8_t>& out, u128 v) {
    while (v > 0x7f) {
        out.push_back(static_cast<std::uint8_t>(v) | 0x80);
        v >>= 7;
    }
    out.push_back(static_cast<std::uint8_t>(v));
}

inline int vbyte_length(u128 v) {
    int n = 1;
    while (v > 0x7f) { v >>= 7; ++n; }
    return n;
}

// Decodes one value from [p, end). Returns the byte position after the value,
// or nullptr on malformed input: truncated, more than 19 bytes, or a 19-byte
// encoding whose payload exceeds 128 bits (the final byte may carry only the
// low 2 of its 7 payload bits). Non-canonical zero-padded encodings are
// accepted; the writer never emits them and sha_body covers integrity.
inline const std::uint8_t* vbyte_decode(const std::uint8_t* p,
                                        const std::uint8_t* end, u128& v) {
    v = 0;
    int shift = 0;
    while (p < end && shift < kVByteMax128 * 7) {
        const std::uint8_t b = *p++;
        if (shift == 126 && (b & 0x7c) != 0) return nullptr;   // bits >= 2^128
        v |= static_cast<u128>(b & 0x7f) << shift;
        if (!(b & 0x80)) return p;
        shift += 7;
    }
    return nullptr;
}

// u64 convenience: rejects values that do not fit.
inline const std::uint8_t* vbyte_decode_u64(const std::uint8_t* p,
                                            const std::uint8_t* end, u64& v) {
    u128 wide;
    p = vbyte_decode(p, end, wide);
    if (p == nullptr || wide > UINT64_MAX) return nullptr;
    v = static_cast<u64>(wide);
    return p;
}

} // namespace cn
