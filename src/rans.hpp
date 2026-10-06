#pragma once
// rans.hpp -- static multi-context rANS (range asymmetric numeral systems),
// the entropy coder of the PWT1 container (docs/PWT_FORMAT.md).
//
// Construction (generalized to many contexts and larger alphabets):
//   * 32-bit state x in [L, L << 8), L = 2^23, byte-wise renormalization
//     (the "rans_byte" construction of Giesen).
//   * one static frequency table per CONTEXT, quantized to kProbBits = 16
//     bits (total 65536) by largest-remainder rounding; every symbol that
//     occurs keeps frequency >= 1.
//   * symbols are coded in BLOCKS: the encoder buffers (context, symbol)
//     pairs, codes a block in reverse (rANS is LIFO) and the decoder reads
//     it forward. A block ends with the state back at L.
// Alphabets: up to 4096 symbols per context (the PWT count alphabet).

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "u128.hpp"
#include "vbyte.hpp"

namespace cn {

using u16 = std::uint16_t;

constexpr u32 kRansProbBits = 16;
constexpr u32 kRansProbScale = 1u << kRansProbBits;
constexpr u32 kRansL = 1u << 23;

// One context's quantized model. `syms` lists the used symbols ascending,
// freq[i] / cum[i] are parallel; cum has one extra trailing entry (= scale).
struct RansContext {
    std::vector<u16> syms;
    std::vector<u32> freq;
    std::vector<u32> cum;
    std::vector<u16> lut;   // lut[b] = index of the symbol whose range holds slot b * 256 (decode speed)

    bool empty() const { return syms.empty(); }

    // Position of symbol s in syms, or -1.
    int index_of(u16 s) const {
        const auto it = std::lower_bound(syms.begin(), syms.end(), s);
        return (it != syms.end() && *it == s) ? int(it - syms.begin()) : -1;
    }

    // Largest-remainder quantization of raw counts (length = alphabet size).
    void build(const u64* counts, std::size_t alphabet) {
        syms.clear(); freq.clear(); cum.clear();
        u64 total = 0;
        for (std::size_t s = 0; s < alphabet; ++s) total += counts[s];
        if (total == 0) return;
        std::vector<std::pair<double, std::size_t>> rema;
        u64 assigned = 0;
        for (std::size_t s = 0; s < alphabet; ++s) {
            if (!counts[s]) continue;
            const double exact = double(counts[s]) / double(total) * kRansProbScale;
            u32 f = static_cast<u32>(exact);
            if (f == 0) f = 1;
            syms.push_back(static_cast<u16>(s));
            freq.push_back(f);
            assigned += f;
            rema.push_back({exact - double(f), freq.size() - 1});
        }
        if (syms.size() > kRansProbScale) throw std::runtime_error("rans: alphabet too large for the scale");
        long long excess = static_cast<long long>(assigned) - static_cast<long long>(kRansProbScale);
        std::sort(rema.begin(), rema.end());
        if (excess < 0)
            for (std::size_t i = rema.size(); i-- > 0 && excess;) { ++freq[rema[i].second]; ++excess; }
        while (excess > 0) {   // shave from the largest frequencies, never below 1
            std::size_t best = 0;
            for (std::size_t i = 1; i < freq.size(); ++i) if (freq[i] > freq[best]) best = i;
            if (freq[best] <= 1) throw std::runtime_error("rans: cannot normalize");
            --freq[best];
            --excess;
        }
        finish_cum_();
    }

    // Serialized form: VByte(m), then m x (VByte(symbol delta), VByte(freq - 1)).
    void serialize(std::vector<std::uint8_t>& out) const {
        vbyte_append(out, syms.size());
        u16 prev = 0;
        for (std::size_t i = 0; i < syms.size(); ++i) {
            vbyte_append(out, static_cast<u64>(syms[i] - prev));
            vbyte_append(out, freq[i] - 1);
            prev = syms[i];
        }
    }
    // Returns bytes consumed; throws on a malformed table. strict: every VByte
    // must be canonical (PWT2 context sets 2/3; set 1 and PWT1 stay lax).
    std::size_t deserialize(const std::uint8_t* p, std::size_t len, std::size_t alphabet, bool strict = false) {
        syms.clear(); freq.clear(); cum.clear();
        std::size_t pos = 0;
        const u64 m = vbyte_read(p, len, pos, strict);
        if (m > alphabet) throw std::runtime_error("rans: table larger than alphabet");
        u64 total = 0, prev = 0;
        for (u64 i = 0; i < m; ++i) {
            const u64 ds = vbyte_read(p, len, pos, strict);
            const u64 fm1 = vbyte_read(p, len, pos, strict);   // freq - 1: must be < scale (guards the +1 wrap)
            const u64 s = (i == 0 ? 0 : prev) + ds;
            if (s >= alphabet || (i > 0 && s <= prev) || fm1 >= kRansProbScale)
                throw std::runtime_error("rans: malformed table");
            const u64 f = fm1 + 1;
            syms.push_back(static_cast<u16>(s));
            freq.push_back(static_cast<u32>(f));
            total += f;
            prev = s;
        }
        if (m > 0 && total != kRansProbScale) throw std::runtime_error("rans: table does not sum to the scale");
        finish_cum_();
        return pos;
    }

    // Symbol for a slot in [0, scale): binary search over cum.
    int slot_index(u32 slot) const {
        // bucket lookup, then a short forward scan (same result as the binary
        // search over cum it replaces: the last k with cum[k] <= slot)
        int k = lut[slot >> 8];
        while (cum[static_cast<std::size_t>(k) + 1] <= slot) ++k;
        return k;
    }

private:
    void finish_cum_() {
        cum.assign(freq.size() + 1, 0);
        for (std::size_t i = 0; i < freq.size(); ++i) cum[i + 1] = cum[i] + freq[i];
        lut.assign(freq.empty() ? 0 : 256, 0);
        std::size_t k = 0;
        for (u32 b = 0; b < lut.size(); ++b) {
            const u32 slot = b << 8;
            while (k + 1 < freq.size() && cum[k + 1] <= slot) ++k;
            lut[b] = static_cast<u16>(k);
        }
    }

    static void vbyte_append(std::vector<std::uint8_t>& out, u64 v) {
        while (v >= 0x80) { out.push_back(static_cast<std::uint8_t>(v | 0x80)); v >>= 7; }
        out.push_back(static_cast<std::uint8_t>(v));
    }
    static u64 vbyte_read(const std::uint8_t* p, std::size_t len, std::size_t& pos, bool strict = false) {
        u64 v = 0;
        for (int shift = 0;; shift += 7) {
            if (pos >= len || shift > 63) throw std::runtime_error("rans: truncated table");
            const std::uint8_t b = p[pos++];
            if (shift == 63 && b > 1) throw std::runtime_error("rans: table VByte above 2^64");
            v |= static_cast<u64>(b & 0x7f) << shift;
            if (!(b & 0x80)) {
                if (strict && shift > 0 && b == 0) throw std::runtime_error("rans: non-canonical VByte in a strict table");
                return v;
            }
        }
    }
};

// Encode one block of (context index, symbol) pairs. Output = bytes to be
// read FORWARD by rans_decode_block (the encoder reverses internally).
// Id = u16 (PWT1, PWT2 context set 1) or u32 (PWT2 context set 2, which can
// exceed 65,535 contexts); the bytes do not depend on the id type.
template <typename Id>
inline void rans_encode_block_ids(const std::vector<std::pair<Id, u16>>& pairs,
                              const std::vector<RansContext>& ctx,
                              std::vector<std::uint8_t>& out) {
    std::vector<std::uint8_t> rev;
    rev.reserve(pairs.size() / 2 + 8);
    u32 x = kRansL;
    for (std::size_t i = pairs.size(); i-- > 0;) {
        const RansContext& c = ctx[pairs[i].first];
        const int k = c.index_of(pairs[i].second);
        if (k < 0) throw std::runtime_error("rans: symbol absent from its context table");
        const u32 f = c.freq[k];
        const u32 xmax = ((kRansL >> kRansProbBits) << 8) * f;
        while (x >= xmax) { rev.push_back(static_cast<std::uint8_t>(x)); x >>= 8; }
        x = ((x / f) << kRansProbBits) + (x % f) + c.cum[k];
    }
    for (int i = 0; i < 4; ++i) { rev.push_back(static_cast<std::uint8_t>(x)); x >>= 8; }
    out.insert(out.end(), rev.rbegin(), rev.rend());
}
inline void rans_encode_block(const std::vector<std::pair<u16, u16>>& pairs,
                              const std::vector<RansContext>& ctx,
                              std::vector<std::uint8_t>& out) {
    rans_encode_block_ids(pairs, ctx, out);
}

// Forward decoder over one block's bytes. The caller supplies the context
// of each symbol as it goes (the PWT tree walk determines it).
class RansBlockDecoder {
public:
    RansBlockDecoder() = default;
    RansBlockDecoder(const std::uint8_t* p, std::size_t len) : p_(p), len_(len) {
        if (len < 4) throw std::runtime_error("rans: block too short");
        x_ = 0;
        for (int i = 0; i < 4; ++i) x_ = (x_ << 8) | p_[pos_++];
        if (x_ < kRansL || x_ >= (kRansL << 8)) throw std::runtime_error("rans: bad initial state");
    }
    u16 decode(const RansContext& c) {
        if (c.empty()) throw std::runtime_error("rans: symbol requested from an empty context");
        const u32 slot = x_ & (kRansProbScale - 1);
        const int k = c.slot_index(slot);
        x_ = c.freq[k] * (x_ >> kRansProbBits) + slot - c.cum[k];
        while (x_ < kRansL) {
            if (pos_ >= len_) throw std::runtime_error("rans: block overrun");
            x_ = (x_ << 8) | p_[pos_++];
        }
        return c.syms[k];
    }
    // True iff every byte was consumed and the state returned to L.
    // Generic decode step for callers with their own model store: find(slot,
    // freq, cum) returns the symbol whose range [cum, cum + freq) holds slot.
    // Same state transition as decode() (used by PWT2 context sets 2 and 3).
    template <typename Find>
    u16 decode_by(Find&& find) {
        const u32 slot = x_ & (kRansProbScale - 1);
        u32 f = 0, c = 0;
        const u16 s = find(slot, f, c);
        x_ = f * (x_ >> kRansProbBits) + slot - c;
        while (x_ < kRansL) {
            if (pos_ >= len_) throw std::runtime_error("rans: block overrun");
            x_ = (x_ << 8) | p_[pos_++];
        }
        return s;
    }
    bool finished() const { return pos_ == len_ && x_ == kRansL; }

private:
    const std::uint8_t* p_ = nullptr;
    std::size_t len_ = 0, pos_ = 0;
    u32 x_ = 0;
};

} // namespace cn
