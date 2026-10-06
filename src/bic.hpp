#pragma once
// bic.hpp -- binary interpolative coding over sorted u128 sequences, exactly
// as specified in docs/CND_FORMAT.md: strict pre-order midpoint recursion
// (LEFT before RIGHT), centered minimal binary spans, MSB-first bitstream.

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "u128.hpp"

namespace cn {

class BitWriter {
public:
    explicit BitWriter(std::vector<std::uint8_t>& out) : out_(out) {}

    // Append the low `width` bits of v, MSB first. width <= 128.
    void put(u128 v, int width) {
        for (int i = width - 1; i >= 0; --i) {
            acc_ = (acc_ << 1) | static_cast<unsigned>((v >> i) & 1);
            if (++nbits_ == 8) { out_.push_back(static_cast<std::uint8_t>(acc_)); acc_ = 0; nbits_ = 0; }
        }
    }

    // Zero-pad to a byte boundary (blob terminator).
    void align() {
        if (nbits_) { out_.push_back(static_cast<std::uint8_t>(acc_ << (8 - nbits_))); acc_ = 0; nbits_ = 0; }
    }

private:
    std::vector<std::uint8_t>& out_;
    unsigned acc_ = 0;
    int nbits_ = 0;
};

class BitReader {
public:
    BitReader(const std::uint8_t* p, std::size_t len) : p_(p), len_(len) {}

    u128 get(int width) {
        u128 v = 0;
        for (int i = 0; i < width; ++i) {
            if (pos_ >= len_ * 8) throw std::runtime_error("bic: bitstream overrun");
            const unsigned bit = (p_[pos_ >> 3] >> (7 - (pos_ & 7))) & 1u;
            v = (v << 1) | bit;
            ++pos_;
        }
        return v;
    }

    std::size_t bit_pos() const { return pos_; }

private:
    const std::uint8_t* p_;
    std::size_t len_;
    std::size_t pos_ = 0;
};

namespace bic_detail {

inline int bits_of(u128 v) {
    int b = 0;
    while (v > 0) { v >>= 1; ++b; }
    return b;
}

// Centered minimal binary, spec section "BIC (normative)".
inline void cmb_put(BitWriter& bw, u128 v, u128 S) {
    if (S <= 1) return;
    const int L = bits_of(S - 1);                  // ceil(log2 S)
    const u128 a = (static_cast<u128>(1) << L) - S;
    const u128 b = S - a;
    u128 vp = v + S - b / 2;                       // (v - floor(b/2)) mod S
    if (vp >= S) vp -= S;
    if (vp >= S) vp -= S;                          // v < S, shift < S: two folds suffice
    if (vp < a) bw.put(vp, L - 1);
    else        bw.put(vp + a, L);
}

inline u128 cmb_get(BitReader& br, u128 S) {
    if (S <= 1) return 0;
    const int L = bits_of(S - 1);
    const u128 a = (static_cast<u128>(1) << L) - S;
    const u128 b = S - a;
    u128 vp = br.get(L - 1);
    if (vp >= a) vp = ((vp << 1) | br.get(1)) - a;
    u128 v = vp + b / 2;                           // undo the rotation
    if (v >= S) v -= S;
    return v;
}

} // namespace bic_detail

// Append the CND1 class blob for the sorted strictly-increasing ks:
// k_first/k_last at end_width bits, pre-order (LEFT first) interior,
// zero-padded to a byte boundary. ks may be empty.
inline void bic_append(std::vector<std::uint8_t>& out,
                       const std::vector<u128>& ks, int end_width) {
    BitWriter bw(out);
    const std::size_t c = ks.size();
    if (c >= 1) bw.put(ks.front(), end_width);
    if (c >= 2) bw.put(ks.back(), end_width);
    if (c >= 3) {
        struct Frame { std::size_t ilo, ihi; u128 vlo, vhi; };
        std::vector<Frame> stack{{1, c - 2, ks.front(), ks.back()}};
        while (!stack.empty()) {
            const Frame f = stack.back();
            stack.pop_back();
            const std::size_t mid = f.ilo + (f.ihi - f.ilo) / 2;
            const u128 lo = f.vlo + 1 + (mid - f.ilo);
            const u128 hi = f.vhi - 1 - (f.ihi - mid);
            bic_detail::cmb_put(bw, ks[mid] - lo, hi - lo + 1);
            // pre-order, LEFT before RIGHT: push RIGHT first so LEFT pops next
            if (mid < f.ihi) stack.push_back({mid + 1, f.ihi, ks[mid], f.vhi});
            if (mid > f.ilo) stack.push_back({f.ilo, mid - 1, f.vlo, ks[mid]});
        }
    }
    bw.align();
}

// Decode a class blob of `count` values. Throws on overrun, a violated
// frame invariant, an under-consumed blob, or nonzero padding (decoder
// obligations 3 and 4 of the spec: the blob must contain exactly the
// demanded bits plus 0-7 zero pad bits).
inline std::vector<u128> bic_decode(const std::uint8_t* blob, std::size_t len,
                                    u64 count, int end_width) {
    std::vector<u128> ks(count);
    if (count == 0) {
        if (len != 0) throw std::runtime_error("bic: nonempty blob for empty class");
        return ks;
    }
    BitReader br(blob, len);
    ks.front() = br.get(end_width);
    if (count >= 2) {
        ks.back() = br.get(end_width);
        if (ks.back() <= ks.front())
            throw std::runtime_error("bic: ends not increasing");
    }
    if (count >= 3) {
        struct Frame { std::size_t ilo, ihi; u128 vlo, vhi; };
        std::vector<Frame> stack{{1, count - 2, ks.front(), ks.back()}};
        while (!stack.empty()) {
            const Frame f = stack.back();
            stack.pop_back();
            const std::size_t mid = f.ilo + (f.ihi - f.ilo) / 2;
            const u128 lo = f.vlo + 1 + (mid - f.ilo);
            const u128 hi = f.vhi - 1 - (f.ihi - mid);
            if (hi < lo) throw std::runtime_error("bic: empty span (corrupt counts)");
            ks[mid] = lo + bic_detail::cmb_get(br, hi - lo + 1);
            if (mid < f.ihi) stack.push_back({mid + 1, f.ihi, ks[mid], f.vhi});
            if (mid > f.ilo) stack.push_back({f.ilo, mid - 1, f.vlo, ks[mid]});
        }
    }
    // exact consumption: only 0-7 zero pad bits may remain
    const std::size_t used = br.bit_pos();
    if ((used + 7) / 8 != len)
        throw std::runtime_error("bic: blob length does not match consumed bits");
    while (br.bit_pos() < len * 8)
        if (br.get(1) != 0) throw std::runtime_error("bic: nonzero padding");
    return ks;
}

} // namespace cn
