#pragma once
// fastbits.hpp -- word-based MSB-first bit I/O, byte-compatible with
// bic.hpp's BitWriter / BitReader (which move one bit at a time). Used by
// the PWT1 and PWT2 raw bit streams (both decoders use the same reader;
// the container bytes are those of the bit-at-a-time writer).

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "u128.hpp"

namespace cn {

class FastBitWriter {
public:
    explicit FastBitWriter(std::vector<std::uint8_t>& out) : out_(out) {}
    void put(u64 v, int w) {   // low w bits of v, 0 <= w <= 64
        if (w <= 0) return;
        if (w > 56) { put(v >> 32, w - 32); put(v & 0xFFFFFFFFull, 32); return; }
        acc_ = (acc_ << w) | (v & ((1ull << w) - 1));
        n_ += w;
        while (n_ >= 8) { out_.push_back(static_cast<std::uint8_t>(acc_ >> (n_ - 8))); n_ -= 8; }
    }
    void put_wide(u128 v, int w) {   // any width up to 128
        while (w > 48) { put(static_cast<u64>(v >> (w - 48)), 48); w -= 48; }
        put(static_cast<u64>(v), w);
    }
    void align() {
        if (n_) { out_.push_back(static_cast<std::uint8_t>(acc_ << (8 - n_))); n_ = 0; }
        acc_ = 0;
    }
    u64 bits_written() const { return 8 * static_cast<u64>(out_.size()) + static_cast<u64>(n_); }

private:
    std::vector<std::uint8_t>& out_;
    u64 acc_ = 0;
    int n_ = 0;
};

class FastBitReader {
public:
    FastBitReader(const std::uint8_t* p, std::size_t len) : p_(p), len_(len), nbits_(8 * static_cast<u64>(len)) {}
    u64 get(int w) {   // 0 <= w <= 64
        if (w <= 0) return 0;
        if (w > 56) { const u64 hi = get(w - 32); return (hi << 32) | get(32); }
        if (pos_ + static_cast<u64>(w) > nbits_) throw std::runtime_error("raw bit stream overrun");
        const u64 byte = pos_ >> 3;
        u64 v;
        if (byte + 8 <= len_) {
            std::memcpy(&v, p_ + byte, 8);
            v = __builtin_bswap64(v);
        } else {
            v = 0;
            for (u64 i = 0; i < 8; ++i) v = (v << 8) | (byte + i < len_ ? p_[byte + i] : 0);
        }
        v <<= (pos_ & 7);
        v >>= (64 - w);
        pos_ += static_cast<u64>(w);
        return v;
    }
    u64 bit_pos() const { return pos_; }
    u64 total_bits() const { return nbits_; }

private:
    const std::uint8_t* p_;
    std::size_t len_;
    u64 nbits_;
    u64 pos_ = 0;
};

} // namespace cn
