#pragma once
// prime_bitmap.hpp -- the PWT decoder's prime table: a bitmap of the odd
// primes below B (bit (p-1)/2 set iff p is an odd prime) with a rank
// directory, rebuilt at load by a sieve of Eratosthenes over the odd
// numbers (B = 1e8: 6.25 MB of bits + 0.4 MB of directory, ~0.3 s).
//
//   index(p) = number of odd primes < p          (rank; 3 -> 0, 5 -> 1, ...)
//   prime(i) = the odd prime with index(p) == i   (select)
//
// The PWT1 "prime-index universe" codes a child prime q as a gap of prime
// INDICES, so the encoder needs rank and the decoder needs select.

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "u128.hpp"

namespace cn {

class PrimeBitmap {
public:
    explicit PrimeBitmap(u64 limit) : limit_(limit) {
        if (limit < 3) throw std::runtime_error("prime bitmap: limit too small");
        nbits_ = limit / 2;                  // odd numbers 1, 3, 5, ... < limit (bit i = 2i+1)
        words_.assign((nbits_ + 63) / 64, ~0ull);
        clear_(0);                           // 1 is not prime
        for (u64 i = 1; ; ++i) {             // i -> odd number 2i+1
            const u64 p = 2 * i + 1;
            if (p * p >= limit) break;
            if (!test_(i)) continue;
            for (u64 m = p * p; m < limit; m += 2 * p) clear_((m - 1) / 2);
        }
        if (nbits_ % 64) words_.back() &= (1ull << (nbits_ % 64)) - 1;
        // rank directory: primes strictly before each 512-bit block
        const std::size_t blocks = (words_.size() + 7) / 8;
        block_rank_.assign(blocks + 1, 0);
        for (std::size_t b = 0; b < blocks; ++b) {
            u64 c = 0;
            for (std::size_t w = b * 8; w < std::min(words_.size(), (b + 1) * 8); ++w) c += pop_(words_[w]);
            block_rank_[b + 1] = block_rank_[b] + c;
        }
        count_ = block_rank_[blocks];
    }

    u64 limit() const { return limit_; }
    u64 count() const { return count_; }     // number of odd primes below limit

    bool is_prime(u64 p) const { return p >= 3 && p < limit_ && (p & 1) && test_((p - 1) / 2); }

    // Number of odd primes < p (p odd, 3 <= p < limit).
    u64 rank(u64 p) const {
        if (p < 3 || p >= limit_ || !(p & 1)) throw std::runtime_error("prime bitmap: rank out of range");
        const u64 i = (p - 1) / 2;
        const std::size_t w = i / 64;
        u64 r = block_rank_[w / 8];
        for (std::size_t k = (w / 8) * 8; k < w; ++k) r += pop_(words_[k]);
        const u64 mask = (i % 64) ? ((1ull << (i % 64)) - 1) : 0;
        return r + pop_(words_[w] & mask);
    }

    // The odd prime with rank i (0 <= i < count()).
    u64 select(u64 i) const {
        if (i >= count_) throw std::runtime_error("prime bitmap: select out of range");
        // binary search the block whose cumulative rank covers i
        std::size_t lo = 0, hi = block_rank_.size() - 1;   // block_rank_[hi] > i
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (block_rank_[mid] <= i) lo = mid; else hi = mid;
        }
        u64 r = block_rank_[lo];
        std::size_t w = lo * 8;
        while (true) {
            const u64 c = pop_(words_[w]);
            if (r + c > i) break;
            r += c;
            ++w;
        }
        return 2 * (w * 64 + nth_set_(words_[w], i - r)) + 1;
    }

    // The g-th odd prime after the odd prime p (g >= 1), whose index the
    // caller knows to be idx: a forward popcount scan when the target is
    // near, select() otherwise. Returns the same value as select(idx).
    u64 advance(u64 p, u64 g, u64 idx) const {
        const u64 i = (p - 1) / 2 + 1;         // first candidate bit after p
        std::size_t w = i / 64;
        if (w >= words_.size()) return select(idx);
        u64 word = words_[w] & (~0ull << (i % 64));
        for (int steps = 0; steps < 64; ++steps) {
            const u64 c = pop_(word);
            if (c >= g) return 2 * (w * 64 + nth_set_(word, g - 1)) + 1;
            g -= c;
            if (++w >= words_.size()) break;
            word = words_[w];
        }
        return select(idx);
    }

    // The next odd prime after p (p < limit), or 0 if none below limit.
    u64 next_prime(u64 p) const {
        u64 i = (p + 1) / 2;   // first odd candidate above p
        if (p % 2 == 0) i = p / 2;   // p even: candidate p+1 -> index p/2
        while (i < nbits_) {
            const std::size_t w = i / 64;
            u64 word = words_[w] & (~0ull << (i % 64));
            if (word) return 2 * (w * 64 + ctz_(word)) + 1;
            i = (w + 1) * 64;
        }
        return 0;
    }

private:
    bool test_(u64 i) const { return (words_[i / 64] >> (i % 64)) & 1ull; }
    void clear_(u64 i) { words_[i / 64] &= ~(1ull << (i % 64)); }
    static unsigned pop_(u64 x) { return static_cast<unsigned>(__builtin_popcountll(x)); }
    static unsigned ctz_(u64 x) { return static_cast<unsigned>(__builtin_ctzll(x)); }
    // position of the k-th (0-based) set bit of x
    static unsigned nth_set_(u64 x, u64 k) {
        for (u64 j = 0; j < k; ++j) x &= x - 1;
        return ctz_(x);
    }

    u64 limit_, nbits_ = 0, count_ = 0;
    std::vector<u64> words_;
    std::vector<u64> block_rank_;
};

} // namespace cn
