#pragma once
// segment.hpp -- windowed segmented sieve of Eratosthenes over [a, b), for
// bounds up to a configured maximum (default 10^12, which covers every prime
// factor appearing in the 10^24 Carmichael table).
//
// Role in the project: this is our replacement for the streaming role that
// Sorenson's Rollsieve plays in Shallue-Webster's tabulation code. Unlike the
// Rollsieve it offers random-access windows (each call is independent, so
// callers can parallelize by window), reports exponents, and certifies its
// cofactor contract: for_factored_in yields cofactor == 1 or a single prime.
//
// Memory: the base-prime table is pi(sqrt(max_bound)) u32s -- 0.3 MB at the
// 10^12 default. Transient window buffers only during calls. Constructing
// with max_bound near 2^64 stores pi(2^32) = 203M primes (~812 MB): legal,
// but that is a deliberate caller choice, not the default.

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "sieve.hpp"

namespace cn {

class SegmentedSieve {
public:
    explicit SegmentedSieve(u64 max_bound = 1'000'000'000'000ull)
        : max_bound_(max_bound) {
        if (max_bound < 4)
            throw std::invalid_argument("SegmentedSieve: max_bound must be >= 4");
        base_limit_ = ceil_sqrt_(max_bound);
        // Simple odd-only Eratosthenes for the base primes <= base_limit_.
        base_primes_.push_back(2);
        const u64 half = base_limit_ / 2;          // odd o <-> integer 2o+1
        std::vector<std::uint8_t> comp(half + 1, 0);
        for (u64 o = 1; o <= half; ++o) {
            const u64 p = 2 * o + 1;
            if (comp[o]) continue;
            base_primes_.push_back(static_cast<u32>(p));
            if (p * p <= base_limit_)
                for (u64 m = (p * p) / 2; m <= half; m += p) comp[m] = 1;
        }
    }

    u64 max_bound() const { return max_bound_; }
    u64 base_limit() const { return base_limit_; }
    const std::vector<u32>& base_primes() const { return base_primes_; }

    // Calls f(p) for every prime p in [a, b), ascending. b <= max_bound().
    template <typename F>
    void for_primes_in(u64 a, u64 b, F&& f) const {
        check_range_(a, b);
        if (a <= 2 && 2 < b) f(static_cast<u64>(2));
        u64 lo = (a < 3) ? 3 : (a | 1);            // first odd >= max(a,3)
        std::vector<std::uint8_t> comp;
        while (lo < b) {
            const u64 hi = (b - lo > 2 * kPrimeSpan) ? lo + 2 * kPrimeSpan : b;
            const u64 count = (hi - lo + 1) / 2;   // odds in [lo, hi)
            comp.assign(count, 0);
            for (std::size_t i = 1; i < base_primes_.size(); ++i) {
                const u64 p = base_primes_[i];
                if (p * p >= hi) break;
                u64 m = ((lo + p - 1) / p) * p;    // first multiple >= lo
                if (m < p * p) m = p * p;
                if ((m & 1) == 0) m += p;          // odd multiples only
                for (; m < hi; m += 2 * p) comp[(m - lo) / 2] = 1;
            }
            for (u64 i = 0; i < count; ++i)
                if (!comp[i]) {
                    const u64 n = lo + 2 * i;
                    if (n != 1) f(n);
                }
            lo = hi | 1;                           // first odd >= hi
        }
    }

    // Calls f(n, factors, cofactor) for every n in [a, b), ascending, where
    // factors lists the PrimePower decomposition over primes <= base_limit()
    // (ascending, exact exponents; the vector is reused between calls) and
    // cofactor is the remaining part: always 1 or a single prime > sqrt(b-1).
    template <typename F>
    void for_factored_in(u64 a, u64 b, F&& f) const {
        check_range_(a, b);
        if (a < 1) a = 1;
        std::vector<u64> rem;
        std::vector<u32> plist;                    // kMaxDistinct slots per n
        std::vector<std::uint8_t> pcount;
        std::vector<PrimePower> factors;
        for (u64 lo = a; lo < b; ) {
            const u64 hi = (b - lo > kFactorSpan) ? lo + kFactorSpan : b;
            const u64 count = hi - lo;
            rem.resize(count);
            for (u64 i = 0; i < count; ++i) rem[i] = lo + i;
            plist.assign(count * kMaxDistinct, 0);
            pcount.assign(count, 0);
            // Every base prime (not just those <= sqrt) can divide entries in
            // the window, so all are scanned; primes above the window span
            // contribute at most one multiple, found by the same first-multiple
            // computation.
            for (const u32 p : base_primes_) {
                for (u64 m = ((lo + p - 1) / p) * p; m < hi; m += p)
                    record_(m - lo, p, plist, pcount);
            }
            for (u64 i = 0; i < count; ++i) {
                const u64 n = lo + i;
                factors.clear();
                u64 r = n;
                for (std::uint8_t j = 0; j < pcount[i]; ++j) {
                    const u32 p = plist[i * kMaxDistinct + j];
                    u32 e = 0;
                    while (r % p == 0) { r /= p; ++e; }
                    factors.push_back({p, e});
                }
                f(n, factors, r);                  // r == 1 or prime
            }
            lo = hi;
        }
    }

private:
    static constexpr u64 kPrimeSpan   = 1u << 19;  // odds per window (span 2^20)
    static constexpr u64 kFactorSpan  = 1u << 16;
    static constexpr int kMaxDistinct = 15;        // max distinct primes of a u64

    static u64 ceil_sqrt_(u64 x) {
        u64 r = static_cast<u64>(std::sqrt(static_cast<double>(x)));
        while (r > 0 && (r - 1) * (r - 1) >= x) --r;
        while (r * r < x) ++r;
        return r;
    }

    void check_range_(u64 a, u64 b) const {
        if (a > b) throw std::invalid_argument("SegmentedSieve: a > b");
        if (b > max_bound_)
            throw std::out_of_range("SegmentedSieve: b exceeds max_bound");
    }

    static void record_(u64 slot, u32 p, std::vector<u32>& plist,
                        std::vector<std::uint8_t>& pcount) {
        plist[slot * kMaxDistinct + pcount[slot]] = p;
        ++pcount[slot];
    }

    u64 max_bound_;
    u64 base_limit_;
    std::vector<u32> base_primes_;
};

} // namespace cn
