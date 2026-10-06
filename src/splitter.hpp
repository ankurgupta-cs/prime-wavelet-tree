#pragma once
// splitter.hpp -- u128 factoring beyond the sieve+MR64 reach: Brent rho plus
// the strong-Fermat-gcd split (Shallue-Webster arXiv:2506.09903 section 3.4),
// completing LinearSieve::factor's honest cofactor. Header-only, C++20, no
// GMP; all arithmetic in unsigned __int128.
//
// Certification domain: 12-base Miller-Rabin is a primality PROOF for every
// n < PSI_12 = 318,665,857,834,031,151,167,461 ~ 3.19e23 (Sorenson-Webster
// 2017). Everything this project must split sits far below that: table
// factors < 10^12, d_min < 2^70, other moduli < 2^80. A value >= PSI_12 that
// merely passes MR-12 is NOT recorded as a prime factor -- it stays in the
// uncertified/cofactor slot, mirroring sieve.hpp's contract that
// Factorization::complete() is a certificate.
//
// Modular multiplication: for moduli <= 2^64 the __int128 product path of
// mr64.hpp applies; above that, b is consumed in 32-bit limbs (Horner), which
// keeps every intermediate below 2^128 provided m < 2^95. The project ceiling
// of 10^24 < 2^80 clears that with room; powmod128 throws on m >= 2^95
// rather than wrapping silently.
//
// Split strategy per composite: peel tiny primes, then try the section-3.4
// gcd trick (one square-and-multiply ladder per base; on Fermat-pseudoprime
// structure -- Carmichael numbers and many of their divisors -- the gcds with
// X[i]+1 usually shatter n into several factors at once), then fall back to
// Brent's rho with batched gcds. Recursion certifies every leaf with MR-12.

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "mr64.hpp"
#include "sieve.hpp"
#include "u128.hpp"

namespace cn {

inline constexpr u128 PSI_12 =
    static_cast<u128>(318665857834ull) * 1000000000000ull + 31151167461ull;

namespace splitter_detail {

inline u128 gcd(u128 a, u128 b) {
    while (b) { const u128 t = a % b; a = b; b = t; }
    return a;
}

// a*b mod m for m < 2^95 (throws above); delegates to the 64-bit path when
// the modulus allows. See header comment for the overflow argument.
inline u128 mulmod(u128 a, u128 b, u128 m) {
    if (m <= UINT64_MAX)
        return mulmod64(static_cast<u64>(a % m), static_cast<u64>(b % m),
                        static_cast<u64>(m));
    a %= m;
    u128 r = 0;
    for (int sh = 64; sh >= 0; sh -= 32) {
        r = (r << 32) % m;
        const u64 limb = static_cast<u64>(b >> sh) & 0xFFFFFFFFull;
        r = (r + a * limb) % m;
    }
    return r;
}

inline u128 powmod(u128 b, u128 e, u128 m) {
    if (m >= static_cast<u128>(1) << 95)
        throw std::invalid_argument("splitter: modulus >= 2^95");
    u128 r = 1 % m;
    b %= m;
    while (e) {
        if (e & 1) r = mulmod(r, b, m);
        b = mulmod(b, b, m);
        e >>= 1;
    }
    return r;
}

inline int ctz128(u128 x) {
    const u64 lo = static_cast<u64>(x);
    if (lo) return std::countr_zero(lo);
    return 64 + std::countr_zero(static_cast<u64>(x >> 64));
}

inline constexpr u64 kBases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

} // namespace splitter_detail

// Strong probable-prime test to the twelve prime bases. A PROOF of primality
// for n < PSI_12; above that a true verdict means "probable prime" only
// (a false verdict is always a proof of compositeness).
inline bool is_prime_u128(u128 n) {
    namespace d = splitter_detail;
    if (n <= UINT64_MAX) return is_prime_u64(static_cast<u64>(n));
    for (u64 p : d::kBases)
        if (n % p == 0) return false;          // n > 37^2 here
    const int s = d::ctz128(n - 1);
    const u128 nd = (n - 1) >> s;
    for (u64 a : d::kBases) {
        u128 x = d::powmod(a, nd, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int i = 1; i < s; ++i) {
            x = d::mulmod(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

// Section-3.4 gcd trick. For odd n > 3 and base b, runs the strong-Fermat
// ladder and returns a nontrivial factor surfaced by
//   b^(n-1) - 1 = (b^n' - 1) * prod_{i<s} (b^(2^i n') + 1),
// i.e. gcd(n, b^n' - 1) and gcd(n, X[i] + 1); returns 0 if none surfaced
// (which says nothing about primality).
inline u128 fermat_gcd_split(u128 n, u64 b) {
    namespace d = splitter_detail;
    const u128 g0 = d::gcd(n, b);
    if (g0 != 1) return g0 < n ? g0 : 0;
    const int s = d::ctz128(n - 1);
    u128 x = d::powmod(b, (n - 1) >> s, n);
    if (x == 1 || x == n - 1) return 0;        // strong pseudoprime to b
    const u128 g1 = d::gcd(n, x - 1);
    if (g1 != 1 && g1 != n) return g1;
    for (int i = 0; i < s; ++i) {
        const u128 g = d::gcd(n, x + 1);
        if (g != 1 && g != n) return g;
        x = d::mulmod(x, x, n);
        if (x == 1) break;                     // earlier gcds had the split if any
    }
    return 0;
}

// Brent's cycle-finding rho with batched gcds. n must be odd, composite and
// coprime to small primes handled by the caller; returns a nontrivial factor.
// Deterministic: parameters sweep a fixed schedule, and for odd composite n
// with a prime factor p the map x^2+c mod n collides mod p within O(p^(1/4))
// steps for all but finitely many c, so the sweep terminates.
inline u128 brent_rho(u128 n) {
    namespace d = splitter_detail;
    for (u64 c = 1;; ++c) {
        if (c % 4 == 2) continue;              // x^2-2-like maps degenerate
        u128 y = 2 + c, q = 1, g = 1, ys = y, x = y;
        const u64 m = 128;
        for (u64 r = 1; g == 1; r <<= 1) {
            x = y;
            for (u64 i = 0; i < r; ++i) y = (d::mulmod(y, y, n) + c) % n;
            for (u64 k = 0; k < r && g == 1; k += m) {
                ys = y;
                const u64 lim = std::min(m, r - k);
                for (u64 i = 0; i < lim; ++i) {
                    y = (d::mulmod(y, y, n) + c) % n;
                    q = d::mulmod(q, x > y ? x - y : y - x, n);
                }
                g = d::gcd(q, n);
            }
        }
        if (g == n) {                          // batch overshot: replay singly
            g = 1;
            while (g == 1) {
                ys = (d::mulmod(ys, ys, n) + c) % n;
                g = d::gcd(x > ys ? x - ys : ys - x, n);
            }
        }
        if (g != n) return g;                  // else cycle degenerate: next c
    }
}

// Certified prime multiset of x, plus the product of any leaves >= PSI_12
// that MR-12 could not certify (always 1 in this project's domain).
struct SplitResult {
    std::vector<u128> primes;   // ascending after factor(); repeated = power
    u128 uncertified = 1;
};

namespace splitter_detail {

inline void split_rec(u128 x, SplitResult& out) {
    if (x == 1) return;
    if (is_prime_u128(x)) {
        if (x < PSI_12) out.primes.push_back(x);
        else out.uncertified *= x;             // probable prime, no proof
        return;
    }
    // The gcd trick fires on Fermat-pseudoprime structure and does so on the
    // first base or two when it fires at all; a wide base sweep only burns
    // ladders on ordinary composites before rho (measured on the d_min set:
    // 12-base sweep 90.5 us/d; the 2-base sweep is the hot path).
    u128 f = fermat_gcd_split(x, 2);
    if (!f) f = fermat_gcd_split(x, 3);
    if (!f) f = brent_rho(x);
    split_rec(f, out);
    split_rec(x / f, out);
}

} // namespace splitter_detail

// Full factorization of any u128 below 2^95. Every returned prime is MR-12
// certified (< PSI_12); anything else lands in uncertified.
inline SplitResult split_u128(u128 x) {
    namespace d = splitter_detail;
    SplitResult out;
    if (x <= 1) return out;
    const int tz = x == 0 ? 0 : d::ctz128(x);
    if (x != 0 && tz) { out.primes.assign(tz, 2); x >>= tz; }
    for (u64 p : d::kBases) {                  // peel tiny odd primes for rho
        if (p == 2) continue;
        while (x % p == 0) { out.primes.push_back(p); x /= p; }
    }
    d::split_rec(x, out);
    std::sort(out.primes.begin(), out.primes.end());
    return out;
}

// Wire-in for LinearSieve::factor's honest cofactor: split it and merge the
// leaves. Every cofactor prime exceeds the primes already recorded (sieve
// contract), so the merged list stays strictly ascending by construction.
// Leaves that do not fit PrimePower's u64 (impossible for table-derived
// values, whose primes are < 10^12) and uncertified leaves remain in
// f.cofactor, keeping complete() honest.
inline void finish_factorization(Factorization& f) {
    if (f.cofactor == 1) return;
    const SplitResult s = split_u128(f.cofactor);
    u128 remain = s.uncertified;
    std::size_t i = 0;
    while (i < s.primes.size()) {
        const u128 p = s.primes[i];
        std::size_t j = i;
        while (j < s.primes.size() && s.primes[j] == p) ++j;
        if (p <= UINT64_MAX)
            f.factors.push_back({static_cast<u64>(p), static_cast<u32>(j - i)});
        else
            for (std::size_t k = i; k < j; ++k) remain *= p;
        i = j;
    }
    f.cofactor = remain;
}

// Convenience: sieve peel + splitter finish in one call.
inline Factorization factor_full(const LinearSieve& sieve, u128 x,
                                 u64 tdiv_limit = UINT64_MAX) {
    Factorization f = sieve.factor(x, tdiv_limit);
    finish_factorization(f);
    return f;
}

} // namespace cn
