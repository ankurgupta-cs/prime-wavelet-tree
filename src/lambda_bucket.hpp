#pragma once
// lambda_bucket.hpp -- the lambda-d bucket map: n = d * (r* + k lambda(d))
// for a divisor d of a Carmichael number n, and the choice of d_min(n).
//
// For a Carmichael n and any divisor d | n: n == 0 (mod d) and, by Korselt,
// n == 1 (mod lambda(d)). Carmichael numbers are cyclic, so gcd(d, lambda(d))
// = 1 and CRT pins n to a single class mod d*lambda(d):
//     n = d * (r* + k*lambda(d)),   r* = d^{-1} mod lambda(d).
// The encoder buckets each n under the divisor d <= 10^8 that MAXIMIZES
// d*lambda(d) (tie-break: larger d), minimizing the bits of the index k.
// Every table entry has smallest prime < B^{1/3} = 10^8, so a bucket exists.
//
// Ranges: d <= 10^8 and lambda(d) < d, so d*lambda(d) < 10^16 fits u64 with
// slack. For REAL table entries k fits u64 (measured max 2^56.9 over the
// full 10^24 table; the encoder enforces the u64 ceiling). Two tail
// mechanisms: (a) lone-prime buckets -- only 3-factor n = pqr with
// q,r > 10^8 (omega >= 4 would force n > 3e24 > B), where Korselt rigidity
// (q-1 | (p-1)(p+s), s <= p) forces p >~ 7e3 and k <= ~2^54; 561 is the
// ONLY Carmichael 3qr and buckets at d=561. (b) composite buckets whose
// lambda is anomalously small (shared p_i - 1 structure), which is how the
// measured max exceeds the lone-prime bound. k stays u128 in the container
// purely for format robustness against non-table inputs.

#include <cstdint>

#include "sieve.hpp"
#include "u128.hpp"

namespace cn {

inline constexpr u64 kBucketCap = 100'000'000ULL;   // divisors d <= 10^8
inline constexpr u128 kCnb1Bound =
    static_cast<u128>(1'000'000'000'000ULL) * 1'000'000'000'000ULL;   // 10^24

inline u64 gcd64(u64 a, u64 b) {
    while (b != 0) { const u64 t = a % b; a = b; b = t; }
    return a;
}
inline u64 lcm64(u64 a, u64 b) { return a / gcd64(a, b) * b; }

struct BucketChoice {
    u64 d = 0;        // chosen divisor (bucket id)
    u64 lambda = 0;   // lambda(d)
    u64 dl = 0;       // d * lambda(d), the maximized quantity
};

// Enumerates subset products of the (ascending, distinct) prime list with
// product <= cap and returns the maximizer of d*lambda(d); ties go to the
// larger d so the encoding is canonical. k <= 14 bounds the recursion.
template <typename Self>
inline void bucket_dfs_(Self&& self, const u64* p, int k, int i, u64 d,
                        u64 lam, u64 cap, BucketChoice& best) {
    if (d > 1) {
        const u64 dl = d * lam;
        if (dl > best.dl || (dl == best.dl && d > best.d)) best = {d, lam, dl};
    }
    for (int j = i; j < k; ++j) {
        if (p[j] > cap / d) break;             // ascending primes: no later j fits
        self(self, p, k, j + 1, d * p[j], lcm64(lam, p[j] - 1), cap, best);
    }
}

inline BucketChoice choose_bucket(const u64* primes, int k,
                                  u64 cap = kBucketCap) {
    BucketChoice best;
    auto rec = [](auto&& self, const u64* p, int kk, int i, u64 d, u64 lam,
                  u64 cap_, BucketChoice& b) -> void {
        bucket_dfs_(self, p, kk, i, d, lam, cap_, b);
    };
    rec(rec, primes, k, 0, 1, 1, cap, best);
    return best;
}

// Windowed variant (experiment): prefer the max-d*lambda divisor
// with d in [lo, hi]; if the window holds no divisor, fall back to the best
// divisor below lo (which for n < lo is n itself, k = 0). Returns the choice
// plus which tier served it (0 = window, 1 = fallback).
struct WindowChoice {
    BucketChoice c;
    int tier = -1;
};

template <typename Self>
inline void window_dfs_(Self&& self, const u64* p, int k, int i, u64 d,
                        u64 lam, u64 lo, u64 hi, BucketChoice& win,
                        BucketChoice& below) {
    if (d > 1) {
        const u64 dl = d * lam;
        if (d >= lo) {
            if (dl > win.dl || (dl == win.dl && d > win.d)) win = {d, lam, dl};
        } else {
            if (dl > below.dl || (dl == below.dl && d > below.d))
                below = {d, lam, dl};
        }
    }
    for (int j = i; j < k; ++j) {
        if (p[j] > hi / d) break;              // ascending primes
        self(self, p, k, j + 1, d * p[j], lcm64(lam, p[j] - 1), lo, hi, win,
             below);
    }
}

inline WindowChoice choose_bucket_window(const u64* primes, int k, u64 lo,
                                         u64 hi) {
    BucketChoice win, below;
    auto rec = [](auto&& self, const u64* p, int kk, int i, u64 d, u64 lam,
                  u64 lo_, u64 hi_, BucketChoice& w, BucketChoice& b) -> void {
        window_dfs_(self, p, kk, i, d, lam, lo_, hi_, w, b);
    };
    rec(rec, primes, k, 0, 1, 1, lo, hi, win, below);
    if (win.d != 0) return {win, 0};
    return {below, 1};
}

// ---- the d-only representation -------------------------------------------
// Among divisors d | n with lambda(d) > n/d (integer division; equivalent to
// d*lambda(d) > n without overflow), the cofactor n/d is forced to be the
// least positive residue d^{-1} mod lambda(d), so n = d * (d^{-1} mod
// lambda(d)) is reconstructible from d ALONE, and n -> d_min is injective.
// d = n is always valid (n == 1 mod lambda(n) by Korselt), so the map is
// total. All arithmetic u128: d can approach n < 2^80, and lambda(d) =
// lcm(p-1) over d's primes stays below n.

inline u128 gcd128(u128 a, u128 b) {
    while (b != 0) { const u128 t = a % b; a = b; b = t; }
    return a;
}
inline u128 lcm128(u128 a, u128 b) { return a / gcd128(a, b) * b; }

struct DminChoice {
    u128 d = 0;
    u128 lambda = 0;
};

// DFS over subset products, ascending primes; a branch is pruned as soon as
// its product reaches the current best (extensions only grow).
template <typename Self>
inline void dmin_dfs_(Self&& self, const u64* p, int k, int i, u128 d,
                      u128 lam, u128 n, DminChoice& best) {
    if (d > 1 && d < best.d && lam > n / d) best = {d, lam};
    for (int j = i; j < k; ++j) {
        const u128 nd = d * p[j];
        if (nd >= best.d) break;               // ascending: no later j helps
        self(self, p, k, j + 1, nd, lcm128(lam, p[j] - 1), n, best);
    }
}

// Smallest divisor of n (given its prime list) with lambda(d) > n/d.
inline DminChoice choose_dmin(const u64* primes, int k, u128 n) {
    DminChoice best;
    best.d = n;                                 // d = n is always valid
    u128 lam = 1;
    for (int i = 0; i < k; ++i) lam = lcm128(lam, primes[i] - 1);
    best.lambda = lam;
    auto rec = [](auto&& self, const u64* p, int kk, int i, u128 d, u128 l,
                  u128 nn, DminChoice& b) -> void {
        dmin_dfs_(self, p, kk, i, d, l, nn, b);
    };
    rec(rec, primes, k, 0, 1, 1, n, best);
    return best;
}

// u128 modular inverse via extended Euclid on signed 128-bit accumulators
// (gcd(a, m) == 1 required; m >= 2).
inline u128 inv_mod128(u128 a, u128 m) {
    using i128 = __int128;
    i128 t = 0, newt = 1;
    i128 r = static_cast<i128>(m), newr = static_cast<i128>(a % m);
    while (newr != 0) {
        const i128 q = r / newr;
        i128 tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += static_cast<i128>(m);
    return static_cast<u128>(t);
}

// Reconstruction from d alone: n = d * (d^{-1} mod lambda(d)).
inline u128 dmin_value(u128 d, u128 lambda) {
    return d * inv_mod128(d, lambda);
}

// Modular inverse of a mod m (gcd(a, m) == 1 required; m >= 1).
inline u64 inv_mod64(u64 a, u64 m) {
    if (m == 1) return 0;
    if (m >> 63) return static_cast<u64>(inv_mod128(a % m, m));   // the signed 64-bit Euclid below overflows at m >= 2^63
    // extended Euclid on (a mod m, m) with signed accumulators
    std::int64_t t = 0, newt = 1;
    std::int64_t r = static_cast<std::int64_t>(m),
                 newr = static_cast<std::int64_t>(a % m);
    while (newr != 0) {
        const std::int64_t q = r / newr;
        std::int64_t tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += static_cast<std::int64_t>(m);
    return static_cast<u64>(t);
}

struct BucketParams {
    u64 lambda;   // lambda(d)
    u64 rstar;    // least positive d^{-1} mod lambda(d), in [0, lambda)
};

// Decode-side recomputation from the bucket id alone: factor d by SPF chase
// (d <= sieve limit), lambda = lcm(p-1), r* = d^{-1} mod lambda. The decoder
// therefore needs nothing but the N = 10^8 sieve.
inline BucketParams bucket_params(u64 d, const LinearSieve& s) {
    const Factorization f = s.factor(d);
    u64 lam = 1;
    for (const PrimePower& pp : f.factors) lam = lcm64(lam, pp.p - 1);
    return {lam, inv_mod64(d, lam)};   // gcd(d, lambda(d)) = 1 for CN divisors
}

// Index of n within bucket d: k = (n/d - r*) / lambda. Returns false if n is
// not in the progression (which would falsify Korselt/cyclicity for a true
// table entry -- the encoder treats that as data corruption).
inline bool bucket_index(u128 n, u64 d, const BucketParams& bp, u128& k_out) {
    k_out = 0;
    if (n % d != 0) return false;
    const u128 m = n / d;
    const u128 r = bp.rstar;
    if (m < r || (m - r) % bp.lambda != 0) return false;
    k_out = (m - r) / bp.lambda;
    return true;
}

// Reconstruction: n = d * (r* + k * lambda).
inline u128 bucket_value(u64 d, const BucketParams& bp, u128 k) {
    return static_cast<u128>(d) * (bp.rstar + k * static_cast<u128>(bp.lambda));
}

} // namespace cn
