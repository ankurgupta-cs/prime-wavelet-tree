#pragma once
// mr64.hpp -- deterministic Miller-Rabin primality for 64-bit integers.
//
// The first twelve primes {2,3,5,7,...,37} as bases are provably sufficient
// for every n < psi_12 = 318,665,857,834,031,151,167,461 ~ 3.19e23
// (Sorenson-Webster, "Strong pseudoprimes to twelve prime bases", Math.
// Comp. 86 (2017)), which covers all of u64 (~1.8e19) with room to spare.
// mulmod goes through unsigned __int128; no GMP.

#include <bit>
#include <cstdint>

#include "u128.hpp"

namespace cn {

inline u64 mulmod64(u64 a, u64 b, u64 m) {
    return static_cast<u64>((static_cast<u128>(a) * b) % m);
}

inline u64 powmod64(u64 b, u64 e, u64 m) {
    u64 r = 1 % m;
    b %= m;
    while (e) {
        if (e & 1) r = mulmod64(r, b, m);
        b = mulmod64(b, b, m);
        e >>= 1;
    }
    return r;
}

// Deterministic for all n < 2^64.
inline bool is_prime_u64(u64 n) {
    constexpr u64 bases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (u64 p : bases)
        if (n % p == 0) return n == p;
    if (n < 41 * 41) return true;   // no prime factor <= 37 and below 41^2

    const int s = std::countr_zero(n - 1);
    const u64 d = (n - 1) >> s;
    for (u64 a : bases) {
        u64 x = powmod64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int i = 1; i < s; ++i) {
            x = mulmod64(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

} // namespace cn
