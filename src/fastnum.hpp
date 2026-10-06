#pragma once
// fastnum.hpp -- division-free helpers for the PWT decoders' per-number
// reconstruction n = d * (d^-1 mod lambda(d)). Results equal
// lambda_bucket.hpp's gcd64 and the exact modular inverse for every even
// modulus below 2^64 (test_pwt2 checks against 128-bit arithmetic;
// inv_mod64 defers to 128 bits for moduli >= 2^63).
//   gcd64_bin      Stein's binary gcd with count-trailing-zeros
//   inv_mod_even   a^-1 mod m for ODD a and EVEN m (the lambda case: d is
//                  odd, lambda(d) is even): y = m^-1 mod a by a binary
//                  extended gcd (a is odd), then x = (1 - m y) / a mod m,
//                  the exact division done as a multiplication by a^-1
//                  mod 2^64 (Newton).

#include <cstdint>
#include <stdexcept>

#include "u128.hpp"

namespace cn {

inline u64 gcd64_bin(u64 a, u64 b) {
    if (a == 0) return b;
    if (b == 0) return a;
    const int s = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) { const u64 t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return a << s;
}

// a^-1 mod p for odd p > 1 and gcd(a, p) = 1 (binary extended gcd).
inline u64 inv_mod_odd(u64 a, u64 p) {
    a %= p;
    if (a == 0) throw std::runtime_error("inv_mod_odd: not invertible");
    u64 u = a, v = p, x1 = 1, x2 = 0;
    auto halve = [p](u64 x) { return (x & 1) ? (x >> 1) + (p >> 1) + 1 : x >> 1; };   // x/2 mod p, p odd, x < p
    while (u != 1 && v != 1) {
        while (!(u & 1)) { u >>= 1; x1 = halve(x1); }
        while (!(v & 1)) { v >>= 1; x2 = halve(x2); }
        if (u >= v) { u -= v; x1 = x1 >= x2 ? x1 - x2 : x1 + (p - x2); }
        else        { v -= u; x2 = x2 >= x1 ? x2 - x1 : x2 + (p - x1); }
        if (u == 0 || v == 0) throw std::runtime_error("inv_mod_odd: not invertible");
    }
    return u == 1 ? x1 : x2;
}

// inverse of an odd a modulo 2^64 (Newton: each step doubles the correct bits)
inline u64 inv_pow2_64(u64 a) {
    u64 x = a;                 // correct to 3 bits for odd a
    for (int i = 0; i < 5; ++i) x *= 2 - a * x;
    return x;
}

// a^-1 mod m for odd a, even m >= 2, gcd(a, m) = 1.
inline u64 inv_mod_even(u64 a, u64 m) {
    a %= m;
    if (a == 1) return 1;
    if (!(a & 1)) throw std::runtime_error("inv_mod_even: a must be odd");
    const u64 y = inv_mod_odd(m % a, a);                  // m y = 1 (mod a)
    const u64 t = (m * y - 1) * inv_pow2_64(a);           // (m y - 1) / a exactly; 0 < t < m
    return m - t;                                         // x = (1 - m y) / a = -t (mod m)
}

// d^-1 mod lambda for the reconstruction (odd a, even m, gcd 1): 32-bit Euclid
// when m < 2^32 (96% of the d_min table; 32-bit division is much cheaper),
// else the division-free even-modulus inverse above (correct for every even
// m < 2^64, unlike a signed 64-bit Euclid, which overflows at m >= 2^63).
inline u64 inv_mod_lambda(u64 a, u64 m) {
    if (m < (u64(1) << 32)) {
        std::uint32_t r0 = static_cast<std::uint32_t>(m), r1 = static_cast<std::uint32_t>(a % m);
        std::int64_t t0 = 0, t1 = 1;
        while (r1) {
            const std::uint32_t q = r0 / r1;
            const std::uint32_t r2 = r0 - q * r1;
            const std::int64_t t2 = t0 - static_cast<std::int64_t>(q) * t1;
            r0 = r1; r1 = r2; t0 = t1; t1 = t2;
        }
        if (r0 != 1) throw std::runtime_error("inv_mod_lambda: not invertible");
        if (t0 < 0) t0 += static_cast<std::int64_t>(m);
        return static_cast<u64>(t0);
    }
    return inv_mod_even(a, m);
}

} // namespace cn
