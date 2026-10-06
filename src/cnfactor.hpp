#pragma once
// cnfactor.hpp -- polynomial-time factorization of a Carmichael number (or
// of any divisor of one): Algorithm 1 of A. Shallue and J. Webster,
// "Algorithms for Carmichael numbers", arXiv:2506.09903.
//
//   For a Carmichael number n and ANY prime q | n, q - 1 divides n - 1, so
//   for every base b coprime to q the strong-Fermat ladder with the
//   exponent n - 1 = 2^s n',
//       X[i] = b^(2^i n') (mod m),   0 <= i < s,
//   computed modulo any divisor m of n, has each prime of m dividing exactly
//   one of the algebraic factors b^n' - 1, X[0] + 1, ..., X[s-1] + 1. The
//   gcds gcd(m, X[i] + 1) therefore split m unless every prime of m lands in
//   the same factor, which for a fixed composite m happens for at most half
//   of the bases. Consecutive prime bases 2, 3, 5, ... split every composite
//   after a few ladders; leaves are proven prime.
//
// The exponent belongs to n, not to the composite being split: that is the
// whole difference from a generic factoring routine (splitter.hpp keeps
// the ladder only for the top-level number and falls back to rho below).
//
// Arithmetic: Montgomery multiplication, one 64-bit limb for moduli below
// 2^64 and two limbs for moduli below 2^126 (the table's n are below 2^80).
// Primality of leaves: an odd-prime bitmap below 1e8, then deterministic
// Miller-Rabin with the smallest proven base sets (Sorenson-Webster 2017
// table: {2,3,5,7,11,13} suffices below 3,474,749,660,383, the twelve
// primes to 37 suffice below psi_12 ~ 3.19e23), all in Montgomery form.

#include <algorithm>
#include <bit>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "prime_bitmap.hpp"
#include "splitter.hpp"
#include "u128.hpp"

namespace cn {

// ---- Montgomery arithmetic, one limb (odd m < 2^64) ------------------------
struct Mont64 {
    u64 m, ninv, r1, r2;   // ninv = -m^-1 mod 2^64; r1 = R mod m; r2 = R^2 mod m, R = 2^64
    explicit Mont64(u64 mod) : m(mod) {
        u64 x = mod;                                   // Newton: x = m^-1 mod 2^64
        for (int i = 0; i < 6; ++i) x *= 2 - mod * x;
        ninv = 0 - x;
        r1 = static_cast<u64>((static_cast<u128>(1) << 64) % mod);
        r2 = static_cast<u64>((static_cast<u128>(r1) * r1) % mod);
    }
    u64 redc(u128 t) const {                           // t < m * 2^64
        const u64 u = static_cast<u64>(t) * ninv;
        const u128 s = t + static_cast<u128>(u) * m;   // low limb becomes 0; may carry out of 2^128
        u64 r = static_cast<u64>(s >> 64);
        if (s < t) r -= m;                             // carry out: true value is 2^64 + r, subtract m
        else if (r >= m) r -= m;
        return r;
    }
    u64 mul(u64 a, u64 b) const { return redc(static_cast<u128>(a) * b); }
    u64 to(u64 a) const { return mul(a % m, r2); }
    u64 from(u64 a) const { return redc(a); }
    u64 pow(u64 a_mont, u128 e) const {                // a in Montgomery form; result in Montgomery form
        u64 r = r1;
        while (e) {
            if (e & 1) r = mul(r, a_mont);
            a_mont = mul(a_mont, a_mont);
            e >>= 1;
        }
        return r;
    }
};

// ---- Montgomery arithmetic, two limbs (odd m < 2^126) ----------------------
struct Mont128 {
    u128 m, r1, r2;
    u64 ninv;
    explicit Mont128(u128 mod) : m(mod) {
        const u64 n0 = static_cast<u64>(mod);
        u64 x = n0;
        for (int i = 0; i < 6; ++i) x *= 2 - n0 * x;
        ninv = 0 - x;
        u128 r = 1 % mod;                              // 2^256 mod m by 256 doublings (r < 2^126 so 2r fits)
        for (int i = 0; i < 128; ++i) { r <<= 1; if (r >= mod) r -= mod; }
        r1 = r;
        for (int i = 0; i < 128; ++i) { r <<= 1; if (r >= mod) r -= mod; }
        r2 = r;
    }
    // a * b * 2^-128 mod m for a, b < m
    u128 mul(u128 a, u128 b) const {
        const u64 a0 = static_cast<u64>(a), a1 = static_cast<u64>(a >> 64);
        const u64 b0 = static_cast<u64>(b), b1 = static_cast<u64>(b >> 64);
        const u64 n0 = static_cast<u64>(m), n1 = static_cast<u64>(m >> 64);
        const u128 p00 = static_cast<u128>(a0) * b0, p01 = static_cast<u128>(a0) * b1;
        const u128 p10 = static_cast<u128>(a1) * b0, p11 = static_cast<u128>(a1) * b1;
        u64 t0 = static_cast<u64>(p00);
        u128 acc = (p00 >> 64) + static_cast<u64>(p01) + static_cast<u64>(p10);
        u64 t1 = static_cast<u64>(acc);
        acc = (acc >> 64) + (p01 >> 64) + (p10 >> 64) + p11;   // < 2^125
        u64 t2 = static_cast<u64>(acc), t3 = static_cast<u64>(acc >> 64), t4 = 0;
        // REDC limb 0
        u64 u = t0 * ninv;
        u128 p = static_cast<u128>(u) * n0;
        u64 lo = static_cast<u64>(p) + t0;                       // == 0 mod 2^64
        u128 c = (p >> 64) + (lo < static_cast<u64>(p) ? 1 : 0);
        c += static_cast<u128>(u) * n1 + t1;
        t1 = static_cast<u64>(c); c >>= 64;
        c += t2; t2 = static_cast<u64>(c); c >>= 64;
        c += t3; t3 = static_cast<u64>(c); t4 = static_cast<u64>(c >> 64);
        // REDC limb 1
        u = t1 * ninv;
        p = static_cast<u128>(u) * n0;
        lo = static_cast<u64>(p) + t1;
        c = (p >> 64) + (lo < static_cast<u64>(p) ? 1 : 0);
        c += static_cast<u128>(u) * n1 + t2;
        t2 = static_cast<u64>(c); c >>= 64;
        c += t3; t3 = static_cast<u64>(c); c >>= 64;
        t4 += static_cast<u64>(c);
        u128 res = (static_cast<u128>(t3) << 64) | t2;
        if (t4 || res >= m) res -= m;
        return res;
    }
    u128 to(u128 a) const { return mul(a % m, r2); }
    u128 from(u128 a) const { return mul(a, 1); }
    u128 pow(u128 a_mont, u128 e) const {
        u128 r = r1;
        while (e) {
            if (e & 1) r = mul(r, a_mont);
            a_mont = mul(a_mont, a_mont);
            e >>= 1;
        }
        return r;
    }
};

// ---- primality of leaves --------------------------------------------------
class LeafPrimality {
public:
    explicit LeafPrimality(u64 bitmap_limit = 100000000ull) : bm_(bitmap_limit) {}
    bool is_prime(u64 x) const {
        if (x < 2) return false;
        if (x == 2) return true;
        if (!(x & 1)) return false;
        if (x < bm_.limit()) return bm_.is_prime(x);
        static constexpr u64 kSmall[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (u64 p : kSmall) if (x % p == 0) return x == p;
        const int nb = x < 3474749660383ull ? 6 : 12;          // proven base counts
        const Mont64 mo(x);
        const int s = std::countr_zero(x - 1);
        const u64 d = (x - 1) >> s;
        const u64 one = mo.r1, minus_one = x - mo.r1;
        for (int k = 0; k < nb; ++k) {
            u64 y = mo.pow(mo.to(kSmall[k]), d);
            if (y == one || y == minus_one) continue;
            bool witness = true;
            for (int i = 1; i < s; ++i) {
                y = mo.mul(y, y);
                if (y == minus_one) { witness = false; break; }
            }
            if (witness) return false;
        }
        return true;
    }
private:
    PrimeBitmap bm_;
};

// ---- Algorithm 1 --------------------------------------------------------------
struct CnFactorResult {
    std::vector<u64> primes;   // ascending, proven
    int bases_used = 0;        // ladders computed (all composites, all bases; trial peels not counted)
    int peels = 0;             // composites split by trial division by the base itself
    bool fell_back = false;    // generic splitter had to finish a composite
};

class CnFactorizer {
public:
    CnFactorizer() : leaf_() {
        // the first 100 odd-prime bases after 2: consecutive primes as in the paper
        bases_.push_back(2);
        for (u64 p = 3; bases_.size() < 100; p += 2) if (leaf_.is_prime(p)) bases_.push_back(p);
    }

    // Factor x, a divisor of the Carmichael number n (x == n allowed), using
    // the exponent n - 1. Throws if the result does not multiply back to x or
    // a leaf cannot be proven prime.
    CnFactorResult factor(u128 x, u128 n) const {
        CnFactorResult out;
        if (x <= 1 || !(x & 1) || n < 3 || x >= (static_cast<u128>(1) << 126))
            throw std::runtime_error("cnfactor: bad input (need odd 1 < x < 2^126, n >= 3)");
        const u128 e = n - 1;
        const int s = splitter_detail::ctz128(e);
        const u128 np = e >> s;
        std::vector<u128> cur, next;
        classify_(x, out, cur);
        for (std::size_t bi = 0; bi < bases_.size() && !cur.empty(); ++bi) {
            const u64 b = bases_[bi];
            next.clear();
            for (const u128 m : cur) {
                std::vector<u128> parts;
                if (m % b == 0) { ++out.peels; parts.push_back(b); parts.push_back(m / b); }
                else {
                    ++out.bases_used;
                    if (m <= UINT64_MAX) ladder_(Mont64(static_cast<u64>(m)), m, b, np, s, parts);
                    else ladder_(Mont128(m), m, b, np, s, parts);
                }
                if (parts.size() <= 1) { next.push_back(m); continue; }   // no split from this base
                for (const u128 q : parts) classify_(q, out, next);
            }
            std::swap(cur, next);
        }
        for (const u128 m : cur) {   // bases exhausted: finish generically (never on table data)
            out.fell_back = true;
            SplitResult sr = split_u128(m);
            if (sr.uncertified != 1) throw std::runtime_error("cnfactor: uncertified leaf");
            for (const u128 p : sr.primes) {
                if (p > UINT64_MAX) throw std::runtime_error("cnfactor: prime above 2^64");
                out.primes.push_back(static_cast<u64>(p));
            }
        }
        std::sort(out.primes.begin(), out.primes.end());
        u128 prod = 1;
        for (const u64 p : out.primes) prod *= p;
        if (prod != x) throw std::runtime_error("cnfactor: product mismatch");
        return out;
    }

    const LeafPrimality& leaves() const { return leaf_; }

private:
    // Put q into the primes (proven) or the composites list.
    void classify_(u128 q, CnFactorResult& out, std::vector<u128>& comps) const {
        if (q == 1) return;
        if (q <= UINT64_MAX && leaf_.is_prime(static_cast<u64>(q))) { out.primes.push_back(static_cast<u64>(q)); return; }
        comps.push_back(q);
    }

    // One ladder mod m with base b: emits the parts of m it separates (at
    // least two entries iff a split happened).
    template <typename Mont, typename T>
    void ladder_impl_(const Mont& mo, u128 m, u64 b, u128 np, int s, std::vector<u128>& parts) const {
        const T minus_one = static_cast<T>(m) - mo.r1;
        T y = mo.pow(mo.to(static_cast<T>(b % m)), np);
        if (y == mo.r1) return;                         // b^n' = 1: every prime in the first factor, no split
        u128 remaining = m;
        std::vector<u128> found;
        for (int i = 0; i < s; ++i) {
            if (y == minus_one) break;                  // X[i] = -1 mod m: all remaining primes at this level
            const u128 g = gcd_(remaining, static_cast<u128>(mo.from(y)) + 1);
            if (g == remaining) break;                  // all remaining primes at this level
            if (g > 1) { found.push_back(g); remaining /= g; }
            y = mo.mul(y, y);
            if (y == mo.r1) break;                      // X[i+1] = 1: nothing more to find
        }
        if (found.empty()) return;
        parts = std::move(found);
        if (remaining > 1) parts.push_back(remaining);
    }
    template <typename Mont>
    void ladder_(const Mont& mo, u128 m, u64 b, u128 np, int s, std::vector<u128>& parts) const {
        if constexpr (std::is_same_v<Mont, Mont64>) ladder_impl_<Mont, u64>(mo, m, b, np, s, parts);
        else ladder_impl_<Mont, u128>(mo, m, b, np, s, parts);
    }
    static u128 gcd_(u128 a, u128 b) {
        if (a <= UINT64_MAX && b <= UINT64_MAX) return std::gcd(static_cast<u64>(a), static_cast<u64>(b));
        while (b) { const u128 t = a % b; a = b; b = t; }
        return a;
    }

    LeafPrimality leaf_;
    std::vector<u64> bases_;
};

} // namespace cn
