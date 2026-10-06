#pragma once
// sieve.hpp -- linear (Euler) sieve with a smallest-prime-factor table, plus
// factorization and a Korselt test, for the Carmichael compression project.
//
// LinearSieve(N):
//   time  O(N)   (every composite is struck exactly once, by its smallest prime)
//   space 4*(N+1) bytes for the spf table (uint32) + ~4*pi(N) bytes of primes
//
// Sizing facts for Carmichael numbers n <= 10^24:
//   * every prime factor p of a Carmichael n satisfies p < sqrt(n), so prime
//     factors range up to ~10^12; no in-RAM spf table covers that, and full
//     factorization of arbitrary table entries is a job for Pollard rho or
//     Miller-Rabin gcd splitting (a later module), not the sieve;
//   * every Carmichael n has at least 3 prime factors, so its smallest prime
//     factor is < n^(1/3) <= 10^8.
// Hence N = 10^8 guarantees the first peel for every entry in the 10^24 table
// and O(log x) factorization via spf-chasing for anything <= N. Cofactors that
// fit in 64 bits are settled definitively by deterministic Miller-Rabin
// (mr64.hpp), so certification extends to 2^64 ~ 1.8e19; only cofactors above
// 2^64 (e.g. q*r for 3-factor entries near 10^24) remain Unknown until the
// rho/Fermat-gcd splitter module lands.
//
// Sizing policy (decided 2026-07): N = 10^8 (400 MB) is the project default.
// Do not chase the uint32 cap of 4.29e9 -- 17 GB of SPF buys certification
// only to 1.8e19, which Miller-Rabin already provides for free.
//
// factor() takes an optional tdiv_limit: trial division stops at the first
// sieved prime exceeding it, returning the rest as an honest cofactor. This
// caps worst-case per-entry cost when processing the full 308M-entry table
// (peel small primes cheaply, hand surviving cofactors to a splitter).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "mr64.hpp"
#include "u128.hpp"

namespace cn {

struct PrimePower {
    u64 p;   // prime
    u32 e;   // exponent
};

// factors lists primes in strictly increasing order. cofactor == 1 means the
// factorization is complete. Otherwise cofactor holds the unfactored part; by
// construction every prime factor of cofactor exceeds the sieve limit N, and
// since cofactor > N^2 its primality is not certified here.
struct Factorization {
    std::vector<PrimePower> factors;
    u128 cofactor = 1;
    bool complete() const { return cofactor == 1; }
};

enum class Tri : std::uint8_t { No, Yes, Unknown };

class LinearSieve {
public:
    explicit LinearSieve(u64 limit) : n_(limit) {
        if (limit < 2)
            throw std::invalid_argument("LinearSieve: limit must be >= 2");
        if (limit > 0xFFFFFFFFull)
            throw std::invalid_argument("LinearSieve: limit must fit in uint32");
        spf_.assign(static_cast<std::size_t>(n_) + 1, 0);
        primes_.reserve(approx_pi(n_));
        for (u64 i = 2; i <= n_; ++i) {
            if (spf_[i] == 0) {
                spf_[i] = static_cast<u32>(i);
                primes_.push_back(static_cast<u32>(i));
            }
            const u32 s = spf_[i];
            for (u32 p : primes_) {
                if (p > s) break;
                const u64 ip = i * p;   // i, p < 2^32 so ip < 2^64: no overflow
                if (ip > n_) break;
                spf_[ip] = p;
            }
        }
    }

    u64 limit() const { return n_; }
    const std::vector<u32>& primes() const { return primes_; }
    std::size_t prime_count() const { return primes_.size(); }

    // Valid for 2 <= x <= limit().
    u32 smallest_prime_factor(u64 x) const {
        if (x < 2 || x > n_)
            throw std::out_of_range("smallest_prime_factor: x outside [2, N]");
        return spf_[x];
    }

    // Valid for x <= limit(); throws above that rather than guessing.
    bool is_prime(u64 x) const {
        if (x > n_) throw std::out_of_range("is_prime: x > N; use factor()");
        return x >= 2 && spf_[x] == x;
    }

    // Factor any u128 value. Trial division scans sieved primes up to
    // min(N, tdiv_limit); a 64-bit remainder is then settled definitively by
    // deterministic Miller-Rabin. cofactor != 1 only when a composite part
    // remains whose prime factors all exceed the scanned bound and which
    // either exceeds 2^64 (MR not applicable) or was proven composite by MR.
    // Values <= N complete via SPF chasing regardless of tdiv_limit.
    Factorization factor(u128 x, u64 tdiv_limit = UINT64_MAX) const {
        Factorization f;
        if (x <= 1) return f;
        if (x <= n_) { chase_(static_cast<u64>(x), f); return f; }

        std::size_t idx = 0;
        const std::size_t np = primes_.size();

        // Phase 1: x exceeds 64 bits, so no early exit by p^2 > x is
        // possible; scan the prime list. x % p is computed by folding the
        // high word (three u64 divisions beat the __umodti3 software path).
        while (x > static_cast<u128>(UINT64_MAX) && idx < np) {
            const u32 p = primes_[idx];
            if (p > tdiv_limit) { f.cofactor = x; return f; }
            const u64 hi = static_cast<u64>(x >> 64);
            const u64 lo = static_cast<u64>(x);
            const u64 b64modp = (UINT64_MAX % p + 1) % p;   // 2^64 mod p
            if (((hi % p) * b64modp + lo % p) % p == 0) {
                u32 e = 0;
                do { x /= p; ++e; } while (x % p == 0);
                f.factors.push_back({p, e});
            }
            ++idx;
        }
        if (x > static_cast<u128>(UINT64_MAX)) {   // prime list exhausted
            f.cofactor = x;
            return f;
        }

        // Phase 2: 64-bit arithmetic from here down.
        u64 y = static_cast<u64>(x);
        while (y > 1) {
            if (y <= n_) { chase_(y, f); return f; }
            if (idx >= np) break;
            const u64 p = primes_[idx];
            if (p > tdiv_limit) break;
            if (p * p > y) {                       // y has no factor <= sqrt(y)
                f.factors.push_back({y, 1});
                return f;
            }
            if (y % p == 0) {
                u32 e = 0;
                do { y /= p; ++e; } while (y % p == 0);
                f.factors.push_back({p, e});
            }
            ++idx;
        }
        if (y > 1) {
            if (is_prime_u64(y))
                f.factors.push_back({y, 1});       // certified prime (MR, u64)
            else
                f.cofactor = y;                    // composite, factors > bound
        }
        return f;
    }

    // Korselt's criterion: n odd, squarefree, with >= 3 prime factors, and
    // p - 1 | n - 1 for every prime p | n. A repeated factor or a Korselt
    // violation among the *found* factors settles No even when the
    // factorization is incomplete; Unknown only when neither disqualifier
    // appears and a composite cofactor remains.
    Tri is_carmichael(u128 n) const {
        if (n < 561 || (n & 1) == 0) return Tri::No;   // smallest CN is 561; CNs are odd
        const Factorization f = factor(n);
        const u128 m = n - 1;
        for (const PrimePower& pp : f.factors) {
            if (pp.e != 1) return Tri::No;             // must be squarefree
            if (m % (pp.p - 1) != 0) return Tri::No;   // Korselt divisibility
        }
        if (!f.complete()) return Tri::Unknown;
        if (f.factors.size() < 3) return Tri::No;
        return Tri::Yes;
    }

private:
    static std::size_t approx_pi(u64 n) {
        if (n < 17) return 8;
        const double nd = static_cast<double>(n);
        return static_cast<std::size_t>(nd / (std::log(nd) - 1.1)) + 16;
    }

    // Requires 2 <= x <= n_. O(log x) via the spf table.
    void chase_(u64 x, Factorization& out) const {
        while (x > 1) {
            const u32 p = spf_[x];
            u32 e = 0;
            do { x /= p; ++e; } while (x % p == 0);
            out.factors.push_back({p, e});
        }
    }

    u64 n_;
    std::vector<u32> spf_;
    std::vector<u32> primes_;
};

} // namespace cn
