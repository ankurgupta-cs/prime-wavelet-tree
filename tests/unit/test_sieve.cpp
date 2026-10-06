#include <algorithm>
#include <cstdio>
#include <random>
#include <vector>

#include "segment.hpp"
#include "sieve.hpp"

using namespace cn;

static int failures = 0;
#define CHECK(cond)                                                     \
    do {                                                                \
        if (!(cond)) {                                                  \
            std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
            ++failures;                                                 \
        }                                                               \
    } while (0)

int main() {
    const u64 N = 10'000'000;
    LinearSieve s(N);

    // ---- prime counting checkpoints ----------------------------------
    const auto& P = s.primes();
    auto pi = [&](u64 x) -> std::size_t {
        return static_cast<std::size_t>(
            std::upper_bound(P.begin(), P.end(), static_cast<u32>(x)) - P.begin());
    };
    CHECK(pi(10) == 4);
    CHECK(pi(100) == 25);
    CHECK(pi(1'000) == 168);
    CHECK(pi(10'000) == 1229);
    CHECK(pi(100'000) == 9592);
    CHECK(pi(1'000'000) == 78498);
    CHECK(pi(10'000'000) == 664579);
    CHECK(P.front() == 2);
    CHECK(P.back() == 9'999'991);

    // ---- spf against brute force -------------------------------------
    for (u64 n = 2; n <= 20'000; ++n) {
        u64 f = n;
        for (u64 d = 2; d * d <= n; ++d)
            if (n % d == 0) { f = d; break; }
        CHECK(s.smallest_prime_factor(n) == f);
    }

    // ---- factorization round trips (all < N^2, so must complete) -----
    std::mt19937_64 rng(0xC0FFEE);
    for (int t = 0; t < 5000; ++t) {
        const u64 x = rng() % 999'999'999'999ULL + 2;   // [2, 1e12)
        const Factorization f = s.factor(x);
        CHECK(f.complete());
        u128 prod = 1;
        u64 prev = 0;
        for (const auto& pp : f.factors) {
            CHECK(pp.p > prev);
            prev = pp.p;
            CHECK(pp.e >= 1);
            if (pp.p <= s.limit()) CHECK(s.is_prime(pp.p));
            for (u32 i = 0; i < pp.e; ++i) prod *= pp.p;
        }
        CHECK(prod == x);
    }

    {   // pure prime power
        const Factorization f = s.factor(static_cast<u128>(1) << 60);
        CHECK(f.complete() && f.factors.size() == 1);
        CHECK(f.factors[0].p == 2 && f.factors[0].e == 60);
    }
    {   // prime above N^2 but within u64: certified by Miller-Rabin now
        const u128 q = parse_u128("999999999999999989");   // prime, ~1e18 > N^2
        const Factorization f = s.factor(q);
        CHECK(f.complete() && f.factors.size() == 1);
        CHECK(f.factors[0].p == 999999999999999989ULL && f.factors[0].e == 1);
    }
    {   // peel a small factor, leave a big prime cofactor uncertified
        const u128 m89 = (static_cast<u128>(1) << 89) - 1;   // Mersenne prime
        const Factorization f = s.factor(3 * m89);
        CHECK(!f.complete());
        CHECK(f.factors.size() == 1 && f.factors[0].p == 3 && f.factors[0].e == 1);
        CHECK(f.cofactor == m89);
    }
    {   // above 2^64 but fully sieve-factorable
        const u64 ps[5] = {1000003, 1000033, 1000037, 1000039, 1000081};
        u128 x = 1;
        for (u64 p : ps) x *= p;
        CHECK(x > static_cast<u128>(UINT64_MAX));
        const Factorization f = s.factor(x);
        CHECK(f.complete() && f.factors.size() == 5);
        for (int i = 0; i < 5; ++i) {
            CHECK(f.factors[i].p == ps[i]);
            CHECK(f.factors[i].e == 1);
        }
    }

    // ---- u128 decimal round trips ------------------------------------
    CHECK(to_string(parse_u128("0")) == "0");
    CHECK(to_string(parse_u128("340282366920938463463374607431768211455")) ==
          "340282366920938463463374607431768211455");
    {
        u128 e24 = 1;
        for (int i = 0; i < 24; ++i) e24 *= 10;
        CHECK(parse_u128("1000000000000000000000000") == e24);
    }
    {
        bool threw = false;
        try { parse_u128("340282366920938463463374607431768211456"); }
        catch (const std::out_of_range&) { threw = true; }
        CHECK(threw);
    }
    {
        bool threw = false;
        try { parse_u128("12a4"); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
    {
        bool threw = false;
        try { parse_u128(""); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // ---- Korselt / Carmichael ----------------------------------------
    const u64 known[] = {561, 1105, 1729, 2465, 2821, 6601, 8911,
                         41041, 62745, 825265, 321197185};
    for (u64 n : known) CHECK(s.is_carmichael(n) == Tri::Yes);

    const u64 nots[] = {2, 6, 91, 100, 341, 1122, 1683, 314721, 1000003};
    for (u64 n : nots) CHECK(s.is_carmichael(n) == Tri::No);

    {   // unknown when the cofactor cannot be certified
        const u128 m89 = (static_cast<u128>(1) << 89) - 1;
        CHECK(s.is_carmichael(3 * m89) == Tri::Unknown);
    }

    // ---- deterministic Miller-Rabin (mr64.hpp) -----------------------
    for (u64 n = 0; n <= 200'000; ++n)
        CHECK(is_prime_u64(n) == (n >= 2 && s.is_prime(n)));
    // Composite traps: strong pseudoprimes to small base sets.
    CHECK(!is_prime_u64(3215031751ULL));            // psp to bases 2,3,5,7
    CHECK(!is_prime_u64(341550071728321ULL));       // psi_7
    CHECK(!is_prime_u64(3825123056546413051ULL));   // psi_9
    CHECK(!is_prime_u64(4294967291ULL * 4294967279ULL));  // near-2^64 semiprime
    // Known primes at the top of the ranges.
    CHECK(is_prime_u64(4294967291ULL));             // 2^32 - 5
    CHECK(is_prime_u64(4294967279ULL));             // 2^32 - 17
    CHECK(is_prime_u64(2305843009213693951ULL));    // Mersenne 2^61 - 1
    CHECK(is_prime_u64(9223372036854775783ULL));    // largest prime < 2^63
    CHECK(is_prime_u64(18446744073709551557ULL));   // largest u64 prime
    CHECK(is_prime_u64(999999999989ULL));           // largest prime < 10^12

    // ---- factor() with tdiv_limit: honest bounded peel ---------------
    {
        const u128 x = static_cast<u128>(3) * 1000003 * 1000033;
        const Factorization f = s.factor(x, 1000);
        CHECK(!f.complete());
        CHECK(f.factors.size() == 1 && f.factors[0].p == 3);
        CHECK(f.cofactor == static_cast<u128>(1000003) * 1000033);
        const Factorization g = s.factor(x);        // unbounded: completes
        CHECK(g.complete() && g.factors.size() == 3);
    }
    {   // MR certifies a u64 prime cofactor far beyond N^2
        const u128 x = static_cast<u128>(3) * 18446744073709551557ULL;
        const Factorization f = s.factor(x);
        CHECK(f.complete() && f.factors.size() == 2);
        CHECK(f.factors[1].p == 18446744073709551557ULL);
    }
    {   // composite u64 cofactor with both factors > N stays honest
        const u64 q1 = 4294967291ULL, q2 = 4294967279ULL;
        const u128 x = static_cast<u128>(3) * q1 * q2;
        const Factorization f = s.factor(x);
        CHECK(!f.complete());
        CHECK(f.factors.size() == 1 && f.factors[0].p == 3);
        CHECK(f.cofactor == static_cast<u128>(q1) * q2);
    }

    // ---- real lines from the 10^24 table (tests/fixtures) --------------
    {   // d = 8: peels through the fold path, MR settles the last factor
        const u128 n = parse_u128("999999913302140456946241");
        CHECK(s.is_carmichael(n) == Tri::Yes);
        const Factorization f = s.factor(n);
        CHECK(f.complete() && f.factors.size() == 8);
        CHECK(f.factors[0].p == 19 && f.factors[7].p == 23116861);
    }
    {   // d = 10 (the last line of the table): all factors below N
        const u128 n = parse_u128("999999999855878641139521");
        const u64 expect[10] = {31, 37, 53, 71, 211, 313, 673, 1093, 2003, 2381};
        const Factorization f = s.factor(n);
        CHECK(f.complete() && f.factors.size() == 10);
        for (int i = 0; i < 10; ++i)
            CHECK(f.factors[i].p == expect[i] && f.factors[i].e == 1);
        CHECK(s.is_carmichael(n) == Tri::Yes);
    }

    // ---- segmented sieve ---------------------------------------------
    {
        SegmentedSieve seg(1'000'000'000'000ULL);
        CHECK(seg.base_limit() == 1'000'000);

        std::vector<u64> low;
        seg.for_primes_in(0, 10, [&](u64 p) { low.push_back(p); });
        CHECK(low.size() == 4 && low[0] == 2 && low[1] == 3 && low[2] == 5 && low[3] == 7);

        std::size_t c1e6 = 0;
        seg.for_primes_in(0, 1'000'000, [&](u64) { ++c1e6; });
        CHECK(c1e6 == 78498);

        // Window straddling 10^7 against LinearSieve ground truth.
        {
            std::vector<u64> got;
            seg.for_primes_in(9'900'000, 10'000'000, [&](u64 p) { got.push_back(p); });
            std::size_t want = 0;
            for (u64 n = 9'900'000; n < 10'000'000; ++n)
                if (s.is_prime(n)) {
                    CHECK(want < got.size() && got[want] == n);
                    ++want;
                }
            CHECK(got.size() == want);
        }

        // Spot window just below 10^12, cross-validated by Miller-Rabin.
        {
            const u64 a = 999'999'900'000ULL, b = 1'000'000'000'000ULL;
            std::vector<u64> got;
            seg.for_primes_in(a, b, [&](u64 p) { got.push_back(p); });
            std::size_t want = 0;
            for (u64 n = a | 1; n < b; n += 2)
                if (is_prime_u64(n)) {
                    CHECK(want < got.size() && got[want] == n);
                    ++want;
                }
            CHECK(got.size() == want);
            CHECK(!got.empty() && got.back() == 999'999'999'989ULL);
        }

        // Factored windows agree with LinearSieve::factor and the cofactor
        // contract (cofactor == 1 or prime).
        {
            std::size_t checked = 0;
            seg.for_factored_in(1'000'000, 1'004'096,
                [&](u64 n, const std::vector<PrimePower>& fs, u64 cof) {
                    u128 prod = cof;
                    u64 prev = 0;
                    for (const auto& pp : fs) {
                        CHECK(pp.p > prev);
                        prev = pp.p;
                        for (u32 i = 0; i < pp.e; ++i) prod *= pp.p;
                    }
                    CHECK(prod == n);
                    CHECK(cof == 1 || is_prime_u64(cof));
                    const Factorization ref = s.factor(n);
                    std::size_t nf = fs.size() + (cof != 1 ? 1 : 0);
                    CHECK(ref.factors.size() == nf);
                    ++checked;
                });
            CHECK(checked == 4096);
        }

        {   // range validation
            bool threw = false;
            try { seg.for_primes_in(0, 2'000'000'000'000ULL, [](u64) {}); }
            catch (const std::out_of_range&) { threw = true; }
            CHECK(threw);
        }
    }

    // ---- census below 10^4 / 10^6 / 10^7 (Pinch: 7 / 43 / 105) -------
    std::size_t c4 = 0, c6 = 0, c7 = 0;
    for (u64 n = 3; n < N; n += 2) {
        if (s.is_carmichael(n) == Tri::Yes) {
            ++c7;
            if (n < 10'000) ++c4;
            if (n < 1'000'000) ++c6;
        }
    }
    CHECK(c4 == 7);
    CHECK(c6 == 43);
    CHECK(c7 == 105);

    if (failures == 0)
        std::printf("all tests passed (sieve N = %llu, %zu primes)\n",
                    static_cast<unsigned long long>(N), s.prime_count());
    else
        std::printf("%d FAILURES\n", failures);
    return failures ? 1 : 0;
}
