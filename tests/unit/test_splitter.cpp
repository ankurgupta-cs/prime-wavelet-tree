// test_splitter.cpp -- tests for splitter.hpp (u128 rho + strong-Fermat-gcd).
// Run from the repo root: fixture lines in tests/fixtures/ are ground truth
// (factor lists straight from the Shallue-Webster table; never test against
// guesses).

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "sieve.hpp"
#include "splitter.hpp"
#include "u128.hpp"

using namespace cn;

static int failures = 0;

#define CHECK(cond, ...)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            std::printf("FAIL %s:%d  ", __FILE__, __LINE__);                 \
            std::printf(__VA_ARGS__);                                        \
            std::printf("\n");                                               \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

static u64 rng_state = 0x9e3779b97f4a7c15ull;
static u64 rnd64() {
    u64 x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return rng_state = x;
}

// Reference modular multiply: binary shift-add, valid for m < 2^127.
static u128 mulmod_ref(u128 a, u128 b, u128 m) {
    a %= m; b %= m;
    u128 r = 0;
    while (b) {
        if (b & 1) { r += a; if (r >= m) r -= m; }
        a += a; if (a >= m) a -= m;
        b >>= 1;
    }
    return r;
}

// First certified prime >= x (MR-12 proof; caller keeps x < PSI_12).
static u128 next_prime(u128 x) {
    if (x <= 2) return 2;
    if ((x & 1) == 0) ++x;
    while (!is_prime_u128(x)) x += 2;
    return x;
}

static void test_mulmod() {
    for (int i = 0; i < 2000; ++i) {
        const u128 m64 = (static_cast<u128>(rnd64()) | 2) + 1;   // odd-ish, <= 2^64
        const u128 hi = rnd64() & ((1ull << 30) - 1);
        const u128 m128 = (hi << 64 | rnd64()) | (static_cast<u128>(1) << 64);
        for (u128 m : {m64, m128}) {
            const u128 a = (static_cast<u128>(rnd64()) << 64 | rnd64()) % m;
            const u128 b = (static_cast<u128>(rnd64()) << 64 | rnd64()) % m;
            CHECK(splitter_detail::mulmod(a, b, m) == mulmod_ref(a, b, m),
                  "mulmod mismatch, m has %d-bit class", m == m64 ? 64 : 94);
        }
    }
    // Boundary modulus just under the 2^95 contract.
    const u128 m = (static_cast<u128>(1) << 95) - 45;
    for (int i = 0; i < 200; ++i) {
        const u128 a = (static_cast<u128>(rnd64()) << 64 | rnd64()) % m;
        const u128 b = (static_cast<u128>(rnd64()) << 64 | rnd64()) % m;
        CHECK(splitter_detail::mulmod(a, b, m) == mulmod_ref(a, b, m),
              "mulmod mismatch at 2^95 boundary");
    }
    bool threw = false;
    try { splitter_detail::powmod(2, 10, static_cast<u128>(1) << 95); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "powmod accepted modulus >= 2^95");
}

static void test_is_prime() {
    // Agreement with the certified u64 tester across magnitudes.
    for (int i = 0; i < 4000; ++i) {
        const u64 x = rnd64() >> (i % 40);
        CHECK(is_prime_u128(x) == is_prime_u64(x), "u64 agreement at x=%llu",
              static_cast<unsigned long long>(x));
    }
    // Primes straddling 2^64: verdicts below PSI_12 are proofs.
    const u128 p = next_prime((static_cast<u128>(1) << 64) + 1);
    CHECK(p > UINT64_MAX, "next_prime failed to cross 2^64");
    CHECK(is_prime_u128(p), "prime above 2^64 rejected");
    CHECK(!is_prime_u128(p * 3), "3p accepted as prime");
    for (u64 b : {2ull, 3ull, 5ull, 7ull})   // Fermat sanity on the same prime
        CHECK(splitter_detail::powmod(b, p - 1, p) == 1,
              "Fermat b^(p-1) != 1 mod p for b=%llu",
              static_cast<unsigned long long>(b));
}

static void expect_split(u128 n, const std::vector<u128>& want,
                         const char* label) {
    const SplitResult got = split_u128(n);
    CHECK(got.uncertified == 1, "%s: uncertified residue %s", label,
          to_string(got.uncertified).c_str());
    if (got.primes != want) {
        CHECK(false, "%s: wrong factorization of %s", label,
              to_string(n).c_str());
        return;
    }
    u128 prod = 1;
    for (u128 q : got.primes) prod *= q;
    CHECK(prod == n && got.uncertified == 1, "%s: product mismatch", label);
}

static void test_semiprimes_and_powers() {
    // Semiprimes with both factors beyond the sieve's 10^8 policy bound,
    // products from ~2^66 up to ~2^79.
    for (int bits = 33; bits <= 39; ++bits) {
        const u128 lo = static_cast<u128>(1) << bits;
        const u128 p = next_prime(lo + (rnd64() & 0xFFFF));
        const u128 q = next_prime(p + 2 + (rnd64() & 0xFFFF));
        expect_split(p * q, {p, q}, "semiprime");
    }
    // Prime square and cube (rho's classic hard case).
    const u128 p = next_prime(3000000019ull);
    expect_split(p * p, {p, p}, "square");
    expect_split(p * p * p, {p, p, p}, "cube");
    // d_min-shaped values: several certified primes < 2^40, product ~2^88
    // (must stay under the 2^95 modulus contract).
    std::vector<u128> ps;
    u128 prod = 1;
    for (u64 seed : {7ull, 11ull, 13ull}) {
        const u128 f = next_prime(15000 * seed);
        ps.push_back(f);
        prod *= f;
    }
    const u128 big = next_prime(static_cast<u128>(1) << 35);
    ps.push_back(big);
    prod *= big;
    std::sort(ps.begin(), ps.end());
    expect_split(prod, ps, "dmin-shaped");
    // Smooth + tiny-prime mix exercises the peel path.
    expect_split(static_cast<u128>(2 * 2 * 3 * 37) * 41,
                 {2, 2, 3, 37, 41}, "smooth");
}

static void test_uncertified_contract() {
    // A probable prime above PSI_12 must land in `uncertified`, not primes.
    const u128 pp = next_prime(PSI_12);
    const SplitResult s = split_u128(pp);
    CHECK(s.primes.empty() && s.uncertified == pp,
          "probable prime above PSI_12 was recorded as certified");
    Factorization f;
    f.cofactor = pp;
    finish_factorization(f);
    CHECK(f.factors.empty() && f.cofactor == pp && !f.complete(),
          "finish_factorization certified an unproven leaf");
}

static int test_fixtures(const std::string& dir, const LinearSieve& sieve) {
    int lines = 0;
    for (const char* name : {"fix_head.txt", "fix_mid.txt", "fix_tail.txt"}) {
        std::ifstream in(dir + "/" + name);
        CHECK(in.good(), "missing fixture %s", name);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string tok;
            ss >> tok;
            const u128 n = parse_u128(tok);
            std::vector<u128> want;
            while (ss >> tok) want.push_back(parse_u128(tok));
            expect_split(n, want, name);

            // Wire-in path: shallow sieve peel, splitter completes it.
            Factorization f = factor_full(sieve, n, /*tdiv_limit=*/1000);
            CHECK(f.complete(), "factor_full left a cofactor on %s",
                  to_string(n).c_str());
            CHECK(f.factors.size() == want.size(), "factor_full k mismatch");
            for (std::size_t i = 0; i < f.factors.size(); ++i) {
                CHECK(f.factors[i].e == 1, "table entry not squarefree?!");
                if (i < want.size())
                    CHECK(f.factors[i].p == want[i], "factor_full factor %zu", i);
            }
            ++lines;
        }
    }
    return lines;
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
    test_mulmod();
    test_is_prime();
    test_semiprimes_and_powers();
    test_uncertified_contract();
    LinearSieve sieve(1000000);
    const int lines = test_fixtures(dir, sieve);
    if (failures) {
        std::printf("%d splitter test failure(s)\n", failures);
        return 1;
    }
    std::printf("all splitter tests passed (%d fixture lines certified)\n",
                lines);
    return 0;
}
