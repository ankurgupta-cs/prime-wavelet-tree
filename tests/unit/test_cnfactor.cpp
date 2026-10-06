// test_cnfactor -- cnfactor.hpp: Montgomery arithmetic against the plain
// modular routines, leaf primality against the 12-base reference, and
// Shallue-Webster Algorithm 1 against the oracle's factor lists on the
// 637 fixture lines (n itself, and every cofactor r* = n / d_min with the
// exponent of n). Run from the repo root (fixtures in tests/fixtures/).

#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "cnfactor.hpp"
#include "lambda_bucket.hpp"
#include "mr64.hpp"
#include "splitter.hpp"
#include "u128.hpp"

using namespace cn;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); ++g_fail; } } while (0)

namespace {

u128 rand_u128(std::mt19937_64& rng, int bits) {
    u128 v = (static_cast<u128>(rng()) << 64) | rng();
    return bits >= 128 ? v : v & ((static_cast<u128>(1) << bits) - 1);
}

// Exact a*b mod m for any m < 2^127 by shift-and-add (the reference for
// moduli where splitter_detail::mulmod's limb loop would overflow).
u128 slow_mulmod(u128 a, u128 b, u128 m) {
    a %= m;
    u128 r = 0;
    for (int i = 127; i >= 0; --i) {
        r <<= 1; if (r >= m) r -= m;
        if ((b >> i) & 1) { r += a; if (r >= m) r -= m; }
    }
    return r;
}

void test_montgomery() {
    std::mt19937_64 rng(1);
    int bad64 = 0, bad128 = 0, badpow = 0;
    for (int t = 0; t < 200000; ++t) {
        const u64 m = (rng() >> (rng() % 40)) | 1;           // odd, up to 64 bits
        if (m < 3) continue;
        const Mont64 mo(m);
        const u64 a = rng() % m, b = rng() % m;
        if (mo.from(mo.mul(mo.to(a), mo.to(b))) != mulmod64(a, b, m)) ++bad64;
        if (t % 50 == 0) {
            const u128 e = rand_u128(rng, 80);
            if (mo.from(mo.pow(mo.to(a), e)) != static_cast<u64>(splitter_detail::powmod(a, e, m))) ++badpow;
        }
    }
    for (int t = 0; t < 200000; ++t) {
        const int bits = 65 + static_cast<int>(rng() % 61);   // 65..125 bits
        u128 m = rand_u128(rng, bits) | 1 | (static_cast<u128>(1) << (bits - 1));
        const Mont128 mo(m);
        const u128 a = rand_u128(rng, 128) % m, b = rand_u128(rng, 128) % m;
        if (mo.from(mo.mul(mo.to(a), mo.to(b))) != slow_mulmod(a, b, m)) ++bad128;
        if (t % 200 == 0 && m < (static_cast<u128>(1) << 95)) {   // the splitter's powmod is valid below 2^95
            const u128 e = rand_u128(rng, 80);
            if (mo.from(mo.pow(mo.to(a), e)) != splitter_detail::powmod(a, e, m)) ++badpow;
        }
    }
    // edge moduli: 2^64 - 1 (odd), 2^80 + 1, 2^126 - 1
    for (u128 m : {static_cast<u128>(~0ull), (static_cast<u128>(1) << 80) + 1, (static_cast<u128>(1) << 126) - 1}) {
        const Mont128 mo(m);
        for (int t = 0; t < 1000; ++t) {
            const u128 a = rand_u128(rng, 128) % m, b = rand_u128(rng, 128) % m;
            if (mo.from(mo.mul(mo.to(a), mo.to(b))) != slow_mulmod(a, b, m)) ++bad128;
        }
    }
    CHECK(bad64 == 0 && bad128 == 0 && badpow == 0, "Montgomery mismatches: 64-bit %d, 128-bit %d, pow %d", bad64, bad128, badpow);
    std::printf("Montgomery: 200k one-limb + 200k two-limb products and 2k powers agree with the reference\n");
}

void test_leaves() {
    const LeafPrimality lp;
    std::mt19937_64 rng(2);
    int bad = 0;
    for (u64 x = 0; x < 200000; ++x) if (lp.is_prime(x) != is_prime_u64(x)) ++bad;
    for (int t = 0; t < 300000; ++t) {
        const u64 x = rng() >> (rng() % 24);                  // up to 64 bits, skewed small
        if (lp.is_prime(x) != is_prime_u64(x)) ++bad;
    }
    // strong pseudoprimes at the base-set boundaries and squares of primes
    for (u64 x : {3215031751ull, 3474749660383ull, 341550071728321ull, 561ull, 1105ull,
                  4294967291ull * 4294967291ull /* largest 32-bit prime squared */, 4295098369ull /* 65537^2 */,
                  1000000007ull * 1000000009ull, 999999999989ull})
        if (lp.is_prime(x) != is_prime_u64(x)) ++bad;
    CHECK(bad == 0, "leaf primality disagreements: %d", bad);
    std::printf("leaf primality: 500k values and the base-set boundary pseudoprimes agree with 12-base MR\n");
}

void test_fixture(const std::string& dir) {
    const CnFactorizer cnf;
    int lines = 0, bad_n = 0, bad_r = 0, fallbacks = 0;
    double ladders_n = 0, ladders_r = 0;
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
            std::vector<u64> primes;
            while (ss >> tok) primes.push_back(static_cast<u64>(parse_u128(tok)));
            ++lines;
            const CnFactorResult fn = cnf.factor(n, n);
            if (fn.primes != primes) ++bad_n;
            if (fn.fell_back) ++fallbacks;
            ladders_n += fn.bases_used;
            const DminChoice c = choose_dmin(primes.data(), static_cast<int>(primes.size()), n);
            const u128 r = n / c.d;
            std::vector<u64> rp;
            for (const u64 p : primes) if (c.d % p != 0) rp.push_back(p);
            const CnFactorResult fr = cnf.factor(r, n);
            if (fr.primes != rp) ++bad_r;
            if (fr.fell_back) ++fallbacks;
            ladders_r += fr.bases_used;
        }
    }
    CHECK(bad_n == 0 && bad_r == 0, "fixture factorizations wrong: n %d, r* %d", bad_n, bad_r);
    std::printf("fixture: %d lines; Algorithm 1 reproduces every factor list for n (%.2f ladders each) and for r* (%.2f), %d fallbacks\n",
                lines, ladders_n / lines, ladders_r / lines, fallbacks);
    // a non-divisor whose primes do not obey p - 1 | n - 1: the ladders
    // happen to split it (base 5), so this exercises the non-divisor path
    {
        const CnFactorResult f = cnf.factor(91, 561);   // 7 - 1 does not divide 560
        const std::vector<u64> expect{7, 13};
        CHECK(f.primes == expect && !f.fell_back, "non-divisor still factored correctly");
    }
    // the generic fallback (the unconditional safety net behind the 100
    // bases): a composite with no Fermat structure at all survives every
    // ladder and is finished by the splitter, correctly
    {
        const CnFactorResult f = cnf.factor(1000003ull * 1000033ull, 561);
        const std::vector<u64> expect{1000003ull, 1000033ull};
        CHECK(f.primes == expect && f.fell_back && f.bases_used == 100,
              "fallback: primes ok %d, fell_back %d, ladders %d", f.primes == expect, f.fell_back, f.bases_used);
    }
    // a prime above 2^64 can never be a proven leaf: must throw, not misreport
    {
        const u128 big = (static_cast<u128>(1) << 64) + 13;   // prime
        bool threw = false;
        try { cnf.factor(big, big); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "prime above 2^64 rejected");
    }
    // an even input and an oversized modulus are refused
    {
        bool threw = false;
        try { cnf.factor(562, 561); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "even input rejected");
        threw = false;
        try { cnf.factor((static_cast<u128>(1) << 126) + 1, (static_cast<u128>(1) << 126) + 1); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "modulus >= 2^126 rejected");
    }
    // the classical small Carmichael numbers
    for (u64 n : {561ull, 1105ull, 1729ull, 2465ull, 2821ull, 6601ull, 8911ull, 41041ull, 62745ull, 63973ull,
                  75361ull, 101101ull, 126217ull, 172081ull, 188461ull, 278545ull, 340561ull, 449065ull}) {
        const CnFactorResult f = cnf.factor(n, n);
        u128 prod = 1;
        for (const u64 p : f.primes) { prod *= p; CHECK(cnf.leaves().is_prime(p), "leaf %llu prime", (unsigned long long)p); }
        CHECK(prod == n && f.primes.size() >= 3, "small Carmichael %llu", (unsigned long long)n);
    }
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
    test_montgomery();
    test_leaves();
    test_fixture(dir);
    if (g_fail) { std::printf("test_cnfactor: %d FAILURE(S)\n", g_fail); return 1; }
    std::printf("test_cnfactor: all OK\n");
    return 0;
}
