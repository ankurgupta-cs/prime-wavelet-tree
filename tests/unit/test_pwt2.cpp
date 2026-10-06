// test_pwt2 -- PWT2 container tests (docs/PWT2_FORMAT.md). Run from the repo
// root (fixtures in tests/fixtures/).
//   1. word-based bit I/O against bic.hpp's bit-serial BitWriter/BitReader
//   2. wheel-210 index/value and the table universe (prime index below B,
//      wheel above) round trip
//   3. the global order: counts, ties to the larger prime, classes, vidx,
//      and the validation of a damaged order
//   4. fixture round trips (637 Carmichael numbers): d_min paths and full
//      factorizations, entropy-coded and byte-aligned; plus alternative
//      valid divisors (the largest prime set that still satisfies
//      d*lambda(d) > n) so that divisor mode is not tied to d_min
//   5. synthetic trees: root child count far above the escape (5,000),
//      primes above the sieve limit in the rank table (both universes),
//      a terminal-internal node, deterministic output
//   6. corruption battery: payload flips (hash must catch them), flips
//      below a recomputed payload hash (structure or targets must catch
//      them), truncation, extension, crafted geometry; never a silent
//      wrong set

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "bic.hpp"
#include "lambda_bucket.hpp"
#include "mr64.hpp"
#include "pwt2.hpp"
#include "u128.hpp"

using namespace cn;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); ++g_fail; } } while (0)

namespace {

struct Fixture {
    std::vector<std::vector<u64>> dmin_paths, full_paths, alt_paths;
    std::vector<u128> ns;   // sorted
};

u128 lambda_of(const std::vector<u64>& ps) {
    u128 lam = 1;
    for (const u64 p : ps) lam = lcm128(lam, p - 1);
    return lam;
}

Fixture load_fixture(const std::string& dir) {
    Fixture fx;
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
            u64 primes[14];
            int k = 0;
            while (ss >> tok) primes[k++] = static_cast<u64>(parse_u128(tok));
            const DminChoice c = choose_dmin(primes, k, n);
            CHECK(dmin_value(c.d, c.lambda) == n, "fixture reconstruction");
            std::vector<u64> dp, fp(primes, primes + k);
            for (int i = 0; i < k; ++i) if (c.d % primes[i] == 0) dp.push_back(primes[i]);
            // an alternative valid divisor: drop the smallest primes of n as
            // long as prod * lambda stays above n (differs from d_min often)
            std::vector<u64> alt(fp);
            for (;;) {
                bool dropped = false;
                for (std::size_t i = 0; i < alt.size() && alt.size() > 1; ++i) {
                    std::vector<u64> t(alt);
                    t.erase(t.begin() + static_cast<std::ptrdiff_t>(i));
                    u128 d = 1;
                    for (const u64 p : t) d *= p;
                    const u128 lam = lambda_of(t);
                    if (d <= n / lam) continue;   // need d * lam > n, i.e. d > floor(n / lam)
                    alt = t;
                    dropped = true;
                    break;
                }
                if (!dropped) break;
            }
            {
                u128 d = 1;
                for (const u64 p : alt) d *= p;
                CHECK(dmin_value(d, lambda_of(alt)) == n, "alternative divisor reconstructs n");
            }
            fx.dmin_paths.push_back(dp);
            fx.full_paths.push_back(fp);
            fx.alt_paths.push_back(alt);
            fx.ns.push_back(n);
        }
    }
    std::sort(fx.ns.begin(), fx.ns.end());
    return fx;
}

PwtTargets targets_of(const std::vector<u128>& ns_sorted) {
    PwtTargets t;
    Sha256 sha;
    u128 sum = 0;
    for (const u128 n : ns_sorted) {
        std::uint8_t le[16];
        for (int b = 0; b < 16; ++b) le[b] = static_cast<std::uint8_t>(n >> (8 * b));
        sha.update(le, 16);
        sum += n;
    }
    t.sha_nset = sha.finish();
    t.total_check = sum;
    t.record_count = ns_sorted.size();
    return t;
}

std::vector<u128> decode_raw(std::vector<std::uint8_t> bytes, bool* targets_ok = nullptr) {
    Pwt2Reader rd(std::move(bytes));
    std::vector<u128> ns;
    u128 sum = 0;
    rd.for_each_n([&](u128 n) { ns.push_back(n); sum += n; });
    std::sort(ns.begin(), ns.end());
    Sha256 sha;
    for (const u128 n : ns) {
        std::uint8_t le[16];
        for (int b = 0; b < 16; ++b) le[b] = static_cast<std::uint8_t>(n >> (8 * b));
        sha.update(le, 16);
    }
    if (targets_ok) *targets_ok = sha.finish() == rd.targets().sha_nset && sum == rd.targets().total_check
                                  && ns.size() == rd.targets().record_count;
    return ns;
}
std::vector<u128> decode_all(std::vector<std::uint8_t> bytes) {
    bool ok = false;
    std::vector<u128> ns = decode_raw(std::move(bytes), &ok);
    if (!ok) throw std::runtime_error("footer targets mismatch");
    return ns;
}

std::vector<std::uint8_t> encode(const std::vector<std::vector<u64>>& paths, const std::vector<u128>& ns_sorted,
                                 bool full, bool naive, u64 B = 100000000ull, Pwt2EncodeStats* st = nullptr, int ctx_set = 1) {
    const Pwt2Order ord = pwt2_order_from_paths(paths);
    const Pwt2VectorSource src(paths, ord);
    Pwt2Encoder enc(full, naive, B, naive ? 1 : ctx_set);
    return enc.encode(src, ord, targets_of(ns_sorted), st);
}

void recompute_payload_sha(std::vector<std::uint8_t>& b) {
    const std::size_t foff = b.size() - kPwt2FooterBytes;
    Sha256 sha;
    sha.update(b.data(), foff);
    const auto h = sha.finish();
    std::copy(h.begin(), h.end(), b.begin() + static_cast<std::ptrdiff_t>(foff));
}

void test_bits() {
    std::mt19937_64 rng(11);
    int bad = 0;
    for (int t = 0; t < 200; ++t) {
        std::vector<std::pair<u64, int>> items;
        std::vector<std::uint8_t> a, b;
        BitWriter slow(a);
        FastBitWriter fast(b);
        const int n = 1 + static_cast<int>(rng() % 3000);
        for (int i = 0; i < n; ++i) {
            const int w = static_cast<int>(rng() % 57);
            const u64 v = w == 0 ? 0 : (rng() & ((w == 64) ? ~0ull : ((1ull << w) - 1)));
            items.emplace_back(v, w);
            slow.put(v, w);
            fast.put(v, w);
        }
        slow.align();
        fast.align();
        if (a != b) ++bad;
        FastBitReader fr(b.data(), b.size());
        BitReader sr(a.data(), a.size());
        for (const auto& [v, w] : items) {
            if (fr.get(w) != v) ++bad;
            if (static_cast<u64>(sr.get(w)) != v) ++bad;
        }
    }
    CHECK(bad == 0, "bit I/O mismatches: %d", bad);
    bool threw = false;
    std::vector<std::uint8_t> one{0xAB};
    FastBitReader r(one.data(), one.size());
    r.get(5);
    try { r.get(4); } catch (const std::exception&) { threw = true; }
    CHECK(threw, "bit reader overrun rejected");
    std::printf("bit I/O: word-based writer/reader byte-identical to the bit-serial ones on 200 random streams\n");
}

void test_fastnum() {
    std::mt19937_64 rng(21);
    int bad = 0;
    for (int t = 0; t < 2000000; ++t) {
        const u64 a = rng() >> (rng() % 64), b = rng() >> (rng() % 64);
        if (gcd64_bin(a, b) != gcd64(a, b)) ++bad;
        const u64 m = ((rng() >> (rng() % 63)) | 2) & ~u64(1);   // even, 2 .. 2^64 - 2
        const u64 x = (rng() % m) | 1;
        if (gcd64(x, m) != 1) continue;
        const u64 y = inv_mod_lambda(x, m);
        if (y >= m || (static_cast<u128>(x) * y) % m != 1 % m) ++bad;
        if (m < (u64(1) << 62) && y != inv_mod64(x, m)) ++bad;
    }
    bool threw = false;
    try { (void)inv_mod_lambda(6, 12); } catch (const std::exception&) { threw = true; }
    CHECK(threw, "non-invertible rejected");
    CHECK(bad == 0, "fastnum mismatches: %d", bad);
    std::puts("fastnum: binary gcd and the lambda inverse agree with exact 128-bit checks on 2M random inputs, all even moduli");
}

void test_bits_quot() {
    int bad = 0;
    for (u64 X = 0; X < 3000; ++X)
        for (u64 R = 1; R < 700; ++R)
            if (Pwt2Keys2::bits_quot(X, R) != static_cast<u64>(pwt_bits_of(X / R))) ++bad;
    std::mt19937_64 rng(31);
    for (int t = 0; t < 3000000; ++t) {
        const u64 X = rng() >> (rng() % 64), R = (rng() >> (rng() % 64)) | 1;
        const u64 R2 = R + (rng() & 1);   // even values too
        if (Pwt2Keys2::bits_quot(X, R) != static_cast<u64>(pwt_bits_of(X / R))) ++bad;
        if (R2 && Pwt2Keys2::bits_quot(X, R2) != static_cast<u64>(pwt_bits_of(X / R2))) ++bad;
    }
    CHECK(bad == 0, "bits_quot mismatches: %d", bad);
    std::puts("bits_quot: bits(floor(X / R)) without division agrees on 2.1M exhaustive + 6M random pairs");
}

void test_universe() {
    const Pwt2Wheel w;
    int bad = 0;
    for (u64 x = 1; x < 100000; ++x) {
        if (w.ridx[x % kPwt2Wheel] < 0) continue;
        if (w.value(w.index(x)) != x) ++bad;
    }
    const u64 B = 1000003;
    const PrimeBitmap bm(B);
    const Pwt2Universe uni(bm, B);
    u64 prev_pos = 0;
    bool first = true;
    for (u64 p = 3; p < 3 * B; p += 2) {
        if (!is_prime_u64(p)) continue;
        const u64 pos = uni.pos(p);
        if (uni.value(pos) != p) ++bad;
        if (!first && pos <= prev_pos) ++bad;
        if ((p < B) != !uni.high(pos)) ++bad;
        prev_pos = pos;
        first = false;
    }
    bool threw = false;
    try { (void)uni.pos(9); } catch (const std::exception&) { threw = true; }
    CHECK(threw, "composite below B rejected by the universe");
    CHECK(bad == 0, "universe mismatches: %d", bad);
    std::printf("universe: wheel-210 round trip; prime index below B, wheel above, strictly increasing across B\n");
}

void test_order() {
    // counts: 19 x3, 31 x3, 7 x2, 13 x2, 101 x1 -> ranks 31 19 | 13 7 | 101
    const std::vector<std::vector<u64>> paths{{7, 19, 31}, {13, 19, 31}, {7, 13, 19, 31}, {101}};
    const Pwt2Order o = pwt2_order_from_paths(paths);
    CHECK(o.m() == 5, "m");
    CHECK(o.prime[1] == 31 && o.prime[2] == 19 && o.prime[3] == 13 && o.prime[4] == 7 && o.prime[5] == 101,
          "rank order (count desc, ties to the larger prime)");
    CHECK(o.class_size.size() == 3 && o.class_size[0] == 2 && o.class_size[1] == 2 && o.class_size[2] == 1, "classes");
    CHECK(o.vidx[4] == 1 && o.vidx[3] == 2 && o.vidx[2] == 3 && o.vidx[1] == 4 && o.vidx[5] == 5, "value index");
    int threw = 0;
    { Pwt2Order d = o; std::swap(d.prime[1], d.prime[2]); try { d.finish(); } catch (const std::exception&) { ++threw; } }
    { Pwt2Order d = o; d.prime[5] = 19; d.class_size = {2, 2, 1}; try { d.finish(); } catch (const std::exception&) { ++threw; } }
    { Pwt2Order d = o; d.class_size = {2, 2}; try { d.finish(); } catch (const std::exception&) { ++threw; } }
    { Pwt2Order d = o; d.prime[5] = 100; try { d.finish(); } catch (const std::exception&) { ++threw; } }
    CHECK(threw == 4, "damaged orders rejected: %d of 4", threw);
    std::printf("order: counts, ties, classes and value index correct; damaged orders rejected\n");
}

void test_fixture(const Fixture& fx) {
    struct Case { const char* name; const std::vector<std::vector<u64>>* paths; bool full; };
    const Case cases[] = {{"d_min", &fx.dmin_paths, false}, {"full", &fx.full_paths, true}, {"alternative divisor", &fx.alt_paths, false}};
    for (const Case& c : cases) {
        int alt_differs = 0;
        if (c.paths == &fx.alt_paths)
            for (std::size_t i = 0; i < fx.alt_paths.size(); ++i) {
                std::vector<u64> a = fx.alt_paths[i], d = fx.dmin_paths[i];
                std::sort(a.begin(), a.end()); std::sort(d.begin(), d.end());
                if (a != d) ++alt_differs;
            }
        for (const bool naive : {false, true}) {
            Pwt2EncodeStats st;
            std::vector<std::uint8_t> bytes;
            try { bytes = encode(*c.paths, fx.ns, c.full, naive, 100000000ull, &st); }
            catch (const std::exception& e) { CHECK(false, "%s %s encode threw: %s", c.name, naive ? "naive" : "rANS", e.what()); continue; }
            std::vector<u128> got;
            try { got = decode_all(bytes); }
            catch (const std::exception& e) { CHECK(false, "%s %s decode threw: %s", c.name, naive ? "naive" : "rANS", e.what()); continue; }
            CHECK(got == fx.ns, "%s %s round trip", c.name, naive ? "naive" : "rANS");
            const auto again = encode(*c.paths, fx.ns, c.full, naive);
            CHECK(again == bytes, "%s %s deterministic", c.name, naive ? "naive" : "rANS");
            std::printf("fixture %-20s %-5s %zu records: %zu bytes, %llu nodes, %llu symbols (table %llu)%s\n", c.name,
                        naive ? "naive" : "rANS", got.size(), bytes.size(), (unsigned long long)st.nodes,
                        (unsigned long long)st.symbols, (unsigned long long)st.table_symbols,
                        c.paths == &fx.alt_paths ? (" -- alternative differs from d_min on " + std::to_string(alt_differs) + " records").c_str() : "");
        }
    }
}

void test_synthetic() {
    // full mode (n = product), so paths need not be Carmichael
    std::mt19937_64 rng(5);
    std::vector<u64> big;   // primes above the sieve limit B = 1000003
    for (u64 x = 1000003 + 2; big.size() < 60; x += 2) if (is_prime_u64(x) && rng() % 7 == 0) big.push_back(x);
    std::vector<u64> small;
    for (u64 x = 3; small.size() < 6000; x += 2) if (is_prime_u64(x)) small.push_back(x);
    std::vector<std::vector<u64>> paths;
    for (int i = 0; i < 5000; ++i) paths.push_back({small[static_cast<std::size_t>(i) + 50]});   // 5,000 single-prime paths: root escape
    paths.push_back({3, 5});                  // 3 and 5 are rare -> deep
    paths.push_back({3, 5, 7});               // terminal-internal node {3,5}
    for (std::size_t i = 0; i < big.size(); ++i) {
        paths.push_back({small[i % 40], big[i]});
        paths.push_back({small[(i + 3) % 40], small[(i + 11) % 40 + 40], big[i]});
    }
    std::sort(paths.begin(), paths.end());
    paths.erase(std::unique(paths.begin(), paths.end()), paths.end());
    std::vector<u128> ns;
    for (const auto& p : paths) { u128 prod = 1; for (const u64 q : p) prod *= q; ns.push_back(prod); }
    std::sort(ns.begin(), ns.end());
    CHECK(std::adjacent_find(ns.begin(), ns.end()) == ns.end(), "synthetic products distinct");
    for (const bool naive : {false, true}) {
        Pwt2EncodeStats st;
        const auto bytes = encode(paths, ns, true, naive, 1000003ull, &st);
        std::vector<u128> got;
        try { got = decode_all(bytes); } catch (const std::exception& e) { CHECK(false, "synthetic decode: %s", e.what()); }
        CHECK(got == ns, "synthetic %s round trip", naive ? "naive" : "rANS");
        if (!naive) CHECK(st.escapes >= 1, "root escape exercised (%llu)", (unsigned long long)st.escapes);
        std::printf("synthetic %-5s %zu records: %zu bytes, %llu escapes, rank table with %zu primes above B\n",
                    naive ? "naive" : "rANS", ns.size(), bytes.size(), (unsigned long long)st.escapes, big.size());
    }
    // a rank outside 1..m or a duplicate path must be refused by the encoder
    {
        const Pwt2Order ord = pwt2_order_from_paths(paths);
        bool threw = false;
        try { Pwt2VectorSource bad({{3, 5}, {5, 3}}, ord); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "duplicate path refused");
    }
}

// Context set 2: fixture and synthetic round trips (both modes),
// determinism, the context-set header byte, and the corruption battery.
void test_ctx2(const Fixture& fx) {
    struct Case { const char* name; const std::vector<std::vector<u64>>* paths; bool full; };
    const Case cases[] = {{"d_min", &fx.dmin_paths, false}, {"full", &fx.full_paths, true}, {"alternative divisor", &fx.alt_paths, false}};
    for (const int cs : {2, 3})
    for (const Case& c : cases) {
        Pwt2EncodeStats st;
        std::vector<std::uint8_t> b;
        try { b = encode(*c.paths, fx.ns, c.full, false, 100000000ull, &st, cs); }
        catch (const std::exception& e) { CHECK(false, "set-2 %s encode threw: %s", c.name, e.what()); continue; }
        CHECK(b[17] == cs, "set-%d header byte", cs);
        std::vector<u128> got;
        try { got = decode_all(b); } catch (const std::exception& e) { CHECK(false, "set-2 %s decode threw: %s", c.name, e.what()); continue; }
        CHECK(got == fx.ns, "set-2 %s round trip", c.name);
        CHECK(encode(*c.paths, fx.ns, c.full, false, 100000000ull, nullptr, cs) == b, "set-%d %s deterministic", cs, c.name);
        std::printf("set %d %-20s %zu records: %zu bytes, %llu models, %llu symbols\n", cs, c.name, got.size(), b.size(),
                    (unsigned long long)st.contexts, (unsigned long long)st.symbols);
    }
    {   // synthetic: 5,000 root children (tail + gamma escape), primes above B, a terminal-internal node
        std::vector<u64> big, small;
        for (u64 x = 1000003 + 2; big.size() < 40; x += 2) if (is_prime_u64(x)) big.push_back(x);
        for (u64 x = 3; small.size() < 6000; x += 2) if (is_prime_u64(x)) small.push_back(x);
        std::vector<std::vector<u64>> paths;
        for (int i = 0; i < 5000; ++i) paths.push_back({small[static_cast<std::size_t>(i) + 50]});
        for (int i = 0; i < 70; ++i) paths.push_back({small[1], small[static_cast<std::size_t>(i) + 3000]});   // a node with 70 children: tail symbol
        paths.push_back({3, 5}); paths.push_back({3, 5, 7});
        for (std::size_t i = 0; i < big.size(); ++i) paths.push_back({small[i % 30], big[i]});
        std::sort(paths.begin(), paths.end());
        paths.erase(std::unique(paths.begin(), paths.end()), paths.end());
        std::vector<u128> ns;
        for (const auto& p : paths) { u128 prod = 1; for (const u64 q : p) prod *= q; ns.push_back(prod); }
        std::sort(ns.begin(), ns.end());
        Pwt2EncodeStats st;
        const auto b = encode(paths, ns, true, false, 1000003ull, &st, 3);
        std::vector<u128> got;
        try { got = decode_all(b); } catch (const std::exception& e) { CHECK(false, "set-2 synthetic decode: %s", e.what()); }
        CHECK(got == ns, "set-2 synthetic round trip");
        CHECK(st.escapes >= 1, "set-2 gamma escape exercised");
        std::printf("set 2 synthetic %zu records: %zu bytes, %llu escapes\n", ns.size(), b.size(), (unsigned long long)st.escapes);
    }
    {   // set 3 with mantissa models actually kept: under parent 3 the children are every third
        // count-1 prime (rank gaps 3 -> L = 2, mantissa always 1), under parent 5 the rest
        // (gaps 1, 2 alternating); full mode, n = product
        std::vector<u64> pr;
        for (u64 x = 1001; pr.size() < 30000; x += 2) if (is_prime_u64(x)) pr.push_back(x);
        std::vector<std::vector<u64>> paths;
        for (std::size_t j = 0; j + 2 < pr.size(); j += 3) {
            paths.push_back({3, pr[j]});
            paths.push_back({5, pr[j + 1]});
            paths.push_back({5, pr[j + 2]});
        }
        std::vector<u128> ns;
        for (const auto& q : paths) ns.push_back(static_cast<u128>(q[0]) * q[1]);
        std::sort(ns.begin(), ns.end());
        for (const int cs : {2, 3}) {
            Pwt2EncodeStats st;
            const auto b = encode(paths, ns, true, false, 100000000ull, &st, cs);
            std::vector<u128> got;
            try { got = decode_all(b); } catch (const std::exception& e) { CHECK(false, "set-%d mantissa synthetic decode: %s", cs, e.what()); }
            CHECK(got == ns, "set-%d mantissa synthetic round trip", cs);
            if (cs == 3) CHECK(st.class_contexts[kC2Mant] > 0, "set 3 keeps mantissa models (%llu)", (unsigned long long)st.class_contexts[kC2Mant]);
            std::printf("set %d mantissa synthetic %zu records: %zu bytes, %llu mantissa models\n", cs, ns.size(), b.size(),
                        (unsigned long long)st.class_contexts[kC2Mant]);
            if (cs == 3) {   // corruption on a container whose mantissa class is populated
                std::mt19937_64 rng(77);
                int wrong = 0, caught = 0, benign = 0;
                for (int t = 0; t < 400; ++t) {
                    auto bb = b;
                    bb[kPwt2HeaderBytes + rng() % (bb.size() - kPwt2FooterBytes - kPwt2HeaderBytes)] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
                    recompute_payload_sha(bb);
                    try { const auto g2 = decode_all(bb); if (g2 != ns) ++wrong; else ++benign; } catch (const std::exception&) { ++caught; }
                }
                CHECK(wrong == 0, "set-3 mantissa corruption: %d wrong sets", wrong);
                std::printf("set 3 mantissa corruption: %d rejected, %d benign, 0 silent wrong sets\n", caught, benign);
            }
        }
    }
    {   // the byte-aligned variant has no context set 2
        bool threw = false;
        try { Pwt2Encoder e(false, true, 100000000ull, 2); (void)e; } catch (const std::exception&) { threw = true; }
        CHECK(threw, "naive + set 2 refused");
    }
    {   // corruption battery on the set-2 d_min fixture container
        std::mt19937_64 rng(1234);
        const auto good = encode(fx.full_paths, fx.ns, true, false, 100000000ull, nullptr, 3);
        int silent = 0, wrong = 0, caught = 0, benign = 0;
        for (int t = 0; t < 400; ++t) {
            auto b = good;
            b[rng() % (b.size() - kPwt2FooterBytes)] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
            try { decode_all(b); ++silent; } catch (const std::exception&) { ++caught; }
        }
        for (int t = 0; t < 400; ++t) {
            auto b = good;
            b[kPwt2HeaderBytes + rng() % (b.size() - kPwt2FooterBytes - kPwt2HeaderBytes)] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
            recompute_payload_sha(b);
            try { const auto got = decode_all(b); if (got != fx.ns) ++wrong; else ++benign; } catch (const std::exception&) { ++caught; }
        }
        for (std::size_t pos = 0; pos < kPwt2HeaderBytes; ++pos)
            for (int bit = 0; bit < 8; ++bit) {
                auto b = good;
                b[pos] ^= static_cast<std::uint8_t>(1u << bit);
                recompute_payload_sha(b);
                try { const auto got = decode_all(b); if (got != fx.ns) ++wrong; else ++benign; } catch (const std::exception&) { ++caught; }
            }
        CHECK(silent == 0 && wrong == 0, "set-2 corruption: %d silent, %d wrong sets", silent, wrong);
        std::printf("set 2 corruption: %d rejected, %d benign, 0 silent wrong sets\n", caught, benign);
    }
}

// Crafted containers (clean-room encoder + mutation), each with its
// expected verdict: ACCEPT = decodes
// and matches its footer targets, REJECT = the reader, the reconstruction
// or the target check refuses it.
void test_crafted(const std::string& dir) {
    std::ifstream man(dir + "/pwt2_craft/MANIFEST.txt");
    CHECK(man.good(), "missing %s/pwt2_craft/MANIFEST.txt", dir.c_str());
    std::string name, want;
    int n = 0, bad = 0;
    while (man >> name >> want) {
        std::ifstream in(dir + "/pwt2_craft/" + name, std::ios::binary);
        std::vector<std::uint8_t> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        bool accepted = false;
        std::string why;
        try {
            bool ok = false;
            const auto ns = decode_raw(b, &ok);
            accepted = ok && std::adjacent_find(ns.begin(), ns.end()) == ns.end();
            if (!accepted) why = ok ? "duplicate values" : "footer targets";
        } catch (const std::exception& e) { why = e.what(); }
        ++n;
        if (accepted != (want == "ACCEPT")) { ++bad; std::printf("  crafted %s: expected %s, got %s (%s)\n", name.c_str(), want.c_str(), accepted ? "ACCEPT" : "REJECT", why.c_str()); }
    }
    CHECK(n > 0 && bad == 0, "crafted containers: %d of %d wrong verdicts", bad, n);
    std::printf("crafted: %d crafted containers (clean-room and mutation cases) all get their expected verdict\n", n);
}

void test_corruption(const Fixture& fx) {
    std::mt19937_64 rng(99);
    for (const bool naive : {false, true}) {
        const auto good = encode(fx.dmin_paths, fx.ns, false, naive);
        int silent = 0, caught = 0;
        // (a) payload flips: the payload hash must reject every one
        for (int t = 0; t < 400; ++t) {
            auto b = good;
            const std::size_t pos = rng() % (b.size() - kPwt2FooterBytes);
            b[pos] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
            try { decode_all(b); ++silent; } catch (const std::exception&) { ++caught; }
        }
        CHECK(silent == 0, "%s payload flips decoded silently: %d", naive ? "naive" : "rANS", silent);
        // (b) flips below a recomputed payload hash: structure or targets must catch them
        int wrong = 0, benign = 0;
        for (int t = 0; t < 400; ++t) {
            auto b = good;
            const std::size_t pos = kPwt2HeaderBytes + rng() % (b.size() - kPwt2FooterBytes - kPwt2HeaderBytes);
            b[pos] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
            recompute_payload_sha(b);
            try {
                const auto got = decode_all(b);
                if (got != fx.ns) ++wrong; else ++benign;
            } catch (const std::exception&) { ++caught; }
        }
        CHECK(wrong == 0, "%s rehashed flips decoded to a WRONG set silently: %d", naive ? "naive" : "rANS", wrong);
        // (c) header flips with the hash recomputed
        for (std::size_t pos = 0; pos < kPwt2HeaderBytes; ++pos)
            for (int bit = 0; bit < 8; ++bit) {
                auto b = good;
                b[pos] ^= static_cast<std::uint8_t>(1u << bit);
                recompute_payload_sha(b);
                try { const auto got = decode_all(b); if (got != fx.ns) ++wrong; else ++benign; }
                catch (const std::exception&) { ++caught; }
            }
        CHECK(wrong == 0, "%s header flips decoded to a WRONG set: %d", naive ? "naive" : "rANS", wrong);
        // (d) truncation and extension
        int trunc_ok = 0;
        for (const std::size_t cut : {std::size_t(1), std::size_t(17), std::size_t(88), std::size_t(200), good.size() / 2}) {
            if (cut >= good.size()) continue;
            auto b = good;
            b.resize(b.size() - cut);
            try { decode_all(b); } catch (const std::exception&) { ++trunc_ok; }
        }
        {
            auto b = good;
            b.push_back(0);
            try { decode_all(b); } catch (const std::exception&) { ++trunc_ok; }
        }
        CHECK(trunc_ok == 6, "%s truncations/extension rejected: %d of 6", naive ? "naive" : "rANS", trunc_ok);
        // (e) crafted geometry: huge lengths that would wrap a u64 addition
        int geo = 0;
        for (const int off : {56, 72, 88, 104}) {
            auto b = good;
            for (int i = 0; i < 8; ++i) b[static_cast<std::size_t>(off + i)] = 0xFF;
            recompute_payload_sha(b);
            try { decode_all(b); } catch (const std::exception&) { ++geo; }
        }
        CHECK(geo == 4, "%s crafted geometry rejected: %d of 4", naive ? "naive" : "rANS", geo);
        std::printf("corruption %-5s: %d rejected, %d benign (header padding-free bits that decode identically), 0 silent wrong sets\n",
                    naive ? "naive" : "rANS", caught, benign);
    }
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
    // --write DIR: dump the four fixture containers (d_min / full x rANS / byte-aligned) and
    // the fixture n-set, for independent reference decoders; no tests run
    if (argc > 3 && std::string(argv[2]) == "--write") {
        const Fixture fx = load_fixture(dir);
        const std::string out = argv[3];
        for (const bool full : {false, true})
            for (const bool naive : {false, true}) {
                const auto b = encode(full ? fx.full_paths : fx.dmin_paths, fx.ns, full, naive);
                std::ofstream o(out + "/fixture_" + (full ? "full" : "dmin") + (naive ? "_naive" : "_rans") + ".pwt2", std::ios::binary);
                o.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
            }
        for (const int cs : {2, 3})
            for (const bool full : {false, true}) {
                const auto b = encode(full ? fx.full_paths : fx.dmin_paths, fx.ns, full, false, 100000000ull, nullptr, cs);
                std::ofstream o(out + "/fixture_" + (full ? "full" : "dmin") + "_set" + std::to_string(cs) + ".pwt2", std::ios::binary);
                o.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
            }
        {   // the synthetic set whose mantissa models are populated (see test_ctx2), sets 2 and 3
            std::vector<u64> pr;
            for (u64 x = 1001; pr.size() < 30000; x += 2) if (is_prime_u64(x)) pr.push_back(x);
            std::vector<std::vector<u64>> paths;
            for (std::size_t j = 0; j + 2 < pr.size(); j += 3) {
                paths.push_back({3, pr[j]}); paths.push_back({5, pr[j + 1]}); paths.push_back({5, pr[j + 2]});
            }
            std::vector<u128> ns;
            for (const auto& q : paths) ns.push_back(static_cast<u128>(q[0]) * q[1]);
            std::sort(ns.begin(), ns.end());
            for (const int cs : {2, 3}) {
                const auto b = encode(paths, ns, true, false, 100000000ull, nullptr, cs);
                std::ofstream o(out + "/synthetic_mant_set" + std::to_string(cs) + ".pwt2", std::ios::binary);
                o.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
            }
            std::ofstream sp(out + "/synthetic_mant_paths.txt");
            for (const auto& q : paths) sp << q[0] << " " << q[1] << "\n";
            std::ofstream sn(out + "/synthetic_mant_n.txt");
            for (const u128 n : ns) sn << to_string(n) << "\n";
        }
        std::ofstream t(out + "/fixture_n.txt");
        for (const u128 n : fx.ns) t << to_string(n) << "\n";
        std::printf("wrote 8 fixture containers (sets 1-3), 2 synthetic mantissa containers, their paths and n-sets to %s\n", out.c_str());
        return 0;
    }
    test_bits();
    test_fastnum();
    test_bits_quot();
    test_universe();
    test_order();
    const Fixture fx = load_fixture(dir);
    std::printf("fixture: %zu Carmichael numbers\n", fx.ns.size());
    test_fixture(fx);
    test_synthetic();
    test_corruption(fx);
    test_crafted(dir);
    test_ctx2(fx);
    if (g_fail) { std::printf("test_pwt2: %d FAILURE(S)\n", g_fail); return 1; }
    std::printf("test_pwt2: all OK\n");
    return 0;
}
