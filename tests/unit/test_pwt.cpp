// test_pwt -- PWT1 container tests (docs/PWT_FORMAT.md). Run from the repo
// root (fixtures in tests/fixtures/).
//   1. prime bitmap: rank/select/advance against a naive sieve
//   2. rANS: multi-context block round trip, skewed and dense alphabets
//   3. fixture round trips, both modes (d_min paths and full factorizations):
//      encode -> decode -> sorted n-set == fixture n-set, footer targets
//   4. synthetic tree with a deep escape (count >= 4095 at the root), a
//      small sieve limit so both universes and the depth-low flags are
//      exercised, and a terminal-internal node
//   5. corruption battery: every kind of damage must be rejected or fail
//      certification, never decode silently to a wrong set

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "lambda_bucket.hpp"
#include "prime_bitmap.hpp"
#include "pwt1.hpp"
#include "pwt_paths.hpp"
#include "rans.hpp"
#include "u128.hpp"

using namespace cn;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); ++g_fail; } } while (0)

namespace {

struct Fixture {
    std::vector<std::vector<u64>> dmin_paths, full_paths;
    std::vector<u128> ns;   // sorted
};

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
            fx.dmin_paths.push_back(dp);
            fx.full_paths.push_back(fp);
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

// Decode a container and return its sorted n-set (throws on structural
// damage only). `targets_ok` reports whether the footer targets match the
// decoded set, which is what pwt_decode checks before certifying.
std::vector<u128> decode_raw(std::vector<std::uint8_t> bytes, bool* targets_ok = nullptr) {
    Pwt1Reader rd(std::move(bytes));
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
    if (targets_ok) *targets_ok = sha.finish() == rd.targets().sha_nset && sum == rd.targets().total_check;
    return ns;
}
// The same, throwing on a target mismatch (pwt_decode's behaviour).
std::vector<u128> decode_all(std::vector<std::uint8_t> bytes) {
    bool ok = false;
    std::vector<u128> ns = decode_raw(std::move(bytes), &ok);
    if (!ok) throw std::runtime_error("footer targets mismatch");
    return ns;
}

// Write paths in trie_paths' on-disk format (K five-byte big-endian
// ASCENDING fields, zero padded) so PathRunSource can be exercised.
template <int K>
void write_path_file(const std::string& path, const std::vector<std::vector<u64>>& paths) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    for (auto p : paths) {
        std::sort(p.begin(), p.end());
        PathRec<K> r{};
        for (std::size_t i = 0; i < p.size(); ++i) path_set_field<K>(r, static_cast<int>(i), p[i]);
        out.write(reinterpret_cast<const char*>(r.b), sizeof r.b);
    }
}

void test_bitmap() {
    const u64 L = 200003;
    PrimeBitmap bm(L);
    std::vector<char> comp(L, 0);
    std::vector<u64> primes;
    for (u64 i = 2; i < L; ++i) {
        if (comp[i]) continue;
        if (i > 2) primes.push_back(i);
        for (u64 j = i * i; j < L; j += i) comp[j] = 1;
    }
    CHECK(bm.count() == primes.size(), "bitmap count %llu vs %zu", (unsigned long long)bm.count(), primes.size());
    for (std::size_t i = 0; i < primes.size(); ++i) {
        CHECK(bm.rank(primes[i]) == i, "rank(%llu)", (unsigned long long)primes[i]);
        CHECK(bm.select(i) == primes[i], "select(%zu)", i);
    }
    std::mt19937_64 rng(7);
    for (int t = 0; t < 20000; ++t) {
        const std::size_t i = rng() % (primes.size() - 1);
        const u64 g = 1 + rng() % std::min<std::size_t>(5000, primes.size() - 1 - i);
        CHECK(bm.advance(primes[i], g, i + g) == primes[i + g], "advance(%zu, %llu)", i, (unsigned long long)g);
    }
    CHECK(bm.is_prime(199999) && !bm.is_prime(200001), "is_prime");
    std::printf("bitmap: %zu primes below %llu, rank/select/advance OK\n", primes.size(), (unsigned long long)L);
}

void test_rans() {
    std::mt19937_64 rng(11);
    std::vector<RansContext> ctx(3);
    std::vector<u64> h0(4096, 0), h1(2, 0), h2(128, 0);
    std::vector<std::pair<u16, u16>> pairs;
    for (int i = 0; i < 300000; ++i) {
        const u16 s0 = static_cast<u16>((rng() % 100 == 0) ? rng() % 4096 : 1);   // skewed, wide
        const u16 s1 = static_cast<u16>(rng() % 1000 == 0);                       // near-deterministic
        const u16 s2 = static_cast<u16>(1 + (rng() % 40));                        // dense
        pairs.emplace_back(0, s0); ++h0[s0];
        pairs.emplace_back(1, s1); ++h1[s1];
        pairs.emplace_back(2, s2); ++h2[s2];
    }
    ctx[0].build(h0.data(), 4096); ctx[1].build(h1.data(), 2); ctx[2].build(h2.data(), 128);
    std::vector<std::uint8_t> tab;
    for (auto& c : ctx) c.serialize(tab);
    std::vector<RansContext> back(3);
    std::size_t pos = 0;
    pos += back[0].deserialize(tab.data() + pos, tab.size() - pos, 4096);
    pos += back[1].deserialize(tab.data() + pos, tab.size() - pos, 2);
    pos += back[2].deserialize(tab.data() + pos, tab.size() - pos, 128);
    CHECK(pos == tab.size(), "table round trip size");
    for (int c = 0; c < 3; ++c) CHECK(back[c].freq == ctx[c].freq && back[c].syms == ctx[c].syms, "table %d", c);
    std::vector<std::uint8_t> enc;
    rans_encode_block(pairs, ctx, enc);
    RansBlockDecoder dec(enc.data(), enc.size());
    bool ok = true;
    for (const auto& [c, s] : pairs) if (dec.decode(back[c]) != s) { ok = false; break; }
    CHECK(ok && dec.finished(), "rANS block round trip");
    double ideal = 0.0;
    for (int c = 0; c < 3; ++c) {
        const std::vector<u64>& h = c == 0 ? h0 : c == 1 ? h1 : h2;
        for (std::size_t i = 0; i < ctx[c].syms.size(); ++i)
            ideal += double(h[ctx[c].syms[i]]) * (16.0 - std::log2(double(ctx[c].freq[i])));
    }
    std::printf("rANS: %zu symbols -> %zu bytes (ideal %.0f bytes, +%.3f%%)\n", pairs.size(), enc.size(),
                ideal / 8, 100.0 * (enc.size() * 8.0 - ideal) / ideal);
    CHECK(enc.size() * 8.0 < ideal * 1.01 + 64, "rANS within 1%% of the model");
    // a symbol absent from a NON-EMPTY model must be refused (the block is
    // coded in reverse, so put the offending pair last)
    {
        std::vector<u64> h1only(4096, 0);
        h1only[1] = 100;
        std::vector<RansContext> one(1);
        one[0].build(h1only.data(), 4096);
        bool threw = false;
        try { rans_encode_block({{0, 1}, {0, 2}}, one, enc); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "absent symbol refused");
        threw = false;
        try { rans_encode_block({{0, 1}}, std::vector<RansContext>{RansContext()}, enc); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "empty context refused");
    }
    // malformed tables: zero frequency via the freq-1 wrap, sum != scale, non-ascending symbols
    {
        auto vb = [](std::vector<std::uint8_t>& o, u64 v) { while (v >= 0x80) { o.push_back(static_cast<std::uint8_t>(v | 0x80)); v >>= 7; } o.push_back(static_cast<std::uint8_t>(v)); };
        std::vector<std::vector<std::uint8_t>> bad(3);
        vb(bad[0], 2); vb(bad[0], 3); vb(bad[0], 65535); vb(bad[0], 1); vb(bad[0], ~0ull);   // freq 0 by wrap
        vb(bad[1], 1); vb(bad[1], 3); vb(bad[1], 65534);                                       // sum 65535
        vb(bad[2], 2); vb(bad[2], 3); vb(bad[2], 32767); vb(bad[2], 0); vb(bad[2], 32767);     // delta 0 = same symbol
        for (int k = 0; k < 3; ++k) {
            RansContext c;
            bool threw = false;
            try { c.deserialize(bad[k].data(), bad[k].size(), 4096); } catch (const std::exception&) { threw = true; }
            CHECK(threw, "malformed table %d rejected", k);
        }
    }
}

std::vector<std::uint8_t> roundtrip(const char* label, const std::vector<std::vector<u64>>& paths, bool full,
                                    u64 B, const std::vector<u128>& ns_sorted, u64 expect_escapes = 0,
                                    bool expect_both_universes = false) {
    PathVectorSource src(paths);
    PwtEncoder enc(full, B, false);
    PwtEncodeStats st;
    const std::vector<std::uint8_t> bytes = enc.encode(src, targets_of(ns_sorted), &st);
    CHECK(st.escapes == expect_escapes, "%s: escapes %llu, expected %llu", label,
          (unsigned long long)st.escapes, (unsigned long long)expect_escapes);
    if (expect_both_universes) CHECK(st.low_edges > 0 && st.high_edges > 0, "%s: both universes used", label);
    std::vector<u128> got;
    {
        Pwt1Reader rd(bytes);
        CHECK(rd.header().full == full && rd.header().n_records == paths.size(), "%s header", label);
        rd.for_each_n([&](u128 n) { got.push_back(n); });
    }
    std::sort(got.begin(), got.end());
    CHECK(got == ns_sorted, "%s: decoded n-set != fixture n-set", label);
    std::printf("%s: %zu paths, %llu nodes, %zu bytes = %.3f bits/el (symbols %.3f vs model %.3f, raw %.3f, tables %llu B), edges low %llu high %llu\n",
                label, paths.size(), (unsigned long long)st.nodes, bytes.size(), 8.0 * bytes.size() / paths.size(),
                8.0 * st.sym_bytes / paths.size(), st.ideal_symbol_bits / paths.size(), double(st.raw_bits) / paths.size(),
                (unsigned long long)st.tables_bytes, (unsigned long long)st.low_edges, (unsigned long long)st.high_edges);
    return bytes;
}

void test_synthetic() {
    // Paths chosen so the tree has: a root with >= 4095 children (escape),
    // chains, a terminal-internal node (a path that is a prefix of another),
    // and primes on both sides of a small sieve limit B = 2000 so depth-low
    // flags and both universes occur. Values are products (full mode), which
    // only requires distinct paths; Korselt structure is not needed here.
    PrimeBitmap bm(1u << 20);
    std::vector<u64> primes;
    for (u64 i = 0; i < bm.count() && primes.size() < 6000; ++i) primes.push_back(bm.select(i));
    std::vector<std::vector<u64>> paths;
    std::mt19937_64 rng(5);
    for (std::size_t i = 1; i < 4200; ++i) paths.push_back({primes[i]});                 // 4199 root children
    for (int t = 0; t < 3000; ++t) {                                                      // deeper paths
        std::vector<u64> p;
        std::size_t top = 400 + rng() % 5000;                                             // depth 1: mostly > B
        p.push_back(primes[top]);
        int len = 2 + static_cast<int>(rng() % 6);
        std::size_t cur = top;
        for (int j = 1; j < len && cur > 1; ++j) {
            cur = rng() % std::min<std::size_t>(cur, 200);                                // depth >= 2: < primes[200] < B
            if (cur == 0 || p.back() % primes[cur] == 1) break;   // keep the Korselt invariant q != 1 (mod p)
            p.push_back(primes[cur]);
        }
        paths.push_back(p);
    }
    std::size_t a = 150, b = 7;                                          // terminal-internal pair obeying the invariant
    while (primes[5000] % primes[a] == 1) ++a;
    while (primes[a] % primes[b] == 1) ++b;
    paths.push_back({primes[5000], primes[a]});                         // terminal-internal ...
    paths.push_back({primes[5000], primes[a], primes[b]});              // ... with a child
    std::sort(paths.begin(), paths.end());
    paths.erase(std::unique(paths.begin(), paths.end()), paths.end());
    std::vector<u128> ns;
    for (const auto& p : paths) { u128 v = 1; for (u64 q : p) v *= q; ns.push_back(v); }
    std::sort(ns.begin(), ns.end());
    CHECK(std::adjacent_find(ns.begin(), ns.end()) == ns.end(), "synthetic products distinct");
    const auto bytes = roundtrip("synthetic (B = 2000, escape, both universes)", paths, true, 2000, ns, 1, true);
    Pwt1Reader rd(bytes);
    CHECK(rd.header().depth_low != 0, "some depth flagged low");
    CHECK(!rd.header().depth_is_low(1), "depth 1 not low");
}

void test_corruption(const std::vector<std::uint8_t>& good, const std::vector<u128>& ns) {
    std::mt19937_64 rng(3);
    int rejected = 0, wrong_set = 0, silent = 0;
    const int trials = 400;
    for (int t = 0; t < trials; ++t) {
        std::vector<std::uint8_t> bad = good;
        const std::size_t i = rng() % bad.size();
        bad[i] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
        try {
            const std::vector<u128> got = decode_all(bad);
            if (got == ns) ++silent; else ++wrong_set;
        } catch (const std::exception&) { ++rejected; }
    }
    // every single-bit flip lands under the payload hash, so all must be rejected
    CHECK(rejected == trials, "bit flips: %d rejected, %d wrong set, %d silent", rejected, wrong_set, silent);
    // truncation and extension
    for (std::size_t cut : {std::size_t(1), std::size_t(17), good.size() / 2}) {
        bool threw = false;
        try { decode_all(std::vector<std::uint8_t>(good.begin(), good.end() - cut)); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "truncation by %zu rejected", cut);
    }
    {
        std::vector<std::uint8_t> ext = good;
        ext.push_back(0);
        bool threw = false;
        try { decode_all(ext); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "extension rejected");
    }
    // damage BELOW the hash: rebuild the payload hash after flipping a symbol
    // byte, so only the walk's own checks can catch it
    int caught = 0, wrong = 0, quiet = 0;
    for (int t = 0; t < 200; ++t) {
        std::vector<std::uint8_t> bad = good;
        const u64 sym_off = pwt_detail::get_le(bad.data() + 48, 8), sym_len = pwt_detail::get_le(bad.data() + 56, 8);
        const u64 raw_off = pwt_detail::get_le(bad.data() + 64, 8), raw_len = pwt_detail::get_le(bad.data() + 72, 8);
        const u64 footer_off = pwt_detail::get_le(bad.data() + 104, 8);
        const bool in_raw = (t % 2) && raw_len >= 16;   // PWT0 has no raw section
        const std::size_t i = in_raw ? raw_off + 8 + rng() % (raw_len - 8) : sym_off + 8 + rng() % (sym_len - 8);
        bad[i] ^= static_cast<std::uint8_t>(1u << (rng() % 8));
        Sha256 sha;
        sha.update(bad.data(), footer_off);
        const auto h = sha.finish();
        std::copy(h.begin(), h.end(), bad.begin() + footer_off);
        try {
            bool targets_ok = false;
            const std::vector<u128> got = decode_raw(bad, &targets_ok);
            if (got == ns) ++quiet;
            else { ++wrong; CHECK(!targets_ok, "wrong set passed the footer targets"); }
        } catch (const std::exception&) { ++caught; }
    }
    std::printf("corruption below the hash: %d caught structurally, %d decoded to a wrong set (footer targets catch), %d unnoticed\n",
                caught, wrong, quiet);
    CHECK(quiet == 0, "a damaged payload decoded to the right set");
}

// Crafted headers whose section lengths wrap a u64 so that naive
// "off + len == next" checks still hold: must be rejected, never read.
void test_crafted_geometry(const std::vector<std::uint8_t>& good) {
    auto put = [](std::vector<std::uint8_t>& b, std::size_t off, u64 v) { for (int i = 0; i < 8; ++i) b[off + i] = static_cast<std::uint8_t>(v >> (8 * i)); };
    auto get = [](const std::vector<std::uint8_t>& b, std::size_t off) { return pwt_detail::get_le(b.data() + off, 8); };
    const u64 sym_off = get(good, 48), footer_off = get(good, 104);
    struct Craft { const char* name; u64 tables_len, sym_len; };
    const Craft crafts[] = {
        {"sym_len = 2^63 with wrapped raw section", 0, u64(1) << 63},
        {"tables_len = 2^64 - 100", u64(0) - 100, 0},
        {"sym_len just past the footer", 0, u64(1) << 40},
    };
    for (const Craft& c : crafts) {
        std::vector<std::uint8_t> bad = good;
        const u64 tables_len = c.tables_len ? c.tables_len : get(good, 40);
        const u64 s_off = 128 + tables_len;
        const u64 s_len = c.sym_len ? c.sym_len : get(good, 56);
        const u64 r_off = s_off + s_len;                 // wraps on purpose
        put(bad, 40, tables_len); put(bad, 48, s_off); put(bad, 56, s_len);
        put(bad, 64, r_off); put(bad, 72, footer_off - r_off);
        if (c.sym_len) { put(bad, sym_off, 1); put(bad, sym_off + 4, 0x7fffffffull); }   // a block header pointing far away
        Sha256 sha;
        sha.update(bad.data(), footer_off);
        const auto h = sha.finish();
        std::copy(h.begin(), h.end(), bad.begin() + footer_off);
        bool threw = false;
        try { decode_raw(bad); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "crafted geometry rejected: %s", c.name);
    }
    std::printf("crafted geometry: wrapped section lengths are rejected before any section is read\n");
}

// The production input path: trie_paths-format file -> external merge sort
// with tiny chunks (many run files, two replays) -> the SAME container as
// the in-memory source, byte for byte.
void test_disk_source(const Fixture& fx) {
    const std::string dir = "build";
    write_path_file<12>(dir + "/test_pwt_dmin.bin", fx.dmin_paths);
    write_path_file<14>(dir + "/test_pwt_full.bin", fx.full_paths);
    for (int mode = 0; mode < 2; ++mode) {
        const bool full = mode == 1;
        PathSortOptions opt;
        opt.scratch = dir; opt.tag = full ? "test_pwt_f" : "test_pwt_d"; opt.chunk_recs = 50; opt.verbose = false;
        std::vector<std::uint8_t> disk_bytes, mem_bytes;
        PwtEncodeStats st;
        if (full) {
            PathRunSource<14> src(dir + "/test_pwt_full.bin", opt);
            CHECK(src.size() == fx.full_paths.size() && src.max_len() == 11, "disk source full: size/max_len");
            disk_bytes = PwtEncoder(true, 100000000ull, false).encode(src, targets_of(fx.ns), &st);
            mem_bytes = PwtEncoder(true, 100000000ull, false).encode(PathVectorSource(fx.full_paths), targets_of(fx.ns));
        } else {
            PathRunSource<12> src(dir + "/test_pwt_dmin.bin", opt);
            CHECK(src.size() == fx.dmin_paths.size() && src.max_len() == 8, "disk source d_min: size/max_len");
            CHECK(src.max_prime_at_depth(1) > 100000000ull && src.max_prime_at_depth(2) < 100000000ull, "max prime per depth");
            disk_bytes = PwtEncoder(false, 100000000ull, false).encode(src, targets_of(fx.ns), &st);
            mem_bytes = PwtEncoder(false, 100000000ull, false).encode(PathVectorSource(fx.dmin_paths), targets_of(fx.ns));
        }
        CHECK(disk_bytes == mem_bytes, "%s: disk-sorted container == vector-sorted container", full ? "full" : "d_min");
        CHECK(decode_all(disk_bytes) == fx.ns, "%s: disk container decodes to the fixture set", full ? "full" : "d_min");
    }
    // a paths file that does not match the oracle targets must be refused by the encoder
    {
        std::vector<u128> wrong = fx.ns;
        wrong[0] += 2;
        bool threw = false;
        try { PwtEncoder(false, 100000000ull, false).encode(PathVectorSource(fx.dmin_paths), targets_of(wrong)); }
        catch (const std::exception&) { threw = true; }
        CHECK(threw, "encoder refuses paths that do not reproduce the oracle's sum");
    }
    std::remove((dir + "/test_pwt_dmin.bin").c_str());
    std::remove((dir + "/test_pwt_full.bin").c_str());
    std::printf("disk source: external sort with 50-record runs reproduces the in-memory container byte for byte\n");
}

// PWT0 (byte-aligned): the container's symbol section must be exactly the
// trie_stats NAIVE byte count of the fixture (d_min 4936, full 8574 B), the
// file must decode to the fixture set, and damage must be rejected.
void test_naive(const Fixture& fx) {
    for (int mode = 0; mode < 2; ++mode) {
        const bool full = mode == 1;
        const auto& paths = full ? fx.full_paths : fx.dmin_paths;
        PwtEncodeStats st;
        const std::vector<std::uint8_t> bytes =
            PwtEncoder(full, 100000000ull, false, true).encode(PathVectorSource(paths), targets_of(fx.ns), &st);
        const u64 expect = full ? 8574 : 4936;
        CHECK(st.sym_bytes == expect, "%s PWT0 stream %llu B, expected %llu (trie_stats NAIVE)",
              full ? "full" : "d_min", (unsigned long long)st.sym_bytes, (unsigned long long)expect);
        CHECK(st.tables_bytes == 0 && st.raw_bytes == 0, "PWT0 has no tables or raw bits");
        CHECK(bytes.size() == st.sym_bytes + kPwt1HeaderBytes + kPwt1FooterBytes, "PWT0 file = stream + framing");
        Pwt1Reader rd(bytes);
        CHECK(rd.header().naive && rd.header().full == full, "PWT0 header");
        CHECK(decode_all(bytes) == fx.ns, "%s PWT0 decodes to the fixture set", full ? "full" : "d_min");
        std::printf("PWT0 %s: %zu bytes = %.3f bits/el (stream %llu B)\n", full ? "full" : "d_min", bytes.size(),
                    8.0 * bytes.size() / fx.ns.size(), (unsigned long long)st.sym_bytes);
        if (!full) test_corruption(bytes, fx.ns);
        // a PWT1 magic on a PWT0 body (and vice versa) must be refused
        std::vector<std::uint8_t> bad = bytes;
        bad[3] = '1';
        bool threw = false;
        try { decode_raw(bad); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "magic/flag disagreement rejected");
    }
}

// A root with more than 65535 children exercises the u16 count sentinel
// and the overflow map (the full-table root has 6.9 M children).
void test_wide_root() {
    PrimeBitmap bm(1u << 21);
    std::vector<std::vector<u64>> paths;
    for (u64 i = 0; i < 70000; ++i) paths.push_back({bm.select(i)});
    std::vector<u128> ns;
    for (const auto& p : paths) ns.push_back(p[0]);
    std::sort(ns.begin(), ns.end());
    roundtrip("wide root (70,000 children)", paths, true, 100000000ull, ns, 1, false);
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
    test_bitmap();
    test_rans();
    const Fixture fx = load_fixture(dir);
    std::printf("fixture: %zu lines\n", fx.ns.size());
    const auto d_bytes = roundtrip("d_min paths, B = 1e8", fx.dmin_paths, false, 100000000ull, fx.ns);
    const auto f_bytes = roundtrip("full paths, B = 1e8", fx.full_paths, true, 100000000ull, fx.ns);
    roundtrip("d_min paths, B = 1e5 (mixed universes)", fx.dmin_paths, false, 100000ull, fx.ns, 0, true);
    roundtrip("full paths, B = 1e4 (mixed universes)", fx.full_paths, true, 10000ull, fx.ns, 0, true);
    test_synthetic();
    test_wide_root();
    test_disk_source(fx);
    test_naive(fx);
    test_corruption(d_bytes, fx.ns);
    test_corruption(f_bytes, fx.ns);
    test_crafted_geometry(d_bytes);
    // a wrong-mode read must not certify: full container decoded as d_min is impossible
    // (mode is in the header) -- instead check the footer targets are the fixture's
    {
        Pwt1Reader rd(d_bytes);
        CHECK(rd.targets().record_count == fx.ns.size(), "footer count");
        CHECK(rd.targets().total_check == targets_of(fx.ns).total_check, "footer sum");
    }
    if (g_fail) { std::printf("test_pwt: %d FAILURE(S)\n", g_fail); return 1; }
    std::printf("test_pwt: all OK\n");
    return 0;
}
