// test_cnd.cpp -- tests for bic.hpp + cnd1.hpp (the CND1 container).
// Run from the repo root (fixtures in tests/fixtures/, temp files in build/).

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "bic.hpp"
#include "cnd1.hpp"
#include "lambda_bucket.hpp"
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

static u64 rng_state = 0xc0ffee123456789ull;
static u64 rnd64() {
    u64 x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return rng_state = x;
}

static void test_cmb() {
    for (u64 S = 1; S <= 300; ++S) {
        for (u64 v = 0; v < S; v += (S > 40 ? 7 : 1)) {
            std::vector<std::uint8_t> buf;
            {
                BitWriter bw(buf);
                bic_detail::cmb_put(bw, v, S);
                bw.align();
            }
            BitReader br(buf.data(), buf.size());
            CHECK(bic_detail::cmb_get(br, S) == v,
                  "cmb round-trip v=%llu S=%llu", (unsigned long long)v,
                  (unsigned long long)S);
        }
    }
    // wide spans
    for (int bits = 33; bits <= 90; bits += 19) {
        const u128 S = (static_cast<u128>(1) << bits) - 3;
        for (int i = 0; i < 50; ++i) {
            const u128 v = ((static_cast<u128>(rnd64()) << 64) | rnd64()) % S;
            std::vector<std::uint8_t> buf;
            {
                BitWriter bw(buf);
                bic_detail::cmb_put(bw, v, S);
                bw.align();
            }
            BitReader br(buf.data(), buf.size());
            CHECK(bic_detail::cmb_get(br, S) == v, "cmb wide %d bits", bits);
        }
    }
}

static std::vector<u128> random_sorted_set(std::size_t n, int val_bits) {
    std::vector<u128> v;
    v.reserve(n);
    const u128 mask = (static_cast<u128>(1) << val_bits) - 1;
    while (v.size() < n)
        v.push_back(((static_cast<u128>(rnd64()) << 64) | rnd64()) & mask);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

static void test_bic() {
    // sizes 0..3 and beyond, sparse and dense
    for (std::size_t n : {0u, 1u, 2u, 3u, 4u, 7u, 100u, 5000u}) {
        std::vector<u128> ks = random_sorted_set(n, 44);
        std::vector<std::uint8_t> blob;
        bic_append(blob, ks, 45);
        const std::vector<u128> back = bic_decode(blob.data(), blob.size(), ks.size(), 45);
        CHECK(back == ks, "bic round-trip n=%zu", ks.size());
    }
    // fully dense run: every span collapses to S=1, zero interior bits
    std::vector<u128> run;
    for (u128 i = 0; i < 1000; ++i) run.push_back(1000 + i);
    std::vector<std::uint8_t> blob;
    bic_append(blob, run, 12);
    CHECK(blob.size() == 3, "dense run should cost only the two ends (got %zu B)",
          blob.size());
    CHECK(bic_decode(blob.data(), blob.size(), run.size(), 12) == run, "dense run");
    // 60-bit values (full-table k scale)
    std::vector<u128> big = random_sorted_set(2000, 59);
    blob.clear();
    bic_append(blob, big, 60);
    CHECK(bic_decode(blob.data(), blob.size(), big.size(), 60) == big, "60-bit ks");
    // truncated blob must throw, not wander
    if (blob.size() > 4) {
        bool threw = false;
        try { bic_decode(blob.data(), blob.size() / 2, big.size(), 60); }
        catch (const std::exception&) { threw = true; }
        CHECK(threw, "truncated bic blob accepted");
    }
}

struct Fixture {
    std::vector<u128> ds, ns;                  // ds sorted; ns sorted
};

static Fixture load_fixture(const std::string& dir) {
    Fixture fx;
    std::vector<std::pair<u128, u128>> dn;
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
            dn.push_back({c.d, n});
        }
    }
    std::sort(dn.begin(), dn.end());
    for (const auto& [d, n] : dn) fx.ds.push_back(d);
    for (std::size_t i = 1; i < fx.ds.size(); ++i)
        CHECK(fx.ds[i] != fx.ds[i - 1], "duplicate fixture d");
    for (const auto& [d, n] : dn) fx.ns.push_back(n);
    std::sort(fx.ns.begin(), fx.ns.end());
    return fx;
}

static Cnd1FooterTargets targets_of(const std::vector<u128>& ns_sorted) {
    Cnd1FooterTargets t;
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
    return t;
}

static std::vector<std::uint8_t> file_bytes(const std::string& p) {
    std::ifstream in(p, std::ios::binary | std::ios::ate);
    std::vector<std::uint8_t> v(static_cast<std::size_t>(in.tellg()));
    in.seekg(0);
    in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size()));
    return v;
}

static void write_bytes(const std::string& p, const std::vector<std::uint8_t>& v) {
    std::ofstream out(p, std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size()));
}

static bool reader_throws(const std::string& p) {
    try { Cnd1Reader rd(p); (void)rd.decode_all(); return false; }
    catch (const std::exception&) { return true; }
}

static void test_container(const std::string& dir) {
    const Fixture fx = load_fixture(dir);
    const std::string tmp = "build/tmp_test_cnd.cnd1";
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        cnd1_encode(fx.ds, kCnd1M, targets_of(fx.ns), out);
    }
    Cnd1Reader rd(tmp);
    CHECK(rd.records() == fx.ds.size(), "record count");
    CHECK(rd.M() == kCnd1M, "M");
    const std::vector<u128> back = rd.decode_all();
    CHECK(back == fx.ds, "container d round-trip");

    // full decode-side reconstruction: factor every d, rebuild every n
    std::vector<u128> ns;
    for (const u128 d : back) {
        const SplitResult f = split_u128(d);
        CHECK(f.uncertified == 1, "uncertified split of d");
        u128 lam = 1;
        for (const u128 p : f.primes) lam = lcm128(lam, p - 1);
        CHECK(gcd128(d, lam) == 1, "gcd(d,lambda) != 1");
        ns.push_back(d * inv_mod128(d, lam));
    }
    std::sort(ns.begin(), ns.end());
    CHECK(ns == fx.ns, "reconstructed n-set != fixture n-set");
    const Cnd1FooterTargets t2 = targets_of(ns);
    CHECK(t2.sha_nset == rd.sha_nset(), "recomputed sha_nset mismatch");
    CHECK(t2.total_check == rd.total_check(), "recomputed total_check mismatch");

    // ---- corruption battery ----
    const std::vector<std::uint8_t> good = file_bytes(tmp);
    const std::string tmp2 = "build/tmp_test_cnd_bad.cnd1";

    auto mutate = [&](auto fn, const char* what) {
        std::vector<std::uint8_t> bad = good;
        fn(bad);
        write_bytes(tmp2, bad);
        CHECK(reader_throws(tmp2), "corruption undetected: %s", what);
    };
    mutate([](std::vector<std::uint8_t>& b) { b[0] = 'X'; }, "bad magic");
    mutate([](std::vector<std::uint8_t>& b) { b.resize(b.size() / 2); }, "truncation");
    mutate([](std::vector<std::uint8_t>& b) { b[b.size() / 2] ^= 0x40; },
           "body bit flip (payload sha)");
    mutate([](std::vector<std::uint8_t>& b) { b[b.size() - 4] ^= 1; },
           "footer count tamper");
    // directory tamper with a RECOMPUTED payload sha: the structural
    // validators (count/body-size consistency), not the sha, must catch it.
    mutate([](std::vector<std::uint8_t>& b) {
        b[64] ^= 3;                            // first class count vbyte
        u64 footer_off = 0;
        for (int i = 7; i >= 0; --i) footer_off = (footer_off << 8) | b[40 + i];
        const auto sha = Sha256::hash(b.data(), static_cast<std::size_t>(footer_off));
        std::memcpy(&b[footer_off], sha.data(), 32);
    }, "directory tamper w/ recomputed sha");

    std::remove(tmp2.c_str());
    std::remove(tmp.c_str());

    // ---- crafted-file battery: files the
    // encoder can never emit, assembled byte-by-byte per the spec with a
    // CORRECT payload sha, so only the targeted validator can reject ----
    struct DirEnt { u64 count; bool has_bl; u64 bl; };
    auto craft = [&](u32 M, u32 end_width, u64 records,
                     const std::vector<std::pair<u32, DirEnt>>& ents,
                     const std::vector<std::uint8_t>& body) {
        std::vector<std::uint8_t> f;
        f.insert(f.end(), {'C', 'N', 'D', '1'});
        cnd1_detail::put32(f, kCnd1Version);
        cnd1_detail::put32(f, M);
        cnd1_detail::put32(f, end_width);
        cnd1_detail::put64(f, records);
        std::vector<std::uint8_t> dir;
        std::size_t ei = 0;
        for (u32 r = 0; r < M; ++r) {
            if (ei < ents.size() && ents[ei].first == r) {
                vbyte_append(dir, ents[ei].second.count);
                if (ents[ei].second.has_bl) vbyte_append(dir, ents[ei].second.bl);
                ++ei;
            } else {
                vbyte_append(dir, 0);
            }
        }
        const u64 body_off = 64 + dir.size();
        cnd1_detail::put64(f, 64);
        cnd1_detail::put64(f, body_off);
        cnd1_detail::put64(f, body_off + body.size());
        f.resize(64, 0);
        f.insert(f.end(), dir.begin(), dir.end());
        f.insert(f.end(), body.begin(), body.end());
        const auto sha = Sha256::hash(f.data(), f.size());
        f.insert(f.end(), sha.begin(), sha.end());
        for (int i = 0; i < 32; ++i) f.push_back(0);   // sha_nset (unreached)
        cnd1_detail::put128(f, 0);                     // total_check
        cnd1_detail::put64(f, records);
        write_bytes(tmp2, f);
        return reader_throws(tmp2);
    };
    // wraparound: two blob lengths whose u64 sum lands exactly on footer_off
    CHECK(craft(kCnd1M, 5, 2,
                {{1, {1, true, 1ull << 60}},
                 {3, {1, true, (0ull - (1ull << 60)) + 2}}},
                {0x08, 0x08}),
          "u64 directory wraparound accepted");
    // blob longer than its bits: trailing garbage byte inside the blob
    CHECK(craft(kCnd1M, 5, 1, {{1, {1, true, 2}}}, {0x08, 0x00}),
          "over-long blob (trailing byte) accepted");
    // nonzero padding bits inside an otherwise exact blob
    CHECK(craft(kCnd1M, 5, 1, {{1, {1, true, 1}}}, {0x09}),
          "nonzero padding accepted");
    // wrong modulus
    CHECK(craft(210, 5, 1, {{1, {1, true, 1}}}, {0x08}),
          "M != 2310 accepted");
    // count exceeding the k-universe of end_width
    CHECK(craft(kCnd1M, 1, 3, {{1, {3, true, 1}}}, {0x40}),
          "class count above 2^end_width accepted");
    // control: the same crafting path with honest fields must be READABLE
    CHECK(!craft(kCnd1M, 5, 1, {{1, {1, true, 1}}}, {0x08}),
          "control crafted file rejected (craft harness broken?)");

    // header tamper on a real file: end_width flip must now hit sha_payload
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        cnd1_encode(fx.ds, kCnd1M, targets_of(fx.ns), out);
    }
    std::vector<std::uint8_t> hbad = file_bytes(tmp);
    hbad[12] ^= 1;                             // end_width low byte
    write_bytes(tmp2, hbad);
    CHECK(reader_throws(tmp2), "header end_width tamper undetected");
    std::remove(tmp2.c_str());
    std::remove(tmp.c_str());

    // encoder domain: a d whose k needs > 64 bits must be refused
    {
        bool threw = false;
        std::ostringstream sink;
        try { cnd1_encode({(static_cast<u128>(1) << 100) + 1}, kCnd1M, {}, sink); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw, "encoder accepted k wider than 64 bits");
    }

    // empty container round-trips
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        cnd1_encode({}, kCnd1M, Cnd1FooterTargets{}, out);
    }
    Cnd1Reader empty(tmp);
    CHECK(empty.records() == 0 && empty.decode_all().empty(), "empty container");
    std::remove(tmp.c_str());
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
    test_cmb();
    test_bic();
    test_container(dir);
    if (failures) {
        std::printf("%d cnd test failure(s)\n", failures);
        return 1;
    }
    std::printf("all cnd tests passed\n");
    return 0;
}
