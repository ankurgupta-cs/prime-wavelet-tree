// test_oracle -- unit tests for vbyte.hpp and orc1.hpp (container round-trip,
// index/rank bookkeeping, input validation). Integration against real table
// samples is exercised separately via oracle_encode/oracle_decode.
//
// Writes a scratch file `tmp_orc1_test.bin` in the current directory.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <vector>

#include "orc1.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

using namespace cn;

static int failures = 0;
#define CHECK(cond)                                                        \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);    \
            ++failures;                                                    \
        }                                                                  \
    } while (0)

static void test_vbyte() {
    std::vector<std::uint8_t> buf;

    // Round-trip small dense range and structured boundary values.
    std::vector<u128> vals;
    for (u64 v = 0; v <= 100'000; ++v) vals.push_back(v);
    for (int s = 0; s < 128; ++s) {
        const u128 p = static_cast<u128>(1) << s;
        vals.push_back(p - 1);
        vals.push_back(p);
        vals.push_back(p + 1);
    }
    vals.push_back(UINT64_MAX);
    vals.push_back(U128_MAX);
    for (u128 v : vals) vbyte_append(buf, v);

    const std::uint8_t* p = buf.data();
    const std::uint8_t* end = buf.data() + buf.size();
    for (u128 v : vals) {
        u128 got;
        p = vbyte_decode(p, end, got);
        CHECK(p != nullptr && got == v);
        if (p == nullptr) return;
    }
    CHECK(p == end);

    // Length function agrees with the encoder.
    for (u128 v : {static_cast<u128>(0), static_cast<u128>(127),
                   static_cast<u128>(128), static_cast<u128>(16383),
                   static_cast<u128>(16384), static_cast<u128>(UINT64_MAX),
                   U128_MAX}) {
        std::vector<std::uint8_t> one;
        vbyte_append(one, v);
        CHECK(static_cast<int>(one.size()) == vbyte_length(v));
    }

    // Truncated input is rejected.
    {
        std::vector<std::uint8_t> t;
        vbyte_append(t, static_cast<u128>(1) << 40);
        u128 got;
        CHECK(vbyte_decode(t.data(), t.data() + t.size() - 1, got) == nullptr);
    }
    // u64 overflow is rejected by the narrow decoder.
    {
        std::vector<std::uint8_t> t;
        vbyte_append(t, static_cast<u128>(UINT64_MAX) + 1);
        u64 got;
        CHECK(vbyte_decode_u64(t.data(), t.data() + t.size(), got) == nullptr);
    }
    // A 19-byte encoding whose payload exceeds 128 bits is rejected, not
    // silently reduced mod 2^128 (18 continuation bytes then 0x04 encodes
    // 4 << 126 = 2^128).
    {
        std::vector<std::uint8_t> t;
        t.reserve(19);                 // single allocation placates -Warray-bounds
        t.assign(18, 0x80);
        t.push_back(0x04);
        u128 got;
        CHECK(vbyte_decode(t.data(), t.data() + t.size(), got) == nullptr);
        // The all-ones maximum still round-trips (last byte 0x03).
        std::vector<std::uint8_t> m;
        vbyte_append(m, U128_MAX);
        CHECK(m.size() == 19 && m.back() == 0x03);
        CHECK(vbyte_decode(m.data(), m.data() + m.size(), got) != nullptr &&
              got == U128_MAX);
    }
}

struct Rec {
    u128 n;
    std::vector<u64> f;
};

static std::vector<Rec> synthetic_records(int count) {
    // Container-level validity only (ascending factors, n == product,
    // ascending n): primality is the encoder's concern, not the container's.
    std::vector<Rec> recs;
    for (int i = 0; i < count; ++i) {
        Rec r;
        r.f = {3, 5, static_cast<u64>(7 + 2 * i)};
        if (i % 5 == 0) r.f.push_back(static_cast<u64>(1'000'000'007 + 2 * i));
        r.n = 1;
        for (u64 x : r.f) r.n *= x;
        recs.push_back(r);
    }
    // Products with the extra factor break global monotonicity of n as
    // generated; sort by n and drop duplicates to restore the contract.
    std::sort(recs.begin(), recs.end(),
              [](const Rec& a, const Rec& b) { return a.n < b.n; });
    recs.erase(std::unique(recs.begin(), recs.end(),
                           [](const Rec& a, const Rec& b) { return a.n == b.n; }),
               recs.end());
    return recs;
}

static void test_orc1() {
    const char* path = "tmp_orc1_test.bin";
    const std::vector<Rec> recs = synthetic_records(1000);

    std::array<std::uint8_t, 32> fake_sha{};
    fake_sha[0] = 0xAB;
    {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        Orc1Writer w(out, /*block_size=*/8);
        for (const Rec& r : recs) w.add(r.n, r.f.data(), static_cast<int>(r.f.size()));
        const Orc1Footer f = w.finish(fake_sha);
        CHECK(w.records_written() == recs.size());
        CHECK(f.n_blocks == (recs.size() + 7) / 8);
    }

    // Scoped so the reader's file handle closes before std::remove below —
    // Windows cannot delete a file that is still open.
    {
    Orc1Reader rd(path);
    CHECK(rd.header().n_records == recs.size());
    CHECK(rd.header().block_size == 8);
    CHECK(rd.footer().sha_src == fake_sha);

    // Full-order decode matches input exactly.
    {
        std::size_t i = 0;
        rd.for_all([&](const Orc1Record& r) {
            CHECK(i < recs.size());
            if (i >= recs.size()) return;
            CHECK(r.n == recs[i].n);
            CHECK(static_cast<std::size_t>(r.k) == recs[i].f.size());
            for (int j = 0; j < r.k; ++j) CHECK(r.factors[j] == recs[i].f[j]);
            ++i;
        });
        CHECK(i == recs.size());
    }

    // total_check.
    {
        u128 sum = 0;
        for (const Rec& r : recs) sum += r.n;
        CHECK(rd.footer().total_check == sum);
    }

    // find_block + in-block scan finds every record, and rejects absentees.
    for (std::size_t i = 0; i < recs.size(); i += 97) {
        const u64 b = rd.find_block(recs[i].n);
        CHECK(b < rd.index().size());
        bool found = false;
        rd.for_block(b, [&](const Orc1Record& r) {
            if (r.n == recs[i].n) found = true;
        });
        CHECK(found);
    }
    CHECK(rd.find_block(recs.front().n - 1) == rd.index().size());
    {
        const u64 b = rd.find_block(recs.back().n + 1);   // beyond last: lands
        CHECK(b == rd.index().size() - 1);                // in final block
    }

    // rank_base bookkeeping.
    for (std::size_t b = 0; b < rd.index().size(); ++b)
        CHECK(rd.index()[b].rank_base == static_cast<u64>(b) * 8);
    }

    // Writer input validation.
    {
        std::ofstream out("tmp_orc1_bad.bin", std::ios::binary | std::ios::trunc);
        Orc1Writer w(out);
        const u64 asc[3] = {3, 5, 7};
        bool threw = false;
        try { w.add(3 * 5 * 7 + 1, asc, 3); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);   // n != product
        const u64 dup[3] = {3, 5, 5};
        threw = false;
        try { w.add(75, dup, 3); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);   // not strictly ascending
        threw = false;
        const u64 two[2] = {3, 5};
        try { w.add(15, two, 2); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);   // k < 3
        w.add(105, asc, 3);
        threw = false;
        try { w.add(105, asc, 3); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);   // n not ascending
        w.finish(fake_sha);
    }

    // Bad magic rejected.
    {
        std::ofstream out("tmp_orc1_bad.bin", std::ios::binary | std::ios::trunc);
        out << "not an orc1 file at all, padded well past sixty-four bytes....";
        out << "................................";
        out.close();
        bool threw = false;
        try { Orc1Reader r2("tmp_orc1_bad.bin"); } catch (const std::exception&) { threw = true; }
        CHECK(threw);
    }

    // Reader hardening: byte-patched containers are rejected.
    {
        auto patched = [&](auto&& mutate) {
            std::ifstream src("tmp_orc1_test.bin", std::ios::binary);
            std::vector<char> bytes((std::istreambuf_iterator<char>(src)),
                                    std::istreambuf_iterator<char>());
            mutate(bytes);
            std::ofstream dst("tmp_orc1_pat.bin",
                              std::ios::binary | std::ios::trunc);
            dst.write(bytes.data(),
                      static_cast<std::streamsize>(bytes.size()));
        };
        auto rejected = [&](const char* what) {
            bool threw = false;
            try {
                Orc1Reader r2("tmp_orc1_pat.bin");
                r2.for_block(0, [](const Orc1Record&) {});
            } catch (const std::exception&) { threw = true; }
            if (!threw) std::printf("NOT REJECTED: %s\n", what);
            CHECK(threw);
        };
        patched([](std::vector<char>& b) { b.push_back(0); });
        rejected("trailing byte after footer");
        patched([](std::vector<char>& b) { b[12] = 1; });
        rejected("nonzero flags");
        patched([](std::vector<char>& b) { b[16] = static_cast<char>(b[16] + 1); });
        rejected("n_records off by one");
        // First record is {3,5,9}: bytes at body start are k=3, 0x03, 0x02,
        // 0x04; zeroing the last delta makes factors non-ascending.
        patched([](std::vector<char>& b) { b[kOrc1HeaderBytes + 3] = 0; });
        rejected("zero factor delta in block");
        std::remove("tmp_orc1_pat.bin");
    }

    std::remove("tmp_orc1_test.bin");
    std::remove("tmp_orc1_bad.bin");
}

int main() {
    test_vbyte();
    test_orc1();
    if (failures == 0) std::printf("all oracle tests passed\n");
    return failures == 0 ? 0 : 1;
}
