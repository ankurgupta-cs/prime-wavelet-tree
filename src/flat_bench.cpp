// flat_bench -- decode-time benchmark for the FLAT cells of the A/B/C table
// (the Decode column of cells 000/001/100/101), on a sorted 16-byte
// little-endian value file (data/n.raw = cell 000/001, data/d.raw = 100/101):
//
//   VByte gaps      encode the sorted gaps as VByte in RAM, time the decode
//                   back to the value array, verify equality (cells 000, 100)
//   per-class BIC   partition by value mod M, code every class with the
//                   CND1 coder (bic.hpp: ends + interpolative, exactly what
//                   cnd_encode writes), account the container exactly as
//                   CND1 would (64 B header + VByte directory of counts and
//                   blob lengths + blobs + 88 B footer), time the blob decode
//                   and the reassembly into one sorted array, verify
//                   (cells 001 with M = 10,810,800; 101 with M = 2310, which
//                   must reproduce the shipped 956,341,716-byte CND1 file)
//
// usage: flat_bench <values.raw> [M ...]      (default M = 2310)
// Times are single-threaded, whole set in RAM, on the calling machine.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#include "bic.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

using namespace cn;

namespace {

double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int bits_of(u128 v) { int b = 0; while (v > 0) { v >>= 1; ++b; } return b; }

void vb_put(std::vector<std::uint8_t>& out, u128 v) {
    while (v >= 0x80) { out.push_back(static_cast<std::uint8_t>(v | 0x80)); v >>= 7; }
    out.push_back(static_cast<std::uint8_t>(v));
}

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    if (argc < 2) { std::fprintf(stderr, "usage: %s <values.raw> [M ...]\n", argv[0]); return 2; }
    std::vector<u64> moduli;
    for (int i = 2; i < argc; ++i) moduli.push_back(std::strtoull(argv[i], nullptr, 10));
    if (moduli.empty()) moduli.push_back(2310);

    // ---- load -----------------------------------------------------------------
    std::vector<u128> vals;
    {
        std::ifstream in(argv[1], std::ios::binary | std::ios::ate);
        if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
        const std::streamsize fsize = in.tellg();
        if (fsize % 16) { std::fprintf(stderr, "truncated record\n"); return 1; }
        vals.resize(static_cast<std::size_t>(fsize / 16));
        in.seekg(0);
        std::vector<std::uint8_t> chunk(16u << 20);
        std::size_t w = 0;
        while (w < vals.size()) {
            in.read(reinterpret_cast<char*>(chunk.data()), static_cast<std::streamsize>(chunk.size()));
            const std::streamsize got = in.gcount();
            if (got <= 0) break;
            for (std::streamsize off = 0; off + 16 <= got; off += 16) {
                u128 v = 0;
                for (int b = 15; b >= 0; --b) v = (v << 8) | chunk[off + b];
                vals[w++] = v;
            }
        }
        if (w != vals.size()) { std::fprintf(stderr, "short read\n"); return 1; }
    }
    const std::size_t N = vals.size();
    const double Nd = double(N);
    for (std::size_t i = 1; i < N; ++i)
        if (vals[i] <= vals[i - 1]) { std::fprintf(stderr, "input not strictly ascending at %zu\n", i); return 1; }
    std::printf("records  %zu   (%s)\n\n", N, argv[1]);

    // ---- VByte gaps -------------------------------------------------------------
    {
        std::vector<std::uint8_t> stream;
        stream.reserve(N * 4);
        vb_put(stream, vals[0]);
        for (std::size_t i = 1; i < N; ++i) vb_put(stream, vals[i] - vals[i - 1]);
        std::vector<u128> back(N);
        const double t0 = now();
        {
            const std::uint8_t* p = stream.data();
            u128 cur = 0;
            for (std::size_t i = 0; i < N; ++i) {
                u128 g = 0;
                int shift = 0;
                std::uint8_t b;
                do { b = *p++; g |= static_cast<u128>(b & 0x7f) << shift; shift += 7; } while (b & 0x80);
                cur += g;
                back[i] = cur;
            }
        }
        const double t1 = now();
        const bool ok = back == vals;
        std::printf("VByte gaps      %llu B = %.3f bits/el   decode to the sorted values %.2f s   %s\n",
                    (unsigned long long)stream.size(), 8.0 * double(stream.size()) / Nd, t1 - t0,
                    ok ? "VERIFIED" : "MISMATCH");
        if (!ok) return 1;
    }

    // ---- per-class BIC, CND1 layout ---------------------------------------------
    for (const u64 M : moduli) {
        std::printf("\nM = %llu\n", (unsigned long long)M);
        std::vector<std::vector<u128>> cls(M);
        {
            std::vector<u64> cnt(M, 0);
            for (const u128 v : vals) ++cnt[static_cast<u64>(v % M)];
            for (u64 r = 0; r < M; ++r) cls[r].reserve(cnt[r]);
        }
        u128 max_k = 0;
        for (const u128 v : vals) {
            const u64 r = static_cast<u64>(v % M);
            const u128 k = (v - r) / M;
            cls[r].push_back(k);
            if (k > max_k) max_k = k;
        }
        const int end_width = std::max(1, bits_of(max_k));   // the CND1 rule (cnd1.hpp)
        std::vector<std::vector<std::uint8_t>> blobs(M);
        std::vector<std::uint8_t> dir;
        u64 body = 0, nonempty = 0;
        const double e0 = now();
        for (u64 r = 0; r < M; ++r) {
            if (!cls[r].empty()) { bic_append(blobs[r], cls[r], end_width); ++nonempty; }
            vb_put(dir, cls[r].size());
            if (!cls[r].empty()) vb_put(dir, blobs[r].size());
            body += blobs[r].size();
        }
        const double e1 = now();
        const u64 container = 64 + dir.size() + body + 88;
        std::printf("  classes %llu nonempty of %llu, end width %d; encode %.1f s\n",
                    (unsigned long long)nonempty, (unsigned long long)M, end_width, e1 - e0);
        std::printf("  container (CND1 layout: header 64 + directory %llu + blobs %llu + footer 88) = %llu B = %.3f bits/el\n",
                    (unsigned long long)dir.size(), (unsigned long long)body, (unsigned long long)container,
                    8.0 * double(container) / Nd);
        // decode: every blob back to its k-set, then values, then one sorted array
        const double d0 = now();
        std::vector<std::vector<u128>> dec(M);
        for (u64 r = 0; r < M; ++r)
            if (!cls[r].empty())
                dec[r] = bic_decode(blobs[r].data(), blobs[r].size(), cls[r].size(), end_width);
        const double d1 = now();
        std::vector<u128> merged;
        merged.reserve(N);
        for (u64 r = 0; r < M; ++r)
            for (const u128 k : dec[r]) merged.push_back(k * M + r);
        std::sort(merged.begin(), merged.end());
        const double d2 = now();
        const bool ok = merged == vals;
        std::printf("  decode: blobs -> k-sets %.1f s; + values and global sort %.1f s = %.1f s to the sorted values   %s\n",
                    d1 - d0, d2 - d1, d2 - d0, ok ? "VERIFIED" : "MISMATCH");
        if (!ok) return 1;
    }
    return 0;
}
