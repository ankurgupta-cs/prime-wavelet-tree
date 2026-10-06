// nstats_orc1 -- gap and trie (POM) measures of the raw sorted n-set,
// streamed from the certified ORC1 oracle (never from the text file).
// Measures what data-aware set encodings cost on the raw n-set before any
// Carmichael structure is used
// (the yardsticks of Gupta, Hon, Shah and Vitter, TCS 2007):
//
//   gap(S)   = sum log2(g_i + 1), the data-aware lower yardstick
//   delta Z  = prefix-code overhead of a decodable delta stream,
//              2 * sum ceil(log2(ceil(log2(g_i + 1)) + 1))
//   trie(S)  = edge count of the binary trie over fixed-width 80-bit keys
//              (POM measure; strie(S) <= gap(S) + 2n - 2 after a best shift)
//   VByte    = actual byte-aligned gap stream (our codec, for scale)
//
// usage: nstats_orc1 <table.orc1> [--emit-n <out.raw> [--emit-only]]

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>

#include "orc1.hpp"
#include "sha256.hpp"
#include "stamps.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

using namespace cn;

namespace {

int bits_of(u128 v) {
    int b = 0;
    while (v > 0) { v >>= 1; ++b; }
    return b;
}

// --emit-only (the extraction step of the flat n cells 000/001): ORC1 ->
// n.raw with NO statistics and NO SHA-256 (hash the output file separately). Kept: the cheap per-record
// invariants (ascending, odd) and the count/sum checks against the oracle.
// Stamps (stamps.hpp): open = header + block index; extract = ORC1 block
// reads + record decode (n = product of the factors) + 16-byte LE emit,
// with the output write calls inside it (their sum printed as "of which");
// close = final write + close; total.
int emit_only(PhaseClock& clk, const char* orc_path, const char* emit_path) {
    Orc1Reader rd(orc_path);
    std::ofstream emit(emit_path, std::ios::binary | std::ios::trunc);
    if (!emit) { std::fprintf(stderr, "cannot write %s\n", emit_path); return 2; }
    clk.stamp("open", "ORC1 header + block index; output created");
    std::vector<std::uint8_t> ebuf;
    ebuf.reserve(1u << 20);
    double t_write = 0.0;
    u128 prev = 0, sum = 0;
    u64 count = 0;
    bool order_ok = true, odd_ok = true;
    rd.for_all([&](const Orc1Record& rec) {
        const u128 n = rec.n;
        if (count > 0 && n <= prev) order_ok = false;
        if ((n & 1) == 0) odd_ok = false;
        std::uint8_t le[16];
        for (int i = 0; i < 16; ++i) le[i] = static_cast<std::uint8_t>(n >> (8 * i));
        ebuf.insert(ebuf.end(), le, le + 16);
        if (ebuf.size() >= (1u << 20)) {
            const double t = wall_now();
            emit.write(reinterpret_cast<const char*>(ebuf.data()), static_cast<std::streamsize>(ebuf.size()));
            t_write += wall_now() - t;
            ebuf.clear();
        }
        sum += n;
        prev = n;
        ++count;
    });
    clk.stamp("extract", "ORC1 block reads + record decode + LE emit (1 MiB write calls included)");
    std::printf("  of which  %10.3f s   output write calls, summed\n", t_write);
    emit.write(reinterpret_cast<const char*>(ebuf.data()), static_cast<std::streamsize>(ebuf.size()));
    emit.close();
    if (!emit) { std::fprintf(stderr, "write failure on %s\n", emit_path); return 1; }
    clk.stamp("close", "last write + close");
    clk.total();
    const bool count_ok = count == rd.header().n_records;
    const bool sum_ok = sum == rd.footer().total_check;
    std::printf("stream check: ascending %s, all odd %s, count == header %s, sum == oracle total_check %s\n",
                order_ok ? "OK" : "FAIL", odd_ok ? "OK" : "FAIL", count_ok ? "OK" : "FAIL", sum_ok ? "OK" : "FAIL");
    std::printf("emitted      %llu records (%llu B) to %s; --emit-only: no statistics, no SHA-256 "
                "(hash the file: full table ae272474...)\n",
                (unsigned long long)count, (unsigned long long)(count * 16), emit_path);
    return (order_ok && odd_ok && count_ok && sum_ok) ? 0 : 1;
}

} // namespace

int main(int argc, char** argv) {
    PhaseClock clk;
    // usage: nstats_orc1 <table.orc1> [--emit-n <out.raw> [--emit-only]]
    //   --emit-n writes the ascending n stream as 16-byte little-endian
    //   records (the format every flat-set tool reads, and byte-for-byte the
    //   stream whose SHA-256 the certified CND1 footer stores as sha_nset).
    //   --emit-only (with --emit-n) skips every statistic and the SHA-256:
    //   the lean extraction step, with phase stamps.
    const char* emit_path = nullptr;
    bool only = false;
    bool bad = argc < 2;
    for (int i = 2; i < argc && !bad; ++i) {
        if (!std::strcmp(argv[i], "--emit-n") && i + 1 < argc && !emit_path) emit_path = argv[++i];
        else if (!std::strcmp(argv[i], "--emit-only") && !only) only = true;
        else bad = true;
    }
    if (bad || (only && !emit_path)) {
        std::fprintf(stderr, "usage: %s <table.orc1> [--emit-n <out.raw> [--emit-only]]\n", argv[0]);
        return 2;
    }
    if (only) {
        std::setvbuf(stdout, nullptr, _IONBF, 0);
        try {
            return emit_only(clk, argv[1], emit_path);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
    }
    try {
        Orc1Reader rd(argv[1]);
        constexpr int kWidth = 80;             // n < 10^24 < 2^80
        u128 prev = 0;
        u64 count = 0, trie_edges = 0, vb_bytes = 0, delta_bits = 0;
        u64 vb_abs_bytes = 0;                  // VByte(n) with NO differencing
        std::ofstream emit;
        std::vector<std::uint8_t> ebuf;
        if (emit_path) {
            emit.open(emit_path, std::ios::binary | std::ios::trunc);
            if (!emit) { std::fprintf(stderr, "cannot write %s\n", emit_path); return 2; }
            ebuf.reserve(1u << 20);
        }
        Sha256 sha;
        u128 sum = 0;
        bool order_ok = true, odd_ok = true;
        double gap_bits = 0.0;
        u128 first_n = 0, last_n = 0;
        // len+raw (order-0 entropy of the gap's bit length + raw bits): the
        // realizable flat-coder family used for every other row of the DCC
        // A/B table (rANS matched it to 0.002 bits/el on the d-set).
        // Once on the gaps as they are, once with the odd-only variant
        // (every Carmichael number is odd, so gaps are even: code g/2).
        u64 len_hist[130] = {0}, len_hist_odd[130] = {0};
        double raw_bits = 0.0, raw_bits_odd = 0.0;
        rd.for_all([&](const Orc1Record& rec) {
            const u128 n = rec.n;
            if (count == 0) {
                first_n = n;
                trie_edges += kWidth;
                vb_bytes += vbyte_length(n);
                gap_bits += std::log2(static_cast<double>(n) + 1.0);
            } else {
                const u128 g = n - prev;       // reader enforces ascending n
                const int gb = bits_of(g);
                gap_bits += std::log2(static_cast<double>(g) + 1.0);
                delta_bits += 2 * bits_of(static_cast<u128>(gb));
                trie_edges += bits_of(n ^ prev);   // = width - lcp(prev, n)
                vb_bytes += vbyte_length(g);
                ++len_hist[gb];
                raw_bits += gb - 1;
                const int gbo = bits_of(g >> 1);   // g even (n odd): g/2 >= 1
                ++len_hist_odd[gbo];
                raw_bits_odd += gbo - 1;
            }
            // the reader checks ascending order only across block boundaries;
            // assert it (and oddness) per record here
            if (count > 0 && n <= prev) order_ok = false;
            if ((n & 1) == 0) odd_ok = false;
            vb_abs_bytes += vbyte_length(n);
            std::uint8_t le[16];
            for (int i = 0; i < 16; ++i) le[i] = static_cast<std::uint8_t>(n >> (8 * i));
            sha.update(le, 16);
            sum += n;
            if (emit_path) {
                ebuf.insert(ebuf.end(), le, le + 16);
                if (ebuf.size() >= (1u << 20)) {
                    emit.write(reinterpret_cast<const char*>(ebuf.data()),
                               static_cast<std::streamsize>(ebuf.size()));
                    ebuf.clear();
                }
            }
            prev = n;
            last_n = n;
            ++count;
        });
        if (emit_path) {
            emit.write(reinterpret_cast<const char*>(ebuf.data()),
                       static_cast<std::streamsize>(ebuf.size()));
            emit.close();
            if (!emit) { std::fprintf(stderr, "write failure on %s\n", emit_path); return 1; }
        }
        const bool count_ok = count == rd.header().n_records;
        const bool sum_ok = sum == rd.footer().total_check;
        std::printf("stream check: ascending %s, all odd %s, count == header %s, sum == oracle total_check %s\n",
                    order_ok ? "OK" : "FAIL", odd_ok ? "OK" : "FAIL",
                    count_ok ? "OK" : "FAIL", sum_ok ? "OK" : "FAIL");
        std::printf("sha256 of the 16-byte-LE n stream: %s  (must equal sha_nset in the CND1 footer)\n",
                    sha256_hex(sha.finish()).c_str());
        if (emit_path) std::printf("emitted      %llu records to %s\n", (unsigned long long)count, emit_path);
        if (!(order_ok && odd_ok && count_ok && sum_ok)) return 1;
        const double N = static_cast<double>(count);
        std::printf("records      %llu\n", (unsigned long long)count);
        std::printf("n span       [%s, %s]\n",
                    to_string(first_n).c_str(), to_string(last_n).c_str());
        std::printf("gap(S)       %.3f bits/el  (%.3f GiB = %.3f GB decimal)\n",
                    gap_bits / N, gap_bits / 8.0 / (1u << 30), gap_bits / 8.0 / 1e9);
        std::printf("delta Z(S)   +%.3f bits/el decodable overhead\n",
                    static_cast<double>(delta_bits) / N);
        std::printf("trie(S)      %.3f edges/el (strie <= gap + 2n - 2 after shift)\n",
                    static_cast<double>(trie_edges) / N);
        std::printf("VByte gaps   %llu B = %.3f bits/el  (%.3f GiB = %.3f GB decimal)\n",
                    (unsigned long long)vb_bytes, vb_bytes * 8.0 / N,
                    static_cast<double>(vb_bytes) / (1u << 30), static_cast<double>(vb_bytes) / 1e9);
        std::printf("VByte absolute (NO differencing) %llu B = %.3f bits/el; fixed 80-bit = %.0f B\n",
                    (unsigned long long)vb_abs_bytes, vb_abs_bytes * 8.0 / N, N * 10.0);
        auto h0 = [](const u64* h) {
            u64 tot = 0;
            for (int i = 0; i < 130; ++i) tot += h[i];
            double b = 0.0;
            for (int i = 0; i < 130; ++i)
                if (h[i]) b -= double(h[i]) * std::log2(double(h[i]) / double(tot));
            return b;
        };
        const double lr = h0(len_hist) + raw_bits, lro = h0(len_hist_odd) + raw_bits_odd;
        std::printf("len+raw n-gaps         %.3f bits/el = %.0f B  (order-0 bit-length entropy + raw bits; realizable)\n",
                    lr / N, lr / 8.0);
        std::printf("len+raw odd-only gaps  %.3f bits/el = %.0f B  (same, coding g/2: the parity lever)\n",
                    lro / N, lro / 8.0);
        std::printf("binomial floor reference: 53.4 bits/el\n");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
