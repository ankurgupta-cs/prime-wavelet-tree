// trie_paths -- prime wavelet tree input, pass 1: from the certified
// ORC1 oracle (n + prime list per record) emit one PRIME PATH per record
// as a fixed-size record of `fields` x 5-byte big-endian primes (factors
// < 10^12 < 2^40), unused fields zero, ascending order. Big-endian packing
// makes memcmp order = lexicographic numeric order with a prefix sorting
// before its extensions. Output is in ORC1 (n) order; pwt2_encode sorts it.
//
// Two path modes (switch A of the A/B/C table):
//   default : d_min path -- the primes of d_min(n), the smallest divisor
//             with d*lambda(d) > n; n == d_min * (d_min^{-1} mod lambda)
//             is verified per record ("both techniques").
//   --full  : the FULL factorization of n, all omega(n) primes; nothing
//             to verify ("prime wavelet tree alone", the general structure).
//
// usage: trie_paths <table.orc1> <out.bin> [--full] [--fields K] [--threads N]
//   --fields defaults to 12 (d_min mode; max observed 12) or 14 (--full;
//   the table's max omega). A path longer than K fields is fatal.
//
// No factoring anywhere. An intermediate file, not a delivery format.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "u128.hpp"

using namespace cn;

namespace {

constexpr int kFieldBytes = 5;

struct Stats {
    u64 records = 0;
    std::vector<u64> pathlen_hist;   // index = path length
    u64 omega_hist[kOrc1MaxFactors + 1] = {0};
    u64 has_pmax = 0, has_pmin = 0;
    int max_pathlen = 0;
    u64 over = 0;
    explicit Stats(int fields) : pathlen_hist(fields + 2, 0) {}
    void add(const Stats& o) {
        records += o.records;
        for (std::size_t i = 0; i < pathlen_hist.size(); ++i) pathlen_hist[i] += o.pathlen_hist[i];
        for (int i = 0; i <= kOrc1MaxFactors; ++i) omega_hist[i] += o.omega_hist[i];
        has_pmax += o.has_pmax;
        has_pmin += o.has_pmin;
        max_pathlen = std::max(max_pathlen, o.max_pathlen);
        over += o.over;
    }
};

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);   // progress survives a killed run
    const char* orc_path = nullptr;
    const char* out_path = nullptr;
    unsigned threads = 0;
    bool full = false;
    int fields = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--threads") && i + 1 < argc)
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!std::strcmp(argv[i], "--fields") && i + 1 < argc)
            fields = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--full")) full = true;
        else if (!orc_path) orc_path = argv[i];
        else if (!out_path) out_path = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!orc_path || !out_path) {
        std::fprintf(stderr, "usage: %s <table.orc1> <out.bin> [--full] [--fields K] [--threads N]\n", argv[0]);
        return 2;
    }
    if (fields == 0) fields = full ? kOrc1MaxFactors : 12;
    if (fields < 1 || fields > 32) { std::fprintf(stderr, "bad --fields\n"); return 2; }
    const int rec_bytes = fields * kFieldBytes;
    if (threads == 0) {
        const unsigned hc = std::thread::hardware_concurrency();
        threads = hc > 2 ? hc - 1 : 1;
    }
    const auto t0 = std::chrono::steady_clock::now();
    u64 n_records = 0, n_blocks = 0;
    try {
        Orc1Reader rd(orc_path);
        n_records = rd.header().n_records;
        n_blocks = rd.index().size();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 2;
    }
    std::printf("records    %llu   blocks %llu   threads %u   mode %s   fields %d (%d B/record)\n",
                (unsigned long long)n_records, (unsigned long long)n_blocks, threads,
                full ? "FULL factorization" : "d_min path", fields, rec_bytes);

    // Pre-size the output so every thread can write its blocks at
    // rank_base * rec_bytes; the file ends up in n order.
    {
        std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
        if (!out) { std::fprintf(stderr, "cannot create %s\n", out_path); return 2; }
        if (n_records > 0) {
            out.seekp(static_cast<std::streamoff>(n_records * rec_bytes - 1));
            out.put('\0');
        }
        out.close();
        if (!out) { std::fprintf(stderr, "cannot size %s\n", out_path); return 2; }
    }

    threads = static_cast<unsigned>(std::min<u64>(threads, std::max<u64>(1, n_blocks)));
    std::vector<Stats> per(threads, Stats(fields));
    std::atomic<bool> failed{false};
    std::mutex err_mtx;
    std::string err;
    auto fail = [&](const std::string& m) {
        failed.store(true);
        std::lock_guard<std::mutex> g(err_mtx);
        if (err.empty()) err = m;
    };

    auto worker = [&](unsigned t, u64 b_lo, u64 b_hi) {
        try {
            Orc1Reader rd(orc_path);
            std::fstream out(out_path, std::ios::binary | std::ios::in | std::ios::out);
            if (!out) { fail("cannot open output for writing"); return; }
            std::vector<std::uint8_t> buf;
            std::vector<std::uint8_t> r(rec_bytes);
            Stats& st = per[t];
            for (u64 b = b_lo; b < b_hi && !failed.load(std::memory_order_relaxed); ++b) {
                buf.clear();
                const u64 rank_base = rd.index()[b].rank_base;
                rd.for_block(b, [&](const Orc1Record& rec) {
                    u128 d = 0;
                    if (!full) {
                        const DminChoice c = choose_dmin(rec.factors, rec.k, rec.n);
                        if (dmin_value(c.d, c.lambda) != rec.n) {
                            fail("RECONSTRUCTION FAILURE at n=" + to_string(rec.n));
                            return;
                        }
                        d = c.d;
                    }
                    std::fill(r.begin(), r.end(), 0);
                    int len = 0;
                    u128 prod = 1;
                    for (int j = 0; j < rec.k; ++j) {
                        if (!full && d % rec.factors[j] != 0) continue;
                        const u64 p = rec.factors[j];
                        if (p >= (1ull << 40)) { fail("prime >= 2^40 at n=" + to_string(rec.n)); return; }
                        if (len >= fields) { ++st.over; ++len; continue; }
                        std::uint8_t* f = r.data() + len * kFieldBytes;
                        for (int q = 0; q < kFieldBytes; ++q)
                            f[q] = static_cast<std::uint8_t>(p >> (8 * (kFieldBytes - 1 - q)));
                        prod *= p;
                        ++len;
                    }
                    if (!full && prod != d) { fail("path product != d at n=" + to_string(rec.n)); return; }
                    if (full && prod != rec.n) { fail("path product != n at n=" + to_string(rec.n)); return; }
                    ++st.records;
                    ++st.pathlen_hist[std::min<std::size_t>(len, st.pathlen_hist.size() - 1)];
                    ++st.omega_hist[rec.k];
                    st.max_pathlen = std::max(st.max_pathlen, len);
                    const u128 dd = full ? rec.n : d;
                    if (dd % rec.factors[rec.k - 1] == 0) ++st.has_pmax;
                    if (dd % rec.factors[0] == 0) ++st.has_pmin;
                    buf.insert(buf.end(), r.begin(), r.end());
                });
                if (failed.load()) return;
                out.seekp(static_cast<std::streamoff>(rank_base * rec_bytes));
                out.write(reinterpret_cast<const char*>(buf.data()),
                          static_cast<std::streamsize>(buf.size()));
                if (!out) { fail("write failure"); return; }
            }
            out.close();
        } catch (const std::exception& e) {
            fail(e.what());
        }
    };

    std::vector<std::thread> pool;
    const u64 per_t = (n_blocks + threads - 1) / threads;
    for (unsigned t = 0; t < threads; ++t) {
        const u64 lo = t * per_t, hi = std::min(n_blocks, lo + per_t);
        if (lo >= hi) break;
        pool.emplace_back(worker, t, lo, hi);
    }
    for (auto& th : pool) th.join();
    if (failed.load()) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }

    Stats all(fields);
    for (const Stats& s : per) all.add(s);
    if (all.records != n_records) {
        std::fprintf(stderr, "record count mismatch: %llu vs %llu\n",
                     (unsigned long long)all.records, (unsigned long long)n_records);
        return 1;
    }
    if (all.over) {
        std::fprintf(stderr, "FATAL: paths longer than --fields=%d (max %d); raise --fields\n",
                     fields, all.max_pathlen);
        return 1;
    }
    const double N = double(all.records);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("verified   %llu records (%s); paths emitted to %s\n",
                (unsigned long long)all.records,
                full ? "product == n" : "n == d_min * r*(d_min), product == d_min", out_path);
    std::printf("path len   max %d; histogram:", all.max_pathlen);
    for (int i = 1; i <= fields; ++i)
        if (all.pathlen_hist[i]) std::printf(" %d:%llu", i, (unsigned long long)all.pathlen_hist[i]);
    std::printf("\nomega(n)   histogram:");
    for (int i = 3; i <= kOrc1MaxFactors; ++i)
        if (all.omega_hist[i]) std::printf(" %d:%llu", i, (unsigned long long)all.omega_hist[i]);
    std::printf("\npath has   p_max(n) %.4f   p_min(n) %.4f\n", all.has_pmax / N, all.has_pmin / N);
    std::printf("elapsed    %.1f s\n", secs);
    return 0;
}
