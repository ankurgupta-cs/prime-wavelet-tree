// dlist_from_paths -- the direct route to the sorted d_min list (cells 100
// and 101 of the A/B/C table), computed from the oracle (decoding cell
// 101's own file with cnd_decode --emit-d would be circular). Two inputs:
//
//   paths mode (default): dpaths.bin from trie_paths in d_min mode (per
//     record `fields` x 5-byte big-endian primes, ascending, zero-padded, in
//     ORC1 = n order) -> d = product of the record's primes
//   --orc1 mode: the certified ORC1 oracle -> d = d_min(n) by choose_dmin
//     (lambda_bucket.hpp), exactly the divisor whose primes trie_paths
//     writes; no 18.5 GB intermediate file
//
// then sort ascending and write d.raw: ascending d as 16-byte little-endian
// records (the layout of cnd_decode --emit-d and of every flat-set tool).
// On the full table the output must hash to 22b9b6ab... (certified d.raw).
//
// Asserted per record (paths mode): every field after the first zero field
// is zero; primes >= 3 and strictly ascending; d nonzero, odd and < 2^80.
// (--orc1 mode: d odd and < 2^80; d | n holds by construction, d being a
// product of a subset of the record's primes.) After the sort: strictly increasing
// (the d_min map is injective).
//
// usage: dlist_from_paths <dpaths.bin> <out d.raw> [--fields K] [--threads T]
//        dlist_from_paths --orc1 <table.orc1> <out d.raw> [--threads T]
//   --fields defaults to 12 (trie_paths' d_min default). --threads defaults
//   to 1 (single-threaded measurement).
//
// Stamps (stamps.hpp):
//   paths mode: read     = the input read calls, summed over 64 MiB chunks
//               multiply = parse + assert + multiply, summed over chunks
//   orc1 mode:  extract  = ORC1 block reads + record decode (n = product)
//                          + choose_dmin, one interleaved loop
//   both:       alloc    = the d array (zero-filled resize); paths mode
//                          also opens the input + 64 MiB read buffer, orc1
//                          mode also reads the ORC1 header + block index
//               sort     = std::sort (T > 1: T run sorts + merge rounds)
//               check    = strictly-increasing pass
//               write    = d.raw written and closed
//               total    = program start to closed output
// UNTIMED, after total: the output file is read back from disk and hashed.
// Parallelizable with --threads: multiply / extract (per chunk or block
// range) and sort (run sorts + merges); read and write are I/O; check is
// one sequential pass (could be split, ~0.3 s).

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "stamps.hpp"
#include "u128.hpp"

using namespace cn;

namespace {

constexpr int kFieldBytes = 5;
constexpr u128 kLimit80 = static_cast<u128>(1) << 80;

// First error wins; workers poll `failed`.
struct ErrorSlot {
    std::atomic<bool> failed{false};
    std::mutex mtx;
    std::string msg;
    void set(const std::string& m) {
        std::lock_guard<std::mutex> g(mtx);
        if (!failed.load()) msg = m;
        failed.store(true);
    }
};

// Product of one dpaths record; returns nullptr on success, else the reason.
const char* path_product(const std::uint8_t* r, int fields, u128& d_out) {
    u128 d = 1;
    u64 prev = 0;
    bool ended = false;
    for (int f = 0; f < fields; ++f, r += kFieldBytes) {
        const u64 p = (static_cast<u64>(r[0]) << 32) | (static_cast<u64>(r[1]) << 24) |
                      (static_cast<u64>(r[2]) << 16) | (static_cast<u64>(r[3]) << 8) | r[4];
        if (p == 0) { ended = true; continue; }
        if (ended) return "nonzero field after a zero field";
        if (p < 3) return "path prime < 3";
        if (p <= prev) return "path primes not strictly ascending";
        prev = p;
        d *= p;                                // d < 2^80, p < 2^40: no u128 wrap
        if (d >= kLimit80) return "d >= 2^80";
    }
    if (prev == 0) return "empty path";
    if ((d & 1) == 0) return "even d";
    d_out = d;
    return nullptr;
}

void multiply_range(const std::uint8_t* buf, std::size_t lo, std::size_t hi, int fields,
                    int rec_bytes, u128* dst, std::size_t rank_base, ErrorSlot& err) {
    for (std::size_t i = lo; i < hi; ++i) {
        u128 d = 0;
        if (const char* why = path_product(buf + i * static_cast<std::size_t>(rec_bytes), fields, d)) {
            err.set(std::string(why) + " at record " + std::to_string(rank_base + i));
            return;
        }
        dst[i] = d;
    }
}

// T run sorts in threads, then pairwise merge rounds (parallel per round).
void parallel_sort(std::vector<u128>& v, unsigned T) {
    const std::size_t N = v.size();
    if (T <= 1 || N < (std::size_t(1) << 16)) { std::sort(v.begin(), v.end()); return; }
    std::vector<std::size_t> b(T + 1);
    for (unsigned i = 0; i <= T; ++i) b[i] = N / T * i + std::min<std::size_t>(i, N % T);
    {
        std::vector<std::thread> pool;
        for (unsigned i = 0; i < T; ++i)
            pool.emplace_back([&, i] { std::sort(v.begin() + b[i], v.begin() + b[i + 1]); });
        for (auto& th : pool) th.join();
    }
    std::vector<u128> tmp(N);
    u128* src = v.data();
    u128* dst = tmp.data();
    while (b.size() > 2) {
        std::vector<std::size_t> nb{0};
        std::vector<std::thread> pool;
        for (std::size_t i = 0; i + 1 < b.size(); i += 2) {
            if (i + 2 < b.size()) {
                const std::size_t lo = b[i], mid = b[i + 1], hi = b[i + 2];
                pool.emplace_back([=] { std::merge(src + lo, src + mid, src + mid, src + hi, dst + lo); });
                nb.push_back(hi);
            } else {
                const std::size_t lo = b[i], hi = b[i + 1];
                pool.emplace_back([=] { std::copy(src + lo, src + hi, dst + lo); });
                nb.push_back(hi);
            }
        }
        for (auto& th : pool) th.join();
        std::swap(src, dst);
        b.swap(nb);
    }
    if (src != v.data()) v.swap(tmp);
}

} // namespace

int main(int argc, char** argv) {
    PhaseClock clk;
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* in_path = nullptr;
    const char* out_path = nullptr;
    bool orc1_mode = false;
    int fields = 12;
    unsigned threads = 1;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--fields") && i + 1 < argc) fields = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc)
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!std::strcmp(argv[i], "--orc1")) orc1_mode = true;
        else if (!in_path) in_path = argv[i];
        else if (!out_path) out_path = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!in_path || !out_path || fields < 1 || fields > 32 || threads < 1 || threads > 256) {
        std::fprintf(stderr,
                     "usage: %s <dpaths.bin> <out d.raw> [--fields K] [--threads T]\n"
                     "       %s --orc1 <table.orc1> <out d.raw> [--threads T]\n", argv[0], argv[0]);
        return 2;
    }
    const int rec_bytes = fields * kFieldBytes;
    try {
        std::vector<u128> ds;
        ErrorSlot err;
        if (!orc1_mode) {
            // ---- paths mode: stream the records, multiply per chunk --------------
            const u64 fsize = file_size_of(in_path);
            if (fsize % static_cast<u64>(rec_bytes)) {
                std::fprintf(stderr, "%s: size %llu is not a multiple of %d (--fields %d)\n", in_path,
                             (unsigned long long)fsize, rec_bytes, fields);
                return 1;
            }
            const std::size_t N = static_cast<std::size_t>(fsize / static_cast<u64>(rec_bytes));
            std::printf("input      %s: %zu records x %d B (fields %d), mode paths, threads %u\n",
                        in_path, N, rec_bytes, fields, threads);
            ds.resize(N);
            std::ifstream in(in_path, std::ios::binary);
            if (!in) { std::fprintf(stderr, "cannot open %s\n", in_path); return 2; }
            const std::size_t chunk_recs = std::max<std::size_t>(1, kIoChunk / static_cast<std::size_t>(rec_bytes));
            std::vector<std::uint8_t> buf(chunk_recs * static_cast<std::size_t>(rec_bytes));
            clk.stamp("alloc", "d array (zero-filled) + 64 MiB read buffer");
            double t_read = 0.0, t_mul = 0.0;
            std::size_t w = 0;
            while (w < N) {
                const std::size_t want = std::min(chunk_recs, N - w);
                double t = wall_now();
                read_exact(in, buf.data(), want * static_cast<std::size_t>(rec_bytes));
                const double t1 = wall_now();
                t_read += t1 - t;
                if (threads == 1) {
                    multiply_range(buf.data(), 0, want, fields, rec_bytes, ds.data() + w, w, err);
                } else {
                    std::vector<std::thread> pool;
                    for (unsigned k = 0; k < threads; ++k) {
                        const std::size_t lo = want * k / threads, hi = want * (k + 1) / threads;
                        if (lo < hi)
                            pool.emplace_back(multiply_range, buf.data(), lo, hi, fields, rec_bytes,
                                              ds.data() + w, w, std::ref(err));
                    }
                    for (auto& th : pool) th.join();
                }
                t_mul += wall_now() - t1;
                if (err.failed.load()) { std::fprintf(stderr, "FATAL: %s\n", err.msg.c_str()); return 1; }
                w += want;
            }
            clk.add("read", t_read, "input read calls, summed over 64 MiB chunks");
            clk.add("multiply", t_mul, "parse 5-byte fields + asserts + product, summed over chunks");
        } else {
            // ---- orc1 mode: d_min per record, straight from the oracle ------------
            u64 n_records = 0, n_blocks = 0;
            {
                Orc1Reader rd(in_path);
                n_records = rd.header().n_records;
                n_blocks = rd.index().size();
            }
            std::printf("input      %s: %llu records, %llu blocks, mode orc1 (choose_dmin), threads %u\n",
                        in_path, (unsigned long long)n_records, (unsigned long long)n_blocks, threads);
            ds.resize(static_cast<std::size_t>(n_records));
            clk.stamp("alloc", "ORC1 header/index + d array, zero-filled");
            const unsigned T = static_cast<unsigned>(std::min<u64>(threads, std::max<u64>(1, n_blocks)));
            std::vector<u64> seen(T, 0);
            auto worker = [&](unsigned t, u64 b_lo, u64 b_hi) {
                try {
                    Orc1Reader rd(in_path);
                    for (u64 b = b_lo; b < b_hi && !err.failed.load(std::memory_order_relaxed); ++b) {
                        u64 idx = rd.index()[b].rank_base;
                        rd.for_block(b, [&](const Orc1Record& rec) {
                            const DminChoice c = choose_dmin(rec.factors, rec.k, rec.n);
                            // d | n by construction (a product of a subset of
                            // rec.factors); only the cheap invariants here
                            if (c.d >= kLimit80 || (c.d & 1) == 0 || idx >= ds.size()) {
                                err.set("bad d_min at n=" + to_string(rec.n));
                                return;
                            }
                            ds[static_cast<std::size_t>(idx++)] = c.d;
                            ++seen[t];
                        });
                    }
                } catch (const std::exception& e) {
                    err.set(e.what());
                }
            };
            if (T == 1) {
                worker(0, 0, n_blocks);
            } else {
                std::vector<std::thread> pool;
                const u64 per = (n_blocks + T - 1) / T;
                for (unsigned t = 0; t < T; ++t) {
                    const u64 lo = t * per, hi = std::min(n_blocks, lo + per);
                    if (lo < hi) pool.emplace_back(worker, t, lo, hi);
                }
                for (auto& th : pool) th.join();
            }
            if (err.failed.load()) { std::fprintf(stderr, "FATAL: %s\n", err.msg.c_str()); return 1; }
            u64 total_seen = 0;
            for (const u64 s : seen) total_seen += s;
            if (total_seen != n_records) {
                std::fprintf(stderr, "FATAL: decoded %llu records, header says %llu\n",
                             (unsigned long long)total_seen, (unsigned long long)n_records);
                return 1;
            }
            clk.stamp("extract", "ORC1 block reads + record decode + choose_dmin");
        }

        // ---- sort, check, write -------------------------------------------------
        parallel_sort(ds, threads);
        clk.stamp("sort", threads == 1 ? "std::sort, one thread" : "run sorts + merge rounds");
        for (std::size_t i = 1; i < ds.size(); ++i)
            if (ds[i] <= ds[i - 1]) {
                std::fprintf(stderr, "FATAL: d not strictly increasing at sorted rank %zu (d = %s)\n", i,
                             to_string(ds[i]).c_str());
                return 1;
            }
        clk.stamp("check", "strictly increasing (d_min injective)");
        raw16_write(out_path, ds);
        clk.stamp("write", "16-byte LE records, closed");
        clk.total();

        // ---- UNTIMED: hash the file as it is on disk ------------------------------
        {
            std::ifstream in(out_path, std::ios::binary);
            if (!in) { std::fprintf(stderr, "cannot reopen %s\n", out_path); return 1; }
            Sha256 sha;
            std::vector<std::uint8_t> buf(kIoChunk);
            u64 bytes = 0;
            while (in) {
                in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
                const std::streamsize got = in.gcount();
                if (got <= 0) break;
                sha.update(buf.data(), static_cast<std::size_t>(got));
                bytes += static_cast<u64>(got);
            }
            std::printf("output     %s: %zu records, %llu B%s\n", out_path, ds.size(),
                        (unsigned long long)bytes, bytes == ds.size() * 16 ? "" : "  SIZE MISMATCH");
            if (!ds.empty())
                std::printf("span       [%s, %s]\n", to_string(ds.front()).c_str(), to_string(ds.back()).c_str());
            std::printf("sha256     %s  (untimed; full table must be 22b9b6ab918bca08...)\n",
                        sha256_hex(sha.finish()).c_str());
            if (bytes != ds.size() * 16) return 1;
        }
        std::printf("parallel   with --threads: multiply/extract and sort; read/write are I/O; check is one pass\n");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
