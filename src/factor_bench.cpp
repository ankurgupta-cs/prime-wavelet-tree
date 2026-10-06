// factor_bench -- what it costs to get from each A/B/C cell's decoded
// values to the WHOLE table content (every n with its complete prime
// factorization), measured with the project's certified splitter
// (splitter.hpp: tiny-prime peel, the strong-Fermat gcd trick, Brent rho,
// MR-12 proofs) and verified record by record against the ORC1 oracle.
//
//   --mode n      the cell stores n (000, 001): factor n from scratch
//   --mode dmin   the cell stores d_min (100, 101): factor d -> lambda ->
//                 n = d * (d^-1 mod lambda) [phase A: to n], then factor
//                 the cofactor r* = n / d [phase B: to the factorization]
//   --mode rstar  the cell stores the tree over d_min (110, 111): the
//                 decoder already has n and the primes of d; factor r* only
//   (010 and 011 need nothing: the tree over full factorizations delivers
//   every prime as it decodes.)
//
// usage: factor_bench <table.orc1> --mode n|dmin|rstar [--threads T]
//                     [--stride S] [--max-blocks K]
//   --stride S processes every S-th oracle block (4096 records each) for a
//   sample; --threads 1 gives clean per-record CPU costs, more threads the
//   wall time on this machine. Every factorization is compared with the
//   oracle's factor list; any mismatch is fatal.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "cnfactor.hpp"
#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "splitter.hpp"
#include "u128.hpp"

using namespace cn;

namespace {

double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

struct Tally {
    u64 records = 0, mismatches = 0;
    double a_secs = 0.0, b_secs = 0.0;   // CPU seconds in phase A / B (per thread, summed)
    u64 ladders = 0, fallbacks = 0;      // --cn statistics
};

bool same_factors(const std::vector<u128>& got, const Orc1Record& rec) {
    if (got.size() != static_cast<std::size_t>(rec.k)) return false;
    for (int i = 0; i < rec.k; ++i) if (got[i] != rec.factors[i]) return false;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* orc = nullptr;
    std::string mode;
    unsigned threads = 1;
    u64 stride = 1, max_blocks = ~0ull;
    bool cn = false;   // --cn: Shallue-Webster Algorithm 1 (cnfactor.hpp) for n and r*; d still needs the generic splitter
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--cn")) cn = true;
        else if (!std::strcmp(argv[i], "--mode") && i + 1 < argc) mode = argv[++i];
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!std::strcmp(argv[i], "--stride") && i + 1 < argc) stride = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--max-blocks") && i + 1 < argc) max_blocks = std::strtoull(argv[++i], nullptr, 10);
        else if (!orc) orc = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!orc || (mode != "n" && mode != "dmin" && mode != "rstar") || threads == 0 || stride == 0) {
        std::fprintf(stderr, "usage: %s <table.orc1> --mode n|dmin|rstar [--threads T] [--stride S] [--max-blocks K]\n", argv[0]);
        return 2;
    }
    u64 n_blocks = 0, n_records = 0;
    try {
        Orc1Reader rd(orc);
        n_blocks = rd.index().size();
        n_records = rd.header().n_records;
    } catch (const std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 2; }
    std::vector<u64> blocks;
    for (u64 b = 0; b < n_blocks && blocks.size() < max_blocks; b += stride) blocks.push_back(b);
    std::printf("mode %s%s   oracle %llu records in %llu blocks   sampling %zu blocks (stride %llu)   threads %u\n",
                mode.c_str(), cn ? " (Shallue-Webster Algorithm 1 for n and r*)" : " (generic splitter)",
                (unsigned long long)n_records, (unsigned long long)n_blocks, blocks.size(),
                (unsigned long long)stride, threads);
    const CnFactorizer cnf;   // shared, read-only after construction (bitmap to 1e8 + 100 bases)

    std::vector<Tally> per(threads);
    std::atomic<bool> failed{false};
    std::mutex err_mtx;
    std::string err;
    auto fail = [&](const std::string& m) { failed.store(true); std::lock_guard<std::mutex> g(err_mtx); if (err.empty()) err = m; };
    const double t0 = now();
    auto worker = [&](unsigned t) {
        try {
            Orc1Reader rd(orc);
            Tally& ty = per[t];
            for (std::size_t bi = t; bi < blocks.size() && !failed.load(std::memory_order_relaxed); bi += threads) {
                rd.for_block(blocks[bi], [&](const Orc1Record& rec) {
                    std::vector<u128> got;
                    auto cn_split = [&](u128 x) {   // Algorithm 1 with the exponent n - 1
                        const CnFactorResult r = cnf.factor(x, rec.n);
                        ty.ladders += static_cast<u64>(r.bases_used);
                        if (r.fell_back) ++ty.fallbacks;
                        got.insert(got.end(), r.primes.begin(), r.primes.end());
                    };
                    if (mode == "n") {
                        const double s0 = now();
                        if (cn) {
                            cn_split(rec.n);
                        } else {
                            SplitResult f = split_u128(rec.n);
                            if (f.uncertified != 1) { fail("uncertified leaf at n=" + to_string(rec.n)); return; }
                            got = std::move(f.primes);
                        }
                        ty.a_secs += now() - s0;
                    } else {
                        const DminChoice c = choose_dmin(rec.factors, rec.k, rec.n);   // the stored value, not timed
                        const u128 r = rec.n / c.d;
                        if (mode == "dmin") {
                            const double s0 = now();
                            SplitResult fd = split_u128(c.d);
                            u128 lam = 1, prod = 1;
                            for (const u128 p : fd.primes) { lam = lcm128(lam, p - 1); prod *= p; }
                            const u128 n2 = prod * inv_mod128(prod, lam);
                            ty.a_secs += now() - s0;
                            if (fd.uncertified != 1 || prod != c.d || n2 != rec.n) { fail("d path failure at n=" + to_string(rec.n)); return; }
                            got = std::move(fd.primes);
                        } else {
                            for (int i = 0; i < rec.k; ++i) if (c.d % rec.factors[i] == 0) got.push_back(rec.factors[i]);
                        }
                        const double s1 = now();
                        if (cn) {
                            cn_split(r);
                        } else {
                            SplitResult fr = split_u128(r);
                            if (fr.uncertified != 1) { fail("uncertified leaf in r* at n=" + to_string(rec.n)); return; }
                            got.insert(got.end(), fr.primes.begin(), fr.primes.end());
                        }
                        ty.b_secs += now() - s1;
                        std::sort(got.begin(), got.end());
                    }
                    ++ty.records;
                    if (!same_factors(got, rec)) { ++ty.mismatches; fail("factorization mismatch at n=" + to_string(rec.n)); }
                });
            }
        } catch (const std::exception& e) { fail(e.what()); }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    const double wall = now() - t0;
    if (failed.load()) { std::fprintf(stderr, "FAILED: %s\n", err.c_str()); return 1; }
    Tally all;
    for (const Tally& t : per) {
        all.records += t.records; all.mismatches += t.mismatches; all.a_secs += t.a_secs; all.b_secs += t.b_secs;
        all.ladders += t.ladders; all.fallbacks += t.fallbacks;
    }
    const double R = double(all.records), N = double(n_records);
    std::printf("records %llu   mismatches %llu   wall %.1f s on %u threads\n",
                (unsigned long long)all.records, (unsigned long long)all.mismatches, wall, threads);
    if (cn)
        std::printf("Algorithm 1: %.2f ladders per record, %llu records needed the generic fallback\n",
                    double(all.ladders) / R, (unsigned long long)all.fallbacks);
    if (mode == "n")
        std::printf("factor n from scratch: %.2f us/record CPU -> full table %.2f h single-thread, %.1f min on %u threads (extrapolated from this sample)\n",
                    1e6 * all.a_secs / R, all.a_secs / R * N / 3600.0, all.a_secs / R * N / 60.0 / threads, threads);
    if (mode == "dmin")
        std::printf("phase A (factor d, lambda, n): %.2f us/record -> %.2f h single-thread; phase B (factor r*): %.2f us/record -> %.2f h; both %.2f h single-thread, %.1f min on %u threads (extrapolated)\n",
                    1e6 * all.a_secs / R, all.a_secs / R * N / 3600.0, 1e6 * all.b_secs / R, all.b_secs / R * N / 3600.0,
                    (all.a_secs + all.b_secs) / R * N / 3600.0, (all.a_secs + all.b_secs) / R * N / 60.0 / threads, threads);
    if (mode == "rstar")
        std::printf("factor r* only: %.2f us/record -> full table %.2f h single-thread, %.1f min on %u threads (extrapolated)\n",
                    1e6 * all.b_secs / R, all.b_secs / R * N / 3600.0, all.b_secs / R * N / 60.0 / threads, threads);
    std::printf("every factorization equals the oracle's factor list\n");
    return 0;
}
