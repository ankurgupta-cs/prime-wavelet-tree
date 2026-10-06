// cnd_decode -- decode and certify the CND1 container (docs/CND_FORMAT.md).
//
// usage: cnd_decode <in.cnd1> [--emit-d <out.raw>]
//                              [--certify <table.orc1>] [--threads N]
//
// Plain decode validates structure + payload sha and rebuilds the sorted
// d-set. --certify additionally factors every d (splitter.hpp; proofs, not
// probables, for d < psi_12), reconstructs n = d * (d^-1 mod lambda(d)),
// sorts the n's, and requires (a) sha of the 16B-LE ascending n-stream ==
// footer sha_nset, (b) elementwise equality with the n-sequence streamed
// from the ORC1 oracle, (c) total_check equality. Factoring parallelizes
// over index ranges (--threads, default hardware_concurrency - 1).

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "cnd1.hpp"
#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "splitter.hpp"
#include "u128.hpp"

using namespace cn;

int main(int argc, char** argv) {
    const char *cnd_path = nullptr, *emit_path = nullptr, *orc_path = nullptr;
    unsigned threads = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--emit-d") && i + 1 < argc) emit_path = argv[++i];
        else if (!std::strcmp(argv[i], "--certify") && i + 1 < argc) orc_path = argv[++i];
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc)
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!cnd_path) cnd_path = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!cnd_path) {
        std::fprintf(stderr,
                     "usage: %s <in.cnd1> [--emit-d <out.raw>] "
                     "[--certify <table.orc1>] [--threads N]\n", argv[0]);
        return 2;
    }
    try {
        Cnd1Reader rd(cnd_path);
        std::printf("records    %llu   M %u   end_width %u\n",
                    (unsigned long long)rd.records(), rd.M(), rd.end_width());
        const std::vector<u128> ds = rd.decode_all();
        if (ds.empty())
            std::printf("decoded    0 d values, payload sha OK (empty container)\n");
        else
            std::printf("decoded    %zu d values, payload sha OK, span [%s, %s]\n",
                        ds.size(), to_string(ds.front()).c_str(),
                        to_string(ds.back()).c_str());

        if (emit_path) {
            std::ofstream out(emit_path, std::ios::binary | std::ios::trunc);
            if (!out) { std::fprintf(stderr, "cannot write %s\n", emit_path); return 2; }
            std::vector<std::uint8_t> buf;
            buf.reserve(1 << 20);
            for (const u128 d : ds) {
                for (int b = 0; b < 16; ++b)
                    buf.push_back(static_cast<std::uint8_t>(d >> (8 * b)));
                if (buf.size() >= (1u << 20)) {
                    out.write(reinterpret_cast<const char*>(buf.data()),
                              static_cast<std::streamsize>(buf.size()));
                    buf.clear();
                }
            }
            out.write(reinterpret_cast<const char*>(buf.data()),
                      static_cast<std::streamsize>(buf.size()));
            out.close();
            if (!out) { std::fprintf(stderr, "write failure\n"); return 1; }
            std::printf("emitted    %s\n", emit_path);
        }

        if (!orc_path) return 0;

        // ---- certification ----
        if (threads == 0) {
            const unsigned hc = std::thread::hardware_concurrency();
            threads = hc > 2 ? hc - 1 : 1;
        }
        threads = static_cast<unsigned>(std::min<u64>(
            {threads, 1024, std::max<u64>(1, ds.size())}));
        std::vector<u128> ns(ds.size());
        std::atomic<bool> failed{false};
        std::mutex err_mtx;
        std::string err;
        auto worker = [&](std::size_t lo, std::size_t hi) {
            for (std::size_t i = lo; i < hi && !failed.load(std::memory_order_relaxed); ++i) {
                const u128 d = ds[i];
                const SplitResult f = split_u128(d);
                u128 lam = 1, prod = 1;
                for (const u128 p : f.primes) { lam = lcm128(lam, p - 1); prod *= p; }
                if (f.uncertified != 1 || prod != d || gcd128(d, lam) != 1) {
                    failed.store(true);
                    std::lock_guard<std::mutex> g(err_mtx);
                    err = "factor/lambda failure at d=" + to_string(d);
                    return;
                }
                ns[i] = d * inv_mod128(d, lam);
            }
        };
        std::vector<std::thread> pool;
        const std::size_t per = (ds.size() + threads - 1) / threads;
        for (unsigned t = 0; t < threads; ++t) {
            const std::size_t lo = static_cast<std::size_t>(t) * per;
            if (lo >= ds.size()) break;
            pool.emplace_back(worker, lo, std::min(ds.size(), lo + per));
        }
        for (auto& th : pool) th.join();
        if (failed) { std::fprintf(stderr, "CERTIFICATION FAILED: %s\n", err.c_str()); return 1; }
        std::printf("factored   all %zu d (threads %u), reconstructing order...\n",
                    ds.size(), threads);

        std::sort(ns.begin(), ns.end());
        Sha256 sha;
        u128 sum = 0;
        {
            std::vector<std::uint8_t> le(16 * 4096);
            std::size_t fill = 0;
            for (const u128 n : ns) {
                for (int b = 0; b < 16; ++b) le[fill + static_cast<std::size_t>(b)] =
                    static_cast<std::uint8_t>(n >> (8 * b));
                fill += 16;
                if (fill == le.size()) { sha.update(le.data(), fill); fill = 0; }
                sum += n;
            }
            if (fill) sha.update(le.data(), fill);
        }
        const auto digest = sha.finish();
        if (digest != rd.sha_nset()) {
            std::fprintf(stderr, "CERTIFICATION FAILED: n-stream sha mismatch\n");
            return 1;
        }
        if (sum != rd.total_check()) {
            std::fprintf(stderr, "CERTIFICATION FAILED: total_check mismatch\n");
            return 1;
        }

        u64 idx = 0;
        bool orc_ok = true;
        {
            Orc1Reader orc(orc_path);
            orc.for_all([&](const Orc1Record& rec) {
                if (!orc_ok) return;
                if (idx >= ns.size() || ns[idx] != rec.n) { orc_ok = false; return; }
                ++idx;
            });
        }
        if (!orc_ok || idx != ns.size()) {
            std::fprintf(stderr, "CERTIFICATION FAILED: mismatch vs oracle at rank %llu\n",
                         (unsigned long long)idx);
            return 1;
        }
        std::printf("sha nset   %s (match)\n", sha256_hex(digest).c_str());
        std::printf("oracle     %llu records elementwise equal\n",
                    (unsigned long long)idx);
        std::printf("round-trip CERTIFIED\n");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
