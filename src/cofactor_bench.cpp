// cofactor_bench -- the whole-file cost of a divisor-mode prime wavelet
// tree with ANY stored divisors (d_min or d*): the tree yields every
// n and the primes of its stored divisor d; the whole factored file also
// needs the primes of the cofactor n / d, factored with Shallue-Webster
// Algorithm 1 (cnfactor.hpp; the exponent n - 1 is known). This tool reads
// the certified ORC1 oracle and a path file in ORC1 order (trie_paths /
// dsel_set2 --emit layout: K five-byte big-endian primes per record,
// ascending, zero padded), and for a sample of blocks times the cofactor
// factorization single-threaded, checking every result against the oracle.
//
// usage: cofactor_bench <table.orc1> <paths.bin> --fields 12|14 [--stride S]
//   (--stride 100 = every 100th ORC1 block, the sample used for the A/B/C
//    table's per-number costs)

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "cnfactor.hpp"
#include "orc1.hpp"

using namespace cn;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* orc = nullptr;
    const char* paths = nullptr;
    int K = 0;
    u64 stride = 100;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--fields") && i + 1 < argc) K = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--stride") && i + 1 < argc) stride = std::strtoull(argv[++i], nullptr, 10);
        else if (!orc) orc = argv[i];
        else if (!paths) paths = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!orc || !paths || (K != 12 && K != 14) || stride == 0) {
        std::fprintf(stderr, "usage: %s <table.orc1> <paths.bin> --fields 12|14 [--stride S]\n", argv[0]);
        return 2;
    }
    try {
        Orc1Reader rd(orc);
        const auto& ix = rd.index();
        const u64 N = rd.header().n_records;
        std::ifstream pin(paths, std::ios::binary | std::ios::ate);
        if (!pin) { std::fprintf(stderr, "cannot open %s\n", paths); return 2; }
        const std::size_t R = static_cast<std::size_t>(K) * 5;
        if (static_cast<u64>(pin.tellg()) != N * R) { std::fprintf(stderr, "path file size != records x %zu\n", R); return 1; }
        const CnFactorizer cnf;
        auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
        double secs = 0.0;
        u64 recs = 0, ladders = 0, fallbacks = 0, trivial = 0, rprimes = 0, dmin_like = 0;
        std::vector<std::uint8_t> buf;
        for (u64 b = 0; b < ix.size(); b += stride) {
            const u64 first = ix[b].rank_base;
            const u64 cnt = (b + 1 < ix.size() ? ix[b + 1].rank_base : N) - first;
            buf.resize(cnt * R);
            pin.seekg(static_cast<std::streamoff>(first * R));
            pin.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
            if (pin.gcount() != static_cast<std::streamsize>(buf.size())) { std::fprintf(stderr, "short read\n"); return 1; }
            u64 i = 0;
            rd.for_block(b, [&](const Orc1Record& rec) {
                const std::uint8_t* pr = buf.data() + i * R;
                ++i;
                u128 d = 1;
                std::vector<u64> dp;
                for (int f = 0; f < K; ++f) {
                    u64 v = 0;
                    for (int q = 0; q < 5; ++q) v = (v << 8) | pr[f * 5 + q];
                    if (!v) break;
                    d *= v;
                    dp.push_back(v);
                }
                if (dp.empty() || rec.n % d != 0) throw std::runtime_error("path is not a divisor of its n (file out of ORC1 order?)");
                std::vector<u64> want;
                for (int f = 0; f < rec.k; ++f) {
                    bool in = false;
                    for (const u64 q : dp) if (q == rec.factors[f]) in = true;
                    if (!in) want.push_back(rec.factors[f]);
                }
                if (want.size() + dp.size() != static_cast<std::size_t>(rec.k)) throw std::runtime_error("path primes not a subset of n's primes");
                const u128 r = rec.n / d;
                std::vector<u64> got;
                const double s0 = now();
                if (r > 1) {
                    const CnFactorResult fr = cnf.factor(r, rec.n);
                    got = fr.primes;
                    ladders += static_cast<u64>(fr.bases_used);
                    if (fr.fell_back) ++fallbacks;
                }
                secs += now() - s0;
                if (r == 1) ++trivial;
                if (got != want) throw std::runtime_error("cofactor factorization differs from the oracle's primes");
                rprimes += want.size();
                ++recs;
                (void)dmin_like;
            });
            if (i != cnt) throw std::runtime_error("record count mismatch in a block");
        }
        const double us = 1e6 * secs / double(recs);
        std::printf("paths %s (fields %d); sampled %llu records (every %llu-th ORC1 block), all cofactor factorizations equal the oracle\n",
                    paths, K, (unsigned long long)recs, (unsigned long long)stride);
        std::printf("cofactor primes %.3f per number, trivial cofactors %llu, Algorithm 1 ladders %.2f per number, generic fallbacks %llu\n",
                    double(rprimes) / double(recs), (unsigned long long)trivial, double(ladders) / double(recs),
                    (unsigned long long)fallbacks);
        std::printf("factor the cofactor: %.3f us/number single-thread -> %.1f min for N = %llu (extrapolated)\n",
                    us, us * double(N) / 60e6, (unsigned long long)N);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
