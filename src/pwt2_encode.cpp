// pwt2_encode -- build a PWT2 container (docs/PWT2_FORMAT.md): the
// frequency-ordered prime wavelet tree. Input = prime paths in the
// trie_paths layout (K five-byte big-endian primes per record, ascending,
// zero padded; any record order); certification targets from the ORC1
// oracle exactly as pwt_encode / cnd_encode.
//
// usage: pwt2_encode <paths.bin> <table.orc1> <out.pwt2> [--fields 12|14]
//                    [--full] [--naive | --ctx-set 1|2|3] [--sieve-limit B]
//                    [--scratch DIR] [--chunk R] [--threads N]
//   --full    the paths are full factorizations (n = product); default: each
//             path is a divisor d of n with d*lambda(d) > n (d_min or any
//             other valid choice; n = d * (d^-1 mod lambda(d)))
//   --naive   byte-aligned variant (every symbol one VByte; cells 010/110)
//   --ctx-set 1|2|3  entropy-coded context set (default 1; 2 = the richer
//             sparse contexts, docs/PWT2_FORMAT.md section 9;
//             3 = set 2 + tree-gap mantissa models, section 10)
// Steps: count path occurrences of every prime -> global order (count
// descending, ties to the larger prime) -> ranked copy of the paths in
// --scratch -> external sort in trie pre-order -> two-pass encode.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "orc1.hpp"
#include "pwt2.hpp"
#include "pwt_paths.hpp"

using namespace cn;

namespace {

constexpr u64 kSmall = u64(1) << 27;   // direct-indexed prime counts / ranks below this

template <class F> void for_chunks(const char* path, int K, F&& f) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error(std::string("cannot open ") + path);
    const std::size_t R = std::size_t(K) * kPathFieldBytes;
    const std::size_t per = (std::size_t(512) << 20) / R;
    std::vector<std::uint8_t> buf(per * R);
    for (;;) {
        in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
        const std::size_t got = static_cast<std::size_t>(in.gcount());
        if (got % R) throw std::runtime_error("partial path record (wrong --fields?)");
        if (got) f(buf.data(), got / R);
        if (got < buf.size()) break;
    }
}
inline u64 getf(const std::uint8_t* r, int i) {
    u64 v = 0;
    for (int q = 0; q < kPathFieldBytes; ++q) v = (v << 8) | r[i * kPathFieldBytes + q];
    return v;
}
inline void putf(std::uint8_t* r, int i, u64 v) {
    for (int q = 0; q < kPathFieldBytes; ++q) r[i * kPathFieldBytes + q] = static_cast<std::uint8_t>(v >> (8 * (kPathFieldBytes - 1 - q)));
}

template <int K>
int run(const char* paths, const char* orc, const char* out_path, bool full, bool naive, u64 B,
        PathSortOptions sopt, int ctx_set) {
    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    std::printf("input       %s (fields %d, %s mode, %s, context set %d)\n", paths, K, full ? "full-factorization" : "divisor (d*lambda(d) > n)",
                naive ? "PWT2 byte-aligned VByte streams" : "PWT2 entropy-coded", naive ? 1 : ctx_set);

    // ---- 1. prime path counts ------------------------------------------------
    std::vector<u64> small(kSmall >> 1, 0);
    std::vector<u64> big;
    u64 recs = 0;
    for_chunks(paths, K, [&](const std::uint8_t* b, std::size_t n) {
        for (std::size_t r = 0; r < n; ++r) {
            const std::uint8_t* rec = b + r * K * kPathFieldBytes;
            u64 prev = 0;
            int len = 0;
            for (; len < K; ++len) {
                const u64 p = getf(rec, len);
                if (!p) break;
                if (p <= prev || !(p & 1) || p < 3) throw std::runtime_error("path not strictly ascending odd primes >= 3");
                prev = p;
                if (p < kSmall) ++small[p >> 1]; else big.push_back(p);
            }
            if (len == 0) throw std::runtime_error("empty path record");
            for (int i = len; i < K; ++i) if (getf(rec, i)) throw std::runtime_error("nonzero field after the padding");
            ++recs;
        }
    });
    std::sort(big.begin(), big.end());
    std::vector<std::pair<u64, u64>> pc;
    for (u64 i = 1; i < (kSmall >> 1); ++i) if (small[i]) pc.emplace_back(2 * i + 1, small[i]);
    for (std::size_t i = 0; i < big.size();) {
        std::size_t j = i;
        while (j < big.size() && big[j] == big[i]) ++j;
        pc.emplace_back(big[i], j - i);
        i = j;
    }
    big.clear(); big.shrink_to_fit();
    const Pwt2Order ord = pwt2_order_from_counts(pc);
    std::printf("order       %llu records, %llu distinct primes in %zu count classes; root-most:",
                (unsigned long long)recs, (unsigned long long)ord.m(), ord.class_size.size());
    for (u64 r = 1; r <= std::min<u64>(8, ord.m()); ++r) std::printf(" %llu", (unsigned long long)ord.prime[r]);
    std::printf(" (%.0f s)\n", secs());

    // ---- 2. ranked copy -------------------------------------------------------
    std::vector<u32> rank_small(kSmall >> 1, 0);
    std::vector<std::pair<u64, u32>> rank_big;
    for (u64 r = 1; r <= ord.m(); ++r) {
        const u64 p = ord.prime[r];
        if (p < kSmall) rank_small[p >> 1] = static_cast<u32>(r);
        else rank_big.emplace_back(p, static_cast<u32>(r));
    }
    std::sort(rank_big.begin(), rank_big.end());
    auto rank_of = [&](u64 p) -> u64 {
        if (p < kSmall) return rank_small[p >> 1];
        const auto it = std::lower_bound(rank_big.begin(), rank_big.end(), std::make_pair(p, u32(0)));
        if (it == rank_big.end() || it->first != p) return 0;
        return it->second;
    };
    const std::string ranked = sopt.scratch + "/" + sopt.tag + "_ranked.bin";
    {
        std::ofstream o(ranked, std::ios::binary | std::ios::trunc);
        if (!o) throw std::runtime_error("cannot create " + ranked);
        const std::size_t R = std::size_t(K) * kPathFieldBytes;
        const unsigned T = std::max(1u, std::min(sopt.threads ? sopt.threads : 4u, 8u));
        std::vector<std::uint8_t> ob;
        bool bad = false;
        for_chunks(paths, K, [&](const std::uint8_t* b, std::size_t n) {
            ob.assign(n * R, 0);
            std::vector<std::thread> pool;
            std::vector<char> bd(T, 0);
            const std::size_t part = (n + T - 1) / T;
            for (unsigned t = 0; t < T; ++t)
                pool.emplace_back([&, t] {
                    const std::size_t lo = t * part, hi = std::min(n, lo + part);
                    u64 rk[kPathMaxFields];
                    for (std::size_t r = lo; r < hi; ++r) {
                        const std::uint8_t* rec = b + r * R;
                        int len = 0;
                        for (; len < K; ++len) {
                            const u64 p = getf(rec, len);
                            if (!p) break;
                            rk[len] = rank_of(p);
                            if (!rk[len]) bd[t] = 1;
                        }
                        std::sort(rk, rk + len);
                        for (int i = 0; i < len; ++i) putf(ob.data() + r * R, i, rk[i]);
                    }
                });
            for (auto& th : pool) th.join();
            for (const char x : bd) if (x) bad = true;
            for (std::size_t w = 0; w < ob.size();) {
                const std::size_t piece = std::min<std::size_t>(ob.size() - w, std::size_t(512) << 20);
                o.write(reinterpret_cast<const char*>(ob.data() + w), static_cast<std::streamsize>(piece));
                w += piece;
            }
        });
        o.close();
        if (bad || !o) throw std::runtime_error(bad ? "prime missing from the order" : "write failure (ranked copy)");
    }
    std::printf("ranked      %s (%.0f s)\n", ranked.c_str(), secs());

    // ---- 3. sort in trie pre-order (ascending ranks) ---------------------------
    sopt.descending = false;
    PathRunSource<K> src(ranked, sopt);
    std::remove(ranked.c_str());
    std::printf("sorted      %zu records, max depth %d (%.0f s)\n", src.size(), src.max_len(), secs());

    // ---- 4. oracle targets ---------------------------------------------------
    PwtTargets targets;
    {
        Sha256 sha;
        u128 sum = 0;
        u64 count = 0;
        Orc1Reader rd(orc);
        rd.for_all([&](const Orc1Record& rec) {
            std::uint8_t le[16];
            for (int i = 0; i < 16; ++i) le[i] = static_cast<std::uint8_t>(rec.n >> (8 * i));
            sha.update(le, 16);
            sum += rec.n;
            ++count;
        });
        targets.sha_nset = sha.finish();
        targets.total_check = sum;
        targets.record_count = count;
        if (count != rd.header().n_records || sum != rd.footer().total_check)
            throw std::runtime_error("oracle stream disagrees with its own header/footer (count or sum)");
        std::printf("oracle      %llu records, sha_nset %s (%.0f s)\n", (unsigned long long)count,
                    sha256_hex(targets.sha_nset).c_str(), secs());
    }
    if (targets.record_count != src.size()) throw std::runtime_error("record count mismatch: oracle vs paths");

    // ---- 5. encode -------------------------------------------------------------
    Pwt2Encoder enc(full, naive, B, naive ? 1 : ctx_set);
    Pwt2EncodeStats st;
    const std::vector<std::uint8_t> bytes = enc.encode(src, ord, targets, &st);
    std::printf("encoded     (%.0f s)\n", secs());
    src.release();
    {
        std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error(std::string("cannot write ") + out_path);
        for (std::size_t w = 0; w < bytes.size();) {
            const std::size_t piece = std::min<std::size_t>(bytes.size() - w, std::size_t(512) << 20);
            out.write(reinterpret_cast<const char*>(bytes.data() + w), static_cast<std::streamsize>(piece));
            w += piece;
        }
        out.close();
        if (!out) throw std::runtime_error("write failure");
    }
    const double N = double(src.size());
    std::printf("nodes       %llu (%.3f per record); symbols %llu (table %llu); escapes %llu\n",
                (unsigned long long)st.nodes, double(st.nodes) / N, (unsigned long long)st.symbols,
                (unsigned long long)st.table_symbols, (unsigned long long)st.escapes);
    std::printf("sections    header %zu + classes %llu + tables %llu + symbols %llu + raw %llu + footer %zu\n",
                kPwt2HeaderBytes, (unsigned long long)st.classes_bytes, (unsigned long long)st.tables_bytes,
                (unsigned long long)st.sym_bytes, (unsigned long long)st.raw_bytes, kPwt2FooterBytes);
    if (naive) {
        std::printf("bits/el     stream %.3f (of which rank table %.3f) + classes %.4f + framing %.4f\n",
                    8.0 * st.sym_bytes / N, 8.0 * st.table_naive_bytes / N, 8.0 * st.classes_bytes / N,
                    8.0 * (kPwt2HeaderBytes + kPwt2FooterBytes) / N);
    } else {
        std::printf("bits/el     symbols %.3f (quantized ideal %.3f, exact entropy %.3f) + raw %.3f + tables %.4f + classes %.4f + framing %.4f\n",
                    8.0 * st.sym_bytes / N, st.ideal_symbol_bits / N, st.entropy_symbol_bits / N, double(st.raw_bits) / N,
                    8.0 * st.tables_bytes / N, 8.0 * st.classes_bytes / N, 8.0 * (kPwt2HeaderBytes + kPwt2FooterBytes) / N);
        if (ctx_set >= 2) {
            static const char* names[kPwt2S3Classes] = {"table", "count", "count tail", "flag", "gap length", "gap mantissa"};
            std::printf("contexts    %llu models; per class: ideal bits/el | models | table bytes\n", (unsigned long long)st.contexts);
            for (int c = 0; c < (ctx_set == 3 ? kPwt2S3Classes : kPwt2S2Classes); ++c)
                std::printf("  %-11s %9.4f | %7llu | %9llu\n", names[c], st.class_ideal_bits[c] / N,
                            (unsigned long long)st.class_contexts[c], (unsigned long long)st.class_table_bytes[c]);
        }
        std::printf("rank table  %.4f bits/el (symbols ideal %.4f + raw %.4f + classes %.4f)\n",
                    (st.table_ideal_bits + double(st.table_raw_bits) + 8.0 * st.classes_bytes) / N,
                    st.table_ideal_bits / N, double(st.table_raw_bits) / N, 8.0 * st.classes_bytes / N);
    }
    std::printf("FILE        %llu bytes = %.3f bits/el   (%s)\n"
                "            paths reproduce the oracle's sum; quote only after pwt2_decode --certify exits 0\n",
                (unsigned long long)st.file_bytes, 8.0 * st.file_bytes / N, out_path);
    std::printf("elapsed     %.0f s\n", secs());
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* paths = nullptr;
    const char* orc = nullptr;
    const char* out = nullptr;
    int fields = 12;
    bool full = false, naive = false;
    int ctx_set = 1;
    bool ctx_given = false;
    u64 B = 100000000ull;
    PathSortOptions sopt;
    sopt.tag = "pwt2_encode_" + std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--full")) full = true;
        else if (!std::strcmp(argv[i], "--naive")) naive = true;
        else if (!std::strcmp(argv[i], "--ctx-set") && i + 1 < argc) { ctx_set = std::atoi(argv[++i]); ctx_given = true; }
        else if (!std::strcmp(argv[i], "--fields") && i + 1 < argc) fields = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--sieve-limit") && i + 1 < argc) B = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--scratch") && i + 1 < argc) sopt.scratch = argv[++i];
        else if (!std::strcmp(argv[i], "--chunk") && i + 1 < argc) sopt.chunk_recs = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) sopt.threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!paths) paths = argv[i];
        else if (!orc) orc = argv[i];
        else if (!out) out = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!paths || !orc || !out) {
        std::fprintf(stderr, "usage: %s <paths.bin> <table.orc1> <out.pwt2> [--fields 12|14] [--full] [--naive]"
                             " [--sieve-limit B] [--scratch DIR] [--chunk R] [--threads N]\n", argv[0]);
        return 2;
    }
    try {
        if (ctx_set < 1 || ctx_set > 3) { std::fprintf(stderr, "--ctx-set must be 1, 2 or 3\n"); return 2; }
        if (naive && ctx_given && ctx_set != 1) { std::fprintf(stderr, "--naive has no context sets (use --ctx-set 1 or omit it)\n"); return 2; }
        if (fields == 12) return run<12>(paths, orc, out, full, naive, B, sopt, ctx_set);
        if (fields == 14) return run<14>(paths, orc, out, full, naive, B, sopt, ctx_set);
        std::fprintf(stderr, "--fields must be 12 or 14\n");
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
