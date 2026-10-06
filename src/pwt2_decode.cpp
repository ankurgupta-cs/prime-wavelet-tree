// pwt2_decode -- decode a PWT2 container (docs/PWT2_FORMAT.md) and certify
// it against the ORC1 oracle, exactly as pwt_decode does for PWT1:
//   decode every stored number (tree order) -> sort -> SHA-256 of the
//   16-byte little-endian stream == footer sha_nset -> sum == total_check
//   -> element by element against the oracle's n stream.
// Time stamps: file read, container checks (payload SHA-256, tables, sieve
// to the sieve limit), rank table, tree, sort; the in-memory decode time of
// the A/B/C table is (tree - loaded).
//
// usage: pwt2_decode <in.pwt2> [--certify <table.orc1>] [--emit-n <out.raw>]
//                    [--emit-paths <out.bin> [--fields 12|14]]
//
// --emit-paths: also write every stored path in the record
// format pwt2_encode reads (pwt_paths.hpp / trie_paths: K five-byte
// big-endian primes per record, ascending, zero padded), one record per
// stored number in TREE order (not trie_paths' n order; pwt2_encode accepts
// any record order). K = --fields, default 14 for a full-factorization file
// and 12 for a divisor file (the d_min inputs used 12; the d* input
// used 14; either re-encodes to the same container when the
// matching --fields is passed to pwt2_encode). Fails if the file's max depth
// exceeds K. The file is written as <out.bin>.partial and renamed only after
// the decode's sha_nset / total_check / count checks pass. With
// --emit-paths the "tree" phase includes writing the paths.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "orc1.hpp"
#include "pwt2.hpp"
#include "pwt_paths.hpp"

using namespace cn;

namespace {

// Streams decoded paths to disk in the pwt2_encode input format.
class PathWriter {
public:
    PathWriter(const std::string& path, int K) : out_(path, std::ios::binary | std::ios::trunc), K_(K) {
        if (!out_) throw std::runtime_error("cannot write " + path);
        buf_.reserve((std::size_t(1) << 20) + 128);
    }
    void add(const u64* primes, int len) {
        if (len < 1 || len > K_ || len > kPathMaxFields) throw std::runtime_error("path length outside [1, fields]");
        u64 p[kPathMaxFields];
        for (int i = 0; i < len; ++i) {      // insertion sort into ascending order (len <= 14)
            const u64 v = primes[i];
            int j = i;
            for (; j > 0 && p[j - 1] > v; --j) p[j] = p[j - 1];
            p[j] = v;
        }
        for (int i = 0; i < K_; ++i) {
            const u64 v = i < len ? p[i] : 0;
            if (i < len && ((i > 0 && v <= p[i - 1]) || v < 3 || v >= (u64(1) << 40)))
                throw std::runtime_error("decoded path is not strictly ascending primes in [3, 2^40)");
            for (int q = 0; q < kPathFieldBytes; ++q)
                buf_.push_back(static_cast<std::uint8_t>(v >> (8 * (kPathFieldBytes - 1 - q))));
        }
        ++records_;
        if (buf_.size() >= (std::size_t(1) << 20)) flush_();
    }
    void close() {
        flush_();
        out_.close();
        if (!out_) throw std::runtime_error("write failure (paths)");
    }
    u64 records() const { return records_; }
    u64 bytes() const { return bytes_; }
    std::array<std::uint8_t, 32> sha() { return sha_.finish(); }

private:
    void flush_() {
        if (buf_.empty()) return;
        sha_.update(buf_.data(), buf_.size());
        out_.write(reinterpret_cast<const char*>(buf_.data()), static_cast<std::streamsize>(buf_.size()));
        if (!out_) throw std::runtime_error("write failure (paths)");
        bytes_ += buf_.size();
        buf_.clear();
    }
    std::ofstream out_;
    int K_;
    std::vector<std::uint8_t> buf_;
    Sha256 sha_;
    u64 records_ = 0, bytes_ = 0;
};

std::string g_partial;   // removed if the run fails after creating it

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* in_path = nullptr;
    const char* orc = nullptr;
    const char* emit = nullptr;
    const char* emit_paths = nullptr;
    int fields = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--certify") && i + 1 < argc) orc = argv[++i];
        else if (!std::strcmp(argv[i], "--emit-n") && i + 1 < argc) emit = argv[++i];
        else if (!std::strcmp(argv[i], "--emit-paths") && i + 1 < argc) emit_paths = argv[++i];
        else if (!std::strcmp(argv[i], "--fields") && i + 1 < argc) fields = std::atoi(argv[++i]);
        else if (!in_path) in_path = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return 2; }
    }
    if (!in_path) {
        std::fprintf(stderr, "usage: %s <in.pwt2> [--certify <table.orc1>] [--emit-n <out.raw>]"
                             " [--emit-paths <out.bin> [--fields 12|14]]\n", argv[0]);
        return 2;
    }
    if (fields && !emit_paths) { std::fprintf(stderr, "--fields applies to --emit-paths only\n"); return 2; }
    if (fields && fields != 12 && fields != 14) { std::fprintf(stderr, "--fields must be 12 or 14 (pwt2_encode's choices)\n"); return 2; }
    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    try {
        std::vector<std::uint8_t> bytes;
        {
            std::ifstream in(in_path, std::ios::binary | std::ios::ate);
            if (!in) { std::fprintf(stderr, "cannot open %s\n", in_path); return 2; }
            const std::streamsize sz = in.tellg();
            in.seekg(0);
            bytes.resize(static_cast<std::size_t>(sz));
            std::size_t got = 0;
            while (got < bytes.size()) {
                const std::size_t piece = std::min<std::size_t>(bytes.size() - got, std::size_t(512) << 20);
                in.read(reinterpret_cast<char*>(bytes.data() + got), static_cast<std::streamsize>(piece));
                if (in.gcount() != static_cast<std::streamsize>(piece)) { std::fprintf(stderr, "short read\n"); return 1; }
                got += piece;
            }
        }
        const double t_read = secs();
        Pwt2Reader rd(std::move(bytes));
        const Pwt2Header& h = rd.header();
        const double t_loaded = secs();
        std::printf("container   %zu bytes, %llu records, %s mode, %s, max depth %d, m %llu primes in %llu classes, context set %d\n",
                    rd.file_bytes(), (unsigned long long)h.n_records, h.full ? "full" : "divisor",
                    h.naive ? "byte-aligned" : "entropy-coded", h.max_depth, (unsigned long long)h.m,
                    (unsigned long long)h.n_classes, h.ctx_set);
        std::printf("structure   OK (payload SHA-256, tables, %llu blocks, %llu symbols); read %.1f s, loaded %.1f s\n",
                    (unsigned long long)h.n_blocks, (unsigned long long)h.total_symbols, t_read, t_loaded);

        std::vector<u128> ns;
        // exact for any real file (every record costs well over a bit); a crafted header can
        // ask for at most 8 records per file byte, the vector grows beyond if needed
        ns.reserve(std::min<u64>(h.n_records, 8 * static_cast<u64>(rd.file_bytes())));
        u128 sum = 0;
        double t_table = 0.0;
        std::unique_ptr<PathWriter> pw;
        int K = 0;
        if (emit_paths) {
            K = fields ? fields : (h.full ? 14 : 12);
            if (h.max_depth > K)
                throw std::runtime_error("max depth " + std::to_string(h.max_depth) + " exceeds " + std::to_string(K) +
                                         " path fields (pass --fields 14)");
            g_partial = std::string(emit_paths) + ".partial";
            pw = std::make_unique<PathWriter>(g_partial, K);
        }
        if (!pw)
            rd.for_each_n([&](u128 n) { ns.push_back(n); sum += n; }, [&] { t_table = secs(); });
        else {
            rd.for_each_path([&](u128 n, const u64* pr, int len) { ns.push_back(n); sum += n; pw->add(pr, len); },
                             [&] { t_table = secs(); });
            pw->close();
        }
        const double t_tree = secs();
        std::printf("rank table  %llu primes decoded (%.1f s cumulative)\n", (unsigned long long)h.m, t_table);
        std::printf("decoded     %zu numbers in tree order (%.1f s cumulative)\n", ns.size(), t_tree);
        std::printf("phases      read %.1f s | checks + sieve %.1f s | rank table %.1f s | tree %.1f s | in-memory decode %.1f s (table + tree)\n",
                    t_read, t_loaded - t_read, t_table - t_loaded, t_tree - t_table, t_tree - t_loaded);
        std::sort(ns.begin(), ns.end());
        for (std::size_t i = 1; i < ns.size(); ++i)
            if (ns[i] <= ns[i - 1]) { std::fprintf(stderr, "DUPLICATE decoded value at %zu\n", i); return 1; }
        Sha256 sha;
        std::vector<std::uint8_t> le(16);
        for (const u128 n : ns) {
            for (int i = 0; i < 16; ++i) le[i] = static_cast<std::uint8_t>(n >> (8 * i));
            sha.update(le.data(), 16);
        }
        const auto got = sha.finish();
        const bool sha_ok = std::memcmp(got.data(), rd.targets().sha_nset.data(), 32) == 0;
        const bool sum_ok = sum == rd.targets().total_check;
        std::printf("sorted      (%.1f s); sha_nset %s; total_check %s; count %s\n", secs(),
                    sha_ok ? "MATCH" : "MISMATCH", sum_ok ? "MATCH" : "MISMATCH",
                    ns.size() == rd.targets().record_count ? "MATCH" : "MISMATCH");
        if (!sha_ok || !sum_ok || ns.size() != rd.targets().record_count) {
            if (!g_partial.empty()) std::remove(g_partial.c_str());
            return 1;
        }

        if (pw) {
            if (pw->records() != ns.size()) throw std::runtime_error("path record count != decoded count");
            std::remove(emit_paths);   // rename does not replace on Windows
            if (std::rename(g_partial.c_str(), emit_paths) != 0)
                throw std::runtime_error(std::string("cannot rename ") + g_partial + " to " + emit_paths);
            g_partial.clear();
            std::printf("paths       %llu records x %d fields = %llu bytes to %s (tree order; tree phase includes this write); sha256 %s\n",
                        (unsigned long long)pw->records(), K, (unsigned long long)pw->bytes(), emit_paths,
                        sha256_hex(pw->sha()).c_str());
        }

        if (emit) {
            std::ofstream out(emit, std::ios::binary | std::ios::trunc);
            if (!out) { std::fprintf(stderr, "cannot write %s\n", emit); return 2; }
            std::vector<std::uint8_t> buf;
            buf.reserve(1u << 20);
            for (const u128 n : ns) {
                for (int i = 0; i < 16; ++i) buf.push_back(static_cast<std::uint8_t>(n >> (8 * i)));
                if (buf.size() >= (1u << 20)) { out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size())); buf.clear(); }
            }
            out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
            out.close();
            if (!out) { std::fprintf(stderr, "write failure\n"); return 1; }
            std::printf("emitted     %zu records to %s\n", ns.size(), emit);
        }
        if (orc) {
            {   // the oracle copy itself: its file hash goes into the log (compare oracle/ORACLE)
                std::ifstream oin(orc, std::ios::binary);
                if (!oin) { std::fprintf(stderr, "cannot open %s\n", orc); return 2; }
                Sha256 fh;
                std::vector<char> ob(std::size_t(1) << 24);
                while (oin) {
                    oin.read(ob.data(), static_cast<std::streamsize>(ob.size()));
                    const std::streamsize g = oin.gcount();
                    if (g > 0) fh.update(reinterpret_cast<const std::uint8_t*>(ob.data()), static_cast<std::size_t>(g));
                }
                std::printf("oracle file SHA-256 %s (%.1f s)\n", sha256_hex(fh.finish()).c_str(), secs());
            }
            Orc1Reader ord(orc);
            std::size_t i = 0;
            u64 mismatches = 0;
            ord.for_all([&](const Orc1Record& rec) {
                if (i >= ns.size() || ns[i] != rec.n) {
                    if (mismatches < 5)
                        std::fprintf(stderr, "MISMATCH at %zu: decoded %s, oracle %s\n", i,
                                     i < ns.size() ? to_string(ns[i]).c_str() : "(none)", to_string(rec.n).c_str());
                    ++mismatches;
                }
                ++i;
            });
            if (i != ns.size()) { std::fprintf(stderr, "count mismatch: oracle %zu vs decoded %zu\n", i, ns.size()); return 1; }
            if (mismatches) { std::fprintf(stderr, "CERTIFICATION FAILED: %llu mismatches\n", (unsigned long long)mismatches); return 1; }
            std::printf("CERTIFIED   every decoded number equals the oracle's, in order (%.1f s); oracle sha_src %s\n",
                        secs(), sha256_hex(ord.footer().sha_src).c_str());
        }
        std::printf("elapsed     %.1f s\n", secs());
        return 0;
    } catch (const std::exception& e) {
        if (!g_partial.empty()) std::remove(g_partial.c_str());
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
