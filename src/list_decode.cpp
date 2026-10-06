// list_decode -- decode and certify the FLAT-cell files of the A/B/C table,
// i.e. the two containers list_build writes:
//   VBG1  VByte gap list (cell 000 = the n list, cell 100 = the d_min list);
//         layout in list_build.cpp's header comment
//   CND   the CND1 layout of docs/CND_FORMAT.md (cnd1.hpp) with ANY modulus
//         M (cell 001 = the n list at M = 10,810,800; cell 101 = the d_min
//         list at M = 2310, the CND1 file carmichael-1e24-cell101-dmin-classes.cnd)
// The format is recognised by its magic. Cnd1Reader (and so cnd_decode) pins
// M = 2310; this tool carries its own relaxed reader with M from the header
// and every structural check of Cnd1Reader, plus a few more (below).
//
// usage: list_decode <in.vbg | in.cnd> [--emit <out.raw>] [--kind n|d]
//                    [--emit-n <out.raw>] [--certify <table.orc1>] [--threads N]
//
//   (no option)  decode + structural checks; prints the SHA-256 of the file
//                and of the sorted value stream (16-byte little-endian
//                records, the layout of n.raw / d.raw), then DECODED.
//   --emit F     write the sorted values as 16-byte LE records (n.raw for
//                cells 000/001, d.raw for 100/101; its SHA-256 is the
//                "values sha256" line).
//   --kind n|d   what the stored values are: n = the Carmichael numbers
//                themselves; d = one divisor per number (the d_min list),
//                n = d * (d^-1 mod lambda(d)). Needed by --emit-n/--certify.
//   --emit-n F   write the sorted n stream (kind d: reconstructed by factoring
//                every d; kind n: the values again).
//   --certify O  compare with the ORC1 oracle: kind n compares the values
//                directly; kind d first factors every d (cnd_decode
//                --certify's code path: split_u128 -> lambda -> inverse ->
//                n; --threads, default hardware_concurrency - 1), sorts the
//                n's, and requires them distinct. Then (a) the n stream equals
//                the oracle's n stream element by element, with the same
//                count, sum mod 2^128 and SHA-256 (the oracle's stream is also
//                checked against its own header count and footer
//                total_check); (b) a CND file's footer targets (sha_nset,
//                total_check) equal that stream's. The oracle file's own
//                SHA-256 is printed (compare oracle/ORACLE).
// The last line is CERTIFIED (with --certify), DECODED (without), or FAILED:
// <reason>. Exit 0 = CERTIFIED / DECODED, 1 = FAILED, 2 = usage.
//
// Checks beyond Cnd1Reader's (CND): M is bounded by the directory size before
// any allocation (every class owns at least one directory byte); every
// directory VByte must be canonical (the writer never emits padded ones);
// every value is odd (cnd1_encode rejects even values); the record count is
// at most 2^40. VBG1 (no checksum in the format): canonical VBytes, exactly N
// values, every gap after the first >= 1, no u128 overflow, the payload
// consumed exactly; at most one record per payload byte. Kind d: every d odd,
// >= 3 and < 2^95 (split_u128's domain) and squarefree (Carmichael numbers
// are); the factorisation must be fully certified (no probable primes), and
// gcd(d, lambda(d)) = 1; n must fit in 128 bits.
//
// Stamps (stamps.hpp, wall seconds): load, file-sha, structure (CND: header,
// directory, payload SHA-256), decode, sort (CND: classes are decoded one by
// one), values-sha, emit, factor + sort-n + n-sha (kind d), emit-n,
// oracle-sha + oracle (--certify). The in-memory decode of the A/B/C table
// is decode (+ sort for CND).

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "bic.hpp"
#include "cnd1.hpp"
#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "splitter.hpp"
#include "stamps.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

using namespace cn;

namespace {

constexpr std::uint8_t kVbgMagic[8] = {'C', 'N', 'V', 'B', 'G', '1', 0, 0};
constexpr std::size_t kVbgHeaderBytes = 24;
constexpr u64 kMaxRecords = u64(1) << 40;
const u128 kSplitDomain = static_cast<u128>(1) << 95;   // split_u128: "any u128 below 2^95"

[[noreturn]] void fail(const std::string& m) { throw std::runtime_error(m); }

u64 get64le(const std::uint8_t* p) {
    u64 v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}

// Strict VByte read: vbyte_decode plus canonical length (no zero padding).
const std::uint8_t* vbyte_strict(const std::uint8_t* p, const std::uint8_t* end, u128& v) {
    const std::uint8_t* q = vbyte_decode(p, end, v);
    if (!q || q - p != vbyte_length(v)) return nullptr;
    return q;
}
const std::uint8_t* vbyte_strict_u64(const std::uint8_t* p, const std::uint8_t* end, u64& v) {
    u128 w;
    const std::uint8_t* q = vbyte_strict(p, end, w);
    if (!q || w > UINT64_MAX) return nullptr;
    v = static_cast<u64>(w);
    return q;
}

// SHA-256 and sum mod 2^128 of a value list as the 16-byte LE stream.
struct StreamDigest {
    std::array<std::uint8_t, 32> sha{};
    u128 sum = 0;
};
StreamDigest digest_of(const std::vector<u128>& v) {
    Sha256 sha;
    StreamDigest d;
    std::vector<std::uint8_t> le(16 * 4096);
    std::size_t fill = 0;
    for (const u128 x : v) {
        for (int b = 0; b < 16; ++b) le[fill + static_cast<std::size_t>(b)] = static_cast<std::uint8_t>(x >> (8 * b));
        fill += 16;
        if (fill == le.size()) { sha.update(le.data(), fill); fill = 0; }
        d.sum += x;
    }
    if (fill) sha.update(le.data(), fill);
    d.sha = sha.finish();
    return d;
}

std::array<std::uint8_t, 32> file_sha(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) fail("cannot open " + path);
    Sha256 h;
    std::vector<char> b(std::size_t(1) << 24);
    while (in) {
        in.read(b.data(), static_cast<std::streamsize>(b.size()));
        const std::streamsize g = in.gcount();
        if (g > 0) h.update(reinterpret_cast<const std::uint8_t*>(b.data()), static_cast<std::size_t>(g));
    }
    if (!in.eof()) fail("read failure on " + path);
    return h.finish();
}

// ---------------------------------------------------------------------------
// VBG1
// ---------------------------------------------------------------------------
std::vector<u128> decode_vbg1(const std::vector<std::uint8_t>& buf, PhaseClock& clk) {
    if (buf.size() < kVbgHeaderBytes) fail("vbg1: file shorter than its 24-byte header");
    const u64 N = get64le(buf.data() + 8), P = get64le(buf.data() + 16);
    if (P != buf.size() - kVbgHeaderBytes)
        fail("vbg1: header payload size " + std::to_string(P) + " != file size - 24 = " +
             std::to_string(buf.size() - kVbgHeaderBytes));
    if (N > P) fail("vbg1: header count " + std::to_string(N) + " exceeds the payload bytes (one byte per value at least)");
    if (N > kMaxRecords) fail("vbg1: record count above the decoder limit 2^40");
    std::printf("structure  VBG1: %llu records, payload %llu B (header consistent with the file size)\n",
                (unsigned long long)N, (unsigned long long)P);
    clk.stamp("structure", "VBG1 header vs file size");

    std::vector<u128> v;
    v.reserve(static_cast<std::size_t>(N));
    const std::uint8_t* p = buf.data() + kVbgHeaderBytes;
    const std::uint8_t* end = buf.data() + buf.size();
    u128 cur = 0;
    for (u64 i = 0; i < N; ++i) {
        u128 g;
        p = vbyte_strict(p, end, g);
        if (!p) fail("vbg1: malformed or non-canonical VByte at record " + std::to_string(i));
        if (i > 0 && g == 0) fail("vbg1: zero gap at record " + std::to_string(i));
        if (g > U128_MAX - cur) fail("vbg1: value overflow at record " + std::to_string(i));
        cur += g;
        v.push_back(cur);
    }
    if (p != end) fail("vbg1: " + std::to_string(end - p) + " payload bytes left after " + std::to_string(N) + " values");
    clk.stamp("decode", "VByte gaps -> values (ascending by construction)");
    return v;
}

// ---------------------------------------------------------------------------
// CND layout, any M (relaxed Cnd1Reader)
// ---------------------------------------------------------------------------
struct CndInfo {
    u32 M = 0, end_width = 0;
    u64 records = 0, nonempty = 0;
    std::array<std::uint8_t, 32> sha_nset{};
    u128 total_check = 0;
};

std::vector<u128> decode_cnd(const std::vector<std::uint8_t>& buf, PhaseClock& clk, CndInfo& info) {
    namespace d = cnd1_detail;
    const u64 fsize = buf.size();
    if (fsize < kCnd1HeaderBytes + kCnd1FooterBytes) fail("cnd: file too small");
    if (d::get32(&buf[4]) != kCnd1Version) fail("cnd: bad version");
    const u32 M = d::get32(&buf[8]);
    const u32 end_width = d::get32(&buf[12]);
    const u64 records = d::get64(&buf[16]);
    const u64 dir_off = d::get64(&buf[24]);
    const u64 body_off = d::get64(&buf[32]);
    const u64 footer_off = d::get64(&buf[40]);
    for (std::size_t i = 48; i < kCnd1HeaderBytes; ++i)
        if (buf[i] != 0) fail("cnd: nonzero reserved header byte");
    if (end_width < 1 || end_width > 64) fail("cnd: bad end_width");
    if (dir_off != kCnd1HeaderBytes || body_off < dir_off || footer_off < body_off ||
        footer_off > fsize || fsize - footer_off != kCnd1FooterBytes)
        fail("cnd: bad offsets");
    if (M == 0) fail("cnd: M = 0");
    if (M > body_off - dir_off) fail("cnd: M exceeds the directory size (every class needs a directory byte)");
    if (records > kMaxRecords) fail("cnd: record count above the decoder limit 2^40");

    std::vector<u64> counts(M, 0), blob_off(M, 0), blob_len(M, 0);
    const std::uint8_t* p = &buf[dir_off];
    const std::uint8_t* dend = buf.data() + body_off;
    u64 running = body_off, total = 0, nonempty = 0;
    for (u32 r = 0; r < M; ++r) {
        u64 c;
        p = vbyte_strict_u64(p, dend, c);
        if (!p) fail("cnd: corrupt directory (count of class " + std::to_string(r) + ")");
        if (c > records - total) fail("cnd: class count exceeds records");
        if (end_width < 64 && c > (static_cast<u64>(1) << end_width)) fail("cnd: class count exceeds universe");
        counts[r] = c;
        total += c;
        if (c) {
            u64 bl;
            p = vbyte_strict_u64(p, dend, bl);
            if (!p) fail("cnd: corrupt directory (blob length of class " + std::to_string(r) + ")");
            if (bl > footer_off - running) fail("cnd: blob overruns body");
            if (bl < (std::min<u64>(c, 2) * end_width + 7) / 8) fail("cnd: blob too short for its ends");
            blob_off[r] = running;
            blob_len[r] = bl;
            running += bl;
            ++nonempty;
        }
    }
    if (p != dend) fail("cnd: directory size mismatch");
    if (running != footer_off) fail("cnd: body size mismatch");
    if (total != records) fail("cnd: record count mismatch");

    const std::uint8_t* f = &buf[footer_off];
    const auto pd = Sha256::hash(buf.data(), static_cast<std::size_t>(footer_off));
    if (std::memcmp(pd.data(), f, 32) != 0) fail("cnd: payload sha mismatch");
    std::memcpy(info.sha_nset.data(), f + 32, 32);
    info.total_check = d::get128(f + 64);
    if (d::get64(f + 80) != records) fail("cnd: footer count mismatch");
    info.M = M;
    info.end_width = end_width;
    info.records = records;
    info.nonempty = nonempty;
    std::printf("structure  CND: M = %u, end_width %u, %llu records in %llu nonempty classes; directory %llu B, "
                "body %llu B; payload SHA-256 OK\n",
                M, end_width, (unsigned long long)records, (unsigned long long)nonempty,
                (unsigned long long)(body_off - dir_off), (unsigned long long)(footer_off - body_off));
    clk.stamp("structure", "header, directory, payload SHA-256 over header||dir||body");

    std::vector<u128> all;
    all.reserve(static_cast<std::size_t>(records));
    for (u32 r = 0; r < M; ++r) {
        if (!counts[r]) continue;
        if ((M % 2 == 0) && (r % 2 == 0)) fail("cnd: nonempty even class (values are odd, M even)");
        const std::vector<u128> ks = bic_decode(&buf[blob_off[r]], static_cast<std::size_t>(blob_len[r]),
                                                counts[r], static_cast<int>(end_width));
        for (std::size_t i = 0; i < ks.size(); ++i) {
            if (i && ks[i] <= ks[i - 1]) fail("cnd: class not increasing");
            const u128 v = ks[i] * M + r;          // k < 2^64, M < 2^32: no overflow
            if ((v & 1) == 0) fail("cnd: even value in class " + std::to_string(r));
            all.push_back(v);
        }
    }
    clk.stamp("decode", "interpolative classes -> values (class by class)");
    std::sort(all.begin(), all.end());
    clk.stamp("sort", "global sort of the decoded classes");
    return all;
}

// ---------------------------------------------------------------------------
// kind d: factor every d -> n (cnd_decode --certify's code path)
// ---------------------------------------------------------------------------
std::vector<u128> reconstruct_n(const std::vector<u128>& ds, unsigned threads) {
    std::vector<u128> ns(ds.size());
    std::atomic<bool> failed{false};
    std::mutex mtx;
    std::size_t bad_at = ds.size();
    std::string err;
    auto worker = [&](std::size_t lo, std::size_t hi) {
        for (std::size_t i = lo; i < hi && !failed.load(std::memory_order_relaxed); ++i) {
            const u128 d = ds[i];
            std::string why;
            if ((d & 1) == 0 || d < 3 || d >= kSplitDomain) {
                why = "d not odd, >= 3 and < 2^95";
            } else {
                const SplitResult f = split_u128(d);
                u128 lam = 1, prod = 1;
                bool squarefree = true;
                for (std::size_t j = 0; j < f.primes.size(); ++j) {
                    const u128 p = f.primes[j];
                    if (j && p <= f.primes[j - 1]) squarefree = false;
                    lam = lcm128(lam, p - 1);
                    prod *= p;
                }
                u128 n = 0;
                if (f.uncertified != 1 || prod != d) why = "factorisation of d not fully certified";
                else if (!squarefree) why = "d is not squarefree";
                else if (gcd128(d, lam) != 1) why = "gcd(d, lambda(d)) != 1";
                else if (__builtin_mul_overflow(d, inv_mod128(d, lam), &n)) why = "n overflows 128 bits";
                else ns[i] = n;
            }
            if (!why.empty()) {
                failed.store(true);
                std::lock_guard<std::mutex> g(mtx);
                if (i < bad_at) { bad_at = i; err = why + " at d = " + to_string(d) + " (sorted index " + std::to_string(i) + ")"; }
                return;
            }
        }
    };
    std::vector<std::thread> pool;
    const std::size_t per = ds.empty() ? 0 : (ds.size() + threads - 1) / threads;
    for (unsigned t = 0; t < threads && per; ++t) {
        const std::size_t lo = static_cast<std::size_t>(t) * per;
        if (lo >= ds.size()) break;
        pool.emplace_back(worker, lo, std::min(ds.size(), lo + per));
    }
    for (auto& th : pool) th.join();
    if (failed) fail(err);
    return ns;
}

void print_span(const char* label, const std::vector<u128>& v) {
    if (v.empty()) std::printf("%-10s 0 values\n", label);
    else std::printf("%-10s %zu values, span [%s, %s], strictly ascending\n", label, v.size(),
                     to_string(v.front()).c_str(), to_string(v.back()).c_str());
}

void check_strict(const std::vector<u128>& v, const char* what) {
    for (std::size_t i = 1; i < v.size(); ++i)
        if (v[i] <= v[i - 1]) fail(std::string("duplicate ") + what + " " + to_string(v[i]) + " at sorted index " + std::to_string(i));
}

} // namespace

int main(int argc, char** argv) {
    PhaseClock clk;
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* in_path = nullptr;
    const char* emit = nullptr;
    const char* emit_n = nullptr;
    const char* orc = nullptr;
    char kind = 0;
    unsigned threads = 0;
    auto usage = [&] {
        std::fprintf(stderr,
                     "usage: %s <in.vbg | in.cnd> [--emit <out.raw>] [--kind n|d]\n"
                     "                  [--emit-n <out.raw>] [--certify <table.orc1>] [--threads N]\n",
                     argv[0]);
        return 2;
    };
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--emit") && i + 1 < argc) emit = argv[++i];
        else if (!std::strcmp(argv[i], "--emit-n") && i + 1 < argc) emit_n = argv[++i];
        else if (!std::strcmp(argv[i], "--certify") && i + 1 < argc) orc = argv[++i];
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc)
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        else if (!std::strcmp(argv[i], "--kind") && i + 1 < argc) {
            const std::string k = argv[++i];
            if (k == "n") kind = 'n';
            else if (k == "d") kind = 'd';
            else { std::fprintf(stderr, "--kind must be n or d\n"); return usage(); }
        } else if (!in_path && argv[i][0] != '-') in_path = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return usage(); }
    }
    if (!in_path) return usage();
    if ((emit_n || orc) && !kind) {
        std::fprintf(stderr, "--emit-n and --certify need --kind n (values are the Carmichael numbers) "
                             "or --kind d (values are divisors, n = d * (d^-1 mod lambda(d)))\n");
        return 2;
    }
    try {
        // ---- load + file hash ----------------------------------------------------
        std::vector<std::uint8_t> buf = file_bytes(in_path);
        clk.stamp("load", "whole file read into memory");
        const bool is_vbg = buf.size() >= 8 && std::memcmp(buf.data(), kVbgMagic, 8) == 0;
        const bool is_cnd = buf.size() >= 4 && buf[0] == 'C' && buf[1] == 'N' && buf[2] == 'D' && buf[3] == '1';
        if (!is_vbg && !is_cnd) fail("unknown format (magic is neither \"CNVBG1\\0\\0\" nor \"CND1\")");
        const u64 file_size = buf.size();
        std::printf("input      %s: %llu B, format %s\n", in_path, (unsigned long long)file_size,
                    is_vbg ? "VBG1 (VByte gaps)" : "CND (classes mod M, interpolative)");
        const auto fsha = Sha256::hash(buf.data(), buf.size());
        std::printf("file       sha256 %s\n", sha256_hex(fsha).c_str());
        clk.stamp("file-sha", "SHA-256 of the whole file");

        // ---- decode --------------------------------------------------------------
        CndInfo info;
        std::vector<u128> vals = is_vbg ? decode_vbg1(buf, clk) : decode_cnd(buf, clk, info);
        buf = std::vector<std::uint8_t>();
        check_strict(vals, "value");
        const double Nd = vals.empty() ? 1.0 : double(vals.size());
        print_span("values", vals);
        std::printf("size       %.3f bits per value (whole file)\n", 8.0 * double(file_size) / Nd);
        const StreamDigest vdig = digest_of(vals);
        std::printf("values     sha256 %s (sorted 16-byte LE stream = the --emit file); sum mod 2^128 %s\n",
                    sha256_hex(vdig.sha).c_str(), to_string(vdig.sum).c_str());
        clk.stamp("values-sha", "strictly-ascending check + SHA-256 and sum of the value stream");
        if (is_cnd)
            std::printf("footer     sha_nset %s; total_check %s\n", sha256_hex(info.sha_nset).c_str(),
                        to_string(info.total_check).c_str());

        if (emit) {
            raw16_write(emit, vals);
            std::printf("emitted    %zu values to %s\n", vals.size(), emit);
            clk.stamp("emit", "sorted values written as 16-byte LE records");
        }
        if (!kind) {
            clk.total();
            std::printf("DECODED    %zu values; structure OK (no --certify: not compared with the oracle)\n", vals.size());
            return 0;
        }

        // ---- the n stream ----------------------------------------------------------
        std::vector<u128> ns_d;                      // kind d only
        const std::vector<u128>* ns = &vals;
        StreamDigest ndig = vdig;
        if (kind == 'd' && (orc || emit_n)) {
            if (threads == 0) {
                const unsigned hc = std::thread::hardware_concurrency();
                threads = hc > 2 ? hc - 1 : 1;
            }
            threads = static_cast<unsigned>(std::min<u64>({threads, 1024, std::max<u64>(1, vals.size())}));
            ns_d = reconstruct_n(vals, threads);
            std::printf("factored   all %zu d (threads %u): n = d * (d^-1 mod lambda(d))\n", vals.size(), threads);
            clk.stamp("factor", "split_u128 every d -> lambda -> inverse -> n");
            std::sort(ns_d.begin(), ns_d.end());
            check_strict(ns_d, "reconstructed n");
            clk.stamp("sort-n", "sort + distinctness of the reconstructed n");
            ndig = digest_of(ns_d);
            ns = &ns_d;
            print_span("n", ns_d);
            std::printf("n          sha256 %s; sum mod 2^128 %s\n", sha256_hex(ndig.sha).c_str(), to_string(ndig.sum).c_str());
            clk.stamp("n-sha", "SHA-256 and sum of the n stream");
        }
        if (is_cnd && (kind == 'n' || ns == &ns_d)) {
            const bool sha_ok = ndig.sha == info.sha_nset, sum_ok = ndig.sum == info.total_check;
            std::printf("footer     sha_nset %s; total_check %s (vs the %s n stream)\n", sha_ok ? "MATCH" : "MISMATCH",
                        sum_ok ? "MATCH" : "MISMATCH", kind == 'n' ? "stored" : "reconstructed");
            if (!sha_ok || !sum_ok) fail("CND footer targets differ from the n stream");
        }
        if (emit_n) {
            raw16_write(emit_n, *ns);
            std::printf("emitted    %zu n to %s\n", ns->size(), emit_n);
            clk.stamp("emit-n", "sorted n written as 16-byte LE records");
        }
        if (!orc) {
            clk.total();
            std::printf("DECODED    %zu values; structure OK (no --certify: not compared with the oracle)\n", vals.size());
            return 0;
        }

        // ---- certification against ORC1 --------------------------------------------
        const auto osha = file_sha(orc);
        std::printf("oracle     %s: file sha256 %s\n", orc, sha256_hex(osha).c_str());
        clk.stamp("oracle-sha", "SHA-256 of the oracle file");
        Orc1Reader ord(orc);
        Sha256 sha;
        u128 osum = 0;
        u64 count = 0, mismatches = 0;
        std::uint8_t le[16];
        ord.for_all([&](const Orc1Record& rec) {
            for (int b = 0; b < 16; ++b) le[b] = static_cast<std::uint8_t>(rec.n >> (8 * b));
            sha.update(le, 16);
            osum += rec.n;
            if (count >= ns->size() || (*ns)[count] != rec.n) {
                if (mismatches < 5)
                    std::fprintf(stderr, "MISMATCH at rank %llu: decoded %s, oracle %s\n", (unsigned long long)count,
                                 count < ns->size() ? to_string((*ns)[count]).c_str() : "(none)", to_string(rec.n).c_str());
                ++mismatches;
            }
            ++count;
        });
        const auto oseq = sha.finish();
        if (count != ord.header().n_records || osum != ord.footer().total_check)
            fail("the oracle's n stream disagrees with its own header count / footer total_check");
        std::printf("oracle     %llu records; n stream sha256 %s; sum %s\n", (unsigned long long)count,
                    sha256_hex(oseq).c_str(), to_string(osum).c_str());
        clk.stamp("oracle", "oracle n stream compared element by element");
        if (count != ns->size())
            fail("count: oracle " + std::to_string(count) + " vs decoded " + std::to_string(ns->size()));
        if (mismatches) fail(std::to_string(mismatches) + " elementwise mismatches with the oracle");
        if (oseq != ndig.sha || osum != ndig.sum) fail("n stream SHA-256 / sum differ from the oracle's");
        clk.total();
        std::printf("CERTIFIED  %llu %s; every n equals the oracle's, in order; count, sum and sha256 %s match%s\n",
                    (unsigned long long)count, kind == 'n' ? "stored n" : "d -> n (factored)",
                    sha256_hex(oseq).c_str(), is_cnd ? "; CND footer targets match" : "");
        return 0;
    } catch (const std::exception& e) {
        std::printf("FAILED     %s\n", e.what());
        std::fprintf(stderr, "FAILED: %s\n", e.what());
        return 1;
    }
}
