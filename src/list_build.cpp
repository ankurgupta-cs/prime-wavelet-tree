// list_build -- timed writer for the FLAT cells of the A/B/C
// table: WRITES the cell's file from a sorted
// 16-byte little-endian value file (data/n.raw -> cells 000/001, data/d.raw
// -> cells 100/101) with phase stamps; then, UNTIMED and after the total,
// reads the file back from disk, decodes it and compares every value with
// the input (prints VERIFIED).
//
// usage: list_build vbyte   <in.raw> <out.vbg>
//        list_build classes <in.raw> <out.cnd> --mod M
//                   [--targets-orc1 <table.orc1> | --targets-raw <n.raw> | --targets-self]
//
// vbyte (cells 000 = n, 100 = d): container VBG1, all integers little-endian
//   offset 0   8 B  magic "CNVBG1\0\0"
//   offset 8   8 B  u64 N, the record count
//   offset 16  8 B  u64 P, the payload byte count
//   offset 24  P B  VByte(v_0), VByte(v_1 - v_0), ..., VByte(v_{N-1} - v_{N-2})
//                   (vbyte.hpp: base-128, low 7 bits first, high bit set on
//                   every byte but the last of a value)
//   file = 24 + P bytes. No checksum: a decoder must produce exactly N
//   values, every gap after the first >= 1, consuming exactly P bytes.
//   P is flat_bench's "VByte gaps" byte count by construction (the same
//   first-absolute-then-gaps stream): 000 = 2,362,758,662; 100 = 1,220,449,369.
//
// classes (cells 001 = n at M = 10,810,800; 101 = d at M = 2310): the CND1
//   layout of docs/CND_FORMAT.md (cnd1.hpp) with modulus M: 64 B header,
//   VByte directory, per-class interpolative blobs (bic.hpp), 88 B footer
//   with sha_payload over header||dir||body plus the certification targets.
//   The TIMED writer is cnd1_encode's own sequence of steps cut at phase
//   boundaries (same primitives, order and allocations); the UNTIMED
//   verification runs the shipped cnd1_encode on the same input, reports its
//   time, and requires the written file to be byte-identical to its output.
//   At M = 2310 on d.raw with the ORC1 (or n.raw) targets the file IS
//   the CND1 file of cell 101 (sha d6607faa...); at M = 10,810,800 on n.raw it is
//   1,545,310,979 B (flat_bench's container figure).
//   Targets (footer sha_nset + total_check): --targets-orc1 streams them from
//   the oracle as cnd_encode does; --targets-raw hashes and sums a 16-byte-LE
//   n stream file; --targets-self hashes and sums the loaded input (cell
//   001: the input IS the n list); none = zero targets (same size; footer and
//   file hash then differ from a certified container).
//   The v1 reader (Cnd1Reader) pins M = 2310, so the verification decodes
//   with a relaxed reader defined below: Cnd1Reader's structural checks,
//   payload sha check and bic_decode, with M taken from the header. At
//   M = 2310 the official Cnd1Reader decode runs as well.
//
// Stamps (stamps.hpp; wall seconds; single thread):
//   load     open + zero-filled allocation + read of the input file
//   targets  classes only, if requested: sha_nset + total_check
//   prepare  vbyte: strictly-ascending input check
//            classes: class split (sorted/odd validation, r = v mod M,
//            k = (v - r) / M appended to class r) + max k -> end width
//   encode   vbyte: the whole VByte gap stream built in RAM
//            classes: interpolative blob of every class + VByte directory
//   seal     classes: header + SHA-256 over header||dir||body + footer
//   write    everything written to <out> and the stream closed
//   total    program start to closed output; the next line gives the total
//            without the targets phase

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ostream>
#include <streambuf>
#include <string>
#include <vector>

#include "bic.hpp"
#include "cnd1.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "stamps.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

using namespace cn;

namespace {

constexpr std::uint8_t kVbgMagic[8] = {'C', 'N', 'V', 'B', 'G', '1', 0, 0};
constexpr std::size_t kVbgHeaderBytes = 24;

// ostream sink appending to a byte vector (the cnd1_encode cross-check).
class VecSink : public std::streambuf {
public:
    explicit VecSink(std::vector<std::uint8_t>& v) : v_(v) {}

protected:
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        v_.insert(v_.end(), reinterpret_cast<const std::uint8_t*>(s),
                  reinterpret_cast<const std::uint8_t*>(s) + n);
        return n;
    }
    int_type overflow(int_type c) override {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
            v_.push_back(static_cast<std::uint8_t>(traits_type::to_char_type(c)));
        return traits_type::not_eof(c);
    }

private:
    std::vector<std::uint8_t>& v_;
};

void put64le(std::vector<std::uint8_t>& o, u64 v) {
    for (int i = 0; i < 8; ++i) o.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
u64 get64le(const std::uint8_t* p) {
    u64 v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}

[[noreturn]] void fail(const std::string& m) { throw std::runtime_error(m); }

// ---------------------------------------------------------------------------
// vbyte mode
// ---------------------------------------------------------------------------
int run_vbyte(PhaseClock& clk, const std::vector<u128>& vals, const std::string& out_path) {
    const std::size_t N = vals.size();
    for (std::size_t i = 1; i < N; ++i)
        if (vals[i] <= vals[i - 1]) fail("input not strictly ascending at record " + std::to_string(i));
    clk.stamp("prepare", "input check: strictly ascending");

    // values < 2^84 give gaps of at most 12 VByte bytes: reserve so the
    // stream never reallocates inside the timed phase
    const std::size_t per = (N && vals.back() >= (static_cast<u128>(1) << 84)) ? kVByteMax128 : 12;
    std::vector<std::uint8_t> payload;
    payload.reserve(N * per + kVByteMax128);
    if (N) {
        vbyte_append(payload, vals[0]);
        for (std::size_t i = 1; i < N; ++i) vbyte_append(payload, vals[i] - vals[i - 1]);
    }
    clk.stamp("encode", "VByte gap stream in RAM (first value absolute)");

    std::vector<std::uint8_t> head(kVbgMagic, kVbgMagic + 8);
    put64le(head, N);
    put64le(head, payload.size());
    {
        std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
        if (!out) fail("cannot write " + out_path);
        write_all(out, head.data(), head.size());
        write_all(out, payload.data(), payload.size());
        out.close();
        if (!out) fail("write failure on " + out_path);
    }
    clk.stamp("write", "24 B header + payload, closed");
    clk.total();
    const double Nd = N ? double(N) : 1.0;
    const u64 P = payload.size(), F = P + kVbgHeaderBytes;
    std::printf("payload    %llu B = %.3f bits/el   (flat_bench \"VByte gaps\" figure)\n",
                (unsigned long long)P, 8.0 * double(P) / Nd);
    std::printf("file       %llu B = %.3f bits/el   (payload + 24 B header)\n",
                (unsigned long long)F, 8.0 * double(F) / Nd);
    payload = std::vector<std::uint8_t>();

    // ---- UNTIMED verification: decode the file from disk ----------------------
    const double v0 = wall_now();
    const std::vector<std::uint8_t> buf = file_bytes(out_path);
    if (buf.size() < kVbgHeaderBytes || std::memcmp(buf.data(), kVbgMagic, 8) != 0) fail("verify: bad magic");
    const u64 n_hdr = get64le(buf.data() + 8), p_hdr = get64le(buf.data() + 16);
    if (n_hdr != N) fail("verify: header count mismatch");
    if (p_hdr != buf.size() - kVbgHeaderBytes) fail("verify: header payload size mismatch");
    const std::uint8_t* p = buf.data() + kVbgHeaderBytes;
    const std::uint8_t* end = buf.data() + buf.size();
    u128 cur = 0;
    for (std::size_t i = 0; i < N; ++i) {
        u128 g;
        p = vbyte_decode(p, end, g);
        if (!p) fail("verify: malformed VByte at record " + std::to_string(i));
        if (i > 0 && g == 0) fail("verify: zero gap at record " + std::to_string(i));
        if (g > U128_MAX - cur) fail("verify: value overflow at record " + std::to_string(i));
        cur += g;
        if (cur != vals[i]) fail("verify: MISMATCH at record " + std::to_string(i));
    }
    if (p != end) fail("verify: trailing payload bytes");
    const double v1 = wall_now();
    u64 brute = 0;                                // independent byte count
    for (std::size_t i = 0; i < N; ++i) brute += static_cast<u64>(vbyte_length(i ? vals[i] - vals[i - 1] : vals[0]));
    if (brute != P) fail("verify: brute-force VByte count " + std::to_string(brute) + " != payload");
    std::printf("verify     (untimed) %zu values decoded from disk in %.3f s, all equal to the input; "
                "brute-force VByte count %llu B == payload: VERIFIED\n",
                N, v1 - v0, (unsigned long long)brute);
    std::printf("sha256     %s  (untimed; the file as read back from disk)\n",
                sha256_hex(Sha256::hash(buf.data(), buf.size())).c_str());
    return 0;
}

// ---------------------------------------------------------------------------
// classes mode
// ---------------------------------------------------------------------------

// Relaxed CND reader for the verification pass: Cnd1Reader's checks with M
// read from the header instead of pinned to 2310. Returns every value,
// class by class (each class ascending).
std::vector<u128> relaxed_cnd_decode(const std::vector<std::uint8_t>& buf, u32 expect_M, u64 expect_records,
                                     std::array<std::uint8_t, 32>& sha_nset, u128& total_check) {
    namespace d = cnd1_detail;
    const u64 fsize = buf.size();
    if (fsize < kCnd1HeaderBytes + kCnd1FooterBytes) fail("cnd: file too small");
    if (buf[0] != 'C' || buf[1] != 'N' || buf[2] != 'D' || buf[3] != '1') fail("cnd: bad magic");
    if (d::get32(&buf[4]) != kCnd1Version) fail("cnd: bad version");
    const u32 M = d::get32(&buf[8]);
    const u32 end_width = d::get32(&buf[12]);
    const u64 records = d::get64(&buf[16]);
    const u64 dir_off = d::get64(&buf[24]);
    const u64 body_off = d::get64(&buf[32]);
    const u64 footer_off = d::get64(&buf[40]);
    if (M != expect_M || M == 0) fail("cnd: header M differs from --mod");
    if (records != expect_records) fail("cnd: header record count differs from the input");
    if (end_width < 1 || end_width > 64) fail("cnd: bad end_width");
    for (std::size_t i = 48; i < kCnd1HeaderBytes; ++i)
        if (buf[i] != 0) fail("cnd: nonzero reserved header byte");
    if (dir_off != kCnd1HeaderBytes || body_off < dir_off || footer_off < body_off ||
        footer_off + kCnd1FooterBytes != fsize)
        fail("cnd: bad offsets");

    std::vector<u64> counts(M, 0), blob_off(M, 0), blob_len(M, 0);
    const std::uint8_t* p = &buf[dir_off];
    const std::uint8_t* dend = buf.data() + body_off;
    u64 running = body_off, total = 0;
    for (u32 r = 0; r < M; ++r) {
        u64 c;
        p = vbyte_decode_u64(p, dend, c);
        if (!p) fail("cnd: corrupt directory");
        if (c > records - total) fail("cnd: class count exceeds records");
        if (end_width < 64 && c > (static_cast<u64>(1) << end_width)) fail("cnd: class count exceeds universe");
        counts[r] = c;
        total += c;
        if (c) {
            u64 bl;
            p = vbyte_decode_u64(p, dend, bl);
            if (!p) fail("cnd: corrupt directory");
            if (bl > footer_off - running) fail("cnd: blob overruns body");
            if (bl < (std::min<u64>(c, 2) * end_width + 7) / 8) fail("cnd: blob too short for its ends");
            blob_off[r] = running;
            blob_len[r] = bl;
            running += bl;
        }
    }
    if (p != dend) fail("cnd: directory size mismatch");
    if (running != footer_off) fail("cnd: body size mismatch");
    if (total != records) fail("cnd: record count mismatch");

    const std::uint8_t* f = &buf[footer_off];
    const auto pd = Sha256::hash(buf.data(), static_cast<std::size_t>(footer_off));
    if (std::memcmp(pd.data(), f, 32) != 0) fail("cnd: payload sha mismatch");
    std::memcpy(sha_nset.data(), f + 32, 32);
    total_check = d::get128(f + 64);
    if (d::get64(f + 80) != records) fail("cnd: footer count mismatch");

    std::vector<u128> all;
    all.reserve(static_cast<std::size_t>(records));
    for (u32 r = 0; r < M; ++r) {
        if (!counts[r]) continue;
        if ((M % 2 == 0) && (r % 2 == 0)) fail("cnd: nonempty even class (values are odd, M even)");
        const std::vector<u128> ks = bic_decode(&buf[blob_off[r]], static_cast<std::size_t>(blob_len[r]),
                                                counts[r], static_cast<int>(end_width));
        for (std::size_t i = 0; i < ks.size(); ++i) {
            if (i && ks[i] <= ks[i - 1]) fail("cnd: class not increasing");
            all.push_back(ks[i] * M + r);
        }
    }
    return all;
}

int run_classes(PhaseClock& clk, const std::vector<u128>& vals, const std::string& out_path, u32 M,
                const Cnd1FooterTargets& targets, double t_targets) {
    namespace d = cnd1_detail;
    const std::size_t N = vals.size();
    const double Nd = N ? double(N) : 1.0;
    u64 file_total = 0;
    u32 end_width = 0;
    {
        // ---- prepare: cnd1_encode's classify, validation and end width ---------
        std::vector<std::vector<u128>> ks(M);
        u128 prev = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (i && vals[i] <= prev) fail("input not sorted/unique at record " + std::to_string(i));
            if ((vals[i] & 1) == 0) fail("even value at record " + std::to_string(i));
            prev = vals[i];
            const u64 r = static_cast<u64>(vals[i] % M);
            ks[r].push_back((vals[i] - r) / M);
        }
        u128 max_k = 0;
        for (const auto& v : ks)
            if (!v.empty() && v.back() > max_k) max_k = v.back();
        end_width = static_cast<u32>(std::max(1, d::bits_of(max_k)));
        if (end_width > 64) fail("k exceeds 64 bits (values outside the CND domain for this M)");
        clk.stamp("prepare", "class split: validate, r = v mod M, k = (v - r)/M per class; end width");

        // ---- encode: blobs, then the directory (it needs their lengths) ------
        std::vector<std::vector<std::uint8_t>> blobs(M);
        for (u32 r = 0; r < M; ++r)
            if (!ks[r].empty()) bic_append(blobs[r], ks[r], static_cast<int>(end_width));
        std::vector<std::uint8_t> dir;
        for (u32 r = 0; r < M; ++r) {
            vbyte_append(dir, ks[r].size());
            if (!ks[r].empty()) vbyte_append(dir, blobs[r].size());
        }
        clk.stamp("encode", "interpolative blob per class + VByte directory");

        // ---- seal: header, payload sha, footer ---------------------------------
        const u64 dir_off = kCnd1HeaderBytes;
        const u64 body_off = dir_off + dir.size();
        u64 body_bytes = 0;
        u64 nonempty = 0;
        for (const auto& b : blobs) { body_bytes += b.size(); nonempty += b.empty() ? 0 : 1; }
        const u64 footer_off = body_off + body_bytes;
        std::vector<std::uint8_t> head;
        head.reserve(kCnd1HeaderBytes);
        head.insert(head.end(), {'C', 'N', 'D', '1'});
        d::put32(head, kCnd1Version);
        d::put32(head, M);
        d::put32(head, end_width);
        d::put64(head, N);
        d::put64(head, dir_off);
        d::put64(head, body_off);
        d::put64(head, footer_off);
        head.resize(kCnd1HeaderBytes, 0);
        Sha256 payload;
        payload.update(head.data(), head.size());
        payload.update(dir.data(), dir.size());
        for (const auto& b : blobs) payload.update(b.data(), b.size());
        std::vector<std::uint8_t> foot;
        const auto pd = payload.finish();
        foot.insert(foot.end(), pd.begin(), pd.end());
        foot.insert(foot.end(), targets.sha_nset.begin(), targets.sha_nset.end());
        d::put128(foot, targets.total_check);
        d::put64(foot, N);
        clk.stamp("seal", "header + SHA-256 over header||dir||body + footer");

        // ---- write: as cnd1_encode writes (header, dir, every blob, footer) ----
        {
            std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
            if (!out) fail("cannot write " + out_path);
            out.write(reinterpret_cast<const char*>(head.data()), static_cast<std::streamsize>(head.size()));
            out.write(reinterpret_cast<const char*>(dir.data()), static_cast<std::streamsize>(dir.size()));
            for (const auto& b : blobs)
                out.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
            out.write(reinterpret_cast<const char*>(foot.data()), static_cast<std::streamsize>(foot.size()));
            out.close();
            if (!out) fail("write failure on " + out_path);
        }
        clk.stamp("write", "header, directory, blobs, footer; closed");
        const double total = clk.total();
        std::printf("stamp  %-14s %10.3f s   (total minus targets)\n", "total-targets", total - t_targets);
        file_total = footer_off + kCnd1FooterBytes;
        std::printf("classes    %llu nonempty of %u, end width %u\n", (unsigned long long)nonempty, M, end_width);
        std::printf("container  header 64 + directory %zu + blobs %llu + footer 88 = %llu B = %.3f bits/el\n",
                    dir.size(), (unsigned long long)body_bytes, (unsigned long long)file_total,
                    8.0 * double(file_total) / Nd);
        std::printf("footer     sha_nset %s\n", sha256_hex(targets.sha_nset).c_str());
    }   // ks, blobs freed before the verification

    // ---- UNTIMED 1: byte identity with the shipped cnd1_encode ------------------
    std::vector<std::uint8_t> disk = file_bytes(out_path);
    if (disk.size() != file_total) fail("verify: file size on disk differs");
    std::printf("sha256     %s  (untimed; the file as read back from disk)\n",
                sha256_hex(Sha256::hash(disk.data(), disk.size())).c_str());
    {
        std::vector<std::uint8_t> ref;
        ref.reserve(static_cast<std::size_t>(file_total));
        VecSink sink(ref);
        std::ostream os(&sink);
        const double c0 = wall_now();
        const u64 ref_total = cnd1_encode(vals, M, targets, os);
        const double c1 = wall_now();
        const bool same = ref_total == file_total && ref == disk;
        std::printf("crosscheck (untimed) shipped cnd1_encode into RAM: %.3f s, %llu B; written file %s\n",
                    c1 - c0, (unsigned long long)ref_total,
                    same ? "BYTE-IDENTICAL to it" : "DIFFERS from it");
        if (!same) return 1;
    }

    // ---- UNTIMED 2: decode the file from disk and compare -----------------------
    {
        const double v0 = wall_now();
        std::array<std::uint8_t, 32> sha_nset{};
        u128 total_check = 0;
        std::vector<u128> all = relaxed_cnd_decode(disk, M, N, sha_nset, total_check);
        const double v1 = wall_now();
        std::sort(all.begin(), all.end());
        const double v2 = wall_now();
        const bool ok = all == vals && sha_nset == targets.sha_nset && total_check == targets.total_check;
        std::printf("verify     (untimed) relaxed CND reader, M = %u from the header: structure + payload sha OK; "
                    "classes decoded %.3f s + global sort %.3f s; %zu values %s\n",
                    M, v1 - v0, v2 - v1, all.size(), ok ? "all equal to the input: VERIFIED" : "MISMATCH");
        if (!ok) return 1;
    }
    disk = std::vector<std::uint8_t>();
    if (M == kCnd1M) {
        const double v0 = wall_now();
        const Cnd1Reader rd(out_path);
        const std::vector<u128> all = rd.decode_all();
        const bool ok = all == vals;
        std::printf("verify     (untimed) official v1 Cnd1Reader decode_all %.3f s: %s\n", wall_now() - v0,
                    ok ? "all equal to the input: VERIFIED" : "MISMATCH");
        if (!ok) return 1;
    } else {
        std::printf("verify     v1 Cnd1Reader not run: it pins M = %u (this file has M = %u)\n", kCnd1M, M);
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    PhaseClock clk;
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    auto usage = [&] {
        std::fprintf(stderr,
                     "usage: %s vbyte   <in.raw> <out.vbg>\n"
                     "       %s classes <in.raw> <out.cnd> --mod M\n"
                     "                  [--targets-orc1 <table.orc1> | --targets-raw <n.raw> | --targets-self]\n",
                     argv[0], argv[0]);
        return 2;
    };
    if (argc < 4) return usage();
    const std::string mode = argv[1];
    const std::string in_path = argv[2], out_path = argv[3];
    u64 M = 0;
    const char* t_orc1 = nullptr;
    const char* t_raw = nullptr;
    bool t_self = false;
    for (int i = 4; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--mod") && i + 1 < argc) M = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--targets-orc1") && i + 1 < argc) t_orc1 = argv[++i];
        else if (!std::strcmp(argv[i], "--targets-raw") && i + 1 < argc) t_raw = argv[++i];
        else if (!std::strcmp(argv[i], "--targets-self")) t_self = true;
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return usage(); }
    }
    const int n_targets = (t_orc1 ? 1 : 0) + (t_raw ? 1 : 0) + (t_self ? 1 : 0);
    if (mode == "vbyte") {
        if (M || n_targets) { std::fprintf(stderr, "vbyte mode takes no --mod / --targets-*\n"); return 2; }
    } else if (mode == "classes") {
        if (M < 1 || M > 0xffffffffull) { std::fprintf(stderr, "classes mode needs --mod M, 1 <= M < 2^32\n"); return 2; }
        if (n_targets > 1) { std::fprintf(stderr, "give at most one --targets-* option\n"); return 2; }
    } else {
        return usage();
    }
    try {
        std::vector<u128> vals;
        raw16_load(in_path, vals);
        clk.stamp("load", "open + zero-filled allocation + read");
        std::printf("input      %s: %zu records; mode %s", in_path.c_str(), vals.size(), mode.c_str());
        if (mode == "classes") std::printf(", M = %llu", (unsigned long long)M);
        std::printf("\n");
        if (mode == "vbyte") return run_vbyte(clk, vals, out_path);

        Cnd1FooterTargets targets;
        double t_targets = 0.0;
        if (n_targets) {
            Sha256 sha;
            u128 sum = 0;
            u64 count = 0;
            if (t_orc1) {
                Orc1Reader rd(t_orc1);
                rd.for_all([&](const Orc1Record& rec) {
                    std::uint8_t le[16];
                    for (int i = 0; i < 16; ++i) le[i] = static_cast<std::uint8_t>(rec.n >> (8 * i));
                    sha.update(le, 16);
                    sum += rec.n;
                    ++count;
                });
            } else if (t_raw) {
                std::ifstream in(t_raw, std::ios::binary | std::ios::ate);
                if (!in) fail(std::string("cannot open ") + t_raw);
                const u64 fsize = static_cast<u64>(in.tellg());
                if (fsize % 16) fail(std::string(t_raw) + ": size is not a multiple of 16");
                in.seekg(0);
                std::vector<u128> chunk(kIoChunk / 16);
                u64 left = fsize / 16;
                while (left) {
                    const std::size_t take = static_cast<std::size_t>(std::min<u64>(left, chunk.size()));
                    read_exact(in, chunk.data(), take * 16);
                    sha.update(chunk.data(), take * 16);
                    for (std::size_t i = 0; i < take; ++i) sum += chunk[i];
                    count += take;
                    left -= take;
                }
            } else {
                sha.update(vals.data(), vals.size() * 16);
                for (const u128 v : vals) sum += v;
                count = vals.size();
            }
            if (count != vals.size())
                fail("targets: " + std::to_string(count) + " n records vs " + std::to_string(vals.size()) + " input records");
            targets.sha_nset = sha.finish();
            targets.total_check = sum;
            t_targets = clk.stamp("targets", t_orc1 ? "sha_nset + total_check streamed from ORC1"
                                            : t_raw ? "sha_nset + total_check of the n.raw file"
                                                    : "sha_nset + total_check of the loaded input");
        } else {
            std::printf("targets    none: footer sha_nset/total_check zero (size unaffected; not a certified container)\n");
        }
        return run_classes(clk, vals, out_path, static_cast<u32>(M), targets, t_targets);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
