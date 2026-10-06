#pragma once
// cnd1.hpp -- writer and reader for the CND1 container (docs/CND_FORMAT.md):
// the sorted d_min-set partitioned by residue class mod M = 2310, each
// class's k-set interpolative-coded (bic.hpp). The footer carries the
// certification targets (sha of the ascending 16B-LE n-stream and the
// mod-2^128 sum), which the encoder takes from the certified ORC1 oracle
// and the certifying decoder must reproduce independently by factoring.

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bic.hpp"
#include "sha256.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

namespace cn {

inline constexpr u32 kCnd1Version = 1;
inline constexpr u32 kCnd1M = 2310;
inline constexpr std::size_t kCnd1HeaderBytes = 64;
inline constexpr std::size_t kCnd1FooterBytes = 88;

struct Cnd1FooterTargets {
    std::array<std::uint8_t, 32> sha_nset{};
    u128 total_check = 0;
};

namespace cnd1_detail {

inline void put32(std::vector<std::uint8_t>& o, u32 v) {
    for (int i = 0; i < 4; ++i) o.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline void put64(std::vector<std::uint8_t>& o, u64 v) {
    for (int i = 0; i < 8; ++i) o.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline void put128(std::vector<std::uint8_t>& o, u128 v) {
    for (int i = 0; i < 16; ++i) o.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline u32 get32(const std::uint8_t* p) {
    u32 v = 0;
    for (int i = 3; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
inline u64 get64(const std::uint8_t* p) {
    u64 v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
inline u128 get128(const std::uint8_t* p) {
    u128 v = 0;
    for (int i = 15; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}

inline int bits_of(u128 v) {
    int b = 0;
    while (v > 0) { v >>= 1; ++b; }
    return b;
}

} // namespace cnd1_detail

// Encode the sorted strictly-increasing d-set. Throws on any violated
// input invariant (unsorted, even d). Returns total bytes written.
inline u64 cnd1_encode(const std::vector<u128>& ds, u32 M,
                       const Cnd1FooterTargets& targets, std::ostream& out) {
    namespace d = cnd1_detail;
    // classify; validate invariants on the way
    std::vector<std::vector<u128>> ks(M);
    u128 prev = 0;
    for (std::size_t i = 0; i < ds.size(); ++i) {
        if (i && ds[i] <= prev) throw std::invalid_argument("cnd1: d-set not sorted/unique");
        if ((ds[i] & 1) == 0) throw std::invalid_argument("cnd1: even d");
        prev = ds[i];
        const u64 r = static_cast<u64>(ds[i] % M);
        ks[r].push_back((ds[i] - r) / M);
    }
    u128 max_k = 0;
    for (const auto& v : ks)
        if (!v.empty() && v.back() > max_k) max_k = v.back();
    const u32 end_width = static_cast<u32>(std::max(1, d::bits_of(max_k)));
    if (end_width > 64)
        throw std::invalid_argument("cnd1: k exceeds 64 bits (d outside domain)");

    // body blobs first (directory needs their lengths)
    std::vector<std::vector<std::uint8_t>> blobs(M);
    for (u32 r = 0; r < M; ++r)
        if (!ks[r].empty()) bic_append(blobs[r], ks[r], static_cast<int>(end_width));

    std::vector<std::uint8_t> dir;
    for (u32 r = 0; r < M; ++r) {
        vbyte_append(dir, ks[r].size());
        if (!ks[r].empty()) vbyte_append(dir, blobs[r].size());
    }

    const u64 dir_off = kCnd1HeaderBytes;
    const u64 body_off = dir_off + dir.size();
    u64 body_bytes = 0;
    for (const auto& b : blobs) body_bytes += b.size();
    const u64 footer_off = body_off + body_bytes;

    std::vector<std::uint8_t> head;
    head.reserve(kCnd1HeaderBytes);
    head.insert(head.end(), {'C', 'N', 'D', '1'});
    d::put32(head, kCnd1Version);
    d::put32(head, M);
    d::put32(head, end_width);
    d::put64(head, ds.size());
    d::put64(head, dir_off);
    d::put64(head, body_off);
    d::put64(head, footer_off);
    head.resize(kCnd1HeaderBytes, 0);

    Sha256 payload;                            // covers header || dir || body
    payload.update(head.data(), head.size());
    payload.update(dir.data(), dir.size());
    for (const auto& b : blobs) payload.update(b.data(), b.size());

    std::vector<std::uint8_t> foot;
    const auto pd = payload.finish();
    foot.insert(foot.end(), pd.begin(), pd.end());
    foot.insert(foot.end(), targets.sha_nset.begin(), targets.sha_nset.end());
    d::put128(foot, targets.total_check);
    d::put64(foot, ds.size());

    out.write(reinterpret_cast<const char*>(head.data()), static_cast<std::streamsize>(head.size()));
    out.write(reinterpret_cast<const char*>(dir.data()), static_cast<std::streamsize>(dir.size()));
    for (const auto& b : blobs)
        out.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
    out.write(reinterpret_cast<const char*>(foot.data()), static_cast<std::streamsize>(foot.size()));
    if (!out) throw std::runtime_error("cnd1: write failure");
    return footer_off + kCnd1FooterBytes;
}

class Cnd1Reader {
public:
    // Loads and validates structure and sha_payload up front (decoder
    // obligations 1, 2 and 4); per-class decode stays lazy.
    explicit Cnd1Reader(const std::string& path) {
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) throw std::runtime_error("cnd1: cannot open " + path);
        const u64 fsize = static_cast<u64>(in.tellg());
        if (fsize < kCnd1HeaderBytes + kCnd1FooterBytes)
            throw std::runtime_error("cnd1: file too small");
        in.seekg(0);
        buf_.resize(fsize);
        in.read(reinterpret_cast<char*>(buf_.data()), static_cast<std::streamsize>(fsize));
        if (static_cast<u64>(in.gcount()) != fsize) throw std::runtime_error("cnd1: short read");

        namespace d = cnd1_detail;
        if (buf_[0] != 'C' || buf_[1] != 'N' || buf_[2] != 'D' || buf_[3] != '1')
            throw std::runtime_error("cnd1: bad magic");
        if (d::get32(&buf_[4]) != kCnd1Version) throw std::runtime_error("cnd1: bad version");
        M_ = d::get32(&buf_[8]);
        end_width_ = d::get32(&buf_[12]);
        records_ = d::get64(&buf_[16]);
        const u64 dir_off = d::get64(&buf_[24]);
        body_off_ = d::get64(&buf_[32]);
        footer_off_ = d::get64(&buf_[40]);
        if (M_ != kCnd1M)                      // version 1 pins M (spec); the
            throw std::runtime_error("cnd1: bad M");   // parity check needs even M
        if (end_width_ < 1 || end_width_ > 64)
            throw std::runtime_error("cnd1: bad end_width");
        if (dir_off != kCnd1HeaderBytes || body_off_ < dir_off || footer_off_ < body_off_ ||
            footer_off_ + kCnd1FooterBytes != fsize)
            throw std::runtime_error("cnd1: bad offsets");

        // directory
        counts_.assign(M_, 0);
        blob_off_.assign(M_, 0);
        blob_len_.assign(M_, 0);
        const std::uint8_t* p = &buf_[dir_off];
        const std::uint8_t* dend = &buf_[body_off_];
        u64 running = body_off_, total = 0;
        for (u32 r = 0; r < M_; ++r) {
            u64 c;
            p = vbyte_decode_u64(p, dend, c);
            if (!p) throw std::runtime_error("cnd1: corrupt directory");
            if (c > records_ - total)          // also blocks u64 wrap of total
                throw std::runtime_error("cnd1: class count exceeds records");
            if (end_width_ < 64 && c > (static_cast<u64>(1) << end_width_))
                throw std::runtime_error("cnd1: class count exceeds universe");
            counts_[r] = c;
            total += c;
            if (c) {
                u64 bl;
                p = vbyte_decode_u64(p, dend, bl);
                if (!p) throw std::runtime_error("cnd1: corrupt directory");
                if (bl > footer_off_ - running)   // per-blob bound; blocks u64
                    throw std::runtime_error("cnd1: blob overruns body");   // wrap of running
                const u64 need = (std::min<u64>(c, 2) * end_width_ + 7) / 8;
                if (bl < need)
                    throw std::runtime_error("cnd1: blob too short for its ends");
                blob_off_[r] = running;
                blob_len_[r] = bl;
                running += bl;
            }
        }
        if (p != dend) throw std::runtime_error("cnd1: directory size mismatch");
        if (running != footer_off_) throw std::runtime_error("cnd1: body size mismatch");
        if (total != records_) throw std::runtime_error("cnd1: record count mismatch");

        // footer + payload sha over header || dir || body (obligation 4)
        const std::uint8_t* f = &buf_[footer_off_];
        Sha256 payload;
        payload.update(buf_.data(), static_cast<std::size_t>(footer_off_));
        const auto pd = payload.finish();
        for (int i = 0; i < 32; ++i)
            if (pd[static_cast<std::size_t>(i)] != f[i])
                throw std::runtime_error("cnd1: payload sha mismatch");
        for (int i = 0; i < 32; ++i) sha_nset_[static_cast<std::size_t>(i)] = f[32 + i];
        total_check_ = d::get128(f + 64);
        if (d::get64(f + 80) != records_) throw std::runtime_error("cnd1: footer count mismatch");
    }

    u32 M() const { return M_; }
    u32 end_width() const { return end_width_; }
    u64 records() const { return records_; }
    u64 class_count(u32 r) const { return counts_[r]; }
    const std::array<std::uint8_t, 32>& sha_nset() const { return sha_nset_; }
    u128 total_check() const { return total_check_; }

    // Decode class r to its d values (ascending). Obligations 3 and 5.
    std::vector<u128> decode_class(u32 r) const {
        std::vector<u128> ds;
        if (counts_[r] == 0) return ds;
        if ((r & 1) == 0) throw std::runtime_error("cnd1: nonempty even class");
        const std::vector<u128> ks =
            bic_decode(&buf_[blob_off_[r]], static_cast<std::size_t>(blob_len_[r]),
                       counts_[r], static_cast<int>(end_width_));
        ds.reserve(ks.size());
        u128 prev = 0;
        for (std::size_t i = 0; i < ks.size(); ++i) {
            if (i && ks[i] <= prev) throw std::runtime_error("cnd1: class not increasing");
            prev = ks[i];
            ds.push_back(static_cast<u128>(ks[i]) * M_ + r);
        }
        return ds;
    }

    // All d values, globally sorted ascending.
    std::vector<u128> decode_all() const {
        std::vector<u128> all;
        all.reserve(records_);
        for (u32 r = 0; r < M_; ++r) {
            const std::vector<u128> c = decode_class(r);
            all.insert(all.end(), c.begin(), c.end());
        }
        std::sort(all.begin(), all.end());
        for (std::size_t i = 1; i < all.size(); ++i)
            if (all[i] == all[i - 1]) throw std::runtime_error("cnd1: duplicate d");
        return all;
    }

private:
    std::vector<std::uint8_t> buf_;
    u32 M_ = 0, end_width_ = 0;
    u64 records_ = 0, body_off_ = 0, footer_off_ = 0;
    std::vector<u64> counts_, blob_off_, blob_len_;
    std::array<std::uint8_t, 32> sha_nset_{};
    u128 total_check_ = 0;
};

} // namespace cn
