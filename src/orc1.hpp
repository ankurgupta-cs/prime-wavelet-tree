#pragma once
// orc1.hpp -- the ORC1 container: writer and reader for the Carmichael table
// oracle. Format spec: docs/ORACLE_FORMAT.md. Records are factor tuples only
// (n = product); the container adds a block index for random access and a
// footer binding the encoding to the source file by SHA-256.

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sha256.hpp"
#include "u128.hpp"
#include "vbyte.hpp"

namespace cn {

inline constexpr std::uint8_t kOrc1Magic[8] = {'C', 'N', 'O', 'R', 'C', '1', 0, 0};
inline constexpr u32 kOrc1Version = 1;
inline constexpr u32 kOrc1DefaultBlockSize = 4096;
inline constexpr std::size_t kOrc1HeaderBytes = 64;
inline constexpr std::size_t kOrc1IndexEntryBytes = 32;
inline constexpr std::size_t kOrc1FooterBytes = 88;
inline constexpr int kOrc1MinFactors = 3;
inline constexpr int kOrc1MaxFactors = 14;

struct Orc1Header {
    u64 n_records = 0;
    u32 block_size = kOrc1DefaultBlockSize;
    u64 index_off = 0;    // absolute file offset of the index
    u64 footer_off = 0;   // absolute file offset of the footer
};

struct Orc1IndexEntry {
    u128 first_n = 0;     // n of the block's first record
    u64 byte_off = 0;     // offset of the block from body start
    u64 rank_base = 0;    // global record index of the block's first record
};

struct Orc1Footer {
    std::array<std::uint8_t, 32> sha_src{};    // SHA-256 of the source text
    std::array<std::uint8_t, 32> sha_body{};   // SHA-256 of body || index
    u64 n_blocks = 0;
    u128 total_check = 0;                      // sum of all n mod 2^128
};

// ---- little-endian serialization helpers ---------------------------------

inline void put_le32(std::vector<std::uint8_t>& out, u32 v) {
    for (int i = 0; i < 4; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline void put_le64(std::vector<std::uint8_t>& out, u64 v) {
    for (int i = 0; i < 8; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline void put_le128(std::vector<std::uint8_t>& out, u128 v) {
    put_le64(out, static_cast<u64>(v));
    put_le64(out, static_cast<u64>(v >> 64));
}
inline u32 get_le32(const std::uint8_t* p) {
    u32 v = 0;
    for (int i = 3; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
inline u64 get_le64(const std::uint8_t* p) {
    u64 v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
inline u128 get_le128(const std::uint8_t* p) {
    return static_cast<u128>(get_le64(p)) | (static_cast<u128>(get_le64(p + 8)) << 64);
}

// ---- writer ---------------------------------------------------------------

class Orc1Writer {
public:
    // The stream must be a fresh binary file positioned at 0 and seekable
    // (the header is rewritten at finish()).
    explicit Orc1Writer(std::ostream& out, u32 block_size = kOrc1DefaultBlockSize)
        : out_(out), block_size_(block_size) {
        if (block_size == 0) throw std::invalid_argument("orc1: block_size 0");
        const std::vector<std::uint8_t> zeros(kOrc1HeaderBytes, 0);
        write_(zeros.data(), zeros.size());
    }

    // Append one record. Factors must be strictly ascending, k in [3,14],
    // n their exact product and strictly greater than the previous n.
    void add(u128 n, const u64* factors, int k) {
        if (k < kOrc1MinFactors || k > kOrc1MaxFactors)
            throw std::invalid_argument("orc1: factor count out of range");
        u128 prod = 1;
        for (int i = 0; i < k; ++i) {
            if (i > 0 && factors[i] <= factors[i - 1])
                throw std::invalid_argument("orc1: factors not strictly ascending");
            prod *= factors[i];
        }
        if (prod != n) throw std::invalid_argument("orc1: n != product of factors");
        if (n_records_ > 0 && n <= prev_n_)
            throw std::invalid_argument("orc1: n not strictly ascending");

        if (records_in_block_ == 0) {
            index_.push_back({n, body_bytes_, n_records_});
        }
        block_buf_.push_back(static_cast<std::uint8_t>(k));
        vbyte_append(block_buf_, factors[0]);
        for (int i = 1; i < k; ++i)
            vbyte_append(block_buf_, factors[i] - factors[i - 1]);

        prev_n_ = n;
        total_check_ += n;
        ++n_records_;
        if (++records_in_block_ == block_size_) flush_block_();
    }

    // Writes index, footer, and the final header. Returns the footer.
    Orc1Footer finish(const std::array<std::uint8_t, 32>& sha_src) {
        flush_block_();
        const u64 index_off = kOrc1HeaderBytes + body_bytes_;

        std::vector<std::uint8_t> ix;
        ix.reserve(index_.size() * kOrc1IndexEntryBytes);
        for (const Orc1IndexEntry& e : index_) {
            put_le128(ix, e.first_n);
            put_le64(ix, e.byte_off);
            put_le64(ix, e.rank_base);
        }
        sha_body_.update(ix.data(), ix.size());
        write_(ix.data(), ix.size());

        Orc1Footer f;
        f.sha_src = sha_src;
        f.sha_body = sha_body_.finish();
        f.n_blocks = index_.size();
        f.total_check = total_check_;

        std::vector<std::uint8_t> fb;
        fb.insert(fb.end(), f.sha_src.begin(), f.sha_src.end());
        fb.insert(fb.end(), f.sha_body.begin(), f.sha_body.end());
        put_le64(fb, f.n_blocks);
        put_le128(fb, f.total_check);
        const u64 footer_off = index_off + ix.size();
        write_(fb.data(), fb.size());

        std::vector<std::uint8_t> hb;
        hb.insert(hb.end(), kOrc1Magic, kOrc1Magic + 8);
        put_le32(hb, kOrc1Version);
        put_le32(hb, 0);                       // flags
        put_le64(hb, n_records_);
        put_le32(hb, block_size_);
        put_le32(hb, 0);                       // reserved
        put_le64(hb, index_off);
        put_le64(hb, footer_off);
        hb.resize(kOrc1HeaderBytes, 0);        // reserved tail
        out_.seekp(0);
        write_(hb.data(), hb.size());
        out_.flush();
        if (!out_) throw std::runtime_error("orc1: final flush failed");
        return f;
    }

    u64 records_written() const { return n_records_; }

private:
    void flush_block_() {
        if (block_buf_.empty()) { records_in_block_ = 0; return; }
        sha_body_.update(block_buf_.data(), block_buf_.size());
        write_(block_buf_.data(), block_buf_.size());
        body_bytes_ += block_buf_.size();
        block_buf_.clear();
        records_in_block_ = 0;
    }

    void write_(const std::uint8_t* p, std::size_t len) {
        out_.write(reinterpret_cast<const char*>(p),
                   static_cast<std::streamsize>(len));
        if (!out_) throw std::runtime_error("orc1: write failed");
    }

    std::ostream& out_;
    u32 block_size_;
    std::vector<std::uint8_t> block_buf_;
    std::vector<Orc1IndexEntry> index_;
    Sha256 sha_body_;
    u128 prev_n_ = 0;
    u128 total_check_ = 0;
    u64 n_records_ = 0;
    u64 body_bytes_ = 0;
    u32 records_in_block_ = 0;
};

// ---- reader ---------------------------------------------------------------

struct Orc1Record {
    u128 n;
    u64 factors[kOrc1MaxFactors];
    int k;
};

class Orc1Reader {
public:
    explicit Orc1Reader(const std::string& path)
        : in_(path, std::ios::binary) {
        if (!in_) throw std::runtime_error("orc1: cannot open " + path);

        in_.seekg(0, std::ios::end);
        const u64 file_size = static_cast<u64>(in_.tellg());

        std::uint8_t hb[kOrc1HeaderBytes];
        read_at_(0, hb, sizeof hb);
        if (std::memcmp(hb, kOrc1Magic, 8) != 0)
            throw std::runtime_error("orc1: bad magic");
        if (get_le32(hb + 8) != kOrc1Version)
            throw std::runtime_error("orc1: unsupported version");
        if (get_le32(hb + 12) != 0)
            throw std::runtime_error("orc1: unsupported flags");
        header_.n_records = get_le64(hb + 16);
        header_.block_size = get_le32(hb + 24);
        header_.index_off = get_le64(hb + 32);
        header_.footer_off = get_le64(hb + 40);
        // Every offset is validated against the actual file size before any
        // allocation is sized from it; the footer must end the file exactly.
        if (header_.block_size == 0 || header_.index_off < kOrc1HeaderBytes ||
            header_.footer_off < header_.index_off ||
            header_.footer_off + kOrc1FooterBytes != file_size)
            throw std::runtime_error("orc1: corrupt header");

        const u64 ix_bytes = header_.footer_off - header_.index_off;
        if (ix_bytes % kOrc1IndexEntryBytes != 0)
            throw std::runtime_error("orc1: corrupt index size");
        const u64 n_blocks = ix_bytes / kOrc1IndexEntryBytes;
        const u64 expect_blocks =
            (header_.n_records + header_.block_size - 1) / header_.block_size;
        if (n_blocks != expect_blocks)
            throw std::runtime_error("orc1: record/block count mismatch");
        std::vector<std::uint8_t> ix(ix_bytes);
        if (ix_bytes > 0) read_at_(header_.index_off, ix.data(), ix_bytes);
        index_.reserve(n_blocks);
        const u64 body_bytes = header_.index_off - kOrc1HeaderBytes;
        for (u64 i = 0; i < n_blocks; ++i) {
            const std::uint8_t* p = ix.data() + i * kOrc1IndexEntryBytes;
            Orc1IndexEntry e{get_le128(p), get_le64(p + 16), get_le64(p + 24)};
            const bool first = index_.empty();
            if ((first && (e.byte_off != 0 || e.rank_base != 0)) ||
                (!first && (e.byte_off <= index_.back().byte_off ||
                            e.rank_base <= index_.back().rank_base ||
                            e.first_n <= index_.back().first_n)) ||
                e.byte_off >= body_bytes || e.rank_base >= header_.n_records)
                throw std::runtime_error("orc1: corrupt index entry");
            index_.push_back(e);
        }

        std::uint8_t fb[kOrc1FooterBytes];
        read_at_(header_.footer_off, fb, sizeof fb);
        std::memcpy(footer_.sha_src.data(), fb, 32);
        std::memcpy(footer_.sha_body.data(), fb + 32, 32);
        footer_.n_blocks = get_le64(fb + 64);
        footer_.total_check = get_le128(fb + 72);
        if (footer_.n_blocks != n_blocks)
            throw std::runtime_error("orc1: index/footer block count mismatch");
    }

    const Orc1Header& header() const { return header_; }
    const Orc1Footer& footer() const { return footer_; }
    const std::vector<Orc1IndexEntry>& index() const { return index_; }

    // Index of the block that could contain n: the last block whose first_n
    // is <= n. Returns n_blocks if n precedes the first block's first_n.
    u64 find_block(u128 n) const {
        std::size_t lo = 0, hi = index_.size();
        while (lo < hi) {                      // first block with first_n > n
            const std::size_t mid = (lo + hi) / 2;
            if (index_[mid].first_n <= n) lo = mid + 1; else hi = mid;
        }
        return lo == 0 ? index_.size() : lo - 1;
    }

    // Decode block b (0-based), calling f(const Orc1Record&) per record.
    template <typename F>
    void for_block(u64 b, F&& f) {
        if (b >= index_.size()) throw std::out_of_range("orc1: block out of range");
        const u64 begin = kOrc1HeaderBytes + index_[b].byte_off;
        const u64 end = (b + 1 < index_.size())
                            ? kOrc1HeaderBytes + index_[b + 1].byte_off
                            : header_.index_off;
        std::vector<std::uint8_t> buf(end - begin);
        read_at_(begin, buf.data(), buf.size());

        const u64 n_here = (b + 1 < index_.size())
                               ? index_[b + 1].rank_base - index_[b].rank_base
                               : header_.n_records - index_[b].rank_base;
        const std::uint8_t* p = buf.data();
        const std::uint8_t* pe = buf.data() + buf.size();
        Orc1Record rec;
        for (u64 i = 0; i < n_here; ++i) {
            if (p >= pe) throw std::runtime_error("orc1: truncated block");
            rec.k = *p++;
            if (rec.k < kOrc1MinFactors || rec.k > kOrc1MaxFactors)
                throw std::runtime_error("orc1: bad factor count in block");
            u64 prev = 0;
            u128 prod = 1;
            for (int j = 0; j < rec.k; ++j) {
                u64 d;
                p = vbyte_decode_u64(p, pe, d);
                if (p == nullptr) throw std::runtime_error("orc1: bad vbyte");
                if ((j > 0 && d == 0) || d > UINT64_MAX - prev)
                    throw std::runtime_error("orc1: corrupt factor delta");
                rec.factors[j] = prev + d;
                if (rec.factors[j] < 2)
                    throw std::runtime_error("orc1: corrupt factor value");
                prev = rec.factors[j];
                if (prod > U128_MAX / rec.factors[j])
                    throw std::runtime_error("orc1: record product overflow");
                prod *= rec.factors[j];
            }
            rec.n = prod;
            f(static_cast<const Orc1Record&>(rec));
        }
        if (p != pe) throw std::runtime_error("orc1: trailing bytes in block");
    }

    // Decode every record in file order.
    template <typename F>
    void for_all(F&& f) {
        for (u64 b = 0; b < index_.size(); ++b) for_block(b, f);
    }

private:
    void read_at_(u64 off, std::uint8_t* dst, std::size_t len) {
        in_.seekg(static_cast<std::streamoff>(off));
        in_.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(len));
        if (in_.gcount() != static_cast<std::streamsize>(len))
            throw std::runtime_error("orc1: short read");
    }

    std::ifstream in_;
    Orc1Header header_;
    Orc1Footer footer_;
    std::vector<Orc1IndexEntry> index_;
};

} // namespace cn
