#pragma once
// pwt1.hpp -- the PWT1 container: a Prime Wavelet Tree over descending
// prime paths (value order), entropy coded with rANS. PWT2 (pwt2.hpp)
// reuses its primitives. Specification: docs/PWT_FORMAT.md.
//
// Tree: one root-to-terminal path per stored number, primes in DESCENDING
// order (largest at the root), children of a node in ascending prime
// order. Two modes: d_min (the path is the primes of d_min(n); the decoder
// reconstructs n = d * (d^-1 mod lambda(d))) and full (the path is the
// full factorization; n = product).
//
// Serialization is a pre-order walk producing three streams:
//   symbols  rANS-coded, one static model per CONTEXT (rans.hpp):
//            count(v)   context (depth(v), bits(prime(v)))   alphabet 4096, 4095 = escape
//            flag(v)    context (depth(v))                   0/1, internal nodes only
//            gaplen(u)  context (depth(u), bucket(c(v)), bits(prime(v)), universe)
//   raw bits the gap below its leading 1 bit; escaped counts as Elias gamma
//   tables   the quantized frequency tables, one per used context
// Universe of a child u of v: PRIME-INDEX ("low") iff prime(v) < B (then
// every child is below the sieve limit B) or the header flags depth(u) as
// entirely below B; otherwise ODD-INTEGER ("high"). Gap of the first child
// is measured from the start (index 0, i.e. prime 3; or the odd integer 1),
// later siblings from the previous sibling. The decoder rebuilds the prime
// bitmap to B at load (prime_bitmap.hpp); nothing about primes is stored.
//
// Layout: [header 128] [tables] [symbol blocks] [raw bits] [footer 88].
// Footer = SHA-256 of everything before it, plus the certification targets
// of the sorted n-set (its SHA-256 over 16-byte little-endian records, the
// sum of all n mod 2^128, the record count) as in CND1.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "bic.hpp"
#include "fastbits.hpp"
#include "fastnum.hpp"
#include "lambda_bucket.hpp"
#include "prime_bitmap.hpp"
#include "rans.hpp"
#include "sha256.hpp"
#include "u128.hpp"

namespace cn {

constexpr char kPwt1Magic[4] = {'P', 'W', 'T', '1'};
constexpr u16 kPwt1Version = 1;
constexpr std::size_t kPwt1HeaderBytes = 128;
constexpr std::size_t kPwt1FooterBytes = 88;
constexpr u16 kPwt1FlagFull = 1;         // path = full factorization (else d_min)
constexpr u16 kPwt1FlagDescending = 2;   // always set in version 1
constexpr u16 kPwt1FlagNaive = 4;        // PWT0: byte-aligned VByte stream, no entropy coding
constexpr char kPwt0Magic[4] = {'P', 'W', 'T', '0'};
constexpr int kPwtPrimeBitsCtx = 42;     // bits(prime) capped to 0..41
constexpr int kPwtBuckets = 4;
constexpr u32 kPwtCountAlphabet = 4096;  // symbol 4095 = escape
constexpr u32 kPwtGapAlphabet = 128;     // bit lengths 1..127
constexpr u32 kPwtBlockSymbols = 1u << 20;

inline int pwt_bits_of(u64 v) { return v ? 64 - __builtin_clzll(v) : 0; }   // bit length (was a shift loop; same values)
inline int pwt_pbits(u64 prime) { return std::min(pwt_bits_of(prime), kPwtPrimeBitsCtx - 1); }
inline int pwt_bucket(u64 c) { return c == 1 ? 0 : c <= 3 ? 1 : c <= 15 ? 2 : 3; }

// Fixed context enumeration for a tree of maximum depth K.
struct PwtContexts {
    int K = 0;
    int count_base = 0, flag_base = 0, gap_base = 0, total = 0;
    explicit PwtContexts(int k = 1) : K(k) {
        count_base = 0;
        flag_base = (K + 1) * kPwtPrimeBitsCtx;
        gap_base = flag_base + K;
        total = gap_base + K * kPwtBuckets * kPwtPrimeBitsCtx * 2;
    }
    int count_ctx(int depth, u64 prime) const { return count_base + depth * kPwtPrimeBitsCtx + pwt_pbits(prime); }
    int flag_ctx(int depth) const { return flag_base + (depth - 1); }
    int gap_ctx(int child_depth, u64 c, u64 parent_prime, bool low) const {
        return gap_base + ((((child_depth - 1) * kPwtBuckets + pwt_bucket(c)) * kPwtPrimeBitsCtx
                            + pwt_pbits(parent_prime)) * 2 + (low ? 0 : 1));
    }
    u32 alphabet(int ctx) const {
        if (ctx < flag_base) return kPwtCountAlphabet;
        if (ctx < gap_base) return 2;
        return kPwtGapAlphabet;
    }
};

struct PwtTargets {
    std::array<std::uint8_t, 32> sha_nset{};
    u128 total_check = 0;
    u64 record_count = 0;
};

struct PwtHeader {
    bool full = false;
    bool naive = false;      // PWT0 variant (magic "PWT0"): symbols are plain VBytes
    u64 n_records = 0;
    int max_depth = 0;
    u32 depth_low = 0;       // bit d set: every prime at depth d is <= sieve_limit
    u64 sieve_limit = 0;
    u64 tables_off = 0, tables_len = 0, sym_off = 0, sym_len = 0, raw_off = 0, raw_len = 0;
    u64 n_blocks = 0, total_symbols = 0, footer_off = 0;
    u32 block_symbols = 0;
    bool depth_is_low(int d) const { return d >= 1 && d <= 31 && ((depth_low >> d) & 1u); }
};

namespace pwt_detail {

inline void put_le(std::vector<std::uint8_t>& out, u64 v, int bytes) {
    for (int i = 0; i < bytes; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline u64 get_le(const std::uint8_t* p, int bytes) {
    u64 v = 0;
    for (int i = bytes - 1; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
inline void put_le128(std::vector<std::uint8_t>& out, u128 v) {
    for (int i = 0; i < 16; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
inline u128 get_le128(const std::uint8_t* p) {
    u128 v = 0;
    for (int i = 15; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}
// base-128 VByte (7 payload bits per byte, least significant group first,
// 0x80 = continuation): the PWT0 symbol code.
inline void put_vbyte(std::vector<std::uint8_t>& out, u64 v) {
    while (v >= 0x80) { out.push_back(static_cast<std::uint8_t>(v | 0x80)); v >>= 7; }
    out.push_back(static_cast<std::uint8_t>(v));
}

// n from a path (primes in any order): d_min / divisor mode reconstructs
// n = d * r*, r* = d^-1 mod lambda(d). lambda(d) = lcm(p - 1) fits 64 bits
// for almost every d, so the lcm chain and the inverse run in u64 (binary
// gcd, 32-bit Euclid below 2^32; fastnum.hpp) and fall back to u128 only
// when the lcm would overflow. CHECKED: the product
// and n must fit 128 bits, d must be invertible mod lambda(d), and lambda
// must stay below 2^127 (the signed 128-bit Euclid's range); a path that
// violates any of these describes no valid n and throws.
inline u128 path_value(const u64* primes, int len, bool full) {
    u128 prod = 1;
    for (int i = 0; i < len; ++i)
        if (__builtin_mul_overflow(prod, static_cast<u128>(primes[i]), &prod))
            throw std::runtime_error("pwt: path product overflows 128 bits");
    if (full) return prod;
    u64 lam = 1;
    int i = 0;
    for (; i < len; ++i) {
        const u64 q = primes[i] - 1;
        const u64 a = lam / gcd64_bin(lam, q);
        if (a > ~u64(0) / q) break;   // lcm would overflow u64
        lam = a * q;
    }
    u128 n;
    if (i == len) {
        const u64 r = inv_mod_lambda(static_cast<u64>(prod % lam), lam);   // throws unless gcd(d, lambda) = 1
        if (__builtin_mul_overflow(prod, static_cast<u128>(r), &n)) throw std::runtime_error("pwt: n overflows 128 bits");
        return n;
    }
    u128 lambda = lam;
    for (; i < len; ++i) {
        const u128 q = primes[i] - 1;
        const u128 a = lambda / gcd128(lambda, q);
        if (__builtin_mul_overflow(a, q, &lambda)) throw std::runtime_error("pwt: lambda overflows 128 bits");
    }
    if (lambda >> 127) throw std::runtime_error("pwt: lambda(d) >= 2^127");
    if (gcd128(prod % lambda, lambda) != 1) throw std::runtime_error("pwt: gcd(d, lambda(d)) != 1");
    if (__builtin_mul_overflow(prod, inv_mod128(prod, lambda), &n)) throw std::runtime_error("pwt: n overflows 128 bits");
    return n;
}

} // namespace pwt_detail

// ---------------------------------------------------------------------------
// Encoder. Two passes over a replayable sorted source (pwt_paths.hpp):
// pass 1 = symbol statistics + the per-node child counts in pre-order;
// pass 2 = emission. Everything is built in RAM.
// ---------------------------------------------------------------------------
struct PwtEncodeStats {
    u64 nodes = 0, symbols = 0, raw_bits = 0, escapes = 0;
    u64 tables_bytes = 0, sym_bytes = 0, raw_bytes = 0, file_bytes = 0;
    double ideal_symbol_bits = 0.0;   // sum of -log2(quantized p): what rANS should achieve
    double entropy_symbol_bits = 0.0; // order-0 entropy with exact frequencies (the cost model's figure)
    u64 low_edges = 0, high_edges = 0;
};

class PwtEncoder {
public:
    // naive = PWT0: every symbol one VByte (node: 2c + terminal; child: the
    // plain integer gap, absolute for a first child), no contexts, no
    // universes, no tables -- the byte-aligned C = 0 cell of the A/B/C table.
    PwtEncoder(bool full, u64 sieve_limit, bool verbose = true, bool naive = false)
        : full_(full), naive_(naive), B_(checked_limit_(sieve_limit)), verbose_(verbose), bitmap_(B_) {}

    template <typename Source>
    std::vector<std::uint8_t> encode(const Source& src, const PwtTargets& targets, PwtEncodeStats* stats = nullptr) {
        if (targets.record_count != src.size()) throw std::runtime_error("pwt: target record count != source size");
        n_ = src.size();
        K_ = src.max_len();
        if (K_ < 1 || K_ > kPathMaxFieldsLimit) throw std::runtime_error("pwt: bad max depth");
        ctx_ = PwtContexts(K_);
        depth_low_ = 0;
        for (int d = 1; d <= K_; ++d)
            if (src.max_prime_at_depth(d) < B_) depth_low_ |= (1u << d);

        // ---- pass 1: histograms + pre-order child counts --------------------
        hist_.assign(ctx_.total, {});
        for (int c = 0; c < ctx_.total; ++c) hist_[c].assign(ctx_.alphabet(c), 0);
        counts_.clear();
        counts_.reserve(n_ * 7 + 16);
        overflow_.clear();
        {
            struct Open { u64 prime; u32 index; u64 c; bool terminal; std::vector<u16> child_len; };
            std::vector<Open> st(K_ + 1);
            int prev_len = 0;
            std::vector<u64> prev(K_ + 1, 0);
            bool first = true;
            auto close = [&](int depth) {
                Open& o = st[depth];
                const u64 c = o.c;
                hist_[ctx_.count_ctx(depth, o.prime)][std::min<u64>(c, kPwtCountAlphabet - 1)]++;
                if (c > 0 && depth >= 1) hist_[ctx_.flag_ctx(depth)][o.terminal ? 1 : 0]++;
                const int bk = pwt_bucket(c);
                for (const u16 packed : o.child_len) {
                    const bool low = packed & 1;
                    const int len = packed >> 1;
                    hist_[ctx_.gap_ctx(depth + 1, c, o.prime, low)][len]++;
                    (void)bk;
                }
                o.child_len.clear();
                if (c >= 65535) { counts_[o.index] = 65535; overflow_[o.index] = c; }
                else counts_[o.index] = static_cast<u16>(c);
            };
            auto open = [&](int depth, u64 prime, bool terminal) {
                Open& o = st[depth];
                o.prime = prime; o.terminal = terminal; o.c = 0;
                o.index = static_cast<u32>(counts_.size());
                if (counts_.size() >= 0xFFFFFFFFull) throw std::runtime_error("pwt: too many nodes");
                counts_.push_back(0);
                o.child_len.clear();
                if (depth >= 1) {
                    Open& p = st[depth - 1];
                    const bool low = universe_low_(depth - 1, p.prime);
                    const u64 g = gap_(low, p.c == 0 ? 0 : prev[depth], prime);
                    p.child_len.push_back(static_cast<u16>((pwt_bits_of(g) << 1) | (low ? 1 : 0)));
                    ++p.c;
                }
            };
            src.for_each([&](const u64* primes, int len) {
                if (len > K_) throw std::runtime_error("pwt: path longer than max depth");
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == primes[lcp]) ++lcp;
                for (int d = prev_len; d > lcp; --d) close(d);
                for (int d = lcp + 1; d <= len; ++d) {
                    check_descent_(d, primes[d - 1], d >= 2 ? primes[d - 2] : 0);
                    open(d, primes[d - 1], d == len);
                }
                for (int d = 0; d < len; ++d) prev[d + 1] = primes[d];
                prev_len = len;
            });
            if (first) throw std::runtime_error("pwt: empty source");
            for (int d = prev_len; d >= 0; --d) close(d);
        }
        nodes_ = counts_.size();

        // ---- models ----------------------------------------------------------
        models_.assign(ctx_.total, RansContext());
        std::vector<std::uint8_t> tables;
        double ideal = 0.0, entropy = 0.0;
        for (int c = 0; c < ctx_.total; ++c) {
            models_[c].build(hist_[c].data(), hist_[c].size());
            models_[c].serialize(tables);
            double total = 0.0;
            for (const u64 h : hist_[c]) total += double(h);
            for (std::size_t i = 0; i < models_[c].syms.size(); ++i) {
                const double h = double(hist_[c][models_[c].syms[i]]);
                ideal += h * (double(kRansProbBits) - std::log2(double(models_[c].freq[i])));
                entropy += h * std::log2(total / h);
            }
        }
        hist_.clear();
        hist_.shrink_to_fit();
        if (naive_) tables.clear();   // PWT0 stores no models

        // ---- pass 2: emission -------------------------------------------------
        std::vector<std::uint8_t> sym_bytes, raw_bytes;
        BitWriter raw(raw_bytes);
        std::vector<std::pair<u16, u16>> block;
        block.reserve(kPwtBlockSymbols);
        u64 total_symbols = 0, n_blocks = 0, raw_bits = 0, escapes = 0, low_edges = 0, high_edges = 0;
        auto flush = [&]() {
            if (block.empty()) return;
            std::vector<std::uint8_t> enc;
            rans_encode_block(block, models_, enc);
            pwt_detail::put_le(sym_bytes, block.size(), 4);
            pwt_detail::put_le(sym_bytes, enc.size(), 4);
            sym_bytes.insert(sym_bytes.end(), enc.begin(), enc.end());
            block.clear();
            ++n_blocks;
        };
        auto emit = [&](int ctx, u32 sym) {
            block.emplace_back(static_cast<u16>(ctx), static_cast<u16>(sym));
            ++total_symbols;
            if (block.size() == kPwtBlockSymbols) flush();
        };
        {
            struct Open { u64 prime; u64 c; u64 seen; };
            std::vector<Open> st(K_ + 1);
            std::vector<u64> prev(K_ + 1, 0);
            int prev_len = 0;
            u32 next_index = 0;
            bool first = true;
            auto node_count = [&](u32 index) -> u64 {
                const u16 c = counts_[index];
                if (c == 65535) return overflow_.at(index);
                return c;
            };
            auto open = [&](int depth, u64 prime, bool terminal) {
                if (naive_) {   // PWT0: plain VBytes, no contexts, no universes
                    if (depth >= 1) {
                        Open& p = st[depth - 1];
                        pwt_detail::put_vbyte(sym_bytes, p.seen == 0 ? prime : prime - prev[depth]);
                        ++total_symbols;
                        ++p.seen;
                        ++high_edges;
                    }
                    Open& o = st[depth];
                    o.prime = prime;
                    o.c = node_count(next_index++);
                    o.seen = 0;
                    pwt_detail::put_vbyte(sym_bytes, 2 * o.c + (terminal ? 1 : 0));
                    ++total_symbols;
                    return;
                }
                if (depth >= 1) {   // this node's gap, in the parent's context
                    Open& p = st[depth - 1];
                    const bool low = universe_low_(depth - 1, p.prime);
                    const u64 g = gap_(low, p.seen == 0 ? 0 : prev[depth], prime);
                    const int len = pwt_bits_of(g);
                    emit(ctx_.gap_ctx(depth, p.c, p.prime, low), static_cast<u32>(len));
                    if (len > 1) raw.put(g & ((static_cast<u128>(1) << (len - 1)) - 1), len - 1);
                    raw_bits += len - 1;
                    ++p.seen;
                    if (low) ++low_edges; else ++high_edges;
                }
                Open& o = st[depth];
                o.prime = prime;
                o.c = node_count(next_index++);
                o.seen = 0;
                if (o.c >= kPwtCountAlphabet - 1) {
                    emit(ctx_.count_ctx(depth, prime), kPwtCountAlphabet - 1);
                    const u64 v = o.c - (kPwtCountAlphabet - 1) + 1;   // >= 1, Elias gamma
                    const int vb = pwt_bits_of(v);
                    raw.put(0, vb - 1);
                    raw.put(v, vb);
                    raw_bits += 2 * vb - 1;
                    ++escapes;
                } else {
                    emit(ctx_.count_ctx(depth, prime), static_cast<u32>(o.c));
                }
                if (o.c > 0 && depth >= 1) emit(ctx_.flag_ctx(depth), terminal ? 1u : 0u);
            };
            u128 value_sum = 0;   // every record is one terminal node: its n must add up to the oracle's sum
            src.for_each([&](const u64* primes, int len) {
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == primes[lcp]) ++lcp;
                for (int d = lcp + 1; d <= len; ++d) open(d, primes[d - 1], d == len);
                for (int d = 0; d < len; ++d) prev[d + 1] = primes[d];
                prev_len = len;
                value_sum += pwt_detail::path_value(primes, len, full_);
            });
            if (next_index != nodes_) throw std::runtime_error("pwt: pass 2 visited a different node count");
            if (value_sum != targets.total_check)
                throw std::runtime_error("pwt: the paths do not reproduce the oracle's total_check "
                                         "(wrong, stale or truncated paths file); nothing written");
            for (int d = 0; d <= K_; ++d) if (st[d].seen > st[d].c) throw std::runtime_error("pwt: child count mismatch");
        }
        flush();
        raw.align();
        counts_.clear();
        counts_.shrink_to_fit();
        overflow_.clear();

        // ---- assemble ----------------------------------------------------------
        PwtHeader h;
        h.full = full_; h.naive = naive_; h.n_records = n_; h.max_depth = K_; h.depth_low = depth_low_; h.sieve_limit = B_;
        h.tables_off = kPwt1HeaderBytes; h.tables_len = tables.size();
        h.sym_off = h.tables_off + h.tables_len; h.sym_len = sym_bytes.size();
        h.raw_off = h.sym_off + h.sym_len; h.raw_len = raw_bytes.size();
        h.n_blocks = n_blocks; h.total_symbols = total_symbols; h.block_symbols = naive_ ? 0 : kPwtBlockSymbols;
        h.footer_off = h.raw_off + h.raw_len;
        std::vector<std::uint8_t> out;
        out.reserve(h.footer_off + kPwt1FooterBytes);
        write_header_(out, h);
        out.insert(out.end(), tables.begin(), tables.end());
        out.insert(out.end(), sym_bytes.begin(), sym_bytes.end());
        out.insert(out.end(), raw_bytes.begin(), raw_bytes.end());
        Sha256 sha;
        sha.update(out.data(), out.size());
        const auto payload = sha.finish();
        out.insert(out.end(), payload.begin(), payload.end());
        out.insert(out.end(), targets.sha_nset.begin(), targets.sha_nset.end());
        pwt_detail::put_le128(out, targets.total_check);
        pwt_detail::put_le(out, targets.record_count, 8);
        if (stats) {
            stats->nodes = nodes_; stats->symbols = total_symbols; stats->raw_bits = raw_bits;
            stats->escapes = escapes; stats->tables_bytes = tables.size();
            stats->sym_bytes = sym_bytes.size(); stats->raw_bytes = raw_bytes.size();
            stats->file_bytes = out.size(); stats->ideal_symbol_bits = ideal;
            stats->entropy_symbol_bits = entropy;
            stats->low_edges = low_edges; stats->high_edges = high_edges;
        }
        return out;
    }

private:
    static constexpr int kPathMaxFieldsLimit = 30;

    // The same bounds the reader enforces on the header field.
    static u64 checked_limit_(u64 B) {
        if (B < 3 || B > (u64(1) << 34)) throw std::runtime_error("pwt: sieve limit must be in [3, 2^34]");
        return B;
    }

    bool universe_low_(int parent_depth, u64 parent_prime) const {
        if (parent_depth >= 1 && parent_prime < B_) return true;
        return ((depth_low_ >> (parent_depth + 1)) & 1u) != 0;
    }
    // gap of child q after previous sibling prev (0 = first child)
    u64 gap_(bool low, u64 prev, u64 q) const {
        if (low) {
            if (q >= B_) throw std::runtime_error("pwt: prime-index universe child above the sieve limit");
            const u64 iq = bitmap_.rank(q);
            return prev == 0 ? iq + 1 : iq - bitmap_.rank(prev);
        }
        if (!(q & 1)) throw std::runtime_error("pwt: even prime");
        return prev == 0 ? (q - 1) / 2 : (q - prev) / 2;
    }
    void check_descent_(int depth, u64 prime, u64 parent) const {
        if (prime < 3) throw std::runtime_error("pwt: prime below 3 on a path");
        if (depth >= 2 && prime >= parent) throw std::runtime_error("pwt: path not strictly descending");
        if (depth >= 2 && parent % prime == 1)   // p > q on a path, so only p = 1 (mod q) can occur
            throw std::runtime_error("pwt: Korselt invariant violated (p = 1 mod q)");
    }
    static void write_header_(std::vector<std::uint8_t>& out, const PwtHeader& h) {
        out.insert(out.end(), h.naive ? kPwt0Magic : kPwt1Magic, (h.naive ? kPwt0Magic : kPwt1Magic) + 4);
        pwt_detail::put_le(out, kPwt1Version, 2);
        pwt_detail::put_le(out, (h.full ? kPwt1FlagFull : 0) | kPwt1FlagDescending | (h.naive ? kPwt1FlagNaive : 0), 2);
        pwt_detail::put_le(out, h.n_records, 8);
        out.push_back(static_cast<std::uint8_t>(h.max_depth));
        out.push_back(0); out.push_back(0); out.push_back(0);
        pwt_detail::put_le(out, h.depth_low, 4);
        pwt_detail::put_le(out, h.sieve_limit, 8);
        pwt_detail::put_le(out, h.tables_off, 8); pwt_detail::put_le(out, h.tables_len, 8);
        pwt_detail::put_le(out, h.sym_off, 8);    pwt_detail::put_le(out, h.sym_len, 8);
        pwt_detail::put_le(out, h.raw_off, 8);    pwt_detail::put_le(out, h.raw_len, 8);
        pwt_detail::put_le(out, h.n_blocks, 8);
        pwt_detail::put_le(out, h.block_symbols, 4);
        pwt_detail::put_le(out, 0, 4);
        pwt_detail::put_le(out, h.total_symbols, 8);
        pwt_detail::put_le(out, h.footer_off, 8);
        while (out.size() < kPwt1HeaderBytes) out.push_back(0);
    }

    bool full_;
    bool naive_;
    u64 B_;
    bool verbose_;
    PrimeBitmap bitmap_;
    std::size_t n_ = 0, nodes_ = 0;
    int K_ = 0;
    u32 depth_low_ = 0;
    PwtContexts ctx_;
    std::vector<std::vector<u64>> hist_;
    std::vector<RansContext> models_;
    std::vector<u16> counts_;
    std::unordered_map<u32, u64> overflow_;
};

// ---------------------------------------------------------------------------
// Reader / decoder. Validates the container structurally (header, sizes,
// payload SHA-256, tables), rebuilds the prime bitmap, then walks the tree
// and hands every stored number to a callback (tree order, NOT sorted).
// Every consumed bit is checked: bad state, overrun, a child that is not a
// prime below its parent, a leftover byte, a wrong block symbol count.
// ---------------------------------------------------------------------------
class Pwt1Reader {
public:
    explicit Pwt1Reader(std::vector<std::uint8_t> bytes) : buf_(std::move(bytes)) {
        if (buf_.size() < kPwt1HeaderBytes + kPwt1FooterBytes) throw std::runtime_error("pwt: file too small");
        const std::uint8_t* p = buf_.data();
        const bool magic1 = std::memcmp(p, kPwt1Magic, 4) == 0, magic0 = std::memcmp(p, kPwt0Magic, 4) == 0;
        if (!magic1 && !magic0) throw std::runtime_error("pwt: bad magic");
        if (pwt_detail::get_le(p + 4, 2) != kPwt1Version) throw std::runtime_error("pwt: unsupported version");
        const u64 flags = pwt_detail::get_le(p + 6, 2);
        if (!(flags & kPwt1FlagDescending) || (flags & ~u64(kPwt1FlagFull | kPwt1FlagDescending | kPwt1FlagNaive)))
            throw std::runtime_error("pwt: bad flags");
        h_.full = (flags & kPwt1FlagFull) != 0;
        h_.naive = (flags & kPwt1FlagNaive) != 0;
        if (h_.naive != magic0) throw std::runtime_error("pwt: magic and naive flag disagree");
        h_.n_records = pwt_detail::get_le(p + 8, 8);
        h_.max_depth = p[16];
        h_.depth_low = static_cast<u32>(pwt_detail::get_le(p + 20, 4));
        h_.sieve_limit = pwt_detail::get_le(p + 24, 8);
        h_.tables_off = pwt_detail::get_le(p + 32, 8); h_.tables_len = pwt_detail::get_le(p + 40, 8);
        h_.sym_off = pwt_detail::get_le(p + 48, 8);    h_.sym_len = pwt_detail::get_le(p + 56, 8);
        h_.raw_off = pwt_detail::get_le(p + 64, 8);    h_.raw_len = pwt_detail::get_le(p + 72, 8);
        h_.n_blocks = pwt_detail::get_le(p + 80, 8);
        h_.block_symbols = static_cast<u32>(pwt_detail::get_le(p + 88, 4));
        h_.total_symbols = pwt_detail::get_le(p + 96, 8);
        h_.footer_off = pwt_detail::get_le(p + 104, 8);
        if (h_.max_depth < 1 || h_.max_depth > 30) throw std::runtime_error("pwt: bad max depth");
        if (h_.sieve_limit < 3 || h_.sieve_limit > (u64(1) << 34)) throw std::runtime_error("pwt: bad sieve limit");
        if (h_.naive) {
            if (h_.block_symbols != 0 || h_.n_blocks != 0 || h_.tables_len != 0 || h_.raw_len != 0)
                throw std::runtime_error("pwt: PWT0 header carries PWT1-only sections");
        } else if (h_.block_symbols != kPwtBlockSymbols) throw std::runtime_error("pwt: block_symbols must be 2^20");
        // section geometry must tile [header, footer) exactly. Every bound is
        // a SUBTRACTION from the real footer offset so no crafted length can
        // wrap a u64 addition and place a section outside the buffer.
        const u64 footer_off = buf_.size() - kPwt1FooterBytes;   // size >= header + footer, checked above
        if (h_.footer_off != footer_off || h_.tables_off != kPwt1HeaderBytes ||
            h_.tables_len > footer_off - h_.tables_off ||
            h_.sym_off != h_.tables_off + h_.tables_len ||
            h_.sym_len > footer_off - h_.sym_off ||
            h_.raw_off != h_.sym_off + h_.sym_len ||
            h_.raw_len != footer_off - h_.raw_off)
            throw std::runtime_error("pwt: section geometry does not tile the file");
        // payload hash
        Sha256 sha;
        sha.update(p, h_.footer_off);
        const auto got = sha.finish();
        const std::uint8_t* f = p + h_.footer_off;
        if (std::memcmp(got.data(), f, 32) != 0) throw std::runtime_error("pwt: payload SHA-256 mismatch");
        std::memcpy(targets_.sha_nset.data(), f + 32, 32);
        targets_.total_check = pwt_detail::get_le128(f + 64);
        targets_.record_count = pwt_detail::get_le(f + 80, 8);
        if (targets_.record_count != h_.n_records) throw std::runtime_error("pwt: footer record count != header");
        // tables
        ctx_ = PwtContexts(h_.max_depth);
        models_.assign(ctx_.total, RansContext());
        std::size_t pos = 0;
        if (!h_.naive)
            for (int c = 0; c < ctx_.total; ++c)
                pos += models_[c].deserialize(p + h_.tables_off + pos, h_.tables_len - pos, ctx_.alphabet(c));
        if (pos != h_.tables_len) throw std::runtime_error("pwt: trailing bytes in the tables section");
        // block directory
        blocks_.clear();
        if (!h_.naive) {
            std::size_t q = 0;
            u64 syms = 0;
            while (q < h_.sym_len) {
                if (q + 8 > h_.sym_len) throw std::runtime_error("pwt: truncated block header");
                const u64 ns = pwt_detail::get_le(p + h_.sym_off + q, 4);
                const u64 nb = pwt_detail::get_le(p + h_.sym_off + q + 4, 4);
                if (ns == 0 || ns > h_.block_symbols || q + 8 + nb > h_.sym_len) throw std::runtime_error("pwt: bad block");
                blocks_.push_back({h_.sym_off + q + 8, nb, ns});
                syms += ns;
                q += 8 + nb;
            }
            if (blocks_.size() != h_.n_blocks || syms != h_.total_symbols)
                throw std::runtime_error("pwt: block directory disagrees with the header");
            for (std::size_t i = 0; i + 1 < blocks_.size(); ++i)
                if (blocks_[i].nsyms != h_.block_symbols) throw std::runtime_error("pwt: short block before the last");
        }
        if (!h_.naive) bitmap_ = std::make_unique<PrimeBitmap>(h_.sieve_limit);
    }

    const PwtHeader& header() const { return h_; }
    const PwtTargets& targets() const { return targets_; }
    std::size_t file_bytes() const { return buf_.size(); }

    // Walk the whole tree; f(u128 n) once per stored number, in tree order.
    template <typename F>
    void for_each_n(F&& f) {
        Walk w(*this);
        w.run(std::forward<F>(f));
    }

private:
    struct Block { u64 off, len, nsyms; };

    class Walk {
    public:
        explicit Walk(Pwt1Reader& r) : r_(r), raw_(r.buf_.data() + r.h_.raw_off, r.h_.raw_len) {}
        template <typename F>
        void run(F&& f) {
            path_.assign(r_.h_.max_depth + 1, 0);
            if (!r_.h_.naive) next_block_();
            node_(0, 0, f);
            if (emitted_ != r_.h_.n_records) throw std::runtime_error("pwt: decoded record count != header");
            if (consumed_ != r_.h_.total_symbols) throw std::runtime_error("pwt: symbol count != header");
            if (r_.h_.naive) {   // PWT0: the VByte stream must be consumed exactly
                if (vb_pos_ != r_.h_.sym_len) throw std::runtime_error("pwt: trailing bytes in the PWT0 stream");
                return;
            }
            if (!dec_.finished()) throw std::runtime_error("pwt: last block did not end cleanly");
            if (block_ != r_.blocks_.size()) throw std::runtime_error("pwt: unread blocks");
            // raw stream: only zero padding may remain
            const std::size_t used = raw_.bit_pos();
            if ((used + 7) / 8 != r_.h_.raw_len) throw std::runtime_error("pwt: raw bit stream length mismatch");
            for (std::size_t b = used; b < r_.h_.raw_len * 8; ++b)
                if (raw_.get(1) != 0) throw std::runtime_error("pwt: nonzero raw padding");
        }

    private:
        // PWT0 symbol: one VByte from the symbol section.
        u64 vb_() {
            const std::uint8_t* s = r_.buf_.data() + r_.h_.sym_off;
            u64 v = 0;
            for (int shift = 0;; shift += 7) {
                if (vb_pos_ >= r_.h_.sym_len || shift > 63) throw std::runtime_error("pwt: PWT0 stream overrun");
                const std::uint8_t b = s[vb_pos_++];
                if (shift == 63 && b > 1) throw std::runtime_error("pwt: PWT0 VByte above 2^64");
                v |= static_cast<u64>(b & 0x7f) << shift;
                if (!(b & 0x80)) break;
            }
            ++consumed_;
            return v;
        }
        u32 sym_(int ctx) {
            if (in_block_ == cur_.nsyms) {
                if (!dec_.finished()) throw std::runtime_error("pwt: block did not end cleanly");
                next_block_();
            }
            ++in_block_; ++consumed_;
            return dec_.decode(r_.models_[ctx]);
        }
        void next_block_() {
            if (block_ >= r_.blocks_.size()) throw std::runtime_error("pwt: symbol stream exhausted");
            cur_ = r_.blocks_[block_++];
            dec_ = RansBlockDecoder(r_.buf_.data() + cur_.off, cur_.len);
            in_block_ = 0;
        }
        template <typename F>
        void node_(int depth, u64 prime, F& f) {
            u64 c;
            bool terminal;
            if (r_.h_.naive) {   // PWT0: one VByte 2c + terminal per node
                const u64 v = vb_();
                c = v >> 1;
                terminal = (v & 1) != 0;
                if (c == 0 && !terminal) throw std::runtime_error("pwt: PWT0 leaf not terminal");
            } else {
                c = sym_(r_.ctx_.count_ctx(depth, prime));
                if (c == kPwtCountAlphabet - 1) {   // escape: Elias gamma of c - 4094
                    int zeros = 0;
                    while (raw_.get(1) == 0) if (++zeros > 40) throw std::runtime_error("pwt: bad escape");
                    u64 v = 1;
                    for (int i = 0; i < zeros; ++i) v = (v << 1) | static_cast<u64>(raw_.get(1));
                    c = v + kPwtCountAlphabet - 2;
                }
                terminal = (c == 0);
                if (c > 0 && depth >= 1) terminal = sym_(r_.ctx_.flag_ctx(depth)) == 1;
            }
            if (depth == 0 && c == 0) throw std::runtime_error("pwt: empty tree");
            if (terminal) {
                if (depth == 0) throw std::runtime_error("pwt: terminal root");
                f(pwt_detail::path_value(path_.data() + 1, depth, r_.h_.full));
                ++emitted_;
            }
            if (depth + 1 > r_.h_.max_depth && c > 0) throw std::runtime_error("pwt: tree deeper than the header");
            const bool low = (depth >= 1 && prime < r_.h_.sieve_limit) || r_.h_.depth_is_low(depth + 1);
            u64 prev = 0, prev_idx = 0;
            for (u64 i = 0; i < c; ++i) {
                u64 q;
                if (r_.h_.naive) {   // PWT0: the plain gap, absolute for the first child
                    const u64 g = vb_();
                    if (g == 0 || (i > 0 && g > ~u64(0) - prev)) throw std::runtime_error("pwt: bad PWT0 gap");
                    q = (i == 0) ? g : prev + g;
                    if (!(q & 1)) throw std::runtime_error("pwt: even child");
                } else {
                    const u32 len = sym_(r_.ctx_.gap_ctx(depth + 1, c, prime, low));
                    if (len == 0 || len > 63) throw std::runtime_error("pwt: bad gap length");
                    u64 g = 1;
                    if (len > 1) g = (u64(1) << (len - 1)) | static_cast<u64>(raw_.get(len - 1));
                    if (low) {
                        const u64 idx = (i == 0) ? g - 1 : prev_idx + g;
                        if (idx >= r_.bitmap_->count()) throw std::runtime_error("pwt: prime index out of range");
                        q = (i == 0) ? r_.bitmap_->select(idx) : r_.bitmap_->advance(prev, g, idx);
                        prev_idx = idx;
                    } else {
                        q = (i == 0) ? 2 * g + 1 : prev + 2 * g;
                    }
                }
                if (depth >= 1 && q >= prime) throw std::runtime_error("pwt: child not below its parent");
                if (q < 3) throw std::runtime_error("pwt: child below 3");
                path_[depth + 1] = q;
                node_(depth + 1, q, f);
                prev = q;
            }
        }

        Pwt1Reader& r_;
        FastBitReader raw_;   // word-based (same bits as bic.hpp's BitReader)
        RansBlockDecoder dec_;
        Block cur_{};
        std::size_t block_ = 0;
        u64 in_block_ = 0, consumed_ = 0, emitted_ = 0;
        std::size_t vb_pos_ = 0;   // PWT0 read position in the symbol section
        std::vector<u64> path_;
    };

    std::vector<std::uint8_t> buf_;
    PwtHeader h_;
    PwtTargets targets_;
    PwtContexts ctx_;
    std::vector<RansContext> models_;
    std::vector<Block> blocks_;
    std::unique_ptr<PrimeBitmap> bitmap_;
};

} // namespace cn
