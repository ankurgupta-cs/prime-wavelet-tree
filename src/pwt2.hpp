#pragma once
// pwt2.hpp -- the PWT2 container: the FREQUENCY-ORDERED Prime Wavelet Tree.
// Specification: docs/PWT2_FORMAT.md.
//
// Global prime order: c(p) = number of stored paths containing p; ranks
// 1..m by c descending, ties to the LARGER prime. A path is its primes
// sorted by ascending rank (most frequent first), so a node near the root
// covers every stored number that contains a frequent prime -- the
// frequency-shaped wavelet tree of Grossi-Gupta-Vitter (SODA 2004), and
// the FP-tree order of Han-Pei-Yin (2000), applied to prime sets.
//
// Serialization (entropy-coded variant):
//   rank table  the primes in rank order: runs of equal count ("classes",
//               sizes in their own section), each class coded ASCENDING as
//               gaps in a prime universe (prime index below the sieve
//               limit B, wheel-210 index above), len+raw with rANS
//   tree        pre-order walk: count(v), flag(v), then for each child the
//               rank gap from the parent (first child) or the previous
//               sibling, len+raw with rANS
// Contexts (context set 1):
//   count(v)   (depth(v), bits(value index of prime(v) among used primes))
//   flag(v)    (depth(v))
//   gaplen(u)  (depth(u), bucket(c(v)), bits(rank(v)))
//   tablelen   (bits(class size), previous member above B?)
// Byte-aligned variant (flag bit 2, the C = 0 cells 010/110): the same walk
// with every symbol one VByte (node 2c + terminal; child rank gap) and the
// table as VByte prime gaps inside each class.
//
// Layout: [header 160] [class sizes] [tables] [symbol blocks] [raw bits]
//         [footer 88: sha_payload, sha_nset, total_check, record_count]

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "pwt1.hpp"   // pwt_detail helpers, PwtTargets, PrimeBitmap, rANS, SHA-256

namespace cn {

constexpr char kPwt2Magic[4] = {'P', 'W', 'T', '2'};
constexpr u16 kPwt2Version = 1;
constexpr std::size_t kPwt2HeaderBytes = 160;
constexpr std::size_t kPwt2FooterBytes = 88;
constexpr u16 kPwt2FlagFull = 1;      // path = full factorization (else a divisor d with d*lambda(d) > n)
constexpr u16 kPwt2FlagNaive = 4;     // byte-aligned VByte streams, no entropy coding
constexpr u16 kPwt2FlagFreq = 8;      // frequency order (always set in version 1)
constexpr u32 kPwt2Wheel = 210;
constexpr u32 kPwt2WheelRes = 48;     // residues coprime to 210
constexpr u32 kPwt2LenAlphabet = 64;  // bit lengths 1..63
constexpr int kPwt2CtxSet = 1;
constexpr int kPwt2MaxDepth = 30;

// ---------------------------------------------------------------------------
// The global order and the table's prime universe.
// ---------------------------------------------------------------------------
struct Pwt2Order {
    std::vector<u64> prime;        // prime[r], r = 1..m (prime[0] = 0)
    std::vector<u64> class_size;   // runs of equal path count, in rank order (sum = m)
    std::vector<u32> vidx;         // vidx[r] = 1-based index of prime[r] in the ascending list of used primes
    u64 m() const { return prime.empty() ? 0 : prime.size() - 1; }

    // Validate (each class strictly DESCENDING by value, primes distinct,
    // odd, >= 3) and derive vidx. Throws on any violation.
    void finish() {
        const u64 M = m();
        if (M == 0) throw std::runtime_error("pwt2: empty rank table");
        u64 tot = 0;
        for (const u64 s : class_size) { if (s == 0) throw std::runtime_error("pwt2: empty class"); tot += s; }
        if (tot != M) throw std::runtime_error("pwt2: class sizes do not sum to m");
        u64 r = 1;
        for (const u64 s : class_size) {
            for (u64 i = 1; i < s; ++i)
                if (prime[r + i] >= prime[r + i - 1]) throw std::runtime_error("pwt2: class not strictly descending");
            r += s;
        }
        std::vector<std::pair<u64, u32>> byv(M);
        for (u64 q = 1; q <= M; ++q) {
            if (prime[q] < 3 || !(prime[q] & 1)) throw std::runtime_error("pwt2: rank table holds a value that is not an odd prime >= 3");
            byv[q - 1] = {prime[q], static_cast<u32>(q)};
        }
        std::sort(byv.begin(), byv.end());
        vidx.assign(M + 1, 0);
        for (u64 i = 0; i < M; ++i) {
            if (i > 0 && byv[i].first == byv[i - 1].first) throw std::runtime_error("pwt2: duplicate prime in the rank table");
            vidx[byv[i].second] = static_cast<u32>(i + 1);
        }
    }
};

// Order from (prime, path count) pairs, any order: count descending, ties
// to the larger prime.
inline Pwt2Order pwt2_order_from_counts(std::vector<std::pair<u64, u64>> pc) {
    std::sort(pc.begin(), pc.end(), [](const auto& a, const auto& b) {
        return a.second != b.second ? a.second > b.second : a.first > b.first;
    });
    Pwt2Order o;
    o.prime.assign(pc.size() + 1, 0);
    for (std::size_t i = 0; i < pc.size(); ++i) {
        if (pc[i].second == 0) throw std::runtime_error("pwt2: prime with count 0");
        o.prime[i + 1] = pc[i].first;
        if (i == 0 || pc[i].second != pc[i - 1].second) o.class_size.push_back(1);
        else ++o.class_size.back();
    }
    o.finish();
    return o;
}

struct Pwt2Wheel {
    int ridx[kPwt2Wheel];
    u32 res[kPwt2WheelRes];
    Pwt2Wheel() {
        u32 k = 0;
        for (u32 r = 0; r < kPwt2Wheel; ++r) {
            u32 a = r, b = kPwt2Wheel;
            while (b) { const u32 t = a % b; a = b; b = t; }
            if (a == 1) { ridx[r] = static_cast<int>(k); res[k++] = r; }
            else ridx[r] = -1;
        }
    }
    u64 index(u64 x) const {
        const int r = ridx[x % kPwt2Wheel];
        if (r < 0) throw std::runtime_error("pwt2: value above the sieve limit shares a factor with 210");
        return kPwt2WheelRes * (x / kPwt2Wheel) + static_cast<u64>(r);
    }
    u64 value(u64 w) const { return kPwt2Wheel * (w / kPwt2WheelRes) + res[w % kPwt2WheelRes]; }
};

// Position of a prime in the table universe: odd-prime index below B, then
// the wheel-210 residues from B upward.
class Pwt2Universe {
public:
    Pwt2Universe(const PrimeBitmap& bm, u64 B) : bm_(bm), B_(B), cntB_(bm.count()) {
        u64 x = B;
        while (wheel_.ridx[x % kPwt2Wheel] < 0) ++x;
        wi0_ = wheel_.index(x);
    }
    u64 pos(u64 p) const {
        if (p < B_) {
            if (!bm_.is_prime(p)) throw std::runtime_error("pwt2: rank-table value below the sieve limit is not prime");
            return bm_.rank(p);
        }
        return cntB_ + wheel_.index(p) - wi0_;
    }
    u64 value(u64 pos) const { return pos < cntB_ ? bm_.select(pos) : wheel_.value(pos - cntB_ + wi0_); }
    bool high(u64 pos) const { return pos >= cntB_; }
    u64 low_count() const { return cntB_; }

private:
    const PrimeBitmap& bm_;
    u64 B_, cntB_, wi0_ = 0;
    Pwt2Wheel wheel_;
};

// ---------------------------------------------------------------------------
// Contexts (context set 1).
// ---------------------------------------------------------------------------
struct Pwt2Contexts {
    int K = 0;
    int table_base = 0, count_base = 0, flag_base = 0, gap_base = 0, total = 0;
    explicit Pwt2Contexts(int k = 1) : K(k) {
        table_base = 0;
        count_base = table_base + 41 * 2;
        flag_base = count_base + (K + 1) * kPwtPrimeBitsCtx;
        gap_base = flag_base + K;
        total = gap_base + K * kPwtBuckets * kPwtPrimeBitsCtx;
    }
    static int cap(u64 v, int hi) { return std::min(pwt_bits_of(v), hi); }
    int table_ctx(u64 class_size, bool prev_high) const { return table_base + cap(class_size, 40) * 2 + (prev_high ? 1 : 0); }
    int count_ctx(int depth, u64 vidx) const { return count_base + depth * kPwtPrimeBitsCtx + cap(vidx, kPwtPrimeBitsCtx - 1); }
    int flag_ctx(int depth) const { return flag_base + (depth - 1); }
    int gap_ctx(int child_depth, u64 c, u64 parent_rank) const {
        return gap_base + ((child_depth - 1) * kPwtBuckets + pwt_bucket(c)) * kPwtPrimeBitsCtx
               + cap(parent_rank, kPwtPrimeBitsCtx - 1);
    }
    u32 alphabet(int ctx) const {
        if (ctx < count_base) return kPwt2LenAlphabet;
        if (ctx < flag_base) return kPwtCountAlphabet;
        if (ctx < gap_base) return 2;
        return kPwt2LenAlphabet;
    }
};

// ---------------------------------------------------------------------------
// Context set 2 (variant gapv 15 / cntv 3 / flagv 1 / tabv 3, no mantissa
// symbols). Every
// symbol class has SPARSE keys; the file stores one model per key that
// occurs, and context ids are 32 bits (the full tree needs more than 2^16).
// All features are known to the decoder before the symbol.
// ---------------------------------------------------------------------------
constexpr int kPwt2S2Classes = 5;   // set 2 uses classes 0..4
constexpr int kPwt2S3Classes = 6;   // set 3 = set 2 + class 5, the tree-gap mantissa
enum { kC2Tab = 0, kC2Cnt = 1, kC2Tail = 2, kC2Flag = 3, kC2Gap = 4, kC2Mant = 5 };
constexpr u64 kPwt2MantR = 1024;   // mantissa context: the exact previous rank a when a < R
constexpr int kPwt2MantLx = 12;    // exact mantissa for L <= Lx, else its top kPwt2MantT bits
constexpr int kPwt2MantT = 4;
inline int pwt2_mant_bits(int L) { return L <= kPwt2MantLx ? L - 1 : kPwt2MantT; }   // coded bits of an L-bit gap
constexpr u32 kPwt2S2Alphabet[kPwt2S2Classes] = {64, 64, 4096, 2, 64};
constexpr u32 kPwt2S2CountTop = 63;     // count symbol min(c, 63); 63 -> tail symbol
constexpr u32 kPwt2S2TailTop = 4095;    // tail symbol min(c, 4095); 4095 -> Elias gamma of c - 4094

struct Pwt2Keys2 {
    u64 m = 0;
    u64 AB = 0;   // bits(m) + 1: the range of bits() of a rank or of an expected gap
    explicit Pwt2Keys2(u64 mm) : m(mm), AB(static_cast<u64>(pwt_bits_of(mm)) + 1) {}
    static u64 bits128(u128 v) {
        const u64 hi = static_cast<u64>(v >> 64);
        return hi ? 64 + static_cast<u64>(pwt_bits_of(hi)) : static_cast<u64>(pwt_bits_of(static_cast<u64>(v)));
    }
    static u64 capb(u64 v, u64 hi) { const u64 b = static_cast<u64>(pwt_bits_of(v)); return b < hi ? b : hi; }
    // bits(floor(X / R)) for R >= 1 without a division (decode speed):
    // with k = bits(X) - bits(R), floor(X / R) lies in [2^(k-1), 2^(k+1)), and it
    // reaches 2^k iff X >= R * 2^k iff (X >> k) >= R. Same values as the division.
    static u64 bits_quot(u64 X, u64 R) {
        const int bx = pwt_bits_of(X), br = pwt_bits_of(R);
        if (bx < br) return 0;
        const int k = bx - br;
        return static_cast<u64>((X >> k) >= R ? k + 1 : k);
    }
    // rank-table gap: class size, previous member at or above B, previous member's position (0 for the first)
    u64 table(u64 s, bool prev_high, u64 prev_pos) const {
        return (capb(s, 40) * 2 + (prev_high ? 1 : 0)) * 42 + capb(prev_pos, 41);
    }
    // count of node u at depth d, rank r, P = product of the primes on u's path (u's included; root: 1);
    // first = u is its parent's first child, after = siblings still to come after u (root: 0, 0)
    u64 count(int d, u64 r, u128 P, bool first, u64 after) const {
        const u64 Pb = std::min<u64>(bits128(P), 81);
        const u64 k = (static_cast<u64>(d) * AB + static_cast<u64>(pwt_bits_of(r))) * 82 + Pb;
        return (k * 2 + (first ? 1 : 0)) * 14 + capb(after, 13);
    }
    u64 tail(int d) const { return static_cast<u64>(d); }
    u64 flag(int d, u64 c) const { return static_cast<u64>(d - 1) * 4 + static_cast<u64>(pwt_bucket(c)); }
    // gap of the i-th child of v (child depth dc): c = child count of v, a = rank of the previous sibling
    // (v's own rank for the first child), rem = c - i children still to come including this one,
    // P = product of the primes on v's path
    // set 3: tree-gap mantissa of an L-bit gap (L >= 2): the exact previous rank a (the parent's
    // rank for a first child) when a < 1024, else (child depth, first, bits(a), bits((m - a) / R));
    // L is the low 6 bits of the key
    u64 mant(int dc, bool first, u64 a, u64 rem, int L) const {
        u64 b;
        if (a < kPwt2MantR) b = (first ? kPwt2MantR : 0) + a;
        else b = 2 * kPwt2MantR + ((static_cast<u64>(dc - 1) * 2 + (first ? 1 : 0)) * AB + static_cast<u64>(pwt_bits_of(a))) * AB
                 + bits_quot(m - a, rem);
        return b * 64 + static_cast<u64>(L);
    }
    u64 gap(int dc, u64 c, bool first, u64 a, u64 rem, u128 P) const {
        const u64 k = ((static_cast<u64>(dc - 1) * 4 + static_cast<u64>(pwt_bucket(c))) * 2 + (first ? 1 : 0)) * AB
                      + static_cast<u64>(pwt_bits_of(a));
        const u64 eb = bits_quot(m - a, rem);
        const u64 Pb = std::min<u64>(bits128(P), 81);
        return ((k * AB + eb) * 14 + capb(rem, 13)) * 41 + Pb / 2;
    }
};

// open-addressing u64 key -> u32 value map (keys are stored + 1, 0 = empty)
class Pwt2KeyMap {
public:
    Pwt2KeyMap() { reset_(1 << 10); }
    std::size_t size() const { return n_; }
    void put(u64 key, u32 v) {
        if (key == ~u64(0)) throw std::runtime_error("pwt2: context key out of range");
        if (2 * (n_ + 1) > keys_.size()) grow_();
        std::size_t i = h_(key) & mask_;
        while (keys_[i] && keys_[i] != key + 1) i = (i + 1) & mask_;
        if (!keys_[i]) ++n_;
        keys_[i] = key + 1;
        vals_[i] = v;
    }
    long long get(u64 key) const {
        std::size_t i = h_(key) & mask_;
        while (keys_[i]) {
            if (keys_[i] == key + 1) return static_cast<long long>(vals_[i]);
            i = (i + 1) & mask_;
        }
        return -1;
    }

private:
    static std::size_t h_(u64 k) {
        k += 1;
        k ^= k >> 33; k *= 0xff51afd7ed558ccdull; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ull; k ^= k >> 33;
        return static_cast<std::size_t>(k);
    }
    void reset_(std::size_t cap) { keys_.assign(cap, 0); vals_.assign(cap, 0); mask_ = cap - 1; n_ = 0; }
    void grow_() {
        std::vector<u64> ok; std::vector<u32> ov;
        ok.swap(keys_); ov.swap(vals_);
        reset_(ok.size() * 2);
        for (std::size_t i = 0; i < ok.size(); ++i) if (ok[i]) put(ok[i] - 1, ov[i]);
    }
    std::vector<u64> keys_;
    std::vector<u32> vals_;
    std::size_t mask_ = 0, n_ = 0;
};

// per-class histograms over sparse keys (encoder pass 1)
struct Pwt2Hist2 {
    u32 A = 0;
    Pwt2KeyMap index;
    std::vector<u64> keys, rows;   // row r: counts rows[r * A .. r * A + A)
    void add(u64 key, u64 sym) {
        if (sym >= A) throw std::runtime_error("pwt2: set-2 symbol above its alphabet");
        long long r = index.get(key);
        if (r < 0) {
            r = static_cast<long long>(keys.size());
            index.put(key, static_cast<u32>(r));
            keys.push_back(key);
            rows.resize(rows.size() + A, 0);
        }
        ++rows[static_cast<std::size_t>(r) * A + sym];
    }
};

// Flat model store for context set 2 decoding: every model's (cum, freq,
// sym) entries in one contiguous array, followed by a sentinel whose cum is
// 65,536; models with more than 16 symbols also get a 256-bucket slot
// table. Built from the deserialized models; decoding is identical to
// RansContext's, only the memory layout differs (58K-94K models stay in
// cache instead of four heap arrays each).
struct Pwt2FlatModels {
    struct E { u32 cum; u16 fm1; u16 sym; };   // fm1 = freq - 1 (a single-symbol model has freq 65,536)
    std::vector<E> e;
    std::vector<u32> off;
    std::vector<u32> lut_off;   // ~0u = no slot table (linear scan)
    std::vector<u16> lut;       // 256 per large model: index within the model
    void build(const std::vector<RansContext>& models) {
        e.clear(); off.clear(); lut_off.clear(); lut.clear();
        off.reserve(models.size()); lut_off.reserve(models.size());
        for (const RansContext& m : models) {
            off.push_back(static_cast<u32>(e.size()));
            for (std::size_t i = 0; i < m.syms.size(); ++i)
                e.push_back({m.cum[i], static_cast<u16>(m.freq[i] - 1), m.syms[i]});
            e.push_back({kRansProbScale, 0, 0});
            if (m.syms.size() > 16) {
                lut_off.push_back(static_cast<u32>(lut.size()));
                std::size_t k = 0;
                for (u32 b = 0; b < 256; ++b) {
                    const u32 slot = b << 8;
                    while (k + 1 < m.syms.size() && m.cum[k + 1] <= slot) ++k;
                    lut.push_back(static_cast<u16>(k));
                }
            } else {
                lut_off.push_back(~u32(0));
            }
            if (e.size() >= 0xFFFFFFF0ull) throw std::runtime_error("pwt2: flat model store too large");
        }
    }
    // (sym, freq, cum) of the entry holding slot in model id
    u16 find(u32 id, u32 slot, u32& f, u32& c) const {
        const E* m = e.data() + off[id];
        std::size_t k = 0;
        const u32 lo = lut_off[id];
        if (lo != ~u32(0)) k = lut[lo + (slot >> 8)];
        while (m[k + 1].cum <= slot) ++k;
        f = static_cast<u32>(m[k].fm1) + 1u;
        c = m[k].cum;
        return m[k].sym;
    }
};

// strict base-128 VByte (value below 2^64, canonical), for the set-2 tables
inline u64 pwt2_vbyte(const std::uint8_t* p, u64 len, u64& pos) {
    u64 v = 0;
    for (int shift = 0;; shift += 7) {
        if (pos >= len || shift > 63) throw std::runtime_error("pwt2: truncated VByte in the set-2 tables");
        const std::uint8_t b = p[pos++];
        if (shift == 63 && b > 1) throw std::runtime_error("pwt2: VByte above 2^64 in the set-2 tables");
        v |= static_cast<u64>(b & 0x7f) << shift;
        if (!(b & 0x80)) {
            if (shift > 0 && b == 0) throw std::runtime_error("pwt2: non-canonical VByte in the set-2 tables");
            return v;
        }
    }
}

struct Pwt2Header {
    bool full = false, naive = false;
    u64 n_records = 0;
    int max_depth = 0;
    int ctx_set = kPwt2CtxSet;
    u64 sieve_limit = 0, m = 0, n_classes = 0;
    u64 classes_off = 0, classes_len = 0, tables_off = 0, tables_len = 0;
    u64 sym_off = 0, sym_len = 0, raw_off = 0, raw_len = 0;
    u64 n_blocks = 0, total_symbols = 0, table_symbols = 0, footer_off = 0;
    u32 block_symbols = 0;
};

struct Pwt2EncodeStats {
    u64 nodes = 0, symbols = 0, table_symbols = 0, raw_bits = 0, table_raw_bits = 0, escapes = 0;
    u64 classes_bytes = 0, tables_bytes = 0, sym_bytes = 0, raw_bytes = 0, file_bytes = 0;
    u64 table_naive_bytes = 0;         // byte-aligned variant: bytes of the rank table inside the stream
    double ideal_symbol_bits = 0.0, entropy_symbol_bits = 0.0;
    double table_ideal_bits = 0.0;     // quantized ideal of the table's symbols alone
    u64 contexts = 0;                  // set 2: models stored
    double class_ideal_bits[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0};
    u64 class_table_bytes[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0}, class_contexts[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0};
};

// In-memory source for tests: paths as prime sets (any order) mapped to
// ranks and sorted in trie pre-order (ascending rank sequences).
class Pwt2VectorSource {
public:
    Pwt2VectorSource(const std::vector<std::vector<u64>>& paths, const Pwt2Order& ord) {
        std::unordered_map<u64, u64> rk;
        for (u64 r = 1; r <= ord.m(); ++r) rk[ord.prime[r]] = r;
        for (const auto& p : paths) {
            std::vector<u64> v;
            for (const u64 q : p) {
                const auto it = rk.find(q);
                if (it == rk.end()) throw std::runtime_error("pwt2 test source: prime not in the order");
                v.push_back(it->second);
            }
            std::sort(v.begin(), v.end());
            if (v.empty()) throw std::runtime_error("empty path");
            max_len_ = std::max(max_len_, static_cast<int>(v.size()));
            paths_.push_back(std::move(v));
        }
        std::sort(paths_.begin(), paths_.end());
        for (std::size_t i = 1; i < paths_.size(); ++i)
            if (paths_[i] == paths_[i - 1]) throw std::runtime_error("duplicate path");
    }
    std::size_t size() const { return paths_.size(); }
    int max_len() const { return max_len_; }
    template <typename F>
    void for_each(F&& f) const { for (const auto& p : paths_) f(p.data(), static_cast<int>(p.size())); }

private:
    std::vector<std::vector<u64>> paths_;
    int max_len_ = 0;
};

// Path counts over prime sets (for tests and small inputs).
inline Pwt2Order pwt2_order_from_paths(const std::vector<std::vector<u64>>& paths) {
    std::unordered_map<u64, u64> c;
    for (const auto& p : paths) for (const u64 q : p) ++c[q];
    std::vector<std::pair<u64, u64>> pc(c.begin(), c.end());
    return pwt2_order_from_counts(std::move(pc));
}

// ---------------------------------------------------------------------------
// Encoder. The source delivers RANK paths (ascending ranks) in trie
// pre-order; two passes (statistics, emission) as in PWT1.
// ---------------------------------------------------------------------------
class Pwt2Encoder {
public:
    Pwt2Encoder(bool full, bool naive, u64 sieve_limit, int ctx_set = kPwt2CtxSet)
        : full_(full), naive_(naive), B_(checked_limit_(sieve_limit)), bitmap_(B_), uni_(bitmap_, B_), ctx_set_(ctx_set) {
        if (ctx_set_ < 1 || ctx_set_ > 3) throw std::runtime_error("pwt2: context set must be 1, 2 or 3");
        if (naive_ && ctx_set_ != 1) throw std::runtime_error("pwt2: the byte-aligned variant has no contexts (set 1 only)");
    }

    template <typename Source>
    std::vector<std::uint8_t> encode(const Source& src, const Pwt2Order& ord, const PwtTargets& targets,
                                     Pwt2EncodeStats* stats = nullptr) {
        if (targets.record_count != src.size()) throw std::runtime_error("pwt2: target record count != source size");
        n_ = src.size();
        K_ = src.max_len();
        if (K_ < 1 || K_ > kPwt2MaxDepth) throw std::runtime_error("pwt2: bad max depth");
        M_ = ord.m();
        if (ord.vidx.size() != M_ + 1) throw std::runtime_error("pwt2: order not finished");
        ctx_ = Pwt2Contexts(K_);
        if (ctx_set_ >= 2) return encode2_(src, ord, targets, stats);

        // ---- pass 1: histograms + pre-order child counts --------------------
        hist_.assign(ctx_.total, {});
        for (int c = 0; c < ctx_.total; ++c) hist_[c].assign(ctx_.alphabet(c), 0);
        table_walk_(ord, [&](int ctx, int len, u64, int) { hist_[ctx][len]++; });
        counts_.clear();
        counts_.reserve(n_ * 2 + 16);
        overflow_.clear();
        {
            struct Open { u64 rank; u32 index; u64 c; bool terminal; std::vector<std::uint8_t> child_len; };
            std::vector<Open> st(K_ + 1);
            std::vector<u64> prev(K_ + 1, 0);
            int prev_len = 0;
            bool first = true;
            auto close = [&](int depth) {
                Open& o = st[depth];
                const u64 c = o.c;
                hist_[ctx_.count_ctx(depth, depth == 0 ? M_ + 1 : ord.vidx[o.rank])][std::min<u64>(c, kPwtCountAlphabet - 1)]++;
                if (c > 0 && depth >= 1) hist_[ctx_.flag_ctx(depth)][o.terminal ? 1 : 0]++;
                for (const std::uint8_t len : o.child_len) hist_[ctx_.gap_ctx(depth + 1, c, o.rank)][len]++;
                o.child_len.clear();
                if (c >= 65535) { counts_[o.index] = 65535; overflow_[o.index] = c; }
                else counts_[o.index] = static_cast<u16>(c);
            };
            auto open = [&](int depth, u64 rank, bool terminal) {
                Open& o = st[depth];
                o.rank = rank; o.terminal = terminal; o.c = 0;
                o.index = static_cast<u32>(counts_.size());
                if (counts_.size() >= 0xFFFFFFFFull) throw std::runtime_error("pwt2: too many nodes");
                counts_.push_back(0);
                o.child_len.clear();
                if (depth >= 1) {
                    Open& p = st[depth - 1];
                    const u64 g = rank - (p.c == 0 ? p.rank : prev[depth]);
                    p.child_len.push_back(static_cast<std::uint8_t>(pwt_bits_of(g)));
                    ++p.c;
                }
            };
            src.for_each([&](const u64* ranks, int len) {
                if (len > K_ || len < 1) throw std::runtime_error("pwt2: bad path length");
                for (int i = 0; i < len; ++i) {
                    if (ranks[i] < 1 || ranks[i] > M_) throw std::runtime_error("pwt2: rank outside 1..m");
                    if (i > 0 && ranks[i] <= ranks[i - 1]) throw std::runtime_error("pwt2: ranks not strictly ascending on a path");
                }
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == ranks[lcp]) ++lcp;
                if (lcp == len) throw std::runtime_error("pwt2: duplicate or unsorted path");   // pre-order puts a prefix before its extensions
                for (int d = prev_len; d > lcp; --d) close(d);
                for (int d = lcp + 1; d <= len; ++d) open(d, ranks[d - 1], d == len);
                for (int d = 0; d < len; ++d) prev[d + 1] = ranks[d];
                prev_len = len;
            });
            if (first) throw std::runtime_error("pwt2: empty source");
            for (int d = prev_len; d >= 0; --d) close(d);
        }
        nodes_ = counts_.size();

        // ---- models ----------------------------------------------------------
        models_.assign(ctx_.total, RansContext());
        std::vector<std::uint8_t> tables;
        double ideal = 0.0, entropy = 0.0, table_ideal = 0.0;
        for (int c = 0; c < ctx_.total; ++c) {
            models_[c].build(hist_[c].data(), hist_[c].size());
            models_[c].serialize(tables);
            double total = 0.0;
            for (const u64 h : hist_[c]) total += double(h);
            for (std::size_t i = 0; i < models_[c].syms.size(); ++i) {
                const double h = double(hist_[c][models_[c].syms[i]]);
                const double bits = h * (double(kRansProbBits) - std::log2(double(models_[c].freq[i])));
                ideal += bits;
                if (c < ctx_.count_base) table_ideal += bits;
                entropy += h * std::log2(total / h);
            }
        }
        hist_.clear();
        hist_.shrink_to_fit();
        if (naive_) tables.clear();

        // ---- class sizes ------------------------------------------------------
        std::vector<std::uint8_t> classes;
        for (const u64 s : ord.class_size) pwt_detail::put_vbyte(classes, s);

        // ---- pass 2: emission -------------------------------------------------
        std::vector<std::uint8_t> sym_bytes, raw_bytes;
        FastBitWriter raw(raw_bytes);
        std::vector<std::pair<u16, u16>> block;
        block.reserve(kPwtBlockSymbols);
        u64 total_symbols = 0, n_blocks = 0, raw_bits = 0, table_raw_bits = 0, escapes = 0;
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
        // rank table first
        table_walk_(ord, [&](int ctx, int len, u64 g, int) {
            if (naive_) return;
            emit(ctx, static_cast<u32>(len));
            raw.put(g, len - 1);
            raw_bits += len - 1;
            table_raw_bits += len - 1;
        });
        u64 table_naive_bytes = 0;
        if (naive_) {   // VByte prime gaps inside each class, ascending
            u64 r = 1;
            for (const u64 s : ord.class_size) {
                u64 prevp = 0;
                for (u64 i = 0; i < s; ++i) {
                    const u64 p = ord.prime[r + s - 1 - i];
                    const std::size_t before = sym_bytes.size();
                    pwt_detail::put_vbyte(sym_bytes, i == 0 ? p : p - prevp);
                    table_naive_bytes += sym_bytes.size() - before;
                    ++total_symbols;
                    prevp = p;
                }
                r += s;
            }
        }
        const u64 table_symbols = total_symbols;
        {
            struct Open { u64 rank; u64 c; u64 seen; };
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
            auto open = [&](int depth, u64 rank, bool terminal) {
                if (depth >= 1) {
                    Open& p = st[depth - 1];
                    const u64 g = rank - (p.seen == 0 ? p.rank : prev[depth]);
                    if (naive_) {
                        pwt_detail::put_vbyte(sym_bytes, g);
                        ++total_symbols;
                    } else {
                        const int len = pwt_bits_of(g);
                        emit(ctx_.gap_ctx(depth, p.c, p.rank), static_cast<u32>(len));
                        raw.put(g, len - 1);
                        raw_bits += len - 1;
                    }
                    ++p.seen;
                }
                Open& o = st[depth];
                o.rank = rank;
                o.c = node_count(next_index++);
                o.seen = 0;
                if (naive_) {
                    pwt_detail::put_vbyte(sym_bytes, 2 * o.c + (terminal ? 1 : 0));
                    ++total_symbols;
                    return;
                }
                const int cctx = ctx_.count_ctx(depth, depth == 0 ? M_ + 1 : ord.vidx[rank]);
                if (o.c >= kPwtCountAlphabet - 1) {
                    emit(cctx, kPwtCountAlphabet - 1);
                    const u64 v = o.c - (kPwtCountAlphabet - 1) + 1;   // >= 1, Elias gamma
                    const int vb = pwt_bits_of(v);
                    raw.put(0, vb - 1);
                    raw.put(v, vb);
                    raw_bits += 2 * vb - 1;
                    ++escapes;
                } else {
                    emit(cctx, static_cast<u32>(o.c));
                }
                if (o.c > 0 && depth >= 1) emit(ctx_.flag_ctx(depth), terminal ? 1u : 0u);
            };
            u128 value_sum = 0;
            u64 pr[kPwt2MaxDepth];
            src.for_each([&](const u64* ranks, int len) {
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == ranks[lcp]) ++lcp;
                for (int d = lcp + 1; d <= len; ++d) open(d, ranks[d - 1], d == len);
                for (int d = 0; d < len; ++d) { prev[d + 1] = ranks[d]; pr[d] = ord.prime[ranks[d]]; }
                prev_len = len;
                value_sum += pwt_detail::path_value(pr, len, full_);
            });
            if (next_index != nodes_) throw std::runtime_error("pwt2: pass 2 visited a different node count");
            if (value_sum != targets.total_check)
                throw std::runtime_error("pwt2: the paths do not reproduce the oracle's total_check "
                                         "(wrong, stale or truncated paths file); nothing written");
            for (int d = 0; d <= K_; ++d) if (st[d].seen > st[d].c) throw std::runtime_error("pwt2: child count mismatch");
        }
        flush();
        raw.align();
        counts_.clear();
        counts_.shrink_to_fit();
        overflow_.clear();

        // ---- assemble ----------------------------------------------------------
        Pwt2Header h;
        h.full = full_; h.naive = naive_; h.n_records = n_; h.max_depth = K_; h.sieve_limit = B_;
        h.m = M_; h.n_classes = ord.class_size.size();
        h.classes_off = kPwt2HeaderBytes; h.classes_len = classes.size();
        h.tables_off = h.classes_off + h.classes_len; h.tables_len = tables.size();
        h.sym_off = h.tables_off + h.tables_len; h.sym_len = sym_bytes.size();
        h.raw_off = h.sym_off + h.sym_len; h.raw_len = raw_bytes.size();
        h.n_blocks = n_blocks; h.total_symbols = total_symbols; h.table_symbols = table_symbols;
        h.block_symbols = naive_ ? 0 : kPwtBlockSymbols;
        h.footer_off = h.raw_off + h.raw_len;
        std::vector<std::uint8_t> out;
        out.reserve(h.footer_off + kPwt2FooterBytes);
        write_header_(out, h);
        out.insert(out.end(), classes.begin(), classes.end());
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
            stats->nodes = nodes_; stats->symbols = total_symbols; stats->table_symbols = table_symbols;
            stats->raw_bits = raw_bits; stats->table_raw_bits = table_raw_bits; stats->escapes = escapes;
            stats->classes_bytes = classes.size(); stats->tables_bytes = tables.size();
            stats->sym_bytes = sym_bytes.size(); stats->raw_bytes = raw_bytes.size();
            stats->file_bytes = out.size(); stats->table_naive_bytes = table_naive_bytes;
            stats->ideal_symbol_bits = ideal; stats->entropy_symbol_bits = entropy; stats->table_ideal_bits = table_ideal;
        }
        return out;
    }

private:
    // ---- context set 2 encoder (same two passes; sparse keyed models) ----
    template <typename F>
    void table_walk2_(const Pwt2Order& ord, const Pwt2Keys2& KY, F&& f) const {
        u64 r = 1;
        for (std::size_t j = 0; j < ord.class_size.size(); ++j) {
            const u64 s = ord.class_size[j];
            u64 prev_pos = 0;
            bool prev_high = false;
            for (u64 i = 0; i < s; ++i) {
                const u64 p = ord.prime[r + s - 1 - i];
                const u64 pos = uni_.pos(p);
                const u64 g = (i == 0) ? pos + 1 : pos - prev_pos;
                if (g == 0 || (i > 0 && pos <= prev_pos)) throw std::runtime_error("pwt2: class not ascending in the universe");
                const int len = pwt_bits_of(g);
                if (len > 56) throw std::runtime_error("pwt2: table gap too long");
                f(KY.table(s, prev_high, i == 0 ? 0 : prev_pos), len, g);
                prev_pos = pos;
                prev_high = uni_.high(pos);
            }
            r += s;
        }
    }

    template <typename Source>
    std::vector<std::uint8_t> encode2_(const Source& src, const Pwt2Order& ord, const PwtTargets& targets,
                                       Pwt2EncodeStats* stats) {
        const Pwt2Keys2 KY(M_);
        // ---- pass 1: histograms + pre-order child counts --------------------
        Pwt2Hist2 H[kPwt2S2Classes];
        for (int c = 0; c < kPwt2S2Classes; ++c) H[c].A = kPwt2S2Alphabet[c];
        const bool mant = (ctx_set_ == 3);
        std::unordered_map<u64, std::vector<u32>> MH;   // set 3: mantissa key -> counts over 2^vb values
        table_walk2_(ord, KY, [&](u64 key, int len, u64) { H[kC2Tab].add(key, static_cast<u64>(len)); });
        counts_.clear();
        counts_.reserve(n_ * 2 + 16);
        overflow_.clear();
        {
            struct Kid { u64 rank; u64 c; bool term; };
            struct Open { u64 rank; u128 P; u32 index; bool terminal; std::vector<Kid> kids; };
            std::vector<Open> st(K_ + 1);
            std::vector<u64> prev(K_ + 1, 0);
            int prev_len = 0;
            bool first = true;
            auto count_sym = [&](int d, u64 key, u64 c) {
                H[kC2Cnt].add(key, std::min<u64>(c, kPwt2S2CountTop));
                if (c >= kPwt2S2CountTop) H[kC2Tail].add(KY.tail(d), std::min<u64>(c, kPwt2S2TailTop));
            };
            auto close = [&](int depth) {
                Open& o = st[depth];
                const u64 c = o.kids.size();
                if (depth == 0) count_sym(0, KY.count(0, 0, 1, false, 0), c);
                else { Kid& me = st[depth - 1].kids.back(); me.c = c; me.term = o.terminal; }
                for (u64 i = 0; i < c; ++i) {
                    const Kid& kd = o.kids[i];
                    const u64 a = (i == 0) ? o.rank : o.kids[i - 1].rank;
                    H[kC2Gap].add(KY.gap(depth + 1, c, i == 0, a, c - i, o.P), static_cast<u64>(pwt_bits_of(kd.rank - a)));
                    if (mant) {
                        const u64 g = kd.rank - a;
                        const int L = pwt_bits_of(g);
                        if (L >= 2) {
                            const int vb = pwt2_mant_bits(L);
                            std::vector<u32>& h = MH[KY.mant(depth + 1, i == 0, a, c - i, L)];
                            if (h.empty()) h.assign(std::size_t(1) << vb, 0);
                            ++h[static_cast<std::size_t>((g - (u64(1) << (L - 1))) >> (L - 1 - vb))];
                        }
                    }
                    const u128 Pk = o.P * ord.prime[kd.rank];
                    count_sym(depth + 1, KY.count(depth + 1, kd.rank, Pk, i == 0, c - i - 1), kd.c);
                    if (kd.c > 0) H[kC2Flag].add(KY.flag(depth + 1, kd.c), kd.term ? 1 : 0);
                }
                o.kids.clear();
                if (c >= 65535) { counts_[o.index] = 65535; overflow_[o.index] = c; }
                else counts_[o.index] = static_cast<u16>(c);
            };
            auto open = [&](int depth, u64 rank, bool terminal) {
                Open& o = st[depth];
                o.rank = rank; o.terminal = terminal;
                o.P = (depth == 0) ? u128(1) : st[depth - 1].P * ord.prime[rank];
                o.index = static_cast<u32>(counts_.size());
                if (counts_.size() >= 0xFFFFFFFFull) throw std::runtime_error("pwt2: too many nodes");
                counts_.push_back(0);
                o.kids.clear();
                if (depth >= 1) st[depth - 1].kids.push_back({rank, 0, false});
            };
            src.for_each([&](const u64* ranks, int len) {
                if (len > K_ || len < 1) throw std::runtime_error("pwt2: bad path length");
                for (int i = 0; i < len; ++i) {
                    if (ranks[i] < 1 || ranks[i] > M_) throw std::runtime_error("pwt2: rank outside 1..m");
                    if (i > 0 && ranks[i] <= ranks[i - 1]) throw std::runtime_error("pwt2: ranks not strictly ascending on a path");
                }
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == ranks[lcp]) ++lcp;
                if (lcp == len) throw std::runtime_error("pwt2: duplicate or unsorted path");
                for (int d = prev_len; d > lcp; --d) close(d);
                for (int d = lcp + 1; d <= len; ++d) open(d, ranks[d - 1], d == len);
                for (int d = 0; d < len; ++d) prev[d + 1] = ranks[d];
                prev_len = len;
            });
            if (first) throw std::runtime_error("pwt2: empty source");
            for (int d = prev_len; d >= 0; --d) close(d);
        }
        nodes_ = counts_.size();

        // ---- models: per class, rows in ascending key order; id = running index
        std::vector<std::uint8_t> tables;
        models_.clear();
        std::vector<Pwt2KeyMap> ids(kPwt2S3Classes);
        double ideal = 0.0, entropy = 0.0, table_ideal = 0.0;
        double cls_ideal[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0};
        u64 cls_tab[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0}, cls_ctx[kPwt2S3Classes] = {0, 0, 0, 0, 0, 0};
        for (int cls = 0; cls < kPwt2S2Classes; ++cls) {
            Pwt2Hist2& hc = H[cls];
            std::vector<std::pair<u64, u64>> order;   // (key, row)
            order.reserve(hc.keys.size());
            for (std::size_t r = 0; r < hc.keys.size(); ++r) order.emplace_back(hc.keys[r], r);
            std::sort(order.begin(), order.end());
            const std::size_t t0 = tables.size();
            pwt_detail::put_vbyte(tables, order.size());
            u64 prevk = 0;
            for (const auto& kr : order) {
                pwt_detail::put_vbyte(tables, kr.first - prevk);
                prevk = kr.first;
                const u64* row = hc.rows.data() + kr.second * hc.A;
                RansContext rc;
                rc.build(row, hc.A);
                rc.serialize(tables);
                double total = 0.0;
                for (u32 s = 0; s < hc.A; ++s) total += double(row[s]);
                for (std::size_t i = 0; i < rc.syms.size(); ++i) {
                    const double h = double(row[rc.syms[i]]);
                    const double bits = h * (double(kRansProbBits) - std::log2(double(rc.freq[i])));
                    ideal += bits;
                    cls_ideal[cls] += bits;
                    if (cls == kC2Tab) table_ideal += bits;
                    entropy += h * std::log2(total / h);
                }
                if (models_.size() >= 0xFFFFFFFFull) throw std::runtime_error("pwt2: too many contexts");
                ids[cls].put(kr.first, static_cast<u32>(models_.size()));
                models_.push_back(std::move(rc));
                ++cls_ctx[cls];
            }
            cls_tab[cls] = tables.size() - t0;
            hc.rows.clear(); hc.rows.shrink_to_fit();
        }
        if (mant) {   // keep a mantissa model only where it beats the raw bits (table bytes charged)
            std::vector<u64> keys;
            keys.reserve(MH.size());
            for (const auto& kv : MH) keys.push_back(kv.first);
            std::sort(keys.begin(), keys.end());
            std::vector<std::uint8_t> sect, tmp;
            u64 kept = 0, prevk = 0;
            std::vector<u64> cnt;
            for (const u64 key : keys) {
                const std::vector<u32>& h = MH[key];
                cnt.assign(h.begin(), h.end());
                u64 N = 0;
                for (const u64 x : cnt) N += x;
                RansContext rc;
                rc.build(cnt.data(), cnt.size());
                double q = 0.0;
                for (std::size_t i = 0; i < rc.syms.size(); ++i)
                    q += double(cnt[rc.syms[i]]) * (double(kRansProbBits) - std::log2(double(rc.freq[i])));
                tmp.clear();
                rc.serialize(tmp);
                const int vb = pwt2_mant_bits(static_cast<int>(key % 64));
                if (double(vb) * double(N) - q - 8.0 * double(tmp.size() + 3) <= 0.0) continue;
                pwt_detail::put_vbyte(sect, key - prevk);
                prevk = key;
                sect.insert(sect.end(), tmp.begin(), tmp.end());
                ideal += q;
                cls_ideal[kC2Mant] += q;
                for (std::size_t i = 0; i < rc.syms.size(); ++i) {
                    const double hh = double(cnt[rc.syms[i]]);
                    entropy += hh * std::log2(double(N) / hh);
                }
                if (models_.size() >= 0xFFFFFFFFull) throw std::runtime_error("pwt2: too many contexts");
                ids[kC2Mant].put(key, static_cast<u32>(models_.size()));
                models_.push_back(std::move(rc));
                ++kept;
            }
            const std::size_t t0 = tables.size();
            pwt_detail::put_vbyte(tables, kept);
            tables.insert(tables.end(), sect.begin(), sect.end());
            cls_tab[kC2Mant] = tables.size() - t0;
            cls_ctx[kC2Mant] = kept;
            MH.clear();
        }

        // ---- class sizes ------------------------------------------------------
        std::vector<std::uint8_t> classes;
        for (const u64 s : ord.class_size) pwt_detail::put_vbyte(classes, s);

        // ---- pass 2: emission -------------------------------------------------
        std::vector<std::uint8_t> sym_bytes, raw_bytes;
        FastBitWriter raw(raw_bytes);
        std::vector<std::pair<u32, u16>> block;
        block.reserve(kPwtBlockSymbols);
        u64 total_symbols = 0, n_blocks = 0, raw_bits = 0, table_raw_bits = 0, escapes = 0;
        auto flush = [&]() {
            if (block.empty()) return;
            std::vector<std::uint8_t> enc;
            rans_encode_block_ids(block, models_, enc);
            pwt_detail::put_le(sym_bytes, block.size(), 4);
            pwt_detail::put_le(sym_bytes, enc.size(), 4);
            sym_bytes.insert(sym_bytes.end(), enc.begin(), enc.end());
            block.clear();
            ++n_blocks;
        };
        auto emit = [&](int cls, u64 key, u64 sym) {
            const long long id = ids[cls].get(key);
            if (id < 0) throw std::runtime_error("pwt2: set-2 context without a model (pass 1 / pass 2 disagree)");
            block.emplace_back(static_cast<u32>(id), static_cast<u16>(sym));
            ++total_symbols;
            if (block.size() == kPwtBlockSymbols) flush();
        };
        table_walk2_(ord, KY, [&](u64 key, int len, u64 g) {
            emit(kC2Tab, key, static_cast<u64>(len));
            raw.put(g, len - 1);
            raw_bits += len - 1;
            table_raw_bits += len - 1;
        });
        const u64 table_symbols = total_symbols;
        {
            struct Open { u64 rank; u128 P; u64 c; u64 seen; };
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
            auto count_emit = [&](int d, u64 key, u64 c) {
                emit(kC2Cnt, key, std::min<u64>(c, kPwt2S2CountTop));
                if (c < kPwt2S2CountTop) return;
                emit(kC2Tail, KY.tail(d), std::min<u64>(c, kPwt2S2TailTop));
                if (c < kPwt2S2TailTop) return;
                const u64 v = c - (kPwt2S2TailTop - 1);   // >= 1, Elias gamma
                const int vb = pwt_bits_of(v);
                raw.put(0, vb - 1);
                raw.put(v, vb);
                raw_bits += 2 * vb - 1;
                ++escapes;
            };
            auto open = [&](int depth, u64 rank, bool terminal) {
                Open& o = st[depth];
                o.rank = rank;
                o.c = node_count(next_index++);
                o.seen = 0;
                if (depth == 0) { o.P = 1; count_emit(0, KY.count(0, 0, 1, false, 0), o.c); return; }
                Open& p = st[depth - 1];
                const u64 i = p.seen;
                const u64 a = (i == 0) ? p.rank : prev[depth];
                const u64 g = rank - a;
                const int len = pwt_bits_of(g);
                emit(kC2Gap, KY.gap(depth, p.c, i == 0, a, p.c - i, p.P), static_cast<u64>(len));
                long long mid = -1;
                if (mant && len >= 2) mid = ids[kC2Mant].get(KY.mant(depth, i == 0, a, p.c - i, len));
                if (mid >= 0) {
                    const int vb = pwt2_mant_bits(len), rest = len - 1 - vb;
                    block.emplace_back(static_cast<u32>(mid), static_cast<u16>((g - (u64(1) << (len - 1))) >> rest));
                    ++total_symbols;
                    if (block.size() == kPwtBlockSymbols) flush();
                    raw.put(g, rest);
                    raw_bits += rest;
                } else {
                    raw.put(g, len - 1);
                    raw_bits += len - 1;
                }
                ++p.seen;
                o.P = p.P * ord.prime[rank];
                count_emit(depth, KY.count(depth, rank, o.P, i == 0, p.c - i - 1), o.c);
                if (o.c > 0) emit(kC2Flag, KY.flag(depth, o.c), terminal ? 1 : 0);
            };
            u128 value_sum = 0;
            u64 pr[kPwt2MaxDepth];
            src.for_each([&](const u64* ranks, int len) {
                if (first) { open(0, 0, false); first = false; }
                int lcp = 0;
                while (lcp < prev_len && lcp < len && prev[lcp + 1] == ranks[lcp]) ++lcp;
                for (int d = lcp + 1; d <= len; ++d) open(d, ranks[d - 1], d == len);
                for (int d = 0; d < len; ++d) { prev[d + 1] = ranks[d]; pr[d] = ord.prime[ranks[d]]; }
                prev_len = len;
                value_sum += pwt_detail::path_value(pr, len, full_);
            });
            if (next_index != nodes_) throw std::runtime_error("pwt2: pass 2 visited a different node count");
            if (value_sum != targets.total_check)
                throw std::runtime_error("pwt2: the paths do not reproduce the oracle's total_check "
                                         "(wrong, stale or truncated paths file); nothing written");
            for (int d = 0; d <= K_; ++d) if (st[d].seen > st[d].c) throw std::runtime_error("pwt2: child count mismatch");
        }
        flush();
        raw.align();
        counts_.clear();
        counts_.shrink_to_fit();
        overflow_.clear();

        Pwt2Header h;
        h.full = full_; h.naive = false; h.n_records = n_; h.max_depth = K_; h.sieve_limit = B_; h.ctx_set = ctx_set_;
        h.m = M_; h.n_classes = ord.class_size.size();
        h.classes_off = kPwt2HeaderBytes; h.classes_len = classes.size();
        h.tables_off = h.classes_off + h.classes_len; h.tables_len = tables.size();
        h.sym_off = h.tables_off + h.tables_len; h.sym_len = sym_bytes.size();
        h.raw_off = h.sym_off + h.sym_len; h.raw_len = raw_bytes.size();
        h.n_blocks = n_blocks; h.total_symbols = total_symbols; h.table_symbols = table_symbols;
        h.block_symbols = kPwtBlockSymbols;
        h.footer_off = h.raw_off + h.raw_len;
        std::vector<std::uint8_t> out;
        out.reserve(h.footer_off + kPwt2FooterBytes);
        write_header_(out, h);
        out.insert(out.end(), classes.begin(), classes.end());
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
            stats->nodes = nodes_; stats->symbols = total_symbols; stats->table_symbols = table_symbols;
            stats->raw_bits = raw_bits; stats->table_raw_bits = table_raw_bits; stats->escapes = escapes;
            stats->classes_bytes = classes.size(); stats->tables_bytes = tables.size();
            stats->sym_bytes = sym_bytes.size(); stats->raw_bytes = raw_bytes.size();
            stats->file_bytes = out.size();
            stats->ideal_symbol_bits = ideal; stats->entropy_symbol_bits = entropy; stats->table_ideal_bits = table_ideal;
            stats->contexts = models_.size();
            for (int c = 0; c < kPwt2S3Classes; ++c) {
                stats->class_ideal_bits[c] = cls_ideal[c]; stats->class_table_bytes[c] = cls_tab[c]; stats->class_contexts[c] = cls_ctx[c];
            }
        }
        return out;
    }

    static u64 checked_limit_(u64 B) {
        if (B < 211 || B > (u64(1) << 34)) throw std::runtime_error("pwt2: sieve limit must be in [211, 2^34]");
        return B;
    }
    // Visit the rank table's symbols in stream order: f(ctx, len, gap, class).
    template <typename F>
    void table_walk_(const Pwt2Order& ord, F&& f) const {
        u64 r = 1;
        for (std::size_t j = 0; j < ord.class_size.size(); ++j) {
            const u64 s = ord.class_size[j];
            u64 prev_pos = 0;
            bool prev_high = false;
            for (u64 i = 0; i < s; ++i) {
                const u64 p = ord.prime[r + s - 1 - i];   // ascending inside the class
                const u64 pos = uni_.pos(p);
                const u64 g = (i == 0) ? pos + 1 : pos - prev_pos;
                if (g == 0 || (i > 0 && pos <= prev_pos)) throw std::runtime_error("pwt2: class not ascending in the universe");
                const int len = pwt_bits_of(g);
                if (len > 56) throw std::runtime_error("pwt2: table gap too long");
                f(ctx_.table_ctx(s, prev_high), len, g, static_cast<int>(j));
                prev_pos = pos;
                prev_high = uni_.high(pos);
            }
            r += s;
        }
    }
    static void write_header_(std::vector<std::uint8_t>& out, const Pwt2Header& h) {
        out.insert(out.end(), kPwt2Magic, kPwt2Magic + 4);
        pwt_detail::put_le(out, kPwt2Version, 2);
        pwt_detail::put_le(out, (h.full ? kPwt2FlagFull : 0) | kPwt2FlagFreq | (h.naive ? kPwt2FlagNaive : 0), 2);
        pwt_detail::put_le(out, h.n_records, 8);                         // 8
        out.push_back(static_cast<std::uint8_t>(h.max_depth));           // 16
        out.push_back(static_cast<std::uint8_t>(h.ctx_set));             // 17
        out.push_back(0); out.push_back(0);                              // 18
        pwt_detail::put_le(out, kPwt2Wheel, 4);                          // 20
        pwt_detail::put_le(out, h.sieve_limit, 8);                       // 24
        pwt_detail::put_le(out, h.m, 8);                                 // 32
        pwt_detail::put_le(out, h.n_classes, 8);                         // 40
        pwt_detail::put_le(out, h.classes_off, 8); pwt_detail::put_le(out, h.classes_len, 8);   // 48, 56
        pwt_detail::put_le(out, h.tables_off, 8);  pwt_detail::put_le(out, h.tables_len, 8);    // 64, 72
        pwt_detail::put_le(out, h.sym_off, 8);     pwt_detail::put_le(out, h.sym_len, 8);       // 80, 88
        pwt_detail::put_le(out, h.raw_off, 8);     pwt_detail::put_le(out, h.raw_len, 8);       // 96, 104
        pwt_detail::put_le(out, h.n_blocks, 8);                          // 112
        pwt_detail::put_le(out, h.block_symbols, 4);                     // 120
        pwt_detail::put_le(out, 0, 4);                                   // 124
        pwt_detail::put_le(out, h.total_symbols, 8);                     // 128
        pwt_detail::put_le(out, h.table_symbols, 8);                     // 136
        pwt_detail::put_le(out, h.footer_off, 8);                        // 144
        while (out.size() < kPwt2HeaderBytes) out.push_back(0);          // 152..159 reserved
    }

    bool full_, naive_;
    u64 B_;
    PrimeBitmap bitmap_;
    Pwt2Universe uni_;
    std::size_t n_ = 0, nodes_ = 0;
    int K_ = 0;
    u64 M_ = 0;
    Pwt2Contexts ctx_;
    int ctx_set_ = kPwt2CtxSet;
    std::vector<std::vector<u64>> hist_;
    std::vector<RansContext> models_;
    std::vector<u16> counts_;
    std::unordered_map<u32, u64> overflow_;
};

// ---------------------------------------------------------------------------
// Reader / decoder.
// ---------------------------------------------------------------------------
class Pwt2Reader {
public:
    explicit Pwt2Reader(std::vector<std::uint8_t> bytes) : buf_(std::move(bytes)) {
        if (buf_.size() < kPwt2HeaderBytes + kPwt2FooterBytes) throw std::runtime_error("pwt2: file too small");
        const std::uint8_t* p = buf_.data();
        if (std::memcmp(p, kPwt2Magic, 4) != 0) throw std::runtime_error("pwt2: bad magic");
        if (pwt_detail::get_le(p + 4, 2) != kPwt2Version) throw std::runtime_error("pwt2: unsupported version");
        const u64 flags = pwt_detail::get_le(p + 6, 2);
        if (!(flags & kPwt2FlagFreq) || (flags & ~u64(kPwt2FlagFull | kPwt2FlagNaive | kPwt2FlagFreq)))
            throw std::runtime_error("pwt2: bad flags");
        h_.full = (flags & kPwt2FlagFull) != 0;
        h_.naive = (flags & kPwt2FlagNaive) != 0;
        h_.n_records = pwt_detail::get_le(p + 8, 8);
        h_.max_depth = p[16];
        h_.ctx_set = p[17];
        if (p[18] != 0 || p[19] != 0) throw std::runtime_error("pwt2: nonzero reserved header bytes");
        if (pwt_detail::get_le(p + 20, 4) != kPwt2Wheel) throw std::runtime_error("pwt2: unsupported wheel");
        h_.sieve_limit = pwt_detail::get_le(p + 24, 8);
        h_.m = pwt_detail::get_le(p + 32, 8);
        h_.n_classes = pwt_detail::get_le(p + 40, 8);
        h_.classes_off = pwt_detail::get_le(p + 48, 8); h_.classes_len = pwt_detail::get_le(p + 56, 8);
        h_.tables_off = pwt_detail::get_le(p + 64, 8);  h_.tables_len = pwt_detail::get_le(p + 72, 8);
        h_.sym_off = pwt_detail::get_le(p + 80, 8);     h_.sym_len = pwt_detail::get_le(p + 88, 8);
        h_.raw_off = pwt_detail::get_le(p + 96, 8);     h_.raw_len = pwt_detail::get_le(p + 104, 8);
        h_.n_blocks = pwt_detail::get_le(p + 112, 8);
        h_.block_symbols = static_cast<u32>(pwt_detail::get_le(p + 120, 4));
        if (pwt_detail::get_le(p + 124, 4) != 0) throw std::runtime_error("pwt2: nonzero reserved header bytes");
        h_.total_symbols = pwt_detail::get_le(p + 128, 8);
        h_.table_symbols = pwt_detail::get_le(p + 136, 8);
        h_.footer_off = pwt_detail::get_le(p + 144, 8);
        for (std::size_t i = 152; i < kPwt2HeaderBytes; ++i) if (p[i]) throw std::runtime_error("pwt2: nonzero reserved header bytes");
        if (h_.max_depth < 1 || h_.max_depth > kPwt2MaxDepth) throw std::runtime_error("pwt2: bad max depth");
        if (h_.ctx_set < 1 || h_.ctx_set > 3) throw std::runtime_error("pwt2: unsupported context set");
        if (h_.naive && h_.ctx_set != 1) throw std::runtime_error("pwt2: byte-aligned file with a context set other than 1");
        if (h_.sieve_limit < 211 || h_.sieve_limit > (u64(1) << 34)) throw std::runtime_error("pwt2: bad sieve limit");
        if (h_.m == 0 || h_.m >= (u64(1) << 32) || h_.n_classes == 0 || h_.n_classes > h_.m)
            throw std::runtime_error("pwt2: bad rank table size");
        if (h_.naive) {
            if (h_.block_symbols != 0 || h_.n_blocks != 0 || h_.tables_len != 0 || h_.raw_len != 0)
                throw std::runtime_error("pwt2: byte-aligned header carries entropy-coded sections");
        } else if (h_.block_symbols != kPwtBlockSymbols) throw std::runtime_error("pwt2: block_symbols must be 2^20");
        // geometry: every bound is a subtraction from the real footer offset
        const u64 footer_off = buf_.size() - kPwt2FooterBytes;
        if (h_.footer_off != footer_off || h_.classes_off != kPwt2HeaderBytes ||
            h_.classes_len > footer_off - h_.classes_off ||
            h_.tables_off != h_.classes_off + h_.classes_len ||
            h_.tables_len > footer_off - h_.tables_off ||
            h_.sym_off != h_.tables_off + h_.tables_len ||
            h_.sym_len > footer_off - h_.sym_off ||
            h_.raw_off != h_.sym_off + h_.sym_len ||
            h_.raw_len != footer_off - h_.raw_off)
            throw std::runtime_error("pwt2: section geometry does not tile the file");
        if (h_.table_symbols > h_.total_symbols) throw std::runtime_error("pwt2: table symbols exceed the total");
        if (h_.table_symbols != h_.m) throw std::runtime_error("pwt2: the rank table must be exactly m symbols");
        if (h_.n_records > h_.total_symbols) throw std::runtime_error("pwt2: more records than symbols");
        if (h_.naive && h_.total_symbols > h_.sym_len) throw std::runtime_error("pwt2: more VBytes than stream bytes");
        Sha256 sha;
        sha.update(p, h_.footer_off);
        const auto got = sha.finish();
        const std::uint8_t* f = p + h_.footer_off;
        if (std::memcmp(got.data(), f, 32) != 0) throw std::runtime_error("pwt2: payload SHA-256 mismatch");
        std::memcpy(targets_.sha_nset.data(), f + 32, 32);
        targets_.total_check = pwt_detail::get_le128(f + 64);
        targets_.record_count = pwt_detail::get_le(f + 80, 8);
        if (targets_.record_count != h_.n_records) throw std::runtime_error("pwt2: footer record count != header");
        // class sizes
        {
            std::size_t pos = 0;
            u64 tot = 0;
            order_.class_size.clear();
            for (u64 j = 0; j < h_.n_classes; ++j) {
                u64 v = 0;
                for (int shift = 0;; shift += 7) {
                    if (pos >= h_.classes_len || shift > 63) throw std::runtime_error("pwt2: truncated class sizes");
                    const std::uint8_t b = p[h_.classes_off + pos++];
                    if (shift == 63 && b > 1) throw std::runtime_error("pwt2: class size VByte above 2^64");
                    v |= static_cast<u64>(b & 0x7f) << shift;
                    if (!(b & 0x80)) { if (shift > 0 && b == 0) throw std::runtime_error("pwt2: non-canonical class size VByte"); break; }
                }
                if (v == 0 || v > h_.m - tot) throw std::runtime_error("pwt2: bad class size");
                tot += v;
                order_.class_size.push_back(v);
            }
            if (tot != h_.m || pos != h_.classes_len) throw std::runtime_error("pwt2: class sizes do not cover m exactly");
        }
        if (h_.ctx_set >= 2) ky2_ = std::make_unique<Pwt2Keys2>(h_.m);
        if (h_.ctx_set >= 2) {   // sets 2 and 3: per class, VByte(rows), then (VByte(key delta), model) in ascending key order
            const std::uint8_t* t = p + h_.tables_off;
            const u64 tl = h_.tables_len;
            u64 tp = 0;
            models_.clear();
            const int ncls = (h_.ctx_set == 3) ? kPwt2S3Classes : kPwt2S2Classes;
            for (int cls = 0; cls < ncls; ++cls) {
                const u64 nrows = pwt2_vbyte(t, tl, tp);
                if (nrows > (tl - tp) / 3) throw std::runtime_error("pwt2: more set-2 rows than the tables can hold");
                u64 prevk = 0;
                for (u64 r = 0; r < nrows; ++r) {
                    const u64 dk = pwt2_vbyte(t, tl, tp);
                    if (r > 0 && dk == 0) throw std::runtime_error("pwt2: set-2 keys not strictly ascending");
                    if (dk > ~u64(0) - 1 - prevk) throw std::runtime_error("pwt2: set-2 key overflow");
                    const u64 key = prevk + dk;
                    u64 alpha = 0;
                    if (cls < kPwt2S2Classes) alpha = kPwt2S2Alphabet[cls];
                    else {   // mantissa: the gap length L is the key's low 6 bits
                        const int L = static_cast<int>(key % 64);
                        if (L < 2) throw std::runtime_error("pwt2: mantissa model for a gap length below 2");
                        alpha = u64(1) << pwt2_mant_bits(L);
                    }
                    RansContext rc;
                    tp += rc.deserialize(t + tp, tl - tp, alpha, true);   // strict: canonical VBytes (spec section 9)
                    if (rc.empty()) throw std::runtime_error("pwt2: empty set-2 model");
                    keys2_[cls].put(key, static_cast<u32>(models_.size()));
                    models_.push_back(std::move(rc));
                    prevk = key;
                }
            }
            if (tp != tl) throw std::runtime_error("pwt2: trailing bytes in the tables section");
            flat2_.build(models_);
        } else {
        ctx_ = Pwt2Contexts(h_.max_depth);
        models_.assign(ctx_.total, RansContext());
        std::size_t pos = 0;
        if (!h_.naive)
            for (int c = 0; c < ctx_.total; ++c)
                pos += models_[c].deserialize(p + h_.tables_off + pos, h_.tables_len - pos, ctx_.alphabet(c));
        if (pos != h_.tables_len) throw std::runtime_error("pwt2: trailing bytes in the tables section");
        }
        blocks_.clear();
        if (!h_.naive) {
            std::size_t q = 0;
            u64 syms = 0;
            while (q < h_.sym_len) {
                if (q + 8 > h_.sym_len) throw std::runtime_error("pwt2: truncated block header");
                const u64 ns = pwt_detail::get_le(p + h_.sym_off + q, 4);
                const u64 nb = pwt_detail::get_le(p + h_.sym_off + q + 4, 4);
                if (ns == 0 || ns > h_.block_symbols || nb > h_.sym_len - q - 8) throw std::runtime_error("pwt2: bad block");
                blocks_.push_back({h_.sym_off + q + 8, nb, ns});
                syms += ns;
                q += 8 + nb;
            }
            if (blocks_.size() != h_.n_blocks || syms != h_.total_symbols)
                throw std::runtime_error("pwt2: block directory disagrees with the header");
            for (std::size_t i = 0; i + 1 < blocks_.size(); ++i)
                if (blocks_[i].nsyms != h_.block_symbols) throw std::runtime_error("pwt2: short block before the last");
        }
        bitmap_ = std::make_unique<PrimeBitmap>(h_.sieve_limit);
    }

    const Pwt2Header& header() const { return h_; }
    const PwtTargets& targets() const { return targets_; }
    std::size_t file_bytes() const { return buf_.size(); }
    const Pwt2Order& order() const { return order_; }   // filled by for_each_n

    // Decode the rank table, then walk the whole tree; f(u128 n) once per
    // stored number in tree order. table_done() is called between the two.
    template <typename F, typename T>
    void for_each_n(F&& f, T&& table_done) {
        Walk w(*this);
        w.run(std::forward<F>(f), std::forward<T>(table_done));
    }
    template <typename F>
    void for_each_n(F&& f) { for_each_n(std::forward<F>(f), [] {}); }

    // As for_each_n, also handing over the stored path (for
    // pwt2_decode --emit-paths): f(u128 n, const u64* primes, int len) with
    // primes in root-to-node order (rank ascending, NOT value order).
    template <typename F, typename T>
    void for_each_path(F&& f, T&& table_done) {
        static_assert(std::is_invocable_v<F&, u128, const u64*, int>, "for_each_path: f(u128, const u64*, int)");
        Walk w(*this);
        w.run(std::forward<F>(f), std::forward<T>(table_done));
    }

private:
    struct Block { u64 off, len, nsyms; };

    class Walk {
    public:
        explicit Walk(Pwt2Reader& r) : r_(r), raw_(r.buf_.data() + r.h_.raw_off, r.h_.raw_len) {}
        template <typename F, typename T>
        void run(F&& f, T&& table_done) {
            const Pwt2Header& h = r_.h_;
            if (!h.naive) next_block_();
            if (h.ctx_set >= 2) read_table2_(); else read_table_();
            if (consumed_ != h.table_symbols) throw std::runtime_error("pwt2: rank table symbol count != header");
            table_done();
            path_.assign(h.max_depth + 1, 0);
            if (h.ctx_set >= 2) node2_(0, 0, 1, false, 0, f); else node_(0, 0, f);
            if (max_seen_ != h.max_depth) throw std::runtime_error("pwt2: max_depth is not the depth of the deepest path");
            if (emitted_ != h.n_records) throw std::runtime_error("pwt2: decoded record count != header");
            if (consumed_ != h.total_symbols) throw std::runtime_error("pwt2: symbol count != header");
            if (h.naive) {
                if (vb_pos_ != h.sym_len) throw std::runtime_error("pwt2: trailing bytes in the byte-aligned stream");
                return;
            }
            if (!dec_.finished()) throw std::runtime_error("pwt2: last block did not end cleanly");
            if (block_ != r_.blocks_.size()) throw std::runtime_error("pwt2: unread blocks");
            const u64 used = raw_.bit_pos();
            if ((used + 7) / 8 != h.raw_len) throw std::runtime_error("pwt2: raw bit stream length mismatch");
            if (used < raw_.total_bits() && raw_.get(static_cast<int>(raw_.total_bits() - used)) != 0)
                throw std::runtime_error("pwt2: nonzero raw padding");
        }

    private:
        // One stored number: f(n) (for_each_n), or f(n, primes, depth) when f
        // takes the path too (for_each_path). Chosen at compile time.
        template <typename F>
        void emit_(F& f, int depth) {
            const u128 n = pwt_detail::path_value(path_.data() + 1, depth, r_.h_.full);
            if constexpr (std::is_invocable_v<F&, u128, const u64*, int>) f(n, static_cast<const u64*>(path_.data() + 1), depth);
            else f(n);
        }
        u64 vb_() {
            const std::uint8_t* s = r_.buf_.data() + r_.h_.sym_off;
            u64 v = 0;
            for (int shift = 0;; shift += 7) {
                if (vb_pos_ >= r_.h_.sym_len || shift > 63) throw std::runtime_error("pwt2: byte-aligned stream overrun");
                const std::uint8_t b = s[vb_pos_++];
                if (shift == 63 && b > 1) throw std::runtime_error("pwt2: VByte above 2^64");
                v |= static_cast<u64>(b & 0x7f) << shift;
                if (!(b & 0x80)) { if (shift > 0 && b == 0) throw std::runtime_error("pwt2: non-canonical VByte"); break; }
            }
            ++consumed_;
            return v;
        }
        u32 sym_(int ctx) {
            if (in_block_ == cur_.nsyms) {
                if (!dec_.finished()) throw std::runtime_error("pwt2: block did not end cleanly");
                next_block_();
            }
            ++in_block_; ++consumed_;
            return dec_.decode(r_.models_[ctx]);
        }
        void next_block_() {
            if (block_ >= r_.blocks_.size()) throw std::runtime_error("pwt2: symbol stream exhausted");
            cur_ = r_.blocks_[block_++];
            dec_ = RansBlockDecoder(r_.buf_.data() + cur_.off, cur_.len);
            in_block_ = 0;
        }
        u32 sym2_(u64 id) {   // set 2: same block handling as sym_, flat model store
            if (in_block_ == cur_.nsyms) {
                if (!dec_.finished()) throw std::runtime_error("pwt2: block did not end cleanly");
                next_block_();
            }
            ++in_block_; ++consumed_;
            const Pwt2FlatModels& F = r_.flat2_;
            const u32 m = static_cast<u32>(id);
            return dec_.decode_by([&F, m](u32 slot, u32& f, u32& c) { return F.find(m, slot, f, c); });
        }
        u64 gap_len2_(u64 id) {
            const u32 len = sym2_(id);
            if (len == 0 || len >= kPwt2LenAlphabet) throw std::runtime_error("pwt2: bad gap length");
            return len == 1 ? 1 : ((u64(1) << (len - 1)) | raw_.get(static_cast<int>(len - 1)));
        }
        u64 gap_len_(int ctx) {
            const u32 len = sym_(ctx);
            if (len == 0 || len >= kPwt2LenAlphabet) throw std::runtime_error("pwt2: bad gap length");
            return len == 1 ? 1 : ((u64(1) << (len - 1)) | raw_.get(static_cast<int>(len - 1)));
        }
        // ---- context set 2 ----
        u64 id2_(int cls, u64 key) const {
            const long long id = r_.keys2_[cls].get(key);
            if (id < 0) throw std::runtime_error("pwt2: set-2 context without a model");
            return static_cast<u64>(id);
        }
        void read_table2_() {
            const Pwt2Header& h = r_.h_;
            Pwt2Order& ord = r_.order_;
            const Pwt2Keys2& KY = *r_.ky2_;
            ord.prime.clear();
            ord.prime.reserve(std::min<u64>(h.m, u64(1) << 24) + 1);
            ord.prime.push_back(0);
            std::vector<u64> cls;
            const Pwt2Universe uni(*r_.bitmap_, h.sieve_limit);
            const u64 max_pos = uni.low_count() + (u64(1) << 40);
            for (const u64 s : ord.class_size) {
                u64 pos = 0, prev_pos = 0;
                bool prev_high = false;
                cls.clear();
                for (u64 i = 0; i < s; ++i) {
                    const u64 g = gap_len2_(id2_(kC2Tab, KY.table(s, prev_high, i == 0 ? 0 : prev_pos)));
                    if (i == 0) pos = g - 1;
                    else { if (g > max_pos - pos) throw std::runtime_error("pwt2: table position overflow"); pos += g; }
                    if (pos >= max_pos) throw std::runtime_error("pwt2: table position out of range");
                    const u64 p = uni.value(pos);
                    if (p >= (u64(1) << 40)) throw std::runtime_error("pwt2: table prime too large");
                    prev_high = uni.high(pos);
                    cls.push_back(p);
                    prev_pos = pos;
                }
                for (std::size_t i = cls.size(); i-- > 0;) ord.prime.push_back(cls[i]);
            }
            ord.finish();
        }
        template <typename F>
        void node2_(int depth, u64 rank, u128 P, bool first, u64 after, F& f) {
            const Pwt2Header& h = r_.h_;
            const Pwt2Keys2& KY = *r_.ky2_;
            u64 c = sym2_(id2_(kC2Cnt, KY.count(depth, rank, P, first, after)));
            if (c == kPwt2S2CountTop) {
                c = sym2_(id2_(kC2Tail, KY.tail(depth)));
                if (c < kPwt2S2CountTop) throw std::runtime_error("pwt2: tail symbol below the count top");
                if (c == kPwt2S2TailTop) {
                    int zeros = 0;
                    while (raw_.get(1) == 0) if (++zeros > 40) throw std::runtime_error("pwt2: bad escape");
                    const u64 v = (u64(1) << zeros) | raw_.get(zeros);
                    c = v + kPwt2S2TailTop - 1;
                }
            }
            bool terminal = (c == 0);
            if (c > 0 && depth >= 1) terminal = sym2_(id2_(kC2Flag, KY.flag(depth, c))) == 1;
            if (depth == 0 && c == 0) throw std::runtime_error("pwt2: empty tree");
            if (depth > max_seen_) max_seen_ = depth;
            if (terminal) {
                if (depth == 0) throw std::runtime_error("pwt2: terminal root");
                emit_(f, depth);
                ++emitted_;
            }
            if (c > 0 && depth + 1 > h.max_depth) throw std::runtime_error("pwt2: tree deeper than the header");
            if (c > h.m) throw std::runtime_error("pwt2: child count above m");
            u64 prev = rank;
            for (u64 i = 0; i < c; ++i) {
                u64 g;
                {
                    const u32 L = sym2_(id2_(kC2Gap, KY.gap(depth + 1, c, i == 0, prev, c - i, P)));
                    if (L == 0 || L >= kPwt2LenAlphabet) throw std::runtime_error("pwt2: bad gap length");
                    if (L == 1) g = 1;
                    else {
                        long long mid = -1;
                        if (h.ctx_set == 3) mid = r_.keys2_[kC2Mant].get(KY.mant(depth + 1, i == 0, prev, c - i, static_cast<int>(L)));
                        if (mid >= 0) {
                            const int vb = pwt2_mant_bits(static_cast<int>(L)), rest = static_cast<int>(L) - 1 - vb;
                            const u64 v = sym2_(static_cast<u64>(mid));
                            g = (u64(1) << (L - 1)) | (v << rest) | raw_.get(rest);
                        } else {
                            g = (u64(1) << (L - 1)) | raw_.get(static_cast<int>(L - 1));
                        }
                    }
                }
                if (g == 0 || g > h.m - prev) throw std::runtime_error("pwt2: child rank outside (parent, m]");
                const u64 q = prev + g;
                const u64 pr = r_.order_.prime[q];
                path_[depth + 1] = pr;
                u128 Pc;
                if (__builtin_mul_overflow(P, static_cast<u128>(pr), &Pc)) throw std::runtime_error("pwt2: path product overflows 128 bits");
                node2_(depth + 1, q, Pc, i == 0, c - i - 1, f);
                prev = q;
            }
        }
        void read_table_() {
            const Pwt2Header& h = r_.h_;
            Pwt2Order& ord = r_.order_;
            ord.prime.clear();
            ord.prime.reserve(std::min<u64>(h.m, u64(1) << 24) + 1);   // grows with the decoded symbols, never from the header alone
            ord.prime.push_back(0);
            std::vector<u64> cls;
            const Pwt2Universe uni(*r_.bitmap_, h.sieve_limit);
            const u64 max_pos = uni.low_count() + (u64(1) << 40);   // far above any 40-bit value
            for (const u64 s : ord.class_size) {
                u64 pos = 0, prevp = 0;
                bool prev_high = false;
                cls.clear();
                for (u64 i = 0; i < s; ++i) {
                    u64 p;
                    if (h.naive) {
                        const u64 g = vb_();
                        if (g == 0 || (i > 0 && g > ~u64(0) - prevp)) throw std::runtime_error("pwt2: bad table gap");
                        p = (i == 0) ? g : prevp + g;
                        if (p >= (u64(1) << 40)) throw std::runtime_error("pwt2: table prime too large");
                        if (p < h.sieve_limit && !r_.bitmap_->is_prime(p)) throw std::runtime_error("pwt2: table value below the sieve limit is not prime");
                    } else {
                        const u64 g = gap_len_(r_.ctx_.table_ctx(s, prev_high));
                        if (i == 0) pos = g - 1;
                        else { if (g > max_pos - pos) throw std::runtime_error("pwt2: table position overflow"); pos += g; }
                        if (pos >= max_pos) throw std::runtime_error("pwt2: table position out of range");
                        p = uni.value(pos);
                        if (p >= (u64(1) << 40)) throw std::runtime_error("pwt2: table prime too large");
                        prev_high = uni.high(pos);
                    }
                    cls.push_back(p);   // ascending inside the class
                    prevp = p;
                }
                for (std::size_t i = cls.size(); i-- > 0;) ord.prime.push_back(cls[i]);   // ranks run descending in value
            }
            ord.finish();   // distinct odd primes, classes strictly descending, vidx
        }
        template <typename F>
        void node_(int depth, u64 rank, F& f) {
            const Pwt2Header& h = r_.h_;
            u64 c;
            bool terminal;
            if (h.naive) {
                const u64 v = vb_();
                c = v >> 1;
                terminal = (v & 1) != 0;
                if (c == 0 && !terminal) throw std::runtime_error("pwt2: leaf not terminal");
            } else {
                c = sym_(r_.ctx_.count_ctx(depth, depth == 0 ? h.m + 1 : r_.order_.vidx[rank]));
                if (c == kPwtCountAlphabet - 1) {
                    int zeros = 0;
                    while (raw_.get(1) == 0) if (++zeros > 40) throw std::runtime_error("pwt2: bad escape");
                    const u64 v = (u64(1) << zeros) | raw_.get(zeros);
                    c = v + kPwtCountAlphabet - 2;
                }
                terminal = (c == 0);
                if (c > 0 && depth >= 1) terminal = sym_(r_.ctx_.flag_ctx(depth)) == 1;
            }
            if (depth == 0 && c == 0) throw std::runtime_error("pwt2: empty tree");
            if (depth > max_seen_) max_seen_ = depth;
            if (terminal) {
                if (depth == 0) throw std::runtime_error("pwt2: terminal root");
                emit_(f, depth);
                ++emitted_;
            }
            if (c > 0 && depth + 1 > h.max_depth) throw std::runtime_error("pwt2: tree deeper than the header");
            if (c > h.m) throw std::runtime_error("pwt2: child count above m");
            u64 prev = rank;
            for (u64 i = 0; i < c; ++i) {
                const u64 g = h.naive ? vb_() : gap_len_(r_.ctx_.gap_ctx(depth + 1, c, rank));
                if (g == 0 || g > h.m - prev) throw std::runtime_error("pwt2: child rank outside (parent, m]");
                const u64 q = prev + g;
                path_[depth + 1] = r_.order_.prime[q];
                node_(depth + 1, q, f);
                prev = q;
            }
        }

        Pwt2Reader& r_;
        FastBitReader raw_;
        RansBlockDecoder dec_;
        Block cur_{};
        std::size_t block_ = 0;
        u64 in_block_ = 0, consumed_ = 0, emitted_ = 0;
        std::size_t vb_pos_ = 0;
        int max_seen_ = 0;
        std::vector<u64> path_;
    };

    std::vector<std::uint8_t> buf_;
    Pwt2Header h_;
    PwtTargets targets_;
    Pwt2Contexts ctx_;
    std::vector<RansContext> models_;
    std::vector<Block> blocks_;
    std::unique_ptr<PrimeBitmap> bitmap_;
    Pwt2Order order_;
    Pwt2KeyMap keys2_[kPwt2S3Classes];
    std::unique_ptr<Pwt2Keys2> ky2_;
    Pwt2FlatModels flat2_;
};

} // namespace cn
