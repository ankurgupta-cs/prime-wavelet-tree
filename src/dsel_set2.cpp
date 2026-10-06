// dsel_set2 -- the d* search: a local search for the global choice of the
// stored divisor per Carmichael number in the frequency-ordered tree,
// priced under the SHIPPING coder: the
// frequency-ordered PWT2 container with CONTEXT SET 2 (docs/PWT2_FORMAT.md
// sections 1-9; every key below is pwt2.hpp Pwt2Keys2, checked symbol by
// symbol at every fresh-rank build). MEASUREMENT ONLY: reads the ORC1 oracle
// (or fixture text); writes its report and, on request, path files and mask
// checkpoints (scratch). Source, universe, CSR tree with overflow nodes,
// external-sort build, DFS and batch machinery, the set-2 model, the rank
// table and the moves are all in this file.
//
// OBJECTIVE = the set-2 code length of the container for the current S_n:
//   * every coded symbol (rank-table gap length, count, count tail, flag,
//     child-gap length) in its set-2 model (class, key): the static empirical
//     entropy f(T) - sum_b f(h_b), f(x) = x log2 x (the encoder's quantized
//     tables track it to ~0.001 bits/el);
//   * raw bits: L - 1 per tree gap and per table gap, Elias gamma escapes;
//   * the tables section, 8 bits per byte: per class VByte(#models); per model
//     VByte(key delta) + VByte(#used symbols) + per used symbol VByte(symbol
//     delta) + VByte(freq - 1), freq = max(1, floor(65536 h / T)) (the
//     real quantizer's +-1 largest-remainder corrections are ignored; the
//     build prints the REAL serialized size next to it);
//   * the class-size VBytes of the rank table.
// Within a sweep the ranks and m (hence AB and every gap key's (m - a) / R)
// are FROZEN at the build (newly used primes keep their build rank > m: their
// bits(rank) is clamped to AB - 1 and m - a to >= 0, a state function); the
// rank table is always priced exactly from the live path counts. Between
// sweeps the ranks are recomputed and everything is rebuilt and repriced.
//
// A MOVE of n = remove its path, insert another valid S (prod(S) lambda(S) > n).
// Its cost is the EXACT change of the objective, symbol by symbol:
//   * removed / new nodes: gap, count (+ tail + escape), flag symbols;
//   * the parent v whose child count changes c -> c +- 1: its count value,
//     its flag (bucket(c) is in the flag key; no flag at c = 0);
//   * v's other children whose key changes (LEMMA below): the gap key has
//     bucket(c), first, bits(a), bits((m - a) / R), cap(R, 13); the count key
//     has first and cap(A, 13) (R = c - i, A = c - i - 1, i = index, a =
//     previous sibling's rank or rank(v));
//   * the successor sibling's gap (split / merged; a and first change);
//   * the terminal flag of an internal end node;
//   * the rank table: a prime whose path count moves a -> a +- 1 leaves class a
//     and joins class a +- 1 (gap split / merge, the successor's key, the class
//     size VByte; when bits(class size) changes, EVERY member is re-keyed via
//     the class's (H, bits(prev pos), L) sub-histogram);
//   * model tables: bins / models appearing and vanishing, freq-byte classes,
//     key deltas against the virtual key set.
// LEMMA (which siblings change). Insert a child at position i0 of v (c -> c+1)
// or delete the child at i0 (c -> c-1), bucket(c) unchanged. A sibling at
// j > i0 (resp. after the deleted one) keeps R = c - j and A because c and j
// shift together; its a and first are unchanged except for the immediate
// successor of the changed position. A sibling at j < i0 keeps first, a and
// its count value and changes only R -> R +- 1, A -> A +- 1. Hence its keys
// change iff one of three bit lengths crosses: cap(R, 13) iff min(R, R') + 1
// is a power of two <= 8192; cap(A, 13) iff min(R, R') is a power of two <=
// 8192; bits(floor(X / R)) with X = m - a iff floor(X / 2^k) = min(R, R') for
// the unique k = bits(X) - bits(min(R, R')) >= 0 (floor(X/(R+1)) < 2^k <=
// floor(X/R) <=> R 2^k <= X < (R+1) 2^k). If bucket(c) changes (c <= 16) all
// children are re-keyed. Nothing else in the tree changes (P(v), depth, ranks
// and m are fixed within a sweep). Small nodes test every j < i0 with the O(1)
// tests; nodes with more than --huge children keep a cache (live children +
// the positions whose quotient sits on a boundary at R or R - 1), so only
// O(#crossings + 42 cap positions) siblings are ever touched.
// Each node sees at most one detach (the old path) and one attach (the new
// path) per move, so a node's virtual child list = live list - D + A.
//
// SEARCH per n: remove its path (virtually), DFS over subsets
// of its primes in rank order: through existing live nodes (only the table
// changes), a divergence = attach, then a chain of new nodes; branch and bound
// on the chain with per-prime table lower bounds and a per-node slack for the
// freq-byte thresholds (--slack; only pruning, never pricing); validity
// reachability pruning. Ties keep S; among equal gains the smaller mask wins.
// T > 1: batches evaluated in parallel on the batch-start state, applied
// sequentially, each re-costed on the live state and committed only if it
// still gains (else re-searched on the live state).
//
// VALIDATION: every build recomputes the whole model from scratch; with fresh
// ranks every symbol's key is compared with pwt2.hpp's Pwt2Keys2; the build
// prints the real serialized table size (rans.hpp) next to the model's;
// --check-rebuild compares the live objective and a fingerprint of every
// live histogram with a same-rank rebuild after every sweep; fixture mode
// with --dump writes every considered candidate's cost, for an independent
// reference implementation that prices every candidate from scratch.
//
// usage: dsel_set2 (--orc1 F [--max-bits B] | --text FIX...) [--sweeps K]
//          [--threads T] [--batch R] [--scratch DIR] [--chunk R]
//          [--max-overflow M] [--check-rebuild] [--masks-in F] [--masks-out F]
//          [--emit F] [--emit-dmin F] [--dump F] [--no-prune] [--min-gain PCT]
//          [--step-cap S] [--dpaths F [--dmin-check K]] [--sample-stride S]
//          [--huge C] [--cb C] [--slack BITS] [--sieve-limit B]

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#endif

#include "lambda_bucket.hpp"
#include "orc1.hpp"
#include "pwt2.hpp"
#include "u128.hpp"

using namespace cn;

namespace {

using i64 = std::int64_t;
using i32 = std::int32_t;
using u8 = std::uint8_t;

[[noreturn]] void die(const char* fmt, ...) {
    std::fflush(stdout);
    va_list ap;
    va_start(ap, fmt);
    std::fprintf(stderr, "FATAL: ");
    std::vfprintf(stderr, fmt, ap);
    std::fprintf(stderr, "\n");
    va_end(ap);
    std::exit(1);
}
#define CHECK(c, ...) do { if (!(c)) die(__VA_ARGS__); } while (0)

double g_t0 = 0;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
double el() { return now() - g_t0; }

struct Mem { double ws_gb, peak_gb, commit_gb; };
#ifdef _WIN32
Mem mem() {
    PROCESS_MEMORY_COUNTERS pmc;
    std::memset(&pmc, 0, sizeof pmc);
    K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc);
    return {pmc.WorkingSetSize / 1e9, pmc.PeakWorkingSetSize / 1e9, pmc.PagefileUsage / 1e9};
}
inline unsigned long proc_id() { return GetCurrentProcessId(); }
#else
// POSIX: only the peak resident set is available (ru_maxrss: KiB on Linux,
// bytes on macOS); the report prints it in all three slots.
Mem mem() {
    struct rusage ru;
    std::memset(&ru, 0, sizeof ru);
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    const double peak = double(ru.ru_maxrss) / 1e9;
#else
    const double peak = double(ru.ru_maxrss) * 1024.0 / 1e9;
#endif
    return {peak, peak, peak};
}
inline unsigned long proc_id() { return static_cast<unsigned long>(getpid()); }
#endif

inline int bits_of(u64 v) { return v ? 64 - __builtin_clzll(v) : 0; }
inline u64 bits128(u128 v) { const u64 hi = static_cast<u64>(v >> 64); return hi ? 64 + bits_of(hi) : bits_of(static_cast<u64>(v)); }
inline u64 capb(u64 v, u64 hi) { const u64 b = static_cast<u64>(bits_of(v)); return b < hi ? b : hi; }
inline int bucket_of(u64 c) { return c == 1 ? 0 : c <= 3 ? 1 : c <= 15 ? 2 : 3; }
inline int vbl(u64 v) { int n = 1; while (v >= 0x80) { v >>= 7; ++n; } return n; }
inline u64 th1(u64 T) { return (129 * T + 65535) >> 16; }     // freq >= 129  <=> h >= th1
inline u64 th2(u64 T) { return (16385 * T + 65535) >> 16; }   // freq >= 16385 <=> h >= th2
inline int blen(u64 h, u64 T) { return h ? 1 + (h >= th1(T) ? 1 : 0) + (h >= th2(T) ? 1 : 0) : 0; }
inline bool pow2(u64 x) { return x && !(x & (x - 1)); }
// bits(floor(X / R)) != bits(floor(X / (R + 1)))
inline bool quot_cross(u64 X, u64 R) {
    if (R == 0) return false;
    const int bx = bits_of(X), br = bits_of(R);
    if (bx < br) return false;
    return (X >> (bx - br)) == R;
}
inline u64 fmix64(u64 k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33; return k;
}
constexpr double kLn2 = 0.69314718055994530942;

std::vector<double> g_df;
inline double df(u64 x) {
    if (x < g_df.size()) return g_df[x];
    const double X = double(x);
    return std::log2(X + 1.0) + X * std::log1p(1.0 / X) / kLn2;
}
void init_df() {
    g_df.resize(1u << 20);
    g_df[0] = 0.0;
    for (u64 x = 1; x < g_df.size(); ++x) { const double X = double(x); g_df[x] = std::log2(X + 1.0) + X * std::log1p(1.0 / X) / kLn2; }
}
inline double fx(u64 x) { return x > 1 ? double(x) * std::log2(double(x)) : 0.0; }

inline u64 gcd64b(u64 a, u64 b) {
    if (a == 0) return b;
    if (b == 0) return a;
    const int s = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do { b >>= __builtin_ctzll(b); if (a > b) std::swap(a, b); b -= a; } while (b);
    return a << s;
}
inline u128 lcm_step(u128 lam, u64 a) {
    const u64 r = (lam >> 64) ? static_cast<u64>(lam % a) : static_cast<u64>(lam) % a;
    return lam * (a / gcd64b(r, a));
}
inline bool valid_pl(u128 prod, u128 lam, u128 n) {
    u128 x;
    if (__builtin_mul_overflow(prod, lam, &x)) return true;
    return x > n;
}

// ------------------------------------------------------------------ options
struct Opt {
    std::string orc;
    std::vector<std::string> text;
    std::string scratch = ".";
    int max_bits = 0;
    int sweeps = 8;
    int threads = 1;
    u64 batch = 65536;
    u64 chunk = u64(1) << 25;
    u64 max_overflow = 150000000;
    bool check_rebuild = false;
    std::string masks_in, masks_out, emit, emit_dmin, dump;
    double min_gain = 0.01;
    u64 step_cap = 4000000;
    u64 stride = 1;
    std::string dpaths;
    u64 dmin_check = 997;
    u32 huge = 128;          // nodes with more live children keep a sibling cache
    u32 cb = 128;            // rank-table classes 1..cb as bitsets, larger as sorted vectors
    double slack = 16.0;     // per future chain node (pruning only)
    bool no_prune = false;
    double tw = 1.0;         // --table-weight: weight of the tables-section bytes in the objective (diagnostic; 1 = the file)
    u64 sieve_limit = 100000000ull;
};
Opt g_opt;

// ------------------------------------------------------------ record source
constexpr int kMaxK = 14;
struct Rec { u128 n; int k; u64 p[kMaxK]; };

struct Source {
    bool is_text = false;
    std::vector<Rec> text;
    std::string path;
    std::vector<std::unique_ptr<Orc1Reader>> rd;
    u64 bs = 4096, N = 0, nblocks = 0;
    u128 X = ~u128(0);

    void open(const Opt& o) {
        if (!o.text.empty()) {
            is_text = true;
            for (const std::string& f : o.text) {
                std::ifstream in(f);
                CHECK(in, "cannot open %s", f.c_str());
                std::string line;
                while (std::getline(in, line)) {
                    if (line.empty()) continue;
                    std::istringstream ss(line);
                    std::string tok;
                    ss >> tok;
                    Rec r{};
                    r.n = parse_u128(tok);
                    r.k = 0;
                    u64 v;
                    u128 prod = 1;
                    while (ss >> v) { CHECK(r.k < kMaxK, "too many factors"); r.p[r.k++] = v; prod *= v; }
                    CHECK(prod == r.n && r.k >= 3, "bad fixture line %s", line.c_str());
                    for (int j = 1; j < r.k; ++j) CHECK(r.p[j] > r.p[j - 1], "factors not ascending");
                    text.push_back(r);
                }
            }
            std::sort(text.begin(), text.end(), [](const Rec& a, const Rec& b) { return a.n < b.n; });
            for (std::size_t i = 1; i < text.size(); ++i) CHECK(text[i].n != text[i - 1].n, "duplicate fixture n");
            N = text.size();
            return;
        }
        path = o.orc;
        for (int t = 0; t < std::max(1, o.threads); ++t) rd.emplace_back(new Orc1Reader(path));
        const auto& ix = rd[0]->index();
        nblocks = ix.size();
        bs = rd[0]->header().block_size;
        if (o.max_bits > 0) X = u128(1) << o.max_bits;
        u64 b = 0;
        while (b < nblocks && ix[b].first_n < X) ++b;
        if (b == 0) { N = 0; return; }
        if (o.max_bits > 0) {
            u64 cnt = 0;
            rd[0]->for_block(b - 1, [&](const Orc1Record& r) { if (r.n < X) ++cnt; });
            N = ix[b - 1].rank_base + cnt;
        } else {
            N = rd[0]->header().n_records;
        }
        for (u64 i = 0; i < nblocks; ++i) CHECK(ix[i].rank_base == i * bs, "ORC1 block %llu not full-size aligned", (unsigned long long)i);
    }

    template <class F> void each(u64 r0, u64 r1, int T, F&& f) {
        if (r1 <= r0) return;
        if (is_text) {
            auto part = [&](int tid, u64 a, u64 b) { for (u64 i = a; i < b; ++i) f(tid, i, text[i]); };
            if (T <= 1) { part(0, r0, r1); return; }
            std::vector<std::thread> th;
            const u64 per = (r1 - r0 + T - 1) / T;
            for (int t = 0; t < T; ++t) {
                const u64 a = r0 + t * per, b = std::min(r1, a + per);
                if (a >= b) break;
                th.emplace_back(part, t, a, b);
            }
            for (auto& x : th) x.join();
            return;
        }
        const u64 b0 = r0 / bs, b1 = (r1 - 1) / bs + 1;
        auto part = [&](int tid, u64 ba, u64 bb) {
            Rec rec;
            for (u64 b = ba; b < bb; ++b) {
                u64 pos = b * bs;
                rd[tid]->for_block(b, [&](const Orc1Record& o) {
                    if (pos >= r0 && pos < r1) {
                        rec.n = o.n; rec.k = o.k;
                        for (int j = 0; j < o.k; ++j) rec.p[j] = o.factors[j];
                        f(tid, pos, static_cast<const Rec&>(rec));
                    }
                    ++pos;
                });
            }
        };
        if (T <= 1) { part(0, b0, b1); return; }
        std::vector<std::thread> th;
        const u64 nb = b1 - b0, per = (nb + T - 1) / T;
        for (int t = 0; t < T; ++t) {
            const u64 a = b0 + t * per, b = std::min(b1, a + per);
            if (a >= b) break;
            th.emplace_back(part, t, a, b);
        }
        for (auto& x : th) x.join();
    }
};

// ---------------------------------------------------------- prime universe
constexpr u64 kSmall = u64(1) << 24;
constexpr u32 kNone = 0xFFFFFFFFu;

struct Universe {
    std::vector<u64> P;          // sorted distinct primes of the (sub-)table
    std::vector<u64> pos;        // PWT2 table-universe position (pwt2.hpp Pwt2Universe)
    std::vector<u8> hi;          // pos >= cnt(B): the prime is at or above B
    std::vector<u32> small;
    std::vector<u64> hk; std::vector<u32> hv; u64 hm = 0;
    u32 idx(u64 p) const {
        if (p < kSmall) { const u32 v = small[p >> 1]; return v ? v - 1 : kNone; }
        u64 i = fmix64(p) & hm;
        while (hk[i] != p) { if (hk[i] == 0) return kNone; i = (i + 1) & hm; }
        return hv[i];
    }
    void finish(u64 B) {
        small.assign(kSmall >> 1, 0);
        u64 nbig = 0;
        for (u64 p : P) if (p >= kSmall) ++nbig;
        u64 cap = 16;
        while (cap < 2 * nbig + 16) cap <<= 1;
        hk.assign(cap, 0); hv.assign(cap, 0); hm = cap - 1;
        for (u32 i = 0; i < P.size(); ++i) {
            const u64 p = P[i];
            if (p < kSmall) small[p >> 1] = i + 1;
            else { u64 j = fmix64(p) & hm; while (hk[j] != 0) j = (j + 1) & hm; hk[j] = p; hv[j] = i; }
        }
        PrimeBitmap bm(B);
        Pwt2Universe uni(bm, B);
        pos.resize(P.size()); hi.resize(P.size());
        for (std::size_t i = 0; i < P.size(); ++i) {
            CHECK(P[i] >= 3 && (P[i] & 1), "prime %llu not odd", (unsigned long long)P[i]);
            pos[i] = uni.pos(P[i]);
            hi[i] = uni.high(pos[i]) ? 1 : 0;
            CHECK(i == 0 || pos[i] > pos[i - 1], "universe positions not ascending");
        }
    }
};

// -------------------------------------------------------------- the ranking
constexpr u32 kRankBits = 23;
constexpr u32 kRankMask = (1u << kRankBits) - 1;

struct Ranking {
    std::vector<u32> rank, r2i;   // idx -> rank (1-based); rank -> idx
    std::vector<u64> pv;          // rank -> prime value (pv[0] = 1)
    u32 m = 0;
};
Ranking make_ranking(const Universe& U, const std::vector<u32>& pc) {
    Ranking R;
    const u32 n = static_cast<u32>(U.P.size());
    CHECK(n < kRankMask, "too many primes for 23-bit ranks");
    std::vector<u32> ord(n);
    for (u32 i = 0; i < n; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](u32 a, u32 b) { return pc[a] != pc[b] ? pc[a] > pc[b] : U.P[a] > U.P[b]; });
    R.rank.assign(n, 0); R.r2i.assign(n + 1, kNone); R.pv.assign(n + 1, 1);
    for (u32 r = 0; r < n; ++r) { R.rank[ord[r]] = r + 1; R.r2i[r + 1] = ord[r]; R.pv[r + 1] = U.P[ord[r]]; }
    u32 used = 0;
    for (u32 i = 0; i < n; ++i) if (pc[i] > 0) ++used;
    R.m = used;
    return R;
}

// ------------------------------------------------------------ set-2 keys
// pwt2.hpp Pwt2Keys2 with m frozen at the build and two clamps for the ranks
// of primes that became used inside a sweep (rank > m): bits(rank) <= AB - 1,
// m - a >= 0. With fresh ranks (all used ranks <= m) the clamps never act and
// the keys equal Pwt2Keys2 (checked at every such build).
struct Keys {
    u64 m = 0, AB = 1;
    void set(u64 mf) { m = mf; AB = static_cast<u64>(bits_of(mf)) + 1; }
    u64 bA(u64 x) const { const u64 b = static_cast<u64>(bits_of(x)); return b < AB ? b : AB - 1; }
    u64 X(u64 a) const { return m >= a ? m - a : 0; }
    u64 count(int d, u64 r, u128 P, bool first, u64 after) const {
        const u64 Pb = std::min<u64>(bits128(P), 81);
        const u64 k = (static_cast<u64>(d) * AB + bA(r)) * 82 + Pb;
        return (k * 2 + (first ? 1 : 0)) * 14 + capb(after, 13);
    }
    u64 flag(int d, u64 c) const { return static_cast<u64>(d - 1) * 4 + static_cast<u64>(bucket_of(c)); }
    u64 gap(int dc, u64 c, bool first, u64 a, u64 rem, u128 P) const {
        const u64 k = ((static_cast<u64>(dc - 1) * 4 + static_cast<u64>(bucket_of(c))) * 2 + (first ? 1 : 0)) * AB + bA(a);
        const u64 eb = bA(X(a) / rem);
        const u64 Pb = std::min<u64>(bits128(P), 81);
        return ((k * AB + eb) * 14 + capb(rem, 13)) * 41 + Pb / 2;
    }
    u64 tablek(u64 s, bool H, u64 pb) const { return (capb(s, 40) * 2 + (H ? 1 : 0)) * 42 + pb; }
};

// ------------------------------------------------------------ model store
constexpr int kNC = 5;
enum { C_TAB = 0, C_CNT = 1, C_TAIL = 2, C_FLAG = 3, C_GAP = 4 };
constexpr u32 kAl[kNC] = {64, 64, 4096, 2, 64};
constexpr u32 kAW[kNC] = {1, 1, 64, 1, 1};
const char* kCName[kNC] = {"table", "count", "count tail", "flag", "gap length"};
inline u64 mkck(int cls, u64 key) { return (static_cast<u64>(cls) << 48) | key; }
inline int ck_cls(u64 ck) { return static_cast<int>(ck >> 48); }
inline u64 ck_key(u64 ck) { return ck & ((u64(1) << 48) - 1); }

// symbol-delta bytes of used symbol s (mask words mk, s's bit clear or set, it
// is ignored): VByte(s - pred) + VByte(succ - s) - VByte(succ - pred), pred = 0
// if none (the serialization's first delta is the symbol itself)
inline int sd_of(int cls, const u64* mk, u32 s) {
    if (cls != C_TAIL) return 1;   // alphabets <= 64: every delta < 128
    // pred
    u64 p = 0; bool hp = false;
    {
        i64 w = s >> 6; u64 x = mk[w] & ((s & 63) ? ((u64(1) << (s & 63)) - 1) : 0);
        for (;;) {
            if (x) { p = static_cast<u64>(w) * 64 + 63 - __builtin_clzll(x); hp = true; break; }
            if (--w < 0) break;
            x = mk[w];
        }
    }
    u64 q = 0; bool hq = false;
    {
        u32 w = s >> 6; u64 x = mk[w] & (((s & 63) == 63) ? 0 : ~((u64(2) << (s & 63)) - 1));
        for (;;) {
            if (x) { q = static_cast<u64>(w) * 64 + __builtin_ctzll(x); hq = true; break; }
            if (++w >= kAW[C_TAIL]) break;
            x = mk[w];
        }
    }
    (void)hp;
    const u64 pp = hp ? p : 0;
    return vbl(s - pp) + (hq ? vbl(q - s) - vbl(q - pp) : 0);
}
inline int fb_all(const u32* h, const u64* mk, u32 W, u64 T) {
    int s = 0;
    for (u32 w = 0; w < W; ++w) {
        u64 x = mk[w];
        while (x) { const int b = __builtin_ctzll(x); x &= x - 1; s += blen(h[w * 64 + b], T); }
    }
    return s;
}

struct Mdl { u64 ck; u64 T; u32 nsym; i32 fb; i32 sd; u32 off; u32 moff; u8 cls; };
struct Store {
    std::vector<Mdl> md;
    std::vector<u32> H;
    std::vector<u64> MK;
    std::vector<u64> hk; std::vector<u32> hv; u64 hm = 0, hn = 0;
    std::vector<u64> keys[kNC];   // sorted keys of the live models (T > 0)
    i64 nmod[kNC] = {0, 0, 0, 0, 0};
    double ent[kNC] = {0, 0, 0, 0, 0};
    i64 by[kNC] = {0, 0, 0, 0, 0};     // tables-section bytes of the class (model approximation)
    u64 syms[kNC] = {0, 0, 0, 0, 0};   // symbols (build only; not maintained in sweeps)
    void clear() {
        std::vector<Mdl>().swap(md); std::vector<u32>().swap(H); std::vector<u64>().swap(MK);
        hk.assign(1 << 12, ~u64(0)); hv.assign(1 << 12, 0); hm = (1 << 12) - 1; hn = 0;
        for (int c = 0; c < kNC; ++c) { std::vector<u64>().swap(keys[c]); nmod[c] = 0; ent[c] = 0; by[c] = 0; syms[c] = 0; }
    }
    i64 find(u64 ck) const {
        u64 i = fmix64(ck) & hm;
        for (;;) { const u64 k = hk[i]; if (k == ck) return hv[i]; if (k == ~u64(0)) return -1; i = (i + 1) & hm; }
    }
    void hput(u64 ck, u32 id) {
        u64 i = fmix64(ck) & hm;
        while (hk[i] != ~u64(0)) i = (i + 1) & hm;
        hk[i] = ck; hv[i] = id; ++hn;
    }
    u32 create(u64 ck) {
        if (2 * (hn + 1) > hk.size()) {
            std::vector<u64> ok; ok.swap(hk);
            hk.assign(ok.size() * 2, ~u64(0)); hv.assign(ok.size() * 2, 0); hm = hk.size() - 1; hn = 0;
            for (u32 i = 0; i < md.size(); ++i) hput(md[i].ck, i);
        }
        const int cls = ck_cls(ck);
        CHECK(cls >= 0 && cls < kNC, "bad class");
        const u32 id = static_cast<u32>(md.size());
        CHECK(H.size() + kAl[cls] < 0xFFFFFFF0ull, "model store full");
        md.push_back(Mdl{ck, 0, 0, 0, 0, static_cast<u32>(H.size()), static_cast<u32>(MK.size()), static_cast<u8>(cls)});
        H.resize(H.size() + kAl[cls], 0);
        MK.resize(MK.size() + kAW[cls], 0);
        hput(ck, id);
        return id;
    }
    // build-time accumulation
    void badd(int cls, u64 key, u32 s, u64 w) {
        CHECK(s < kAl[cls], "symbol %u above the alphabet of class %d", s, cls);
        const u64 ck = mkck(cls, key);
        i64 id = find(ck);
        if (id < 0) id = create(ck);
        H[md[id].off + s] += static_cast<u32>(w);
        md[id].T += w;
    }
    // after the build: masks, per-model byte state, entropies, key sets
    void finish() {
        for (int c = 0; c < kNC; ++c) { keys[c].clear(); nmod[c] = 0; ent[c] = 0; by[c] = 0; syms[c] = 0; }
        i64 mb[kNC] = {0, 0, 0, 0, 0};
        for (Mdl& m : md) {
            const u32 A = kAl[m.cls], W = kAW[m.cls];
            u64* mk = MK.data() + m.moff;
            const u32* h = H.data() + m.off;
            for (u32 w = 0; w < W; ++w) mk[w] = 0;
            m.nsym = 0; m.sd = 0;
            double e = fx(m.T);
            for (u32 s = 0; s < A; ++s) if (h[s]) { mk[s >> 6] |= u64(1) << (s & 63); ++m.nsym; e -= fx(h[s]); }
            if (!m.T) continue;
            // symbol deltas
            u64 prev = 0;
            for (u32 s = 0; s < A; ++s) if (h[s]) { m.sd += vbl(s - prev); prev = s; }
            m.fb = fb_all(h, mk, W, m.T);
            ent[m.cls] += e;
            keys[m.cls].push_back(ck_key(m.ck));
            ++nmod[m.cls];
            mb[m.cls] += vbl(m.nsym) + m.sd + m.fb;
            syms[m.cls] += m.T;
        }
        for (int c = 0; c < kNC; ++c) {
            std::sort(keys[c].begin(), keys[c].end());
            i64 kd = 0; u64 prev = 0;
            for (u64 k : keys[c]) { kd += vbl(k - prev); prev = k; }
            by[c] = vbl(static_cast<u64>(nmod[c])) + mb[c] + kd;
        }
    }
    // order-independent fingerprint of every live histogram
    u64 fingerprint() const {
        u64 f = 0;
        for (const Mdl& m : md) {
            if (!m.T) continue;
            const u32* h = H.data() + m.off;
            for (u32 s = 0; s < kAl[m.cls]; ++s) if (h[s]) f += fmix64(m.ck * 0x9E3779B97F4A7C15ull + s * 0x632BE59BD9B4E019ull + h[s]);
        }
        return f;
    }
};

// ------------------------------------------------------------- rank table
// Live classes by path count: members (idx order = prime order = pos order),
// class sizes, class-size VByte bytes, and per class the sub-histogram of its
// members' symbols by (H, bits(prev pos), L) for re-keying on a size change.
inline u32 triple(bool H, u64 pb, int L) { return (H ? 1u << 12 : 0u) | (static_cast<u32>(pb) << 6) | static_cast<u32>(L); }
constexpr u32 kDenseC = 1u << 22;

struct TabLive {
    std::vector<u32> pc;
    u32 CB = 128;
    u64 W = 0;
    std::vector<u64> bits;                                  // classes 1..CB: bitset over idx
    std::unordered_map<u32, std::vector<u32>> sv;           // classes > CB: sorted idx
    std::vector<u32> dn; std::unordered_map<u32, u32> sp;   // class sizes
    std::unordered_map<u32, std::vector<std::pair<u32, u32>>> sub;
    u32 csize(u32 c) const {
        if (c < kDenseC) return dn[c];
        auto it = sp.find(c); return it == sp.end() ? 0 : it->second;
    }
    void set_csize(u32 c, u32 s) {
        if (c < kDenseC) { dn[c] = s; return; }
        if (s == 0) sp.erase(c); else sp[c] = s;
    }
    // strictly below / above idx in live class c, -1 if none
    i64 pred(u32 c, u32 idx) const {
        if (csize(c) == 0) return -1;
        if (c <= CB) {
            const u64* b = bits.data() + static_cast<u64>(c - 1) * W;
            i64 w = idx >> 6;
            u64 x = b[w] & ((idx & 63) ? ((u64(1) << (idx & 63)) - 1) : 0);
            for (;;) {
                if (x) return w * 64 + 63 - __builtin_clzll(x);
                if (--w < 0) return -1;
                x = b[w];
            }
        }
        auto it = sv.find(c);
        if (it == sv.end()) return -1;
        const auto& v = it->second;
        auto p = std::lower_bound(v.begin(), v.end(), idx);
        if (p == v.begin()) return -1;
        return *(p - 1);
    }
    i64 succ(u32 c, u32 idx) const {
        if (csize(c) == 0) return -1;
        if (c <= CB) {
            const u64* b = bits.data() + static_cast<u64>(c - 1) * W;
            u64 w = idx >> 6;
            u64 x = b[w] & (((idx & 63) == 63) ? 0 : ~((u64(2) << (idx & 63)) - 1));
            for (;;) {
                if (x) return static_cast<i64>(w * 64 + __builtin_ctzll(x));
                if (++w >= W) return -1;
                x = b[w];
            }
        }
        auto it = sv.find(c);
        if (it == sv.end()) return -1;
        const auto& v = it->second;
        auto p = std::upper_bound(v.begin(), v.end(), idx);
        if (p == v.end()) return -1;
        return *p;
    }
    void mem_add(u32 c, u32 idx) {
        if (c <= CB) { bits[static_cast<u64>(c - 1) * W + (idx >> 6)] |= u64(1) << (idx & 63); return; }
        auto& v = sv[c];
        v.insert(std::lower_bound(v.begin(), v.end(), idx), idx);
    }
    void mem_del(u32 c, u32 idx) {
        if (c <= CB) { bits[static_cast<u64>(c - 1) * W + (idx >> 6)] &= ~(u64(1) << (idx & 63)); return; }
        auto it = sv.find(c);
        CHECK(it != sv.end(), "class %u missing", c);
        auto p = std::lower_bound(it->second.begin(), it->second.end(), idx);
        CHECK(p != it->second.end() && *p == idx, "member %u missing from class %u", idx, c);
        it->second.erase(p);
        if (it->second.empty()) sv.erase(it);
    }
    u32 subget(u32 c, u32 tr) const {
        auto it = sub.find(c);
        if (it == sub.end()) return 0;
        auto p = std::lower_bound(it->second.begin(), it->second.end(), std::make_pair(tr, 0u));
        return (p != it->second.end() && p->first == tr) ? p->second : 0;
    }
    void subadd(u32 c, u32 tr, i64 w) {
        auto& v = sub[c];
        auto p = std::lower_bound(v.begin(), v.end(), std::make_pair(tr, 0u));
        if (p != v.end() && p->first == tr) {
            CHECK(i64(p->second) + w >= 0, "sub-histogram underflow");
            p->second = static_cast<u32>(i64(p->second) + w);
            if (p->second == 0) v.erase(p);
        } else {
            CHECK(w > 0, "sub-histogram underflow (absent)");
            v.insert(p, std::make_pair(tr, static_cast<u32>(w)));
        }
        if (v.empty()) sub.erase(c);
    }
};

// ------------------------------------------------------------------- tree
constexpr int kMaxD = 14;
constexpr u32 kOV = 0x80000000u;
constexpr u32 kTermBit = 0x80000000u;

struct Tree {
    std::vector<u32> rk[kMaxD + 2], ct[kMaxD + 2], nc[kMaxD + 2], fc[kMaxD + 2], ol[kMaxD + 2];
    u32 off[kMaxD + 3] = {0};
    std::vector<u32> ork, oct, onc, ool;
    std::vector<u32> pool;
    u64 snap = 0;
    void clear() {
        for (int t = 0; t <= kMaxD + 1; ++t) { for (auto* v : {&rk[t], &ct[t], &nc[t], &fc[t], &ol[t]}) { std::vector<u32>().swap(*v); } }
        std::vector<u32>().swap(ork); std::vector<u32>().swap(oct); std::vector<u32>().swap(onc); std::vector<u32>().swap(ool);
        std::vector<u32>().swap(pool);
        snap = 0;
    }
    u64 overflow() const { return ork.size(); }
    u32& RK(u32 ref, int t) { return (ref & kOV) ? ork[ref ^ kOV] : rk[t][ref - off[t]]; }
    u32& CT(u32 ref, int t) { return (ref & kOV) ? oct[ref ^ kOV] : ct[t][ref - off[t]]; }
    u32& NC(u32 ref, int t) { return (ref & kOV) ? onc[ref ^ kOV] : nc[t][ref - off[t]]; }
    u32& OL(u32 ref, int t) { return (ref & kOV) ? ool[ref ^ kOV] : ol[t][ref - off[t]]; }
    u32 RK(u32 ref, int t) const { return (ref & kOV) ? ork[ref ^ kOV] : rk[t][ref - off[t]]; }
    u32 CT(u32 ref, int t) const { return (ref & kOV) ? oct[ref ^ kOV] : ct[t][ref - off[t]]; }
    u32 NC(u32 ref, int t) const { return (ref & kOV) ? onc[ref ^ kOV] : nc[t][ref - off[t]]; }
    u32 OL(u32 ref, int t) const { return (ref & kOV) ? ool[ref ^ kOV] : ol[t][ref - off[t]]; }
    int depth_of(u32 ref) const {   // snapshot refs only
        for (int t = 0; t <= kMaxD + 1; ++t) if (ref >= off[t] && ref < off[t + 1]) return t;
        return -1;
    }
    void range(u32 P, int t, u32& a, u32& b) const {
        if (P & kOV) { a = b = 0; return; }
        const u32 i = P - off[t];
        a = fc[t][i]; b = fc[t][i + 1];
    }
    static u32 lower_rank(const std::vector<u32>& v, u32 a, u32 b, u32 r) {
        while (a < b) { const u32 m = a + (b - a) / 2; if ((v[m] & kRankMask) < r) a = m + 1; else b = m; }
        return a;
    }
    const u32* olist(u32 P, int t) const { const u32 h = OL(P, t); return h ? pool.data() + h : nullptr; }
    static u32 ol_lower(const u32* L, u32 r) {
        u32 a = 0, b = L[0];
        while (a < b) { const u32 m = a + (b - a) / 2; if (L[2 + 2 * m] < r) a = m + 1; else b = m; }
        return a;
    }
    u32 child_ref(u32 P, int t, u32 r) const {
        if (t >= kMaxD) return kNone;
        u32 a, b; range(P, t, a, b);
        if (a < b) {
            const u32 i = lower_rank(rk[t + 1], a, b, r);
            if (i < b && (rk[t + 1][i] & kRankMask) == r) return off[t + 1] + i;
        }
        if (const u32* L = olist(P, t)) {
            const u32 i = ol_lower(L, r);
            if (i < L[0] && L[2 + 2 * i] == r) return kOV | L[3 + 2 * i];
        }
        return kNone;
    }
    // live (CT > 0) children of P in rank order: f(rank, ref)
    template <class F> void live_children(u32 P, int t, F&& f) const {
        if (t >= kMaxD) return;
        u32 a, b; range(P, t, a, b);
        const u32* L = olist(P, t);
        u32 i = a, j = 0;
        const u32 nl = L ? L[0] : 0;
        while (i < b || j < nl) {
            const u32 ra = i < b ? (rk[t + 1][i] & kRankMask) : 0xFFFFFFFFu;
            const u32 rb = j < nl ? L[2 + 2 * j] : 0xFFFFFFFFu;
            if (ra < rb) { const u32 x = off[t + 1] + i; if (ct[t + 1][i] > 0) f(ra, x); ++i; }
            else { const u32 x = kOV | L[3 + 2 * j]; if (oct[L[3 + 2 * j]] > 0) f(rb, x); ++j; }
        }
    }
    void ol_insert(u32 P, int t, u32 r, u32 j) {
        u32 h = OL(P, t);
        if (h == 0) {
            if (pool.empty()) pool.push_back(0);
            h = static_cast<u32>(pool.size());
            CHECK(pool.size() + 4 < 0xFFFFFFF0ull, "overflow pool full");
            pool.push_back(1); pool.push_back(1); pool.push_back(r); pool.push_back(j);
            OL(P, t) = h;
            return;
        }
        u32 sz = pool[h], cap = pool[h + 1];
        if (sz == cap) {
            const u32 nh = static_cast<u32>(pool.size());
            CHECK(pool.size() + 2 + 4ull * cap < 0xFFFFFFF0ull, "overflow pool full");
            pool.resize(pool.size() + 2 + 4ull * cap);
            pool[nh] = sz; pool[nh + 1] = 2 * cap;
            std::memcpy(pool.data() + nh + 2, pool.data() + h + 2, sizeof(u32) * 2 * sz);
            h = nh; OL(P, t) = h; cap *= 2;
        }
        const u32 i = ol_lower(pool.data() + h, r);
        u32* e = pool.data() + h + 2;
        std::memmove(e + 2 * (i + 1), e + 2 * i, sizeof(u32) * 2 * (sz - i));
        e[2 * i] = r; e[2 * i + 1] = j;
        pool[h] = sz + 1;
    }
};

// sibling cache of a node with many children (live state)
struct HugeCache {
    u32 ref = 0; int t = 0; bool dirty = true;
    std::vector<u32> rk, rf, q;   // live children (rank, ref); positions whose quotient crosses at R or R - 1
};

// ------------------------------------------------------------ global state
struct State {
    Source src;
    Universe U;
    Ranking R;
    Keys K;
    std::vector<u16> masks, dmask;
    Tree T;
    Store M;
    TabLive TL;
    std::vector<HugeCache> huge;
    std::unordered_map<u32, u32> huge_id;
    double live = 0;                       // live objective (bits)
    double p_ent[kNC] = {0, 0, 0, 0, 0}, p_rawt = 0, p_rawb = 0, p_esc = 0;
    i64 p_by[kNC] = {0, 0, 0, 0, 0}, p_csb = 0;
    u64 N = 0;
    double total_from_parts() const {
        double t = p_rawt + p_rawb + p_esc + 8.0 * double(p_csb);
        for (int c = 0; c < kNC; ++c) t += p_ent[c] + 8.0 * g_opt.tw * double(p_by[c]);
        return t;
    }
};
State* G = nullptr;

void build_cache(HugeCache& h) {
    const State& S = *G;
    h.rk.clear(); h.rf.clear(); h.q.clear();
    S.T.live_children(h.ref, h.t, [&](u32 r, u32 x) { h.rk.push_back(r); h.rf.push_back(x); });
    const u32 prank = h.ref == 0 ? 0 : (S.T.RK(h.ref, h.t) & kRankMask);
    const u32 n = static_cast<u32>(h.rk.size());
    for (u32 j = 0; j < n; ++j) {
        const u32 a = j ? h.rk[j - 1] : prank;
        const u64 X = S.K.X(a), R = n - j;
        if (quot_cross(X, R) || (R >= 2 && quot_cross(X, R - 1))) h.q.push_back(j);
    }
    h.dirty = false;
}
void register_huge(u32 ref, int t) {
    State& S = *G;
    if (S.huge_id.count(ref)) return;
    S.huge_id[ref] = static_cast<u32>(S.huge.size());
    S.huge.emplace_back();
    S.huge.back().ref = ref; S.huge.back().t = t;
    build_cache(S.huge.back());
}
const HugeCache* huge_of(u32 ref) {
    auto it = G->huge_id.find(ref);
    if (it == G->huge_id.end()) return nullptr;
    const HugeCache* h = &G->huge[it->second];
    CHECK(!h->dirty, "stale sibling cache");
    return h;
}

// ------------------------------------------------ overlay (virtual context)
struct Slot { u64 key; i64 a; i64 b; };
struct SMap {
    u32 cap = 0;
    std::vector<Slot> s;
    std::vector<u32> touched;
    explicit SMap(u32 c = 1u << 12) : cap(c), s(c, Slot{0, 0, 0}) {}
    u32 find(u64 key) const {
        u32 i = static_cast<u32>(fmix64(key)) & (cap - 1);
        while (s[i].key != 0) { if (s[i].key == key) return i; i = (i + 1) & (cap - 1); }
        return cap;
    }
    u32 slot(u64 key, bool& fresh) {
        u32 i = static_cast<u32>(fmix64(key)) & (cap - 1);
        u32 probes = 0;
        while (s[i].key != 0) {
            if (s[i].key == key) { fresh = false; return i; }
            i = (i + 1) & (cap - 1);
            CHECK(++probes < cap, "overlay map full");
        }
        fresh = true;
        return i;
    }
    void clear() { for (u32 i : touched) s[i] = Slot{0, 0, 0}; touched.clear(); }
};

struct VM { u64 ck; i64 mid; u64 T; u32 nsym; i32 fb; i32 sd; u32 hoff; u32 moff; u8 cls; };
struct VUndo { u32 vm; u32 bin; u32 oh; u32 onsym; i32 ofb; i32 osd; u64 oT; u64 omask; };
struct Sc { double cost; double ent[kNC]; double rawt, rawb, esc; i64 by[kNC], nm[kNC], csb; };

struct View {
    const u32* rk = nullptr; const u32* rf = nullptr; u32 n = 0;
    i32 D = -1; u32 A = 0; u32 Apos = 0; u32 prank = 0;
    const HugeCache* hc = nullptr;
    u32 vc() const { return n - (D >= 0 ? 1 : 0) + (A ? 1 : 0); }
    u32 vidx(u32 j) const { return j - ((D >= 0 && static_cast<u32>(D) < j) ? 1 : 0) + ((A && Apos <= j) ? 1 : 0); }
    u32 vidxA() const { return Apos - ((D >= 0 && static_cast<u32>(D) < Apos) ? 1 : 0); }
    u32 vprev(u32 j) const {
        i64 pj = static_cast<i64>(j) - 1;
        if (pj == D) --pj;
        const u32 lp = pj >= 0 ? rk[pj] : prank;
        if (A && Apos <= j && A > lp) return A;
        return lp;
    }
    u32 vprevA() const {
        i64 pj = static_cast<i64>(Apos) - 1;
        if (pj == D) --pj;
        return pj >= 0 ? rk[pj] : prank;
    }
    u32 lpos(u32 r) const { return static_cast<u32>(std::lower_bound(rk, rk + n, r) - rk); }
};

// node context: ref, depth, rank, P = product of the path's primes (root 1),
// parent ref / rank, live index in the parent's live children (lazy)
struct NCtx { u32 ref; int t; u32 rank; u128 P; u32 pref; u32 prank; i32 pj; };

struct VirtualCx {
    const State& S;
    enum { MNODE, MTERM, MVIEW, MTP, MCS, MSUB, MVMAP, NMAPS };
    SMap mp[NMAPS];
    struct Undo { u8 map; u32 slot; Slot old; };
    std::vector<Undo> undo;
    std::vector<VM> vms; std::vector<u32> vcnt; std::vector<u64> vmsk; std::vector<VUndo> vundo;
    Sc sc;
    SMap kmemo; std::vector<u32> kr, kf;
    struct Mark { std::size_t u, vu, nv, nc, nm; Sc sc; };

    explicit VirtualCx(const State& s)
        : S(s), mp{SMap(1u << 13), SMap(1u << 12), SMap(1u << 12), SMap(1u << 10), SMap(1u << 10), SMap(1u << 13), SMap(1u << 15)},
          kmemo(1u << 13) {
        undo.reserve(1 << 16); vundo.reserve(1 << 16); vms.reserve(1 << 12);
        vcnt.reserve(1u << 22); vmsk.reserve(1u << 16);
        kr.reserve(1u << 21); kf.reserve(1u << 21);
        reset();
    }
    void reset() {
        for (auto& x : mp) x.clear();
        undo.clear(); vundo.clear(); vms.clear(); vcnt.clear(); vmsk.clear();
        std::memset(&sc, 0, sizeof sc);
        kmemo.clear(); kr.clear(); kf.clear();
    }
    Mark mark() const { return {undo.size(), vundo.size(), vms.size(), vcnt.size(), vmsk.size(), sc}; }
    void restore(const Mark& k) {
        while (vundo.size() > k.vu) {
            const VUndo& e = vundo.back();
            VM& v = vms[e.vm];
            vcnt[v.hoff + e.bin] = e.oh; vmsk[v.moff + (e.bin >> 6)] = e.omask;
            v.T = e.oT; v.nsym = e.onsym; v.fb = e.ofb; v.sd = e.osd;
            vundo.pop_back();
        }
        vms.resize(k.nv); vcnt.resize(k.nc); vmsk.resize(k.nm);
        while (undo.size() > k.u) {
            const Undo& e = undo.back();
            mp[e.map].s[e.slot] = e.old;
            undo.pop_back();
        }
        sc = k.sc;
    }
    i64 geta(int map, u64 key) const { const u32 i = mp[map].find(key); return i == mp[map].cap ? 0 : mp[map].s[i].a; }
    const Slot* getslot(int map, u64 key) const { const u32 i = mp[map].find(key); return i == mp[map].cap ? nullptr : &mp[map].s[i]; }
    Slot& put(int map, u64 key) {
        bool fresh;
        SMap& m = mp[map];
        const u32 i = m.slot(key, fresh);
        undo.push_back({static_cast<u8>(map), i, m.s[i]});
        if (fresh) { m.s[i] = Slot{key, 0, 0}; m.touched.push_back(i); }
        return m.s[i];
    }

    // ---------------------------------------------------------- structure
    u32 eff_ct(u32 ref, int t) const { return static_cast<u32>(i64(S.T.CT(ref, t)) + geta(MNODE, u64(ref) + 1)); }
    void add_ct(u32 ref, int t, int d) { (void)t; put(MNODE, u64(ref) + 1).a += d; }
    u32 eff_nc(u32 ref, int t) const {
        const Slot* s = getslot(MNODE, u64(ref) + 1);
        return static_cast<u32>(i64(S.T.NC(ref, t)) + (s ? s->b : 0));
    }
    void set_nc(u32 ref, int t, u32 v) { Slot& s = put(MNODE, u64(ref) + 1); s.b = i64(v) - i64(S.T.NC(ref, t)); }
    bool eff_term(u32 ref, int t) const {
        const Slot* s = getslot(MTERM, u64(ref) + 1);
        return s ? s->a != 0 : (S.T.RK(ref, t) & kTermBit) != 0;
    }
    void set_term(u32 ref, int t, bool v) { (void)t; put(MTERM, u64(ref) + 1).a = v ? 1 : 0; }
    u32 child_ref(u32 P, int t, u32 r) const { return S.T.child_ref(P, t, r); }

    // live children of v: sibling cache (huge) or a per-record memoized scan
    void live_kids(u32 v, int t, View& V) {
        if (const HugeCache* h = huge_of(v)) { V.rk = h->rk.data(); V.rf = h->rf.data(); V.n = static_cast<u32>(h->rk.size()); V.hc = h; return; }
        V.hc = nullptr;
        bool fresh;
        const u32 sl = kmemo.slot(u64(v) + 1, fresh);
        if (fresh) {
            kmemo.s[sl] = Slot{u64(v) + 1, static_cast<i64>(kr.size()), 0};
            kmemo.touched.push_back(sl);
            const std::size_t c0 = kr.size();
            S.T.live_children(v, t, [&](u32 r, u32 x) {
                CHECK(kr.size() < kr.capacity(), "child scan buffer full");
                kr.push_back(r); kf.push_back(x);
            });
            kmemo.s[sl].b = static_cast<i64>(kr.size() - c0);
        }
        V.rk = kr.data() + kmemo.s[sl].a; V.rf = kf.data() + kmemo.s[sl].a; V.n = static_cast<u32>(kmemo.s[sl].b);
    }
    View view(u32 v, int t, u32 prank) {
        View V;
        live_kids(v, t, V);
        V.prank = prank;
        if (const Slot* s = getslot(MVIEW, u64(v) + 1)) {
            V.D = static_cast<i32>(s->a & 0xFFFFFFFF) - 1;
            V.A = static_cast<u32>(s->b & 0xFFFFFFFF); V.Apos = static_cast<u32>(s->b >> 32);
        }
        return V;
    }
    void set_view(u32 v, int t, i32 D, u32 A, u32 Apos) {
        Slot& s = put(MVIEW, u64(v) + 1);
        s.a = static_cast<i64>(D + 1) | (static_cast<i64>(t) << 40);
        s.b = static_cast<i64>(A) | (static_cast<i64>(Apos) << 32);
    }

    // ------------------------------------------------------------- models
    u32 vm_get(u64 ck) {
        bool fresh;
        SMap& m = mp[MVMAP];
        const u32 sl = m.slot(ck + 1, fresh);
        if (!fresh) return static_cast<u32>(m.s[sl].a);
        undo.push_back({static_cast<u8>(MVMAP), sl, m.s[sl]});
        m.s[sl] = Slot{ck + 1, static_cast<i64>(vms.size()), 0};
        m.touched.push_back(sl);
        VM v;
        v.ck = ck; v.cls = static_cast<u8>(ck_cls(ck)); v.mid = S.M.find(ck);
        const u32 A = kAl[v.cls], W = kAW[v.cls];
        v.hoff = static_cast<u32>(vcnt.size()); v.moff = static_cast<u32>(vmsk.size());
        vcnt.resize(vcnt.size() + A); vmsk.resize(vmsk.size() + W);
        if (v.mid >= 0) {
            const Mdl& L = S.M.md[v.mid];
            std::memcpy(vcnt.data() + v.hoff, S.M.H.data() + L.off, sizeof(u32) * A);
            std::memcpy(vmsk.data() + v.moff, S.M.MK.data() + L.moff, sizeof(u64) * W);
            v.T = L.T; v.nsym = L.nsym; v.fb = L.fb; v.sd = L.sd;
        } else {
            std::memset(vcnt.data() + v.hoff, 0, sizeof(u32) * A);
            std::memset(vmsk.data() + v.moff, 0, sizeof(u64) * W);
            v.T = 0; v.nsym = 0; v.fb = 0; v.sd = 0;
        }
        vms.push_back(v);
        return static_cast<u32>(vms.size() - 1);
    }
    bool live_exists(const VM& v) const { return v.mid >= 0 && S.M.md[v.mid].T > 0; }
    bool vkey_deleted(int cls, u64 key) const {
        const u32 i = mp[MVMAP].find(mkck(cls, key) + 1);
        if (i == mp[MVMAP].cap) return false;
        return vms[static_cast<std::size_t>(mp[MVMAP].s[i].a)].T == 0;
    }
    // VByte cost of inserting key into the virtual key set of cls (key itself excluded)
    int key_ins(int cls, u64 key) const {
        const std::vector<u64>& K = S.M.keys[cls];
        u64 p = 0; bool hp = false, hq = false; u64 q = 0;
        auto it = std::lower_bound(K.begin(), K.end(), key);
        for (auto jt = it; jt != K.begin();) { --jt; if (*jt == key) continue; if (!vkey_deleted(cls, *jt)) { p = *jt; hp = true; break; } }
        for (auto jt = it; jt != K.end(); ++jt) { if (*jt == key) continue; if (!vkey_deleted(cls, *jt)) { q = *jt; hq = true; break; } }
        for (const VM& v : vms) {
            if (v.cls != cls || v.T == 0 || live_exists(v)) continue;
            const u64 k2 = ck_key(v.ck);
            if (k2 < key && (!hp || k2 > p)) { p = k2; hp = true; }
            if (k2 > key && (!hq || k2 < q)) { q = k2; hq = true; }
        }
        const u64 pp = hp ? p : 0;
        return vbl(key - pp) + (hq ? vbl(q - key) - vbl(q - pp) : 0);
    }
    void sym(int cls, u64 key, u32 s, i64 w) {
        const u64 ck = mkck(cls, key);
        const u32 vi = vm_get(ck);
        VM& v = vms[vi];
        u32* h = vcnt.data() + v.hoff;
        u64* mk = vmsk.data() + v.moff;
        const u64 h0 = h[s], T0 = v.T;
        CHECK(w > 0 || h0 >= static_cast<u64>(-w), "virtual symbol underflow: class %d key %llu symbol %u (h %llu, w %lld)",
              cls, (unsigned long long)key, s, (unsigned long long)h0, (long long)w);
        vundo.push_back({vi, s, static_cast<u32>(h0), v.nsym, v.fb, v.sd, T0, mk[s >> 6]});
        const u64 h1 = static_cast<u64>(i64(h0) + w), T1 = static_cast<u64>(i64(T0) + w);
        double de;
        if (w == 1) de = df(T0) - df(h0);
        else if (w == -1) de = df(h0 - 1) - df(T0 - 1);
        else de = (fx(T1) - fx(T0)) - (fx(h1) - fx(h0));
        const i64 mb0 = T0 ? vbl(v.nsym) + v.sd + v.fb : 0;
        if (h0 == 0 && h1 > 0) { v.sd += sd_of(cls, mk, s); mk[s >> 6] |= u64(1) << (s & 63); ++v.nsym; }
        else if (h0 > 0 && h1 == 0) { mk[s >> 6] &= ~(u64(1) << (s & 63)); v.sd -= sd_of(cls, mk, s); --v.nsym; }
        h[s] = static_cast<u32>(h1);
        v.T = T1;
        if (th1(T0) == th1(T1) && th2(T0) == th2(T1)) v.fb += blen(h1, T1) - blen(h0, T0);
        else v.fb = fb_all(h, mk, kAW[cls], T1);
        const i64 mb1 = T1 ? vbl(v.nsym) + v.sd + v.fb : 0;
        i64 dby = mb1 - mb0;
        if (T0 == 0 && T1 > 0) {
            const u64 nm = static_cast<u64>(S.M.nmod[cls] + sc.nm[cls]);
            dby += key_ins(cls, ck_key(ck)) + vbl(nm + 1) - vbl(nm);
            ++sc.nm[cls];
        } else if (T0 > 0 && T1 == 0) {
            const u64 nm = static_cast<u64>(S.M.nmod[cls] + sc.nm[cls]);
            dby += -key_ins(cls, ck_key(ck)) + vbl(nm - 1) - vbl(nm);
            --sc.nm[cls];
        }
        sc.ent[cls] += de; sc.by[cls] += dby;
        sc.cost += de + 8.0 * g_opt.tw * double(dby);
    }
    void rawt(i64 d) { sc.rawt += double(d); sc.cost += double(d); }
    void rawb(i64 d) { sc.rawb += double(d); sc.cost += double(d); }
    void esc(i64 d) { sc.esc += double(d); sc.cost += double(d); }
    void csb(i64 d) { sc.csb += d; sc.cost += 8.0 * double(d); }
    double cost() const { return sc.cost; }

    // count symbol (+ tail + escape) of a node with c children
    void count_sym(int d, u64 key, u64 c, int w) {
        sym(C_CNT, key, static_cast<u32>(std::min<u64>(c, 63)), w);
        if (c >= 63) sym(C_TAIL, static_cast<u64>(d), static_cast<u32>(std::min<u64>(c, 4095)), w);
        if (c >= 4095) { const u64 v = c - 4094; esc(i64(w) * (2 * bits_of(v) - 1)); }
    }

    // ---------------------------------------------------------- tree keys
    u64 count_key_of(NCtx& v) {
        if (v.t == 0) return S.K.count(0, 0, 1, false, 0);
        const View Vp = view(v.pref, v.t - 1, v.prank);
        if (v.pj < 0) {
            v.pj = static_cast<i32>(Vp.lpos(v.rank));
            CHECK(static_cast<u32>(v.pj) < Vp.n && Vp.rk[v.pj] == v.rank, "node rank %u not among its parent's live children", v.rank);
        }
        CHECK(Vp.D != v.pj, "count key of a detached node");
        const u32 i = Vp.vidx(static_cast<u32>(v.pj)), c = Vp.vc();
        return S.K.count(v.t, v.rank, v.P, i == 0, c - i - 1);
    }
    // v's own count value (and flag) move from c0 to c1 children
    void node_count_change(NCtx& v, u32 c0, u32 c1) {
        const u64 key = count_key_of(v);
        if (std::min<u64>(c0, 63) != std::min<u64>(c1, 63) || c0 >= 63 || c1 >= 63) {
            count_sym(v.t, key, c0, -1);
            count_sym(v.t, key, c1, +1);
        }
        if (v.t >= 1) {
            const bool term = eff_term(v.ref, v.t);
            const bool f0 = c0 > 0, f1 = c1 > 0;
            const u64 k0 = f0 ? S.K.flag(v.t, c0) : 0, k1 = f1 ? S.K.flag(v.t, c1) : 0;
            if (!(f0 && f1 && k0 == k1)) {
                if (f0) sym(C_FLAG, k0, term ? 1 : 0, -1);
                if (f1) sym(C_FLAG, k1, term ? 1 : 0, +1);
            }
        }
    }
    // gap and count keys of the surviving live child j of v in view V
    void child_keys(const View& V, u32 j, int t, u128 Pv, u64& gk, int& L, u64& kk) const {
        const u32 c = V.vc(), i = V.vidx(j), a = V.vprev(j), r = V.rk[j];
        gk = S.K.gap(t + 1, c, i == 0, a, c - i, Pv);
        L = bits_of(r - a);
        kk = S.K.count(t + 1, r, Pv * S.R.pv[r], i == 0, c - i - 1);
    }
    void rekey_child(const NCtx& v, const View& V0, const View& V1, u32 j) {
        u64 g0, g1, k0, k1; int L0, L1;
        child_keys(V0, j, v.t, v.P, g0, L0, k0);
        child_keys(V1, j, v.t, v.P, g1, L1, k1);
        if (g0 != g1 || L0 != L1) {
            sym(C_GAP, g0, static_cast<u32>(L0), -1); rawt(-(L0 - 1));
            sym(C_GAP, g1, static_cast<u32>(L1), +1); rawt(L1 - 1);
        }
        if (k0 != k1) {
            const u32 cx = eff_nc(V0.rf[j], v.t + 1);
            const u32 s = static_cast<u32>(std::min<u64>(cx, 63));
            sym(C_CNT, k0, s, -1);
            sym(C_CNT, k1, s, +1);
        }
    }
    // cap(R, 13) / cap(A, 13) / bits((m - a) / R) boundary between Rm and Rm + 1
    bool crosses(u64 Rm, u32 a) const {
        if (Rm == 0) return false;
        if (pow2(Rm + 1) && Rm + 1 <= 8192) return true;
        if (pow2(Rm) && Rm <= 8192) return true;
        return quot_cross(S.K.X(a), Rm);
    }
    std::vector<u32> cand;
    void rekey(const NCtx& v, const View& V0, const View& V1, u32 lim, bool all) {
        const u32 n = V0.n;
        auto surv = [&](u32 j) { return static_cast<i32>(j) != V0.D && static_cast<i32>(j) != V1.D; };
        if (all) { for (u32 j = 0; j < n; ++j) if (surv(j)) rekey_child(v, V0, V1, j); return; }
        cand.clear();
        auto next_surv = [&](u32 j) { while (j < n && !surv(j)) ++j; return j; };
        if (V1.D >= 0) { const u32 j = next_surv(static_cast<u32>(V1.D) + 1); if (j < n) cand.push_back(j); }
        if (V0.D >= 0) { const u32 j = next_surv(static_cast<u32>(V0.D) + 1); if (j < n) cand.push_back(j); }
        if (V1.A) { const u32 j = next_surv(V1.Apos); if (j < n) cand.push_back(j); }
        if (V0.A) { const u32 j = next_surv(V0.Apos); if (j < n) cand.push_back(j); }
        if (V0.hc) {
            for (u32 q : V0.hc->q) { if (q >= lim) break; if (surv(q)) cand.push_back(q); }
            for (int k = 0; k <= 13; ++k) {
                const u64 p = u64(1) << k;
                for (u64 R : {p - 1, p, p + 1}) {
                    if (R == 0 || R > n) continue;
                    const u32 j = static_cast<u32>(n - R);
                    if (j < lim && surv(j)) cand.push_back(j);
                }
            }
        } else {
            const u32 c0 = V0.vc(), c1 = V1.vc();
            for (u32 j = 0; j < lim && j < n; ++j) {
                if (!surv(j)) continue;
                const u64 R0 = c0 - V0.vidx(j), R1 = c1 - V1.vidx(j);
                if (R0 == R1) continue;
                if (crosses(std::min(R0, R1), V0.vprev(j))) cand.push_back(j);
            }
        }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        // copy: rekey_child may recurse into nothing that touches cand, but be safe
        const std::vector<u32> cc = cand;
        for (u32 j : cc) rekey_child(v, V0, V1, j);
    }
    // attach a new child of rank ru under v; returns its virtual index and c(v) after
    void attach(NCtx& v, u32 ru, u32& ia, u32& c1out) {
        const View V0 = view(v.ref, v.t, v.rank);
        CHECK(V0.A == 0, "second attach at one node");
        View V1 = V0;
        V1.A = ru; V1.Apos = V0.lpos(ru);
        const u32 c0 = V0.vc(), c1 = V1.vc();
        CHECK(c0 == eff_nc(v.ref, v.t), "view child count %u != nc %u", c0, eff_nc(v.ref, v.t));
        CHECK(ru > v.rank, "child rank not above parent");
        node_count_change(v, c0, c1);
        rekey(v, V0, V1, V1.Apos, c0 > 0 && bucket_of(c0) != bucket_of(c1));
        ia = V1.vidxA();
        const u32 a = V1.vprevA();
        const int L = bits_of(ru - a);
        sym(C_GAP, S.K.gap(v.t + 1, c1, ia == 0, a, c1 - ia, v.P), static_cast<u32>(L), +1);
        rawt(L - 1);
        set_view(v.ref, v.t, V1.D, ru, V1.Apos);
        set_nc(v.ref, v.t, c1);
        c1out = c1;
    }
    // detach v's child x (rank xr; its own children already gone)
    void detach(NCtx& v, u32 xr, u32 xref) {
        const View V0 = view(v.ref, v.t, v.rank);
        CHECK(V0.A == 0 && V0.D < 0, "detach at a node with a virtual change");
        const u32 D = V0.lpos(xr);
        CHECK(D < V0.n && V0.rk[D] == xr && V0.rf[D] == xref, "detached child not live");
        View V1 = V0; V1.D = static_cast<i32>(D);
        const u32 c0 = V0.vc(), c1 = V1.vc();
        CHECK(c0 == eff_nc(v.ref, v.t), "view child count mismatch (detach)");
        const u32 i = V0.vidx(D), a = V0.vprev(D);
        const int L = bits_of(xr - a);
        sym(C_GAP, S.K.gap(v.t + 1, c0, i == 0, a, c0 - i, v.P), static_cast<u32>(L), -1);
        rawt(-(L - 1));
        const u32 cx = eff_nc(xref, v.t + 1);
        CHECK(cx == 0, "detached child still has children");
        count_sym(v.t + 1, S.K.count(v.t + 1, xr, v.P * S.R.pv[xr], i == 0, c0 - i - 1), cx, -1);
        node_count_change(v, c0, c1);
        rekey(v, V0, V1, D, c1 > 0 && bucket_of(c0) != bucket_of(c1));
        set_view(v.ref, v.t, static_cast<i32>(D), 0, 0);
        set_nc(v.ref, v.t, c1);
    }

    // --------------------------------------------------------- rank table
    u32 vpc(u32 i) const { return static_cast<u32>(i64(S.TL.pc[i]) + geta(MTP, u64(i) + 1)); }
    void set_pc(u32 i, u32 v) { put(MTP, u64(i) + 1).a = i64(v) - i64(S.TL.pc[i]); }
    u32 vcs(u32 c) const { return static_cast<u32>(i64(S.TL.csize(c)) + geta(MCS, u64(c) + 1)); }
    void add_cs(u32 c, int d) { put(MCS, u64(c) + 1).a += d; }
    void sub_add(u32 c, u32 tr, i64 w) { put(MSUB, ((u64(c) << 16) | tr) + 1).a += w; }
    i64 vpred(u32 c, u32 idx) const {
        i64 best = -1;
        u32 x = idx;
        for (;;) { const i64 q = S.TL.pred(c, x); if (q < 0) break; if (vpc(static_cast<u32>(q)) == c) { best = q; break; } x = static_cast<u32>(q); }
        const SMap& m = mp[MTP];
        for (u32 sl : m.touched) {
            const Slot& s = m.s[sl];
            if (s.key == 0 || s.a == 0) continue;
            const u32 i = static_cast<u32>(s.key - 1);
            if (i < idx && static_cast<i64>(i) > best && S.TL.pc[i] != c && vpc(i) == c) best = i;
        }
        return best;
    }
    i64 vsucc(u32 c, u32 idx) const {
        i64 best = -1;
        u32 x = idx;
        for (;;) { const i64 q = S.TL.succ(c, x); if (q < 0) break; if (vpc(static_cast<u32>(q)) == c) { best = q; break; } x = static_cast<u32>(q); }
        const SMap& m = mp[MTP];
        for (u32 sl : m.touched) {
            const Slot& s = m.s[sl];
            if (s.key == 0 || s.a == 0) continue;
            const u32 i = static_cast<u32>(s.key - 1);
            if (i > idx && (best < 0 || static_cast<i64>(i) < best) && S.TL.pc[i] != c && vpc(i) == c) best = i;
        }
        return best;
    }
    // member p's gap symbol (previous member q, -1 = first) in class c of size s
    void tab_sym(u32 c, u32 s, i64 q, u32 p, i64 w) {
        const u64 pp = q >= 0 ? S.U.pos[q] : 0;
        const bool H = q >= 0 && S.U.hi[q];
        const u64 g = q >= 0 ? S.U.pos[p] - pp : S.U.pos[p] + 1;
        const int L = bits_of(g);
        const u64 pb = capb(pp, 41);
        sym(C_TAB, S.K.tablek(s, H, pb), static_cast<u32>(L), w);
        rawb(w * (L - 1));
        sub_add(c, triple(H, pb, L), w);
    }
    void rekey_class(u32 c, u32 s0, u32 s1) {
        if (capb(s0, 40) == capb(s1, 40)) return;
        std::vector<std::pair<u32, i64>> ent;
        auto it = S.TL.sub.find(c);
        if (it != S.TL.sub.end()) for (const auto& e : it->second) ent.emplace_back(e.first, i64(e.second) + geta(MSUB, ((u64(c) << 16) | e.first) + 1));
        const SMap& m = mp[MSUB];
        for (u32 sl : m.touched) {
            const Slot& sx = m.s[sl];
            if (sx.key == 0 || sx.a == 0) continue;
            const u64 k = sx.key - 1;
            if ((k >> 16) != c) continue;
            const u32 tr = static_cast<u32>(k & 0xFFFF);
            if (S.TL.subget(c, tr) == 0) ent.emplace_back(tr, sx.a);
        }
        std::sort(ent.begin(), ent.end());
        ent.erase(std::unique(ent.begin(), ent.end()), ent.end());
        for (const auto& e : ent) {
            CHECK(e.second >= 0, "negative class sub-histogram");
            if (e.second == 0) continue;
            const bool H = (e.first >> 12) & 1; const u64 pb = (e.first >> 6) & 63; const u32 L = e.first & 63;
            sym(C_TAB, S.K.tablek(s0, H, pb), L, -e.second);
            sym(C_TAB, S.K.tablek(s1, H, pb), L, +e.second);
        }
    }
    void class_remove(u32 c, u32 p) {
        const u32 s = vcs(c);
        CHECK(s >= 1, "remove from empty class %u", c);
        const i64 q = vpred(c, p), t = vsucc(c, p);
        tab_sym(c, s, q, p, -1);
        if (t >= 0) tab_sym(c, s, static_cast<i64>(p), static_cast<u32>(t), -1);
        add_cs(c, -1);
        if (s - 1 > 0) rekey_class(c, s, s - 1);
        if (t >= 0) tab_sym(c, s - 1, q, static_cast<u32>(t), +1);
        csb(-vbl(s) + (s - 1 > 0 ? vbl(s - 1) : 0));
    }
    void class_insert(u32 c, u32 p) {
        const u32 s = vcs(c);
        const i64 q = vpred(c, p), t = vsucc(c, p);
        if (t >= 0) tab_sym(c, s, q, static_cast<u32>(t), -1);
        add_cs(c, +1);
        if (s > 0) rekey_class(c, s, s + 1);
        tab_sym(c, s + 1, q, p, +1);
        if (t >= 0) tab_sym(c, s + 1, static_cast<i64>(p), static_cast<u32>(t), +1);
        csb((s > 0 ? -vbl(s) : 0) + vbl(s + 1));
    }
    void tbl(u32 idx, int d) {
        const u32 a = vpc(idx);
        CHECK(d > 0 || a > 0, "prime count underflow (idx %u)", idx);
        const u32 b = static_cast<u32>(i64(a) + d);
        if (a > 0) class_remove(a, idx);
        set_pc(idx, b);
        if (b > 0) class_insert(b, idx);
    }
};

// ---------------------------------------------------------------- per record
struct Prep {
    u128 n; int k;
    u64 q[kMaxK]; u32 idx[kMaxK], rnk[kMaxK];
    int ord[kMaxK];
    u128 sp[kMaxK + 1], sl[kMaxK + 1], lam_all;
    int L; u32 prk[kMaxD + 1], pidx[kMaxD + 1], ref[kMaxD + 1];
    u128 pp[kMaxD + 1];
};
void prep(const Rec& r, u16 mask, Prep& P) {
    if (r.k < 1 || r.k > kMaxK) __builtin_unreachable();
    P.n = r.n; P.k = r.k;
    for (int j = 0; j < r.k; ++j) {
        P.q[j] = r.p[j];
        const u32 i = G->U.idx(r.p[j]);
        CHECK(i != kNone, "prime %llu not in universe", (unsigned long long)r.p[j]);
        P.idx[j] = i; P.rnk[j] = G->R.rank[i];
        P.ord[j] = j;
    }
    std::sort(P.ord, P.ord + r.k, [&](int a, int b) { return P.rnk[a] < P.rnk[b]; });
    P.sp[r.k] = 1; P.sl[r.k] = 1;
    for (int jj = r.k - 1; jj >= 0; --jj) {
        const u64 q = P.q[P.ord[jj]];
        P.sp[jj] = P.sp[jj + 1] * q;
        P.sl[jj] = lcm_step(P.sl[jj + 1], q - 1);
    }
    P.lam_all = P.sl[0];
    P.L = 0; P.ref[0] = 0; P.pp[0] = 1; P.prk[0] = 0;
    for (int jj = 0; jj < r.k; ++jj) {
        const int pos = P.ord[jj];
        if (!(mask >> pos & 1)) continue;
        ++P.L;
        P.prk[P.L] = P.rnk[pos]; P.pidx[P.L - 1] = P.idx[pos];
        P.pp[P.L] = P.pp[P.L - 1] * P.q[pos];
    }
    CHECK(P.L >= 1, "empty path");
}
void walk_path(const Tree& T, Prep& P) {
    if (P.L < 1 || P.L > kMaxD) __builtin_unreachable();
    for (int t = 1; t <= P.L; ++t) {
        P.ref[t] = T.child_ref(P.ref[t - 1], t - 1, P.prk[t]);
        CHECK(P.ref[t] != kNone && T.CT(P.ref[t], t) > 0, "current path missing from the tree (depth %d)", t);
    }
}
inline NCtx path_ctx(const Prep& P, int t) {
    NCtx c;
    c.ref = P.ref[t]; c.t = t; c.rank = t ? P.prk[t] : 0; c.P = P.pp[t];
    c.pref = t ? P.ref[t - 1] : 0; c.prank = t >= 2 ? P.prk[t - 1] : 0; c.pj = -1;
    return c;
}

void remove_path(VirtualCx& vx, const Prep& P) {
    const int L = P.L;
    const u32 E = P.ref[L];
    CHECK(vx.eff_term(E, L), "path end not terminal");
    vx.set_term(E, L, false);
    const u32 cE = vx.eff_nc(E, L);
    if (cE > 0) { const u64 fk = vx.S.K.flag(L, cE); vx.sym(C_FLAG, fk, 1, -1); vx.sym(C_FLAG, fk, 0, +1); }
    else CHECK(vx.eff_ct(E, L) == 1, "terminal leaf shared by two paths");
    for (int i = 0; i < L; ++i) vx.tbl(P.pidx[i], -1);
    for (int t = L; t >= 1; --t) {
        const u32 x = P.ref[t];
        const u32 c = vx.eff_ct(x, t);
        CHECK(c >= 1, "path count underflow");
        vx.add_ct(x, t, -1);
        if (c == 1) {
            CHECK(vx.eff_nc(x, t) == 0, "dying node has live children");
            NCtx v = path_ctx(P, t - 1);
            vx.detach(v, P.prk[t], x);
        }
    }
    vx.add_ct(0, 0, -1);
}

struct Tip { u32 rank; int t; u128 P; bool first; u32 after; };
inline u64 tip_key(const Keys& K, const Tip& w) { return K.count(w.t, w.rank, w.P, w.first, w.after); }

// straight-line insertion of a given path (apply / re-cost), same ops as the DFS
void insert_virtual(VirtualCx& vx, const Prep& P) {
    const Keys& K = vx.S.K;
    const int L = P.L;
    for (int i = 0; i < L; ++i) vx.tbl(P.pidx[i], +1);
    vx.add_ct(0, 0, +1);
    NCtx cur = path_ctx(P, 0);
    bool alive = true;
    Tip tip{};
    for (int t = 1; t <= L; ++t) {
        if (alive) {
            const u32 ex = vx.child_ref(cur.ref, t - 1, P.prk[t]);
            if (ex != kNone && vx.eff_ct(ex, t) > 0) {
                vx.add_ct(ex, t, +1);
                NCtx nx; nx.ref = ex; nx.t = t; nx.rank = P.prk[t]; nx.P = P.pp[t]; nx.pref = cur.ref; nx.prank = cur.rank; nx.pj = -1;
                cur = nx;
                continue;
            }
            u32 ia, c1;
            vx.attach(cur, P.prk[t], ia, c1);
            tip = Tip{P.prk[t], t, P.pp[t], ia == 0, c1 - ia - 1};
            alive = false;
            continue;
        }
        vx.count_sym(tip.t, tip_key(K, tip), 1, +1);
        vx.sym(C_FLAG, K.flag(tip.t, 1), 0, +1);
        const int Lg = bits_of(P.prk[t] - tip.rank);
        vx.sym(C_GAP, K.gap(t, 1, true, tip.rank, 1, tip.P), static_cast<u32>(Lg), +1);
        vx.rawt(Lg - 1);
        tip = Tip{P.prk[t], t, P.pp[t], true, 0};
    }
    if (alive) {
        CHECK(!vx.eff_term(cur.ref, L), "insert ends at another number's terminal");
        const u32 c = vx.eff_nc(cur.ref, L);
        CHECK(c > 0, "live non-terminal leaf");
        const u64 fk = K.flag(L, c);
        vx.sym(C_FLAG, fk, 0, -1); vx.sym(C_FLAG, fk, 1, +1);
        vx.set_term(cur.ref, L, true);
    } else {
        vx.count_sym(tip.t, tip_key(K, tip), 0, +1);
    }
}

constexpr u32 kNoMask = 0xFFFFFFFFu;

struct Evaluator {
    VirtualCx vx;
    Prep P;
    double best = 0; u32 best_mask = kNoMask; u32 old_mask = 0;
    bool old_seen = false; double old_total = 0;
    double tlb[kMaxK + 1];
    u64 steps = 0, steps0 = 0, capped = 0, old_bad = 0;
    bool cap_hit = false;
    std::vector<std::pair<u32, double>>* dumpc = nullptr;
    Evaluator() : vx(*G) {}

    inline bool reachable(u128 prod, u128 lam, int jj) const {
        u128 Lm;
        if (__builtin_mul_overflow(lam, P.sl[jj], &Lm) || Lm > P.lam_all) Lm = P.lam_all;
        u128 x, y;
        if (__builtin_mul_overflow(prod, P.sp[jj], &x)) return true;
        if (__builtin_mul_overflow(x, Lm, &y)) return true;
        return y > P.n;
    }
    inline void consider(double total, u32 mask) {
        if (dbg) std::printf("    DBG consider mask %u total %.6f\n", mask, total);
        if (dumpc) dumpc->emplace_back(mask, total);
        if (mask == old_mask) { old_seen = true; old_total = total; return; }
        if (total < best - 1e-9) { best = total; best_mask = mask; }
        else if (total <= best + 1e-9 && best_mask != kNoMask && mask < best_mask) { best_mask = mask; if (total < best) best = total; }
    }
    bool dbg = false;
    inline bool prune(double cost, int jn) const {
        const bool pr = !g_opt.no_prune && cost + tlb[jn] > best + 1e-9;
        if (dbg) std::printf("    DBG prune? cost %.6f tlb[%d] %.6f best %.6f -> %s\n", cost, jn, tlb[jn], best, pr ? "PRUNE" : "go");
        return pr;
    }

    void dfs_chain(const Tip& tip, int j0, u128 prod, u128 lam, u32 mask) {
        const Keys& K = G->K;
        if (tip.t >= kMaxD) return;
        const u64 tk = tip_key(K, tip);
        for (int jj = j0; jj < P.k; ++jj) {
            if (!reachable(prod, lam, jj)) break;
            if (++steps - steps0 > g_opt.step_cap) { cap_hit = true; return; }
            const int pos = P.ord[jj];
            const u64 q = P.q[pos];
            const u128 prod1 = prod * q, lam1 = lcm_step(lam, q - 1);
            const u32 mask1 = mask | (1u << pos);
            const u32 rn = P.rnk[pos];
            const VirtualCx::Mark m = vx.mark();
            vx.tbl(P.idx[pos], +1);
            vx.count_sym(tip.t, tk, 1, +1);
            vx.sym(C_FLAG, K.flag(tip.t, 1), 0, +1);
            const int Lg = bits_of(rn - tip.rank);
            vx.sym(C_GAP, K.gap(tip.t + 1, 1, true, tip.rank, 1, tip.P), static_cast<u32>(Lg), +1);
            vx.rawt(Lg - 1);
            const Tip w{rn, tip.t + 1, prod1, true, 0};
            if (valid_pl(prod1, lam1, P.n)) {
                const VirtualCx::Mark m2 = vx.mark();
                vx.count_sym(w.t, tip_key(K, w), 0, +1);
                consider(vx.cost(), mask1);
                vx.restore(m2);
            }
            if (!prune(vx.cost(), jj + 1)) dfs_chain(w, jj + 1, prod1, lam1, mask1);
            vx.restore(m);
            if (cap_hit) return;
        }
    }
    void dfs_alive(NCtx& Pn, int j0, u128 prod, u128 lam, u32 mask) {
        const Keys& K = G->K;
        const int t = Pn.t;
        for (int jj = j0; jj < P.k; ++jj) {
            if (!reachable(prod, lam, jj)) break;
            if (++steps - steps0 > g_opt.step_cap) { cap_hit = true; return; }
            const int pos = P.ord[jj];
            const u64 q = P.q[pos];
            const u128 prod1 = prod * q, lam1 = lcm_step(lam, q - 1);
            const u32 mask1 = mask | (1u << pos);
            const u32 rn = P.rnk[pos];
            const bool v1 = valid_pl(prod1, lam1, P.n);
            if (t >= kMaxD) break;
            const VirtualCx::Mark m = vx.mark();
            vx.tbl(P.idx[pos], +1);
            const u32 ex = vx.child_ref(Pn.ref, t, rn);
            if (ex != kNone && vx.eff_ct(ex, t + 1) > 0) {
                vx.add_ct(ex, t + 1, +1);
                if (v1) {
                    const VirtualCx::Mark m2 = vx.mark();
                    CHECK(!vx.eff_term(ex, t + 1), "a valid S is another number's path (depth %d)", t + 1);
                    const u32 c = vx.eff_nc(ex, t + 1);
                    CHECK(c > 0, "live non-terminal leaf in search");
                    const u64 fk = K.flag(t + 1, c);
                    vx.sym(C_FLAG, fk, 0, -1); vx.sym(C_FLAG, fk, 1, +1);
                    consider(vx.cost(), mask1);
                    vx.restore(m2);
                }
                NCtx X; X.ref = ex; X.t = t + 1; X.rank = rn; X.P = prod1; X.pref = Pn.ref; X.prank = Pn.rank; X.pj = -1;
                dfs_alive(X, jj + 1, prod1, lam1, mask1);
            } else {
                u32 ia, c1;
                vx.attach(Pn, rn, ia, c1);
                const Tip tip{rn, t + 1, prod1, ia == 0, c1 - ia - 1};
                if (v1) {
                    const VirtualCx::Mark m2 = vx.mark();
                    vx.count_sym(tip.t, tip_key(K, tip), 0, +1);
                    consider(vx.cost(), mask1);
                    vx.restore(m2);
                }
                if (!prune(vx.cost(), jj + 1)) dfs_chain(tip, jj + 1, prod1, lam1, mask1);
            }
            vx.restore(m);
            if (cap_hit) return;
        }
    }
    void run(const Rec& r, u16 mask) {
        prep(r, mask, P);
        walk_path(G->T, P);
        old_mask = mask;
        vx.reset();
        remove_path(vx, P);
        // per-prime table lower bounds after the removal; primes of n with
        // path counts within 1 of each other may interact (shared class,
        // neighbours, class-size bits): slack 128 bits
        double dl[kMaxK];
        u32 cc[kMaxK];
        for (int jj = 0; jj < P.k; ++jj) cc[jj] = vx.vpc(P.idx[P.ord[jj]]);
        for (int jj = 0; jj < P.k; ++jj) {
            const VirtualCx::Mark m = vx.mark();
            const double c0 = vx.cost();
            vx.tbl(P.idx[P.ord[jj]], +1);
            double d = vx.cost() - c0 - 1.0;
            vx.restore(m);
            for (int o = 0; o < P.k; ++o)
                if (o != jj && (cc[o] + 1 >= cc[jj] && cc[o] <= cc[jj] + 1)) { d -= 128.0; break; }
            dl[jj] = std::min(0.0, d) - g_opt.slack;
        }
        tlb[P.k] = 0;
        for (int jj = P.k - 1; jj >= 0; --jj) tlb[jj] = tlb[jj + 1] + dl[jj];
        {
            static const char* e = std::getenv("DSEL_DEBUG_N");
            dbg = e && to_string(r.n) == e;
            if (dbg) { std::printf("  DBG n=%s mask %u k %d; after removal cost %.6f; dl:", e, mask, P.k, vx.cost()); for (int jj = 0; jj < P.k; ++jj) std::printf(" [rank %u pc %u] %.3f", P.rnk[P.ord[jj]], cc[jj], dl[jj]); std::printf("\n"); }
        }
        best = 0; best_mask = kNoMask; old_seen = false; old_total = 0; cap_hit = false; steps0 = steps;
        NCtx root = path_ctx(P, 0);
        dfs_alive(root, 0, 1, 1, 0);
        if (cap_hit) ++capped;
        if (old_seen && std::fabs(old_total) > 1e-6) {
            ++old_bad;
            if (old_bad <= 5) std::printf("  WARNING n=%s mask %u: re-insertion of the current S costs %.9f bits (must be 0)\n", to_string(r.n).c_str(), mask, old_total);
        }
    }
};

// ------------------------------------------------------------------ build
struct Key { u8 b[kMaxK * 3]; };
inline bool key_less(const Key& a, const Key& b) { return std::memcmp(a.b, b.b, sizeof a.b) < 0; }
inline u32 key_field(const Key& k, int i) { return (u32(k.b[3 * i]) << 16) | (u32(k.b[3 * i + 1]) << 8) | k.b[3 * i + 2]; }
inline int key_len(const Key& k) { int n = 0; while (n < kMaxK && key_field(k, n)) ++n; return n; }

struct BuildInfo {
    u64 nodes = 0; u64 level[kMaxD + 2] = {0};
    double total = 0;
    double ent[kNC] = {0, 0, 0, 0, 0}, rawt = 0, rawb = 0, esc = 0;
    i64 by[kNC] = {0, 0, 0, 0, 0}, csb = 0;
    i64 real_by[kNC] = {0, 0, 0, 0, 0};   // real serialized tables (rans.hpp RansContext)
    double qideal[kNC] = {0, 0, 0, 0, 0}; // quantized ideal bits
    u64 syms[kNC] = {0, 0, 0, 0, 0}, nmod[kNC] = {0, 0, 0, 0, 0};
    u64 classes = 0, m = 0, total_syms = 0;
    u64 path_len_hist[kMaxK + 1] = {0};
    u64 at_dmin = 0, at_all = 0, order_changed = 0; double sum_len = 0;
    u64 keychecks = 0;
    u64 fp = 0;
};

u64 g_prev_level[kMaxD + 2] = {0};
int g_build_no = 0;

// Rebuilds tree, models and table from the masks under ranking G->R and keys
// G->K; fresh = the ranks were just recomputed (keys checked against
// Pwt2Keys2). Checks every S (valid, rebuilds n) and the live prime counts.
BuildInfo build(const std::vector<u32>* prev_rank, bool fresh) {
    State& S = *G;
    const double tb0 = now();
    ++g_build_no;
    S.T.clear();
    S.M.clear();
    const int T = std::max(1, g_opt.threads);
    const u64 N = S.N;
    BuildInfo bi;
    std::vector<std::string> runs;
    std::vector<std::vector<u32>> cnt(T, std::vector<u32>(S.U.P.size(), 0));
    struct TS { u64 lh[kMaxK + 1] = {0}; u64 at_dmin = 0, at_all = 0, oc = 0; double sl = 0; };
    std::vector<TS> ts(T);
    {
        std::vector<Key> buf;
        for (u64 c0 = 0; c0 < N; c0 += g_opt.chunk) {
            const u64 c1 = std::min(N, c0 + g_opt.chunk);
            buf.resize(c1 - c0);
            S.src.each(c0, c1, T, [&](int tid, u64 pos, const Rec& r) {
                const u16 mask = S.masks[pos];
                u32 rk[kMaxK], ro[kMaxK]; int L = 0;
                u128 prod = 1, lam = 1;
                for (int j = 0; j < r.k; ++j) {
                    if (!(mask >> j & 1)) continue;
                    const u32 i = S.U.idx(r.p[j]);
                    CHECK(i != kNone, "prime missing");
                    ++cnt[tid][i];
                    rk[L] = S.R.rank[i];
                    ro[L] = prev_rank ? (*prev_rank)[i] : 0;
                    ++L;
                    prod *= r.p[j]; lam = lcm_step(lam, r.p[j] - 1);
                }
                if (L > kMaxK) __builtin_unreachable();
                CHECK(L >= 1 && valid_pl(prod, lam, r.n) && dmin_value(prod, lam) == r.n,
                      "record %llu: stored S does not rebuild n", (unsigned long long)pos);
                if (prev_rank) {
                    u32 a[kMaxK], b[kMaxK];
                    int ia[kMaxK];
                    for (int x = 0; x < L; ++x) ia[x] = x;
                    std::sort(ia, ia + L, [&](int x, int y) { return rk[x] < rk[y]; });
                    for (int x = 0; x < L; ++x) a[x] = ia[x];
                    for (int x = 0; x < L; ++x) ia[x] = x;
                    std::sort(ia, ia + L, [&](int x, int y) { return ro[x] < ro[y]; });
                    for (int x = 0; x < L; ++x) b[x] = ia[x];
                    if (std::memcmp(a, b, sizeof(u32) * L) != 0) ++ts[tid].oc;
                }
                for (int x = 1; x < L; ++x) { const u32 v = rk[x]; int y = x; while (y > 0 && rk[y - 1] > v) { rk[y] = rk[y - 1]; --y; } rk[y] = v; }
                Key& k = buf[pos - c0];
                std::memset(k.b, 0, sizeof k.b);
                for (int x = 0; x < L; ++x) { k.b[3 * x] = u8(rk[x] >> 16); k.b[3 * x + 1] = u8(rk[x] >> 8); k.b[3 * x + 2] = u8(rk[x]); }
                ++ts[tid].lh[L]; ts[tid].sl += L;
                if (mask == (1u << r.k) - 1) ++ts[tid].at_all;
                if (S.dmask[pos] == mask) ++ts[tid].at_dmin;
            });
            const u64 n = c1 - c0, per = (n + T - 1) / T;
            std::vector<std::thread> th;
            for (int t = 0; t < T; ++t) {
                const u64 a = t * per, b = std::min(n, a + per);
                if (a >= b) break;
                th.emplace_back([&, a, b] { std::sort(buf.begin() + a, buf.begin() + b, key_less); });
            }
            for (auto& x : th) x.join();
            for (int t = 0; t < T; ++t) {
                const u64 a = t * per, b = std::min(n, a + per);
                if (a >= b) break;
                const std::string rp = g_opt.scratch + "/dsel_set2_run_" + std::to_string(proc_id()) + "_" + std::to_string(runs.size()) + ".bin";
                std::ofstream o(rp, std::ios::binary | std::ios::trunc);
                CHECK(o, "cannot create %s", rp.c_str());
                for (u64 w = a; w < b;) {
                    const u64 piece = std::min<u64>(b - w, (u64(256) << 20) / sizeof(Key));
                    o.write(reinterpret_cast<const char*>(buf.data() + w), std::streamsize(piece * sizeof(Key)));
                    w += piece;
                }
                o.close();
                CHECK(o, "write failure %s", rp.c_str());
                runs.push_back(rp);
            }
        }
    }
    {
        u64 bad = 0;
        for (std::size_t i = 0; i < S.U.P.size(); ++i) {
            u64 c = 0;
            for (int t = 0; t < T; ++t) c += cnt[t][i];
            if (c != S.TL.pc[i]) ++bad;
        }
        CHECK(bad == 0, "%llu prime counts differ from the recount", (unsigned long long)bad);
        std::vector<std::vector<u32>>().swap(cnt);
    }
    for (const TS& x : ts) {
        for (int i = 0; i <= kMaxK; ++i) bi.path_len_hist[i] += x.lh[i];
        bi.at_dmin += x.at_dmin; bi.at_all += x.at_all; bi.order_changed += x.oc; bi.sum_len += x.sl;
    }
    const double t_runs = now() - tb0;
    // --- merge + scan into the CSR snapshot; every set-2 symbol of the tree
    struct RR {
        std::ifstream in; std::vector<Key> b; std::size_t p = 0, n = 0; bool eof = false;
        bool fill() { if (eof) return false; b.resize(1 << 16); in.read(reinterpret_cast<char*>(b.data()), std::streamsize(b.size() * sizeof(Key))); n = std::size_t(in.gcount()) / sizeof(Key); p = 0; if (!n) { eof = true; return false; } return true; }
        bool ready() { return p < n || fill(); }
    };
    std::vector<RR> rr(runs.size());
    std::vector<std::size_t> heap;
    for (std::size_t i = 0; i < runs.size(); ++i) {
        rr[i].in.open(runs[i], std::ios::binary);
        CHECK(rr[i].in, "cannot reopen %s", runs[i].c_str());
        if (rr[i].ready()) heap.push_back(i);
    }
    auto hg = [&](std::size_t a, std::size_t b) { return key_less(rr[b].b[rr[b].p], rr[a].b[rr[a].p]); };
    std::make_heap(heap.begin(), heap.end(), hg);
    Tree& Tr = S.T;
    Store& M = S.M;
    const Keys& K = S.K;
    const Pwt2Keys2 KY(S.R.m);
    const u32 mfr = S.R.m;
    for (int t = 0; t <= kMaxD + 1; ++t) {
        const u64 rs = g_prev_level[t] ? g_prev_level[t] + g_prev_level[t] / 16 + 16 : 0;
        if (rs) { Tr.rk[t].reserve(rs); Tr.ct[t].reserve(rs); Tr.nc[t].reserve(rs); Tr.fc[t].reserve(rs + 1); }
    }
    u32 stk[kMaxD + 2];
    u128 Pst[kMaxD + 2];
    Pst[0] = 1;
    double rawt = 0, esc = 0;
    auto bcount = [&](int d, u64 key, u64 c) {
        M.badd(C_CNT, key, static_cast<u32>(std::min<u64>(c, 63)), 1);
        if (c >= 63) M.badd(C_TAIL, static_cast<u64>(d), static_cast<u32>(std::min<u64>(c, 4095)), 1);
        if (c >= 4095) { const u64 v = c - 4094; esc += 2 * bits_of(v) - 1; }
    };
    auto open = [&](int d, u32 r) {
        const u32 i = static_cast<u32>(Tr.rk[d].size());
        Tr.rk[d].push_back(r); Tr.ct[d].push_back(0); Tr.nc[d].push_back(0);
        Tr.fc[d].push_back(static_cast<u32>(Tr.rk[d + 1].size()));
        stk[d] = i;
        if (d > 0) Pst[d] = Pst[d - 1] * S.R.pv[r];
    };
    auto close = [&](int d) {
        const u32 i = stk[d];
        const u32 a = Tr.fc[d][i];
        const u32 c = static_cast<u32>(Tr.rk[d + 1].size()) - a;
        Tr.nc[d][i] = c;
        const u32 rp = Tr.rk[d][i] & kRankMask;
        const u128 Pv = Pst[d];
        if (d == 0) {
            bcount(0, K.count(0, 0, 1, false, 0), c);
            if (fresh) { CHECK(K.count(0, 0, 1, false, 0) == KY.count(0, 0, 1, false, 0), "root key differs from Pwt2Keys2"); ++bi.keychecks; }
        }
        u32 prev = rp;
        for (u32 j = 0; j < c; ++j) {
            const u32 xr = Tr.rk[d + 1][a + j];
            const u32 x = xr & kRankMask;
            CHECK(x > prev, "children not ascending");
            const int L = bits_of(x - prev);
            const u64 gk = K.gap(d + 1, c, j == 0, prev, c - j, Pv);
            M.badd(C_GAP, gk, static_cast<u32>(L), 1);
            rawt += L - 1;
            const u32 cx = Tr.nc[d + 1][a + j];
            const u128 Px = Pv * S.R.pv[x];
            const u64 kk = K.count(d + 1, x, Px, j == 0, c - j - 1);
            bcount(d + 1, kk, cx);
            if (cx > 0) M.badd(C_FLAG, K.flag(d + 1, cx), (xr & kTermBit) ? 1 : 0, 1);
            if (fresh) {
                CHECK(x <= mfr && prev <= mfr, "fresh build with a rank above m");
                CHECK(gk == KY.gap(d + 1, c, j == 0, prev, c - j, Pv), "gap key differs from Pwt2Keys2");
                CHECK(kk == KY.count(d + 1, x, Px, j == 0, c - j - 1), "count key differs from Pwt2Keys2");
                if (cx > 0) CHECK(K.flag(d + 1, cx) == KY.flag(d + 1, cx), "flag key differs from Pwt2Keys2");
                bi.keychecks += 2 + (cx > 0);
            }
            prev = x;
        }
    };
    open(0, 0);
    Key prevk{}; std::memset(prevk.b, 0, sizeof prevk.b);
    int plen = 0;
    u64 scanned = 0;
    while (!heap.empty()) {
        std::pop_heap(heap.begin(), heap.end(), hg);
        const std::size_t ri = heap.back(); heap.pop_back();
        const Key k = rr[ri].b[rr[ri].p++];
        if (rr[ri].ready()) { heap.push_back(ri); std::push_heap(heap.begin(), heap.end(), hg); }
        const int len = key_len(k);
        CHECK(scanned == 0 || key_less(prevk, k), "duplicate or unsorted path (two numbers share an S?)");
        int lcp = 0;
        while (lcp < plen && lcp < len && key_field(prevk, lcp) == key_field(k, lcp)) ++lcp;
        for (int d = plen; d > lcp; --d) close(d);
        for (int d = lcp + 1; d <= len; ++d) open(d, key_field(k, d - 1));
        for (int d = 0; d <= len; ++d) ++Tr.ct[d][stk[d]];
        CHECK(!(Tr.rk[len][stk[len]] & kTermBit), "terminal twice");
        Tr.rk[len][stk[len]] |= kTermBit;
        prevk = k; plen = len; ++scanned;
    }
    for (int d = plen; d >= 0; --d) close(d);
    CHECK(scanned == N, "scanned %llu of %llu", (unsigned long long)scanned, (unsigned long long)N);
    for (auto& x : rr) x.in.close();
    for (const std::string& rp : runs) std::remove(rp.c_str());
    u64 tot = 0;
    for (int d = 0; d <= kMaxD + 1; ++d) {
        Tr.fc[d].push_back(static_cast<u32>(d + 1 <= kMaxD + 1 ? Tr.rk[d + 1].size() : 0));
        Tr.ol[d].assign(Tr.rk[d].size(), 0);
        Tr.off[d] = static_cast<u32>(tot);
        bi.level[d] = Tr.rk[d].size();
        g_prev_level[d] = Tr.rk[d].size();
        tot += Tr.rk[d].size();
        CHECK(tot < 0x7FFFFF00ull, "tree too large for 31-bit refs");
    }
    Tr.off[kMaxD + 2] = static_cast<u32>(tot);
    Tr.snap = tot;
    for (int d = 0; d <= kMaxD + 1; ++d) { Tr.rk[d].shrink_to_fit(); Tr.ct[d].shrink_to_fit(); Tr.nc[d].shrink_to_fit(); Tr.fc[d].shrink_to_fit(); }
    {
        const u64 mo = std::min<u64>(g_opt.max_overflow, G->N);
        const u64 ro = std::min<u64>(mo + mo / 8 + 1024, u64(1) << 31);
        Tr.ork.reserve(ro); Tr.oct.reserve(ro); Tr.onc.reserve(ro); Tr.ool.reserve(ro);
        Tr.pool.reserve(std::min<u64>(5 * ro + (u64(1) << 20), u64(0xFFFFFF00)));
    }
    bi.nodes = tot;
    const double t_tree = now() - tb0 - t_runs;
    // --- the rank table (live counts; classes ascending by position)
    TabLive& TL = S.TL;
    TL.CB = g_opt.cb;
    TL.W = (S.U.P.size() + 63) / 64 + 1;
    TL.bits.assign(static_cast<u64>(TL.CB) * TL.W, 0);
    TL.sv.clear(); TL.sub.clear(); TL.sp.clear();
    TL.dn.assign(kDenseC, 0);
    double rawb = 0; i64 csb = 0;
    {
        std::unordered_map<u32, std::vector<u32>> cls;
        for (u32 i = 0; i < S.U.P.size(); ++i) if (TL.pc[i]) cls[TL.pc[i]].push_back(i);
        bi.classes = cls.size();
        for (auto& kv : cls) {
            const u32 c = kv.first;
            const std::vector<u32>& mem = kv.second;
            const u32 s = static_cast<u32>(mem.size());
            TL.set_csize(c, s);
            csb += vbl(s);
            bi.m += s;
            i64 q = -1;
            for (u32 p : mem) {
                const u64 pp = q >= 0 ? S.U.pos[q] : 0;
                const bool H = q >= 0 && S.U.hi[q];
                const u64 g = q >= 0 ? S.U.pos[p] - pp : S.U.pos[p] + 1;
                const int L = bits_of(g);
                CHECK(L <= 56, "table gap too long");
                const u64 pb = capb(pp, 41);
                const u64 key = K.tablek(s, H, pb);
                if (fresh) { CHECK(key == KY.table(s, H, q >= 0 ? pp : 0), "table key differs from Pwt2Keys2"); ++bi.keychecks; }
                M.badd(C_TAB, key, static_cast<u32>(L), 1);
                rawb += L - 1;
                TL.subadd(c, triple(H, pb, L), 1);
                TL.mem_add(c, p);
                q = p;
            }
        }
    }
    M.finish();
    // --- real serialized tables (for the report) and quantized ideal bits
    {
        std::vector<u64> row(4096);
        std::vector<std::uint8_t> ser;
        for (int c = 0; c < kNC; ++c) {
            bi.real_by[c] = vbl(static_cast<u64>(M.nmod[c]));
            bi.nmod[c] = static_cast<u64>(M.nmod[c]);
            bi.syms[c] = M.syms[c];
        }
        u64 prevk2[kNC] = {0, 0, 0, 0, 0};
        std::vector<std::pair<u64, u32>> ord;
        for (u32 i = 0; i < M.md.size(); ++i) if (M.md[i].T) ord.emplace_back(M.md[i].ck, i);
        std::sort(ord.begin(), ord.end());
        for (const auto& e : ord) {
            const Mdl& m = M.md[e.second];
            const u32 A = kAl[m.cls];
            for (u32 s = 0; s < A; ++s) row[s] = M.H[m.off + s];
            RansContext rc;
            rc.build(row.data(), A);
            ser.clear();
            rc.serialize(ser);
            const u64 key = ck_key(m.ck);
            bi.real_by[m.cls] += vbl(key - prevk2[m.cls]) + static_cast<i64>(ser.size());
            prevk2[m.cls] = key;
            for (std::size_t i = 0; i < rc.syms.size(); ++i)
                bi.qideal[m.cls] += double(row[rc.syms[i]]) * (16.0 - std::log2(double(rc.freq[i])));
        }
    }
    for (int c = 0; c < kNC; ++c) { bi.ent[c] = M.ent[c]; bi.by[c] = M.by[c]; bi.total_syms += M.syms[c]; }
    bi.rawt = rawt; bi.rawb = rawb; bi.esc = esc; bi.csb = csb;
    S.p_rawt = rawt; S.p_rawb = rawb; S.p_esc = esc; S.p_csb = csb;
    for (int c = 0; c < kNC; ++c) { S.p_ent[c] = M.ent[c]; S.p_by[c] = M.by[c]; }
    bi.total = S.total_from_parts();
    S.live = bi.total;
    bi.fp = M.fingerprint();
    // --- sibling caches
    S.huge.clear(); S.huge_id.clear();
    for (int d = 0; d <= kMaxD; ++d)
        for (u32 i = 0; i < Tr.rk[d].size(); ++i)
            if (Tr.nc[d][i] > g_opt.huge) register_huge(Tr.off[d] + i, d);
    std::printf("  build %d: %llu runs %.0f s, tree+symbols %.0f s, table+models %.0f s, nodes %llu (%.4f /el), %zu sibling caches, keys checked vs Pwt2Keys2 %llu, mem ws %.2f GB (peak %.2f)\n",
                g_build_no, (unsigned long long)runs.size(), t_runs, t_tree, now() - tb0 - t_runs - t_tree, (unsigned long long)tot,
                double(tot) / double(N), S.huge.size(), (unsigned long long)bi.keychecks, mem().ws_gb, mem().peak_gb);
    return bi;
}

void print_model(const char* tag, const BuildInfo& bi) {
    const double Nd = double(G->N);
    double ent = 0, qid = 0; i64 by = 0, rby = 0;
    for (int c = 0; c < kNC; ++c) { ent += bi.ent[c]; qid += bi.qideal[c]; by += bi.by[c]; rby += bi.real_by[c]; }
    // predicted file: header 160 + classes + tables (real) + symbols (quantized ideal + 8 B per block + 4 B rANS flush) + raw + footer 88
    const u64 blocks = (bi.total_syms + (1u << 20) - 1) >> 20;
    const double file = 160.0 + double(bi.csb) + double(rby) + std::ceil(qid / 8.0) + 12.0 * double(blocks) + std::ceil((bi.rawt + bi.rawb + bi.esc) / 8.0) + 88.0;
    std::printf("MODEL %-18s nodes %llu = %.4f /el | OBJECTIVE %.4f bits/el = %.0f bytes | predicted file %.0f bytes = %.4f bits/el\n",
                tag, (unsigned long long)bi.nodes, double(bi.nodes) / Nd, bi.total / Nd, bi.total / 8.0, file, 8.0 * file / Nd);
    std::printf("      entropy (bits/el): table %.4f + count %.4f + tail %.4f + flag %.4f + gap %.4f = %.4f (quantized ideal %.4f); raw tree %.4f + raw table %.4f + escapes %.4f\n",
                bi.ent[0] / Nd, bi.ent[1] / Nd, bi.ent[2] / Nd, bi.ent[3] / Nd, bi.ent[4] / Nd, ent / Nd, qid / Nd, bi.rawt / Nd, bi.rawb / Nd, bi.esc / Nd);
    std::printf("      tables section: model %lld B (%.4f bits/el) vs REAL serialized %lld B (%+.3f%%); per class model/real/models:",
                (long long)by, 8.0 * by / Nd, (long long)rby, 100.0 * (double(by) - double(rby)) / std::max(1.0, double(rby)));
    for (int c = 0; c < kNC; ++c) std::printf(" %s %lld/%lld/%llu", kCName[c], (long long)bi.by[c], (long long)bi.real_by[c], (unsigned long long)bi.nmod[c]);
    std::printf("\n      rank table: m %llu primes in %llu classes, class sizes %lld B, symbols %.4f + raw %.4f bits/el; symbols %llu (blocks %llu); fingerprint %016llx\n",
                (unsigned long long)bi.m, (unsigned long long)bi.classes, (long long)bi.csb, bi.ent[0] / Nd, bi.rawb / Nd,
                (unsigned long long)bi.total_syms, (unsigned long long)blocks, (unsigned long long)bi.fp);
    std::printf("      nodes by depth:");
    for (int d = 0; d <= kMaxD; ++d) if (bi.level[d]) std::printf(" %d:%llu", d, (unsigned long long)bi.level[d]);
    std::printf("\n      paths: mean |S| %.4f, S=d_min %.3f%%, S=all primes %.3f%%, |S| histogram:", bi.sum_len / Nd,
                100.0 * bi.at_dmin / Nd, 100.0 * bi.at_all / Nd);
    for (int i = 0; i <= kMaxK; ++i) if (bi.path_len_hist[i]) std::printf(" %d:%llu", i, (unsigned long long)bi.path_len_hist[i]);
    std::printf("\n");
}

// --------------------------------------------------------------- commit
// Applies an evaluated move (the overlay of vx on the live state) to the
// live model store, rank table, tree and caches; new chain nodes materialized.
void commit(VirtualCx& vx, const Prep& Pn) {
    State& S = *G;
    // models
    for (const VM& v : vx.vms) {
        const bool lex = v.mid >= 0 && S.M.md[v.mid].T > 0;
        u32 id;
        if (v.mid >= 0) id = static_cast<u32>(v.mid);
        else { if (v.T == 0) continue; id = S.M.create(v.ck); }
        Mdl& L = S.M.md[id];
        std::memcpy(S.M.H.data() + L.off, vx.vcnt.data() + v.hoff, sizeof(u32) * kAl[v.cls]);
        std::memcpy(S.M.MK.data() + L.moff, vx.vmsk.data() + v.moff, sizeof(u64) * kAW[v.cls]);
        L.T = v.T; L.nsym = v.nsym; L.fb = v.fb; L.sd = v.sd;
        std::vector<u64>& K = S.M.keys[v.cls];
        const u64 key = ck_key(v.ck);
        if (!lex && v.T > 0) K.insert(std::lower_bound(K.begin(), K.end(), key), key);
        else if (lex && v.T == 0) { auto it = std::lower_bound(K.begin(), K.end(), key); CHECK(it != K.end() && *it == key, "key missing"); K.erase(it); }
    }
    for (int c = 0; c < kNC; ++c) {
        S.M.nmod[c] += vx.sc.nm[c];
        CHECK(S.M.nmod[c] == static_cast<i64>(S.M.keys[c].size()), "model count mismatch in class %d", c);
        S.p_ent[c] += vx.sc.ent[c]; S.p_by[c] += vx.sc.by[c];
    }
    S.p_rawt += vx.sc.rawt; S.p_rawb += vx.sc.rawb; S.p_esc += vx.sc.esc; S.p_csb += vx.sc.csb;
    S.live += vx.sc.cost;
    // rank table
    TabLive& TL = S.TL;
    {
        SMap& m = vx.mp[VirtualCx::MTP];
        for (u32 sl : m.touched) {
            Slot& s = m.s[sl];
            if (s.key == 0 || s.a == 0) continue;
            const u32 i = static_cast<u32>(s.key - 1);
            const u32 a = TL.pc[i], b = static_cast<u32>(i64(a) + s.a);
            if (a > 0) TL.mem_del(a, i);
            if (b > 0) TL.mem_add(b, i);
            TL.pc[i] = b;
            s.a = 0;
        }
        SMap& mc = vx.mp[VirtualCx::MCS];
        for (u32 sl : mc.touched) {
            Slot& s = mc.s[sl];
            if (s.key == 0 || s.a == 0) continue;
            const u32 c = static_cast<u32>(s.key - 1);
            TL.set_csize(c, static_cast<u32>(i64(TL.csize(c)) + s.a));
            s.a = 0;
        }
        SMap& ms = vx.mp[VirtualCx::MSUB];
        for (u32 sl : ms.touched) {
            Slot& s = ms.s[sl];
            if (s.key == 0 || s.a == 0) continue;
            const u64 k = s.key - 1;
            TL.subadd(static_cast<u32>(k >> 16), static_cast<u32>(k & 0xFFFF), s.a);
            s.a = 0;
        }
    }
    // tree: counts, child counts, terminal bits of existing nodes
    Tree& T = S.T;
    std::vector<std::pair<u32, int>> dirty;
    {
        SMap& mn = vx.mp[VirtualCx::MNODE];
        for (u32 sl : mn.touched) {
            Slot& s = mn.s[sl];
            if (s.key == 0 || (s.a == 0 && s.b == 0)) continue;
            const u32 ref = static_cast<u32>(s.key - 1);
            const int t = (ref & kOV) ? 0 : T.depth_of(ref);
            CHECK(t >= 0, "bad node ref");
            T.CT(ref, t) = static_cast<u32>(i64(T.CT(ref, t)) + s.a);
            T.NC(ref, t) = static_cast<u32>(i64(T.NC(ref, t)) + s.b);
            s.a = 0; s.b = 0;
        }
        SMap& mt = vx.mp[VirtualCx::MTERM];
        for (u32 sl : mt.touched) {
            Slot& s = mt.s[sl];
            if (s.key == 0) continue;
            const u32 ref = static_cast<u32>(s.key - 1);
            const int t = (ref & kOV) ? 0 : T.depth_of(ref);
            u32& x = T.RK(ref, t);
            x = s.a ? (x | kTermBit) : (x & ~kTermBit);
            s.key = 0;   // processed (duplicates in touched)
        }
        SMap& mv = vx.mp[VirtualCx::MVIEW];
        for (u32 sl : mv.touched) {
            Slot& s = mv.s[sl];
            if (s.key == 0) continue;
            dirty.emplace_back(static_cast<u32>(s.key - 1), static_cast<int>(s.a >> 40));
            s.key = 0;
        }
    }
    // materialize the new path's chain (every chain parent's child list changes)
    {
        u32 cur = 0;
        for (int t = 1; t <= Pn.L; ++t) {
            const u32 r = Pn.prk[t];
            const u32 ex = T.child_ref(cur, t - 1, r);
            if (ex != kNone && T.CT(ex, t) > 0) { cur = ex; continue; }
            dirty.emplace_back(cur, t - 1);
            u32 nr;
            if (ex != kNone) {
                CHECK(T.NC(ex, t) == 0, "revive of a node with children");
                nr = ex;
            } else {
                const u32 j = static_cast<u32>(T.ork.size());
                CHECK(j < 0x7FFFFF00u, "overflow ids exhausted");
                T.ork.push_back(r); T.oct.push_back(0); T.onc.push_back(0); T.ool.push_back(0);
                T.ol_insert(cur, t - 1, r, j);
                nr = kOV | j;
            }
            T.CT(nr, t) = 1;
            T.NC(nr, t) = t < Pn.L ? 1 : 0;
            u32& x = T.RK(nr, t);
            x = t == Pn.L ? (x | kTermBit) : (x & ~kTermBit);
            cur = nr;
        }
    }
    // sibling caches of every node whose child list changed
    std::sort(dirty.begin(), dirty.end());
    dirty.erase(std::unique(dirty.begin(), dirty.end()), dirty.end());
    for (const auto& d : dirty) {
        auto it = S.huge_id.find(d.first);
        if (it != S.huge_id.end()) build_cache(S.huge[it->second]);
        else if (T.NC(d.first, d.second) > g_opt.huge) register_huge(d.first, d.second);
    }
}

// --------------------------------------------------------------- apply
struct ApplyStats { u64 moves = 0, rejected = 0, researched = 0, eval_mismatch = 0; double gain = 0; };

void apply_move(const Rec& r, u64 pos, u32 new_mask, double eval_total, bool immediate, VirtualCx& vx, ApplyStats& as) {
    State& S = *G;
    const u16 old_mask = S.masks[pos];
    Prep Po, Pn;
    prep(r, old_mask, Po);
    walk_path(S.T, Po);
    prep(r, static_cast<u16>(new_mask), Pn);
    {
        u128 prod = 1, lam = 1;
        for (int j = 0; j < r.k; ++j) if (new_mask >> j & 1) { prod *= r.p[j]; lam = lcm_step(lam, r.p[j] - 1); }
        CHECK(valid_pl(prod, lam, r.n), "proposed S invalid");
    }
    vx.reset();
    remove_path(vx, Po);
    insert_virtual(vx, Pn);
    const double dv = vx.cost();
    if (immediate && std::fabs(dv - eval_total) > 1e-6) {
        if (++as.eval_mismatch <= 10)
            std::printf("  WARNING record %llu: search total %.9f vs re-cost %.9f\n", (unsigned long long)pos, eval_total, dv);
    }
    if (!(dv < -1e-9)) { ++as.rejected; return; }
    commit(vx, Pn);
    S.masks[pos] = static_cast<u16>(new_mask);
    ++as.moves; as.gain += dv;
}

// ------------------------------------------------------------------ search
struct SweepStats { ApplyStats a; u64 steps = 0, capped = 0, old_bad = 0, rebuilds = 0, evaluated = 0; double secs = 0, t_eval = 0, t_apply = 0, t_rebuild = 0; };

std::ofstream g_dump;

void search(int sweep, SweepStats& ss) {
    State& S = *G;
    const double ts = now();
    const int T = std::max(1, g_opt.threads);
    std::vector<std::unique_ptr<Evaluator>> ev;
    for (int t = 0; t < T; ++t) ev.emplace_back(new Evaluator());
    std::unique_ptr<VirtualCx> avx(new VirtualCx(S));
    const u64 N = S.N;
    u64 next_report = N / 20 + 1;
    std::vector<Rec> batch;
    std::vector<u64> evn(T, 0);
    struct Prop { u32 mask; double total; };
    std::vector<Prop> prop;
    const double obj0 = S.live;
    std::vector<std::pair<u32, double>> dc;
    for (u64 r0 = 0; r0 < N;) {
        const u64 B = T == 1 ? std::max<u64>(g_opt.batch, 4096) : g_opt.batch;
        const u64 r1 = std::min(N, r0 + B);
        if (T == 1) {
            S.src.each(r0, r1, 1, [&](int, u64 pos, const Rec& r) {
                if (pos % g_opt.stride) return;
                ++evn[0];
                if (g_dump) { dc.clear(); ev[0]->dumpc = &dc; }
                const u16 om = S.masks[pos];
                ev[0]->run(r, om);
                ev[0]->dumpc = nullptr;
                if (ev[0]->best_mask != kNoMask) apply_move(r, pos, ev[0]->best_mask, ev[0]->best, true, *avx, ss.a);
                if (g_dump) {
                    g_dump << "REC " << pos << " " << om << " " << S.masks[pos] << "\n";
                    char buf[96];
                    for (const auto& c : dc) { std::snprintf(buf, sizeof buf, "C %u %.9f\n", c.first, c.second); g_dump << buf; }
                }
            });
        } else {
            batch.resize(r1 - r0); prop.resize(r1 - r0);
            const double te0 = now();
            S.src.each(r0, r1, T, [&](int tid, u64 pos, const Rec& r) {
                batch[pos - r0] = r;
                if (pos % g_opt.stride) { prop[pos - r0] = {kNoMask, 0.0}; return; }
                ++evn[tid];
                ev[tid]->run(r, S.masks[pos]);
                prop[pos - r0] = {ev[tid]->best_mask, ev[tid]->best};
            });
            const double te1 = now();
            ss.t_eval += te1 - te0;
            for (u64 i = 0; i < r1 - r0; ++i)
                if (prop[i].mask != kNoMask) {
                    const u64 rej = ss.a.rejected;
                    apply_move(batch[i], r0 + i, prop[i].mask, prop[i].total, false, *avx, ss.a);
                    if (ss.a.rejected != rej) {
                        ++ss.a.researched;
                        ev[0]->run(batch[i], S.masks[r0 + i]);
                        if (ev[0]->best_mask != kNoMask) apply_move(batch[i], r0 + i, ev[0]->best_mask, ev[0]->best, true, *avx, ss.a);
                    }
                }
            ss.t_apply += now() - te1;
        }
        r0 = r1;
        if (r0 >= next_report || r0 == N) {
            next_report += N / 20 + 1;
            u64 st = 0, en = 0; for (auto& e : ev) st += e->steps; for (u64 x : evn) en += x;
            const double sec = now() - ts;
            const Mem mm = mem();
            std::printf("  sweep %d: %6.2f%% records, moves %llu (%.3f%%), rejected %llu, live %.4f bits/el (%+.3f%%), overflow nodes %llu, caches %zu, steps/eval %.1f, %.0f s (ETA %.0f s), ws %.2f GB\n",
                        sweep, 100.0 * r0 / N, (unsigned long long)ss.a.moves, 100.0 * ss.a.moves / r0, (unsigned long long)ss.a.rejected,
                        S.live / N, 100.0 * (S.live - obj0) / obj0, (unsigned long long)S.T.overflow(), S.huge.size(),
                        double(st) / double(std::max<u64>(1, en)), sec, sec / r0 * (N - r0), mm.ws_gb);
        }
        if (r0 < N && S.T.overflow() > g_opt.max_overflow) {
            const double lt = S.live;
            std::printf("  overflow %llu nodes > %llu: rebuilding with the same ranks\n", (unsigned long long)S.T.overflow(), (unsigned long long)g_opt.max_overflow);
            for (auto& e : ev) e.reset();
            avx.reset();
            const double tr0 = now();
            BuildInfo bi = build(nullptr, false);
            ss.t_rebuild += now() - tr0;
            for (int t = 0; t < T; ++t) ev[t].reset(new Evaluator());
            avx.reset(new VirtualCx(S));
            std::printf("  CHECK live %.6f vs rebuilt %.6f bits (diff %.3g)\n", lt, bi.total, bi.total - lt);
            ++ss.rebuilds;
        }
    }
    for (auto& e : ev) { ss.steps += e->steps; ss.capped += e->capped; ss.old_bad += e->old_bad; }
    for (u64 x : evn) ss.evaluated += x;
    ss.secs = now() - ts;
}

// ------------------------------------------------------------------ emits
void emit_paths(const std::string& fn, bool dmin) {
    State& S = *G;
    std::ofstream o(fn, std::ios::binary | std::ios::trunc);
    CHECK(o, "cannot create %s", fn.c_str());
    const u64 CH = u64(1) << 22;
    std::vector<u8> buf;
    for (u64 c0 = 0; c0 < S.N; c0 += CH) {
        const u64 c1 = std::min(S.N, c0 + CH);
        buf.assign((c1 - c0) * 70, 0);
        S.src.each(c0, c1, std::max(1, g_opt.threads), [&](int, u64 pos, const Rec& r) {
            u16 mask = dmin ? S.dmask[pos] : S.masks[pos];
            u128 prod = 1, lam = 1;
            u8* rec = buf.data() + (pos - c0) * 70;
            int len = 0;
            for (int j = 0; j < r.k; ++j) {
                if (!(mask >> j & 1)) continue;
                prod *= r.p[j]; lam = lcm_step(lam, r.p[j] - 1);
                for (int q = 0; q < 5; ++q) rec[len * 5 + q] = u8(r.p[j] >> (8 * (4 - q)));
                ++len;
            }
            CHECK(len >= 1 && valid_pl(prod, lam, r.n) && dmin_value(prod, lam) == r.n, "emit: S does not rebuild n at %llu", (unsigned long long)pos);
        });
        o.write(reinterpret_cast<const char*>(buf.data()), std::streamsize(buf.size()));
    }
    o.close();
    CHECK(o, "write failure %s", fn.c_str());
    std::printf("emitted    %s (%llu records x 14 fields, every S rebuilds its n)\n", fn.c_str(), (unsigned long long)S.N);
}
void save_masks(const std::string& fn) {
    std::ofstream o(fn, std::ios::binary | std::ios::trunc);
    CHECK(o, "cannot create %s", fn.c_str());
    const u64 n = G->N;
    o.write(reinterpret_cast<const char*>(&n), 8);
    for (u64 w = 0; w < n;) {
        const u64 piece = std::min<u64>(n - w, u64(1) << 27);
        o.write(reinterpret_cast<const char*>(G->masks.data() + w), std::streamsize(piece * 2));
        w += piece;
    }
    o.close();
    CHECK(o, "write failure %s", fn.c_str());
}

void rank_movement(const Universe& U, const Ranking& A, const Ranking& B) {
    const u32 n = static_cast<u32>(U.P.size());
    u64 changed = 0, top_changed = 0, newly = 0, gone = 0, used_both = 0;
    for (u32 i = 0; i < n; ++i) {
        const bool ua = A.rank[i] <= A.m, ub = B.rank[i] <= B.m;
        if (!ua && ub) ++newly;
        if (ua && !ub) ++gone;
        if (ua && ub) { ++used_both; if (A.rank[i] != B.rank[i]) ++changed; }
    }
    for (u32 r = 1; r <= std::min<u32>(1000, A.m); ++r) if (B.r2i[r] != A.r2i[r]) ++top_changed;
    std::printf("RANKS  used %u -> %u (newly used %llu, dropped %llu); of %llu used in both, rank changed %llu; top-1000 positions changed %llu; top 10:",
                A.m, B.m, (unsigned long long)newly, (unsigned long long)gone, (unsigned long long)used_both,
                (unsigned long long)changed, (unsigned long long)top_changed);
    for (u32 r = 1; r <= 10 && r <= B.m; ++r) std::printf(" %llu", (unsigned long long)U.P[B.r2i[r]]);
    std::printf("\n");
}

void dump_state(const char* tag, int sweep, double total, u64 nodes) {
    if (!g_dump) return;
    char buf[256];
    std::snprintf(buf, sizeof buf, "SWEEP %d %s\nMODEL %.9f %llu\n", sweep, tag, total, (unsigned long long)nodes);
    g_dump << buf << "MASKS";
    for (u64 i = 0; i < G->N; ++i) g_dump << " " << G->masks[i];
    g_dump << "\n";
    g_dump.flush();
}

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    g_t0 = now();
    Opt& o = g_opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() -> std::string { CHECK(i + 1 < argc, "missing value for %s", a.c_str()); return argv[++i]; };
        if (a == "--orc1") o.orc = nx();
        else if (a == "--text") { while (i + 1 < argc && argv[i + 1][0] != '-') o.text.push_back(argv[++i]); }
        else if (a == "--max-bits") o.max_bits = std::atoi(nx().c_str());
        else if (a == "--sweeps") o.sweeps = std::atoi(nx().c_str());
        else if (a == "--threads") o.threads = std::atoi(nx().c_str());
        else if (a == "--batch") o.batch = std::strtoull(nx().c_str(), nullptr, 10);
        else if (a == "--chunk") o.chunk = std::strtoull(nx().c_str(), nullptr, 10);
        else if (a == "--scratch") o.scratch = nx();
        else if (a == "--max-overflow") o.max_overflow = std::strtoull(nx().c_str(), nullptr, 10);
        else if (a == "--check-rebuild") o.check_rebuild = true;
        else if (a == "--masks-in") o.masks_in = nx();
        else if (a == "--masks-out") o.masks_out = nx();
        else if (a == "--emit") o.emit = nx();
        else if (a == "--emit-dmin") o.emit_dmin = nx();
        else if (a == "--dump") o.dump = nx();
        else if (a == "--no-prune") o.no_prune = true;
        else if (a == "--table-weight") o.tw = std::atof(nx().c_str());
        else if (a == "--min-gain") o.min_gain = std::atof(nx().c_str());
        else if (a == "--step-cap") o.step_cap = std::strtoull(nx().c_str(), nullptr, 10);
        else if (a == "--sample-stride") o.stride = std::max<u64>(1, std::strtoull(nx().c_str(), nullptr, 10));
        else if (a == "--dpaths") o.dpaths = nx();
        else if (a == "--dmin-check") o.dmin_check = std::max<u64>(1, std::strtoull(nx().c_str(), nullptr, 10));
        else if (a == "--huge") o.huge = static_cast<u32>(std::strtoul(nx().c_str(), nullptr, 10));
        else if (a == "--cb") o.cb = static_cast<u32>(std::strtoul(nx().c_str(), nullptr, 10));
        else if (a == "--slack") o.slack = std::atof(nx().c_str());
        else if (a == "--sieve-limit") o.sieve_limit = std::strtoull(nx().c_str(), nullptr, 10);
        else die("unexpected arg %s", a.c_str());
    }
    CHECK(!o.orc.empty() || !o.text.empty(), "usage: dsel_set2 (--orc1 F [--max-bits B] | --text FIX...) [options]; see the head comment");
    CHECK(o.threads >= 1 && o.threads <= 6, "--threads must be 1..6");
    CHECK(o.batch >= 1 && o.chunk >= 1024, "bad --batch / --chunk");
    CHECK(o.cb >= 1 && o.huge >= 1, "bad --cb / --huge");
    init_df();
    State S;
    G = &S;
    S.src.open(o);
    S.N = S.src.N;
    CHECK(S.N > 0, "no records");
    std::printf("dsel_set2  %s%s: %llu records; threads %d, batch %llu, max overflow %llu, huge > %u children, bitset classes <= %u, slack %.1f, table weight %.4f%s, scratch %s\n",
                o.text.empty() ? o.orc.c_str() : "fixture text", o.max_bits ? (" below 2^" + std::to_string(o.max_bits)).c_str() : "",
                (unsigned long long)S.N, o.threads, (unsigned long long)o.batch, (unsigned long long)o.max_overflow, o.huge, o.cb, o.slack, o.tw,
                o.no_prune ? ", NO PRUNING" : "", o.scratch.c_str());
    {
        const int T = std::max(1, o.threads);
        std::vector<u8> seen(kSmall >> 1, 0);
        std::vector<std::vector<u64>> big(T);
        S.src.each(0, S.N, T, [&](int tid, u64, const Rec& r) {
            for (int j = 0; j < r.k; ++j) {
                if (r.p[j] < kSmall) seen[r.p[j] >> 1] = 1;
                else {
                    big[tid].push_back(r.p[j]);
                    if (big[tid].size() >= (u64(1) << 25)) { std::sort(big[tid].begin(), big[tid].end()); big[tid].erase(std::unique(big[tid].begin(), big[tid].end()), big[tid].end()); if (big[tid].capacity() > 2 * big[tid].size() + (u64(1) << 25)) big[tid].shrink_to_fit(); }
                }
            }
        });
        for (u64 i = 1; i < (kSmall >> 1); ++i) if (seen[i]) S.U.P.push_back(2 * i + 1);
        std::vector<u64> all;
        for (auto& b : big) { all.insert(all.end(), b.begin(), b.end()); std::vector<u64>().swap(b); }
        std::sort(all.begin(), all.end());
        all.erase(std::unique(all.begin(), all.end()), all.end());
        S.U.P.insert(S.U.P.end(), all.begin(), all.end());
        S.U.finish(o.sieve_limit);
        std::printf("universe   %zu distinct primes (%.0f s)\n", S.U.P.size(), el());
    }
    S.masks.assign(S.N, 0);
    S.dmask.assign(S.N, 0);
    S.TL.pc.assign(S.U.P.size(), 0);
    if (!o.masks_in.empty()) {
        std::ifstream in(o.masks_in, std::ios::binary);
        CHECK(in, "cannot open %s", o.masks_in.c_str());
        u64 n = 0;
        in.read(reinterpret_cast<char*>(&n), 8);
        CHECK(n == S.N, "masks file has %llu records, table %llu", (unsigned long long)n, (unsigned long long)S.N);
        for (u64 w = 0; w < n;) {
            const u64 piece = std::min<u64>(n - w, u64(1) << 27);
            in.read(reinterpret_cast<char*>(S.masks.data() + w), std::streamsize(piece * 2));
            w += piece;
        }
        CHECK(in, "short masks file");
    }
    {
        const int T = std::max(1, o.threads);
        std::vector<std::vector<u32>> cnt(T, std::vector<u32>(S.U.P.size(), 0));
        const bool from_file = !o.masks_in.empty();
        std::vector<std::ifstream> dp;
        if (!o.dpaths.empty()) {
            for (int t = 0; t < T; ++t) { dp.emplace_back(o.dpaths, std::ios::binary); CHECK(dp.back(), "cannot open %s", o.dpaths.c_str()); }
            dp[0].seekg(0, std::ios::end);
            CHECK(u64(dp[0].tellg()) >= S.N * 60, "%s shorter than %llu records", o.dpaths.c_str(), (unsigned long long)S.N);
        }
        std::vector<u64> dp_next(T, ~u64(0));
        std::atomic<u64> dchecked{0};
        S.src.each(0, S.N, T, [&](int tid, u64 pos, const Rec& r) {
            {
                u16 dm = 0;
                if (!dp.empty()) {
                    if (dp_next[tid] != pos) dp[tid].seekg(std::streamoff(pos * 60));
                    u8 rec[60];
                    dp[tid].read(reinterpret_cast<char*>(rec), 60);
                    CHECK(dp[tid].gcount() == 60, "short read in %s", o.dpaths.c_str());
                    dp_next[tid] = pos + 1;
                    u128 prod = 1, lam = 1;
                    int j = 0;
                    for (int f = 0; f < 12; ++f) {
                        u64 v = 0;
                        for (int q = 0; q < 5; ++q) v = (v << 8) | rec[5 * f + q];
                        if (!v) break;
                        while (j < r.k && r.p[j] < v) ++j;
                        CHECK(j < r.k && r.p[j] == v, "dpaths record %llu: prime %llu does not divide n", (unsigned long long)pos, (unsigned long long)v);
                        dm |= u16(1u << j);
                        prod *= v; lam = lcm_step(lam, v - 1);
                    }
                    CHECK(dm && valid_pl(prod, lam, r.n), "dpaths record %llu: invalid divisor", (unsigned long long)pos);
                    if (pos % o.dmin_check == 0) {
                        const DminChoice dc = choose_dmin(r.p, r.k, r.n);
                        u16 d2 = 0;
                        for (int jj = 0; jj < r.k; ++jj) if (dc.d % r.p[jj] == 0) d2 |= u16(1u << jj);
                        CHECK(d2 == dm, "dpaths record %llu differs from choose_dmin", (unsigned long long)pos);
                        ++dchecked;
                    }
                } else {
                    const DminChoice dc = choose_dmin(r.p, r.k, r.n);
                    CHECK(dmin_value(dc.d, dc.lambda) == r.n, "d_min failure at %llu", (unsigned long long)pos);
                    for (int j = 0; j < r.k; ++j) if (dc.d % r.p[j] == 0) dm |= u16(1u << j);
                }
                S.dmask[pos] = dm;
                if (!from_file) S.masks[pos] = dm;
            }
            const u16 mask = S.masks[pos];
            CHECK(mask != 0 && mask < (1u << r.k), "bad mask at %llu", (unsigned long long)pos);
            for (int j = 0; j < r.k; ++j) if (mask >> j & 1) ++cnt[tid][S.U.idx(r.p[j])];
        });
        for (std::size_t i = 0; i < S.U.P.size(); ++i) { u64 c = 0; for (int t = 0; t < T; ++t) c += cnt[t][i]; S.TL.pc[i] = static_cast<u32>(c); }
        std::printf("start      %s masks (%.0f s)%s\n", from_file ? o.masks_in.c_str() : "d_min", el(),
                    dp.empty() ? "" : (" [d_min from " + o.dpaths + ", " + std::to_string(dchecked.load()) + " records re-derived by choose_dmin, all equal]").c_str());
    }
    if (!o.dump.empty()) { g_dump.open(o.dump, std::ios::trunc); CHECK(g_dump, "cannot create %s", o.dump.c_str()); }
    if (!o.emit_dmin.empty()) emit_paths(o.emit_dmin, true);

    S.R = make_ranking(S.U, S.TL.pc);
    S.K.set(S.R.m);
    BuildInfo bi = build(nullptr, true);
    print_model("start", bi);
    dump_state("start", 0, bi.total, bi.nodes);
    const double start_total = bi.total;
    double best_total = bi.total; int best_sweep = 0;
    std::vector<u16> best_masks = S.masks;
    double prev_total = start_total;
    for (int s = 1; s <= o.sweeps; ++s) {
        SweepStats ss;
        search(s, ss);
        const double Nd = double(S.N);
        const Mem mm = mem();
        std::printf("SWEEP %d  moves %llu (%.3f%%), rejected %llu (re-searched %llu), gain %.4f bits/el, live (old ranks) %.4f bits/el; steps/eval %.1f, evaluated %llu (%.2f us each incl. apply, wall), capped %llu, re-insertion violations %llu, eval/apply mismatches %llu, rebuilds %llu; %.0f s (parallel search %.0f s, sequential apply %.0f s, rebuilds %.0f s); ws %.2f GB, peak %.2f GB, commit %.2f GB\n",
                    s, (unsigned long long)ss.a.moves, 100.0 * ss.a.moves / Nd, (unsigned long long)ss.a.rejected, (unsigned long long)ss.a.researched,
                    -ss.a.gain / Nd, S.live / Nd, double(ss.steps) / double(std::max<u64>(1, ss.evaluated)),
                    (unsigned long long)ss.evaluated, 1e6 * ss.secs / double(std::max<u64>(1, ss.evaluated)),
                    (unsigned long long)ss.capped, (unsigned long long)ss.old_bad, (unsigned long long)ss.a.eval_mismatch,
                    (unsigned long long)ss.rebuilds, ss.secs, ss.t_eval, ss.t_apply, ss.t_rebuild, mm.ws_gb, mm.peak_gb, mm.commit_gb);
        dump_state("search", s, S.live, 0);
        const double live_parts = S.total_from_parts();
        if (o.check_rebuild) {
            const double lt = S.live;
            const u64 lfp = S.M.fingerprint();
            BuildInfo cb = build(nullptr, false);
            std::printf("CHECK  same-rank rebuild: live %.6f vs exact %.6f bits (diff %.3g; live parts %.6f); histogram fingerprint live %016llx vs rebuilt %016llx %s\n",
                        lt, cb.total, cb.total - lt, live_parts, (unsigned long long)lfp, (unsigned long long)cb.fp, lfp == cb.fp ? "EQUAL" : "DIFFER");
            CHECK(std::fabs(cb.total - lt) < 1e-3 * std::max(1.0, lt / 1e6) && lfp == cb.fp, "live model differs from the from-scratch recount");
        }
        const double same_rank = S.live;
        Ranking R2 = make_ranking(S.U, S.TL.pc);
        rank_movement(S.U, S.R, R2);
        std::vector<u32> old_rank = std::move(S.R.rank);
        S.R = std::move(R2);
        S.K.set(S.R.m);
        bi = build(&old_rank, true);
        std::printf("      paths whose prime order changed under the new ranks: %llu (%.3f%%)\n",
                    (unsigned long long)bi.order_changed, 100.0 * bi.order_changed / Nd);
        char tag[64];
        std::snprintf(tag, sizeof tag, "after sweep %d", s);
        print_model(tag, bi);
        const double tot = bi.total;
        std::printf("      re-rank effect %+.4f bits/el; vs start %+.3f%%\n", (tot - same_rank) / Nd, 100.0 * (tot - start_total) / start_total);
        dump_state(tag, s, bi.total, bi.nodes);
        const double g = 100.0 * (prev_total - tot) / prev_total;
        prev_total = tot;
        if (tot < best_total) { best_total = tot; best_sweep = s; best_masks = S.masks; }
        if (!o.masks_out.empty()) { std::vector<u16> cur; cur.swap(S.masks); S.masks = best_masks; save_masks(o.masks_out); S.masks.swap(cur); }
        if (ss.a.moves == 0 || g < o.min_gain) { std::printf("stop       sweep gain %.4f%% < %.4f%% (or no moves)\n", g, o.min_gain); break; }
    }
    std::printf("best       state after sweep %d: %.4f bits/el (model, fresh ranks) = %+.3f%% vs start (emitted / saved)\n", best_sweep,
                best_total / double(S.N), 100.0 * (best_total - start_total) / start_total);
    S.masks = best_masks;
    if (!o.emit.empty()) emit_paths(o.emit, false);
    const Mem mm = mem();
    std::printf("elapsed    %.1f s, peak working set %.2f GB\n", el(), mm.peak_gb);
    return 0;
}
