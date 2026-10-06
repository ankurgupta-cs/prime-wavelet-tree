// paths_cmp -- compare two prime-path files as SETS, whatever their record
// order (e.g. a pwt2_decode --emit-paths export against the encoder's input).
//
// Input format: the pwt2_encode / trie_paths record (pwt_paths.hpp): K
// five-byte big-endian primes per record, strictly ascending, zero padded.
// pwt2_decode --emit-paths writes records in TREE order, trie_paths and
// dsel_set2 --emit in n order, so two files holding the same paths differ
// as byte streams. This tool sorts each file externally (pwt_paths.hpp's
// PathRunSource, ascending record layout, memcmp order of the record bytes;
// a duplicate record fails the run), then walks the two sorted streams in
// lockstep (merge-join) and classifies every record as common, only in A or
// only in B. For each file it also prints the record count, the path-length
// histogram and the SHA-256 of its SORTED record stream: two files hold the
// same set of paths exactly when the merge-join finds no one-sided record,
// and then the two sorted-stream hashes are equal too.
//
// usage: paths_cmp <a.bin> <b.bin> [--fields 12|14] [--scratch DIR]
//                  [--chunk R] [--threads N]
//   --fields   record width K (default 14); both files use the same K
//   --scratch  directory for the sorted run files (default data/scratch);
//              run files are named pathscmp_<pid>_{a,b}_run_<i>.bin and
//              removed at exit. Needs about the two input sizes in space.
//   --chunk    records per sorted run (default 2^25 = 2.35 GB at K = 14)
//   --threads  sort threads (default: hardware threads - 1)
// Last line: SAME SET, DIFFERENT or FAILED. Exit 0 = SAME SET, 1 = DIFFERENT
// or FAILED, 2 = usage.

#include <array>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <mutex>
#include <process.h>
#include <deque>
#include <string>
#include <thread>
#include <vector>

#include "pwt_paths.hpp"
#include "sha256.hpp"
#include "stamps.hpp"

using namespace cn;

namespace {

constexpr std::size_t kBatchRecs = std::size_t(1) << 16;

struct Cancelled {};

// Bounded queue of record batches from the producer (file A's merge) to the
// consumer (file B's merge). An empty batch marks the end of A.
class BatchQueue {
public:
    bool push(std::vector<std::uint8_t>&& b) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return q_.size() < kDepth || closed_; });
        if (closed_) return false;
        q_.push_back(std::move(b));
        cv_.notify_all();
        return true;
    }
    std::vector<std::uint8_t> pop() {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return !q_.empty(); });
        std::vector<std::uint8_t> b = std::move(q_.front());
        q_.pop_front();
        cv_.notify_all();
        return b;
    }
    void close() {   // consumer gives up: unblock and stop the producer
        std::lock_guard<std::mutex> lk(m_);
        closed_ = true;
        cv_.notify_all();
    }

private:
    static constexpr std::size_t kDepth = 16;
    std::mutex m_;
    std::condition_variable cv_;
    std::deque<std::vector<std::uint8_t>> q_;
    bool closed_ = false;
};

template <int K>
struct Side {
    static constexpr int kRecBytes = K * kPathFieldBytes;
    Sha256 sha;
    u64 records = 0;
    u64 by_len[kPathMaxFields + 1] = {0};
    // Serialize (ascending primes, len) back into the record bytes, checking
    // the record is strictly ascending and nonzero.
    void encode(const u64* primes, int len, std::uint8_t* out) {
        if (len < 1 || len > K) throw std::runtime_error("path length outside [1, fields]");
        for (int i = 0; i < K; ++i) {
            const u64 v = i < len ? primes[i] : 0;
            if (i < len && (v == 0 || (i > 0 && v <= primes[i - 1]) || v >= (u64(1) << 40)))
                throw std::runtime_error("record is not strictly ascending nonzero 40-bit fields");
            for (int q = 0; q < kPathFieldBytes; ++q)
                out[i * kPathFieldBytes + q] = static_cast<std::uint8_t>(v >> (8 * (kPathFieldBytes - 1 - q)));
        }
        sha.update(out, kRecBytes);
        ++records;
        ++by_len[len];
    }
};

template <int K>
std::string rec_str(const std::uint8_t* r) {
    std::string s = "[";
    for (int i = 0; i < K; ++i) {
        u64 v = 0;
        for (int q = 0; q < kPathFieldBytes; ++q) v = (v << 8) | r[i * kPathFieldBytes + q];
        if (v == 0) break;
        if (i) s += ' ';
        s += std::to_string(v);
    }
    return s + "]";
}

template <int K>
int run(const char* fa, const char* fb, const PathSortOptions& base) {
    constexpr int RB = K * kPathFieldBytes;
    PhaseClock clk;
    const std::string pid = std::to_string(_getpid());

    PathSortOptions oa = base, ob = base;
    oa.descending = ob.descending = false;     // keep the ascending record layout
    oa.tag = "pathscmp_" + pid + "_a";
    ob.tag = "pathscmp_" + pid + "_b";

    std::printf("A          %s\n", fa);
    PathRunSource<K> a(fa, oa);
    clk.stamp("sort_a", "external sort of A into runs");
    std::printf("B          %s\n", fb);
    PathRunSource<K> b(fb, ob);
    clk.stamp("sort_b", "external sort of B into runs");

    Side<K> sa, sb;
    BatchQueue q;
    std::exception_ptr perr;
    std::thread producer([&] {
        try {
            std::vector<std::uint8_t> batch;
            batch.reserve(kBatchRecs * RB);
            a.for_each([&](const u64* primes, int len) {
                const std::size_t off = batch.size();
                batch.resize(off + RB);
                sa.encode(primes, len, batch.data() + off);
                if (batch.size() == kBatchRecs * RB) {
                    if (!q.push(std::move(batch))) throw Cancelled{};
                    batch = std::vector<std::uint8_t>();
                    batch.reserve(kBatchRecs * RB);
                }
            });
            if (!batch.empty() && !q.push(std::move(batch))) throw Cancelled{};
            q.push(std::vector<std::uint8_t>());   // end of A
        } catch (const Cancelled&) {
        } catch (...) {
            perr = std::current_exception();
            q.push(std::vector<std::uint8_t>());
        }
    });

    // Consumer: B's sorted stream against A's, merge-join.
    std::vector<std::uint8_t> cur;
    std::size_t pos = 0;
    bool a_done = false;
    auto a_ready = [&]() -> bool {
        while (!a_done && pos >= cur.size()) {
            cur = q.pop();
            pos = 0;
            if (cur.empty()) a_done = true;
        }
        return !a_done;
    };
    u64 common = 0, only_a = 0, only_b = 0;
    std::vector<std::string> ex_a, ex_b;
    try {
        std::uint8_t rb[RB];
        b.for_each([&](const u64* primes, int len) {
            sb.encode(primes, len, rb);
            for (;;) {
                if (!a_ready()) { ++only_b; if (ex_b.size() < 5) ex_b.push_back(rec_str<K>(rb)); return; }
                const int c = std::memcmp(cur.data() + pos, rb, RB);
                if (c < 0) {
                    ++only_a;
                    if (ex_a.size() < 5) ex_a.push_back(rec_str<K>(cur.data() + pos));
                    pos += RB;
                    continue;
                }
                if (c == 0) { ++common; pos += RB; }
                else { ++only_b; if (ex_b.size() < 5) ex_b.push_back(rec_str<K>(rb)); }
                return;
            }
        });
        while (a_ready()) {
            ++only_a;
            if (ex_a.size() < 5) ex_a.push_back(rec_str<K>(cur.data() + pos));
            pos += RB;
        }
    } catch (...) {
        q.close();
        producer.join();
        throw;
    }
    producer.join();
    if (perr) std::rethrow_exception(perr);
    clk.stamp("compare", "lockstep merge of the two sorted streams");

    auto report = [&](const char* tag, Side<K>& s, std::size_t n_file, int max_len) {
        std::printf("%s          %llu records (file %zu), max length %d; sorted-stream sha256 %s\n", tag,
                    (unsigned long long)s.records, n_file, max_len, sha256_hex(s.sha.finish()).c_str());
        std::printf("%s lengths ", tag);
        for (int l = 1; l <= K; ++l)
            if (s.by_len[l]) std::printf(" %d:%llu", l, (unsigned long long)s.by_len[l]);
        std::printf("\n");
    };
    report("A", sa, a.size(), a.max_len());
    report("B", sb, b.size(), b.max_len());
    std::printf("merge      common %llu, only in A %llu, only in B %llu\n", (unsigned long long)common,
                (unsigned long long)only_a, (unsigned long long)only_b);
    for (const auto& s : ex_a) std::printf("  only A   %s\n", s.c_str());
    for (const auto& s : ex_b) std::printf("  only B   %s\n", s.c_str());
    clk.total();
    if (only_a == 0 && only_b == 0 && sa.records == sb.records && common == sa.records) {
        std::printf("SAME SET   %llu paths; every record of A is in B and vice versa (no duplicates in either)\n",
                    (unsigned long long)common);
        return 0;
    }
    std::printf("DIFFERENT  %llu common, %llu only in A, %llu only in B\n", (unsigned long long)common,
                (unsigned long long)only_a, (unsigned long long)only_b);
    return 1;
}

int usage(const char* argv0) {
    std::fprintf(stderr, "usage: %s <a.bin> <b.bin> [--fields 12|14] [--scratch DIR] [--chunk R] [--threads N]\n",
                 argv0);
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const char* fa = nullptr;
    const char* fb = nullptr;
    int fields = 14;
    PathSortOptions opt;
    opt.chunk_recs = std::size_t(1) << 25;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--fields") && i + 1 < argc) fields = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--scratch") && i + 1 < argc) opt.scratch = argv[++i];
        else if (!std::strcmp(argv[i], "--chunk") && i + 1 < argc) opt.chunk_recs = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) opt.threads = static_cast<unsigned>(std::atoi(argv[++i]));
        else if (argv[i][0] == '-' && argv[i][1] == '-') { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return usage(argv[0]); }
        else if (!fa) fa = argv[i];
        else if (!fb) fb = argv[i];
        else { std::fprintf(stderr, "unexpected arg %s\n", argv[i]); return usage(argv[0]); }
    }
    if (!fa || !fb) return usage(argv[0]);
    if (fields != 12 && fields != 14) { std::fprintf(stderr, "--fields must be 12 or 14\n"); return 2; }
    if (opt.chunk_recs == 0) { std::fprintf(stderr, "--chunk must be positive\n"); return 2; }
    try {
        return fields == 14 ? run<14>(fa, fb, opt) : run<12>(fa, fb, opt);
    } catch (const std::exception& e) {
        std::printf("FAILED     %s\n", e.what());
        std::fprintf(stderr, "FAILED: %s\n", e.what());
        return 1;
    }
}
