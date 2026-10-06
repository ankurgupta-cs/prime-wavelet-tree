#pragma once
// pwt_paths.hpp -- prime-path records (the trie_paths output: K five-byte
// big-endian prime fields per record, ascending, zero padded) and an
// external merge sort that delivers them in DESCENDING trie pre-order
// (largest prime first; a prefix sorts before its extensions). The same
// sorted stream can be replayed any number of times (the PWT encoder needs
// two passes), because the sorted runs stay on disk until release().

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "u128.hpp"

namespace cn {

constexpr int kPathFieldBytes = 5;
constexpr int kPathMaxFields = 14;

template <int K> struct PathRec { std::uint8_t b[K * kPathFieldBytes]; };
template <int K> inline bool path_less(const PathRec<K>& a, const PathRec<K>& b) {
    return std::memcmp(a.b, b.b, K * kPathFieldBytes) < 0;
}
template <int K> inline u64 path_field(const PathRec<K>& r, int i) {
    u64 v = 0;
    for (int q = 0; q < kPathFieldBytes; ++q) v = (v << 8) | r.b[i * kPathFieldBytes + q];
    return v;
}
template <int K> inline void path_set_field(PathRec<K>& r, int i, u64 v) {
    for (int q = 0; q < kPathFieldBytes; ++q)
        r.b[i * kPathFieldBytes + q] = static_cast<std::uint8_t>(v >> (8 * (kPathFieldBytes - 1 - q)));
}
template <int K> inline int path_len(const PathRec<K>& r) {
    int n = 0;
    while (n < K && path_field<K>(r, n) != 0) ++n;
    return n;
}

struct PathSortOptions {
    std::string scratch = "data/scratch";
    std::string tag = "pwt";                 // run-file name prefix (make it unique per run)
    std::size_t chunk_recs = std::size_t(1) << 26;
    unsigned threads = 0;
    bool descending = true;
    bool verbose = true;
};

// Sorted, replayable stream of paths. Records are delivered as
// (primes, len) with primes[0] the LARGEST prime when descending.
template <int K>
class PathRunSource {
public:
    using Rec = PathRec<K>;
    static constexpr int kRecBytes = K * kPathFieldBytes;

    PathRunSource(const std::string& input, PathSortOptions opt) : opt_(std::move(opt)) {
        if (opt_.threads == 0) {
            const unsigned hc = std::thread::hardware_concurrency();
            opt_.threads = hc > 2 ? hc - 1 : 1;
        }
        std::ifstream in(input, std::ios::binary | std::ios::ate);
        if (!in) throw std::runtime_error("cannot open " + input);
        const std::streamsize fsize = in.tellg();
        if (fsize % kRecBytes) throw std::runtime_error("truncated path record (wrong --fields?)");
        n_ = static_cast<std::size_t>(fsize / kRecBytes);
        in.seekg(0);
        std::vector<Rec> chunk;
        chunk.reserve(std::min(opt_.chunk_recs, n_));
        std::size_t done = 0;
        while (done < n_) {
            const std::size_t take = std::min(opt_.chunk_recs, n_ - done);
            chunk.resize(take);
            std::size_t got = 0;
            while (got < take) {   // <= 512 MB per read (multi-GB reads fail on this CRT)
                const std::size_t piece = std::min<std::size_t>(take - got, (std::size_t(512) << 20) / kRecBytes);
                in.read(reinterpret_cast<char*>(chunk.data() + got), static_cast<std::streamsize>(piece * kRecBytes));
                if (in.gcount() != static_cast<std::streamsize>(piece * kRecBytes)) throw std::runtime_error("short read");
                got += piece;
            }
            for (Rec& r : chunk) {
                const int len = path_len<K>(r);
                if (len == 0) throw std::runtime_error("empty path record");
                if (opt_.descending)
                    for (int i = 0; i < len / 2; ++i) {
                        const u64 a = path_field<K>(r, i), b = path_field<K>(r, len - 1 - i);
                        path_set_field<K>(r, i, b);
                        path_set_field<K>(r, len - 1 - i, a);
                    }
                for (int i = 0; i < len; ++i) {          // per-depth maximum prime (depth i+1)
                    const u64 p = path_field<K>(r, i);
                    if (p > max_at_depth_[i + 1]) max_at_depth_[i + 1] = p;
                }
                if (len > max_len_) max_len_ = len;
            }
            parallel_sort_(chunk);
            const std::string rp = opt_.scratch + "/" + opt_.tag + "_run_" + std::to_string(runs_.size()) + ".bin";
            if (std::ifstream(rp, std::ios::binary))   // never clobber another process's run file
                throw std::runtime_error("run file already exists (another run sharing --scratch and tag?): " + rp);
            std::ofstream out(rp, std::ios::binary | std::ios::trunc);
            if (!out) throw std::runtime_error("cannot create " + rp);
            for (std::size_t w = 0; w < take;) {
                const std::size_t piece = std::min<std::size_t>(take - w, (std::size_t(512) << 20) / kRecBytes);
                out.write(reinterpret_cast<const char*>(chunk.data() + w), static_cast<std::streamsize>(piece * kRecBytes));
                w += piece;
            }
            out.close();
            if (!out) throw std::runtime_error("write failure " + rp);
            runs_.push_back(rp);
            done += take;
            if (opt_.verbose) std::printf("  run %zu: %zu records sorted\n", runs_.size(), take);
        }
    }
    ~PathRunSource() { release(); }
    PathRunSource(const PathRunSource&) = delete;
    PathRunSource& operator=(const PathRunSource&) = delete;

    std::size_t size() const { return n_; }
    int max_len() const { return max_len_; }
    u64 max_prime_at_depth(int d) const { return (d >= 1 && d <= K) ? max_at_depth_[d] : 0; }

    // Replay the sorted stream. f(const u64* primes, int len). Throws on a
    // duplicate or mis-sorted record.
    template <typename F>
    void for_each(F&& f) const {
        std::vector<RunReader> rd(runs_.size());
        std::vector<std::size_t> heap;
        for (std::size_t i = 0; i < runs_.size(); ++i) {
            rd[i].in.open(runs_[i], std::ios::binary);
            if (!rd[i].in) throw std::runtime_error("cannot reopen " + runs_[i]);
            if (rd[i].ready()) heap.push_back(i);
        }
        auto greater = [&](std::size_t a, std::size_t b) { return path_less<K>(rd[b].cur(), rd[a].cur()); };
        std::make_heap(heap.begin(), heap.end(), greater);
        Rec prev{};
        std::memset(prev.b, 0, kRecBytes);
        std::size_t scanned = 0;
        u64 primes[K];
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end(), greater);
            const std::size_t ri = heap.back();
            heap.pop_back();
            const Rec r = rd[ri].cur();
            rd[ri].advance();
            if (rd[ri].ready()) { heap.push_back(ri); std::push_heap(heap.begin(), heap.end(), greater); }
            if (scanned > 0 && !path_less<K>(prev, r)) throw std::runtime_error("duplicate or unsorted path");
            const int len = path_len<K>(r);
            for (int i = 0; i < len; ++i) primes[i] = path_field<K>(r, i);
            f(static_cast<const u64*>(primes), len);
            prev = r;
            ++scanned;
        }
        if (scanned != n_) throw std::runtime_error("merge delivered a wrong record count");
    }

    void release() {
        for (const std::string& rp : runs_) std::remove(rp.c_str());
        runs_.clear();
    }

private:
    struct RunReader {
        std::ifstream in;
        std::vector<Rec> buf;
        std::size_t pos = 0, len = 0;
        bool eof = false;
        bool fill() {
            if (eof) return false;
            buf.resize(std::size_t(1) << 18);
            in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size() * sizeof(Rec)));
            len = static_cast<std::size_t>(in.gcount()) / sizeof(Rec);
            pos = 0;
            if (len == 0) { eof = true; return false; }
            return true;
        }
        bool ready() { return pos < len || fill(); }
        const Rec& cur() const { return buf[pos]; }
        void advance() { ++pos; }
    };

    void parallel_sort_(std::vector<Rec>& chunk) const {
        const std::size_t take = chunk.size();
        const std::size_t T = std::max<std::size_t>(1, std::min<std::size_t>(opt_.threads, take));
        const std::size_t part = (take + T - 1) / T;
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t lo = t * part, hi = std::min(take, lo + part);
            if (lo >= hi) break;
            pool.emplace_back([&chunk, lo, hi] { std::sort(chunk.begin() + lo, chunk.begin() + hi, path_less<K>); });
        }
        for (auto& th : pool) th.join();
        for (std::size_t width = part; width < take; width *= 2)
            for (std::size_t lo = 0; lo + width < take; lo += 2 * width)
                std::inplace_merge(chunk.begin() + lo, chunk.begin() + lo + width,
                                   chunk.begin() + std::min(take, lo + 2 * width), path_less<K>);
    }

    PathSortOptions opt_;
    std::size_t n_ = 0;
    int max_len_ = 0;
    u64 max_at_depth_[kPathMaxFields + 2] = {0};
    std::vector<std::string> runs_;
};

// In-memory source for tests: records as vectors of primes in ANY order;
// the source sorts them descending like the disk version.
class PathVectorSource {
public:
    explicit PathVectorSource(std::vector<std::vector<u64>> paths) : paths_(std::move(paths)) {
        for (auto& p : paths_) {
            std::sort(p.begin(), p.end(), std::greater<u64>());
            if (p.empty()) throw std::runtime_error("empty path");
            for (std::size_t i = 0; i < p.size(); ++i)
                if (p[i] > max_at_depth_[i + 1]) max_at_depth_[i + 1] = p[i];
            if (static_cast<int>(p.size()) > max_len_) max_len_ = static_cast<int>(p.size());
        }
        std::sort(paths_.begin(), paths_.end());   // lexicographic = trie pre-order
        for (std::size_t i = 1; i < paths_.size(); ++i)
            if (paths_[i] == paths_[i - 1]) throw std::runtime_error("duplicate path");
    }
    std::size_t size() const { return paths_.size(); }
    int max_len() const { return max_len_; }
    u64 max_prime_at_depth(int d) const { return (d >= 1 && d <= kPathMaxFields) ? max_at_depth_[d] : 0; }
    template <typename F>
    void for_each(F&& f) const {
        for (const auto& p : paths_) f(p.data(), static_cast<int>(p.size()));
    }

private:
    std::vector<std::vector<u64>> paths_;
    int max_len_ = 0;
    u64 max_at_depth_[kPathMaxFields + 2] = {0};
};

} // namespace cn
