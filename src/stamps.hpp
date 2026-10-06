#pragma once
// stamps.hpp -- phase stamps and bulk 16-byte-record I/O for the timed
// construction tools (dlist_from_paths, list_build, list_decode, paths_cmp,
// nstats_orc1 --emit-only).
//
// Stamp lines are machine-readable, one per phase:
//   stamp  <phase>  <seconds> s   cum <seconds> s   <note>
// <seconds> is the phase's own wall time (steady_clock); cum is the wall
// time since the PhaseClock was constructed (first statement of main), read
// when the line is printed. A phase measured as a sum over interleaved
// chunks (e.g. read vs multiply in one streaming loop) is printed with
// add(); its seconds are that sum, so consecutive add() lines need not
// differ in cum by exactly their seconds. total() prints
//   stamp  total  <seconds> s   (phases sum <s> s, unaccounted <s> s)
// where "unaccounted" is loop/allocation overhead outside every phase.
// Writes go to the OS (stream closed); no fsync, as in every other tool.

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "u128.hpp"

namespace cn {

static_assert(std::endian::native == std::endian::little,
              "16-byte LE record I/O copies u128 memory directly: little-endian host required");
static_assert(sizeof(u128) == 16, "u128 must be 16 bytes");

inline double wall_now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

class PhaseClock {
public:
    PhaseClock() : t0_(wall_now()), last_(t0_) {}

    // Phase from the previous stamp (or construction) to now.
    double stamp(const char* phase, const char* note = "") {
        const double t = wall_now();
        const double dt = t - last_;
        last_ = t;
        sum_ += dt;
        print_(phase, dt, t - t0_, note);
        return dt;
    }

    // Phase measured elsewhere (a sum over chunks); resets the interval start
    // so the next stamp() does not count this time again.
    void add(const char* phase, double secs, const char* note = "") {
        const double t = wall_now();
        last_ = t;
        sum_ += secs;
        print_(phase, secs, t - t0_, note);
    }

    // Restart the interval without printing (time outside every phase).
    void skip() { last_ = wall_now(); }

    double elapsed() const { return wall_now() - t0_; }

    double total() {
        const double t = wall_now() - t0_;
        std::printf("stamp  %-14s %10.3f s   (phases sum %.3f s, unaccounted %.3f s)\n",
                    "total", t, sum_, t - sum_);
        return t;
    }

private:
    static void print_(const char* phase, double dt, double cum, const char* note) {
        std::printf("stamp  %-14s %10.3f s   cum %10.3f s   %s\n", phase, dt, cum, note);
    }

    double t0_, last_;
    double sum_ = 0.0;
};

inline constexpr std::size_t kIoChunk = std::size_t(64) << 20;   // 64 MiB per read/write call

// Chunked write of len bytes (single huge write calls are avoided on Windows).
inline void write_all(std::ofstream& out, const void* data, std::size_t len) {
    const char* p = static_cast<const char*>(data);
    while (len > 0) {
        const std::size_t take = std::min(len, kIoChunk);
        out.write(p, static_cast<std::streamsize>(take));
        if (!out) throw std::runtime_error("write failure");
        p += take;
        len -= take;
    }
}

// Chunked read of exactly len bytes.
inline void read_exact(std::ifstream& in, void* data, std::size_t len) {
    char* p = static_cast<char*>(data);
    while (len > 0) {
        const std::size_t take = std::min(len, kIoChunk);
        in.read(p, static_cast<std::streamsize>(take));
        if (static_cast<std::size_t>(in.gcount()) != take) throw std::runtime_error("short read");
        p += take;
        len -= take;
    }
}

inline std::uint64_t file_size_of(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("cannot open " + path);
    return static_cast<std::uint64_t>(in.tellg());
}

// Load a file of 16-byte little-endian records into v (resized): chunked
// reads straight into v's storage (the host is little-endian).
inline void raw16_load(const std::string& path, std::vector<u128>& v) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("cannot open " + path);
    const std::uint64_t fsize = static_cast<std::uint64_t>(in.tellg());
    if (fsize % 16) throw std::runtime_error(path + ": size is not a multiple of 16");
    in.seekg(0);
    v.resize(static_cast<std::size_t>(fsize / 16));
    read_exact(in, v.data(), static_cast<std::size_t>(fsize));
}

// Write v as 16-byte little-endian records (truncating), then close.
inline void raw16_write(const std::string& path, const std::vector<u128>& v) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot write " + path);
    write_all(out, v.data(), v.size() * 16);
    out.close();
    if (!out) throw std::runtime_error("write failure on " + path);
}

// Whole file into a byte vector (verification passes).
inline std::vector<std::uint8_t> file_bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("cannot open " + path);
    const std::uint64_t fsize = static_cast<std::uint64_t>(in.tellg());
    in.seekg(0);
    std::vector<std::uint8_t> buf(static_cast<std::size_t>(fsize));
    read_exact(in, buf.data(), buf.size());
    return buf;
}

} // namespace cn
