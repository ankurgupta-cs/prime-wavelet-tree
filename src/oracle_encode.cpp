// oracle_encode -- new_table.txt (or any prefix/fragment with complete lines)
// -> ORC1. Streams the input once in fixed-size chunks: SHA-256 runs over the
// raw bytes while lines are parsed and verified. Verification per line:
// strictly ascending n, k in [3,14], strictly ascending factors, n == product,
// Korselt lambda(n) | n-1, factors < 10^12, p1 < 10^8. --prove-primes adds a
// deterministic Miller-Rabin proof for every factor.
//
// usage: oracle_encode <input.txt> <output.orc1> [--prove-primes]

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "mr64.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "u128.hpp"

using namespace cn;

namespace {

u128 gcd128(u128 a, u128 b) {
    while (b != 0) { const u128 t = a % b; a = b; b = t; }
    return a;
}

struct LineStats {
    u64 lines = 0;
    long long text_bytes = 0;
};

// Parses and verifies one line (without trailing LF); appends to writer.
void handle_line(const char* s, std::size_t len, u64 lineno, bool prove,
                 Orc1Writer& writer) {
    auto fail = [&](const char* why) {
        std::fprintf(stderr, "line %llu: %s: %.*s\n",
                     static_cast<unsigned long long>(lineno), why,
                     static_cast<int>(len > 80 ? 80 : len), s);
        throw std::runtime_error("input verification failed");
    };
    if (len == 0) fail("empty line");

    const char* p = s;
    const char* const e = s + len;
    auto read_number = [&](u128& v) {
        if (p >= e || *p < '0' || *p > '9') fail("expected digit");
        if (*p == '0' && p + 1 < e && p[1] >= '0' && p[1] <= '9')
            fail("leading zero (non-canonical decimal)");
        v = 0;
        while (p < e && *p >= '0' && *p <= '9') {
            const unsigned d = static_cast<unsigned>(*p - '0');
            if (v > (U128_MAX - d) / 10) fail("number overflows u128");
            v = v * 10 + d;
            ++p;
        }
    };

    static const u128 kBound = parse_u128("1000000000000000000000000");  // 10^24
    u128 n;
    read_number(n);
    if (n >= kBound) fail("n >= 10^24");
    u64 factors[kOrc1MaxFactors];
    int k = 0;
    u128 prod = 1;
    while (p < e) {
        if (*p != ' ') fail("expected single space");
        ++p;
        if (k >= kOrc1MaxFactors) fail("more than 14 factors");
        u128 f;
        read_number(f);
        if (f < 3) fail("factor < 3");
        if (f >= 1'000'000'000'000ULL) fail("factor >= 10^12");
        factors[k] = static_cast<u64>(f);
        if (k > 0 && factors[k] <= factors[k - 1]) fail("factors not ascending");
        if (prod > kBound / factors[k]) fail("factor product exceeds 10^24");
        prod *= factors[k];
        ++k;
    }
    if (k < kOrc1MinFactors) fail("fewer than 3 factors");
    if (factors[0] >= 100'000'000ULL) fail("smallest factor >= 10^8");
    if (prod != n) fail("n != product of factors");

    // Korselt: lcm(p_i - 1) divides n - 1.
    u128 lam = 1;
    for (int i = 0; i < k; ++i) {
        const u128 pm1 = factors[i] - 1;
        lam = lam / gcd128(lam, pm1) * pm1;
    }
    if ((n - 1) % lam != 0) fail("Korselt violation");

    if (prove)
        for (int i = 0; i < k; ++i)
            if (!is_prime_u64(factors[i])) fail("factor not prime");

    writer.add(n, factors, k);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <input.txt> <output.orc1> [--prove-primes]\n",
                     argv[0]);
        return 2;
    }
    const bool prove = argc > 3 && std::strcmp(argv[3], "--prove-primes") == 0;

    std::ifstream in(argv[1], std::ios::binary);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    std::ofstream out(argv[2], std::ios::binary | std::ios::trunc);
    if (!out) { std::fprintf(stderr, "cannot open %s\n", argv[2]); return 2; }

    try {
        Orc1Writer writer(out);
        Sha256 sha_src;
        LineStats st;

        constexpr std::size_t kChunk = 8u << 20;
        std::vector<char> chunk(kChunk);
        std::string carry;
        for (;;) {
            in.read(chunk.data(), static_cast<std::streamsize>(kChunk));
            const std::size_t got = static_cast<std::size_t>(in.gcount());
            if (got == 0) break;
            sha_src.update(chunk.data(), got);
            st.text_bytes += static_cast<long long>(got);

            std::size_t start = 0;
            for (std::size_t i = 0; i < got; ++i) {
                if (chunk[i] == '\r') {
                    std::fprintf(stderr, "CR byte at offset ~%lld: not LF-only\n",
                                 st.text_bytes - static_cast<long long>(got - i));
                    return 1;
                }
                if (chunk[i] != '\n') continue;
                ++st.lines;
                if (!carry.empty()) {
                    carry.append(chunk.data() + start, i - start);
                    handle_line(carry.data(), carry.size(), st.lines, prove, writer);
                    carry.clear();
                } else {
                    handle_line(chunk.data() + start, i - start, st.lines, prove,
                                writer);
                }
                start = i + 1;
            }
            carry.append(chunk.data() + start, got - start);
            if (carry.size() > 4096) {
                std::fprintf(stderr,
                             "line exceeds 4096 bytes: not the table format\n");
                return 1;
            }
        }
        if (!carry.empty()) {
            std::fprintf(stderr, "input does not end with LF\n");
            return 1;
        }

        const Orc1Footer f = writer.finish(sha_src.finish());
        out.seekp(0, std::ios::end);           // finish() leaves tellp at the header
        const long long orc_bytes =
            static_cast<long long>(out.tellp());
        out.close();
        if (!out || orc_bytes <= 0) {
            std::fprintf(stderr, "output write/close failed\n");
            return 1;
        }
        std::printf("records    %llu\n",
                    static_cast<unsigned long long>(writer.records_written()));
        std::printf("blocks     %llu\n", static_cast<unsigned long long>(f.n_blocks));
        std::printf("text bytes %lld\n", st.text_bytes);
        std::printf("orc1 bytes %lld  (%.2fx)\n", orc_bytes,
                    orc_bytes > 0 ? double(st.text_bytes) / double(orc_bytes) : 0.0);
        std::printf("sha256 src  %s\n", sha256_hex(f.sha_src).c_str());
        std::printf("sha256 body %s\n", sha256_hex(f.sha_body).c_str());
        std::printf("primes     %s\n", prove ? "proven (MR-12)" : "not re-proven");
        return 0;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 1;
    }
}
