// oracle_decode -- ORC1 -> exact text reconstruction, with the round-trip
// certificate: the SHA-256 of the regenerated stream must equal the sha_src
// recorded at encode time. Also re-verifies sha_body (raw body||index bytes)
// and the total_check sum before declaring success.
//
// usage: oracle_decode <input.orc1> <output.txt|->     ("-" = stdout, no cert file write)

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "orc1.hpp"
#include "sha256.hpp"
#include "u128.hpp"

using namespace cn;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <input.orc1> <output.txt|->\n", argv[0]);
        return 2;
    }
#ifdef _WIN32
    // The reconstruction must be byte-exact: without this, the CRT's text
    // mode would rewrite every LF as CRLF on the "-" (stdout) path while the
    // certificate hash is computed over the untranslated bytes.
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    try {
        // Pass 1: sha_body over raw [header end, index end).
        {
            std::ifstream raw(argv[1], std::ios::binary);
            if (!raw) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
            std::uint8_t hb[kOrc1HeaderBytes];
            raw.read(reinterpret_cast<char*>(hb), sizeof hb);
            if (raw.gcount() != static_cast<std::streamsize>(sizeof hb)) {
                std::fprintf(stderr, "truncated header\n");
                return 1;
            }
            const u64 footer_off = get_le64(hb + 40);
            Sha256 body;
            std::vector<char> buf(8u << 20);
            u64 remaining = footer_off - kOrc1HeaderBytes;
            while (remaining > 0) {
                const std::size_t want =
                    remaining < buf.size() ? static_cast<std::size_t>(remaining)
                                           : buf.size();
                raw.read(buf.data(), static_cast<std::streamsize>(want));
                if (static_cast<std::size_t>(raw.gcount()) != want) {
                    std::fprintf(stderr, "truncated body/index\n");
                    return 1;
                }
                body.update(buf.data(), want);
                remaining -= want;
            }
            Orc1Reader probe(argv[1]);
            if (body.finish() != probe.footer().sha_body) {
                std::fprintf(stderr, "sha_body MISMATCH: container corrupt\n");
                return 1;
            }
        }

        Orc1Reader reader(argv[1]);
        const bool to_stdout = std::strcmp(argv[2], "-") == 0;
        std::ofstream file_out;
        if (!to_stdout) {
            file_out.open(argv[2], std::ios::binary | std::ios::trunc);
            if (!file_out) { std::fprintf(stderr, "cannot open %s\n", argv[2]); return 2; }
        }

        Sha256 sha_out;
        u128 total = 0;
        u64 records = 0;
        bool write_ok = true;
        std::string buf;
        buf.reserve(16u << 20);
        auto flush = [&](bool force) {
            if (!force && buf.size() < (8u << 20)) return;
            sha_out.update(buf.data(), buf.size());
            if (to_stdout) {
                if (std::fwrite(buf.data(), 1, buf.size(), stdout) != buf.size())
                    write_ok = false;
            } else {
                file_out.write(buf.data(),
                               static_cast<std::streamsize>(buf.size()));
                if (!file_out) write_ok = false;
            }
            buf.clear();
        };

        reader.for_all([&](const Orc1Record& r) {
            buf += to_string(r.n);
            for (int i = 0; i < r.k; ++i) {
                buf += ' ';
                buf += to_string(static_cast<u128>(r.factors[i]));
            }
            buf += '\n';
            total += r.n;
            ++records;
            flush(false);
        });
        flush(true);
        if (to_stdout) {
            if (std::fflush(stdout) != 0) write_ok = false;
        } else {
            file_out.close();
            if (!file_out) write_ok = false;
        }

        const auto digest = sha_out.finish();
        const bool sha_ok = digest == reader.footer().sha_src;
        const bool cnt_ok = records == reader.header().n_records;
        const bool sum_ok = total == reader.footer().total_check;
        const bool ok = sha_ok && cnt_ok && sum_ok && write_ok;
        std::fprintf(stderr, "records   %llu (%s)\n",
                     static_cast<unsigned long long>(records),
                     cnt_ok ? "ok" : "MISMATCH");
        std::fprintf(stderr, "total_chk %s\n", sum_ok ? "ok" : "MISMATCH");
        std::fprintf(stderr, "output    %s\n",
                     write_ok ? "ok" : "WRITE FAILED (truncated?)");
        std::fprintf(stderr, "sha256    %s\n", sha256_hex(digest).c_str());
        std::fprintf(stderr, "expected  %s\n",
                     sha256_hex(reader.footer().sha_src).c_str());
        std::fprintf(stderr, "round-trip %s\n", ok ? "CERTIFIED" : "FAILED");
        return ok ? 0 : 1;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 1;
    }
}
