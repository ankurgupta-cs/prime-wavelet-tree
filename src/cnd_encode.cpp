// cnd_encode -- build the CND1 container (docs/CND_FORMAT.md) from the
// emitted sorted d-set, taking the certification targets (n-stream sha,
// mod-2^128 sum) from the certified ORC1 oracle -- never from the text.
//
// usage: cnd_encode <d.raw> <table.orc1> <out.cnd1>

#include <cstdio>
#include <fstream>
#include <vector>

#include "cnd1.hpp"
#include "orc1.hpp"
#include "sha256.hpp"
#include "u128.hpp"

using namespace cn;

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <d.raw> <table.orc1> <out.cnd1>\n", argv[0]);
        return 2;
    }
    try {
        std::ifstream in(argv[1], std::ios::binary | std::ios::ate);
        if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
        const u64 fsize = static_cast<u64>(in.tellg());
        if (fsize % 16) { std::fprintf(stderr, "truncated d record\n"); return 1; }
        const std::size_t nrec = static_cast<std::size_t>(fsize / 16);
        in.seekg(0);
        std::vector<u128> ds(nrec);
        std::vector<std::uint8_t> chunk(16 * 65536);
        std::size_t w = 0;
        while (in && w < nrec) {
            in.read(reinterpret_cast<char*>(chunk.data()),
                    static_cast<std::streamsize>(chunk.size()));
            const std::streamsize got = in.gcount();
            for (std::streamsize off = 0; off < got; off += 16) {
                u128 d = 0;
                for (int b = 15; b >= 0; --b) d = (d << 8) | chunk[off + b];
                ds[w++] = d;
            }
        }
        if (w != nrec) { std::fprintf(stderr, "short read\n"); return 1; }

        // certification targets from the oracle
        Cnd1FooterTargets targets;
        Sha256 sha;
        u128 sum = 0;
        u64 orc_records = 0;
        {
            Orc1Reader rd(argv[2]);
            rd.for_all([&](const Orc1Record& rec) {
                std::uint8_t le[16];
                for (int i = 0; i < 16; ++i)
                    le[i] = static_cast<std::uint8_t>(rec.n >> (8 * i));
                sha.update(le, 16);
                sum += rec.n;
                ++orc_records;
            });
        }
        if (orc_records != nrec) {
            std::fprintf(stderr, "record mismatch: oracle %llu vs d-set %zu\n",
                         (unsigned long long)orc_records, nrec);
            return 1;
        }
        targets.sha_nset = sha.finish();
        targets.total_check = sum;

        std::ofstream out(argv[3], std::ios::binary | std::ios::trunc);
        if (!out) { std::fprintf(stderr, "cannot write %s\n", argv[3]); return 2; }
        const u64 total = cnd1_encode(ds, kCnd1M, targets, out);
        out.close();
        if (!out) { std::fprintf(stderr, "write failure\n"); return 1; }

        std::printf("records    %zu\n", nrec);
        std::printf("cnd1 bytes %llu  (%.3f bits/el, %.2fx vs 18,358,310,309 B text)\n",
                    (unsigned long long)total,
                    double(total) * 8.0 / double(nrec),
                    18358310309.0 / double(total));
        std::printf("sha nset   %s\n", sha256_hex(targets.sha_nset).c_str());
        std::printf("total_chk  %s\n", to_string(targets.total_check).c_str());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
