# CND: sorted lists coded by residue class (interpolative coding)

Two cells use this layout:

- **cell 101 (CND1 proper):** the d_min list, M = 2310 = 2 * 3 * 5 * 7 * 11;
- **cell 001:** the n list, M = 10,810,800, in the same layout.

**The modulus is read from the header.** Version 1 of this format was first
specified with M = 2310 fixed, and `cnd_decode` (the original CND1 decoder)
still accepts only that value. `list_decode` reads any M >= 1 from the
header and applies the even-class rule (below) when M is even; `list_build
classes --mod M` writes any M. The bytes of both files are exactly as
described here; the only field that differs is M (and the values stored).

Implementation: `src/bic.hpp`, `src/cnd1.hpp`, `src/cnd_encode.cpp`,
`src/cnd_decode.cpp`, `src/list_build.cpp`, `src/list_decode.cpp`.

## What it stores

A sorted set of distinct odd values v: the Carmichael numbers n themselves
(cell 001), or one divisor per number (cell 101): d_min(n), the smallest
divisor d | n with lambda(d) > n / d, from which n = d * (d^{-1} mod
lambda(d)); the map n -> d_min(n) is injective. The set is partitioned by
residue class r = v mod M; class membership is recomputable from the value,
so the partition costs nothing at decode time and prices the residue bias of
the values (all odd; for d_min, H(d mod 3) = 0.995; squarefree exclusions).

Within class r the stored values are k = (v - r) / M, a strictly increasing
sequence coded with binary interpolative coding (BIC).

## Number formats

- Fixed-width integers: little-endian.
- VByte: base-128, low 7 bits first, high bit = continuation (`src/vbyte.hpp`).
- Bitstreams: bits are consumed MSB-first within each byte; each class blob
  starts byte-aligned and is zero-padded to a byte boundary at its end.

## Layout

    [header 64 B] [directory] [body] [footer 88 B]

### Header (64 bytes)

    off size field
      0    4 magic "CND1"
      4    4 version        u32 = 1
      8    4 M              u32   the modulus (2310 for cell 101,
                                  10,810,800 for cell 001)
     12    4 end_width      u32   bit width of the two absolute ends of
                                  every class blob (global: bits(max k))
     16    8 record_count   u64
     24    8 dir_off        u64   = 64
     32    8 body_off       u64
     40    8 footer_off     u64
     48   16 reserved       zero

### Directory

For every class r = 0 .. M-1 in order:

    VByte count_r                  number of values with v mod M = r
    VByte blob_bytes_r             ONLY present when count_r > 0: byte
                                   length of the class's body blob

Byte offsets of blobs are the running sum of blob_bytes_r from body_off;
classes are therefore independently seekable and decodable in parallel.

### Body

Per nonempty class, in class order, one byte-aligned blob:

    k_first   end_width bits       (absent when count_r = 0)
    k_last    end_width bits       (absent when count_r < 2)
    interior  BIC bits             (absent when count_r < 3)
    padding   0-7 zero bits

### BIC (normative)

The interior elements ks[1..c-2] of a class (c = count_r) are coded by
midpoint recursion over frames (ilo, ihi, vlo, vhi), where every element
with index in [ilo, ihi] is strictly between the values vlo and vhi:

    initial frame: (1, c-2, ks[0], ks[c-1])
    mid  = ilo + (ihi - ilo) / 2
    lo   = vlo + 1 + (mid - ilo)         smallest value ks[mid] could take
    hi   = vhi - 1 - (ihi - mid)         largest
    S    = hi - lo + 1
    emit centered minimal binary of (ks[mid] - lo) within span S
    recurse LEFT  (ilo, mid-1, vlo, ks[mid])   if mid > ilo
    recurse RIGHT (mid+1, ihi, ks[mid], vhi)   if mid < ihi

Traversal order is strict pre-order, LEFT before RIGHT. Frames with
ilo > ihi never arise (guarded by the two recursion conditions).

Centered minimal binary of offset v in span S (v in [0, S)):

    S = 1: emit nothing.
    L = ceil(log2 S); a = 2^L - S (count of short codewords);
    b = S - a (count of long codewords)
    v' = (v - floor(b/2)) mod S          rotate: middle offsets get short codes
    if v' < a:  emit v' in L-1 bits
    else:       emit v' + a in L bits

Prefix-freeness: L-bit codewords occupy [2a, 2^L); the decoder reads L-1
bits, accepts if the value is < a, otherwise appends one more bit and
subtracts a. The decoder recomputes lo/hi/S from its identical frame
arithmetic, then inverts the rotation.

### Footer (88 bytes)

    off size field
      0   32 sha_payload    SHA-256 of header || directory || body (the
                            whole file up to footer_off; a tampered header
                            field such as end_width must not decode)
     32   32 sha_nset       SHA-256 of the certification stream: the
                            Carmichael numbers n sorted ascending, each as
                            16 bytes little-endian (for cell 101 the n
                            rebuilt from the stored d)
     64   16 total_check    sum of all n mod 2^128 (equals the ORC1
                            footer's total_check for the same table)
     80    8 record_count   u64, repeated

## Decoder obligations

1. Validate the header before allocating from any of it: magic, version,
   M (cnd_decode: M = 2310 exactly; list_decode: 1 <= M <= the directory
   size in bytes, since every class owns at least one directory byte),
   end_width in [1, 64], offsets consistent with the file size.
2. Directory bounds, checked per class so no u64 accumulator can wrap:
   count_r <= record_count - (sum so far); count_r <= 2^end_width when
   end_width < 64 (a class cannot hold more values than its universe);
   blob_bytes_r <= footer_off - (running offset); blob_bytes_r >=
   ceil(min(count_r, 2) * end_width / 8). The sum of count_r must equal
   record_count and the last blob must end exactly at footer_off.
3. Within a class: k_first < k_last when count_r >= 2, decoded ks strictly
   increasing (the frame arithmetic enforces this when every span S >= 1;
   a span < 1 or a read past blob_bytes_r is a container error).
4. Exact consumption: a blob must contain exactly the bits the recursion
   demands; ceil(bits / 8) must equal blob_bytes_r and every pad bit must
   be zero. sha_payload must match.
5. v = M * k + r; require v odd. With even M this is exactly: every
   nonempty class has odd r.

cnd_decode accepts non-canonical (zero-padded) VByte encodings in the
directory, as ORC1 does; the writer never emits them, and sha_payload
distinguishes the byte streams. list_decode rejects them.

## Certification (the accepted proof of correctness)

Decode every class. For cell 101, for each d: factor d (`src/splitter.hpp`;
d < 2^70 < psi_12, so the 12-base Miller-Rabin verdicts are proofs),
lambda(d) = lcm(p - 1), require gcd(d, lambda(d)) = 1, and set
n = d * (d^{-1} mod lambda(d)). Sort the n ascending; then

- the SHA-256 of the 16-byte-LE n stream must equal the footer's sha_nset,
- the count and every element must match the n sequence streamed from the
  oracle (`docs/ORACLE_FORMAT.md`), element by element (never the text),
- total_check must equal the oracle's.

`cnd_decode --certify` (M = 2310) and `list_decode --kind n|d --certify`
(any M) perform all three; the factoring parallelizes by index ranges.
