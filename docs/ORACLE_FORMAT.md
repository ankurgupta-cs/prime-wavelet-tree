# ORC1: the oracle format

The *oracle* is an exact, self-indexing binary encoding of the Shallue-Webster
text table `new_table.txt`. It has three jobs: regenerate the text byte for
byte, serve as the ground truth that every other file in this repository is
certified against, and answer simple queries by block. It is deliberately
simple; it is not one of the compressed representations.

Implementation: `src/orc1.hpp` (writer and reader), `src/oracle_encode.cpp`
(text to ORC1), `src/oracle_decode.cpp` (ORC1 to text, with the round-trip
certificate). Size, hash and commands: `oracle/ORACLE`, `oracle/TEXT`.

## 1. Source file

- `https://blue.butler.edu/~jewebste/new_table.txt` (size and SHA-256 in
  `oracle/TEXT`); the server supports HTTP byte ranges.
- One line per Carmichael number: `n SP p1 SP ... SP pk LF` (LF only, the file
  ends in LF, no header). n ascends across lines; the factors ascend strictly
  within a line; n = p1 * ... * pk exactly; k in [3, 14]; every pi < 10^12;
  p1 < 10^8.

## 2. Record layout (factors only)

The factor tuple is the whole information content; n is reconstructed as the
product. Per Carmichael number, in file order:

```
record := k:u8 | VByte(p1) | VByte(p2 - p1) | ... | VByte(pk - p(k-1))
```

- `k` = factor count, valid range [3, 14] (values 0-2 and 15-255 are
  reserved; readers reject them).
- VByte: little-endian base-128, low 7 bits per byte, high bit =
  continuation (`src/vbyte.hpp`).
- Factor deltas are >= 1 (the encoder asserts strict ascent).
- No per-record n and no cross-record state: any record decodes given its
  start offset, and n = product(pi) is recomputed in 128-bit arithmetic
  (n < 2^80).

## 3. Container layout

```
file := header | body | index | footer
header (64 B, all fixed-width little-endian):
    magic       8 B   "CNORC1\0\0"
    version     u32   = 1
    flags       u32   = 0 (reserved)
    n_records   u64
    block_size  u32   = 4096 (records per block)
    reserved    u32
    index_off   u64   byte offset of the index
    footer_off  u64   byte offset of the footer
    reserved    16 B
body := concatenated records, grouped in blocks of block_size records
    (the last block ragged). Blocks are byte-contiguous, without padding.
    Every block starts a fresh record, so a block decodes independently.
index := n_blocks entries of
    first_n     u128  (16 B LE)  n of the block's first record
    byte_off    u64   offset of the block's first record from the body start
    rank_base   u64   global record index of the block's first record
footer (88 B):
    sha256_src  32 B  SHA-256 of the original text file
    sha256_body 32 B  SHA-256 of body||index as written
    n_blocks    u64
    total_check u128  sum of all n mod 2^128
```

Structural rules (readers enforce them): flags == 0; the footer ends the file
exactly (footer_off + 88 == file size); n_blocks == ceil(n_records /
block_size); index entries strictly increase in first_n, byte_off and
rank_base, with entry 0 at (0, 0); decoded records have factors >= 2,
strictly ascending, with no overflow in the reconstruction.

## 4. Queries

The index (32 B per block of 4096 records) fits in memory. A binary search on
`first_n` finds the covering block, and at most 4096 records are decoded
sequentially. This gives membership, factorization, predecessor and
successor, range enumeration, rank (count of CNs <= n, via rank_base plus the
position in the block) and select.

## 5. Verification

Encoding (`oracle_encode`, streaming, one pass over the text):
1. SHA-256 of the source bytes while parsing.
2. Per line: n ascending, k in [3, 14], factors strictly ascending,
   n = product(pi), and Korselt's criterion lambda(n) = lcm(pi - 1) divides
   n - 1.
3. `--prove-primes`: deterministic Miller-Rabin (`src/mr64.hpp`) on every
   factor; all factors are below 10^12, so the test is a proof.
4. Header, body, index and footer are written with both SHA-256 values.

Decoding (`oracle_decode`, the certificate): the full text is regenerated
from the ORC1 file alone, and its SHA-256 must equal sha256_src
("round-trip CERTIFIED"). `oracle_decode carmichael-1e24-oracle.orc1 -` does
this without writing the text (`make verify-oracle`).

Every other file in this repository is certified against the oracle, never
against the text.
