# PWT1: the Prime Wavelet Tree container

Version 1. PWT1 is the value-ordered prime wavelet tree (largest prime
nearest the root), the predecessor of the frequency-ordered PWT2
(docs/PWT2_FORMAT.md) used by every tree file in this repository. PWT2
reuses its rANS symbol blocks, raw bit stream and table format (section 4),
which is why this specification is included. Implementation: src/pwt1.hpp
(format, encoder, reader), src/rans.hpp (entropy coder),
src/prime_bitmap.hpp (prime table rebuilt at load), src/pwt_paths.hpp
(sorted path input); tests test_pwt (in `make test`). The per-depth
universe flags let the same format hold the tree over d_min paths and the
tree over full factorizations.

## 1. What is stored

A set of positive integers, each given by a PATH of distinct odd primes.
Two modes, chosen by a header flag:

- d_min mode: the path is the prime set of d_min(n), the smallest divisor
  d of n with d * lambda(d) > n. The decoder recovers
  n = d * (d^-1 mod lambda(d)), lambda(d) = lcm(p - 1 : p | d)
  (requires gcd(d, lambda(d)) = 1, true for every Carmichael number).
- full mode: the path is the complete factorization; n = product.

Tree: primes in DESCENDING order along a path (largest at the root's child);
a node is identified by the multiset of primes above it, so equal prefixes
are shared. Children of a node are ordered by ascending prime. A node is
TERMINAL if some stored path ends there; leaves are always terminal and
an internal node may be terminal too (its number divides no other stored
number's path, it is merely a prefix of one).

Invariants the encoder asserts on every path (data invariants, never
assumed): primes >= 3 and strictly descending; for consecutive primes
p > q on a path, p != 1 (mod q) (Korselt: q | n and p - 1 | n - 1 forbid
q | p - 1). The reader re-checks what it can (children below their parent,
primes below the sieve limit where claimed).

## 2. Layout

    [header 128 B] [tables] [symbol blocks] [raw bits] [footer 88 B]

All integers little-endian. Offsets in the header must tile the file
exactly (tables at 128, then symbols, then raw bits, then the footer at
footer_off = file size - 88); a reader rejects anything else.

### Header (128 bytes)

    off  size  field
      0    4   magic "PWT1"
      4    2   version = 1
      6    2   flags: bit 0 = full mode (0 = d_min), bit 1 = descending (always 1)
      8    8   n_records
     16    1   max_depth K (1..30; 12 for the d_min tree, 14 for the full tree)
     17    3   reserved (0)
     20    4   depth_low: bit d (1 <= d <= K) set iff EVERY prime at depth d
               is below sieve_limit
     24    8   sieve_limit B (3 <= B <= 2^34; 100,000,000 for both tables)
     32    8   tables_off (= 128)      40  8  tables_len
     48    8   sym_off                 56  8  sym_len
     64    8   raw_off                 72  8  raw_len
     80    8   n_blocks
     88    4   block_symbols (symbols per block, 2^20; the last block is shorter)
     92    4   reserved (0)
     96    8   total_symbols
    104    8   footer_off
    112   16   reserved (0)

### Footer (88 bytes, the CND1 convention)

    off  size  field
      0   32   sha_payload: SHA-256 of bytes [0, footer_off)
     32   32   sha_nset: SHA-256 of the sorted n-set as 16-byte LE records
     64   16   total_check: sum of all n mod 2^128
     80    8   record_count (= n_records)

sha_nset, total_check and record_count are taken from the certified ORC1
oracle at encode time. A decoder verifies them after decoding and sorting;
the elementwise comparison against the oracle stream is the final gate.

## 3. Contexts and models

Every coded symbol belongs to a CONTEXT with its own static frequency
table (order-0 model). With K = max_depth, bits(p) = bit length of p
capped at 41 (42 values), bucket(c) in {0: c = 1, 1: c in 2..3,
2: c in 4..15, 3: c >= 16}, the contexts are enumerated:

    count ctx   id = depth * 42 + bits(prime(v))          depth 0..K     alphabet 4096
    flag ctx    id = (K + 1) * 42 + (depth - 1)            depth 1..K     alphabet 2
    gap ctx     id = (K + 1) * 42 + K
                   + ((((child_depth - 1) * 4 + bucket(c)) * 42 + bits(prime(parent))) * 2 + U)
                                                            child_depth 1..K, U in {0 low, 1 high}
                                                                           alphabet 128
    total       (K + 1) * 42 + K + K * 4 * 42 * 2   (= 4,590 for K = 12; 5,348 for K = 14)

The root has prime 0 (bits 0) and depth 0. Gap-length symbols are bounded
by the data: a prime below 2^40 gives an odd-integer gap of at most 39
bits and a prime index below 2^34 gives an index gap of at most 34 bits;
the reader rejects any length symbol above 63.

Tables section: for each context id in order, VByte(m) then m pairs
(VByte(symbol - previous symbol), VByte(freq - 1)) over the USED symbols
in ascending order; frequencies are quantized to 16 bits (sum = 65536,
every used symbol >= 1, largest-remainder rounding; freq - 1 must be
below 65536). An unused context is the single byte 0. VByte here and in
the block headers is the project's base-128 code: 7 payload bits per byte,
least significant group first, 0x80 = continuation. The full-table tables
are 64-68 KB.

## 4. Symbol stream (pre-order walk)

Starting at the root (depth 0, prime 0), for a node v with child count c:

1. count symbol in count ctx(depth(v), prime(v)): min(c, 4095). If
   c >= 4095 the raw stream carries Elias gamma of (c - 4094): that many
   zero bits, then the value's bits from its leading 1.
2. if c > 0 and depth(v) >= 1: flag symbol in flag ctx(depth(v)):
   1 iff v is terminal. (Leaves are terminal without a flag; the root is
   never terminal.)
3. for each child u of v in ascending prime order:
   - universe U: LOW (prime-index) iff depth(v) >= 1 and prime(v) < B, or
     depth_low has bit depth(u) set; else HIGH (odd-integer).
   - gap g >= 1:
       LOW : first child  g = index(prime(u)) + 1
             later child  g = index(prime(u)) - index(prime(previous sibling))
             where index(p) = number of odd primes below p (3 -> 0)
       HIGH: first child  g = (prime(u) - 1) / 2
             later child  g = (prime(u) - prime(previous sibling)) / 2
   - gap-length symbol L = bit length of g in gap ctx(depth(u), c, prime(v), U);
     the raw stream carries the L - 1 bits of g below its leading 1.
   - then u's own count, flag and children (recursion).

Why the universe rule is decodable: prime(v) and depth_low are known
before v's children are read. Why it is valid: in a descending path every
child is smaller than its parent, so a parent below B has all children
below B; a depth flagged low has, by the encoder's scan, no prime at or
above B. On the d_min tree every depth >= 2 is below 1e8; on the full
tree depth 2 reaches 3.8e9, so depth-2
children use the prime-index universe only under parents below B, and
depths >= 3 (max 9.5e7) are flagged low.

Symbol blocks: the symbol sequence is cut into blocks of block_symbols;
each block is stored as u32 symbol count, u32 byte length, then the rANS
bytes. rANS (rans.hpp, the "rans_byte" construction): 32-bit state x in
[2^23, 2^31), 16-bit probabilities (freq f, cumulative cum per symbol
within its context). The block's bytes begin with the initial state as 4
bytes, most significant first, followed by the renormalization bytes in
decode order. Decode step for a symbol in context C:
slot = x mod 2^16; find the symbol s of C with cum(s) <= slot < cum(s) + f(s);
x = f(s) * floor(x / 2^16) + slot - cum(s); while x < 2^23: x = 256 x + next byte.
The encoder produces this by coding the block's symbols in REVERSE order
(x_max = (2^23 / 2^16) * 256 * f; while x >= x_max emit the low byte and
shift; x = floor(x / f) * 2^16 + x mod f + cum). A block must end with the
state back at exactly 2^23 and no byte left; the reader rejects otherwise.

Raw bits: one MSB-first bit stream (bic.hpp BitWriter), zero-padded to a
byte at the very end; a reader checks the padding is zero and that the
stream's length matches.

## 5. Reader obligations

A conforming reader rejects the file unless: magic, version and flags are
as above; the sections tile [0, footer_off) and the footer ends the file;
sha_payload matches; every table deserializes (ascending symbols within
the alphabet, frequencies summing to 65536); the block directory matches
n_blocks and total_symbols and every block but the last holds
block_symbols symbols; during the walk every rANS block ends cleanly,
every child prime is >= 3 and below its parent, every gap-length symbol
is at most 63, every prime-index is within the sieve, the escape
terminates, the decoded record count equals
n_records, the symbol count equals total_symbols, and only zero padding
remains in the raw stream. test_pwt flips every kind of byte (400 random
single-bit flips per container, truncations, an extension, and 200 flips
below the payload hash with the hash recomputed) and requires each to be
rejected or to fail the footer targets; none may decode to the right set
silently.

## 6. PWT0: the byte-aligned variant (the C = 0 tree cells)

Same tree, same pre-order walk, same header (magic "PWT0", flags bit 2
set, tables_len = raw_len = n_blocks = block_symbols = 0) and footer, but
the symbol section is a plain VByte stream and nothing is entropy coded:

- per node one VByte holding 2c + terminal (leaves therefore carry 1);
- per child one VByte holding the plain integer gap: the prime itself for
  a first child, prime minus the previous sibling's prime otherwise. No
  universes, no prime table, no differencing tricks beyond that sibling
  delta.

This is exactly the "every symbol one VByte" serialization (test_pwt
checks the stream length, 4,936 and 8,574 bytes, on the fixture). The
reader applies the same structural checks (child odd and below its parent,
record and symbol counts, stream consumed exactly, payload SHA-256, footer
targets); pwt1.hpp writes and reads either variant.

## 7. Decoding cost

Prime table: a bitmap of the odd numbers below B (B/16 bytes, 6.25 MB at
1e8) with a rank directory, built by a sieve at load. select() (first
child) is a binary search over the directory plus popcounts; later
siblings advance from the previous sibling by a popcount scan. Per stored
number the decoder then does one product (full mode) or one product, one
lcm chain and one modular inverse (d_min mode). No factoring anywhere.
