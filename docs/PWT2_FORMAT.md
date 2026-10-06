# PWT2: the frequency-ordered Prime Wavelet Tree container

Version 1. Implementation: src/pwt2.hpp (format, encoder, reader),
src/rans.hpp (entropy coder), src/prime_bitmap.hpp (sieve rebuilt at load);
programs pwt2_encode, pwt2_decode (`--certify`, `--emit-n`, `--emit-paths`);
tests test_pwt2 (in `make test`). PWT1 (docs/PWT_FORMAT.md), the
value-ordered tree, is a separate magic; PWT2 reuses its symbol blocks.
Cells 010, 011, 110, 111 and the PWT + d* file are PWT2 files.

## 1. What is stored

A set of positive integers, each given by a PATH: a set of distinct odd
primes. Two modes (header flag):

- divisor mode: the path is the prime set of a divisor d of n with
  d * lambda(d) > n (d_min(n) in cells 110 and 111; any such divisor is
  allowed, as in the PWT + d* file). The decoder recovers
  n = d * (d^-1 mod lambda(d)), lambda(d) = lcm(p - 1 : p | d).
- full mode: the path is the complete factorization; n = product.

## 2. The order and the tree

Every prime that occurs on some path gets a path count c(p) (the number of
stored paths containing p). Ranks 1..m order the primes by c DESCENDING,
ties to the LARGER prime. A path is written as its ranks in ascending order
(most frequent prime first); the tree is the trie of those rank sequences,
children in ascending rank. A node is TERMINAL if a stored path ends there;
leaves are terminal and an internal node may be terminal too.

Why this order: a node covers every number below it, so putting the most
frequent primes nearest the root makes each node cover as many numbers as
possible (the frequency-shaped wavelet tree of Grossi, Gupta and Vitter,
SODA 2004; the FP-tree order of Han, Pei and Yin, 2000). On the full table
the d_min tree has 1.63 nodes per number instead of 2.74 in the descending
PWT1 tree. In the full-factorization tree every multiple of the most
frequent prime (19) lies in ONE subtree, and the multiples of the k-th most
frequent prime in at most 2^(k-1) (in the d_min tree the 19-subtree holds
the 26.3M numbers whose d_min contains 19; the other multiples carry 19 in
the cofactor).

The runs of equal count are the CLASSES: consecutive ranks with the same
c(p), each class in strictly descending prime order. The decoder never
needs the counts, only the class sizes.

## 3. Layout

    [header 160 B] [class sizes] [tables] [symbol blocks] [raw bits] [footer 88 B]

All integers little-endian. The sections tile [0, footer_off) exactly and
the footer ends the file; a reader rejects anything else.

### Header (160 bytes)

    off  size  field
      0    4   magic "PWT2"
      4    2   version = 1
      6    2   flags: bit 0 = full mode (0 = divisor mode), bit 2 = byte-aligned
               variant, bit 3 = frequency order (always set); other bits 0
      8    8   n_records
     16    1   max_depth K = the length of the longest path (1..30; the
               reader checks that some path reaches it)
     17    1   context set: 1 (sections 4-5), 2 (section 9) or 3 (section 10);
               the byte-aligned variant always uses 1
     18    2   reserved (0)
     20    4   wheel modulus of the table's high universe (210)
     24    8   sieve limit B (211 <= B <= 2^34; 100,000,000 for the tables)
     32    8   m = number of distinct primes (the rank universe), 1 <= m < 2^32
     40    8   n_classes (1..m)
     48    8   classes_off (= 160)     56  8  classes_len
     64    8   tables_off              72  8  tables_len
     80    8   sym_off                 88  8  sym_len
     96    8   raw_off                104  8  raw_len
    112    8   n_blocks
    120    4   block_symbols: exactly 2^20 (entropy-coded), 0 (byte-aligned)
    124    4   reserved (0)
    128    8   total_symbols
    136    8   table_symbols = m (the first m symbols are the rank table;
               redundant, kept so a reader can check the split)
    144    8   footer_off
    152    8   reserved (0)

### Footer (88 bytes, the CND1 / PWT1 convention)

    off  size  field
      0   32   sha_payload: SHA-256 of bytes [0, footer_off)
     32   32   sha_nset: SHA-256 of the sorted n-set as 16-byte LE records
     64   16   total_check: sum of all n mod 2^128
     80    8   record_count (= n_records)

sha_nset, total_check and record_count come from the certified ORC1 oracle
at encode time; pwt2_decode --certify verifies them after decoding and
sorting, then compares element by element with the oracle stream.

### Class sizes

n_classes VBytes, each >= 1, summing to m, filling the section exactly.
VBytes here and in section 6 are the project's base-128 code (7 payload
bits per byte, least significant group first, 0x80 = continuation), with
two strictness rules: the value must be below 2^64 (a tenth byte may only
be 0x00 or 0x01), and the encoding must be canonical (a multi-byte VByte
does not end in 0x00). The set of stored numbers is never empty.

## 4. Contexts and models (context set 1)

Every coded symbol has a CONTEXT with its own static quantized frequency
table (order-0 model per context). bits(v) = bit length of v; cap(v, h) =
min(bits(v), h). bucket(c) in {0: c = 1, 1: c in 2..3, 2: c in 4..15,
3: c >= 16}. bits(0) = 0 (the root's rank). vidx(r) = 1-based position of
prime(rank r) in the ascending list of the m used primes; the root counts
as vidx = m + 1.

    table ctx  id = cap(class size, 40) * 2 + H          H = 1 iff the previous member is at or above B   alphabet 64
    count ctx  id = 82 + depth * 42 + cap(vidx(node), 41)           depth 0..K                    alphabet 4096
    flag ctx   id = 82 + (K + 1) * 42 + (depth - 1)                 depth 1..K                    alphabet 2
    gap ctx    id = 82 + (K + 1) * 42 + K
                  + ((child_depth - 1) * 4 + bucket(c)) * 42 + cap(rank(parent), 41)
                                                                    child_depth 1..K              alphabet 64
    total      82 + (K + 1) * 42 + K + K * 4 * 42     (2,656 for K = 12; 3,078 for K = 14)

Tables section: for each context id in order, VByte(number of used
symbols) then (VByte(symbol - previous symbol), VByte(freq - 1)) pairs over
the used symbols ascending; frequencies quantized to 16 bits (sum 65,536,
every used symbol >= 1). An unused context is the single byte 0. Absent in
the byte-aligned variant.

Quantization (needed only for byte-identical encoders; a decoder accepts
any table that sums to 65,536). The reference encoder (rans.hpp
RansContext::build) counts each context's symbols over the whole file,
then for every used symbol s with count c_s out of total T: exact_s =
c_s / T * 65536 in IEEE double, f_s = max(1, floor(exact_s)); the pairs
(exact_s - f_s, index of s among the used symbols) are sorted ascending
(std::sort on pairs, so equal remainders order by index); if the f_s sum
to less than 65,536 the deficit is paid by adding 1 to the entries from
the END of that sorted list backwards (largest remainder first, ties to
the higher index); if they sum to more, 1 is taken repeatedly from the
first maximal f_s (lowest index among the largest) until the sum is
65,536.

## 5. Symbol stream

The stream is the rank table followed by the tree.

### Rank table

For each class in rank order (sizes from the class-size section), its
primes in ASCENDING order (the reverse of their rank order) as gaps in the
table universe:

    pos(p) = index(p)                         if p < B   (odd primes below p; 3 -> 0)
           = cnt(B) + w(p) - w(B')            if p >= B  (w = wheel-210 index, B' = least integer >= B coprime to 210)
    w(x)   = 48 * floor(x / 210) + (position of x mod 210 among the 48 residues coprime to 210)
    first member: g = pos + 1; later members: g = pos - pos(previous member)

Each gap is a length symbol L = bits(g) in table ctx(class size, previous
member at or above B; false for the first member) followed by the L - 1 bits of
g below its leading 1 in the raw stream. The decoder maps pos back
(select in the prime bitmap below cnt(B), the wheel above) and stores the
class in reverse (descending) order at its ranks.

The table is valid only if every class is strictly descending in value,
all m primes are distinct, odd, >= 3 and below 2^40 (both variants), and
every value below B is a prime of the sieve. Values at or above B are not
rechecked for primality by the reader (a composite there can only decode
a wrong n, which the footer targets and the oracle comparison reject).

### Tree (pre-order walk)

Starting at the root (depth 0, rank 0), for a node v with child count c:

1. count symbol min(c, 4095) in count ctx(depth(v), vidx(v)); if c >= 4095
   the raw stream carries the Elias gamma code of v = c - 4094 >= 1:
   bits(v) - 1 zero bits, then the bits(v) bits of v from its leading 1
   (at most 40 zero bits are accepted).
2. if c > 0 and depth(v) >= 1: flag symbol (1 = terminal) in flag ctx.
3. for each child u of v in ascending rank: gap g = rank(u) - rank(v) for
   the first child, rank(u) - rank(previous sibling) for the others
   (g >= 1, rank(u) <= m); length symbol L = bits(g) in gap
   ctx(depth(u), c, rank(v)), then the L - 1 raw bits; then u's own
   count, flag and children.

A terminal node at depth t yields n from the primes of the ranks on its
path (any order): the product (full mode) or d * (d^-1 mod lambda(d))
(divisor mode). The file is invalid if a product or n reaches 2^128, if
gcd(d, lambda(d)) != 1 (no inverse), if lambda(d) >= 2^127, or if two
terminals yield the same n (the stored numbers are a set).

Symbol blocks: exactly as PWT1 (docs/PWT_FORMAT.md section 4): blocks of
block_symbols symbols, each stored as u32 symbol count, u32 byte length,
then the rANS bytes (32-bit state, 16-bit probabilities, the rans_byte
construction); every block but the last is full, there is no empty block,
so n_blocks = ceil(total_symbols / 2^20); a block must end with the state
back at 2^23 and no byte left. Table and tree symbols share the blocks.

Raw bits: one MSB-first bit stream, zero-padded to the next byte boundary
at the very end: raw_len = ceil(bits used / 8) and the padding bits are 0.

## 6. Byte-aligned variant (flag bit 2; the C = 0 cells 010 and 110)

Same order, classes and walk, no tables, no raw bits, no blocks
(tables_len = raw_len = n_blocks = block_symbols = 0). The symbol section
is a VByte stream:

- rank table: for each class, ascending: VByte(first prime), then
  VByte(prime - previous prime);
- tree: per node VByte(2c + terminal); per child VByte(rank gap) as in
  section 5.

total_symbols counts VBytes (so total_symbols <= sym_len); table_symbols =
m. The rank-table checks of section 5 apply unchanged (values below 2^40,
primes below B checked against the sieve); values at or above B need not
be coprime to 210 here, since no wheel is involved.

## 7. Reader obligations

A conforming reader rejects the file unless: magic, version, flags,
context set, wheel, block_symbols and reserved bytes are as above; m <
2^32, table_symbols = m, n_records <= total_symbols (and total_symbols <=
sym_len in the byte-aligned variant); the sections tile the file and the
footer ends it; sha_payload matches; the class sizes are canonical VBytes,
positive and sum to m; every table deserializes; the block directory
matches n_blocks and total_symbols; the rank table decodes to exactly m
symbols and passes section 5's checks; during the walk every child rank
lies in (rank(previous), m], every gap-length symbol is in 1..63, the
escape terminates, the depth never exceeds K and reaches it; every path's
n is valid (section 5); the decoded record count equals n_records, the
symbol count equals total_symbols, the last block ends cleanly, the
byte-aligned stream is consumed exactly, and only zero padding remains in
the raw stream. Context sets 2 and 3 add: every VByte of the tables section
canonical and below 2^64 (the model bytes included), keys strictly
ascending within a class and below 2^64 - 1, no empty model, symbols
within each class's alphabet, class-5 keys with L >= 2, a model present
for every key the walk uses (class 5: presence decides coded or raw), and
a count-tail symbol >= 63. Then (pwt2_decode): the decoded values, sorted, are
distinct and match sha_nset, total_check and record_count. The reader
never sizes an allocation by a header field alone. test_pwt2 flips every
header bit and 400 random payload bits (with and without a recomputed
payload hash), truncates, extends and crafts wrapping section lengths, and
runs the crafted containers in tests/fixtures/pwt2_craft/ (expected
verdicts in MANIFEST.txt); none may decode silently to a wrong set.

## 8. Decoding cost

Load: payload SHA-256, the prime bitmap to B (sieve, 6.25 MB at 1e8), the
context tables. Rank table: m gap decodes plus one select (below B) or
wheel step (above) each, then a sort of the m primes for vidx. Tree: one
count symbol, at most one flag and one gap per node; a child's prime is an
array lookup by rank. Per stored number one product (full mode) or one
product, one lcm chain and one modular inverse (divisor mode). No
factoring anywhere; the divisor-mode cofactor n / d is factored with
Shallue-Webster Algorithm 1 (src/cnfactor.hpp) when the whole factored
file is wanted.

## 9. Context set 2 (header byte 17 = 2; entropy-coded variant only)

The contexts were chosen by measurement on the set-1 files (variant
gapv 15 / cntv 3 / flagv 1 / tabv 3, without mantissa symbols). Same layout, order, rank
table, walk, blocks, raw stream and footer as set 1; only the symbol
classes, their contexts and the tables section differ.

Symbol classes and alphabets (index: name, alphabet):
0 table-gap length (64), 1 count (64), 2 count tail (4096), 3 flag (2),
4 child-gap length (64).

Counts: the count symbol is min(c, 63). If c >= 63 a tail symbol follows,
min(c, 4095) (a tail symbol below 63 makes the file invalid), in the tail
context of the node's depth; if c >= 4095 the raw stream carries the Elias
gamma code of c - 4094 (as in section 5). Compared with set 1 this adds one
tail symbol per node with 63 or more children (+656,490 symbols on the
d_min table, +568,657 on the full table).

Keys (all known before the symbol; bits() and cap as in section 4, AB =
bits(m) + 1, bits128(P) = bit length of a 128-bit value, P(v) = product of
the primes on the path to v inclusive, P(root) = 1):

    table  key = (cap(class size, 40) * 2 + H) * 42 + cap(pos(previous member), 41)
                 H as in set 1; the previous position is 0 for a class's first member
    count  k   = (d * AB + bits(rank(u))) * 82 + min(bits128(P(u)), 81)
           key = (k * 2 + F) * 14 + cap(A, 13)
                 u at depth d; F = 1 iff u is its parent's first child; A = siblings still
                 to come after u; the root has d = 0, rank 0, P = 1, F = 0, A = 0
    tail   key = d
    flag   key = (d - 1) * 4 + bucket(c(u))
    gap    k   = (((dc - 1) * 4 + bucket(c(v))) * 2 + F) * AB + bits(a)
           key = ((k * AB + bits((m - a) / R)) * 14 + cap(R, 13)) * 41 + floor(min(bits128(P(v)), 81) / 2)
                 the i-th child of v at child depth dc; F = 1 iff i = 0; a = rank of the
                 previous sibling (rank(v) for i = 0); R = c(v) - i children still to come
                 including this one; (m - a) / R is floor division

Tables section: for each class 0..4 in order, VByte(number of models),
then for each model in ascending key order VByte(key - previous key; the
first model's delta is its key), followed by the model in the section-4
table format (used symbols, frequencies summing to 65,536; an empty model
is invalid). Model (context) ids run 0, 1, 2, ... over the whole section in
this order and are 32-bit; a symbol whose key has no model makes the file
invalid. Every VByte of this section, the ones inside the models
included, follows the strictness rules of section 3; keys stay below
2^64 - 1. A model whose key the walk never uses is permitted (a reader
accepts it; the reference encoder writes only used models).

Decoding cost: one keyed lookup (hash) per symbol on top of set 1, and the
extra tail symbols.

## 10. Context set 3 (header byte 17 = 3): set 2 + tree-gap mantissa models

Set 3 models exact gap values given the previous rank: an order-1 model
of the child rank below the root and the 1,023 most frequent primes, and
bucketed contexts for rarer previous ranks. Everything of section 9 applies; the tables section
has a sixth class, and tree gaps may carry a mantissa symbol.

Class 5, tree-gap mantissa: for a tree gap g of bit length L >= 2 (table
gaps never have one), vb = L - 1 if L <= 12, else 4. Key:

    b   = a + 1024 * F                                   if a < 1024
        = 2048 + (((dc - 1) * 2 + F) * AB + bits(a)) * AB + bits((m - a) / R)   otherwise
    key = b * 64 + L

(a, F, R, dc as in the section-9 gap key). The model's alphabet is 2^vb.
If the tables hold a model for the key, the symbol stream carries
v = (g - 2^(L-1)) >> (L - 1 - vb) right after the gap-length symbol and
the raw stream carries only the low L - 1 - vb bits of g; otherwise the
raw stream carries all L - 1 bits as in set 2. The tables list, for class
5, only the keys whose model pays for itself. Reference encoder rule (an
encoder choice: a decoder only needs the presence or absence of each key):
build the quantized model (section 4) from the key's value counts; q = the
sum over its used symbols, in ascending symbol order, of count * (16 -
log2 freq) in IEEE double; S = the byte length of the model's own
serialization (the key-delta VByte is not included; the constant 3 below
stands in for it); keep the model iff vb * N - q - 8 * (S + 3) > 0, N =
the key's gap count. The model's alphabet is implied by L = key mod 64
(L >= 2 required).
