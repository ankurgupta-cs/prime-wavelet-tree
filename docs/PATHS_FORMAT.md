# Path files: the encoder input of the prime wavelet tree

A path file holds one record per Carmichael number: the primes of the
stored divisor (or of the full factorization), which `pwt2_encode` turns
into a tree. Code: `src/pwt_paths.hpp` (record layout, external sort).

## Record layout

    record := K fields of 5 bytes; field = one prime, BIG-endian (primes < 2^40)
              primes strictly ascending, unused fields zero (at the end)

- K = 12 for the d_min paths (`dpaths.bin`; a d_min has at most 12 primes),
  K = 14 for the full factorizations (`fpaths.bin`; a CN has at most 14) and
  for the d* paths (`dstar_paths.bin`). Pass K to the tools with `--fields K`.
- Big-endian fields make `memcmp` order equal the lexicographic numeric
  order, with a prefix sorting before its extensions.
- There is no header and no checksum; the file size is N x 5K bytes and the
  files are identified by their SHA-256 (`results/INTERMEDIATES`).

## Writers and order

| writer | content | record order |
|---|---|---|
| `trie_paths ORC1 out --fields 12` | primes of d_min(n) | n order (the oracle's) |
| `trie_paths ORC1 out --full --fields 14` | every prime of n | n order |
| `dsel_set2 ... --emit out` | primes of the chosen divisor (d* or d_min), always 14 fields | n order |
| `pwt2_decode file --emit-paths out --fields K` | the stored paths of a PWT2 file | tree order |

`pwt2_encode` accepts any record order (it sorts), so a file exported from a
tree re-encodes to the same container. Two path files hold the same set of
paths exactly when `paths_cmp a b --fields K` prints SAME SET; their raw
SHA-256 differ when the record orders differ.

## Validity

A divisor path is valid for n when d = product of its primes satisfies
d * lambda(d) > n, lambda(d) = lcm(p - 1); then n = d * (d^{-1} mod
lambda(d)). A full path is valid when the product is n. The writers check
this per record; `pwt2_encode` checks the sum of the rebuilt n against the
oracle's total_check, and `pwt2_decode --certify` checks every n.
