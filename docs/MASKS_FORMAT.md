# d* masks: the divisor search checkpoint

`dsel_set2` (the d* search) chooses, for every Carmichael number n, the set S
of its primes whose product d is stored in the tree. Its checkpoint records
that choice as one mask per number. This is the format the search writes.

The masks of the published search, after its fifth sweep, are not hosted.
Their size and SHA-256 are `masks_bytes` and `masks_sha256` in
`cells/111x-pwt-dstar/CELL`. The d* file stores the chosen divisors, so it
can be rebuilt without repeating the search, from its own exported paths
(`cells/111x-pwt-dstar/build.sh`, `ROUTE=paths`, the default).

Code: `src/dsel_set2.cpp` (`--masks-out` writes, `--masks-in` reads).

## Layout

All integers little-endian.

    u64 N                  the record count (must equal the oracle's)
    N x u16 mask           in oracle (n) order

Bit j of a mask (j = 0 for the smallest prime) is set when the j-th prime of
n, in the ascending order of the oracle record, belongs to S. A mask is
nonzero and below 2^k for a CN with k primes, and the divisor must satisfy
d * lambda(d) > n. The file size is 8 + 2N bytes.

## Use

    dsel_set2 --orc1 ORC1 --dpaths dpaths.bin --sweeps 0 --masks-in MASKS --emit dstar_paths.bin ...

With `--sweeps 0` no search runs: the loaded masks are emitted as a path
file (`docs/PATHS_FORMAT.md`, 14 fields, n order), which `pwt2_encode
--fields 14 --ctx-set 3` turns into the d* file. With `--sweeps K` the
search resumes from the masks (the published run: four sweeps, then one more
from the saved masks). The exact flags are in
`cells/111x-pwt-dstar/build.sh`.
