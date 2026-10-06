# Tests

| what | run with | needs |
|---|---|---|
| `unit/` C++ unit tests of the headers and containers (sieve, oracle, splitter, CND, PWT1, PWT2, Algorithm 1 factoring), including corruption and crafted-file batteries | `make test` | `fixtures/` |
| `prefix/` the self-test: every tool of the pipeline on the Carmichael numbers below 10^12, each output checked against `prefix/golden_prefix.tsv` (frozen sizes and SHA-256), every decoder certified against the prefix oracle | `make check` | the tools in `build/` |
| `test_release_tools.sh` fixture checks of `list_decode`, `pwt2_decode --emit-paths` (encode, export, re-encode byte for byte) and `paths_cmp`, with corruption cases | `make test-release-tools` | `fixtures/`, perl |

The folder `fixtures/` holds these files.
- `fix_head.txt`, `fix_mid.txt` and `fix_tail.txt` are lines from the head, middle and tail of the table. They hold 637 Carmichael numbers with their factorizations, the ground truth of the unit tests.
- `pwt2_craft/` holds crafted PWT2 containers. Some were written by an independent encoder and some are mutations. `MANIFEST.txt` gives the verdict a reader must reach for each one.

The folder `prefix/` holds these files.
- `prefix_1e12.txt` is the first part of the table, every Carmichael number below 10^12.
- `prefix.orc1` is its oracle.
- `golden_files/` holds the nine cell files of the prefix.
- `golden_prefix.tsv` gives the frozen size and SHA-256 of every file the self-test builds. A mismatch means a tool writes different bytes.
- `selftest_runs.txt` and `runs/` record the runs that froze and confirmed the goldens.
- `SHA256SUMS` lists the hashes of the files in `prefix/`.

The table-derived files here (`fixtures/fix_*.txt`, `prefix/prefix_1e12.txt`, `prefix/prefix.orc1`, `prefix/golden_files/`) are excerpts of the Shallue-Webster table. They are data, dedicated to the public domain under CC0 1.0 like the hosted files. See the README's Data section.
