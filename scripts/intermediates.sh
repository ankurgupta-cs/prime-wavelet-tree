#!/usr/bin/env bash
# intermediates.sh NAME... -- build the shared intermediate files from the oracle.
#
#   n_list       $WORK/n.raw       nstats_orc1 ORC1 --emit-n n.raw --emit-only
#   dmin_list    $WORK/d.raw       dlist_from_paths --orc1 ORC1 d.raw --threads T
#   dmin_paths   $WORK/dpaths.bin  trie_paths ORC1 dpaths.bin --fields 12 --threads T
#   full_paths   $WORK/fpaths.bin  trie_paths ORC1 fpaths.bin --full --fields 14 --threads T
#
# Sizes and SHA-256 come from results/INTERMEDIATES (rendered from
# results/cells.tsv). A file that is already present with the right size and
# hash is kept; otherwise it is written to <file>.partial, checked, then
# renamed. The outputs do not depend on the thread count (the prefix self-test
# checks this). n.raw and d.raw are 16-byte little-endian values, ascending;
# the path files hold 5-byte big-endian primes per field (docs/PATHS_FORMAT.md).
ROOT=$(cd "$(dirname "$0")/.." && pwd)
. "$ROOT/scripts/common.sh"
M=$ROOT/results/INTERMEDIATES
[ $# -ge 1 ] || die "usage: intermediates.sh n_list|dmin_list|dmin_paths|full_paths ..."
mkdir -p "$WORK"

build_one() {
    local name=$1 file bytes sha out part
    file=$(mget "$M" "$name.file"); bytes=$(mget "$M" "$name.bytes"); sha=$(mget "$M" "$name.sha256")
    [ -n "$file" ] && [ -n "$sha" ] || die "unknown intermediate $name"
    out=$WORK/$file
    if [ -f "$out" ] && [ "$(bytes_of "$out")" = "$bytes" ]; then
        if [ "$(sha256_of "$out")" = "$sha" ]; then
            printf 'present   %s (%s, sha256 %.16s...)\n' "$name" "$out" "$sha"; return 0
        fi
    fi
    need_oracle
    part=$out.partial
    case $name in
        n_list)     run "$(tool nstats_orc1)" "$ORC1" --emit-n "$part" --emit-only ;;
        dmin_list)  run "$(tool dlist_from_paths)" --orc1 "$ORC1" "$part" --threads "$THREADS" ;;
        dmin_paths) run "$(tool trie_paths)" "$ORC1" "$part" --fields 12 --threads "$THREADS" ;;
        full_paths) run "$(tool trie_paths)" "$ORC1" "$part" --full --fields 14 --threads "$THREADS" ;;
        *) die "no build rule for intermediate $name" ;;
    esac || { echo "FAILED    $name" >&2; return 1; }
    check_file "$part" "$bytes" "$sha" "$name" || return 1
    mv -f "$part" "$out"
}

rc=0
for n in "$@"; do build_one "$n" || rc=1; done
exit $rc
