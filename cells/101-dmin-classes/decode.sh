#!/usr/bin/env bash
# cells/101-dmin-classes/decode.sh -- verify, certify or decode cell 101
# decode.sh [verify|certify|decode] -- check and decode cell 101's file
# (default: $DATA/<file>, else $WORK/<file>; override with FILE=...).
#   verify   size + SHA-256, then list_decode: structure, every value odd and
#            strictly ascending, and the SHA-256 of the sorted value stream
#            (16-byte LE) must equal the manifest's values_sha256
#   certify  cnd_decode --certify (or list_decode --kind d --certify with
#            CERTIFY_TOOL=list_decode): every d is factored, n = d (d^-1 mod
#            lambda(d)) rebuilt, and every n compared with the oracle, in order,
#            plus the footer targets (THREADS threads)
#   decode   write the sorted values to $WORK/c101.values.raw and the CNs to $WORK/c101.n.raw
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/../.." && pwd)
. "$ROOT/scripts/common.sh"
CELL=$here/CELL
tag=$(mget "$CELL" tag)
mode=${1:-verify}
F=$(cell_file "$CELL")
case $mode in
    verify)
        check_file "$F" "$(mget "$CELL" bytes)" "$(mget "$CELL" sha256)" "cell $tag file" || exit 1
        want=$(mget "$CELL" values_sha256)
        out=$(run "$(tool list_decode)" "$F") || { printf '%s\n' "$out"; exit 1; }
        printf '%s\n' "$out"
        got=$(printf '%s\n' "$out" | awk '$1 == "values" && $2 == "sha256" { print $3 }')
        if [ "$got" = "$want" ]; then echo "VERIFIED  cell $tag: file hash and sorted values hash $want"
        else echo "FAILED    cell $tag: values sha256 $got, expected $want" >&2; exit 1; fi ;;
    certify)
        need_oracle
        # cnd_decode is the CND1 decoder the file was certified with; list_decode
        # gives the same verdict (CERTIFY_TOOL=list_decode).
        if [ "${CERTIFY_TOOL:-cnd_decode}" = list_decode ]; then
            run "$(tool list_decode)" "$F" --kind d --certify "$ORC1" --threads "$THREADS"
        else
            run "$(tool cnd_decode)" "$F" --certify "$ORC1" --threads "$THREADS"
        fi ;;
    decode)
        mkdir -p "$WORK"
        run "$(tool list_decode)" "$F" --emit "$WORK/c101.values.raw" --kind d --emit-n "$WORK/c101.n.raw" --threads "$THREADS" ;;
    *) die "usage: decode.sh verify|certify|decode" ;;
esac
