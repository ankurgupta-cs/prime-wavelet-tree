#!/usr/bin/env bash
# cells/111x-pwt-dstar/decode.sh -- verify, certify or decode cell 111x
# decode.sh [verify|certify|decode|paths] -- check and decode cell 111x's file
# (default: $DATA/<file>, else $WORK/<file>; override with FILE=...).
#   verify   size + SHA-256, then pwt2_decode: payload SHA-256, tables, every
#            number rebuilt from its path, sorted, and the footer's sha_nset,
#            total_check and count must match (the footer pins the CN set)
#   certify  pwt2_decode --certify: every decoded CN compared with the oracle
#   decode   write the CNs (16-byte LE, tree order) to $WORK/c111x.n.raw
#   paths    export the stored paths (14 fields, tree order) to $WORK/c111x.paths.bin;
#            they re-encode to the same file with the build's pwt2_encode flags
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
        run "$(tool pwt2_decode)" "$F" || exit 1
        echo "VERIFIED  cell $tag: file hash, structure and footer targets" ;;
    certify)
        need_oracle
        run "$(tool pwt2_decode)" "$F" --certify "$ORC1" ;;
    decode)
        mkdir -p "$WORK"
        run "$(tool pwt2_decode)" "$F" --emit-n "$WORK/c$tag.n.raw" ;;
    paths)
        mkdir -p "$WORK"
        run "$(tool pwt2_decode)" "$F" --emit-paths "$WORK/c$tag.paths.bin" --fields 14 ;;
    *) die "usage: decode.sh verify|certify|decode|paths" ;;
esac
