#!/usr/bin/env bash
# oracle/decode.sh [verify|certify|decode] -- check or decode the oracle (ORC1).
#   verify   size + SHA-256 of the oracle, then oracle_decode regenerates the
#            text in memory (nothing written) and checks the round trip: the
#            regenerated text's SHA-256 must equal the sha_src recorded in the
#            oracle, which is the text's (oracle/TEXT)
#   certify  the same: the oracle is the reference every cell is certified
#            against, and its own reference is the text
#   decode   write the text to $WORK/new_table.txt (the full table) and hash it
# Env: BIN, DATA, WORK, ORC1 (scripts/common.sh).
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/.." && pwd)
. "$ROOT/scripts/common.sh"
O=$here/ORACLE
T=$here/TEXT
mode=${1:-verify}
case $mode in
    verify|certify)
        need_oracle
        check_file "$ORC1" "$(mget "$O" bytes)" "$(mget "$O" sha256)" "oracle" || exit 1
        run "$(tool oracle_decode)" "$ORC1" - > /dev/null || exit 1
        echo "VERIFIED  oracle: file hash and round trip to the text's SHA-256 $(mget "$T" sha256)" ;;
    decode)
        need_oracle
        mkdir -p "$WORK"
        run "$(tool oracle_decode)" "$ORC1" "$WORK/$(mget "$T" file)" || exit 1
        check_file "$WORK/$(mget "$T" file)" "$(mget "$T" bytes)" "$(mget "$T" sha256)" "decoded text" ;;
    *) die "usage: decode.sh verify|certify|decode" ;;
esac
