#!/usr/bin/env bash
# oracle/build.sh -- build the oracle (ORC1) from the Shallue-Webster text table.
#
#   text (oracle/TEXT) --oracle_encode --prove-primes--> carmichael-1e24-oracle.orc1
#
# Every prime factor is proven prime while encoding (deterministic
# Miller-Rabin; every factor is below 10^12). The result must have the
# published size and SHA-256 (oracle/ORACLE); then oracle_decode regenerates
# the text in memory and checks its SHA-256 against the text's (round trip).
# Env: BIN, DATA, WORK (scripts/common.sh); TEXT_FILE (default $DATA/new_table.txt,
# fetched with `make fetch-text`), OUT (default $WORK/carmichael-1e24-oracle.orc1).
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/.." && pwd)
. "$ROOT/scripts/common.sh"
T=$here/TEXT
O=$here/ORACLE
TEXT_FILE=${TEXT_FILE:-$DATA/$(mget "$T" file)}
OUT=${OUT:-$WORK/$(mget "$O" file)}
[ -f "$TEXT_FILE" ] || die "text table not found: $TEXT_FILE (make fetch-text, or set TEXT_FILE=...)"
mkdir -p "$(dirname "$OUT")"
check_file "$TEXT_FILE" "$(mget "$T" bytes)" "$(mget "$T" sha256)" "text table" || exit 1
run "$(tool oracle_encode)" "$TEXT_FILE" "$OUT" --prove-primes || exit 1
check_file "$OUT" "$(mget "$O" bytes)" "$(mget "$O" sha256)" "oracle ($OUT)" || exit 1
run "$(tool oracle_decode)" "$OUT" - > /dev/null
