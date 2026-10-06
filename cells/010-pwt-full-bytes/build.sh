#!/usr/bin/env bash
# cells/010-pwt-full-bytes/build.sh -- rebuild cell 010
# Rebuild cell 010 from the oracle and require the published bytes.
# Route: ORC1 --trie_paths --full--> fpaths.bin --pwt2_encode --full --naive--> PWT2 byte-aligned file
# Usage: build.sh   (env: BIN, DATA, WORK, ORC1, THREADS, OUT; see scripts/common.sh)
# The extraction step builds a shared intermediate in $WORK (kept for the other
# cells); the encoder writes $OUT (default $WORK/<file>). Outputs do not depend
# on THREADS or ENC_THREADS.
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/../.." && pwd)
. "$ROOT/scripts/common.sh"
CELL=$here/CELL
tag=$(mget "$CELL" tag)
need_oracle
mkdir -p "$WORK/scratch"
OUT=${OUT:-$WORK/$(mget "$CELL" file)}

# 1. extraction (ORC1 -> fpaths.bin)
bash "$ROOT/scripts/intermediates.sh" full_paths || exit 1

# 2. encode
run "$(tool pwt2_encode)" "$WORK/fpaths.bin" "$ORC1" "$OUT" --fields 14 --full --naive \
    --scratch "$WORK/scratch" --chunk 33554432 --threads "${ENC_THREADS:-4}" || exit 1
check_file "$OUT" "$(mget "$CELL" bytes)" "$(mget "$CELL" sha256)" "cell $tag ($OUT)"
