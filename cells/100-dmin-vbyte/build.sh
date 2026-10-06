#!/usr/bin/env bash
# cells/100-dmin-vbyte/build.sh -- rebuild cell 100
# Rebuild cell 100 from the oracle and require the published bytes.
# Route: ORC1 --dlist_from_paths --orc1 (d_min search)--> d.raw --list_build vbyte--> VBG1 file
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

# 1. extraction (ORC1 -> d.raw)
bash "$ROOT/scripts/intermediates.sh" dmin_list || exit 1

# 2. encode
run "$(tool list_build)" vbyte "$WORK/d.raw" "$OUT" || exit 1
check_file "$OUT" "$(mget "$CELL" bytes)" "$(mget "$CELL" sha256)" "cell $tag ($OUT)"
