#!/usr/bin/env bash
# cells/111x-pwt-dstar/build.sh -- rebuild the PWT + d* file and require the published bytes.
#
# The d* divisors come from a 33-hour local search (dsel_set2). The file
# stores the chosen divisors, so it can be rebuilt without repeating the
# search. ROUTE selects how (default: paths).
#
#   ROUTE=paths   from the hosted d* file itself: pwt2_decode --emit-paths
#                 exports the stored divisors' paths (14 fields, tree order),
#                 and pwt2_encode re-encodes them to the same bytes. Needs the
#                 d* file and the oracle.
#   ROUTE=search  rerun the search with the logged command lines (6 threads,
#                 4 sweeps, then one more sweep from the saved masks; about
#                 33 h and 16 GB RAM), then encode. NOT rerun: the published
#                 search (2026-09-24/26) ran once, with an earlier build of
#                 dsel_set2 than src/ gives, and it prices candidates in
#                 floating point, so a rerun is not known to reproduce the
#                 masks that search saved byte for byte (masks_sha256 in CELL).
#   ROUTE=masks   from masks the search writes (--masks-out, one mask per CN,
#                 docs/MASKS_FORMAT.md; default $WORK/dstar_masks5.bin, the
#                 output of ROUTE=search): dsel_set2 --sweeps 0 skips the
#                 search and emits the paths (their SHA-256 must equal the
#                 dstar_paths intermediate), then pwt2_encode. Needs the
#                 oracle, the masks and the d_min paths. The masks of the
#                 published search are not hosted.
#
# Env: BIN, DATA, WORK, ORC1, THREADS, OUT (scripts/common.sh); DSTAR (the
# d* file for ROUTE=paths), MASKS (the masks for ROUTE=masks).
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/../.." && pwd)
. "$ROOT/scripts/common.sh"
CELL=$here/CELL
M=$ROOT/results/INTERMEDIATES
need_oracle
mkdir -p "$WORK/scratch"
OUT=${OUT:-$WORK/$(mget "$CELL" file)}
ROUTE=${ROUTE:-paths}
# masks route: the settings of the 2026-10-06 rebuild (batch, huge, cb and
# slack spelled out at their defaults)
DS="--threads 6 --batch 65536 --max-overflow 60000000 --huge 128 --cb 128 --slack 16 --min-gain 0.01 --check-rebuild --chunk 33554432"
# search route: the flags of the logged search runs, verbatim
DS_LOGGED="--threads 6 --min-gain 0.01 --check-rebuild --max-overflow 60000000 --chunk 33554432"

case $ROUTE in
    paths)
        DSTAR=${DSTAR:-$DATA/$(mget "$CELL" file)}
        [ -f "$DSTAR" ] || die "d* file not found: $DSTAR (make fetch-111x, or set DSTAR=...)"
        [ "$DSTAR" != "$OUT" ] || die "OUT must differ from the input d* file"
        P=$WORK/c111x.paths.bin
        run "$(tool pwt2_decode)" "$DSTAR" --emit-paths "$P" --fields 14 || exit 1 ;;
    masks)
        MASKS=${MASKS:-$WORK/dstar_masks5.bin}
        [ -f "$MASKS" ] || die "masks file not found: $MASKS (ROUTE=search writes it; or set MASKS=...)"
        check_file "$MASKS" "$(mget "$CELL" masks_bytes)" "$(mget "$CELL" masks_sha256)" "d* masks" || exit 1
        # --dpaths supplies d_min (dsel_set2 re-derives a sample and compares).
        bash "$ROOT/scripts/intermediates.sh" dmin_paths || exit 1
        P=$WORK/$(mget "$M" dstar_paths.file)
        # shellcheck disable=SC2086
        run "$(tool dsel_set2)" --orc1 "$ORC1" --dpaths "$WORK/dpaths.bin" --sweeps 0 $DS \
            --scratch "$WORK/scratch" --masks-in "$MASKS" --emit "$P" || exit 1
        check_file "$P" "$(mget "$M" dstar_paths.bytes)" "$(mget "$M" dstar_paths.sha256)" "d* paths" || exit 1 ;;
    search)
        bash "$ROOT/scripts/intermediates.sh" dmin_paths || exit 1
        P=$WORK/$(mget "$M" dstar_paths.file)
        # The logged first run also wrote the 4-sweep paths (--emit); they are
        # not needed here and are left out (the masks do not depend on it).
        # shellcheck disable=SC2086
        run "$(tool dsel_set2)" --orc1 "$ORC1" --dpaths "$WORK/dpaths.bin" --sweeps 4 $DS_LOGGED \
            --scratch "$WORK/scratch" --masks-out "$WORK/dstar_masks4.bin" || exit 1
        # shellcheck disable=SC2086
        run "$(tool dsel_set2)" --orc1 "$ORC1" --dpaths "$WORK/dpaths.bin" --sweeps 1 $DS_LOGGED \
            --scratch "$WORK/scratch" --masks-in "$WORK/dstar_masks4.bin" \
            --masks-out "$WORK/dstar_masks5.bin" --emit "$P" || exit 1
        check_file "$WORK/dstar_masks5.bin" "$(mget "$CELL" masks_bytes)" "$(mget "$CELL" masks_sha256)" "d* masks" || exit 1
        check_file "$P" "$(mget "$M" dstar_paths.bytes)" "$(mget "$M" dstar_paths.sha256)" "d* paths" || exit 1 ;;
    *) die "ROUTE must be paths, masks or search" ;;
esac

# encode (the d* paths use 14 fields, as the published input did)
run "$(tool pwt2_encode)" "$P" "$ORC1" "$OUT" --fields 14 --ctx-set 3 \
    --scratch "$WORK/scratch" --chunk 33554432 --threads "${ENC_THREADS:-4}" || exit 1
check_file "$OUT" "$(mget "$CELL" bytes)" "$(mget "$CELL" sha256)" "cell 111x ($OUT)"
