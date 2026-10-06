#!/usr/bin/env bash
# selftest_prefix.sh -- the 10^12-prefix self-test behind `make check`.
#
# Rebuilds every file of the pipeline from the 8,241 Carmichael numbers below
# 10^12 (the first 8,241 lines of Shallue and Webster's table) and checks each
# output's size and SHA-256 against golden_prefix.tsv. Every output that has a
# decoder is also decoded and certified against the prefix oracle.
#
#   text --oracle_encode--> prefix.orc1 --oracle_decode--> text (byte for byte)
#   prefix.orc1 --nstats_orc1--> n.raw        --list_build--> c000 (VByte), c001 (classes M=10,810,800)
#   prefix.orc1 --trie_paths---> dpaths.bin   --dlist_from_paths--> d.raw
#                                d.raw        --list_build--> c100 (VByte), c101 (CND1, M=2310)
#                                d.raw        --cnd_encode--> c101 again (must be byte-identical)
#   prefix.orc1 --trie_paths --full--> fpaths.bin --pwt2_encode--> c010 (naive), c011 (set 3)
#                                dpaths.bin   --pwt2_encode--> c110 (naive), c111 (set 3)
#   prefix.orc1 --dsel_set2 (4 sweeps, then 1 more from the masks)--> d* paths --pwt2_encode--> c111x
#
# Decoders, all required: pwt2_decode --certify / --emit-n (all five trees),
# cnd_decode --certify / --emit-d (101), list_decode --emit / --certify /
# --emit-n (000, 001, 100, 101), and pwt2_decode --emit-paths (all five trees).
# Each exported path file must (a) re-encode to its tree file byte for byte,
# (b) hold the same records as the encoder's input as a multiset (for c111x:
# the d* paths of the search), (c) match the SHA-256 the decoder prints, and
# (d) match its golden (tree order). Also checked: thread and sort-run
# invariance of trie_paths, pwt2_encode and dsel_set2 (batch mode), the
# masks-only d* route (--sweeps 0), cnd_encode == list_build (101), and the
# stored copies in golden_files/ (when present) against the rebuilt files.
#
# Usage:
#   tests/prefix/selftest_prefix.sh [--bin DIR] [--work DIR] [--keep]
#                                   [--write-goldens FILE]
#   --bin DIR    directory with the built tools (default: $BIN, else ../../build
#                relative to this script, else ./build)
#   --work DIR   working directory (default: a fresh mktemp directory, removed
#                at the end unless --keep or a check failed)
#   --strict     accepted for compatibility; every check is required
#   --write-goldens FILE  write the size/hash list of this run to FILE (still
#                compares against the existing goldens when they are present)
#
# Inputs next to this script: prefix_1e12.txt (or prefix_1e12.txt.gz) and
# golden_prefix.tsv (frozen 2026-10-06); optionally prefix.orc1 (then it must
# equal the rebuilt oracle) and golden_files/. The text must keep LF line ends
# (no CRLF conversion on checkout: mark tests/prefix/* -text in
# .gitattributes). Exit status: 0 = every check passed, 1 = a check failed,
# 2 = usage or setup error (a missing tool included).
#
# Portable bash (3.2+: macOS, Linux, MSYS2/Git Bash). Needs cmp, awk, wc, od,
# sort, one of sha256sum / shasum / openssl, and gzip for .gz inputs.

set -u
export LC_ALL=C

here=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-}
WORK=
KEEP=0
WRITE_GOLDENS=
while [ $# -gt 0 ]; do
    case "$1" in
        --bin) [ $# -ge 2 ] || { echo "--bin needs a directory" >&2; exit 2; }; BIN=$2; shift 2 ;;
        --work) [ $# -ge 2 ] || { echo "--work needs a directory" >&2; exit 2; }; WORK=$2; shift 2 ;;
        --keep) KEEP=1; shift ;;
        --strict) shift ;;
        --write-goldens) [ $# -ge 2 ] || { echo "--write-goldens needs a file" >&2; exit 2; }; WRITE_GOLDENS=$2; shift 2 ;;
        -h|--help) awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$WRITE_GOLDENS" in ''|/*|[A-Za-z]:*) ;; *) WRITE_GOLDENS="$PWD/$WRITE_GOLDENS" ;; esac

# ------------------------------------------------------------------ setup
if [ -z "$BIN" ]; then
    if [ -d "$here/../../build" ]; then BIN="$here/../../build"; else BIN=./build; fi
fi
[ -d "$BIN" ] || { echo "tool directory not found: $BIN (use --bin DIR)" >&2; exit 2; }
BIN=$(cd "$BIN" && pwd)

exe() { if [ -f "$BIN/$1.exe" ]; then printf '%s' "$BIN/$1.exe"; else printf '%s' "$BIN/$1"; fi; }

if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum < "$1" | awk '{print $1}'; }
elif command -v shasum >/dev/null 2>&1; then
    sha256() { shasum -a 256 < "$1" | awk '{print $1}'; }
elif command -v openssl >/dev/null 2>&1; then
    sha256() { openssl dgst -sha256 < "$1" | awk '{print $NF}'; }
else
    echo "no sha256sum, shasum or openssl found" >&2; exit 2
fi
fbytes() { wc -c < "$1" | tr -d ' \t\r\n'; }

if [ -n "${EPOCHREALTIME:-}" ]; then now() { printf '%s' "$EPOCHREALTIME"; }; else now() { date +%s; }; fi
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.2f", b - a }'; }

GOLDEN="$here/golden_prefix.tsv"
HAVE_GOLDEN=0
[ -f "$GOLDEN" ] && HAVE_GOLDEN=1
if [ $HAVE_GOLDEN -eq 0 ] && [ -z "$WRITE_GOLDENS" ]; then
    echo "golden list not found: $GOLDEN" >&2; exit 2
fi

TOOLS="oracle_encode oracle_decode nstats_orc1 trie_paths dlist_from_paths list_build
       cnd_encode cnd_decode list_decode pwt2_encode pwt2_decode dsel_set2"
for t in $TOOLS; do
    [ -f "$(exe $t)" ] || { echo "missing tool: $(exe $t)" >&2; exit 2; }
done

CREATED_WORK=0
if [ -z "$WORK" ]; then
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/cn-prefix-selftest.XXXXXX") || exit 2
    CREATED_WORK=1
fi
mkdir -p "$WORK/logs" "$WORK/scratch" || exit 2
WORK=$(cd "$WORK" && pwd)
cd "$WORK" || exit 2

NOK=0; NFAIL=0
FAILS=
ok()   { NOK=$((NOK + 1)); printf '  ok       %s\n' "$1"; }
fail() { NFAIL=$((NFAIL + 1)); FAILS="$FAILS
  $1"; printf '  FAIL     %s\n' "$1"; }

T_START=$(now)
step() { printf '\n== %s\n' "$1"; }

# run NAME CMD... : log to logs/NAME.log, report time and exit status
run() {
    local name=$1; shift
    local log="$WORK/logs/$name.log" t0 t1 rc
    printf '%s\n' "$*" > "$log"
    t0=$(now)
    "$@" >> "$log" 2>&1
    rc=$?
    t1=$(now)
    if [ $rc -ne 0 ]; then
        fail "$name: exit $rc ($(elapsed "$t0" "$t1") s; log $log)"
        tail -n 5 "$log" | sed 's/^/           | /'
        return 1
    fi
    printf '  ran      %-34s %6s s\n' "$name" "$(elapsed "$t0" "$t1")"
    return 0
}

# same A B WHAT : byte identity of two files
same() {
    if [ ! -f "$1" ] || [ ! -f "$2" ]; then fail "$3: missing file ($1 / $2)"; return 1; fi
    if cmp -s "$1" "$2"; then ok "$3"; return 0; fi
    fail "$3: files differ ($1 vs $2)"; return 1
}

# grep_log NAME PATTERN WHAT : a tool's own verdict line
grep_log() {
    if grep -q -- "$2" "$WORK/logs/$1.log" 2>/dev/null; then ok "$3"; else fail "$3: '$2' not in logs/$1.log"; fi
}

N=0
GOLD_OUT="$WORK/golden_this_run.tsv"
: > "$GOLD_OUT"
# golden NAME FILE WHAT : size + SHA-256 against golden_prefix.tsv
golden() {
    local name=$1 file=$2 what=$3 b h gb gh bits
    if [ ! -f "$file" ]; then fail "golden $name: not produced"; return 1; fi
    b=$(fbytes "$file"); h=$(sha256 "$file")
    printf '%s\t%s\t%s\t%s\n' "$name" "$b" "$h" "$what" >> "$GOLD_OUT"
    bits=$(awk -v b="$b" -v n="$N" 'BEGIN { if (n > 0) printf "%.3f", 8 * b / n; else printf "-" }')
    if [ $HAVE_GOLDEN -eq 0 ]; then
        printf '  new      %-20s %10s B %8s bits/CN  %s\n' "$name" "$b" "$bits" "$h"; return 0
    fi
    gb=$(awk -F '\t' -v n="$name" '{ sub(/\r$/, "") } $1 == n { print $2 }' "$GOLDEN")
    gh=$(awk -F '\t' -v n="$name" '{ sub(/\r$/, "") } $1 == n { print $3 }' "$GOLDEN")
    if [ -z "$gh" ]; then fail "golden $name: no entry in golden_prefix.tsv (got $b B $h)"; return 1; fi
    if [ "$b" = "$gb" ] && [ "$h" = "$gh" ]; then
        ok "$(printf '%-20s %10s B %8s bits/CN  %.16s' "$name" "$b" "$bits" "$h")"; return 0
    fi
    fail "golden $name: got $b B $h, expected $gb B $gh"; return 1
}

O=prefix.orc1
S=(--scratch "$WORK/scratch")

printf 'prefix self-test\n  tools    %s\n  work     %s\n' "$BIN" "$WORK"
if [ $HAVE_GOLDEN -eq 1 ]; then
    printf '  goldens  %s (sha256 %s, %s entries)\n' "$GOLDEN" "$(sha256 "$GOLDEN")" \
        "$(awk '{ sub(/\r$/, "") } $0 !~ /^#/ && NF >= 3 { k++ } END { print k + 0 }' "$GOLDEN")"
fi
{
    for t in $TOOLS; do
        printf '%s  %s\n' "$(sha256 "$(exe $t)")" "$(basename "$(exe $t)")"
    done
} > "$WORK/binaries.sha256"

# ------------------------------------------------------------------ 0. input text
step "0. input: the 10^12 prefix of the table"
if [ -f "$here/prefix_1e12.txt" ]; then
    cp "$here/prefix_1e12.txt" prefix_1e12.txt
elif [ -f "$here/prefix_1e12.txt.gz" ]; then
    gzip -dc "$here/prefix_1e12.txt.gz" > prefix_1e12.txt
else
    echo "prefix_1e12.txt(.gz) not found next to the script" >&2; exit 2
fi
N=$(wc -l < prefix_1e12.txt | tr -d ' ')
if [ "$(tr -dc '\r' < prefix_1e12.txt | wc -c | tr -d ' ')" != 0 ]; then
    fail "prefix_1e12.txt has CR bytes (CRLF checkout?); the table is LF-only"
fi
golden prefix_1e12.txt prefix_1e12.txt "input: first 8,241 lines (bytes 0..249,923) of the 18,358,310,309-byte table"
if [ "$N" = 8241 ]; then ok "8,241 records"; else fail "record count $N, expected 8241"; fi

# ------------------------------------------------------------------ 1. oracle
step "1. oracle (ORC1) and its byte-for-byte round trip"
run oracle_encode "$(exe oracle_encode)" prefix_1e12.txt $O --prove-primes
golden prefix.orc1 $O "ORC1 oracle, every factor MR-proven"
[ -f "$here/prefix.orc1" ] && same "$here/prefix.orc1" $O "stored prefix.orc1 == rebuilt oracle"
run oracle_decode "$(exe oracle_decode)" $O oracle_decoded.txt && grep_log oracle_decode "round-trip CERTIFIED" "oracle_decode certificate"
same oracle_decoded.txt prefix_1e12.txt "oracle decodes to the input text byte for byte"

# ------------------------------------------------------------------ 2. intermediates
step "2. intermediates: n list, d_min paths, full paths, d_min list"
run nstats_orc1 "$(exe nstats_orc1)" $O --emit-n n.raw --emit-only
golden n.raw n.raw "n list, 16-byte LE, ascending"
run trie_paths_dmin "$(exe trie_paths)" $O dpaths.bin --fields 12 --threads 1
golden dpaths.bin dpaths.bin "d_min paths, 12 x 5-byte fields per CN, n order"
run trie_paths_dmin_t4 "$(exe trie_paths)" $O dpaths_t4.bin --fields 12 --threads 4
same dpaths_t4.bin dpaths.bin "trie_paths d_min: --threads 4 == --threads 1"
run trie_paths_full "$(exe trie_paths)" $O fpaths.bin --full --fields 14 --threads 1
golden fpaths.bin fpaths.bin "full factorizations, 14 x 5-byte fields per CN, n order"
run dlist_from_paths "$(exe dlist_from_paths)" dpaths.bin d.raw --fields 12
golden d.raw d.raw "d_min list, 16-byte LE, ascending"
run dlist_from_orc1 "$(exe dlist_from_paths)" --orc1 $O d_orc1.raw --threads 2
same d_orc1.raw d.raw "d_min list: --orc1 route (choose_dmin) == paths route"

# ------------------------------------------------------------------ 3. list cells
step "3. list cells 000, 001, 100, 101"
run list_000 "$(exe list_build)" vbyte n.raw c000.vbg && grep_log list_000 "VERIFIED" "list_build 000 read-back"
golden c000.vbg c000.vbg "cell 000: n list, VByte gaps (VBG1)"
run list_001 "$(exe list_build)" classes n.raw c001.cnd --mod 10810800 --targets-self && grep_log list_001 "VERIFIED" "list_build 001 read-back"
golden c001.cnd c001.cnd "cell 001: n list, classes mod 10,810,800 (CND1 layout)"
run list_100 "$(exe list_build)" vbyte d.raw c100.vbg && grep_log list_100 "VERIFIED" "list_build 100 read-back"
golden c100.vbg c100.vbg "cell 100: d_min list, VByte gaps (VBG1)"
run list_101 "$(exe list_build)" classes d.raw c101.cnd1 --mod 2310 --targets-orc1 $O && grep_log list_101 "official v1 Cnd1Reader decode_all.*VERIFIED" "list_build 101 read-back"
golden c101.cnd1 c101.cnd1 "cell 101: d_min list, CND1 mod 2310"
run cnd_encode "$(exe cnd_encode)" d.raw $O c101_cnd_encode.cnd1
same c101_cnd_encode.cnd1 c101.cnd1 "cnd_encode == list_build classes (cell 101)"

# ------------------------------------------------------------------ 4. tree cells
step "4. tree cells 010, 011, 110, 111"
run pwt2_010 "$(exe pwt2_encode)" fpaths.bin $O c010.pwt2n --fields 14 --full --naive --threads 1 "${S[@]}"
golden c010.pwt2n c010.pwt2n "cell 010: PWT2 byte-aligned, full factorizations"
run pwt2_011 "$(exe pwt2_encode)" fpaths.bin $O c011.pwt2 --fields 14 --full --ctx-set 3 --threads 1 "${S[@]}"
golden c011.pwt2 c011.pwt2 "cell 011: PWT2 context set 3, full factorizations"
run pwt2_110 "$(exe pwt2_encode)" dpaths.bin $O c110.pwt2n --fields 12 --naive --threads 1 "${S[@]}"
golden c110.pwt2n c110.pwt2n "cell 110: PWT2 byte-aligned, d_min"
run pwt2_111 "$(exe pwt2_encode)" dpaths.bin $O c111.pwt2 --fields 12 --ctx-set 3 --threads 1 "${S[@]}"
golden c111.pwt2 c111.pwt2 "cell 111: PWT2 context set 3, d_min"
# the external sort with many runs and threads must not change a byte
run pwt2_010_runs "$(exe pwt2_encode)" fpaths.bin $O c010_runs.pwt2n --fields 14 --full --naive --threads 3 --chunk 777 "${S[@]}"
same c010_runs.pwt2n c010.pwt2n "pwt2_encode 010: 11 sort runs, 3 threads == 1 run"
run pwt2_111_runs "$(exe pwt2_encode)" dpaths.bin $O c111_runs.pwt2 --fields 12 --ctx-set 3 --threads 4 --chunk 1000 "${S[@]}"
same c111_runs.pwt2 c111.pwt2 "pwt2_encode 111: 9 sort runs, 4 threads == 1 run"

# ------------------------------------------------------------------ 5. d* search
step "5. d* search (dsel_set2, the published settings) and cell 111 with d*"
DS="--threads 6 --batch 65536 --max-overflow 60000000 --huge 128 --cb 128 --slack 16"
run dsel_sweeps4 "$(exe dsel_set2)" --orc1 $O --dpaths dpaths.bin --sweeps 4 $DS "${S[@]}" \
    --masks-out dstar_masks4.bin --emit dstar_paths4.bin
golden dstar_masks4.bin dstar_masks4.bin "d* masks after 4 sweeps (u64 count + u16 mask per CN)"
run dsel_sweep5 "$(exe dsel_set2)" --orc1 $O --dpaths dpaths.bin --masks-in dstar_masks4.bin --sweeps 1 $DS "${S[@]}" \
    --masks-out dstar_masks5.bin --emit dstar_paths5.bin
golden dstar_masks5.bin dstar_masks5.bin "d* masks after a 5th sweep resumed from the masks (as masks_full5.bin)"
golden dstar_paths5.bin dstar_paths5.bin "d* paths, 14 x 5-byte fields per CN (as s2_full5.bin)"
run dsel_masks_route "$(exe dsel_set2)" --orc1 $O --masks-in dstar_masks5.bin --sweeps 0 $DS "${S[@]}" --emit dstar_paths5_rt.bin
same dstar_paths5_rt.bin dstar_paths5.bin "d* paths from the saved masks (--sweeps 0) == search output"
run dsel_threads2 "$(exe dsel_set2)" --orc1 $O --dpaths dpaths.bin --sweeps 4 --threads 2 --batch 65536 \
    --max-overflow 60000000 --huge 128 --cb 128 --slack 16 "${S[@]}" --masks-out dstar_masks4_t2.bin
same dstar_masks4_t2.bin dstar_masks4.bin "dsel_set2 batch mode: --threads 2 == --threads 6"
run pwt2_111x "$(exe pwt2_encode)" dstar_paths5.bin $O c111x.pwt2 --fields 14 --ctx-set 3 --threads 1 "${S[@]}"
golden c111x.pwt2 c111x.pwt2 "cell 111 with d*: PWT2 context set 3 over the d* divisors"

# ------------------------------------------------------------------ 6. decode + certify
step "6. decode and certify against the prefix oracle"
for f in c010.pwt2n c011.pwt2 c110.pwt2n c111.pwt2 c111x.pwt2; do
    tag=${f%%.*}
    run "pwt2_decode_$tag" "$(exe pwt2_decode)" $f --certify $O --emit-n "dec_n_$tag.raw" \
        && grep_log "pwt2_decode_$tag" "^CERTIFIED" "pwt2_decode --certify $f"
    same "dec_n_$tag.raw" n.raw "pwt2_decode --emit-n $f == n.raw"
done
run cnd_decode_101 "$(exe cnd_decode)" c101.cnd1 --emit-d dec_d_c101.raw --certify $O --threads 1 \
    && grep_log cnd_decode_101 "round-trip CERTIFIED" "cnd_decode --certify c101.cnd1"
same dec_d_c101.raw d.raw "cnd_decode --emit-d c101.cnd1 == d.raw"

# ------------------------------------------------------------------ 7. list_decode + emit-paths
step "7. list_decode (every list cell) and pwt2_decode --emit-paths (every tree)"
# list_decode decodes and certifies the flat-cell files (VBG1, and the CND
# layout at any M):
#     list_decode <in.vbg|in.cnd> [--emit <out.raw>] [--kind n|d]
#                 [--emit-n <out.raw>] [--certify <table.orc1>] [--threads N]
# Every list decodes to its input stream (n.raw or d.raw); the n lists (000,
# 001) certify directly; the d lists (100, 101) certify after factoring every
# d (n = d * (d^-1 mod lambda(d))), and every --emit-n stream equals n.raw.
check_list() {         # $1 = cell tag, $2 = file, $3 = stored stream, $4 = kind n|d
    run "list_decode_$1" "$(exe list_decode)" "$2" --emit "dec_$1.raw" \
        && grep_log "list_decode_$1" "^DECODED" "list_decode $2 structure"
    same "dec_$1.raw" "$3" "list_decode --emit $2 == $3"
    run "list_decode_${1}_certify" "$(exe list_decode)" "$2" --kind "$4" --certify $O --threads 1 \
        --emit-n "decn_$1.raw" && grep_log "list_decode_${1}_certify" "^CERTIFIED" "list_decode --kind $4 --certify $2"
    same "decn_$1.raw" n.raw "list_decode --emit-n $2 == n.raw"
}
check_list c000 c000.vbg  n.raw n
check_list c001 c001.cnd  n.raw n
check_list c100 c100.vbg  d.raw d
check_list c101 c101.cnd1 d.raw d

# pwt2_decode --emit-paths exports a tree's stored paths in the trie_paths
# record layout (5-byte big-endian fields, zero-padded), in tree order:
#     pwt2_decode <in.pwt2> --emit-paths <out.bin> [--fields 12|14]
# For each tree: the export re-encodes, with the cell's own flags, to the tree
# file byte for byte; its records are the encoder's input records as a
# multiset (for c111x, the d* paths of the search, 14 fields as in the
# published s2_full5.bin); the SHA-256 the decoder prints is the file's; and
# the export itself (tree order) matches its golden.

# records FILE WIDTH : one hex line per WIDTH-byte record, sorted (cached in FILE.rec)
records() {
    if [ ! -f "$1.rec" ]; then
        od -A n -v -t x1 "$1" | awk -v r="$2" '
            { for (i = 1; i <= NF; i++) { s = s $i; if (++k == r) { print s; s = ""; k = 0 } } }
            END { if (k) print "partial-record " s }' | sort > "$1.rec"
    fi
    printf '%s' "$1.rec"
}
check_paths() {        # $1 = file, $2 = fields, $3 = re-encode flags, $4 = encoder input
    local tag=${1%%.*} sz h
    run "emit_paths_$tag" "$(exe pwt2_decode)" "$1" --emit-paths "paths_$tag.bin" --fields "$2" || return 0
    sz=$(fbytes "paths_$tag.bin")
    if [ "$sz" -ne $((N * 5 * $2)) ]; then fail "emit-paths $1: $sz B, expected N x 5 x $2"; return 0; fi
    h=$(sha256 "paths_$tag.bin")
    grep_log "emit_paths_$tag" "^paths .* sha256 $h" "pwt2_decode --emit-paths $1: printed SHA-256 == file"
    same "$(records "paths_$tag.bin" $((5 * $2)))" "$(records "$4" $((5 * $2)))" \
        "paths exported from $1 == $4 as a multiset of records"
    golden "paths_$tag.bin" "paths_$tag.bin" "pwt2_decode --emit-paths $1 --fields $2 (tree order)"
    # shellcheck disable=SC2086
    run "reencode_$tag" "$(exe pwt2_encode)" "paths_$tag.bin" $O "re_$1" --fields "$2" $3 --threads 1 "${S[@]}" &&
        same "re_$1" "$1" "re-encode of the paths exported from $1 == $1"
}
if "$(exe pwt2_decode)" 2>&1 | grep -q -- '--emit-paths'; then
    check_paths c010.pwt2n 14 "--full --naive"     fpaths.bin
    check_paths c011.pwt2  14 "--full --ctx-set 3" fpaths.bin
    check_paths c110.pwt2n 12 "--naive"            dpaths.bin
    check_paths c111.pwt2  12 "--ctx-set 3"        dpaths.bin
    check_paths c111x.pwt2 14 "--ctx-set 3"        dstar_paths5.bin
else
    fail "pwt2_decode has no --emit-paths option (a build older than 2026-10-06)"
fi

# ------------------------------------------------------------------ 8. stored copies, completeness, summary
step "8. stored copies and summary"
# golden_files/ (optional) holds the nine prefix cell files; each must equal
# the file rebuilt above (.gz copies are compared after gzip -dc).
if [ -d "$here/golden_files" ]; then
    for g in "$here"/golden_files/*; do
        [ -f "$g" ] || continue
        base=$(basename "$g"); name=${base%.gz}
        if [ "$base" != "$name" ]; then gzip -dc "$g" > "stored_$name"; else cp "$g" "stored_$name"; fi
        same "stored_$name" "$name" "stored golden_files/$base == rebuilt $name"
    done
fi
if [ $HAVE_GOLDEN -eq 1 ]; then
    missing=$(awk -F '\t' '{ sub(/\r$/, "") } NR == FNR { if ($0 !~ /^#/ && NF >= 3) seen[$1] = 1; next }
                           $0 !~ /^#/ && NF >= 3 && !($1 in seen) { print $1 }' "$GOLD_OUT" "$GOLDEN")
    if [ -n "$missing" ]; then fail "golden entries not produced: $(echo $missing)"; else ok "every golden entry was produced and checked"; fi
fi
if [ -n "$WRITE_GOLDENS" ]; then
    {
        printf '# golden_prefix.tsv -- size and SHA-256 of every file the 10^12-prefix self-test builds\n'
        printf '# (tests/prefix/selftest_prefix.sh). Input: the first 8,241 lines = bytes 0..249,923 of\n'
        printf '# the 18,358,310,309-byte table (sha256 c934af2a...), all CNs below 10^12.\n'
        printf '# Outputs depend only on the tools and the flags in the script (thread counts\n'
        printf '# included for dsel_set2: --threads 1 runs a different, sequential search).\n'
        printf '# The paths_*.bin entries are pwt2_decode --emit-paths exports, in tree order.\n'
        printf '# name<TAB>bytes<TAB>sha256<TAB>what\n'
        cat "$GOLD_OUT"
    } > "$WRITE_GOLDENS"
    printf '  wrote    %s\n' "$WRITE_GOLDENS"
fi

T_END=$(now)
printf '\n%d ok, %d failed; %s s total; logs in %s/logs\n' \
    "$NOK" "$NFAIL" "$(elapsed "$T_START" "$T_END")" "$WORK"
[ $NFAIL -gt 0 ] && printf 'failures:%s\n' "$FAILS"
rc=0
[ $NFAIL -gt 0 ] && rc=1
if [ $rc -eq 0 ] && [ $CREATED_WORK -eq 1 ] && [ $KEEP -eq 0 ]; then
    cd / && rm -rf "$WORK"
fi
exit $rc
