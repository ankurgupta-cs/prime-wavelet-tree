#!/usr/bin/env bash
# test_release_tools.sh -- fixture checks for list_decode, pwt2_decode --emit-paths
# and paths_cmp (run by `make test-release-tools`; usage:
# bash tests/test_release_tools.sh [build-dir], from the repository root):
#   list_decode            VBG1 and CND-at-any-M files, kinds n and d, --emit,
#                          --emit-n, --certify; cross-checked with cnd_decode;
#                          plain and crafted corruption cases (CND cases with a
#                          recomputed payload SHA-256 reach the deeper checks)
#   pwt2_decode --emit-paths  encode -> decode --emit-paths -> re-encode must be
#                          byte-identical (d_min, full and a d*-like divisor
#                          tree; context sets 1-3 and byte-aligned); the
#                          emitted records are the input records as a multiset;
#                          failure cases leave no output behind
#   paths_cmp              path files compared as sets (emitted vs input,
#                          multi-run sorts, exact common / one-sided counts
#                          vs comm, duplicate / truncated / unsorted records)
# Input: tests/fixtures/fix_{head,mid,tail}.txt (637 table lines) plus the two
# smallest Carmichael numbers with 13 and 14 prime factors (OEIS A006931; the
# lines as the certified oracle stores them), so full-factorization paths fill
# all 14 fields. Everything is built from that text by the project's own
# tools (oracle_encode --prove-primes, nstats_orc1, trie_paths,
# dlist_from_paths, list_build, pwt2_encode). Scratch: a fresh temp directory,
# removed at the end unless KEEP=1. Needs bash, perl (core Digest::SHA), cmp,
# od, sort, sha256sum, stat. Exit 0 iff every check passes.
set -u
B=${1:-build}
EXE=""
[ -f "$B/list_decode.exe" ] && EXE=.exe
for t in list_decode oracle_encode nstats_orc1 trie_paths dlist_from_paths list_build pwt2_encode pwt2_decode cnd_decode paths_cmp; do
    [ -x "$B/$t$EXE" ] || { echo "missing $B/$t$EXE (make test-release-tools builds it)"; exit 2; }
done
FIX=${FIX:-tests/fixtures}
[ -f "$FIX/fix_head.txt" ] || { echo "run from the repository root ($FIX/ not found)"; exit 2; }
T=$(mktemp -d "${TMPDIR:-/tmp}/cn_reltools.XXXXXX") || exit 2
mkdir -p "$T/s"
cleanup() { if [ "${KEEP:-0}" = 1 ]; then echo "kept $T"; else rm -rf "$T"; fi; }
trap cleanup EXIT

npass=0
nfail=0
ok()  { echo "PASS  $1"; npass=$((npass + 1)); }
bad() {
    echo "FAIL  $1"
    nfail=$((nfail + 1))
    if [ -n "${2:-}" ] && [ -f "$2" ]; then tail -n 4 "$2" | sed 's/^/      | /'; fi
}
die() { echo "setup failed: $1"; [ -n "${2:-}" ] && [ -f "$2" ] && tail -n 5 "$2"; exit 1; }
setup() { local log=$1; shift; "$@" > "$log" 2>&1 || die "$*" "$log"; }
sz() { stat -c %s "$1"; }
cmpf() { if cmp -s "$2" "$3"; then ok "$1"; else bad "$1"; fi; }

# perl helpers: xor one byte; patch a little-endian integer field; insert
# bytes; recompute a CND file's payload SHA-256 (footer bytes 0-31 over
# header || directory || body).
xorbyte() { perl -e 'my ($f, $o, $x) = @ARGV; open(my $h, "+<:raw", $f) or die; seek($h, $o, 0);
    read($h, my $b, 1) == 1 or die "offset"; seek($h, $o, 0); print $h chr(ord($b) ^ hex($x)); close $h or die;' "$1" "$2" "$3"; }
setle() { perl -e 'my ($f, $o, $w, $v) = @ARGV; open(my $h, "+<:raw", $f) or die; seek($h, $o, 0);
    print $h substr(pack("Q<", $v), 0, $w); close $h or die;' "$1" "$2" "$3" "$4"; }
addle() { perl -e 'my ($f, $o, $w, $d) = @ARGV; open(my $h, "+<:raw", $f) or die; seek($h, $o, 0);
    read($h, my $b, $w) == $w or die; my $v = unpack("Q<", $b . ("\0" x (8 - $w))); seek($h, $o, 0);
    print $h substr(pack("Q<", $v + $d), 0, $w); close $h or die;' "$1" "$2" "$3" "$4"; }
insbytes() { perl -e 'my ($f, $o, $hex) = @ARGV; local $/; open(my $h, "<:raw", $f) or die; my $s = <$h>; close $h;
    substr($s, $o, 0) = pack("H*", $hex); open($h, ">:raw", $f) or die; print $h $s; close $h or die;' "$1" "$2" "$3"; }
reseal() { perl -MDigest::SHA=sha256 -e 'my ($f) = @ARGV; local $/; open(my $h, "<:raw", $f) or die; my $s = <$h>; close $h;
    my $fo = length($s) - 88; substr($s, $fo, 32) = sha256(substr($s, 0, $fo));
    open($h, ">:raw", $f) or die; print $h $s; close $h or die;' "$1"; }
raw16() { local out=$1; shift; perl -e 'binmode STDOUT; print pack("Q<Q<", $_, 0) for @ARGV' "$@" > "$out"; }

# ---------------------------------------------------------------------------
echo "== setup: fixture oracle and intermediates in $T"
{
    cat "$FIX/fix_head.txt"
    echo "1791562810662585767521 11 13 17 19 31 37 43 71 73 97 109 113 127"
    echo "87674969936234821377601 7 13 17 19 23 31 37 41 61 67 89 163 193 241"
    cat "$FIX/fix_mid.txt" "$FIX/fix_tail.txt"
} > "$T/fix.txt"
NREC=$(wc -l < "$T/fix.txt" | tr -d ' ')
setup "$T/oe.log" "$B/oracle_encode$EXE" "$T/fix.txt" "$T/fix.orc1" --prove-primes
setup "$T/oe2.log" "$B/oracle_encode$EXE" "$FIX/fix_head.txt" "$T/other.orc1" --prove-primes
setup "$T/ns.log" "$B/nstats_orc1$EXE" "$T/fix.orc1" --emit-n "$T/n.raw" --emit-only
setup "$T/tp1.log" "$B/trie_paths$EXE" "$T/fix.orc1" "$T/dpaths.bin" --threads 1
setup "$T/tp2.log" "$B/trie_paths$EXE" "$T/fix.orc1" "$T/fpaths.bin" --full --threads 1
setup "$T/dl.log" "$B/dlist_from_paths$EXE" "$T/dpaths.bin" "$T/d.raw"
[ "$(sz "$T/n.raw")" = $((NREC * 16)) ] && [ "$(sz "$T/d.raw")" = $((NREC * 16)) ] || die "n.raw / d.raw sizes"
# a d*-like divisor path file (14 fields): every third record stores d = n
# (its full factorization; n * (n^-1 mod lambda(n)) = n), the rest d_min
perl -e 'open(D, "<:raw", $ARGV[0]) or die; open(F, "<:raw", $ARGV[1]) or die; open(O, ">:raw", $ARGV[2]) or die;
    my $i = 0; while (read(D, my $d, 60) == 60) { read(F, my $f, 70) == 70 or die "fpaths short";
    print O ($i % 3 == 0 ? $f : $d . ("\0" x 10)); $i++ } read(F, my $x, 1) and die "fpaths long"; close O or die;' \
    "$T/dpaths.bin" "$T/fpaths.bin" "$T/xpaths.bin" || die "xpaths"
LB="$B/list_build$EXE"
setup "$T/lb1.log" "$LB" vbyte "$T/n.raw" "$T/c000.vbg"
setup "$T/lb2.log" "$LB" vbyte "$T/d.raw" "$T/c100.vbg"
setup "$T/lb3.log" "$LB" classes "$T/n.raw" "$T/c001.cnd" --mod 10810800 --targets-self
setup "$T/lb4.log" "$LB" classes "$T/d.raw" "$T/c101.cnd" --mod 2310 --targets-orc1 "$T/fix.orc1"
setup "$T/lb5.log" "$LB" classes "$T/n.raw" "$T/c001_m510510.cnd" --mod 510510 --targets-self
setup "$T/lb6.log" "$LB" classes "$T/d.raw" "$T/c101_m1.cnd" --mod 1 --targets-raw "$T/n.raw"
setup "$T/lb7.log" "$LB" classes "$T/d.raw" "$T/c101_m3.cnd" --mod 3 --targets-orc1 "$T/fix.orc1"
setup "$T/lb8.log" "$LB" classes "$T/n.raw" "$T/c001_none.cnd" --mod 510510
echo "   $NREC records; flat files: c000 $(sz "$T/c000.vbg") B, c100 $(sz "$T/c100.vbg") B, c001 $(sz "$T/c001.cnd") B, c101 $(sz "$T/c101.cnd") B"

# ---------------------------------------------------------------------------
# ld NAME RC WANT SUBSTR -- args: run list_decode; require exit RC, and unless
# WANT is "-", a final CERTIFIED/DECODED/FAILED line starting with WANT and
# containing SUBSTR.
ld() {
    local name=$1 want_rc=$2 want=$3 sub=$4
    shift 4
    local log="$T/ld_$name.log"
    "$B/list_decode$EXE" "$@" > "$log" 2>&1
    local rc=$?
    local last
    last=$(grep -E '^(CERTIFIED|DECODED|FAILED) ' "$log" | tail -n 1)
    if [ "$rc" -eq "$want_rc" ] && { [ "$want" = "-" ] || { [[ "$last" == "$want"* ]] && [[ "$last" == *"$sub"* ]]; }; }; then
        ok "list_decode $name"
    else
        bad "list_decode $name (exit $rc, want $want_rc $want '$sub')" "$log"
    fi
}

echo "== list_decode: certified round trips"
O="$T/fix.orc1"
ld c000_n 0 CERTIFIED "" "$T/c000.vbg" --kind n --certify "$O" --emit "$T/o000.raw" --emit-n "$T/on000.raw"
cmpf "c000 --emit == n.raw" "$T/o000.raw" "$T/n.raw"
cmpf "c000 --emit-n == n.raw" "$T/on000.raw" "$T/n.raw"
ld c001_n 0 CERTIFIED "footer targets match" "$T/c001.cnd" --kind n --certify "$O" --emit "$T/o001.raw"
cmpf "c001 (M = 10,810,800) --emit == n.raw" "$T/o001.raw" "$T/n.raw"
ld c001_m510510_n 0 CERTIFIED "footer targets match" "$T/c001_m510510.cnd" --kind n --certify "$O"
ld c100_d 0 CERTIFIED "factored" "$T/c100.vbg" --kind d --certify "$O" --threads 2 --emit "$T/o100.raw" --emit-n "$T/on100.raw"
cmpf "c100 --emit == d.raw" "$T/o100.raw" "$T/d.raw"
cmpf "c100 --emit-n == n.raw" "$T/on100.raw" "$T/n.raw"
ld c101_d 0 CERTIFIED "footer targets match" "$T/c101.cnd" --kind d --certify "$O" --threads 3 --emit "$T/o101.raw" --emit-n "$T/on101.raw"
cmpf "c101 --emit == d.raw" "$T/o101.raw" "$T/d.raw"
cmpf "c101 --emit-n == n.raw" "$T/on101.raw" "$T/n.raw"
ld c101_m1_d 0 CERTIFIED "footer targets match" "$T/c101_m1.cnd" --kind d --certify "$O" --threads 1
ld c101_m3_d 0 CERTIFIED "footer targets match" "$T/c101_m3.cnd" --kind d --certify "$O"
ld c101_decode_only 0 DECODED "" "$T/c101.cnd" --emit "$T/o101b.raw"
cmpf "c101 decode-only --emit == d.raw" "$T/o101b.raw" "$T/d.raw"
vs=$(sed -n 's/^values  *sha256 \([0-9a-f]*\).*/\1/p' "$T/ld_c101_decode_only.log")
[ "$vs" = "$(sha256sum "$T/d.raw" | cut -c1-64)" ] && ok "values sha256 line == sha256 of d.raw" || bad "values sha256 line == sha256 of d.raw"
fs=$(sed -n 's/^file  *sha256 \([0-9a-f]*\).*/\1/p' "$T/ld_c101_decode_only.log")
[ "$fs" = "$(sha256sum "$T/c101.cnd" | cut -c1-64)" ] && ok "file sha256 line == sha256 of the file" || bad "file sha256 line == sha256 of the file"
# an n list read as divisors certifies too: d = n is a valid divisor (Korselt)
ld c000_as_d 0 CERTIFIED "" "$T/c000.vbg" --kind d --certify "$O" --threads 2
# cross-check with the shipped v1 decoder on the M = 2310 file
if "$B/cnd_decode$EXE" "$T/c101.cnd" --emit-d "$T/cd101.raw" --certify "$O" --threads 2 > "$T/cnd_decode.log" 2>&1 &&
    grep -q "round-trip CERTIFIED" "$T/cnd_decode.log"; then ok "cnd_decode certifies the same c101"; else bad "cnd_decode certifies the same c101" "$T/cnd_decode.log"; fi
cmpf "cnd_decode --emit-d == list_decode --emit" "$T/cd101.raw" "$T/o101.raw"

echo "== list_decode: wrong claims and usage"
ld d_list_as_n_cnd 1 FAILED "footer targets" "$T/c101.cnd" --kind n --certify "$O"
ld d_list_as_n_vbg 1 FAILED "mismatches" "$T/c100.vbg" --kind n --certify "$O"
ld wrong_oracle 1 FAILED "count" "$T/c000.vbg" --kind n --certify "$T/other.orc1"
ld no_targets 1 FAILED "footer targets" "$T/c001_none.cnd" --kind n --certify "$O"
ld usage_no_kind 2 - "" "$T/c000.vbg" --certify "$O"
ld usage_bad_kind 2 - "" "$T/c000.vbg" --kind x
printf 'not a list file, just text....' > "$T/junk.bin"
ld bad_magic 1 FAILED "unknown format" "$T/junk.bin"
raw16 "$T/sq.raw" 45 51
setup "$T/lb_sq.log" "$LB" vbyte "$T/sq.raw" "$T/sq.vbg"
ld d_not_squarefree 1 FAILED "not squarefree" "$T/sq.vbg" --kind d --emit-n "$T/sq_n.raw" --threads 1
raw16 "$T/ev.raw" 50 51
setup "$T/lb_ev.log" "$LB" vbyte "$T/ev.raw" "$T/ev.vbg"
ld d_even 1 FAILED "not odd" "$T/ev.vbg" --kind d --emit-n "$T/ev_n.raw" --threads 1

echo "== list_decode: corrupted VBG1 (no checksum in the format: certification catches content changes)"
V="$T/c000.vbg"
P=$(( $(sz "$V") - 24 ))
cp "$V" "$T/v1.vbg"; xorbyte "$T/v1.vbg" $((24 + P / 2)) 01
ld vbg_flip_certify 1 FAILED "" "$T/v1.vbg" --kind n --certify "$O"
head -c $(( $(sz "$V") - 1 )) "$V" > "$T/v2.vbg"
ld vbg_truncated 1 FAILED "payload size" "$T/v2.vbg"
cp "$V" "$T/v3.vbg"; printf '\001' >> "$T/v3.vbg"
ld vbg_trailing_byte 1 FAILED "payload size" "$T/v3.vbg"
cp "$V" "$T/v4.vbg"; addle "$T/v4.vbg" 8 8 1
ld vbg_count_plus1 1 FAILED "malformed" "$T/v4.vbg"
cp "$V" "$T/v5.vbg"; addle "$T/v5.vbg" 8 8 -1
ld vbg_count_minus1 1 FAILED "left after" "$T/v5.vbg"
# 561 = b1 04; non-canonical b1 84 00 (one more payload byte)
cp "$V" "$T/v6.vbg"; [ "$(od -An -tx1 -j24 -N2 "$V" | tr -d ' ')" = "b104" ] || die "c000 does not start with 561"
xorbyte "$T/v6.vbg" 25 80; insbytes "$T/v6.vbg" 26 00; addle "$T/v6.vbg" 16 8 1
ld vbg_noncanonical 1 FAILED "non-canonical" "$T/v6.vbg"
# a zero gap: one extra record whose gap is 0
cp "$V" "$T/v7.vbg"; insbytes "$T/v7.vbg" 26 00; addle "$T/v7.vbg" 8 8 1; addle "$T/v7.vbg" 16 8 1
ld vbg_zero_gap 1 FAILED "zero gap at record 1" "$T/v7.vbg"

echo "== list_decode: corrupted CND (plain, then crafted with a valid payload SHA-256)"
C="$T/c101.cnd"
CS=$(sz "$C")
cp "$C" "$T/k1.cnd"; xorbyte "$T/k1.cnd" $(( CS - 88 - 10 )) 04
ld cnd_body_flip 1 FAILED "payload sha mismatch" "$T/k1.cnd"
head -c $(( CS - 1 )) "$C" > "$T/k2.cnd"
ld cnd_truncated 1 FAILED "" "$T/k2.cnd"
cp "$C" "$T/k3.cnd"; xorbyte "$T/k3.cnd" $(( CS - 88 + 40 )) 01      # footer sha_nset: outside the payload hash
ld cnd_footer_target_decode_only 0 DECODED "" "$T/k3.cnd"
ld cnd_footer_target_certify 1 FAILED "footer targets" "$T/k3.cnd" --kind d --certify "$O" --threads 2
cp "$C" "$T/k4.cnd"; setle "$T/k4.cnd" 8 4 4294967295; reseal "$T/k4.cnd"
ld cnd_huge_M 1 FAILED "M exceeds the directory size" "$T/k4.cnd"
cp "$C" "$T/k5.cnd"; setle "$T/k5.cnd" 8 4 2311; reseal "$T/k5.cnd"
ld cnd_M_plus1 1 FAILED "corrupt directory (count of class 2310)" "$T/k5.cnd"
cp "$C" "$T/k6.cnd"; addle "$T/k6.cnd" 16 8 1; reseal "$T/k6.cnd"
ld cnd_records_plus1 1 FAILED "record count mismatch" "$T/k6.cnd"
# class 0 is empty (count byte 00): write it non-canonically as 80 00
cp "$C" "$T/k7.cnd"; [ "$(od -An -tx1 -j64 -N1 "$C" | tr -d ' ')" = "00" ] || die "class 0 of c101 not empty"
insbytes "$T/k7.cnd" 64 80; addle "$T/k7.cnd" 32 8 1; addle "$T/k7.cnd" 40 8 1; reseal "$T/k7.cnd"
ld cnd_noncanonical_dir 1 FAILED "corrupt directory (count of class 0)" "$T/k7.cnd"
cp "$C" "$T/k8.cnd"; setle "$T/k8.cnd" 48 1 1; reseal "$T/k8.cnd"
ld cnd_reserved_byte 1 FAILED "reserved" "$T/k8.cnd"
# a body bit flipped and the payload hash recomputed: structure passes,
# the content change is caught by decode or certification
cp "$C" "$T/k9.cnd"; xorbyte "$T/k9.cnd" $(( CS - 88 - 700 )) 10; reseal "$T/k9.cnd"
ld cnd_crafted_body_certify 1 FAILED "" "$T/k9.cnd" --kind d --certify "$O" --threads 2

# ---------------------------------------------------------------------------
echo "== pwt2_decode --emit-paths: encode -> emit paths -> re-encode"
# pw NAME PATHS K ENCODER-FLAGS...
pw() {
    local name=$1 paths=$2 K=$3
    shift 3
    if ! "$B/pwt2_encode$EXE" "$T/$paths" "$O" "$T/$name.a.pwt2" --fields "$K" "$@" --scratch "$T/s" > "$T/pe_$name.a.log" 2>&1; then
        bad "pwt2 $name: encode" "$T/pe_$name.a.log"; return
    fi
    if ! "$B/pwt2_decode$EXE" "$T/$name.a.pwt2" --certify "$O" --emit-paths "$T/$name.paths" --fields "$K" > "$T/pd_$name.log" 2>&1 ||
        ! grep -q '^CERTIFIED' "$T/pd_$name.log"; then
        bad "pwt2 $name: decode --emit-paths --certify" "$T/pd_$name.log"; return
    fi
    local w=$((K * 5))
    od -An -v -tx1 -w$w "$T/$paths" | sort > "$T/$name.in.sorted"
    od -An -v -tx1 -w$w "$T/$name.paths" | sort > "$T/$name.out.sorted"
    if [ "$(sz "$T/$name.paths")" = $((NREC * w)) ] && [ "$(wc -l < "$T/$name.in.sorted")" -eq "$NREC" ] &&
        cmp -s "$T/$name.in.sorted" "$T/$name.out.sorted"; then
        ok "pwt2 $name: emitted paths == encoder input (as a multiset of $NREC records x $K fields)"
    else
        bad "pwt2 $name: emitted paths == encoder input"
    fi
    if ! "$B/pwt2_encode$EXE" "$T/$name.paths" "$O" "$T/$name.b.pwt2" --fields "$K" "$@" --scratch "$T/s" > "$T/pe_$name.b.log" 2>&1; then
        bad "pwt2 $name: re-encode" "$T/pe_$name.b.log"; return
    fi
    cmpf "pwt2 $name: re-encoded container byte-identical ($(sz "$T/$name.a.pwt2") B)" "$T/$name.a.pwt2" "$T/$name.b.pwt2"
}
pw dmin_s1 dpaths.bin 12 --ctx-set 1
pw dmin_s2 dpaths.bin 12 --ctx-set 2
pw dmin_s3 dpaths.bin 12 --ctx-set 3
pw dmin_naive dpaths.bin 12 --naive
pw full_s1 fpaths.bin 14 --full --ctx-set 1
pw full_s2 fpaths.bin 14 --full --ctx-set 2
pw full_s3 fpaths.bin 14 --full --ctx-set 3
pw full_naive fpaths.bin 14 --full --naive
pw dstar_like_s3 xpaths.bin 14 --ctx-set 3
pw dstar_like_naive xpaths.bin 14 --naive

echo "== pwt2_decode --emit-paths: defaults, field choice, failures"
grep -q "max depth 14" "$T/pd_full_s3.log" && ok "full tree uses all 14 fields (max depth 14)" || bad "full tree max depth 14" "$T/pd_full_s3.log"
"$B/pwt2_decode$EXE" "$T/dmin_s3.a.pwt2" --emit-paths "$T/def12.paths" > "$T/pd_def12.log" 2>&1 &&
    cmp -s "$T/def12.paths" "$T/dmin_s3.paths" && ok "divisor file: default --fields 12" || bad "divisor file: default --fields 12" "$T/pd_def12.log"
"$B/pwt2_decode$EXE" "$T/full_s3.a.pwt2" --emit-paths "$T/def14.paths" > "$T/pd_def14.log" 2>&1 &&
    cmp -s "$T/def14.paths" "$T/full_s3.paths" && ok "full file: default --fields 14" || bad "full file: default --fields 14" "$T/pd_def14.log"
# the container does not depend on the field count: 14-field d_min paths re-encode identically
if "$B/pwt2_decode$EXE" "$T/dmin_s3.a.pwt2" --emit-paths "$T/d14.paths" --fields 14 > "$T/pd_d14.log" 2>&1 &&
    "$B/pwt2_encode$EXE" "$T/d14.paths" "$O" "$T/d14.pwt2" --fields 14 --ctx-set 3 --scratch "$T/s" > "$T/pe_d14.log" 2>&1 &&
    cmp -s "$T/d14.pwt2" "$T/dmin_s3.a.pwt2"; then ok "d_min paths at 14 fields re-encode to the same container"
else bad "d_min paths at 14 fields re-encode to the same container" "$T/pe_d14.log"; fi
"$B/pwt2_decode$EXE" "$T/full_s3.a.pwt2" --emit-paths "$T/f12.paths" --fields 12 > "$T/pd_f12.log" 2>&1
rc=$?
if [ $rc -eq 1 ] && grep -q "exceeds 12 path fields" "$T/pd_f12.log" && [ ! -e "$T/f12.paths" ] && [ ! -e "$T/f12.paths.partial" ]; then
    ok "max depth 14 > --fields 12: refused, nothing written"
else bad "max depth 14 > --fields 12: refused, nothing written (exit $rc)" "$T/pd_f12.log"; fi
"$B/pwt2_decode$EXE" "$T/full_s3.a.pwt2" --emit-paths "$T/f13.paths" --fields 13 > "$T/pd_f13.log" 2>&1
[ $? -eq 2 ] && ok "--fields 13 rejected (usage)" || bad "--fields 13 rejected (usage)" "$T/pd_f13.log"
"$B/pwt2_decode$EXE" "$T/full_s3.a.pwt2" --fields 14 > "$T/pd_fo.log" 2>&1
[ $? -eq 2 ] && ok "--fields without --emit-paths rejected (usage)" || bad "--fields without --emit-paths rejected" "$T/pd_fo.log"
# footer sha_nset changed (outside the payload hash): the walk completes and
# writes the partial file, the target check fails, the partial is removed
cp "$T/dmin_s3.a.pwt2" "$T/bad.pwt2"; xorbyte "$T/bad.pwt2" $(( $(sz "$T/bad.pwt2") - 88 + 40 )) 01
"$B/pwt2_decode$EXE" "$T/bad.pwt2" --emit-paths "$T/bad.paths" > "$T/pd_bad.log" 2>&1
rc=$?
if [ $rc -eq 1 ] && grep -q "sha_nset MISMATCH" "$T/pd_bad.log" && [ ! -e "$T/bad.paths" ] && [ ! -e "$T/bad.paths.partial" ]; then
    ok "target mismatch: exit 1, no paths file or partial left"
else bad "target mismatch: exit 1, no paths file or partial left (exit $rc)" "$T/pd_bad.log"; fi
cp "$T/dmin_s3.a.pwt2" "$T/bad2.pwt2"; xorbyte "$T/bad2.pwt2" 300 01
"$B/pwt2_decode$EXE" "$T/bad2.pwt2" --emit-paths "$T/bad2.paths" > "$T/pd_bad2.log" 2>&1
rc=$?
[ $rc -eq 1 ] && [ ! -e "$T/bad2.paths" ] && [ ! -e "$T/bad2.paths.partial" ] && ok "payload byte flipped: exit 1, nothing written" ||
    bad "payload byte flipped: exit 1, nothing written (exit $rc)" "$T/pd_bad2.log"
# default behaviour (no --emit-paths) unchanged: certify + --emit-n
"$B/pwt2_decode$EXE" "$T/dstar_like_s3.a.pwt2" --certify "$O" --emit-n "$T/xn.raw" > "$T/pd_plain.log" 2>&1 &&
    grep -q '^CERTIFIED' "$T/pd_plain.log" && ! grep -q '^paths' "$T/pd_plain.log" && cmp -s "$T/xn.raw" "$T/n.raw" &&
    ok "plain pwt2_decode --certify --emit-n (no paths line, n.raw identical)" || bad "plain pwt2_decode --certify --emit-n" "$T/pd_plain.log"

# ---------------------------------------------------------------------------
echo "== paths_cmp: path files compared as sets (record order ignored)"
# pc NAME EXPECT-EXIT EXPECT-LAST-WORD ARGS...   (also: no run files left behind)
pc() {
    local name=$1 want=$2 last=$3
    shift 3
    "$B/paths_cmp$EXE" "$@" --scratch "$T/s" > "$T/pc_$name.log" 2>&1
    local rc=$?
    if [ $rc -eq "$want" ] && tail -n 1 "$T/pc_$name.log" | grep -q "^$last" &&
        [ -z "$(ls "$T/s" | grep '^pathscmp_')" ]; then
        ok "paths_cmp $name: exit $rc, $last, no run files left"
    else
        bad "paths_cmp $name (exit $rc, want $want $last)" "$T/pc_$name.log"
    fi
}
pc emitted_vs_input_dstar 0 "SAME SET" "$T/xpaths.bin" "$T/dstar_like_s3.paths" --fields 14
pc emitted_vs_input_dmin_7runs 0 "SAME SET" "$T/dpaths.bin" "$T/dmin_s3.paths" --fields 12 --chunk 100 --threads 3
pc emitted_vs_input_full 0 "SAME SET" "$T/full_s3.paths" "$T/fpaths.bin" --fields 14 --chunk 64
# d*-like vs full paths: the d = n records coincide, the rest differ; the
# counts must equal comm's on the od-sorted records
pc dstar_vs_full 1 DIFFERENT "$T/xpaths.bin" "$T/fpaths.bin" --fields 14 --chunk 50
want_c=$(comm -12 "$T/dstar_like_s3.in.sorted" "$T/full_s3.in.sorted" | wc -l | tr -d ' ')
got=$(grep '^merge' "$T/pc_dstar_vs_full.log")
if [ "$got" = "merge      common $want_c, only in A $((NREC - want_c)), only in B $((NREC - want_c))" ] &&
    [ "$want_c" -gt 0 ] && [ "$want_c" -lt "$NREC" ]; then
    ok "paths_cmp counts equal comm's ($want_c common of $NREC)"
else
    bad "paths_cmp counts vs comm (want common $want_c): $got"
fi
head -c $(( (NREC - 1) * 70 )) "$T/xpaths.bin" > "$T/xminus.bin"
pc one_record_missing 1 DIFFERENT "$T/xminus.bin" "$T/xpaths.bin" --chunk 100
grep -q '^merge      common '$((NREC - 1))', only in A 0, only in B 1$' "$T/pc_one_record_missing.log" &&
    ok "paths_cmp one missing record: only in B 1" || bad "paths_cmp one missing record" "$T/pc_one_record_missing.log"
{ cat "$T/xpaths.bin"; head -c 70 "$T/xpaths.bin"; } > "$T/xdup.bin"
pc duplicate_record_in_B 1 FAILED "$T/xpaths.bin" "$T/xdup.bin" --chunk 100
pc duplicate_record_in_A 1 FAILED "$T/xdup.bin" "$T/xpaths.bin" --chunk 100
{ cat "$T/xpaths.bin"; printf 'abc'; } > "$T/xtrunc.bin"
pc truncated_record 1 FAILED "$T/xpaths.bin" "$T/xtrunc.bin"
# a record with descending fields (7, 5) is not a path record
{ cat "$T/xpaths.bin"; perl -e 'print pack("C5", 0,0,0,0,7), pack("C5", 0,0,0,0,5), "\0" x 60'; } > "$T/xdesc.bin"
pc unsorted_fields 1 FAILED "$T/xdesc.bin" "$T/xpaths.bin"
"$B/paths_cmp$EXE" "$T/xpaths.bin" > /dev/null 2>&1
[ $? -eq 2 ] && ok "paths_cmp one file: usage (exit 2)" || bad "paths_cmp one file: usage"
"$B/paths_cmp$EXE" "$T/xpaths.bin" "$T/xpaths.bin" --fields 13 > /dev/null 2>&1
[ $? -eq 2 ] && ok "paths_cmp --fields 13: usage (exit 2)" || bad "paths_cmp --fields 13: usage"

echo "== $npass passed, $nfail failed"
[ "$nfail" -eq 0 ]
