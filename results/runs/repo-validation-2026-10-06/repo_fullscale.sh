#!/usr/bin/env bash
# Full-scale exercise of the repository's make targets (clean-copy validation).
# Reads the hosted copies in release_data and the rebuilt C = 0 files in
# release_build; writes only under release_build/repo_work. Certified
# originals are only read (copied into repo_work where a rebuild needs an
# intermediate).
export PATH=<toolchain>/bin:$PATH
RC=<repo>
export DATA=<data>
export WORK=<work>
L=<work>/logs
mkdir -p "$L"
cd "$RC" || exit 2
S=$L/summary.tsv
: > "$S"

job() {   # job NAME CMD...
    local name=$1; shift
    local t0 t1 rc
    t0=$(date +%s.%N)
    { echo "cmd:   $*"; echo "start: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "----"; } > "$L/$name.txt"
    "$@" >> "$L/$name.txt" 2>&1
    rc=$?
    t1=$(date +%s.%N)
    { echo "----"; echo "end:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "exit:  $rc"; printf 'seconds: %.1f\n' "$(awk -v a="$t0" -v b="$t1" 'BEGIN { print b - a }')"; } >> "$L/$name.txt"
    printf '%s\t%s\t%.1f\t%s\n' "$name" "$rc" "$(awk -v a="$t0" -v b="$t1" 'BEGIN { print b - a }')" "$*" >> "$S"
    echo "$name exit $rc"
}

step=${1:-all}
if [ "$step" = all ] || [ "$step" = verify ]; then
  for c in 001 101 011 111 111x oracle; do job "verify-$c" make "verify-$c"; done
  RB=<scratch>
  FILE=$RB/c000.vbg   job verify-000 make verify-000
  FILE=$RB/c100.vbg   job verify-100 make verify-100
  FILE=$RB/c010.pwt2n job verify-010 make verify-010
  FILE=$RB/c110.pwt2n job verify-110 make verify-110
fi
if [ "$step" = all ] || [ "$step" = rebuild ]; then
  mkdir -p "$WORK"
  [ -f "$WORK/n.raw" ] || job copy-n-raw cp <certified>/n.raw "$WORK/n.raw"
  [ -f "$WORK/d.raw" ] || job copy-d-raw cp <certified>/d.raw "$WORK/d.raw"
  [ -f "$WORK/dpaths.bin" ] || job copy-dpaths cp <certified>/dpaths.bin "$WORK/dpaths.bin"
  job rebuild-001 make rebuild-001
  job rebuild-101 make rebuild-101
  ROUTE=masks job rebuild-111x-masks make rebuild-111x
fi
if [ "$step" = all ] || [ "$step" = certify ]; then
  job certify-001 make certify-001
  job certify-111 make certify-111
fi
echo "fullscale $step done"
