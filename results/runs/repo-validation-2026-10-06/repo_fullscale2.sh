#!/usr/bin/env bash
# Second stage: wait for repo_fullscale.sh, then certify cell 101 with the default tool (cnd_decode).
until grep -q 'fullscale all done' <scratch>/repo_fullscale_console.txt 2>/dev/null; do sleep 20; done
export PATH=<toolchain>/bin:$PATH
RC=<repo>
export DATA=<data>
export WORK=<work>
L=<work>/logs
cd "$RC" || exit 2
t0=$(date +%s.%N)
{ echo "cmd:   make certify-101"; echo "start: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "----"; } > "$L/certify-101.txt"
make certify-101 >> "$L/certify-101.txt" 2>&1
rc=$?
t1=$(date +%s.%N)
s=$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.1f", b - a }')
{ echo "----"; echo "end:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "exit:  $rc"; echo "seconds: $s"; } >> "$L/certify-101.txt"
printf 'certify-101\t%s\t%s\tmake certify-101\n' "$rc" "$s" >> "$L/summary.tsv"
echo "certify-101 exit $rc; stage2 done"
