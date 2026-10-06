#!/bin/bash
# run_job_b1.sh TAG OUT EXPECTED_SHA -- command args...
# Runs one job, logs stdout to logs_b1/TAG.txt and stderr to logs_b1/TAG.err, then
# records exit code, wall seconds, output bytes and SHA-256 (vs EXPECTED) in
# logs_b1/summary.tsv. OUT = "-" for a job without one output file.
export PATH=/usr/bin:$PATH
tag=$1; out=$2; exp=$3; shift 3
[ "$1" = "--" ] && shift
L=<work>/release_build/logs_b1
mkdir -p "$L"
{
  echo "tag:     $tag"
  echo "cmd:     $*"
  echo "cwd:     $(pwd)"
  echo "start:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "----"
} > "$L/$tag.txt"
t0=$(date +%s.%N)
"$@" >> "$L/$tag.txt" 2> "$L/$tag.err"
rc=$?
t1=$(date +%s.%N)
secs=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
bytes="-"; sha="-"; status="NO-OUTPUT"
if [ "$out" != "-" ] && [ -f "$out" ]; then
  bytes=$(stat -c %s "$out")
  h0=$(date +%s.%N)
  sha=$(sha256sum "$out" | cut -d' ' -f1)
  h1=$(date +%s.%N)
  hsecs=$(awk -v a="$h0" -v b="$h1" 'BEGIN{printf "%.1f", b-a}')
  if [ "$exp" = "-" ]; then status="NO-EXPECTED"
  elif [ "$sha" = "$exp" ]; then status="MATCH"
  else status="MISMATCH"; fi
fi
{
  echo "----"
  echo "end:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "exit:    $rc"
  echo "seconds: $secs (process wall, measured by run_job.sh)"
  echo "output:  $out"
  echo "bytes:   $bytes"
  echo "sha256:  $sha (hashing ${hsecs:-0} s)"
  echo "expect:  $exp"
  echo "status:  $status"
} >> "$L/$tag.txt"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tag" "$rc" "$secs" "$bytes" "$sha" "$exp" "$status" "$*" >> "$L/summary.tsv"
echo "$tag exit=$rc secs=$secs bytes=$bytes sha=$sha status=$status"
exit $rc
