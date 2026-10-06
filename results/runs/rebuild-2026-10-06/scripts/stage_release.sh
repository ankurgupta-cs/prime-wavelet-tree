#!/bin/bash
# Note added for publication (2026-10-06): this staging predates the carmichael-1e24-*
# file names (results/cells.tsv), and c111x-dstar-masks.bin was later dropped from the
# release (the search masks are not hosted). <repo> is the project tree, not this repository.
# Step 6: copy every VERIFIED hosted file into <work>/release_data,
# re-hash each copy, and write SHA256SUMS + MANIFEST.tsv from the copies.
# Sources are read-only; a copy that does not hash to its expected value is
# removed again and reported (never "fixed").
export PATH=/usr/bin:$PATH
R=<work>/release_data
L=<work>/release_build/logs
mkdir -p "$R"
cd "$R" || exit 2
: > "$L/stage_release.tsv"
stage() { # name cell expected source
  name=$1; cell=$2; exp=$3; src=$4
  t0=$(date +%s.%N)
  cp --no-preserve=mode "$src" "$R/$name.part" && mv -f "$R/$name.part" "$R/$name"
  rc=$?
  t1=$(date +%s.%N)
  s=$(sha256sum "$R/$name" | cut -d' ' -f1)
  b=$(stat -c %s "$R/$name")
  st=MISMATCH; [ "$s" = "$exp" ] && st=MATCH
  [ $st = MISMATCH ] && { mv -f "$R/$name" "$R/$name.MISMATCH"; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$cell" "$b" "$s" "$exp" "$st" "$src" \
    "$(awk -v a="$t0" -v c="$t1" 'BEGIN{printf "%.1f", c-a}')" | tee -a "$L/stage_release.tsv"
}
stage c001-n-classes.cnd     001   0a623eaeab8e77c115e35177e0fe154b5fd244758b29ac8999ba72ea28d954cd <work>/release_build/c001.cnd
stage c101-dmin-cnd1.cnd1    101   d6607faa1f0196c041b608d41812720b921cce90dfd2b0d94e509fc6a73b0acb <repo>/artifacts/new_table.cnd1
stage c011-pwt-full.pwt2     011   dd69b283217bb41c96b26751591d5bd04f5fb261d29abbbc4fb10313d666adad <work>/pwt2/new_table_full.s3.pwt2
stage c111-pwt-dmin.pwt2     111   0b387083f508818bb42dda65b3e0bc1bdc7161a009788f8036346f475247e317 <repo>/artifacts/new_table.pwt2
stage c111x-pwt-dstar.pwt2   111x  3bc8369bd280892da4517dcdc8573a409081fb937a7bd568abaa13f6ab8cfed6 <work>/pwt2/new_table_opt5.s3.pwt2
stage new_table.orc1         oracle a970d7d665e3ce6e11b5017371833a2ae9fafc12ef5509f0cba12ce102a5815f <work>/data/new_table.orc1
stage c111x-dstar-masks.bin  111x  e599e76119b7bfcad9e981d67c98ec34a5fcd6e2437825874fc077bc50dd08fc <work>/scratch/dsel_set2/full5/masks_full5.bin
# SHA256SUMS from the copies themselves (only files that matched)
{
  for f in $(awk -F'\t' '$6=="MATCH"{print $1}' "$L/stage_release.tsv" | sort); do
    sha256sum "$f" | sed 's/ \*/  /'
  done
} > SHA256SUMS
{
  printf 'file\tcell\tbytes\tsha256\tcopied from\n'
  awk -F'\t' '$6=="MATCH"{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$7}' "$L/stage_release.tsv" | sort
} > MANIFEST.tsv
sha256sum -c SHA256SUMS
echo "stage_release done"
