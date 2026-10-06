#!/bin/bash
# Rebuild run (2026-10-06), steps 3-5, one job at a time. Inputs are read-only certified files;
# every output and scratch file goes under <work>/release_build.
export PATH=/usr/bin:$PATH
W=<work>/release_build
B=<repo>/build
DATA=<work>/data
ORC=$DATA/new_table.orc1
J=<work>/release_build/run_job.sh
cd <work>/release_build || exit 2
mkdir -p sc010 sc011 sc110 sc111 sc111x dsel_scratch

step=${1:-all}

if [ "$step" = all ] || [ "$step" = lists ]; then
# ---- step 3: flat cells from n.raw / d.raw (list_build, single-threaded)
$J c001 $W/c001.cnd 0a623eaeab8e77c115e35177e0fe154b5fd244758b29ac8999ba72ea28d954cd -- \
  $B/list_build.exe classes $DATA/n.raw $W/c001.cnd --mod 10810800 --targets-self
$J c000 $W/c000.vbg 5778a94fc44cfe9136ead0eb4d8563745597aa277861dc79bc82ff1c83af07f8 -- \
  $B/list_build.exe vbyte $DATA/n.raw $W/c000.vbg
$J c100 $W/c100.vbg 57e5be328b64de9a00ea0bbf03de95a4141f68bf5d70f67ba9aab12cfd417575 -- \
  $B/list_build.exe vbyte $DATA/d.raw $W/c100.vbg
fi

if [ "$step" = all ] || [ "$step" = trees ]; then
# ---- step 4: byte-identity gate, current pwt2_encode, logged flags, 4 threads
$J c010 $W/c010.pwt2n b7e2e68993e5e09b6c36a642eb01f9629551b5f86e219dbabe644e2fc43ab9e6 -- \
  $B/pwt2_encode.exe $DATA/fpaths.bin $ORC $W/c010.pwt2n --fields 14 --full --naive --scratch $W/sc010 --chunk 33554432 --threads 4
$J c011 $W/c011.pwt2 dd69b283217bb41c96b26751591d5bd04f5fb261d29abbbc4fb10313d666adad -- \
  $B/pwt2_encode.exe $DATA/fpaths.bin $ORC $W/c011.pwt2 --fields 14 --full --ctx-set 3 --scratch $W/sc011 --chunk 33554432 --threads 4
$J c110 $W/c110.pwt2n d2ba4ade81f4e972eadd573bef11e309a7999f9122789be56dbfdbed6a54a1bf -- \
  $B/pwt2_encode.exe $DATA/dpaths.bin $ORC $W/c110.pwt2n --fields 12 --naive --scratch $W/sc110 --chunk 33554432 --threads 4
$J c111 $W/c111.pwt2 0b387083f508818bb42dda65b3e0bc1bdc7161a009788f8036346f475247e317 -- \
  $B/pwt2_encode.exe $DATA/dpaths.bin $ORC $W/c111.pwt2 --fields 12 --ctx-set 3 --scratch $W/sc111 --chunk 33554432 --threads 4
fi

if [ "$step" = all ] || [ "$step" = dstar ]; then
# ---- step 5: d* masks route (search skipped with --sweeps 0), logged settings
$J dsel_emit $W/s2_from_masks.bin 92302de4e0a0d3d7ab4ef10dab431e4fd21805607166ea87b7fccfdb57fc8e8e -- \
  $B/dsel_set2.exe --orc1 $ORC --dpaths $DATA/dpaths.bin --threads 6 --sweeps 0 --min-gain 0.01 --check-rebuild \
  --max-overflow 60000000 --chunk 33554432 --scratch $W/dsel_scratch \
  --masks-in <work>/scratch/dsel_set2/full5/masks_full5.bin --emit $W/s2_from_masks.bin
if grep -q "^dsel_emit	0	.*	MATCH	" logs/summary.tsv; then
  $J c111x $W/c111x.pwt2 3bc8369bd280892da4517dcdc8573a409081fb937a7bd568abaa13f6ab8cfed6 -- \
    $B/pwt2_encode.exe $W/s2_from_masks.bin $ORC $W/c111x.pwt2 --fields 14 --ctx-set 3 --scratch $W/sc111x --chunk 33554432 --threads 4
else
  echo "dsel_emit did not MATCH 92302de4...; d* encode skipped" | tee -a logs/summary_notes.txt
fi
fi
echo "chain $step done"
