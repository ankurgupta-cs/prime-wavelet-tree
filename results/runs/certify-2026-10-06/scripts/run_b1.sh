#!/bin/bash
# Certification run (2026-10-06): certify the flat cells with list_decode, rebuild d*
# from its own decoded paths (pwt2_decode --emit-paths -> pwt2_encode), and
# certify the hosted tree files with the current decoder. One job at a time.
# Inputs are read-only; every output and scratch file goes under
# <work>/release_build.
export PATH=/usr/bin:$PATH
W=<work>/release_build
R=<work>/release_data
B=<repo>/build
ORC=<work>/data/new_table.orc1
J=<work>/release_build/run_job_b1.sh
cd <work>/release_build || exit 2
step=${1:-all}
if [ "$step" = all ] || [ "$step" = lists ]; then
$J b1_c001 - - -- $B/list_decode.exe $R/c001-n-classes.cnd --kind n --certify $ORC
$J b1_c000 - - -- $B/list_decode.exe $W/c000.vbg --kind n --certify $ORC
$J b1_c100 - - -- $B/list_decode.exe $W/c100.vbg --kind d --threads 13 --certify $ORC
$J b1_c101 - - -- $B/list_decode.exe $R/c101-dmin-cnd1.cnd1 --kind d --threads 13 --certify $ORC
fi
if [ "$step" = all ] || [ "$step" = dstar ]; then
$J b1_dstar_emit $W/dstar_paths.bin - -- $B/pwt2_decode.exe <work>/pwt2/new_table_opt5.s3.pwt2 --emit-paths $W/dstar_paths.bin --fields 14
$J b1_dstar_encode $W/dstar_reenc.pwt2 3bc8369bd280892da4517dcdc8573a409081fb937a7bd568abaa13f6ab8cfed6 -- \
  $B/pwt2_encode.exe $W/dstar_paths.bin $ORC $W/dstar_reenc.pwt2 --fields 14 --ctx-set 3 --scratch $W/sc_dstar_b1 --chunk 33554432 --threads 4
fi
if [ "$step" = all ] || [ "$step" = certify ]; then
$J b1_cert_c011 - - -- $B/pwt2_decode.exe $R/c011-pwt-full.pwt2 --certify $ORC
$J b1_cert_c111 - - -- $B/pwt2_decode.exe $R/c111-pwt-dmin.pwt2 --certify $ORC
$J b1_cert_dstar - - -- $B/pwt2_decode.exe $R/c111x-pwt-dstar.pwt2 --certify $ORC
fi
echo "chain $step done"
