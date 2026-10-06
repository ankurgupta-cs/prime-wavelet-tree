#!/bin/bash
# Step 2: hash-verify the intermediates and the oracle copy in <work> (read-only).
export PATH=/usr/bin:$PATH
L=<work>/release_build/logs
D=<work>/data
out=$L/hash_inputs.tsv
: > "$out"
check() { # file expected
  t0=$(date +%s.%N)
  s=$(sha256sum "$1" | cut -d' ' -f1)
  t1=$(date +%s.%N)
  b=$(stat -c %s "$1")
  st=MISMATCH; [ "$s" = "$2" ] && st=MATCH
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$b" "$s" "$2" "$st" "$(awk -v a="$t0" -v c="$t1" 'BEGIN{printf "%.1f", c-a}')" | tee -a "$out"
}
check $D/n.raw         ae272474e7bb0459034f867c21b4f11d8bd2f96b45d2b7b7f0bf69c5653e2e62
check $D/d.raw         22b9b6ab918bca087418dbe67fbbb1f49b4766cf6a4298a293a57d4d7412dcc8
check $D/new_table.orc1 a970d7d665e3ce6e11b5017371833a2ae9fafc12ef5509f0cba12ce102a5815f
check $D/dpaths.bin    391ed7caa7315833c9144c34ff86dcb9f45a0143b270806833c480c654a3dd39
check $D/fpaths.bin    86828c5e8d78a661394acc72533fbded080b2e107ce9defc0f31f364e4b13136
check <work>/scratch/dsel_set2/full5/masks_full5.bin e599e76119b7bfcad9e981d67c98ec34a5fcd6e2437825874fc077bc50dd08fc
check <work>/scratch/dsel_set2/full5/s2_full5.bin    92302de4e0a0d3d7ab4ef10dab431e4fd21805607166ea87b7fccfdb57fc8e8e
echo "hash_inputs done"
