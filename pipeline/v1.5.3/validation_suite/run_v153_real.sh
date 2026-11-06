#!/bin/bash
# Release 1.5.3 on the five complete Colab runs of 15 October 2026 (real Orekit):
# Cells 0-5 compared row by row with Colab's refined events, then all figures drawn.
P=/home/claude/arf/arf_v153.py
D=/home/claude/drive152
E=2026-10-15T08:36:44.671104
O=/home/claude/eval1015/v153
cd /home/claude/eval1015
: > $O/progress.txt
for spec in "openPC earth-resources+planet conservative" "openPO earth-resources+planet optimistic" \
            "PC_v152 planet conservative" "PO_v152 planet optimistic" "open earth-resources optimistic"; do
  set -- $spec; r=$1; fl=$2; case=$3
  timeout 1500 python3 real_orekit_run.py --pipe $P --run $r --out $O/c05_$r --gp-dir $D/$r/ARF_cache \
    --tiles $D/$r/arf_mosaic_tiles.csv --colab-refined $D/$r/arf_refined_events.csv --roi ckgr_okavango \
    --fleet $fl --case $case --ref-epoch $E > $O/c05_$r.log 2>&1
  echo "$r cells0-5 exit $?" >> $O/progress.txt
  timeout 1700 python3 render_real.py $P $D/$r $O/render_$r ckgr_okavango $fl $case $E > $O/render_$r.log 2>&1
  echo "$r render exit $?" >> $O/progress.txt
done
echo DONE >> $O/progress.txt
