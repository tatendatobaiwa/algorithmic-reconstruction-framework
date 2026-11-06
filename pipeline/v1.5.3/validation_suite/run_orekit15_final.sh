#!/bin/bash
# Cells 0-5 of the FINAL release file with real Orekit on the saved elements
# (15 runs, 1.5.2 strip rule); each refined-events file is compared byte for
# byte with the earlier 1.5.2 run (runs/orekit15/<run>_v152), the input of the
# exact replays (all runs but ckgr_open, replayed from the Colab output).
PIPE=/home/claude/arf/arf_v152.py
cd /home/claude/validate
mkdir -p /home/claude/runs/orekit15_final logs/final
: > logs/final/orekit15_final.txt
sha256sum $PIPE >> logs/final/orekit15_final.txt
for r in ckgr_open ckgr_openPO ckgr_openPC ckgr_PO ckgr_PC lz_open lz_openPO lz_openPC lz_PO lz_PC karoo_open karoo_openPO karoo_openPC karoo_PO karoo_PC; do
  timeout 3000 python3 real_orekit_c05.py $PIPE $r /home/claude/runs/orekit15_final/$r 0 0 > logs/final/orekit15_final_$r.log 2>&1
  rc=$?
  a=$(sha256sum /home/claude/runs/orekit15_final/$r/arf_refined_events.csv | cut -c1-64)
  b=$(sha256sum /home/claude/runs/orekit15/${r}_v152/arf_refined_events.csv | cut -c1-64)
  [ "$a" = "$b" ] && same=identical || same=DIFFERENT
  echo "$r exit $rc refined events $same to the earlier 1.5.2 run (orekit15/<run>_v152)" >> logs/final/orekit15_final.txt
done
echo DONE >> logs/final/orekit15_final.txt
