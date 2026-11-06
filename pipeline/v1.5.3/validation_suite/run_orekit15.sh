#!/bin/bash
# Real-Orekit reproduction of Cells 0-5 for the 15 runs of 12 Oct 2026:
# legacy end caps (must reproduce Colab's refined events) and the v1.5.2 rule.
PIPE=/home/claude/arf/arf_v152.py
for r in ckgr_open ckgr_openPO ckgr_openPC ckgr_PO ckgr_PC lz_open lz_openPO lz_openPC lz_PO lz_PC karoo_open karoo_openPO karoo_openPC karoo_PO karoo_PC; do
  for mode in legacy v152; do
    L=0; [ "$mode" = legacy ] && L=1
    timeout 3000 python3 real_orekit_c05.py $PIPE $r /home/claude/runs/orekit15/${r}_${mode} $L 0 > logs/orekit15_${r}_${mode}.log 2>&1
    echo "$r $mode exit $?" >> logs/orekit15_progress.log
  done
done
echo DONE >> logs/orekit15_progress.log
