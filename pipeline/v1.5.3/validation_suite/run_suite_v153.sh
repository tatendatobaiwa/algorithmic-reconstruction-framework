#!/bin/bash
# The whole validation suite on release 1.5.3 (the 1.5.2 suite, unchanged tests,
# plus test_presentation_v153.py; the confirmatory-analysis test reads the
# release of record from the 1.5.3 protocol).
P=/home/claude/arf/arf_v153.py
S=/home/claude/arf/ship153
cd /home/claude/validate
L=logs/v153/suite
mkdir -p $L
sha256sum $P > $L/pipeline.sha256
run() { name=$1; shift; t0=$(date +%s); timeout 3000 python3 "$@" > $L/$name.log 2>&1; rc=$?; echo "$name rc=$rc $(( $(date +%s) - t0 ))s $(grep -c '^PASS' $L/$name.log) pass $(grep -c '^FAIL' $L/$name.log) fail" >> $L/summary.txt; }
: > $L/summary.txt
run test_cell0 test_cell0.py $P
run test_cell0_outdir test_cell0_outdir.py $P
run test_cell1 test_cell1.py $P
run test_cell2 test_cell2.py $P
run test_cell3 test_cell3.py $P
run test_cell45 test_cell45.py $P /home/claude/arf/arf_v150.py
run test_cell5_exclusions test_cell5_exclusions.py $P
run test_cell6_7c test_cell6_7c.py $P /home/claude/arf/arf_v140.py
run test_cell7d test_cell7d.py $P
run test_cells7to9 test_cells7to9.py $P
run test_presentation_v153 test_presentation_v153.py $P
run test_confirmatory_analysis test_confirmatory_analysis.py $S/arf_confirmatory_analysis.py $S/arf_confirmatory_protocol.json
echo SUITE_DONE >> $L/summary.txt
