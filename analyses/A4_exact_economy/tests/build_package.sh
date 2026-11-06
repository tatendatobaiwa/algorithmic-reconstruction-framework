#!/bin/bash
# Assemble ARF_paper1_A4.zip: the scripts, the A4 protocol (made from the A3 protocol with the
# final script checksums), the drafts, the tests, README and SHA256SUMS.
set -euo pipefail
P=/home/claude/arf/paper1
S=/home/claude/arf/ship_p1
A3=/home/claude/arf/ship154/arf_confirmatory_protocol.json
mkdir -p "$S"
cp "$P"/arf_paper1_evidence.py "$P"/arf_exact_economy_a4.py "$P"/arf_sensitivity_runner_a4.py \
   "$P"/arf_sensitivity_summary_a4.py "$P"/arf_colab_install_a4.py "$S"/
cp "$P"/manuscript/p1_positioning_and_framing.tex "$S"/
python3 "$P"/make_protocol_amendment_a4.py "$A3" "$S"/arf_sensitivity_runner_a4.py "$S"/arf_sensitivity_summary_a4.py \
  "$S"/arf_paper1_evidence.py "$S"/arf_exact_economy_a4.py "$S"/arf_confirmatory_protocol.json
( cd "$P" && rm -f "$S"/arf_paper1_tests.zip && zip -q -r "$S"/arf_paper1_tests.zip tests/*.py tests/*.txt tests/*.log \
    make_protocol_amendment_a4.py build_package.sh fixture_real_w1 )
cd "$S"
rm -f SHA256SUMS.txt
sha256sum arf_paper1_evidence.py arf_exact_economy_a4.py arf_sensitivity_runner_a4.py arf_sensitivity_summary_a4.py \
  arf_colab_install_a4.py arf_confirmatory_protocol.json p1_positioning_and_framing.tex arf_paper1_tests.zip \
  README_paper1_A4.txt > SHA256SUMS.txt
rm -f ARF_paper1_A4.zip
zip -q ARF_paper1_A4.zip arf_paper1_evidence.py arf_exact_economy_a4.py arf_sensitivity_runner_a4.py \
  arf_sensitivity_summary_a4.py arf_colab_install_a4.py arf_confirmatory_protocol.json p1_positioning_and_framing.tex \
  arf_paper1_tests.zip README_paper1_A4.txt SHA256SUMS.txt
echo "built $S/ARF_paper1_A4.zip"; cat SHA256SUMS.txt
