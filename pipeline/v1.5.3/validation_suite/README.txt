ARF v1.5.3 validation suite (18 October 2026)
=============================================

What this is
  The scripts and logs behind release 1.5.3 (algorithmic_reconstruction_framework.py,
  SHA-256 14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950), which changes
  figures and report labels only. The 1.5.2 audit's scripts are included, updated where noted;
  its logs are in arf_validation_suite_v152.zip.

How 1.5.3 was made
  make_v153.py              builds 1.5.3 from the shipped 1.5.2 file by exact, counted text
                            replacements (it refuses any other input file)
  ship_cells.py             splits the release into its 14 cells (header + cells == file)

Evidence that no computed value changed
  compare_v152_v153.py      Cells 6-8 of 1.5.2 and 1.5.3 replayed on the same saved inputs (15
                            runs of 12 October through replay150.py, 5 Colab runs of 15 October
                            through arf_colab_parity.replay); every CSV, JSON value, report line
                            and figure compared. --recheck re-applies the comparison to recorded
                            output folders: wall-clock columns (ilp_seconds) are excluded and
                            listed, and Figures 1-4, 8 join the redrawn set (whole-number k ticks,
                            Figure 8 width)
  logs/v153/compare_recheck.jsonl, arf_v152_v153_identity.csv
                            20 of 20 runs identical
  run_v153_real.sh          Cells 0-5 of 1.5.3 with real Orekit on the five 15 October runs
  real_orekit_run.py        (row-by-row comparison with Colab: real_runs_15Oct/*_comparison.json,
                            135,629 rows identical) and full runs that draw every figure
  render_real.py            (real_runs_15Oct/*_render.log: figure sizes)

Tests
  run_suite_v153.sh         the 1.5.2 suite on 1.5.3 plus test_presentation_v153.py:
                            logs/v153/suite/summary.txt, 109 of 109 checks pass
  test_presentation_v153.py Figures 3, 5, 7, 9 (axes, widths, labels), the superseded labels, the
                            per-run statuses and the pooled label (12 checks)
  test_confirmatory_analysis.py
                            now reads the release of record from the protocol file (passes on
                            the 1.5.2 and the 1.5.3 protocol)

Parity reference values and the command-line check
  gen_parity153.py          arf_parity_expected_v153.json (five 15 October runs; each value equal
                            to what the Colab run recorded) and arf_parity_expected_v153_12oct.json
                            (15 runs of 12 October; equal to the 1.5.2 reference)
  logs/v153/parity_cli/     arf_colab_parity.py run as in Colab: 5 of 5 MATCH on a copy of the
                            15 October folders; 15 of 15 MATCH on ARF_12Oct2026_saved_inputs.zip
                            unzipped into a separate folder (--drive)

Protocol
  make_protocol_amendment_a2.py
                            writes amendment A2 from the 1.5.2 protocol, the 1.5.3 file and the
                            comparison (it refuses unless every run is identical)

Not included (by design)
  - WDPA boundary files: the WDPA licence forbids redistribution.
  - Colab outputs (refined events, tile lists, element caches): they are on your Drive. Scripts
    that need them expect /home/claude/realdata/ (12 October) and /home/claude/drive152/
    (15 October); edit the paths at the top of a script to point elsewhere.
