ARF PAPER 1 PACKAGE, PROTOCOL AMENDMENT A4 (26 October 2026)
==============================================================

What this is
------------
Post hoc evidence for the Earth-observation manuscript (Paper 1), added after
the confirmatory results were read and registered as such in amendment A4 of
arf_confirmatory_protocol.json. Release 1.5.3 (the release of record) is not
changed; no confirmatory decision can change.

Files
-----
arf_colab_install_a4.py          install cell (run once; checks everything)
arf_confirmatory_protocol.json   the A3 protocol with amendment A4 appended
arf_paper1_evidence.py           E1-E5, read-only, from the confirmatory runs
arf_exact_economy_a4.py          exact minimum images at every onset (H1b ties)
arf_sensitivity_runner_a4.py     48 post hoc runs: k = 5/10/15/20, EarthEnv cloud,
                                 full-window exact optimum
arf_sensitivity_summary_a4.py    tables of those runs
p1_positioning_and_framing.tex   positioning table, H1a/H1b framing, data
                                 statement (drafts; citations checked)
arf_paper1_tests.zip             the tests and their logs
SHA256SUMS.txt

Colab, in this order
--------------------
1. Upload ARF_paper1_A4.zip to MyDrive (do not unzip it there).
2. Cell: paste arf_colab_install_a4.py (the one in THIS zip) and run it; it
   replaces an earlier A4 install and keeps it as arf_confirmatory_protocol.a4_prev.json.
3. !python /content/drive/MyDrive/ARF_code/arf_paper1_evidence.py
4. !python /content/drive/MyDrive/ARF_code/arf_exact_economy_a4.py
   (about 10 min per window, 1.5 h; replays Cells 6-8 of each primary run and
   stops if the replay differs from the run by one byte)
5. !python /content/drive/MyDrive/ARF_code/arf_sensitivity_runner_a4.py
   (48 runs, about 3.5 to 4 h; resumable; no Space-Track secrets needed, the
   element sets come from ARF_batch/shared_cache; run it with !python, not
   runpy, so the notebook kernel is not affected; the earthenv runs download the
   EarthEnv monthly cloud files for August and September once, into
   MyDrive/ARF_climatology)
6. !python /content/drive/MyDrive/ARF_code/arf_sensitivity_summary_a4.py
Send me the printed output of steps 3, 4 and 6.

Safety
------
- The sensitivity runs write ONLY to MyDrive/ARF_sensitivity_A4/<variant>/<window>.
  Release 1.5.3 names its output folder from ROI, fleet and window only, so a
  changed budget run would otherwise land in the confirmatory folder; the runner
  prevents this and stops if any confirmatory file changes.
- Do not run the 1.5.3 or 1.5.4 launch cells again: they would put an older
  protocol back over A4. The A3 protocol is kept as
  ARF_code/arf_confirmatory_protocol.a3.json.
- Every A4 number is post hoc and is reported beside, never instead of, the
  confirmatory decisions.

Tested here (sandbox, 26 October 2026)
--------------------------------------
test_evidence.py           8/8   real Cell 6-8 outputs + the real Cell 7D files of window 1
test_exact_economy_v2.py   9/9   real saved inputs of 12 October (CKGR+Okavango), "run" made
                                 with Cell 2's exact altitude and period; negative control: the
                                 parity replay differs in decay_tau_min (4 rows, rel 7.6e-16),
                                 the Colab failure of 26 October; the revised replay is
                                 byte-identical; 19 exact solves, all optimal. No tie occurs in
                                 that data: 14 September is the tie branch's first real test.
test_sensitivity_runner.py 13/13 real arf_batch_runner.py; Cell 0's folder code verbatim
test_sensitivity_summary.py 7/7  real replays at k = 30, k = 10 and the full-window ILP
                                 (the ILP run's 7 result files equal the k = 30 run's)
test_install_a4.py         12/12 real release files, including the upgrade from the first package
