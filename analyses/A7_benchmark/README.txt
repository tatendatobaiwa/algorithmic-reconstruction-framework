ARF amendment A7: benchmark metrics (POST HOC, read-only; no runs, no replay)

Reads the files already written by the 56 real-fleet runs of release 1.5.3 (and the 16
A5 wet-season runs) and reports, per arm and window:
  - the physical floor t* in hours (delivery time after the onset to 95% of the tiles
    reachable in the 24-h horizon), the reachable share of the region, and the share of
    onsets at which all passes reach 95% of the WHOLE region;
  - images ECC, coverage-first greedy, time-aware greedy and FCFS each need for that
    level, and the paired savings;
  - the cloud floor per realisation;
  - window-block bootstrap 95% intervals (10,000 resamples, seed 20261030).
  Revision 30 Oct: A5 earthenv runs are skipped (their clear-sky outputs equal the
  p = 0.35 run of the same window), so no window counts twice. Re-installing over the
  first A7 is supported (the earlier protocol is kept as arf_confirmatory_protocol.a7_prev.json).

Install (Colab)
  1. Upload ARF_benchmark_A7.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a7.py (mounts Drive if needed). It checks
     the checksums, release 1.5.3, arf_confirmatory_analysis.py and the A6 protocol;
     installs the A7 protocol (A6 kept as arf_confirmatory_protocol.a6.json) and the script.
  3. !python /content/drive/MyDrive/ARF_code/arf_benchmark_a7.py      (a few minutes)
  Output: MyDrive/ARF_batch/benchmark_a7/. Send the console log and a7_summary.json.
  Do NOT run any earlier install or launch cell again.
