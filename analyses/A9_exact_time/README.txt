ARF amendment A9: exact time at tight budgets, and the Fail-Signal in hours
(POST HOC; replays of saved inputs only, no new pipeline runs)

Part 1. For k = 5, 10, 15, 20 and onsets every 80 min of the 16 primary windows (8 dry,
8 wet A5 p = 0.35), ECC is compared with the exact earliest time to 95% of reachable
tiles: HiGHS decides whether k images delivered by t can cover the level
(solve_ilp_min_images of release 1.5.3), bisection over delivery times. Verdicts: ECC at
floor / at exact time / delay avoidable / later than exact / missed a cover / no k-cover
exists / unknown (time limit). An unknown is never reported as a result.
Part 2. For the 20 cloud realisations of each window, the Fail-Signal plan and the
open-loop plan of the same size are rebuilt (they must equal the run's
arf_cloud_robustness.csv) and compared in hours: open-loop's time to its own final
coverage minus the Fail-Signal's time to the same coverage; times to 50% and 80%.

Install (Colab)
  1. Upload ARF_exacttime_A9.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a9.py. It checks the checksums, release
     1.5.3, arf_colab_parity.py, arf_exact_economy_a4.py and the A8 protocol; installs the
     A9 protocol (A8 kept) and the script.
  3. !pip install -q h3==4.5.0
  4. !python /content/drive/MyDrive/ARF_code/arf_exact_time_a9.py
     Several hours (a replay of about 5 min per window, then the exact solves and the
     20 realisations). Resumable: rerun the same cell after a disconnection; finished
     windows come from ARF_batch/exact_time_a9/cache. A first look at one window:
       !python .../arf_exact_time_a9.py --only 2026-08-03
     (its cache is reused by the full run).
  Output: MyDrive/ARF_batch/exact_time_a9/. Send the console log and a9_summary.json.
  Do NOT run any earlier install or launch cell again.

Revision 30 Oct: in the Fail-Signal summary, a realisation in which the adaptive plan ends
with LESS coverage than the open-loop plan (so never reaches its level) now counts as the
worst outcome (later), instead of being dropped; the first version dropped these, which
flattered the Fail-Signal. Re-install over the first A9 is supported; all windows come
from the cache, so the rerun takes seconds.
