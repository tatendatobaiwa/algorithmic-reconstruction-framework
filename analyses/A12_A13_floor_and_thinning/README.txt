ARF post hoc amendments A12 and A13 (4 November 2026), for Paper 1's evidence base
==================================================================================
A12  Floor decomposition: what sets the earliest possible mosaic (wait for the first usable
     pass, coverage geometry, or delivery). Script: arf_floor_decomposition.py
A13  PlanetScope acquisition fraction: the floor, reachable share and H1a/H1b shares when
     only 75, 50 or 25% of PlanetScope's predicted passes acquire (20 seeds, nested).
     Script: arf_planet_thinning.py
Both are READ-ONLY replays of the 8 dry-season primary runs and the 8 wet-season A5 runs.
Both STOP unless they reproduce the run's own floors (A12) or floors, reachable counts,
ECC times and images and FCFS images (A13) exactly at every onset.

How to run
  1. Upload ARF_posthoc_A12A13.zip to MyDrive. Do not unzip it there.
  2. Paste arf_colab_run_a12_a13.py into one Colab cell and run it once.
  3. Upload the ARF_A12A13_results_<time>.zip it writes to MyDrive.

Tested here (sandbox): both scripts on a copy of the 3 August 2026 primary run (73/73 onsets
reproduced); the Colab cell end to end on a mock Drive. Not tested: your real Drive layout of
the wet-season runs (the cell looks for exactly one run folder under each
ARF_wetseason_A5/p035/<date> and stops if it finds none or several).
