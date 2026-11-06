ARF WET-SEASON ARM, PROTOCOL AMENDMENT A5 (EXPLORATORY), 28 October 2026
=========================================================================
Release 1.5.3 unchanged. Primary design (CKGR + Okavango, open-data + PlanetScope,
conservative swath), eight 48-h windows starting Mondays 5 January to 23 February
2026, run retrospectively; each with the protocol's constant cloud probability
(p035) and with the EarthEnv monthly climatology (earthenv): 16 runs.
Exploratory: no hypothesis decision; reported beside the confirmatory results.

Colab, in this order
1. Upload ARF_wetseason_A5.zip to MyDrive (do not unzip it there); paste
   arf_colab_install_a5.py into a cell and run it.
2. Download MODCF_monthlymean_01.tif and MODCF_monthlymean_02.tif yourself from
   https://www.earthenv.org/cloud into MyDrive/ARF_climatology (names unchanged).
3. With the Space-Track secrets set and the tab open, paste
   arf_colab_prefetch_wet_a5.py into a cell and run it (8 windows, one query each).
4. !pip install -q h3==4.5.0
   !python /content/drive/MyDrive/ARF_code/arf_wetseason_runner_a5.py   (about 1 h)
5. !python /content/drive/MyDrive/ARF_code/arf_wetseason_summary_a5.py
Send me the output of steps 3 and 5.

Safety: runs write only to MyDrive/ARF_wetseason_A5/<variant>/<window>; the
confirmatory and A4 folders are checked after every run and the batch stops at the
first change. Do not run the A4 install cell or any launch cell again.

Tested here (sandbox, 28 October 2026)
test_prefetch_a5.py   6/6   real 1.5.3 fetch functions (AST); real CelesTrak lists of
                            12 October; cache names identical to the real Drive files
test_wet_runner_a5.py 10/10 real batch runner, real A4 runner; Cell 0 folder code verbatim
test_wet_summary_a5.py 7/7  real Cell 6-8 outputs; statistics equal scipy; table compiles
test_install_a5.py    8/8   real release files, real A4 protocol and A4 runner
Stated limit: the fleet is built from the current CelesTrak list, so a satellite
retired since January 2026 drops out (Cell 2 prints the fleet).
