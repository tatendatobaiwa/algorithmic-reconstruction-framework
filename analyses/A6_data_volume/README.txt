ARF amendment A6: data volume (POST HOC, read-only; no new pipeline runs)

What it does
  For every onset of the eight confirmatory primary windows, it replays Cells 6-8 of
  release 1.5.3 on the run's own saved inputs (the A4 exact-economy procedure, which
  must reproduce the run's arf_onset_robustness.csv), takes the images ECC (k = 30) and
  first-come-first-served use to reach the H1a level (95% of reachable tiles), and
  converts them to data volume: usable H3 tiles x 5.16 km2 x the sensor's
  distributed-product rate (MB per km2, low and high):
    Sentinel-2 L2A   0.066-0.080   Copernicus Collection-1 Level-2A product page
    Landsat 8/9 L2SP 0.020-0.035   USGS LSDS-1619 v6.0, Table A-1; scene size LSDS-2082
    PlanetScope      0.78-2.43     Planet PSScene product specification (June 2026)
  Sentinel-2 and Landsat acquire on fixed plans, so the saving is the user's download,
  storage and processing, not on-board energy or downlink.

Install (Colab)
  1. Upload ARF_datavolume_A6.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a6.py (it mounts Drive if needed).
     It checks the checksums, release 1.5.3, arf_colab_parity.py, arf_exact_economy_a4.py
     and the A5 protocol; installs the A6 protocol (A5 kept as
     arf_confirmatory_protocol.a5.json) and arf_data_volume_a6.py.
  3. !pip install -q h3==4.5.0
  4. !python /content/drive/MyDrive/ARF_code/arf_data_volume_a6.py      (about 30 min)
  Output: MyDrive/ARF_batch/data_volume_a6/ (a6_onsets.csv, a6_windows.csv,
  a6_summary.json, a6_table.tex). Send the console log and a6_summary.json.
  Do NOT run any earlier install or launch cell again.

Tests (sandbox, real 12 October inputs): tests/test_data_volume_a6.py 5/5,
tests/test_install_a6.py (see SHA256SUMS.txt for the shipped files).
