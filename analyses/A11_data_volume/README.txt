ARF amendment A11: data volume with corrected sensor rates (corrects A10 and A6)
(POST HOC; a read-only replay of the 8 primary windows, no new pipeline runs)

Why. Checking the A6 rate sources against their documents found two faults.
  PlanetScope: the upper rate (2.43 MB/km2) used the frame of PS2 (Dove Classic,
  25 x 11.5 km), whose imagery ended in April 2022 (Planet, PlanetScope Product
  Specifications, June 2026). For the current PS2.SD and PSB.SD instruments the
  rate is 500-700 MB per 4-band scene over 575-637 km2: 0.78-1.22 MB/km2.
  Landsat: the 185 x 180 km scene area is not in the source cited. A11 uses the
  USGS EROS scene size (about 170 x 183 km) with LSDS-1619 Table A-1: 0.022-0.038.
  Sentinel-2 is unchanged (0.066-0.080). A11 also writes each selector's usable
  footprint area per sensor, so any later rate change is a recalculation, not a
  replay. Tile areas as A10. A11 supersedes A10 and A6.

Install (Colab)
  1. Upload ARF_datavolume_A11.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a11.py. It checks the checksums, release
     1.5.3, arf_colab_parity.py, arf_exact_economy_a4.py and the A10 protocol; installs the
     A11 protocol (A10 kept as arf_confirmatory_protocol.a10.json) and the script.
  3. !pip install -q h3==4.5.0
  4. !python /content/drive/MyDrive/ARF_code/arf_data_volume_a11.py
     About 30 min (a replay of about 2 min per window).
  Output: MyDrive/ARF_batch/data_volume_a11/ (a11_onsets.csv, a11_windows.csv,
  a11_summary.json, a11_table.tex). Send the console log.
  Do NOT run any earlier install or launch cell again.
