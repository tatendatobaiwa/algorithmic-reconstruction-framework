ARF amendment A10: data volume with each tile's own area (corrects A6)
(POST HOC; a read-only replay of the 8 primary windows, no new pipeline runs)

Why. A6 gave every usable tile the global mean area of an H3 resolution-7 cell
(5.16 km2). Over CKGR+Okavango the cells average 5.96 km2 (12,231 cells, 72,853 km2
by h3.cell_area), so every A6 volume in GB was about 15% too low. A10 sums each
image's usable tiles by their own h3.cell_area. Everything else in A6 is unchanged
(rates, sensor mapping, FCFS comparison). A10 supersedes the A6 results.

Install (Colab)
  1. Upload ARF_datavolume_A10.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a10.py. It checks the checksums, release
     1.5.3, arf_colab_parity.py, arf_exact_economy_a4.py and the A9 protocol; installs the
     A10 protocol (A9 kept as arf_confirmatory_protocol.a9.json) and the script.
  3. !pip install -q h3==4.5.0
  4. !python /content/drive/MyDrive/ARF_code/arf_data_volume_a10.py
     About 30 min (a replay of about 4 min per window).
  Output: MyDrive/ARF_batch/data_volume_a10/ (a10_onsets.csv, a10_windows.csv,
  a10_summary.json, a10_table.tex). Send the console log and a10_summary.json.
  Do NOT run any earlier install or launch cell again.
