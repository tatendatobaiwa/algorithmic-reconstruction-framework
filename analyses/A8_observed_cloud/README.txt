ARF amendment A8: real Sentinel-2 cloud (POST HOC; no pipeline runs)

For the 16 primary-design windows (8 dry: Aug-Sep 2026; 8 wet: Jan-Feb 2026), every
Sentinel-2 L2A scene over CKGR+Okavango in the 48-h window is found in Earth Search
(the release's own stac_search) and its Scene Classification Layer is read at 80 m and
sampled at 7 points per H3 tile. Reported per window: overpasses, the share of tile
observations that are cloud or cloud shadow (against the model's p = 0.35), the share of
overpasses that would raise the Fail-Signal, and the time to 95% Sentinel-2 coverage with
cloud ignored and with real cloud (the delay cloud adds, in hours).
Not measured: PlanetScope (no open per-pixel masks), Landsat (requester-pays assets).

Install (Colab)
  1. Upload ARF_realcloud_A8.zip to MyDrive. Do NOT unzip it there.
  2. New cell: paste and run arf_colab_install_a8.py. It checks the checksums, release
     1.5.3 and the A7 protocol; installs the A8 protocol (A7 kept) and the script.
  3. !pip install -q rasterio h3==4.5.0
  4. !python /content/drive/MyDrive/ARF_code/arf_real_cloud_a8.py
     It reads scenes over the internet (about 20-40 per window; expect tens of minutes).
     Resumable: finished windows are cached in ARF_batch/real_cloud_a8/cache.
  Output: MyDrive/ARF_batch/real_cloud_a8/. Send the full console log, a8_summary.json and
  a8_unreadable.csv. This is the first live test of the code against Earth Search: if a
  window fails or many scenes are unreadable, stop and send the log.
  Do NOT run any earlier install or launch cell again.

Revision 30 Oct: the summary now pools over tile observations (each overpass weighted
by the tiles it observes) and by month, counts "full" passes (>= 5% of the region), and
writes tab_a8_pooled.tex. Re-install over the first A8 is supported; finished windows come
from the cache, so the rerun takes seconds and makes no new catalogue queries.
