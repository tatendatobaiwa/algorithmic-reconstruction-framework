Figure scripts for the manuscript (matplotlib 3.11, numpy, pandas, scipy, h3 4.5.0). Output: PNG, 600 dpi, RGB.
  figstyle.py               shared style (Okabe-Ito palette, serif, MDPI 600 dpi); save() checks dpi and mode
  fig_reachability_map.py   RUN_DIR OUT.png  -> Figure "reach" (needs arf_refined_events.csv, arf_onset_robustness.csv,
                            arf_mosaic_tiles.csv of one run); asserts per-onset reachable counts equal the run's
  fig_results.py            OUT_DIR -> fig_decisions, fig_images, fig_budget, fig_failsignal, fig_cloud; every figure
                            asserts its inputs against the stored summaries (lower bounds against
                            arf_confirmatory_results.json, verdict counts, pooled cloud, correlations)
Tested on 2 Nov 2026: reachability map on the 3 Aug 2026 primary run and on the 12 Oct 2026 run; all checks passed.
Colab: copy this folder to MyDrive and run e.g.
  !python fig_reachability_map.py /content/drive/MyDrive/ARF_Outputs_..._from20260810 fig_reach_0810.png --window "10 August 2026"
