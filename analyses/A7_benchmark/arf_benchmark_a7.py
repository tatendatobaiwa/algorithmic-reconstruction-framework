# =============================================================================
# ARF BENCHMARK METRICS (protocol amendment A7, POST HOC; read-only)
# =============================================================================
# Answers four reviewer questions from the saved outputs of every real-fleet
# run of release 1.5.3 (the 56 confirmatory runs: primary design, fleet
# subsets, PlanetScope swath cases, Sentinel-2A exclusion, Lower Zambezi and
# Karoo) and, when present, the 16 wet-season runs of amendment A5. No
# pipeline run, no replay: only the files the runs wrote.
#   1. The physical floor in hours. t* = delivery time, after the onset, at
#      which all passes together first cover 95% of the tiles reachable in the
#      24-h horizon (best_possible_clear_t95r_min of Cell 7C), with the share of
#      the region reachable at all and the time to 95% and 99% of the WHOLE
#      region where reached. Times are delivery times (pass time plus the
#      modelled downlink wait), not product availability.
#   2. Images each selector needs to reach that level: ECC against
#      coverage-first greedy, time-aware greedy and first-come-first-served
#      (all from Cell 7C), with the time each reaches it.
#   3. Window-clustered uncertainty: onsets share orbits and cloud within a
#      window, so every interval is a window-block bootstrap (windows resampled
#      with replacement, 10,000 resamples, fixed seed); the window is the unit.
#      Runs with another cloud source (A5 earthenv) are skipped: their passes and
#      clear-sky outputs equal those of the p = 0.35 run of the same window.
#   4. Under the realised cloud field: the cloud floor and ECC's gap to it per
#      cloud realisation (Cell 8, arf_cloud_robustness.csv).
# Writes ARF_batch/benchmark_a7/ (a7_onsets.csv, a7_windows.csv,
# a7_summary.json, tab_a7_floor.tex, tab_a7_images.tex).
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_benchmark_a7.py
# =============================================================================
import os
import json
import runpy
import argparse

import numpy as np
import pandas as pd

LEVEL = 95.0
B = 10000
SEED = 20261030
SELECTORS = (("earliest_completion", "ECC"), ("coverage_first", "Coverage-first greedy"),
             ("time_aware", "Time-aware greedy"), ("fcfs_unbudgeted", "FCFS"))
INDEXES = ("confirmatory_v153_primary_index.csv", "confirmatory_v153_secondary_index.csv",
           "confirmatory_v153_exploratory_index.csv", "wetseason_v153_a5_index.csv")
ARM_ORDER = [  # (label, roi, fleet, swath, exclude, season)
    ("Primary: CKGR+Okavango, all 96, conservative", "ckgr_okavango", "earth-resources+planet", "conservative", "", "dry"),
    ("PlanetScope optimistic swath", "ckgr_okavango", "earth-resources+planet", "optimistic", "", "dry"),
    ("Without Sentinel-2A", "ckgr_okavango", "earth-resources+planet", "conservative", "SENTINEL-2A", "dry"),
    ("Earth-resources only (no PlanetScope)", "ckgr_okavango", "earth-resources", "conservative", "", "dry"),
    ("Earth-resources only, without Sentinel-2A", "ckgr_okavango", "earth-resources", "conservative", "SENTINEL-2A", "dry"),
    ("Lower Zambezi", "lower_zambezi", "earth-resources+planet", "conservative", "", "dry"),
    ("Karoo Nature Reserve", "karoo_nature_reserve", "earth-resources+planet", "conservative", "", "dry"),
    ("Primary design, wet season (A5)", "ckgr_okavango", "earth-resources+planet", "conservative", "", "wet"),
]


def boot_ci(values, stat=np.mean, b=B, seed=SEED):
    """Window-block bootstrap: percentile 95% interval of `stat` over windows resampled with replacement."""
    x = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if x.size < 2:
        return None, None
    rng = np.random.default_rng(seed)
    s = np.array([stat(x[rng.integers(0, x.size, x.size)]) for _ in range(b)])
    return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))


def run_dirs(batch):
    out = []
    for name in INDEXES:
        p = os.path.join(batch, name)
        if not os.path.exists(p):
            continue
        idx = pd.read_csv(p, dtype=str).fillna("")
        for _, r in idx[idx.status == "done"].iterrows():
            d = r["out_dir"]
            if os.path.exists(os.path.join(d, "arf_reproducibility_manifest.json")):
                out.append((d, "wet" if name.startswith("wetseason") else "dry", r.get("config", "")))
    seen, uniq = set(), []
    for d, s, c in out:
        if d not in seen:
            seen.add(d); uniq.append((d, s, c))
    return uniq


def onset_rows(d):
    on = pd.read_csv(os.path.join(d, "arf_onset_robustness.csv"))
    n_tiles = len(pd.read_csv(os.path.join(d, "arf_mosaic_tiles.csv"), usecols=["tile"]))
    on = on[on.reachable_clear_tiles > 0]
    r = pd.DataFrame({"onset_shift_min": on.onset_shift_min,
                      "reach_pct_of_region": 100.0 * on.reachable_clear_tiles / n_tiles,
                      "floor_h": on.best_possible_clear_t95r_min / 60.0,
                      "region_t95_h": on.best_possible_clear_t95_min / 60.0,
                      "region_t99_h": on.best_possible_clear_t99_min / 60.0})
    for key, _ in SELECTORS:
        r[f"{key}_images"] = on[f"{key}_clear_images95r"]
        r[f"{key}_gap_min"] = on[f"{key}_clear_t95r_min"] - on.best_possible_clear_t95r_min
    return r, n_tiles


def cloud_rows(d):
    p = os.path.join(d, "arf_cloud_robustness.csv")
    if not os.path.exists(p):
        return None
    c = pd.read_csv(p)
    if not {"floor_t95r_min", "ecc_t95r_min"} <= set(c.columns):
        return None
    return {"cloud_floor_h_median": float((c.floor_t95r_min / 60).median()),
            "cloud_ecc_gap_min_median": float((c.ecc_t95r_min - c.floor_t95r_min).median()),
            "cloud_ecc_images_median": float(c.ecc_images95r.median()) if "ecc_images95r" in c else None,
            "cloud_fcfs_images_median": float(c.fcfs_images95r.median()) if "fcfs_images95r" in c else None}


def arm_of(w):
    for lab, roi, fleet, sw, ex, season in ARM_ORDER:
        if (w.roi, w.fleet, (w.planet_swath_case or "conservative"), w.exclude_satellites or "", w.season) == \
                (roi, fleet, sw, ex, season):
            return lab
    return None


def summarise(wins):
    summ = {}
    for lab, *_ in ARM_ORDER:
        g = wins[wins.arm == lab]
        if g.empty:
            continue
        s = {"windows": int(len(g)), "onsets": int(g.onsets.sum())}
        for col in ("floor_h_median", "reach_pct_median", "region_t95_h_median", "region_t95_reached_pct",
                    "ecc_images_median", "coverage_first_images_median", "time_aware_images_median",
                    "fcfs_unbudgeted_images_median", "ecc_saving_vs_coverage_first_mean",
                    "ecc_saving_vs_fcfs_mean", "ecc_fewer_than_coverage_first_pct", "ecc_gap_min_max",
                    "cloud_floor_h_median"):
            if col in g:
                v = g[col].astype(float)
                lo, hi = boot_ci(v.values)
                s[col] = {"mean_over_windows": float(v.mean()) if v.notna().any() else None,
                          "ci95_window_bootstrap": [lo, hi]}
        summ[lab] = s
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description="Benchmark metrics from the saved real-fleet runs (A7, post hoc).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    a, _ = ap.parse_known_args(argv)
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    an = runpy.run_path(os.path.join(a.code, "arf_confirmatory_analysis.py"), run_name="arf_confirmatory_analysis")
    frames, wins = [], []
    skipped_cloud_variant = 0
    for d, season, _ in run_dirs(a.batch_dir):
        w = an["window_outcomes"](d, int(proto["min_evaluable_onsets"]))
        if w.get("release") != proto["release"] or w.get("versions_consistent") is not True:
            continue
        cfg = json.load(open(os.path.join(d, "arf_reproducibility_manifest.json"))).get("config_overrides") or {}
        if cfg.get("cloud_probability_source", "constant") != "constant":
            skipped_cloud_variant += 1      # same passes and clear-sky outputs as the p = 0.35 run of that window:
            continue                        # counting it again would double the window in the bootstrap
        on, n_tiles = onset_rows(d)
        win = str(w["analysis_epoch"])[:10]
        meta = {"window": win, "season": season, "roi": w["roi"], "fleet": w["fleet"],
                "planet_swath_case": w["planet_swath_case"] or "conservative",
                "exclude_satellites": w["exclude_satellites"] or "", "run_dir": d}
        lab = arm_of(pd.Series(meta))
        if lab is None:
            continue
        on.insert(0, "arm", lab)
        for k, v in meta.items():
            on[k] = v
        frames.append(on)
        row = {"arm": lab, **meta, "region_tiles": n_tiles, "onsets": int(len(on)),
               "floor_h_median": float(on.floor_h.median()),
               "floor_h_p90": float(on.floor_h.quantile(0.9)),
               "reach_pct_median": float(on.reach_pct_of_region.median()),
               "region_t95_reached_pct": float(100.0 * on.region_t95_h.notna().mean()),
               "region_t95_h_median": float(on.region_t95_h.median()) if on.region_t95_h.notna().any() else None,
               "ecc_images_median": float(on.earliest_completion_images.median()),
               "coverage_first_images_median": float(on.coverage_first_images.median()),
               "time_aware_images_median": float(on.time_aware_images.median()),
               "fcfs_unbudgeted_images_median": float(on.fcfs_unbudgeted_images.median()),
               "ecc_saving_vs_coverage_first_mean": float((on.coverage_first_images - on.earliest_completion_images).mean()),
               "ecc_saving_vs_fcfs_mean": float((on.fcfs_unbudgeted_images - on.earliest_completion_images).mean()),
               "ecc_fewer_than_coverage_first_pct": float(100.0 * (on.earliest_completion_images < on.coverage_first_images).mean()),
               "ecc_more_than_coverage_first_onsets": int((on.earliest_completion_images > on.coverage_first_images).sum()),
               "ecc_gap_min_max": float(on.earliest_completion_gap_min.max()),
               "coverage_first_gap_min_median": float(on.coverage_first_gap_min.median()),
               "h1b_pct": w.get("h1b_pct"), "h2_diff": w.get("h2_diff")}
        row.update(cloud_rows(d) or {})
        wins.append(row)
        print(json.dumps({k: row[k] for k in ("arm", "window", "floor_h_median", "reach_pct_median",
                                               "ecc_images_median", "coverage_first_images_median",
                                               "fcfs_unbudgeted_images_median")}), flush=True)
    assert wins, "no run of the release of record found in the indexes"
    dup = pd.DataFrame(wins).duplicated(["arm", "window"])
    assert not dup.any(), "a window appears twice in one arm"
    print(f"cloud-source variants skipped (clear-sky outputs identical to the p = 0.35 run): {skipped_cloud_variant}")
    out = os.path.join(a.batch_dir, "benchmark_a7")
    os.makedirs(out, exist_ok=True)
    all_on = pd.concat(frames, ignore_index=True)
    all_on.to_csv(os.path.join(out, "a7_onsets.csv"), index=False)
    W = pd.DataFrame(wins).sort_values(["arm", "window"])
    W.to_csv(os.path.join(out, "a7_windows.csv"), index=False)
    S = summarise(W)
    for lab in S:     # window-level H1b and H2 with bootstrap intervals, alongside the protocol's t bounds
        g = W[W.arm == lab]
        for col in ("h1b_pct", "h2_diff"):
            v = g[col].astype(float)
            lo, hi = boot_ci(v.values)
            S[lab][col] = {"mean_over_windows": float(v.mean()) if v.notna().any() else None,
                           "ci95_window_bootstrap": [lo, hi]}
    summ = {"label": "POST HOC (amendment A7)", "level_pct_of_reachable": LEVEL, "bootstrap": {"B": B, "seed": SEED,
            "unit": "window"}, "times": "delivery times after the onset (pass time plus modelled downlink wait)",
            "cloud_source_variants_skipped": skipped_cloud_variant, "arms": S}
    json.dump(summ, open(os.path.join(out, "a7_summary.json"), "w"), indent=1, default=str)
    write_tables(W, S, out)
    pd.set_option("display.width", 250)
    print(W.groupby("arm", sort=False)[["floor_h_median", "reach_pct_median", "region_t95_reached_pct",
                                        "ecc_images_median", "coverage_first_images_median",
                                        "fcfs_unbudgeted_images_median", "ecc_gap_min_max"]].mean().round(2).to_string())
    print(f"Files: {out}/a7_onsets.csv, a7_windows.csv, a7_summary.json, tab_a7_floor.tex, tab_a7_images.tex")
    return W, summ


def _ci(s, col, fmt="{:.1f}"):
    v = s.get(col) or {}
    m, (lo, hi) = v.get("mean_over_windows"), (v.get("ci95_window_bootstrap") or [None, None])
    if m is None:
        return "--"
    return fmt.format(m) + ("" if lo is None else " [" + fmt.format(lo) + ", " + fmt.format(hi) + "]")


def write_tables(W, S, out):
    head = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
            "\\setlength{\\tabcolsep}{3.5pt}"]
    L = head + ["\\caption{Post hoc (amendment A7): physical floor $t^{\\ast}$, the delivery time after the onset at "
                "which all passes together cover 95\\% of the tiles reachable in the 24-h horizon, in hours. Values "
                "are means over windows of the per-window median over onsets, with a 95\\% window-block bootstrap "
                "interval (10,000 resamples). Reachable: share of the region seen by any pass in the horizon. "
                "Whole region: share of onsets at which all passes reach 95\\% of the whole region.}",
                "\\label{tab:a7_floor}", "\\begin{tabular}{|l|c|c|c|c|}", "\\hline",
                "Arm & Windows & $t^{\\ast}$ (h) & Reachable (\\%) & Whole region 95\\% (\\% of onsets) \\\\",
                "\\hline"]
    for lab, s in S.items():
        L += [f"{lab} & {s['windows']} & {_ci(s, 'floor_h_median')} & {_ci(s, 'reach_pct_median')} & "
              f"{_ci(s, 'region_t95_reached_pct', '{:.0f}')} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "tab_a7_floor.tex"), "w").write("\n".join(L) + "\n")
    L = head + ["\\caption{Post hoc (amendment A7): images each selector uses to reach 95\\% of the reachable tiles "
                "(clear sky); means over windows of the per-window median, with 95\\% window-block bootstrap "
                "intervals for the paired savings. All selectors reach the level at the floor time unless "
                "stated in the text.}",
                "\\label{tab:a7_images}", "\\begin{tabular}{|l|c|c|c|c|c|c|}", "\\hline",
                "Arm & ECC & Cov.-first & Time-aware & FCFS & Saving vs cov.-first & Saving vs FCFS \\\\",
                "\\hline"]
    for lab, s in S.items():
        L += [f"{lab} & {_ci(s, 'ecc_images_median').split(' [')[0]} & "
              f"{_ci(s, 'coverage_first_images_median').split(' [')[0]} & "
              f"{_ci(s, 'time_aware_images_median').split(' [')[0]} & "
              f"{_ci(s, 'fcfs_unbudgeted_images_median').split(' [')[0]} & "
              f"{_ci(s, 'ecc_saving_vs_coverage_first_mean')} & {_ci(s, 'ecc_saving_vs_fcfs_mean')} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "tab_a7_images.tex"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
