# =============================================================================
# ARF DATA VOLUME (protocol amendment A6, POST HOC; read-only)
# =============================================================================
# Converts the images each selector uses into data volume. For every onset of
# every primary window (clear sky, the H1a level: 95% of reachable tiles), the
# images ECC uses (k = 30) and the images first-come-first-served (FCFS) uses to
# reach the same level are taken from a replay of Cells 6-8 of release 1.5.3
# on the run's own saved inputs (the procedure of arf_exact_economy_a4.py:
# the run's own altitude and period; the replay must reproduce the run's
# arf_onset_robustness.csv, exactly in the columns used here). Each image's
# footprint inside the region is its usable H3 tiles (utility > 0) x 5.16 km2
# (mean area of an H3 resolution-7 cell). Volume = footprint area x the
# distributed-product rate of its sensor (MB per km2), as a range:
#   Sentinel-2 Level-2A       0.066-0.080  (800 MB per 110 x 110 km tile on a
#                             100 km step; Copernicus Collection-1 L2A page)
#   Landsat 8/9 Level-2 SP    0.020-0.035  (file sizes of LSDS-1619 v6.0, Table
#                             A-1, over a 185 x 180 km scene)
#   PlanetScope 4-band        0.78-2.43    (500-700 MB per scene, Planet
#                             PSScene product specification, over the 287.5,
#                             575 and 637 km2 frames of the three generations)
# It is the volume of data covering the region; whole products downloaded are
# larger (tiles and scenes extend beyond it). Sentinel-2 and Landsat acquire on
# fixed plans, so a smaller selection saves the user's download, storage and
# processing, not on-board energy or downlink; this is stated in the output.
# Writes ARF_batch/data_volume_a6/ (per-onset CSV, per-window CSV, JSON, LaTeX).
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_data_volume_a6.py
# =============================================================================
import os
import re
import json
import time
import runpy
import argparse

import numpy as np
import pandas as pd

LEVEL = 95.0
TILE_KM2 = 5.16
RATES = {   # MB per km2 (low, high); sources in the header
    "Sentinel-2": (0.066, 0.080),
    "Landsat": (0.020, 0.035),
    "PlanetScope": (0.78, 2.43),
}
PATTERNS = (("Sentinel-2", r"^SENTINEL[- ]?2[A-Z]?$"), ("Landsat", r"^LANDSAT[- ]?[89]$"),
            ("PlanetScope", r"^FLOCK\b"))


def sensor(name):
    for s, p in PATTERNS:
        if re.match(p, str(name).strip(), flags=re.IGNORECASE):
            return s
    raise ValueError(f"satellite {name!r} has no sensor rate")


def image_rows(ns, sub, keys):
    """Per image: sensor and usable tiles inside the region."""
    util = ns["_effective_util"](sub, None)
    out = []
    for key in keys:
        s = sub.rows(sub.index[key])
        n = int((util[s] > 0.0).sum())
        out.append({"sat": key[0], "t_min": key[1], "sensor": sensor(key[0]), "tiles": n})
    return out


def volume(images):
    lo = sum(i["tiles"] * TILE_KM2 * RATES[i["sensor"]][0] for i in images)
    hi = sum(i["tiles"] * TILE_KM2 * RATES[i["sensor"]][1] for i in images)
    by = {}
    for i in images:
        by.setdefault(i["sensor"], [0, 0.0, 0.0])
        by[i["sensor"]][0] += 1
        by[i["sensor"]][1] += i["tiles"] * TILE_KM2 * RATES[i["sensor"]][0]
        by[i["sensor"]][2] += i["tiles"] * TILE_KM2 * RATES[i["sensor"]][1]
    return lo, hi, by


def onsets(ns):
    H = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
    k = int(ns["MAX_SELECTION_STEPS"])
    rows = []
    for d in range(0, ns["ONSET_SHIFT_HOURS"] * 60 + 1, ns["ONSET_SHIFT_STEP_MIN"]):
        sub, _ = ns["subset_by_time"](ns["event_data"], None, d, d + H)
        R = ns["reachable_tile_count"](sub, None) if len(sub) else 0
        if not R:
            continue
        need = int(np.ceil(LEVEL / 100.0 * R - 1e-9))
        t_e, keys_e, n_e = ns["earliest_completion"](sub, None, k, LEVEL, denominator=R)
        sel_f = ns["select_fcfs"](sub, None, len(sub))[0]
        t_f, n_f, _ = ns["level_stats"](sub, None, sel_f, need)
        if t_e is None or n_f is None:
            continue
        keys_f = sorted(sel_f, key=ns["delivery_sort_key"](sub))[:n_f]    # the images FCFS used to reach the level
        ie, iff = image_rows(ns, sub, keys_e), image_rows(ns, sub, keys_f)
        (el, eh, eb), (fl, fh, fb) = volume(ie), volume(iff)
        row = {"onset_shift_min": d, "ecc_images": n_e, "fcfs_images": n_f,
               "ecc_mb_low": el, "ecc_mb_high": eh, "fcfs_mb_low": fl, "fcfs_mb_high": fh,
               "saving_mb_low": fl - el, "saving_mb_high": fh - eh}
        for s in RATES:
            row[f"ecc_{s}_images"] = eb.get(s, [0])[0]
            row[f"fcfs_{s}_images"] = fb.get(s, [0])[0]
            row[f"saving_{s}_mb_low"] = fb.get(s, [0, 0, 0])[1] - eb.get(s, [0, 0, 0])[1]
            row[f"saving_{s}_mb_high"] = fb.get(s, [0, 0, 0])[2] - eb.get(s, [0, 0, 0])[2]
        rows.append(row)
    return pd.DataFrame(rows)


def gb(x):
    return f"{x / 1000.0:,.1f}" if x is not None and np.isfinite(x) else "--"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Data volume of the images each selector uses (A6, post hoc).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    ap.add_argument("--index", default="confirmatory_v153_primary_index.csv")
    a, _ = ap.parse_known_args(argv)
    ex = runpy.run_path(os.path.join(a.code, "arf_exact_economy_a4.py"), run_name="arf_exact_economy_a4")
    par = ex["load"](os.path.join(a.code, "arf_colab_parity.py"), "arf_colab_parity")
    pipe = os.path.join(a.code, "algorithmic_reconstruction_framework.py")
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    assert par.sha256_file(pipe) == proto["pipeline_sha256"], "ARF_code does not hold the release of record"
    idx = pd.read_csv(os.path.join(a.batch_dir, a.index), dtype=str).fillna("")
    frames, wins = [], []
    for _, r in idx[(idx.status == "done") & (idx.release == proto["release"])].iterrows():
        d = r["out_dir"]
        man = json.load(open(os.path.join(d, "arf_reproducibility_manifest.json")))
        cfg, inp = man.get("config_overrides") or {}, man["inputs"]
        if cfg.get("analysis_epoch") not in proto["windows"]:
            continue
        ov = man["output_validation"]["files"]
        refined, tiles_csv = os.path.join(d, "arf_refined_events.csv"), os.path.join(d, "arf_mosaic_tiles.csv")
        for p in (refined, tiles_csv):
            assert par.sha256_file(p) == ov[os.path.basename(p)]["sha256"], f"{p} differs from what Cell 9 recorded"
        key = list(man["scenario_results"])[0]
        orbit = man["scenario_results"][key]["orbital_period_min"]
        c2 = json.load(open(os.path.join(d, "arf_cell2_summary.json")))
        alt, torb = float(c2["median_alt_km"]), float(c2["t_max_wait_min"])
        assert torb == float(orbit), f"{d}: Cell 2's period differs from the manifest's"
        tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
        win = str(cfg["analysis_epoch"])[:10]
        t0 = time.time()
        ns = ex["replay_exact"](par, alt, torb, pipe, refined, tiles, orbit, inp["active_roi"], inp["fleet"],
                                inp["config"].get("planet_swath_case", "conservative"), inp["n_satellites"])
        ident, ok, rep = ex["compare_onsets"](os.path.join(ns["_out_dir"], "arf_onset_robustness.csv"),
                                              os.path.join(d, "arf_onset_robustness.csv"))
        assert ok, f"{d}: the replay does not reproduce the run's onset results: {rep[:3]}"
        df = onsets(ns)
        on = pd.read_csv(os.path.join(d, "arf_onset_robustness.csv")).set_index("onset_shift_min")
        chk = df.set_index("onset_shift_min")
        assert (chk.ecc_images == on.loc[chk.index, "earliest_completion_clear_images95r"]).all() and \
            (chk.fcfs_images == on.loc[chk.index, "fcfs_unbudgeted_clear_images95r"]).all(), \
            f"{d}: image counts differ from the run's onset results"
        df.insert(0, "window", win)
        frames.append(df)
        med = df.median(numeric_only=True)
        wins.append({"window": win, "onsets": len(df), "replay_identical": ident,
                     "ecc_images_median": med.ecc_images, "fcfs_images_median": med.fcfs_images,
                     **{c: float(med[c]) for c in ("ecc_mb_low", "ecc_mb_high", "fcfs_mb_low", "fcfs_mb_high",
                                                   "saving_mb_low", "saving_mb_high")},
                     "onsets_saving_positive_low": int((df.saving_mb_low > 0).sum()),
                     "onsets_saving_negative_high": int((df.saving_mb_high < 0).sum()),
                     "planetscope_share_of_saving_high": (float(df["saving_PlanetScope_mb_high"].sum() /
                                                                df.saving_mb_high.sum())
                                                          if df.saving_mb_high.sum() else None),
                     "seconds": round(time.time() - t0, 1)})
        print(json.dumps(wins[-1]), flush=True)
    out = os.path.join(a.batch_dir, "data_volume_a6")
    os.makedirs(out, exist_ok=True)
    all_df = pd.concat(frames)
    all_df.to_csv(os.path.join(out, "a6_onsets.csv"), index=False)
    w = pd.DataFrame(wins).sort_values("window")
    w.to_csv(os.path.join(out, "a6_windows.csv"), index=False)
    summ = {"label": "POST HOC (amendment A6)", "rates_mb_per_km2": RATES, "tile_km2": TILE_KM2,
            "scope": ("volume of distributed products covering the region (footprint inside the region); whole "
                      "products are larger; Sentinel-2 and Landsat acquire on fixed plans, so the saving is in the "
                      "user's download, storage and processing, not on-board energy or downlink"),
            "onsets": int(len(all_df)),
            "saving_gb_median_over_onsets": [float(all_df.saving_mb_low.median() / 1000),
                                             float(all_df.saving_mb_high.median() / 1000)],
            "onsets_where_ecc_needs_more_data_high_rate": int((all_df.saving_mb_high < 0).sum()),
            "onsets_where_ecc_needs_more_data_low_rate": int((all_df.saving_mb_low < 0).sum()),
            "windows": wins}
    json.dump(summ, open(os.path.join(out, "a6_summary.json"), "w"), indent=1, default=str)
    L = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         "\\setlength{\\tabcolsep}{3.5pt}",
         "\\caption{Post hoc (amendment A6): data volume of the images each selector uses to reach 95\\% of the "
         "reachable tiles (clear sky), median over onsets, in GB, as a range over the distributed-product rates "
         "(Sentinel-2 L2A 0.066--0.080, Landsat L2SP 0.020--0.035, PlanetScope 4-band 0.78--2.43 MB\\,km$^{-2}$). "
         "Volume covering the region only; whole products are larger.}",
         "\\label{tab:a6_volume}", "\\begin{tabular}{|l|r|r|r|r|r|}", "\\hline",
         "Window & Img.\\ ECC & Img.\\ FCFS & ECC (GB) & FCFS (GB) & Saving (GB) \\\\", "\\hline"]
    for r in w.itertuples():
        L += [f"{r.window} & {r.ecc_images_median:g} & {r.fcfs_images_median:g} & "
              f"{gb(r.ecc_mb_low)}--{gb(r.ecc_mb_high)} & {gb(r.fcfs_mb_low)}--{gb(r.fcfs_mb_high)} & "
              f"{gb(r.saving_mb_low)}--{gb(r.saving_mb_high)} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "a6_table.tex"), "w").write("\n".join(L) + "\n")
    pd.set_option("display.width", 220)
    print(w.drop(columns=["seconds"]).round(1).to_string(index=False))
    print(json.dumps({k: v for k, v in summ.items() if k != "windows"}, indent=1, default=str))
    print(f"Files: {out}/a6_onsets.csv, a6_windows.csv, a6_summary.json, a6_table.tex")
    return all_df, summ


if __name__ == "__main__":
    main()
