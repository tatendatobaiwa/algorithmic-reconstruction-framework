"""Decompose the physical floor of every onset into the parts that set it.

For an onset at minute d of a window, the floor F (Methods, Eq. floor) is the
delivery time, relative to d, at which ALL passes acquired in [d, d + H) first
cover a share q of the cells they can reach. It splits exactly into

    F = W + G + L
    W  wait:      onset -> acquisition of the first pass with a usable cell
    G  geometry:  that first acquisition -> acquisition of the pass that completes level q
    L  latency:   acquisition -> delivery of the completing pass (downlink wait)

and each onset is labelled with the sun elevation at the region's centroid at
the onset (from the pipeline's own NOAA routine, copied verbatim below), so a
floor set by waiting for daylight is told apart from one set by geometry.

The floor is recomputed from the refined events with the pipeline's rules
(events = (satellite, t_min) passes, usable = utility > 0 with
utility = max(0, 1 - A/theta_max) * D, ordered by delivery time t_min + wait)
and, when the run's arf_onset_robustness.csv is present, checked onset by onset
against the pipeline's best_possible_clear_t95r_min. Any disagreement stops the
script: the decomposition is only reported for floors it reproduces exactly.

Clear sky only (the floor of the hypotheses). Release-independent: reads the
CSV outputs of 1.5.x runs.

Usage: arf_floor_decomposition.py RUN_DIR WINDOW_START_UTC OUT_CSV
         [--level 95] [--horizon-min 1440] [--shift-hours 24] [--step-min 20]
"""
import argparse
import json
import math
import os
import sys

import h3
import numpy as np
import pandas as pd


def solar_elevation_deg(t_utc, lat_deg, lon_deg):
    """Verbatim from the ARF pipeline (Cell 5, releases 1.5.2-1.5.4): NOAA GML
    general solar position equations, geometric elevation, no refraction."""
    t = np.atleast_1d(np.asarray(t_utc, dtype="datetime64[s]"))
    years = t.astype("datetime64[Y]")
    doy = (t.astype("datetime64[D]") - years).astype(int) + 1
    yr = years.astype(int) + 1970
    ndays = np.where((yr % 4 == 0) & ((yr % 100 != 0) | (yr % 400 == 0)), 366, 365)
    sec = (t - t.astype("datetime64[D]")).astype(int)
    hour = sec / 3600.0
    g = 2.0 * np.pi / ndays * (doy - 1 + (hour - 12.0) / 24.0)
    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                       - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = hour * 60.0 + eqtime + 4.0 * np.asarray(lon_deg)
    ha = np.radians(tst / 4.0 - 180.0)
    lat = np.radians(lat_deg)
    cos_zen = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(ha)
    return 90.0 - np.degrees(np.arccos(np.clip(cos_zen, -1.0, 1.0)))


def load_events(run_dir):
    """Usable passes: one row per (sat, t_min) with its delivery time and the
    set of cells with utility > 0 (the pipeline's rule)."""
    ev = pd.read_csv(os.path.join(run_dir, "arf_refined_events.csv"))
    theta = ev["theta_max_deg"].to_numpy(float)
    util = np.clip(1.0 - ev["A_nadir_deg"].to_numpy(float) / theta, 0.0, None) * ev["D_downlink"].to_numpy(float)
    ev = ev.loc[util > 0.0, ["sat", "t_min", "tile", "T_wait_min"]]
    ev = ev[np.isfinite(ev["T_wait_min"].to_numpy(float))]
    tiles = pd.Index(pd.unique(ev["tile"]))
    ev = ev.assign(tid=tiles.get_indexer(ev["tile"]))
    passes = []
    for (sat, t), g in ev.groupby(["sat", "t_min"], sort=True):
        w = g["T_wait_min"].to_numpy(float)
        passes.append((int(t), int(t) + float(w[0]), sat, np.unique(g["tid"].to_numpy())))
    return passes, tiles


def decompose(passes, n_tiles, d, horizon, level):
    """Floor and its split for onset d, or None when nothing is reachable."""
    sub = [p for p in passes if d <= p[0] < d + horizon]
    if not sub:
        return None
    reach = np.unique(np.concatenate([p[3] for p in sub])).size
    need = level / 100.0 * reach
    sub.sort(key=lambda p: p[1])                                # by delivery time
    seen = np.zeros(n_tiles, bool)
    n_seen = 0
    first_acq = min(p[0] for p in sub)
    for acq, dlv, sat, cells in sub:
        new = cells[~seen[cells]]
        if new.size:
            seen[new] = True
            n_seen += new.size
            if 100.0 * n_seen / reach >= level - 1e-9:
                return {"reachable_cells": int(reach), "floor_min": dlv - d,
                        "W_wait_min": first_acq - d, "G_geometry_min": acq - first_acq,
                        "L_latency_min": dlv - acq, "completing_sat": sat,
                        "passes_in_horizon": len(sub)}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("window_start_utc")
    ap.add_argument("out_csv")
    ap.add_argument("--level", type=float, default=95.0)
    ap.add_argument("--horizon-min", type=int, default=1440)
    ap.add_argument("--shift-hours", type=int, default=24)
    ap.add_argument("--step-min", type=int, default=20)
    a = ap.parse_args()

    passes, tiles = load_events(a.run_dir)
    mt = pd.read_csv(os.path.join(a.run_dir, "arf_mosaic_tiles.csv"), usecols=["tile"])
    ll = np.array([h3.cell_to_latlng(t) for t in mt["tile"]])
    # Centroid of the region's cells on the unit sphere (cells are nearly equal in area).
    x = np.cos(np.radians(ll[:, 0])) * np.cos(np.radians(ll[:, 1]))
    y = np.cos(np.radians(ll[:, 0])) * np.sin(np.radians(ll[:, 1]))
    z = np.sin(np.radians(ll[:, 0]))
    clat = math.degrees(math.atan2(z.mean(), math.hypot(x.mean(), y.mean())))
    clon = math.degrees(math.atan2(y.mean(), x.mean()))
    t0 = np.datetime64(a.window_start_utc.replace("Z", ""), "s")

    rows = []
    for d in range(0, a.shift_hours * 60 + 1, a.step_min):
        r = decompose(passes, len(tiles), d, a.horizon_min, a.level)
        onset = t0 + np.timedelta64(d * 60, "s")
        sun = float(solar_elevation_deg(onset, clat, clon)[0])
        lst_h = ((d / 60.0 + clon / 15.0) + pd.Timestamp(str(t0)).hour) % 24.0   # mean solar time, h
        row = {"onset_min": d, "onset_utc": str(onset), "sun_elev_centroid_deg": round(sun, 2),
               "local_mean_solar_time_h": round(lst_h, 2)}
        if r:
            row.update(r)
            assert abs(r["W_wait_min"] + r["G_geometry_min"] + r["L_latency_min"] - r["floor_min"]) < 1e-9
        rows.append(row)
    out = pd.DataFrame(rows)

    check = {"checked": False}
    rob = os.path.join(a.run_dir, "arf_onset_robustness.csv")
    if os.path.exists(rob) and a.level == 95.0:
        ref = pd.read_csv(rob)[["onset_shift_min", "best_possible_clear_t95r_min", "reachable_clear_tiles"]]
        m = out.merge(ref, left_on="onset_min", right_on="onset_shift_min", how="inner")
        bad_f = m[~np.isclose(m["floor_min"], m["best_possible_clear_t95r_min"], equal_nan=True)]
        bad_r = m[m["reachable_cells"] != m["reachable_clear_tiles"]]
        if len(bad_f) or len(bad_r) or len(m) != len(out):
            print("MISMATCH with the pipeline's floors:")
            print(bad_f[["onset_min", "floor_min", "best_possible_clear_t95r_min"]].head(10).to_string())
            print(bad_r[["onset_min", "reachable_cells", "reachable_clear_tiles"]].head(10).to_string())
            sys.exit(1)
        check = {"checked": True, "onsets": int(len(m)), "floors_equal": True, "reachable_equal": True}
    out.to_csv(a.out_csv, index=False)
    ok = out.dropna(subset=["floor_min"])
    night = ok["sun_elev_centroid_deg"] < 0
    summ = {"run_dir": os.path.basename(os.path.normpath(a.run_dir)), "centroid": [round(clat, 3), round(clon, 3)],
            "onsets": int(len(out)), "onsets_with_floor": int(len(ok)), "pipeline_check": check,
            "median_floor_h": round(float(ok["floor_min"].median()) / 60, 2),
            "median_W_h": round(float(ok["W_wait_min"].median()) / 60, 2),
            "median_G_h": round(float(ok["G_geometry_min"].median()) / 60, 2),
            "median_L_h": round(float(ok["L_latency_min"].median()) / 60, 2),
            "median_floor_h_onset_sun_up": (round(float(ok.loc[~night, "floor_min"].median()) / 60, 2)
                                            if (~night).any() else None),
            "median_floor_h_onset_night": (round(float(ok.loc[night, "floor_min"].median()) / 60, 2)
                                           if night.any() else None)}
    print(json.dumps(summ, indent=1))
    json.dump(summ, open(os.path.splitext(a.out_csv)[0] + "_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
