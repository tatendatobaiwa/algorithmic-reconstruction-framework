# =============================================================================
# ARF REAL-CLOUD CHECK, SENTINEL-2 (protocol amendment A8, POST HOC; read-only)
# =============================================================================
# The pipeline's cloud is modelled (p = 0.35 correlated field, or the EarthEnv
# climatology). A8 measures real cloud for the one open part of the fleet whose
# per-pixel cloud masks are freely readable: Sentinel-2. For each primary window
# (8 dry, 8 wet), every Sentinel-2 L2A scene over the region in the 48-h window
# is found in the Earth Search catalogue (the release's own stac_search, Cell
# 7D), and its Scene Classification Layer (SCL, 20 m) is read at reduced
# resolution and sampled at 7 points per H3 resolution-7 tile (centre and the
# midpoints towards the 6 vertices).
#   SCL classes (ESA Level-2A algorithm): 0 no data, 1 saturated/defective,
#   2 topographic/dark shadows, 3 cloud shadows, 4 vegetation, 5 not vegetated,
#   6 water, 7 unclassified, 8 cloud medium, 9 cloud high, 10 thin cirrus,
#   11 snow/ice.
#   A tile is OBSERVED by a scene if >= 4 of its 7 points hold data (class != 0);
#   CLEAR if >= half of its valid points are clear (4, 5, 6, 11); CLOUDY if
#   >= half are cloud or cloud shadow (3, 8, 9, 10). Sensitivity: class 7
#   counted as clear.
# Reported per window: Sentinel-2 overpasses; the observed share of cloudy tile
# observations (to set against the model's p = 0.35 and the EarthEnv month);
# the share of overpasses that would raise the Fail-Signal (fewer than half of
# their observed tiles clear); and the REALISED clear coverage of the region by
# Sentinel-2 over time against the same scenes with cloud ignored, i.e. the
# delay real cloud adds, in hours (acquisition times from the catalogue).
# PlanetScope (no open per-pixel masks without a licence) and Landsat (Earth
# Search links its Collection 2 assets to a requester-pays bucket) are not
# measured; this is stated in the output.
# Writes ARF_batch/real_cloud_a8/ (a8_overpasses.csv, a8_windows.csv,
# a8_summary.json, tab_a8_cloud.tex, tab_a8_pooled.tex). Resumable: finished windows are cached in
# real_cloud_a8/cache and not read again (delete that folder to start over).
#
# Use (Colab):  !pip install -q rasterio h3==4.5.0
#               !python /content/drive/MyDrive/ARF_code/arf_real_cloud_a8.py
# =============================================================================
import os
import ast
import json
import argparse
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

WINDOW_H = 48.0
DECIMATE = 4                     # 20 m SCL read at 80 m (internal overviews)
CLEAR = (4, 5, 6, 11)
CLOUD = (3, 8, 9, 10)
LEVELS = (50.0, 80.0, 95.0)
COLLECTIONS = ("sentinel-2-c1-l2a", "sentinel-2-l2a")      # c1 preferred when a scene is in both
INDEXES = (("confirmatory_v153_primary_index.csv", "dry"), ("wetseason_v153_a5_index.csv", "wet"))
MODEL_P = 0.35
FULL_PASS_PCT = 5.0               # an overpass observing >= 5% of the region's tiles is a "full" pass


def release_stac_search(pipeline):
    """stac_search and VALIDATION_STAC_URL taken from the release file (Cell 7D), unchanged."""
    src = open(pipeline).read()
    tree = ast.parse(src)
    ns = {}
    exec("import json, urllib.request", ns)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "VALIDATION_STAC_URL" for t in node.targets):
            exec(compile(ast.Module([node], []), pipeline, "exec"), ns)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "stac_search":
            exec(compile(ast.Module([node], []), pipeline, "exec"), ns)
    return ns["stac_search"]


def tile_points(tiles):
    """(n_tiles, 7) lon and lat arrays: centre and midpoints towards the 6 vertices."""
    import h3
    lon = np.zeros((len(tiles), 7)); lat = np.zeros((len(tiles), 7))
    for i, t in enumerate(tiles):
        clat, clon = h3.cell_to_latlng(t)
        v = h3.cell_to_boundary(t)[:6]
        lat[i, 0], lon[i, 0] = clat, clon
        for j, (vlat, vlon) in enumerate(v, start=1):
            lat[i, j], lon[i, j] = (clat + vlat) / 2.0, (clon + vlon) / 2.0
    return lon, lat


def sample_scl(href, lon, lat, decimate=DECIMATE):
    """SCL class at each point (0 where outside the raster)."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import transform
    from affine import Affine
    with rasterio.Env(AWS_NO_SIGN_REQUEST="YES", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        with rasterio.open(href) as src:
            h, w = max(1, src.height // decimate), max(1, src.width // decimate)
            arr = src.read(1, out_shape=(h, w), resampling=Resampling.nearest)
            T = src.transform * Affine.scale(src.width / w, src.height / h)
            xs, ys = transform("EPSG:4326", src.crs, lon.ravel().tolist(), lat.ravel().tolist())
    cols, rows = ~T * (np.asarray(xs), np.asarray(ys))
    rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
    ok = (rows >= 0) & (rows < h) & (cols >= 0) & (cols < w)
    out = np.zeros(rows.size, dtype=np.uint8)
    out[ok] = arr[rows[ok], cols[ok]]
    return out.reshape(lon.shape)


def classify(scl, clear_classes=CLEAR):
    """Per tile: observed, clear, cloudy (boolean arrays) from (n_tiles, 7) SCL classes."""
    valid = scl != 0
    nv = valid.sum(1)
    nclear = np.isin(scl, clear_classes).sum(1)
    ncloud = np.isin(scl, CLOUD).sum(1)
    observed = nv >= 4
    clear = observed & (nclear * 2 >= nv)
    cloudy = observed & ~clear & (ncloud * 2 >= nv)
    return observed, clear, cloudy


def scene_key(it):
    p = it.get("properties", {})
    mgrs = p.get("s2:mgrs_tile") or p.get("grid:code") or it.get("id", "")[:-8]
    return (str(p.get("platform", "")).lower(), str(mgrs).replace("MGRS-", ""), str(p.get("datetime", ""))[:16])


def scl_href(it):
    a = it.get("assets", {})
    for k in ("scl", "SCL"):
        if k in a and a[k].get("href"):
            return a[k]["href"]
    return None


def time_to(levels, times, sets, denom):
    """Hours at which the union of `sets` (ordered by `times`) first reaches each level (% of denom)."""
    out, cov, seen = {}, 0, set()
    order = np.argsort(times)
    curve = []
    for i in order:
        seen |= sets[i]
        curve.append((times[i], 100.0 * len(seen) / denom if denom else 0.0))
    for lv in levels:
        out[lv] = next((t for t, c in curve if c >= lv - 1e-9), None)
    return out, (curve[-1][1] if curve else 0.0)


def window_analysis(items, tiles, start, lon, lat, sampler=sample_scl, clear_classes=CLEAR):
    """items: STAC items (both collections). Returns (overpass rows, window row, clear-with-7 row)."""
    best = {}
    for it in items:
        if not str(it.get("properties", {}).get("platform", "")).lower().startswith("sentinel-2"):
            continue
        k = scene_key(it)
        c = it.get("collection", "")
        if k not in best or (c == COLLECTIONS[0] and best[k].get("collection") != COLLECTIONS[0]):
            best[k] = it
    scenes, fails = [], []
    for k, it in sorted(best.items(), key=lambda kv: kv[0][2]):
        dt = datetime.fromisoformat(str(it["properties"]["datetime"]).replace("Z", "+00:00")).astimezone(timezone.utc).replace(tzinfo=None)
        t_h = (dt - start).total_seconds() / 3600.0
        if not (0.0 <= t_h < WINDOW_H):
            continue
        href = scl_href(it)
        if href is None:
            fails.append({"scene": it.get("id"), "error": "no scl asset"}); continue
        try:
            scl = sampler(href, lon, lat)
        except Exception as e:                      # unreadable asset: recorded, not hidden
            fails.append({"scene": it.get("id"), "error": f"{type(e).__name__}: {e}"[:200]}); continue
        scenes.append({"platform": k[0], "time_h": t_h, "id": it.get("id"), "scl": scl,
                       "eo_cloud_cover": it["properties"].get("eo:cloud_cover")})
    res = {}
    for label, cc in (("main", clear_classes), ("with7", tuple(clear_classes) + (7,))):
        ov_rows, obs_sets, clr_sets, times = [], [], [], []
        # one overpass = one platform within 10 min (as Cell 7D)
        groups = []
        for s in sorted(scenes, key=lambda s: (s["platform"], s["time_h"])):
            if groups and groups[-1][0]["platform"] == s["platform"] and s["time_h"] - groups[-1][-1]["time_h"] <= 10 / 60:
                groups[-1].append(s)
            else:
                groups.append([s])
        n_obs = n_cloudy = n_clear = 0
        for g in groups:
            obs = np.zeros(len(tiles), bool); clr = np.zeros(len(tiles), bool); cld = np.zeros(len(tiles), bool)
            for s in g:
                o, c, k_ = classify(s["scl"], cc)
                clr |= c; cld |= k_ & ~clr; obs |= o
            cld &= ~clr
            n_o, n_c, n_k = int(obs.sum()), int(clr.sum()), int(cld.sum())
            n_obs += n_o; n_clear += n_c; n_cloudy += n_k
            t = float(np.median([s["time_h"] for s in g]))
            ov_rows.append({"platform": g[0]["platform"], "time_h": t, "scenes": len(g), "tiles_observed": n_o,
                            "tiles_clear": n_c, "tiles_cloudy": n_k,
                            "cloudy_pct_of_observed": 100.0 * n_k / n_o if n_o else None,
                            "fail_signal": bool(n_o and n_c * 2 < n_o),
                            "eo_cloud_cover_mean": float(np.mean([s["eo_cloud_cover"] for s in g
                                                                  if s["eo_cloud_cover"] is not None]))
                            if any(s["eo_cloud_cover"] is not None for s in g) else None})
            obs_sets.append(set(np.flatnonzero(obs))); clr_sets.append(set(np.flatnonzero(clr))); times.append(t)
        n = len(tiles)
        reach = set().union(*obs_sets) if obs_sets else set()
        t_clr, cov_clr = time_to(LEVELS, times, clr_sets, n)
        t_sky, cov_sky = time_to(LEVELS, times, obs_sets, n)
        t_clr_r, _ = time_to(LEVELS, times, clr_sets, len(reach))
        t_sky_r, _ = time_to(LEVELS, times, obs_sets, len(reach))
        res[label] = {
            "overpasses": len(groups), "scenes": len(scenes), "scenes_unreadable": len(fails),
            "s2_reachable_pct_of_region": 100.0 * len(reach) / n if n else None,
            "cloudy_obs_pct": 100.0 * n_cloudy / n_obs if n_obs else None,
            "clear_obs_pct": 100.0 * n_clear / n_obs if n_obs else None,
            "fail_signal_overpass_pct": 100.0 * np.mean([r["fail_signal"] for r in ov_rows]) if ov_rows else None,
            "clear_cov_48h_pct_of_region": cov_clr, "observed_cov_48h_pct_of_region": cov_sky,
            **{f"t{lv:g}_reach_clear_h": t_clr_r[lv] for lv in LEVELS},
            **{f"t{lv:g}_reach_nocloud_h": t_sky_r[lv] for lv in LEVELS},
            **{f"t{lv:g}_region_clear_h": t_clr[lv] for lv in LEVELS},
            **{f"t{lv:g}_region_nocloud_h": t_sky[lv] for lv in LEVELS},
            "_overpasses": ov_rows, "_fails": fails}
    return res


def windows_to_run(batch):
    out = []
    for name, season in INDEXES:
        p = os.path.join(batch, name)
        if not os.path.exists(p):
            continue
        idx = pd.read_csv(p, dtype=str).fillna("")
        for _, r in idx[idx.status == "done"].iterrows():
            d = r["out_dir"]
            mp = os.path.join(d, "arf_reproducibility_manifest.json")
            if not os.path.exists(mp):
                continue
            cfg = json.load(open(mp)).get("config_overrides") or {}
            if (cfg.get("active_roi"), cfg.get("fleet"), cfg.get("planet_swath_case", "conservative")) != \
                    ("ckgr_okavango", "earth-resources+planet", "conservative") or cfg.get("exclude_satellites") \
                    or cfg.get("cloud_probability_source", "constant") != "constant":
                continue
            out.append({"season": season, "window": str(cfg["analysis_epoch"])[:10],
                        "start": datetime.fromisoformat(str(cfg["analysis_epoch"])[:19]), "run_dir": d})
    df = pd.DataFrame(out).drop_duplicates("window")
    return df.sort_values("window").to_dict("records")


def main(argv=None, search_fn=None, sampler=sample_scl):
    ap = argparse.ArgumentParser(description="Real Sentinel-2 cloud over the primary windows (A8, post hoc).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    a, _ = ap.parse_known_args(argv)
    search = search_fn or release_stac_search(os.path.join(a.code, "algorithmic_reconstruction_framework.py"))
    out = os.path.join(a.batch_dir, "real_cloud_a8")
    os.makedirs(out, exist_ok=True)
    wins, ovs, fails = [], [], []
    cache_dir = os.path.join(out, "cache")
    os.makedirs(cache_dir, exist_ok=True)
    for w in windows_to_run(a.batch_dir):
        cp = os.path.join(cache_dir, f"{w['season']}_{w['window']}.json")
        if os.path.exists(cp):                      # resumable: a finished window is not read again
            c = json.load(open(cp))
            wins.append(c["row"]); ovs += c["overpasses"]; fails += c["fails"]
            print(f"{w['window']}: from cache", flush=True)
            continue
        tiles = list(pd.read_csv(os.path.join(w["run_dir"], "arf_mosaic_tiles.csv"), usecols=["tile"])["tile"].astype(str))
        lon, lat = tile_points(tiles)
        bbox = [float(lon.min()) - 0.1, float(lat.min()) - 0.1, float(lon.max()) + 0.1, float(lat.max()) + 0.1]
        items = search(list(COLLECTIONS), bbox, w["start"], w["start"] + timedelta(hours=WINDOW_H))
        r = window_analysis(items, tiles, w["start"], lon, lat, sampler=sampler)
        m = r["main"]
        for o in m["_overpasses"]:
            ovs.append({"season": w["season"], "window": w["window"], **o})
        fails += [{"window": w["window"], **f} for f in m["_fails"]]
        row = {"season": w["season"], "window": w["window"], "region_tiles": len(tiles),
               **{k: v for k, v in m.items() if not k.startswith("_")},
               "cloudy_obs_pct_with7clear": r["with7"]["cloudy_obs_pct"],
               "t95_reach_clear_h_with7clear": r["with7"]["t95_reach_clear_h"]}
        row["cloud_delay_t95_reach_h"] = (None if row["t95_reach_clear_h"] is None or row["t95_reach_nocloud_h"] is None
                                          else row["t95_reach_clear_h"] - row["t95_reach_nocloud_h"])
        wins.append(row)
        json.dump({"row": row, "overpasses": [o for o in ovs if o["window"] == w["window"]],
                   "fails": [x for x in fails if x["window"] == w["window"]]}, open(cp, "w"), default=str)
        print(json.dumps({k: row[k] for k in ("season", "window", "overpasses", "scenes", "scenes_unreadable",
                                               "s2_reachable_pct_of_region", "cloudy_obs_pct",
                                               "fail_signal_overpass_pct", "t95_reach_nocloud_h",
                                               "t95_reach_clear_h", "clear_cov_48h_pct_of_region")}, default=str), flush=True)
    assert wins, "no primary-design window found in the indexes"
    W = pd.DataFrame(wins)
    W.to_csv(os.path.join(out, "a8_windows.csv"), index=False)
    pd.DataFrame(ovs).to_csv(os.path.join(out, "a8_overpasses.csv"), index=False)
    pd.DataFrame(fails, columns=["window", "scene", "error"]).to_csv(os.path.join(out, "a8_unreadable.csv"), index=False)
    seasons = {}
    for s, g in W.groupby("season"):
        seasons[s] = {"windows": int(len(g)), "overpasses": int(g.overpasses.sum()), "scenes": int(g.scenes.sum()),
                      "scenes_unreadable": int(g.scenes_unreadable.sum()),
                      "cloudy_obs_pct_mean": float(g.cloudy_obs_pct.mean()),
                      "cloudy_obs_pct_range": [float(g.cloudy_obs_pct.min()), float(g.cloudy_obs_pct.max())],
                      "fail_signal_overpass_pct_mean": float(g.fail_signal_overpass_pct.mean()),
                      "s2_reachable_pct_mean": float(g.s2_reachable_pct_of_region.mean()),
                      "clear_cov_48h_pct_mean": float(g.clear_cov_48h_pct_of_region.mean()),
                      "windows_t95_reach_clear": int(g.t95_reach_clear_h.notna().sum()),
                      "cloud_delay_t95_reach_h_median": (float(g.cloud_delay_t95_reach_h.median())
                                                         if g.cloud_delay_t95_reach_h.notna().any() else None)}
    # Pooled over tile observations (each overpass weighted by the tiles it observes): a pass clipping a corner
    # of the region (a few dozen tiles) must not weigh as much as a full pass. "Full" passes observe >= 5% of it.
    O = pd.DataFrame(ovs)
    if not O.empty:
        O["month"] = O.window.str[5:7]
        O["full_pass"] = O.tiles_observed >= FULL_PASS_PCT / 100.0 * W.region_tiles.iloc[0]
        def pooled(g):
            n = int(g.tiles_observed.sum())
            ok = g.eo_cloud_cover_mean.notna()
            return {"overpasses": int(len(g)), "full_passes": int(g.full_pass.sum()), "tile_observations": n,
                    "cloudy_obs_pct": 100.0 * g.tiles_cloudy.sum() / n if n else None,
                    "fail_signal_overpasses": int(g.fail_signal.sum()),
                    "fail_signal_full_passes": int((g.fail_signal & g.full_pass).sum()),
                    "passes_le5pct_cloudy": int((g.cloudy_pct_of_observed <= 5).sum()),
                    "passes_ge50pct_cloudy": int((g.cloudy_pct_of_observed >= 50).sum()),
                    "corr_with_eo_cloud_cover": (float(np.corrcoef(g.cloudy_pct_of_observed[ok], g.eo_cloud_cover_mean[ok])[0, 1])
                                                 if ok.sum() >= 3 else None)}
        for s, g in O.groupby("season"):
            seasons[s]["pooled"] = pooled(g)
            seasons[s]["by_month"] = {m: pooled(h) for m, h in g.groupby("month")}
    summ = {"label": "POST HOC (amendment A8)", "sensor": "Sentinel-2 L2A SCL only",
            "model_cloud_probability_constant": MODEL_P, "clear_classes": CLEAR, "cloud_classes": CLOUD,
            "points_per_tile": 7, "read_resolution_m": 20 * DECIMATE, "times": "acquisition times (catalogue)",
            "not_measured": "PlanetScope (no open per-pixel masks without a licence); Landsat (requester-pays assets)",
            "seasons": seasons}
    json.dump(summ, open(os.path.join(out, "a8_summary.json"), "w"), indent=1, default=str)
    f = lambda x, d=1: "--" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"
    L = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         "\\setlength{\\tabcolsep}{3.5pt}",
         "\\caption{Post hoc (amendment A8): real cloud in the Sentinel-2 scenes of each primary window (L2A Scene "
         "Classification, 7 points per H3 tile). Cloudy: share of tile observations classed cloud or cloud shadow "
         "(model: $p = 0.35$). Fail-Signal: share of overpasses with fewer than half of their observed tiles clear. "
         "$t_{95}$: hours after the window start until Sentinel-2 covers 95\\% of the tiles it observes in the window, "
         "ignoring cloud and with real cloud. PlanetScope and Landsat are not measured.}",
         "\\label{tab:a8_cloud}", "\\begin{tabular}{|l|c|c|c|c|c|c|c|}", "\\hline",
         "Window & Overpasses & Observed (\\%) & Cloudy (\\%) & Fail-Signal (\\%) & $t_{95}$ no cloud (h) & "
         "$t_{95}$ real (h) & Clear at 48 h (\\%) \\\\", "\\hline"]
    for r in W.itertuples():
        L += [f"{r.window} ({r.season}) & {r.overpasses} & {f(r.s2_reachable_pct_of_region)} & {f(r.cloudy_obs_pct)} & "
              f"{f(r.fail_signal_overpass_pct, 0)} & {f(r.t95_reach_nocloud_h)} & {f(r.t95_reach_clear_h)} & "
              f"{f(r.clear_cov_48h_pct_of_region)} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "tab_a8_cloud.tex"), "w").write("\n".join(L) + "\n")
    if not O.empty:
        L = ["\\begin{table}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
             "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Post hoc (amendment A8): real Sentinel-2 cloud, pooled over tile observations (each overpass "
             "weighted by the tiles it observes). Fail: overpasses with fewer than half of their observed tiles clear "
             f"(full passes: $\\geq$ {FULL_PASS_PCT:g}\\% of the region). Bimodal: passes $\\leq$ 5\\% / $\\geq$ 50\\% cloudy.}}",
             "\\label{tab:a8_pooled}", "\\begin{tabular}{|l|c|c|c|c|}", "\\hline",
             "Season, month & Passes (full) & Cloudy (\\%) & Fail (full) & Bimodal \\\\", "\\hline"]
        for s in sorted(seasons):
            rows = [(s, seasons[s]["pooled"])] + [(f"{s}, {m}", v) for m, v in sorted(seasons[s].get("by_month", {}).items())]
            for lab, v in rows:
                L += [f"{lab} & {v['overpasses']} ({v['full_passes']}) & {f(v['cloudy_obs_pct'])} & "
                      f"{v['fail_signal_overpasses']} ({v['fail_signal_full_passes']}) & "
                      f"{v['passes_le5pct_cloudy']} / {v['passes_ge50pct_cloudy']} \\\\", "\\hline"]
        L += ["\\end{tabular}", "\\end{table}"]
        open(os.path.join(out, "tab_a8_pooled.tex"), "w").write("\n".join(L) + "\n")
    print(json.dumps(summ["seasons"], indent=1))
    print(f"Files: {out}/a8_windows.csv, a8_overpasses.csv, a8_unreadable.csv, a8_summary.json, tab_a8_cloud.tex")
    return W, summ


if __name__ == "__main__":
    main()
