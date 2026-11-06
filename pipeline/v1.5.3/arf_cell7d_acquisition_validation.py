# =============================================================================
# CELL 7D: VALIDATION AGAINST REAL ACQUISITIONS (Sentinel-2 and Landsat archives)
# =============================================================================
# Do the passes ARF predicts correspond to images the satellites actually took?
# For every open-data satellite in the fleet (Sentinel-2, Landsat 8/9), the
# predicted overpasses of the ROI are matched against the scenes in the public
# Earth Search STAC catalogue (Element 84; https://earth-search.aws.element84.com/v1)
# for the same window: same platform, overpass times within
# VALIDATION_MATCH_TOL_MIN. Reported: precision (predicted overpasses that were
# really acquired), recall (real acquisitions that were predicted), how much of
# each matched pass's predicted ROI coverage lies inside the real scenes, and
# the real scene cloud cover (eo:cloud_cover, %) as an observed counterpart of
# the cloud model. Item fields used (checked on live items, 2026-10-09):
# properties.datetime, properties.platform ("sentinel-2b", "landsat-9"),
# properties["eo:cloud_cover"], and the item geometry.
#
# Only the part of the window older than each platform's archive lag is
# evaluated (scenes take time to reach the archive; v1.5.2: per platform,
# VALIDATION_ARCHIVE_LAG_H_BY_PLATFORM, with precision and recall reported per
# platform as well as pooled). PlanetScope has no open archive, so its passes
# are not validated here. Needs internet access; if the catalogue cannot be
# reached the cell records "skipped" and the run goes on.
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 7D is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["7D"] = _ARF_CELL_VERSION

import re
import json
import urllib.request
from datetime import datetime, timedelta, timezone
import numpy as np
import pandas as pd
import h3

logger = logging.getLogger("arf.cell7d")
_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()

VALIDATION_STAC_URL = "https://earth-search.aws.element84.com/v1"
# Both Sentinel-2 L2A collections are queried; the same acquisition appearing
# in both is counted once (grouped into one overpass by platform and time).
VALIDATION_COLLECTIONS = {"Sentinel-2": ["sentinel-2-c1-l2a", "sentinel-2-l2a"], "Landsat 8/9": ["landsat-c2-l2"]}
VALIDATION_MATCH_TOL_MIN = 15.0      # same overpass if times differ by at most this
VALIDATION_GROUP_GAP_MIN = 10.0      # scenes/minutes closer than this belong to one overpass
VALIDATION_ARCHIVE_LAG_H = 48.0      # ignore the most recent part of the window (Sentinel-2; default)
# v1.5.2: per platform. The Landsat collection queried is Collection 2 Level-2.
# USGS (https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-
# science-products, accessed 2026-10-14): Level-2 products "are typically
# available within 24-72 hours" after a scene is processed into Tier 1 or 2;
# "For Landsat 8 OLI/TIRS, it takes approximately 14-16 days to process to a
# Tier 1 or Tier 2 product"; Landsat 9 takes "approximately 3 days to process
# to a Level-2 product". (USGS's Level-1 page gives 3 to 10 days for Landsat 8
# tiering and its generation timeline 4 to 11 days for Level-2; the longest
# statement is used, so a pass is never counted as missed while it may still be
# in processing.) With the single 48 h lag, recent Landsat passes that were
# acquired but not yet published counted as false predictions. Lags: Landsat 9,
# 3 days + 1 day margin; Landsat 8, 16 + 3 days. Sentinel-2 keeps 48 h (the
# Sentinel-2 passes of the October 2026 runs were in the Earth Search catalogue
# within that time).
VALIDATION_ARCHIVE_LAG_H_BY_PLATFORM = {"landsat-9": 96.0, "landsat-8": 456.0}


def _archive_lag_h(platform):
    return float(VALIDATION_ARCHIVE_LAG_H_BY_PLATFORM.get(platform, VALIDATION_ARCHIVE_LAG_H))


def _stac_platform(sat):
    """STAC platform string for a fleet satellite, or None if not validated."""
    name, mission = str(sat.get("id", "")).upper(), str(sat.get("mission", ""))
    m = re.match(r"^SENTINEL[- ]?2([A-Z])", name)
    if mission.startswith("Sentinel-2") and m:
        return f"sentinel-2{m.group(1).lower()}"
    m = re.match(r"^LANDSAT[- ]?([89])", name)
    if mission.startswith("Landsat") and m:
        return f"landsat-{m.group(1)}"
    return None


def stac_search(collections, bbox, start, end, url=VALIDATION_STAC_URL, limit=100, max_pages=100):
    """All items of `collections` intersecting bbox [W, S, E, N] with datetime in
    [start, end] (UTC datetimes), following the STAC API "next" links."""
    body = {"collections": list(collections), "bbox": list(bbox), "limit": limit,
            "datetime": f"{start:%Y-%m-%dT%H:%M:%SZ}/{end:%Y-%m-%dT%H:%M:%SZ}"}
    href, method, items = url.rstrip("/") + "/search", "POST", []
    for _ in range(max_pages):
        data = json.dumps(body).encode() if method == "POST" else None
        req = urllib.request.Request(href, data=data, method=method,
                                     headers={"Content-Type": "application/json", "User-Agent": "ARF-Research-Pipeline/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            page = json.loads(resp.read().decode("utf-8"))
        items.extend(page.get("features", []))
        nxt = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt:
            break
        href, method = nxt["href"], str(nxt.get("method", "GET")).upper()
        if method == "POST":
            body = nxt.get("body", body) if nxt.get("merge") is not True else {**body, **nxt.get("body", {})}
    return items


def _item_tiles(geometry, roi_tiles, res):
    """ROI tiles inside a GeoJSON Polygon / MultiPolygon footprint."""
    polys = ([geometry["coordinates"]] if geometry.get("type") == "Polygon"
             else geometry.get("coordinates", []) if geometry.get("type") == "MultiPolygon" else [])
    cells = set()
    for poly in polys:
        outer = [(lat, lon) for lon, lat in poly[0]]
        holes = [[(lat, lon) for lon, lat in ring] for ring in poly[1:]]
        try:
            cells |= set(h3.polygon_to_cells(h3.LatLngPoly(outer, *holes), res))
        except Exception:
            continue
    return cells & roi_tiles


def _group_overpasses(records, gap_min):
    """records: (platform, datetime, tiles, cloud or None) -> overpasses (one per
    platform and time cluster): platform, time (median), tiles (union), cloud
    (mean of the scenes' eo:cloud_cover), n_scenes."""
    out = []
    for plat in sorted({r[0] for r in records}):
        rs = sorted((r for r in records if r[0] == plat), key=lambda r: r[1])
        cluster = []
        for r in rs + [None]:
            if cluster and (r is None or (r[1] - cluster[-1][1]).total_seconds() / 60.0 > gap_min):
                times = sorted(c[1] for c in cluster)
                clouds = [c[3] for c in cluster if c[3] is not None]
                out.append({"platform": plat, "time": times[len(times) // 2],
                            "tiles": set().union(*(c[2] for c in cluster)),
                            "cloud_pct": float(np.mean(clouds)) if clouds else None, "n_scenes": len(cluster)})
                cluster = []
            if r is not None:
                cluster.append(r)
    return out


def validate_acquisitions(refined_df, fleet, sim_start, window_h, roi_tiles, res,
                          search_fn=stac_search, platform_fn=_stac_platform, now=None):
    """Predicted vs real overpasses for the open-data satellites. Returns (summary, rows)."""
    now = now or datetime.now(timezone.utc).replace(tzinfo=None)
    plat_of = {s["id"]: platform_fn(s) for s in fleet}
    plat_of = {k: v for k, v in plat_of.items() if v}
    w_end = sim_start + timedelta(hours=window_h)
    # v1.5.2: each platform is evaluated up to (now - its archive lag)
    t_end_of = {p: min(w_end, now - timedelta(hours=_archive_lag_h(p))) for p in set(plat_of.values())}
    t_end = max(t_end_of.values()) if t_end_of else min(w_end, now - timedelta(hours=VALIDATION_ARCHIVE_LAG_H))
    base = {"stac_url": VALIDATION_STAC_URL, "match_tolerance_min": VALIDATION_MATCH_TOL_MIN,
            "archive_lag_h": VALIDATION_ARCHIVE_LAG_H, "window_start_utc": f"{sim_start:%Y-%m-%dT%H:%M:%S}",
            "archive_lag_h_by_platform": {p: _archive_lag_h(p) for p in sorted(t_end_of)},
            "evaluated_until_utc_by_platform": {p: f"{t:%Y-%m-%dT%H:%M:%S}" for p, t in sorted(t_end_of.items())},
            "evaluated_until_utc": f"{t_end:%Y-%m-%dT%H:%M:%S}", "platforms": sorted(set(plat_of.values()))}
    if not plat_of:
        return {**base, "status": "skipped: no Sentinel-2 or Landsat 8/9 satellite in this fleet"}, []
    live = {p for p, t in t_end_of.items() if t > sim_start}       # platforms with any archive-ready time
    if not live:
        _lag = min(_archive_lag_h(p) for p in t_end_of)
        return {**base, "status": f"skipped: window not yet older than {_lag:g} h "
                                  f"(re-run Cell 7D after {sim_start + timedelta(hours=_lag):%Y-%m-%d %H:%M} UTC)"}, []

    # Predicted overpasses: refined events of the validated satellites, grouped by time.
    pred_recs = []
    r = refined_df[refined_df["sat"].astype(str).isin(plat_of)]
    for (sat, t_min), g in r.groupby(["sat", "t_min"], observed=True):
        t = sim_start + timedelta(minutes=float(t_min))
        p = plat_of[str(sat)]
        if p in live and t <= t_end_of[p]:
            pred_recs.append((p, t, set(g["tile"].astype(str)), None))
    predicted = _group_overpasses(pred_recs, VALIDATION_GROUP_GAP_MIN)

    # Real overpasses: catalogue scenes touching the ROI tiles.
    lat_lon = np.array([h3.cell_to_latlng(t) for t in roi_tiles])
    bbox = [float(lat_lon[:, 1].min()) - 0.1, float(lat_lon[:, 0].min()) - 0.1,
            float(lat_lon[:, 1].max()) + 0.1, float(lat_lon[:, 0].max()) + 0.1]
    want = set(plat_of.values())
    real_recs, n_items = [], 0
    try:
        for mission, cols in VALIDATION_COLLECTIONS.items():
            if not any(p.startswith("sentinel-2" if mission == "Sentinel-2" else "landsat") for p in want):
                continue
            for it in search_fn(cols, bbox, sim_start, t_end):
                n_items += 1
                p = it.get("properties", {})
                plat = str(p.get("platform", "")).lower()
                if plat not in want:
                    continue
                tiles = _item_tiles(it.get("geometry") or {}, roi_tiles, res)
                if not tiles:
                    continue
                dt = datetime.fromisoformat(str(p["datetime"]).replace("Z", "+00:00")).astimezone(timezone.utc).replace(tzinfo=None)
                if plat not in live or dt > t_end_of[plat]:
                    continue                       # v1.5.2: beyond this platform's archive-ready time
                cc = p.get("eo:cloud_cover")
                real_recs.append((plat, dt, tiles, None if cc is None else float(cc)))
    except Exception as e:                         # no network, catalogue down, ...
        return {**base, "status": f"skipped: catalogue query failed ({type(e).__name__}: {e})"}, []
    real = _group_overpasses(real_recs, VALIDATION_GROUP_GAP_MIN)

    # One-to-one matching, nearest in time first.
    pairs = sorted(((abs((a["time"] - b["time"]).total_seconds()) / 60.0, i, j)
                    for i, a in enumerate(predicted) for j, b in enumerate(real) if a["platform"] == b["platform"]))
    used_p, used_r, match = set(), set(), []
    for dtm, i, j in pairs:
        if dtm <= VALIDATION_MATCH_TOL_MIN and i not in used_p and j not in used_r:
            used_p.add(i); used_r.add(j); match.append((i, j, dtm))
    rows = []
    for i, j, dtm in match:
        a, b = predicted[i], real[j]
        ov = len(a["tiles"] & b["tiles"])
        rows.append({"kind": "matched", "platform": a["platform"], "predicted_utc": f"{a['time']:%Y-%m-%dT%H:%M}",
                     "actual_utc": f"{b['time']:%Y-%m-%dT%H:%M:%S}", "dt_min": round(dtm, 2),
                     "predicted_tiles": len(a["tiles"]), "actual_tiles": len(b["tiles"]), "overlap_tiles": ov,
                     "cloud_pct": b["cloud_pct"], "scenes": b["n_scenes"]})
    for i, a in enumerate(predicted):
        if i not in used_p:
            rows.append({"kind": "predicted_only", "platform": a["platform"], "predicted_utc": f"{a['time']:%Y-%m-%dT%H:%M}",
                         "actual_utc": None, "dt_min": None, "predicted_tiles": len(a["tiles"]), "actual_tiles": None,
                         "overlap_tiles": None, "cloud_pct": None, "scenes": 0})
    for j, b in enumerate(real):
        if j not in used_r:
            rows.append({"kind": "actual_only", "platform": b["platform"], "predicted_utc": None,
                         "actual_utc": f"{b['time']:%Y-%m-%dT%H:%M:%S}", "dt_min": None, "predicted_tiles": None,
                         "actual_tiles": len(b["tiles"]), "overlap_tiles": None, "cloud_pct": b["cloud_pct"],
                         "scenes": b["n_scenes"]})
    m = [x for x in rows if x["kind"] == "matched"]
    pt = sum(x["predicted_tiles"] for x in m)
    at = sum(x["actual_tiles"] for x in m)
    ov = sum(x["overlap_tiles"] for x in m)
    clouds = [b["cloud_pct"] for b in real if b["cloud_pct"] is not None]
    # v1.5.2: per platform -- a satellite that is tracked but not acquiring
    # (e.g. Sentinel-2A on its extension-campaign plan) shows up as predicted
    # overpasses with no match, which a pooled precision hides.
    by_platform = {}
    for plat in sorted(t_end_of):
        n_p = sum(1 for a in predicted if a["platform"] == plat)
        n_r = sum(1 for b in real if b["platform"] == plat)
        n_m = sum(1 for x in m if x["platform"] == plat)
        by_platform[plat] = {"evaluated": plat in live, "archive_lag_h": _archive_lag_h(plat),
                             "predicted_overpasses": n_p, "actual_overpasses": n_r, "matched": n_m,
                             "precision": (n_m / n_p) if n_p else None, "recall": (n_m / n_r) if n_r else None,
                             "evaluated_until_utc": f"{t_end_of[plat]:%Y-%m-%dT%H:%M:%S}" if plat in live else None}
    summary = {**base, "status": "ok", "catalogue_items_returned": n_items, "by_platform": by_platform,
               "predicted_overpasses": len(predicted), "actual_overpasses": len(real), "matched": len(m),
               "precision": (len(m) / len(predicted)) if predicted else None,
               "recall": (len(m) / len(real)) if real else None,
               "median_time_offset_min": float(np.median([x["dt_min"] for x in m])) if m else None,
               "predicted_tiles_confirmed_pct": (100.0 * ov / pt) if pt else None,
               "actual_tiles_predicted_pct": (100.0 * ov / at) if at else None,
               "actual_cloud_cover_mean_pct": float(np.mean(clouds)) if clouds else None,
               "actual_cloud_cover_median_pct": float(np.median(clouds)) if clouds else None,
               "model_cloud_probability_pct": 100.0 * CLOUD_PROBABILITY,
               "note": "PlanetScope passes are not validated (no open archive)."}
    return summary, rows


_val_cols = ["kind", "platform", "predicted_utc", "actual_utc", "dt_min", "predicted_tiles", "actual_tiles",
             "overlap_tiles", "cloud_pct", "scenes"]
if globals().get("FLEET", "starlink") == "starlink":
    acquisition_validation = {"status": "skipped: the Starlink proxy has no imaging archive"}
    _val_rows = []
else:
    _sim_start = datetime.fromisoformat(str(_sim_epoch_datetime64())[:26])
    acquisition_validation, _val_rows = validate_acquisitions(
        refined_events_df, constellation, _sim_start, WINDOW_HOURS, set(target_tiles), H3_RESOLUTION)
pd.DataFrame(_val_rows, columns=_val_cols).to_csv(os.path.join(_out_dir, "arf_acquisition_validation.csv"), index=False)
with open(os.path.join(_out_dir, "arf_acquisition_validation_summary.json"), "w") as f:
    json.dump(acquisition_validation, f, indent=2, default=str)
scenario_results.setdefault(_scenario_key(), {})["acquisition_validation"] = acquisition_validation
logger.info(f"Acquisition validation: {acquisition_validation['status']}"
            + ("" if acquisition_validation["status"] != "ok" else
               f" | precision {acquisition_validation['precision']}, recall {acquisition_validation['recall']}"))
print("\n=== CELL 7D COMPLETE ===")

