"""Cell 7D (validation against the Sentinel-2 / Landsat archives), offline:
the catalogue search is replaced by a mock that returns STAC items built here.
Usage: test_cell7d.py PIPELINE
A. Matching: a predicted overpass with a real scene of the same platform
   within 15 min is matched; one without is predicted-only; a real one with no
   prediction is actual-only; pooled precision/recall follow.
B. Per-platform archive lags (v1.5.2): with 'now' 6 days after the start of
   a 4-day window, Sentinel-2 passes are evaluated over the whole window (48 h
   lag), Landsat 9 up to day 2 (96 h), Landsat 8 not at all (456 h: USGS
   gives up to 14-16 days to Tier 1/2 plus 24-72 h to Level-2), so its acquired-but-unpublished
   pass is neither a false prediction nor a miss.
C. Per-platform precision/recall: a tracked satellite that does not acquire
   (Sentinel-2A on its extension-campaign plan) shows precision 0 for its
   platform while the pooled precision stays high.
D. Platform names, grouping of scenes within 10 min into one overpass, the
   'window not yet old enough' message."""
import sys, os, ast, json, logging
from datetime import datetime, timedelta, timezone
import numpy as np, pandas as pd, h3
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


C = e2e15.cells_of(open(P, encoding="utf-8").read())
tree = ast.parse(C["7D"])
keep = []
for n in tree.body:
    if isinstance(n, ast.FunctionDef):
        keep.append(n)
    elif isinstance(n, ast.Assign) and all(isinstance(t, ast.Name) and t.id.startswith("VALIDATION_") for t in n.targets):
        keep.append(n)
tree.body = keep
ns = {"np": np, "pd": pd, "h3": h3, "re": __import__("re"), "json": json, "datetime": datetime, "timedelta": timedelta,
      "timezone": timezone, "urllib": __import__("urllib.request"), "CLOUD_PROBABILITY": 0.35,
      "logger": logging.getLogger("t7d")}
exec(compile(tree, "c7d", "exec"), ns)

tiles = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
roi = set(tiles)
ll = np.array([h3.cell_to_latlng(t) for t in tiles])
lat0, lat1, lon0, lon1 = ll[:, 0].min(), ll[:, 0].max(), ll[:, 1].min(), ll[:, 1].max()
poly = {"type": "Polygon", "coordinates": [[(lon0 - 0.2, lat0 - 0.2), (lon1 + 0.2, lat0 - 0.2), (lon1 + 0.2, lat1 + 0.2),
                                            (lon0 - 0.2, lat1 + 0.2), (lon0 - 0.2, lat0 - 0.2)]]}
start = datetime(2026, 9, 1, 0, 0, 0)
fleet = [{"id": "SENTINEL-2A", "mission": "Sentinel-2"}, {"id": "SENTINEL-2B", "mission": "Sentinel-2"},
         {"id": "LANDSAT 8", "mission": "Landsat 8/9"}, {"id": "LANDSAT 9", "mission": "Landsat 8/9"},
         {"id": "FLOCK 4X-1", "mission": "PlanetScope"}]
check([ns["_stac_platform"](s) for s in fleet] == ["sentinel-2a", "sentinel-2b", "landsat-8", "landsat-9", None],
      "platform names: sentinel-2a/2b, landsat-8/9; PlanetScope not validated")
# predicted passes (minutes from start): S2B day 1 08:30, S2A day 1 08:45 (not acquired), L9 day 1 08:20,
# L8 day 2 08:25 (too recent for its 288 h lag)
pred = [("SENTINEL-2B", 24 * 60 + 510), ("SENTINEL-2A", 24 * 60 + 525), ("LANDSAT 9", 24 * 60 + 500),
        ("LANDSAT 8", 2 * 24 * 60 + 505)]
rows = [{"sat": s, "t_min": t, "tile": tl} for s, t in pred for tl in tiles[:50]]
refined = pd.DataFrame(rows)
refined["sat"] = refined["sat"].astype("category")


def item(platform, when, cloud):
    return {"properties": {"platform": platform, "datetime": when.strftime("%Y-%m-%dT%H:%M:%S.000Z"),
                           "eo:cloud_cover": cloud}, "geometry": poly}


real = [item("sentinel-2b", start + timedelta(minutes=24 * 60 + 512), 10.0),
        item("sentinel-2b", start + timedelta(minutes=24 * 60 + 513), 20.0),     # same overpass (2 scenes)
        item("landsat-9", start + timedelta(minutes=24 * 60 + 499), 5.0),
        item("sentinel-2b", start + timedelta(minutes=3 * 24 * 60 + 510), 30.0),  # real, not predicted (actual-only)
        item("landsat-8", start + timedelta(minutes=2 * 24 * 60 + 506), 0.0)]     # within L8's 456 h lag: not evaluated
calls = []


def search(cols, bbox, s, e):
    calls.append((tuple(cols), s, e))
    return [it for it in real if s <= datetime.fromisoformat(it["properties"]["datetime"].replace("Z", "")) <= e]


now = start + timedelta(days=6)                   # 6 days after the start of a 4-day window
summ, out_rows = ns["validate_acquisitions"](refined, fleet, start, 96, roi, 7, search_fn=search, now=now)
bp = summ["by_platform"]
check(summ["status"] == "ok" and summ["matched"] == 2 and summ["predicted_overpasses"] == 3 and summ["actual_overpasses"] == 3,
      f"pooled: {summ['matched']} matched of {summ['predicted_overpasses']} predicted and {summ['actual_overpasses']} real "
      f"(two S2B scenes one minute apart are one overpass)")
check(bp["landsat-8"]["evaluated"] is False and bp["landsat-8"]["predicted_overpasses"] == 0
      and bp["landsat-8"]["actual_overpasses"] == 0 and summ["archive_lag_h_by_platform"]["landsat-8"] == 456.0
      and summ["archive_lag_h_by_platform"]["landsat-9"] == 96.0 and summ["archive_lag_h_by_platform"]["sentinel-2a"] == 48.0
      and bp["landsat-9"]["evaluated_until_utc"] == "2026-09-03T00:00:00",
      "per-platform lags: Landsat 8 (456 h) not evaluated 6 days on, so its acquired-but-unpublished pass is "
      "neither a false prediction nor a miss; Landsat 9 (96 h) evaluated to day 2; Sentinel-2 48 h")
check(bp["sentinel-2a"]["precision"] == 0.0 and bp["sentinel-2a"]["predicted_overpasses"] == 1
      and bp["sentinel-2b"]["precision"] == 1.0 and abs(bp["sentinel-2b"]["recall"] - 0.5) < 1e-12
      and bp["landsat-9"]["precision"] == 1.0,
      f"per platform: S2A precision {bp['sentinel-2a']['precision']} (tracked, not acquiring), S2B precision "
      f"{bp['sentinel-2b']['precision']} recall {bp['sentinel-2b']['recall']}, L9 precision {bp['landsat-9']['precision']}")
kinds = sorted(r["kind"] for r in out_rows)
check(kinds == ["actual_only", "matched", "matched", "predicted_only"],
      f"rows: {kinds}")
summ2, _ = ns["validate_acquisitions"](refined, fleet, start, 96, roi, 7, search_fn=search, now=start + timedelta(hours=30))
check(summ2["status"].startswith("skipped: window not yet older than 48 h"),
      f"30 h after the start: {summ2['status'][:60]}...")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell7d.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
