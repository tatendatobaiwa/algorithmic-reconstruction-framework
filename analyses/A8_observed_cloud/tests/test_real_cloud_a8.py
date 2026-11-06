"""Test of arf_real_cloud_a8.py with the REAL tile list of the CKGR+Okavango region (12,231 H3 res-7 tiles, from
the 12 October run) and SYNTHETIC Sentinel-2 SCL rasters written as tiled GeoTIFFs with overviews in UTM 34S
(EPSG:32734), the CRS of real Sentinel-2 tiles over the region. Earth Search cannot be reached from this workspace,
so the catalogue is replaced by items with the real Earth Search field layout (collection, properties.platform,
properties.datetime, properties["s2:mgrs_tile"], properties["eo:cloud_cover"], assets.scl.href).
Checks:
  - the release's own stac_search is extracted unchanged from release 1.5.3
  - per-tile classes at 80 m read agree with an independent full-resolution point sample (>= 99% of tiles)
  - a scene in both collections counts once (c1 kept); a Landsat item, a scene outside the window and an
    unreadable asset are excluded, the last recorded in a8_unreadable.csv
  - cloudy share, Fail-Signal share and times to 50/80/95% recomputed independently
  - one window per season from the indexes; earthenv duplicate skipped; LaTeX table compiles at IEEE width
Usage: test_real_cloud_a8.py SCRIPT PIPELINE_V153 TILES"""
import sys, os, json, shutil, tempfile, subprocess, importlib.util, io, contextlib, hashlib, ast
from datetime import datetime, timedelta
import numpy as np, pandas as pd
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import transform

S, P, TILES = sys.argv[1:4]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


spec = importlib.util.spec_from_file_location("a8", S); a8 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a8)
ss = a8.release_stac_search(P)
src = open(P).read()
orig = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "stac_search")
check(ss.__code__.co_code == next(c for c in compile(ast.Module([orig], []), P, "exec").co_consts if hasattr(c, "co_code")).co_code
      and ss.__defaults__[0] == "https://earth-search.aws.element84.com/v1",
      "stac_search extracted unchanged from release 1.5.3 (Earth Search v1 URL)")

root = tempfile.mkdtemp(prefix="a8_")
tiles = list(pd.read_csv(TILES, usecols=["tile"])["tile"].astype(str))
lon, lat = a8.tile_points(tiles)


def raster(path, kind):
    """40 m UTM 34S raster over lon 22.6-24.4, lat -21.2..-19.2: data domain on the west 85%, clear west /
    cloud (9) east / shadow (3) strip / cirrus (10) block / class 7 block, laid out in 4 km blocks."""
    x0, y0 = transform("EPSG:4326", "EPSG:32734", [22.6], [-19.2])
    x1, y1 = transform("EPSG:4326", "EPSG:32734", [24.4], [-21.2])
    px = 40.0
    w, h = int((x1[0] - x0[0]) / px), int((y0[0] - y1[0]) / px)
    yy, xx = np.mgrid[0:h, 0:w]
    bx, by = xx // 100, yy // 100                      # 4 km blocks (multiples of the 80 m read)
    a = np.full((h, w), 4, np.uint8)
    if kind == "mixed":
        a[(bx + by) % 3 == 0] = 9
        a[(bx % 5 == 0)] = 3
        a[(by % 7 == 0) & (bx % 2 == 0)] = 10
        a[(by % 11 == 0)] = 7
    elif kind == "cloudy":
        a[:] = 9; a[(bx % 4 == 0)] = 5
    a[:, int(w * 0.85) // 100 * 100:] = 0                # no data east (swath edge)
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=1, dtype="uint8", crs="EPSG:32734",
                       transform=from_origin(x0[0], y0[0], px, px), tiled=True, blockxsize=512, blockysize=512,
                       nodata=0) as d:
        d.write(a, 1)
    with rasterio.open(path, "r+") as d:
        d.build_overviews([2, 4, 8], rasterio.enums.Resampling.nearest)
    return path


r_mixed = raster(os.path.join(root, "mixed.tif"), "mixed")
r_cloudy = raster(os.path.join(root, "cloudy.tif"), "cloudy")
# independent: full-resolution sample at the same points
with rasterio.open(r_mixed) as d:
    xs, ys = transform("EPSG:4326", d.crs, lon.ravel().tolist(), lat.ravel().tolist())
    full = np.array([v[0] for v in d.sample(zip(xs, ys), masked=False)], dtype=np.uint8).reshape(lon.shape)
fast = a8.sample_scl(r_mixed, lon, lat)
o1, c1, k1 = a8.classify(full); o2, c2, k2 = a8.classify(fast)
agree = np.mean((o1 == o2) & (c1 == c2) & (k1 == k2))
check(agree >= 0.99 and o2.sum() > 1000 and c2.sum() > 100 and k2.sum() > 100,
      f"80 m read agrees with a full-resolution sample on {100 * agree:.2f}% of tiles "
      f"(observed {o2.sum()}, clear {c2.sum()}, cloudy {k2.sum()})")

start = datetime(2026, 8, 3)


def item(coll, plat, t_h, mgrs, href, cc=10.0, iid=None):
    return {"id": iid or f"S2X_{mgrs}_{coll}_{t_h}", "collection": coll,
            "properties": {"platform": plat, "datetime": (start + timedelta(hours=t_h)).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
                           "s2:mgrs_tile": mgrs, "eo:cloud_cover": cc},
            "assets": {"scl": {"href": href}}}


ITEMS = [item("sentinel-2-c1-l2a", "sentinel-2b", 8.5, "34KFD", r_mixed),
         item("sentinel-2-l2a", "sentinel-2b", 8.5, "34KFD", r_cloudy),           # same scene, old collection
         item("sentinel-2-c1-l2a", "sentinel-2c", 32.5, "34KFD", r_cloudy),
         item("sentinel-2-c1-l2a", "sentinel-2a", 20.0, "34KGD", os.path.join(root, "missing.tif")),
         item("sentinel-2-c1-l2a", "sentinel-2b", 50.0, "34KFD", r_mixed),         # after the window
         item("landsat-c2-l2", "landsat-9", 10.0, "x", r_mixed)]


def search(cols, bbox, a, b):
    assert list(cols) == list(a8.COLLECTIONS) and len(bbox) == 4
    return ITEMS


# batch with one dry and one wet primary window (+ an earthenv duplicate that must be skipped)
batch = os.path.join(root, "ARF_batch"); os.makedirs(batch)
def run(name, epoch, extra=None):
    d = os.path.join(root, name); os.makedirs(d)
    shutil.copy(TILES, os.path.join(d, "arf_mosaic_tiles.csv"))
    json.dump({"framework_version": "1.5.3", "config_overrides": {"active_roi": "ckgr_okavango",
               "fleet": "earth-resources+planet", "planet_swath_case": "conservative", "analysis_epoch": epoch, **(extra or {})}},
              open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    return {"run_id": name, "status": "done", "release": "1.5.3", "out_dir": d}
pd.DataFrame([run("dry", "2026-08-03T00:00:00")]).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
pd.DataFrame([run("wet", "2026-01-05T00:00:00"), run("wete", "2026-01-05T00:00:00", {"cloud_probability_source": "earthenv"})]
             ).to_csv(os.path.join(batch, "wetseason_v153_a5_index.csv"), index=False)
code = os.path.join(root, "ARF_code"); os.makedirs(code); shutil.copy(P, os.path.join(code, "algorithmic_reconstruction_framework.py"))


def search2(cols, bbox, a, b):     # same items, shifted to each window's start
    out = []
    for it in ITEMS:
        it = json.loads(json.dumps(it))
        t = datetime.fromisoformat(it["properties"]["datetime"][:19]) - start + a
        it["properties"]["datetime"] = t.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        out.append(it)
    return out


with contextlib.redirect_stdout(io.StringIO()):
    W, summ = a8.main(["--code", code, "--batch-dir", batch], search_fn=search2)
calls = []
def search3(*a):
    calls.append(1); return search2(*a)
with contextlib.redirect_stdout(io.StringIO()):
    W2, _ = a8.main(["--code", code, "--batch-dir", batch], search_fn=search3)
check(not calls and W2.round(9).equals(W.round(9)), "second run resumes from the cache: no catalogue query, same results")
check(len(W) == 2 and set(W.season) == {"dry", "wet"}, "one dry and one wet window; earthenv duplicate skipped")
r = W.iloc[0]
check(r.scenes == 2 and r.overpasses == 2 and r.scenes_unreadable == 1, "duplicate scene counted once, Landsat and late "
      "scene excluded, unreadable asset recorded")
unr = pd.read_csv(os.path.join(batch, "real_cloud_a8", "a8_unreadable.csv"))
check(len(unr) == 2 and unr.scene.str.contains("34KGD").all(), "unreadable asset listed per window in a8_unreadable.csv")
# independent recomputation for the window: overpass 1 = mixed (c1 kept), overpass 2 = cloudy
s_m, s_c = full, None
with rasterio.open(r_cloudy) as d:
    xs, ys = transform("EPSG:4326", d.crs, lon.ravel().tolist(), lat.ravel().tolist())
    s_c = np.array([v[0] for v in d.sample(zip(xs, ys), masked=False)], dtype=np.uint8).reshape(lon.shape)
oa, ca, ka = a8.classify(s_m); ob, cb, kb = a8.classify(s_c)
cloudy = 100 * ((ka & ~ca).sum() + (kb & ~cb).sum()) / (oa.sum() + ob.sum())
fs = 100 * np.mean([ca.sum() * 2 < oa.sum(), cb.sum() * 2 < ob.sum()])
reach = oa | ob
t95_sky = 8.5 if oa.sum() >= 0.95 * reach.sum() else (32.5 if reach.sum() else None)
clr_union = ca | cb
t95_clear = 8.5 if ca.sum() >= 0.95 * reach.sum() else (32.5 if clr_union.sum() >= 0.95 * reach.sum() else None)
ok = abs(r.cloudy_obs_pct - cloudy) < 1.0 and abs(r.fail_signal_overpass_pct - fs) < 1e-9
ok &= abs(r.s2_reachable_pct_of_region - 100 * reach.sum() / len(tiles)) < 1.0
ok &= (r.t95_reach_nocloud_h == t95_sky or (pd.isna(r.t95_reach_nocloud_h) and t95_sky is None))
ok &= (pd.isna(r.t95_reach_clear_h) and t95_clear is None) or abs(r.t95_reach_clear_h - t95_clear) < 1e-9
check(ok, f"cloudy {r.cloudy_obs_pct:.1f}% (indep. {cloudy:.1f}), Fail-Signal {r.fail_signal_overpass_pct:.0f}% "
          f"(indep. {fs:.0f}), reachable and t95 recomputed independently")
ps = summ["seasons"]["dry"]["pooled"]
n_obs = oa.sum() + ob.sum(); n_cld = (ka & ~ca).sum() + (kb & ~cb).sum()
full = [x >= 0.05 * len(tiles) for x in (oa.sum(), ob.sum())]
ok = ps["overpasses"] == 2 and ps["tile_observations"] == n_obs and abs(ps["cloudy_obs_pct"] - 100 * n_cld / n_obs) < 1e-9
ok &= ps["full_passes"] == sum(full) and ps["fail_signal_overpasses"] == int(ca.sum() * 2 < oa.sum()) + int(cb.sum() * 2 < ob.sum())
ok &= set(summ["seasons"]["dry"]["by_month"]) == {"08"} and set(summ["seasons"]["wet"]["by_month"]) == {"01"}
check(ok, f"pooled cloudy share {ps['cloudy_obs_pct']:.1f}% weighted by tile observations, full passes and months recomputed")
lv, c = a8.time_to((50.0, 100.0), [2.0, 1.0], [{1, 2}, {3}], 4)
check(lv == {50.0: 2.0, 100.0: None} and abs(c - 75.0) < 1e-12, "time_to orders by time and reports unreached levels")
okt = True
for tb in ("tab_a8_cloud.tex", "tab_a8_pooled.tex"):
    d = tempfile.mkdtemp(dir=root)
    open(os.path.join(d, "t.tex"), "w").write(
        "\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,textwidth=7.16in,"
        "columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
        + open(os.path.join(batch, "real_cloud_a8", tb)).read() + "x\n\\end{document}\n")
    rr = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=d, capture_output=True)
    okt &= rr.returncode == 0 and "Overfull \\hbox" not in open(os.path.join(d, "t.log"), errors="ignore").read()
check(okt, "tab_a8_cloud.tex and tab_a8_pooled.tex compile at IEEE width without overfull boxes")
print(W[["season", "window", "overpasses", "scenes", "s2_reachable_pct_of_region", "cloudy_obs_pct",
         "fail_signal_overpass_pct", "t95_reach_nocloud_h", "t95_reach_clear_h"]].round(2).to_string(index=False))
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
