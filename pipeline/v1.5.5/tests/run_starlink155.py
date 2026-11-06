"""Full run of the 1.5.4 Starlink proxy at a past window start (dense-constellation
arm, amendment A3) with real Orekit, on SYNTHETIC element sets
(make_synthetic_starlink.py): Cells 0-9, Cell 7D included (for the proxy it records 'skipped' only).
Substitutions as in render_real.py: Orekit data from /home/claude/orekit-data, and
the ROI tiles injected from a run's tile list (the WDPA boundary files are not
used here). Records each cell's wall-clock time.

Usage: run_starlink154.py PIPELINE OUT_DIR N WINDOW_START CACHE_DIR TILES_CSV [SUN_DEG|none] [STOP_AFTER_CELL]
"""
import sys, os, io, re, json, time, shutil, contextlib, traceback
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIPE, OUT, N, EPOCH, CACHE, TILES = sys.argv[1:7]
SUN = None if len(sys.argv) < 8 or sys.argv[7] == "none" else float(sys.argv[7])
STOP = sys.argv[8] if len(sys.argv) > 8 else None


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


C = cells_of(open(PIPE, encoding="utf-8").read())
os.makedirs(OUT, exist_ok=True)
cache = os.path.join(OUT, "ARF_cache")
os.makedirs(cache, exist_ok=True)
for f in os.listdir(CACHE):
    if f.startswith("arf_gp_history_starlink_"):
        shutil.copy2(os.path.join(CACHE, f), cache)
c1 = C["1"]
a = c1.index("from orekit_jpype.pyhelpers import setup_orekit_curdir")
b = c1.index("# =============================================================================\n# 3. SPATIAL")
c1 = (c1[:a] + "from java.io import File as _File\nfrom org.orekit.data import DataContext as _DC, DirectoryCrawler as _DCr\n"
      "_DC.getDefault().getDataProvidersManager().addProvider(_DCr(_File('/home/claude/orekit-data')))\n"
      "_orekit_data_loaded = True\n\n" + c1[b:])
call = re.search(r"target_tiles, roi_boundary_report = discretise_wdpa_roi\(.*?\)\n", c1, re.S)
c1 = c1[:call.start()] + "target_tiles, roi_boundary_report = set(_INJECTED_TILES), None\n" + c1[call.end():]
tiles = list(pd.read_csv(TILES, usecols=["tile"])["tile"])
ov = {"active_roi": "ckgr_okavango", "fleet": "starlink", "active_n": int(N), "analysis_epoch": EPOCH}
if SUN is not None:
    ov["proxy_min_sun_elev_deg"] = SUN
if os.environ.get("ARF_PROXY_SWATH_KM"):            # v1.5.5 test: proxy ground swath
    ov["proxy_swath_km"] = float(os.environ["ARF_PROXY_SWATH_KM"])
if os.environ.get("ARF_NO_CLOUD_SWEEP") == "1":       # dense-arm option: skip the cloud-sensitivity sweep
    ov["cloud_sweep_lengths_km"], ov["cloud_sweep_times_h"] = [], []
ns = {"__name__": "__main__", "display": lambda *a, **k: None, "ARF_SKIP_INSTALL": True, "ARF_ALLOW_MIXED_VERSIONS": True,
      "ARF_SHARED_CACHE_DIR": cache, "_INJECTED_TILES": tiles, "ARF_BOUNDARY_DIR": "/nonexistent",
      "ARF_CONFIG_OVERRIDES": ov}
plt.show = lambda *a, **k: None
sink = io.StringIO()
times, errors = {}, {}
old = os.getcwd()
os.chdir(OUT)
try:
    with contextlib.redirect_stdout(sink):
        # Cell 7D included: for the proxy it only records "skipped" (no network); as in Colab
        for name in ("0", "1", "2", "3", "4", "5", "6", "6B", "7", "7B", "7C", "7D", "8", "9"):
            code = c1 if name == "1" else C[name]
            t0 = time.time()
            try:
                exec(compile(code, f"cell{name}", "exec"), ns)
            except Exception:
                errors[name] = traceback.format_exc()[-3000:]
                times[name] = round(time.time() - t0, 1)
                break
            times[name] = round(time.time() - t0, 1)
            plt.close("all")
            if name == STOP:
                break
finally:
    os.chdir(old)
open(os.path.join(OUT, "run_log.txt"), "w").write(sink.getvalue())
if "refined_events_df" in ns:                      # v1.5.5 test: keep Cell 5 output for comparison
    ns["refined_events_df"].to_csv(os.path.join(OUT, "refined_events_cell5.csv"), index=False)
if "sensor_config" in ns:
    json.dump({"footprint_radius_km": ns["sensor_config"].footprint_radius_km, "min_altitude_km": ns["sensor_config"].min_altitude_km}, open(os.path.join(OUT, "sensor.json"), "w"))
if "constellation" in ns:
    json.dump([{"id": c["id"], "swath_km": c.get("swath_km"), "min_sun_elev_deg": c.get("min_sun_elev_deg"), "alt_km": float((c["tle"].getMeanMotion() and ((398600.4418 / (c["tle"].getMeanMotion()) ** 2) ** (1 / 3) - 6378.137)))} for c in ns["constellation"]], open(os.path.join(OUT, "constellation.json"), "w"))
res = {"N": int(N), "epoch": EPOCH, "sun": SUN, "cell_seconds": times, "total_seconds": round(sum(times.values()), 1),
       "errors": errors, "n_constellation": len(ns.get("constellation", [])),
       "out_dir": ns.get("OUT_DIR"), "fleet_tag": ns.get("FLEET_TAG"),
       "op_stats": ns.get("_op_stats"), "reference_epoch": ns.get("_reference_epoch"),
       "n_candidate_events": int(len(ns["refined_events_df"])) if "refined_events_df" in ns else None,
       "peak_rss_gb": round(__import__("resource").getrusage(__import__("resource").RUSAGE_SELF).ru_maxrss / 1024 ** 2, 2)}
json.dump(res, open(os.path.join(OUT, "timing.json"), "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
