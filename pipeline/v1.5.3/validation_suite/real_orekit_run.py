"""Real-Orekit run of Cells 0-5 of a pipeline file on the CACHED orbital
elements of a Colab run (15 October 2026 runs: ARF_cache of the run folder),
then a comparison of the refined events with the ones Colab wrote.

What is real here: Cell 0, 2, 4, 5 run unmodified (orekit-jpype 13.1.9.0 and
jdk4py 25.0.2.1, the versions in the Colab manifest of that run); Cell 1 runs
with two substitutions only: (1) Orekit data from /home/claude/orekit-data
(hand-written tai-utc.dat + IERS finals2000A.all from astropy-iers-data,
because gitlab.orekit.org is not reachable from here), and (2) the ROI tiles
injected from the run's own tile list (the Okavango World Heritage Site
boundary, WDPA 555577555, is not available here; the tile list's SHA-256 is
checked against the reference). Cell 2 reads the run's CelesTrak GP JSON
cache (re-stamped as fresh so it is used instead of a live fetch).

Usage: real_orekit_run.py --pipe PIPELINE --run NAME --out DIR --gp-dir RUN/ARF_cache --tiles RUN/arf_mosaic_tiles.csv --colab-refined RUN/arf_refined_events.csv --roi ROI --fleet FLEET [--case CASE] --ref-epoch ISO
"""
import sys, os, io, re, json, time, shutil, hashlib, contextlib, typing
import numpy as np, pandas as pd

import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--pipe", required=True)
_ap.add_argument("--run", required=True)
_ap.add_argument("--out", required=True)
_ap.add_argument("--gp-dir", required=True, help="the run's ARF_cache folder (arf_gp_resource.json, arf_gp_planet.json)")
_ap.add_argument("--tiles", required=True, help="the run's arf_mosaic_tiles.csv")
_ap.add_argument("--colab-refined", required=True, help="the run's arf_refined_events.csv from Colab")
_ap.add_argument("--roi", required=True)
_ap.add_argument("--fleet", required=True)
_ap.add_argument("--case", default="conservative")
_ap.add_argument("--ref-epoch", required=True)
_a = _ap.parse_args()
PIPE, RUN, OUT = _a.pipe, _a.run, _a.out
LEGACY, WITH3, EXTRA_OV = False, False, {}
roi, fleet, case, ref_epoch, tiles_csv = _a.roi, _a.fleet, _a.case, _a.ref_epoch, _a.tiles
gp_files = {g: os.path.join(_a.gp_dir, f"arf_gp_{g}.json") for g in ("resource", "planet")}
COLAB_REFINED = _a.colab_refined


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


src = open(PIPE, encoding="utf-8").read()
C = cells_of(src)
os.makedirs(OUT, exist_ok=True)
cache = os.path.join(OUT, "ARF_cache")
os.makedirs(cache, exist_ok=True)
for g, p in gp_files.items():
    d = json.load(open(p))
    d["fetched_at"] = time.time()            # re-stamp so Cell 2 uses it (no live fetch possible here)
    json.dump(d, open(os.path.join(cache, f"arf_gp_{g}.json"), "w"))

# --- Cell 1 with the two substitutions ------------------------------------------------
c1 = C["1"]
orekit_block_start = c1.index("from orekit_jpype.pyhelpers import setup_orekit_curdir")
orekit_block_end = c1.index("# =============================================================================\n# 3. SPATIAL")
c1 = (c1[:orekit_block_start]
      + "from java.io import File as _File\nfrom org.orekit.data import DataContext as _DC, DirectoryCrawler as _DCr\n"
        "_DC.getDefault().getDataProvidersManager().addProvider(_DCr(_File('/home/claude/orekit-data')))\n"
        "_orekit_data_loaded = True\n\n"
      + c1[orekit_block_end:])
call = re.search(r"target_tiles, roi_boundary_report = discretise_wdpa_roi\(.*?\)\n", c1, re.S)
assert call, "Cell 1 discretisation call not found"
c1 = c1[:call.start()] + "target_tiles, roi_boundary_report = set(_INJECTED_TILES), None\n" + c1[call.end():]

tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
ns = {"__name__": "__main__", "display": lambda *a, **k: None, "ARF_SKIP_INSTALL": True,
      "ARF_SHARED_CACHE_DIR": cache, "_INJECTED_TILES": tiles, "ARF_BOUNDARY_DIR": "/nonexistent",
      "ARF_CONFIG_OVERRIDES": {"active_roi": roi, "fleet": fleet, "planet_swath_case": case,
                               "reference_epoch_utc": ref_epoch, "legacy_strip_end_caps": LEGACY, **EXTRA_OV}}
sink = io.StringIO()
times = {}
old = os.getcwd()
os.chdir(OUT)
try:
    with contextlib.redirect_stdout(sink):
        for name, code in (("0", C["0"]), ("1", c1), ("2", C["2"])) + ((("3", C["3"]),) if WITH3 else ()) + \
                          (("4", C["4"]), ("5", C["5"])):
            t0 = time.time()
            exec(compile(code, f"cell{name}", "exec"), ns)
            times[name] = round(time.time() - t0, 1)
finally:
    os.chdir(old)
log = sink.getvalue()
open(os.path.join(OUT, "run_log.txt"), "w").write(log)
print("cell seconds", times)
print("sim_epoch", ns["sim_epoch"], "| satellites", len(ns["constellation"]), "| median alt", ns["CONSTELLATION_ALTITUDE_KM"])

# --- comparison with the Colab refined events ------------------------------------------
mine = ns["refined_events_df"].copy()
mine["sat"], mine["tile"] = mine["sat"].astype(str), mine["tile"].astype(str)
ref = pd.read_csv(COLAB_REFINED)
key = ["sat", "tile", "t_min"]
m = ref.merge(mine, on=key, how="outer", suffixes=("_colab", "_here"), indicator=True)
both = m[m["_merge"] == "both"]
res = {"run": RUN, "legacy_end_caps": LEGACY, "extra_overrides": EXTRA_OV, "rows_colab": len(ref), "rows_here": len(mine),
       "rows_both": len(both), "rows_only_colab": int((m["_merge"] == "left_only").sum()),
       "rows_only_here": int((m["_merge"] == "right_only").sum())}
for col in ("A_nadir_deg", "theta_max_deg", "sun_elev_deg", "D_downlink", "T_wait_min"):
    dd = (both[f"{col}_colab"] - both[f"{col}_here"]).abs()
    res[f"max_abs_diff_{col}"] = float(dd.max()) if len(dd) else None
    res[f"p99_abs_diff_{col}"] = float(dd.quantile(0.99)) if len(dd) else None
only = m[m["_merge"] != "both"]
if len(only):
    a = only.get("A_nadir_deg_colab").fillna(only.get("A_nadir_deg_here"))
    th = only.get("theta_max_deg_colab").fillna(only.get("theta_max_deg_here"))
    res["unmatched_rows_margin_to_edge_deg_max"] = float((th - a).abs().max())
    res["unmatched_rows_margin_to_edge_deg_median"] = float((th - a).abs().median())
res["cell5_summary"] = json.load(open(os.path.join(OUT, "arf_cell5_summary.json")))
res["cell4_summary"] = json.load(open(os.path.join(OUT, "arf_cell4_summary.json")))
res["cell2_summary"] = {k: v for k, v in json.load(open(os.path.join(OUT, "arf_cell2_summary.json"))).items()
                        if k in ("n_satellites", "median_alt_km", "t_max_wait_min", "sim_epoch", "reference_epoch",
                                 "reference_epoch_source")}
res["tile_set_sha256"] = hashlib.sha256("\n".join(sorted(tiles)).encode()).hexdigest()
json.dump(res, open(os.path.join(OUT, "comparison.json"), "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
