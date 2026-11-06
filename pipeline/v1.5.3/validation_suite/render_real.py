"""Full real-data run of a pipeline file on a Colab run's saved inputs, to draw
its figures: Cells 0-5 with real Orekit on the run's cached CelesTrak elements
(as real_orekit_run.py: Orekit data from /home/claude/orekit-data, the run's
own tile list injected for Cell 1's WDPA discretisation), Cell 3 for the map
helpers, then Cells 6, 6B, 7, 7B, 7C, 8 and 9. Cell 7D (online archive) is not
run. Without the WDPA boundary the maps draw the outline of the ROI's cells.

Usage: render_real.py PIPELINE RUN_DIR OUT_DIR ROI FLEET CASE REF_EPOCH
"""
import sys, os, io, re, json, time, contextlib, typing, hashlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIPE, RUN_DIR, OUT, roi, fleet, case, ref_epoch = sys.argv[1:8]


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


src = open(PIPE, encoding="utf-8").read()
C = cells_of(src)
os.makedirs(OUT, exist_ok=True)
cache = os.path.join(OUT, "ARF_cache")
os.makedirs(cache, exist_ok=True)
for g in ("resource", "planet"):
    d = json.load(open(os.path.join(RUN_DIR, "ARF_cache", f"arf_gp_{g}.json")))
    d["fetched_at"] = time.time()             # re-stamped so Cell 2 uses the saved elements (no live fetch here)
    json.dump(d, open(os.path.join(cache, f"arf_gp_{g}.json"), "w"))

c1 = C["1"]
a = c1.index("from orekit_jpype.pyhelpers import setup_orekit_curdir")
b = c1.index("# =============================================================================\n# 3. SPATIAL")
c1 = (c1[:a] + "from java.io import File as _File\nfrom org.orekit.data import DataContext as _DC, DirectoryCrawler as _DCr\n"
      "_DC.getDefault().getDataProvidersManager().addProvider(_DCr(_File('/home/claude/orekit-data')))\n"
      "_orekit_data_loaded = True\n\n" + c1[b:])
call = re.search(r"target_tiles, roi_boundary_report = discretise_wdpa_roi\(.*?\)\n", c1, re.S)
c1 = c1[:call.start()] + "target_tiles, roi_boundary_report = set(_INJECTED_TILES), None\n" + c1[call.end():]

tiles = list(pd.read_csv(os.path.join(RUN_DIR, "arf_mosaic_tiles.csv"), usecols=["tile"])["tile"])
ns = {"__name__": "__main__", "display": lambda *a, **k: None, "ARF_SKIP_INSTALL": True, "ARF_ALLOW_MIXED_VERSIONS": True,
      "ARF_SHARED_CACHE_DIR": cache, "_INJECTED_TILES": tiles, "ARF_BOUNDARY_DIR": "/nonexistent",
      "ARF_CONFIG_OVERRIDES": {"active_roi": roi, "fleet": fleet, "planet_swath_case": case,
                               "reference_epoch_utc": ref_epoch}}
plt.show = lambda *a, **k: None
sink = io.StringIO()
times, errors = {}, {}
old = os.getcwd()
os.chdir(OUT)
try:
    with contextlib.redirect_stdout(sink):
        for name in ("0", "1", "2", "3", "4", "5", "6", "6B", "7", "7B", "7C", "8", "9"):
            code = c1 if name == "1" else C[name]
            t0 = time.time()
            try:
                exec(compile(code, f"cell{name}", "exec"), ns)
            except Exception as e:                   # Cell 3's study-area map needs the WDPA boundary
                if name != "3":
                    raise
                errors[name] = repr(e)[:300]
            times[name] = round(time.time() - t0, 1)
            plt.close("all")
finally:
    os.chdir(old)
open(os.path.join(OUT, "render_log.txt"), "w").write(sink.getvalue())
print("release", ns.get("__version__"), "| cell seconds", times, "| errors", errors)
from PIL import Image
for f in sorted(os.listdir(OUT)):
    if f.endswith(".png"):
        w, h = Image.open(os.path.join(OUT, f)).size
        print(f"{f}: {w / 300:.2f} x {h / 300:.2f} in")
