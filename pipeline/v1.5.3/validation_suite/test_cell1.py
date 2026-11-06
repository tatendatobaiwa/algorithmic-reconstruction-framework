"""Cell 1 (ROI loading and discretisation), on PRIVATE copies of the WDPA
boundaries (WDPA licence: never copied to outputs; nothing of their geometry
is printed or saved here -- only tile counts, hashes and pass/fail).
Usage: test_cell1.py PIPELINE
A. Lower Zambezi (7962) and Karoo (555563860): the tiles equal the tile lists
   of the 12 October 2026 Colab runs (SHA-256 of the sorted IDs, as the parity
   file records them); area check within tolerance; report fields.
B. CKGR (7510): its tiles are a subset of the CKGR + Okavango ROI tiles of the
   Colab runs (the Okavango WHS file, 555577555, is not held here).
C. Determinism: loading twice, and from a zip (nested in a folder), gives the
   same tiles; split parts are de-duplicated by content.
D. Identity checks: the Okavango Ramsar file (145516) offered as a site that
   must be a World Heritage Site is refused; a wrong REP_AREA is refused; a
   missing site and a point-only record give their instructions.
E. Sample tile (arf_cell1_validation.json): v1.5.2 uses min(tiles), the same
   in every Python process; next(iter(set)) of v1.5.1 depends on the
   per-process string-hash salt (shown by running both under two
   PYTHONHASHSEED values)."""
import sys, os, ast, json, shutil, hashlib, zipfile, tempfile, subprocess, logging
import numpy as np, pandas as pd, h3
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
PRIV = "/home/claude/wdpa_private"
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


def sha(tiles):
    return hashlib.sha256("\n".join(sorted(tiles)).encode()).hexdigest()


C = e2e15.cells_of(open(P, encoding="utf-8").read())
tree = ast.parse(C["1"])
keep = []
for n in tree.body:
    if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
        keep.append(n)
    elif isinstance(n, ast.Assign) and any(getattr(t, "id", "") in ("WDPA_AREA_TOLERANCE", "WDPA_TILE_TOLERANCE",
                                                                     "WDPA_REP_AREA_TOLERANCE", "_MONTHS",
                                                                     "ROI_DEFINITIONS") for t in n.targets):
        keep.append(n)
tree.body = keep


class ConfigurationError(Exception):
    pass


ns = {"os": os, "h3": h3, "np": np, "logger": logging.getLogger("t1"), "ConfigurationError": ConfigurationError,
      "List": list, "Tuple": tuple, "Dict": dict, "Set": set, "Optional": object, "Any": object}
exec(compile(tree, "cell1_defs", "exec"), ns)
RD = ns["ROI_DEFINITIONS"]
exp = json.load(open("/mnt/user-data/outputs/arf_parity_expected_v150.json"))
exp_sha = {}
for name, r in exp["runs"].items():
    exp_sha.setdefault(r["roi"], r["tile_set_sha256"])

# A
for roi, sid in (("lower_zambezi", 7962), ("karoo_nature_reserve", 555563860)):
    tiles, reps = ns["discretise_wdpa_roi"](RD[roi], 7, boundary_dir=PRIV)
    rep = reps[0]
    check(sha(tiles) == exp_sha.get(roi) and rep["site_id"] == sid and abs(rep["area_deviation_pct"]) <= rep["area_tolerance_pct"],
          f"{roi}: {len(tiles)} tiles, identical to the Colab runs' tile list (SHA-256); tiled area "
          f"{rep['area_deviation_pct']:+.2f}% of GIS_AREA (tolerance {rep['area_tolerance_pct']:.1f}%)")
# B
ckgr_roi = pd.read_csv("/home/claude/realdata/ckgr_openPO_mosaic_tiles.csv", usecols=["tile"])["tile"]
spec = RD["ckgr_okavango"]["wdpa_sites"][7510]
t7510, rep = ns["load_wdpa_site"](7510, spec, PRIV, 7)
check(t7510 <= set(ckgr_roi) and abs(rep["area_deviation_pct"]) <= rep["area_tolerance_pct"] and "Game Reserve" in rep["designation"][0],
      f"CKGR (7510): {len(t7510)} tiles, all inside the {len(ckgr_roi)}-tile CKGR + Okavango ROI of the Colab runs; "
      f"designation {rep['designation']}")
# C
t2, _ = ns["load_wdpa_site"](7510, spec, PRIV, 7)
tmp = tempfile.mkdtemp(prefix="wdpa_zip_", dir=PRIV)          # private area only
try:
    os.makedirs(os.path.join(tmp, "nested", "deeper"))
    zp = os.path.join(tmp, "nested", "deeper", "WDPA_WDOECM_Oct2026_Public_7510_shp.zip")
    with zipfile.ZipFile(zp, "w") as z:
        for f in os.listdir(PRIV):
            if "_7510_shp-polygons." in f:
                z.write(os.path.join(PRIV, f), arcname=f"WDPA_WDOECM_Oct2026_Public_7510_shp_0/{f}")
    t_zip, rep_zip = ns["load_wdpa_site"](7510, spec, os.path.join(tmp, "nested"), 7)
    # extracted copy AND zip side by side: same content, read once (no doubled area)
    for f in os.listdir(PRIV):
        if "_7510_shp-polygons." in f:
            shutil.copy(os.path.join(PRIV, f), os.path.join(tmp, "nested", f))
    t_both, rep_both = ns["load_wdpa_site"](7510, spec, os.path.join(tmp, "nested"), 7)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
check(t2 == t7510 and t_zip == t7510 and t_both == t7510 and rep_both["n_parts"] == 1
      and abs(rep_both["wdpa_gis_area_km2"] - rep["wdpa_gis_area_km2"]) < 1e-6,
      "deterministic; read from a zip two folders deep; zip + extracted copy read once (area not doubled)")
# D
errs = []
try:
    ns["load_wdpa_site"](145516, {"name": "Okavango Delta (World Heritage Site)", "designation_contains": "World Heritage"},
                         PRIV, 7)
    errs.append(("ramsar", False))
except ConfigurationError as e:
    errs.append(("ramsar", "Ramsar" in str(e) or "designated" in str(e)))
try:
    ns["load_wdpa_site"](7962, {"name": "Lower Zambezi", "designation_contains": None, "expected_rep_area_km2": 5000.0}, PRIV, 7)
    errs.append(("rep_area", False))
except ConfigurationError as e:
    errs.append(("rep_area", "REP_AREA" in str(e)))
try:
    ns["load_wdpa_site"](999999999, {"name": "Nowhere"}, PRIV, 7)
    errs.append(("missing", False))
except ConfigurationError as e:
    errs.append(("missing", "protectedplanet.net/999999999" in str(e)))
tmp = tempfile.mkdtemp(prefix="wdpa_pt_", dir=PRIV)
try:
    open(os.path.join(tmp, "WDPA_WDOECM_Oct2026_Public_555547578_shp-points.shp"), "wb").close()
    try:
        ns["load_wdpa_site"](555547578, {"name": "Middle Zambezi"}, tmp, 7)
        errs.append(("point", False))
    except ConfigurationError as e:
        errs.append(("point", "NO boundary polygon" in str(e)))
finally:
    shutil.rmtree(tmp, ignore_errors=True)
check(all(ok for _, ok in errs), f"identity checks refuse: Ramsar file for a World Heritage Site, wrong REP_AREA, "
                                 f"missing site, point-only record ({errs})")
# E
snip = ("import pandas as pd; t = set(pd.read_csv('/home/claude/realdata/lz_tiles.csv')['tile']); "
        "print(next(iter(t)), min(t))")
outs = [subprocess.run([sys.executable, "-c", snip], capture_output=True, text=True,
                       env={**os.environ, "PYTHONHASHSEED": str(seed)}).stdout.split() for seed in (1, 2, 3)]
v151 = {o[0] for o in outs}
v152 = {o[1] for o in outs}
src1 = C["1"]
check(len(v152) == 1 and "_sample_tile = min(target_tiles)" in src1,
      f"sample tile: min(tiles) identical under three hash seeds ({v152}); next(iter(set)) gave {len(v151)} different "
      f"tiles (v1.5.1 rule)")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell1.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
