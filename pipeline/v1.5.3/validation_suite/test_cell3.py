"""Cell 3 (study-area map): no vector copy of the WDPA boundary in the PDF.
Usage: test_cell3.py PIPELINE
Cells 0-3 run with real Orekit for the Lower Zambezi ROI, loading its boundary
from the PRIVATE WDPA copy (WDPA licence: the boundary must not be
redistributed, so a saved PDF may carry it only as pixels). The map is drawn
three times, changing only how draw_outline draws (test-only switches
inserted into the cell's source):
  A  as released (WDPA outline rasterised)
  B  outline forced to vector
  C  outline not drawn
PyMuPDF lists every vector path in each PDF (page.get_drawings()):
  A must equal C exactly (the outline adds no vector path) and hold more
  embedded images than C (the outline as pixels); B must hold more vector
  paths than C (so the comparison can see an outline when there is one).
The PDFs are written to a private temporary folder and deleted; only counts
are reported. Also checks the caption names the stations as notional."""
import sys, os, io, re, json, time, shutil, tempfile, contextlib
import numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
RD = "/home/claude/realdata/"
PRIV = "/home/claude/wdpa_private"
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


C = e2e15.cells_of(open(P, encoding="utf-8").read())
c1 = C["1"]
a = c1.index("from orekit_jpype.pyhelpers import setup_orekit_curdir")
b = c1.index("# =============================================================================\n# 3. SPATIAL")
c1 = (c1[:a] + "from java.io import File as _File\nfrom org.orekit.data import DataContext as _DC, DirectoryCrawler as _DCr\n"
      "_DC.getDefault().getDataProvidersManager().addProvider(_DCr(_File('/home/claude/orekit-data')))\n"
      "_orekit_data_loaded = True\n\n" + c1[b:])
c3 = C["3"]
anchor = '    lines = geo_lines(geoms)\n    if lines:\n        lc = LineCollection(lines, colors=color or MAP_COLORS["boundary"]'
assert anchor in c3, "draw_outline body changed: update the test switches"
c3 = c3.replace(anchor, '    if globals().get("_TEST_SKIP_OUTLINE"):\n        return\n'
                        '    if globals().get("_TEST_VECTOR_OUTLINE"):\n        rasterized = False\n' + anchor)

work = tempfile.mkdtemp(prefix="c3_", dir=PRIV)
try:
    cache = os.path.join(work, "ARF_cache")
    os.makedirs(cache)
    for g in ("resource", "planet"):
        d = json.load(open(RD + f"ckgr_openPO_arf_gp_{g}.json"))
        d["fetched_at"] = time.time()
        json.dump(d, open(os.path.join(cache, f"arf_gp_{g}.json"), "w"))
    ns = {"__name__": "__main__", "display": lambda *a, **k: None, "ARF_SKIP_INSTALL": True,
          "ARF_SHARED_CACHE_DIR": cache, "ARF_BOUNDARY_DIR": PRIV,
          "ARF_CONFIG_OVERRIDES": {"active_roi": "lower_zambezi", "fleet": "earth-resources+planet",
                                   "planet_swath_case": "conservative",
                                   "reference_epoch_utc": "2026-10-12T08:27:09.082080"}}
    import matplotlib
    matplotlib.use("Agg")
    results = {}
    old = os.getcwd()
    os.chdir(work)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            for name, code in (("0", C["0"]), ("1", c1), ("2", C["2"])):
                exec(compile(code, f"cell{name}", "exec"), ns)
            for variant in ("A", "B", "C"):
                ns["_TEST_SKIP_OUTLINE"] = variant == "C"
                ns["_TEST_VECTOR_OUTLINE"] = variant == "B"
                exec(compile(c3, "cell3", "exec"), ns)
                shutil.copy(os.path.join(ns["OUT_DIR"], "arf_coverage_map.pdf"), os.path.join(work, f"map_{variant}.pdf"))
                results[variant] = json.load(open(os.path.join(ns["OUT_DIR"], "arf_cell3_metadata.json")))
    finally:
        os.chdir(old)
    import pymupdf
    info = {}
    for v in ("A", "B", "C"):
        doc = pymupdf.open(os.path.join(work, f"map_{v}.pdf"))
        page = doc[0]
        dr = page.get_drawings()
        info[v] = {"paths": len(dr), "items": sum(len(d["items"]) for d in dr),
                   "images": len(page.get_images(full=True)),
                   "sig": [(d.get("type"), round(d["rect"].x0, 2), round(d["rect"].y0, 2), round(d["rect"].x1, 2),
                            round(d["rect"].y1, 2), len(d["items"])) for d in dr]}
        doc.close()
    check(info["A"]["sig"] == info["C"]["sig"],
          f"released map: identical vector content with and without the WDPA outline ({info['A']['paths']} paths, "
          f"{info['A']['items']} items each): no vector copy of the boundary")
    check(info["A"]["images"] > info["C"]["images"],
          f"the outline is in the PDF as pixels: {info['A']['images']} embedded images vs {info['C']['images']} without it")
    check(info["B"]["items"] > info["C"]["items"],
          f"control: a vector outline adds {info['B']['items'] - info['C']['items']} path items, so the comparison sees one")
    cap = results["A"]["draft_caption"]
    check("notional ground stations" in cap and "not the operators' ground networks" in cap
          and results["A"]["outline_source"] == "WDPA",
          "caption: stations described as notional (approximate coordinates, not the operators' networks); outline "
          "source WDPA; WDPA citation present" if results["A"].get("wdpa_citation") else "caption check")
finally:
    shutil.rmtree(work, ignore_errors=True)

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell3.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
