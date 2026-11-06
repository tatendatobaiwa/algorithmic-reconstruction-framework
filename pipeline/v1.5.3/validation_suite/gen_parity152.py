"""Reference values for arf_colab_parity.py (v1.5.2): Cells 6-8 of the release
on the SAVED inputs of the 12 Oct 2026 Colab runs (arf_refined_events.csv as
written by Colab, unmodified; tile lists), by the parity script's own
procedure (replay + collect). Also compares with the v1.5.0 reference values.
Usage: gen_parity152.py PIPELINE PARITY_SCRIPT EXPECTED_V150 OUT_JSON [run folder ...]"""
import sys, os, json, time, hashlib, importlib.util, datetime
import pandas as pd
P, PAR, EXP150, OUT = sys.argv[1:5]
only = sys.argv[5:]
spec = importlib.util.spec_from_file_location("arf_colab_parity", PAR)
par = importlib.util.module_from_spec(spec); spec.loader.exec_module(par)
exp = json.load(open(EXP150))
RD = "/home/claude/realdata/"
TILES = {"ckgr_okavango": RD + "ckgr_openPO_mosaic_tiles.csv", "lower_zambezi": RD + "lz_tiles.csv",
         "karoo_nature_reserve": RD + "karoo_tiles.csv"}
out = json.load(open(OUT)) if os.path.exists(OUT) else {
    "release": None, "pipeline_sha256": par.sha256_file(P),
    "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "procedure": "arf_colab_parity.replay(...) then collect(...)", "runs": {}}
src = open(P, encoding="utf-8").read()
import re
out["release"] = re.search(r'^__version__ = "([^"]+)"', src, re.M).group(1)
out["pipeline_sha256"] = par.sha256_file(P)
for folder, r in exp["runs"].items():
    if only and folder not in only:
        continue
    csv = RD + f"{r['name']}_refined_events.csv"
    assert par.sha256_file(csv) == r["refined_events_sha256"], f"{csv}: not the saved Colab input"
    tiles = list(pd.read_csv(TILES[r["roi"]], usecols=["tile"])["tile"])
    assert par.tile_set_sha(tiles) == r["tile_set_sha256"]
    t0 = time.time()
    ns = par.replay(P, csv, tiles, r["orbital_period_min"], r["roi"], r["fleet"], r["planet_swath_case"], r["n_satellites"])
    vals = par.collect(ns)
    rec = {k: r[k] for k in ("name", "roi", "fleet", "planet_swath_case", "n_satellites", "orbital_period_min",
                             "refined_events_sha256", "refined_events_bytes", "tile_set_sha256")}
    rec["values"] = vals
    diffs = par.compare(r["values"], vals)
    rec["_vs_v150"] = diffs
    rec["_seconds"] = round(time.time() - t0, 1)
    out["runs"][folder] = rec
    json.dump(out, open(OUT, "w"), indent=1)
    print(folder, r["name"], f"{rec['_seconds']} s", "same as v1.5.0" if not diffs else f"{len(diffs)} differ: "
          + "; ".join(f"{d['field']} {d['expected']} -> {d['got']}" for d in diffs[:6]), flush=True)
