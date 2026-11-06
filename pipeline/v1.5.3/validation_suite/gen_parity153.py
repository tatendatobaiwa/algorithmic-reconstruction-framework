"""Reference values for arf_colab_parity.py (v1.5.3), by the parity script's own
procedure (replay + collect) on saved inputs, with release 1.5.3:

  oct15  the five complete Colab runs of 15 October 2026 (inputs downloaded from
         their Drive folders, checked against each manifest's SHA-256); keyed by
         the Drive folder name; each value is also compared with the value the
         Colab run itself recorded (scenario_results and statuses in its manifest)
  oct12  the 15 runs of 12 October 2026 (inputs as in arf_parity_expected_v152.json);
         each value is compared with the 1.5.2 reference

Usage: gen_parity153.py PIPELINE PARITY_SCRIPT oct15 OUT_JSON RUN_DIR...
       gen_parity153.py PIPELINE PARITY_SCRIPT oct12 OUT_JSON EXPECTED_V152
"""
import sys, os, re, json, time, datetime, importlib.util
import pandas as pd

P, PAR, mode, OUT, *rest = sys.argv[1:]
spec = importlib.util.spec_from_file_location("arf_colab_parity", PAR)
par = importlib.util.module_from_spec(spec)
spec.loader.exec_module(par)
src = open(P, encoding="utf-8").read()
out = {"release": re.search(r'^__version__ = "([^"]+)"', src, re.M).group(1), "pipeline_sha256": par.sha256_file(P),
       "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
       "procedure": "arf_colab_parity.replay(...) then collect(...)", "inputs": None, "runs": {}}


def colab_values(man, key):
    """The compared values as the Colab run recorded them (manifest)."""
    d = man["scenario_results"][key]
    vals = {p: par._get(d, p) for p in par.METRICS}
    hv = man.get("hypothesis_validation") or {}
    lab = None
    for h in par.STATUSES:
        a = (hv.get(h) or {}).get("assessments") or {}
        lab = lab or next((l for l in a if l.endswith("fleet=" + key.split("__N", 1)[1])), None)
        vals[f"status/{h}"] = (a.get(lab) or {}).get("status") if lab else None
    return json.loads(json.dumps(vals, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


if mode == "oct15":
    out["inputs"] = ("the five complete runs of 15 October 2026 in their Drive folders (window start "
                     "2026-10-15T08:36:44.671104 UTC); nothing to restore")
    for rd in rest:
        man = json.load(open(os.path.join(rd, "arf_reproducibility_manifest.json")))
        inp = man["inputs"]
        fleet, case, roi = inp["fleet"], inp["config"]["planet_swath_case"], inp["active_roi"]
        key = f"{roi}__N{par.fleet_tag(fleet, case)}"
        folder = os.path.basename(man["environment"]["out_dir"])
        refined = os.path.join(rd, "arf_refined_events.csv")
        tiles_csv = os.path.join(rd, "arf_mosaic_tiles.csv")
        ov = man["output_validation"]["files"]
        assert par.sha256_file(refined) == ov["arf_refined_events.csv"]["sha256"], refined
        assert par.sha256_file(tiles_csv) == ov["arf_mosaic_tiles.csv"]["sha256"], tiles_csv
        tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
        orbit = man["scenario_results"][key]["orbital_period_min"]
        t0 = time.time()
        ns = par.replay(P, refined, tiles, orbit, roi, fleet, case, inp["n_satellites"])
        vals = par.collect(ns)
        rec = {"name": os.path.basename(rd.rstrip("/")), "roi": roi, "fleet": fleet, "planet_swath_case": case,
               "n_satellites": inp["n_satellites"], "orbital_period_min": orbit,
               "refined_events_sha256": par.sha256_file(refined), "refined_events_bytes": os.path.getsize(refined),
               "tile_set_sha256": par.tile_set_sha(tiles), "values": vals}
        rec["_vs_colab_run"] = par.compare(colab_values(man, key), vals)
        rec["_seconds"] = round(time.time() - t0, 1)
        out["runs"][folder] = rec
        json.dump(out, open(OUT, "w"), indent=1)
        print(folder, f"{rec['_seconds']} s", "same as the Colab run" if not rec["_vs_colab_run"] else
              f"{len(rec['_vs_colab_run'])} differ: " + "; ".join(f"{d['field']} {d['expected']} -> {d['got']}"
                                                                for d in rec["_vs_colab_run"][:6]), flush=True)
else:
    exp = json.load(open(rest[0]))
    out["inputs"] = ("the 15 runs of 12 October 2026: unzip ARF_12Oct2026_saved_inputs.zip into a new folder and "
                     "pass --drive pointing at it (never into MyDrive itself: two folder names are now used by "
                     "the runs of 15 October)")
    RD = "/home/claude/realdata/"
    TILES = {"ckgr_okavango": RD + "ckgr_openPO_mosaic_tiles.csv", "lower_zambezi": RD + "lz_tiles.csv",
             "karoo_nature_reserve": RD + "karoo_tiles.csv"}
    for folder, r in exp["runs"].items():
        csv = RD + f"{r['name']}_refined_events.csv"
        assert par.sha256_file(csv) == r["refined_events_sha256"], f"{csv}: not the saved Colab input"
        tiles = list(pd.read_csv(TILES[r["roi"]], usecols=["tile"])["tile"])
        assert par.tile_set_sha(tiles) == r["tile_set_sha256"]
        t0 = time.time()
        ns = par.replay(P, csv, tiles, r["orbital_period_min"], r["roi"], r["fleet"], r["planet_swath_case"],
                        r["n_satellites"])
        vals = par.collect(ns)
        rec = {k: r[k] for k in ("name", "roi", "fleet", "planet_swath_case", "n_satellites", "orbital_period_min",
                                 "refined_events_sha256", "refined_events_bytes", "tile_set_sha256")}
        rec["values"] = vals
        rec["_vs_v152"] = par.compare(r["values"], vals, rtol=0.0)
        rec["_seconds"] = round(time.time() - t0, 1)
        out["runs"][folder] = rec
        json.dump(out, open(OUT, "w"), indent=1)
        print(folder, r["name"], f"{rec['_seconds']} s", "same as v1.5.2" if not rec["_vs_v152"] else
              f"{len(rec['_vs_v152'])} differ", flush=True)
print(f"{OUT}: {len(out['runs'])} runs")
