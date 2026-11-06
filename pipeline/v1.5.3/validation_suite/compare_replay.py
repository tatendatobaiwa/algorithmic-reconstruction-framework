"""Parity of the 15 October 2026 Colab runs (release 1.5.2) in the sandbox:
Cells 6-8 of the shipped release, replayed on each run's saved inputs
(arf_refined_events.csv and arf_mosaic_tiles.csv from its Drive folder,
checked against the manifest's SHA-256), compared field by field with the
scenario_results and hypothesis statuses the run recorded in Colab.

Usage: compare_replay.py PIPELINE PARITY_SCRIPT RUN_DIR OUT_JSON
The run folder holds the downloaded files; they are only read."""
import sys, os, json, math, time, hashlib, importlib.util

P, PAR, RUN_DIR, OUT = sys.argv[1:5]
spec = importlib.util.spec_from_file_location("arf_colab_parity", PAR)
par = importlib.util.module_from_spec(spec)
spec.loader.exec_module(par)
import pandas as pd

man = json.load(open(os.path.join(RUN_DIR, "arf_reproducibility_manifest.json")))
inp = man["inputs"]
fleet, case, roi = inp["fleet"], inp["config"]["planet_swath_case"], inp["active_roi"]
tag = par.fleet_tag(fleet, case)
key = f"{roi}__N{tag}"
colab = man["scenario_results"][key]
refined = os.path.join(RUN_DIR, "arf_refined_events.csv")
tiles_csv = os.path.join(RUN_DIR, "arf_mosaic_tiles.csv")
ov = man["output_validation"]["files"]
for f, p in (("arf_refined_events.csv", refined), ("arf_mosaic_tiles.csv", tiles_csv)):
    assert par.sha256_file(p) == ov[f]["sha256"], f"{p}: not the file the run wrote"
tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])

t0 = time.time()
ns = par.replay(P, refined, tiles, colab["orbital_period_min"], roi, fleet, case, inp["n_satellites"])
mine = ns["scenario_results"][ns["_scenario_key"]()]
secs = round(time.time() - t0, 1)


def flat(d, pre=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{pre}/{k}" if pre else str(k)))
    elif isinstance(d, (list, tuple)):
        for i, v in enumerate(d):
            out.update(flat(v, f"{pre}[{i}]"))
    else:
        out[pre] = d
    return out


def norm(v):
    if hasattr(v, "item"):
        v = v.item()
    if isinstance(v, float) and math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return v


a = {k: norm(v) for k, v in flat(json.loads(json.dumps(colab, default=str))).items()}
b = {k: norm(v) for k, v in flat(json.loads(json.dumps(mine, default=lambda o: o.item() if hasattr(o, "item") else str(o)))).items()}
SKIP = ("seconds", "time_s", "runtime", "elapsed")
common = sorted(set(a) & set(b))
diffs, solver = [], []
for k in common:
    x, y = a[k], b[k]
    if any(s in k for s in SKIP):
        continue
    same = x == y
    if not same and isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        same = math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-12)
    if not same and isinstance(x, str) and isinstance(y, str):
        try:
            same = math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-12)
        except ValueError:
            pass
    if not same:
        (solver if k.startswith("ilp_gap") else diffs).append({"field": k, "colab": x, "sandbox": y})

hv_c = man["hypothesis_validation"]
hv_s = ns.get("hypothesis_validation") or {}
label = [l for l in hv_c["H1_revised"]["assessments"] if l.endswith(f"fleet={tag}")][0]
st = {}
for h in par.STATUSES:
    c = (hv_c.get(h, {}).get("assessments", {}).get(label) or {}).get("status")
    s = list((hv_s.get(h, {}).get("assessments") or {}).values())
    st[h] = {"colab": c, "sandbox": s[0].get("status") if s else None}
res = {"run": os.path.basename(RUN_DIR.rstrip("/")), "scenario": key, "pipeline_sha256": par.sha256_file(P),
       "seconds": secs, "fields_compared": len(common), "fields_only_colab": len(set(a) - set(b)),
       "fields_only_sandbox": len(set(b) - set(a)), "differences": diffs, "solver_dependent_differences": solver,
       "statuses": st, "statuses_equal": all(v["colab"] == v["sandbox"] for v in st.values())}
json.dump(res, open(OUT, "w"), indent=1, default=str)
print(f"{res['run']}: {res['fields_compared']} fields compared in {secs} s; {len(diffs)} differ (+{len(solver)} ILP); "
      f"statuses equal: {res['statuses_equal']}")
for d in diffs[:25]:
    print("   ", d)
for d in solver[:10]:
    print("    [ILP]", d)
