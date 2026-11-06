"""Test of arf_exact_economy_a4.py (revised 26 October 2026) on the REAL saved inputs of the 12 October 2026
CKGR+Okavango run. The "Colab run" is made as Cell 2 makes it: a median altitude ALT0 and the period
T0 = 2 pi sqrt((Re + ALT0)^3 / mu) / 60 computed by the pipeline's own expression; Cells 6-8 are replayed
with exactly these, and their onset CSV is the run's.
  1 negative control: arf_colab_parity.replay (period rebuilt from the rounded value) does NOT reproduce
    the run's onset CSV byte for byte; which columns differ, and by how much (the Colab failure)
  2 the revised script, which reads ALT0 and T0 from arf_cell2_summary.json, reproduces it byte for byte
  3 compare_onsets: a 1e-12 relative change in a time-weighted column is accepted and reported; a
    change of one image in ECC's count is refused
  4 per onset, ECC images, FCFS images, floor and ECC time equal the CSV; exact <= ECC and <= FCFS
Usage: test_exact_economy_v2.py SCRIPT PIPELINE_V153 PARITY_V153 PROTOCOL_A3 REFINED TILES ORBIT_MIN"""
import sys, os, json, shutil, tempfile, hashlib, importlib.util, math
import numpy as np, pandas as pd

S, P, PAR, PROT, REF, TILES, ORB = sys.argv[1:8]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
spec = importlib.util.spec_from_file_location("ex", S); ex = importlib.util.module_from_spec(spec); spec.loader.exec_module(ex)
spec = importlib.util.spec_from_file_location("par", PAR); par = importlib.util.module_from_spec(spec); spec.loader.exec_module(par)
mu_earth, Re_km = 398600.4418, 6378.137                      # the pipeline's constants (Cell 2)
ALT0 = (mu_earth * (float(ORB) * 60.0 / (2 * math.pi)) ** 2) ** (1.0 / 3.0) - Re_km + 0.123456789
T0 = 2 * np.pi * np.sqrt((Re_km + ALT0) ** 3 / mu_earth) / 60.0   # Cell 2, verbatim
T0 = float(T0)
tiles = list(pd.read_csv(TILES, usecols=["tile"])["tile"])
args = (P, REF, tiles, T0, "ckgr_okavango", "earth-resources+planet", "optimistic", 97)
run_ns = ex.replay_exact(par, ALT0, T0, *args)               # the "Colab run"
check(run_ns["T_MAX_WAIT_MINUTES"] == T0 and run_ns["CONSTELLATION_ALTITUDE_KM"] == ALT0,
      "run made with the exact altitude and period")

root = tempfile.mkdtemp(prefix="exact2_")
code, batch, run = (os.path.join(root, x) for x in ("ARF_code", "ARF_batch", "run_w1"))
for x in (code, batch, run):
    os.makedirs(x)
shutil.copy(P, os.path.join(code, "algorithmic_reconstruction_framework.py"))
shutil.copy(PAR, os.path.join(code, "arf_colab_parity.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
proto = json.load(open(PROT))
shutil.copy(REF, os.path.join(run, "arf_refined_events.csv"))
shutil.copy(TILES, os.path.join(run, "arf_mosaic_tiles.csv"))
shutil.copy(os.path.join(run_ns["_out_dir"], "arf_onset_robustness.csv"), run)
json.dump({"median_alt_km": ALT0, "t_max_wait_min": T0}, open(os.path.join(run, "arf_cell2_summary.json"), "w"))
man = {"framework_version": "1.5.3",
       "config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                            "planet_swath_case": "optimistic", "analysis_epoch": proto["windows"][0]},
       "inputs": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "n_satellites": 97,
                  "config": {"planet_swath_case": "optimistic"}},
       "scenario_results": {"k": {"orbital_period_min": T0}},
       "output_validation": {"files": {f: {"sha256": sha(os.path.join(run, f))}
                                       for f in ("arf_refined_events.csv", "arf_mosaic_tiles.csv")}}}
json.dump(man, open(os.path.join(run, "arf_reproducibility_manifest.json"), "w"))
pd.DataFrame([{"run_id": "w1", "status": "done", "release": "1.5.3", "out_dir": run}]).to_csv(
    os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
run_csv = os.path.join(run, "arf_onset_robustness.csv")

# 1 negative control: the parity replay as the first version used it
old = par.replay(*args)
ident, ok, rep = ex.compare_onsets(os.path.join(old["_out_dir"], "arf_onset_robustness.csv"), run_csv)
print("   parity replay: period", repr(old["T_MAX_WAIT_MINUTES"]), "vs run", repr(T0))
print("   " + "\n   ".join(rep))
check(not ident, "negative control: the parity replay's onset CSV is NOT byte-identical (the Colab failure)")

# 2 the revised script
df, summ = ex.main(["--code", code, "--batch-dir", batch])
check(summ["windows"][0]["replay_identical"] is True, "revised script: replay byte-identical to the run")

# 3 compare_onsets tolerance rules
a = pd.read_csv(run_csv)
num = [c for c in a.columns if c not in ex.EXACT_COLUMNS and a[c].dtype.kind == "f" and a[c].abs().max() > 0]
p1 = os.path.join(root, "p1.csv"); b = a.copy(); b[num[0]] = b[num[0]] * (1 + 1e-12); b.to_csv(p1, index=False)
i1, ok1, r1 = ex.compare_onsets(p1, run_csv)
check(not i1 and ok1 and any(num[0] in x for x in r1), f"1e-12 relative change in {num[0]} accepted and reported")
p2 = os.path.join(root, "p2.csv"); b = a.copy()
b.loc[0, "earliest_completion_clear_images95r"] += 1; b.to_csv(p2, index=False)
i2, ok2, r2 = ex.compare_onsets(p2, run_csv)
check(not ok2 and any("must be exact" in x and "FAIL" in x for x in r2), "one extra ECC image refused")
# a missing value counts as a difference
p3 = os.path.join(root, "p3.csv"); b = a.copy(); b.loc[0, num[0]] = np.nan; b.to_csv(p3, index=False)
check(not ex.compare_onsets(p3, run_csv)[1], "a value turned missing is refused")

# 4 per-onset values and invariants
on = a.set_index("onset_shift_min"); d = df.set_index("onset_shift_min"); d = d[d.reachable_tiles > 0]
eq = lambda x, y: np.allclose(pd.to_numeric(x, errors="coerce").values, pd.to_numeric(y, errors="coerce").values, equal_nan=True)
check(eq(d.ecc_images, on.loc[d.index, "earliest_completion_clear_images95r"]) and
      eq(d.fcfs_images, on.loc[d.index, "fcfs_unbudgeted_clear_images95r"]) and
      eq(d.floor_min, on.loc[d.index, "best_possible_clear_t95r_min"]) and
      eq(d.ecc_min, on.loc[d.index, "earliest_completion_clear_t95r_min"]), "per-onset values equal the run's CSV")
s = d[d.exact_status == "optimal"]
check(len(s) == int(d.exact_status.notna().sum()) and (s.exact_images <= s.ecc_images).all()
      and (s.exact_images <= s.fcfs_images).all(), f"{len(s)} exact solves optimal; exact <= ECC and <= FCFS")
# period mismatch between Cell 2 and the manifest is refused
man["scenario_results"]["k"]["orbital_period_min"] = T0 + 1e-9
json.dump(man, open(os.path.join(run, "arf_reproducibility_manifest.json"), "w"))
try:
    ex.main(["--code", code, "--batch-dir", batch]); refused = False
except AssertionError as e:
    refused = "differs from the manifest" in str(e)
check(refused, "a Cell 2 period that differs from the manifest is refused")
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
