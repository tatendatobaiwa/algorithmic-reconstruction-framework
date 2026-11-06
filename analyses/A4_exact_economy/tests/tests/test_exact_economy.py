"""Test of arf_exact_economy_a4.py on the REAL saved inputs of the 12 October 2026 CKGR+Okavango run
(earth-resources+planet, optimistic; 97 satellites), presented as one primary run folder of window 1:
  - the replay reproduces the folder's arf_onset_robustness.csv (made by the same replay procedure)
  - per onset, ECC images, FCFS images and the floor equal the onset CSV
  - invariants: exact minimum <= ECC images; exact minimum <= FCFS images (FCFS reaches the level by the
    floor, ECC finishes at or after it); every solve optimal
  - a tampered refined-events file is refused (SHA-256 check against the manifest)
Usage: test_exact_economy.py SCRIPT PIPELINE_V153 PARITY_V153 PROTOCOL_A3 REFINED TILES ONSET_CSV ORBIT_MIN"""
import sys, os, json, shutil, tempfile, hashlib, importlib.util
import numpy as np, pandas as pd

S, P, PAR, PROT, REF, TILES, ONSET, ORB = sys.argv[1:9]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
root = tempfile.mkdtemp(prefix="exact_")
code, batch, run = (os.path.join(root, x) for x in ("ARF_code", "ARF_batch", "run_w1"))
for x in (code, batch, run):
    os.makedirs(x)
shutil.copy(P, os.path.join(code, "algorithmic_reconstruction_framework.py"))
shutil.copy(PAR, os.path.join(code, "arf_colab_parity.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
proto = json.load(open(PROT))
shutil.copy(REF, os.path.join(run, "arf_refined_events.csv"))
shutil.copy(TILES, os.path.join(run, "arf_mosaic_tiles.csv"))
shutil.copy(ONSET, os.path.join(run, "arf_onset_robustness.csv"))
man = {"framework_version": "1.5.3",
       "config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                            "planet_swath_case": "optimistic", "analysis_epoch": proto["windows"][0]},
       "inputs": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "n_satellites": 97,
                  "config": {"planet_swath_case": "optimistic"}},
       "scenario_results": {"ckgr_okavango__Nearth-resources_planet-optimistic": {"orbital_period_min": float(ORB)}},
       "output_validation": {"files": {f: {"sha256": sha(os.path.join(run, f))}
                                       for f in ("arf_refined_events.csv", "arf_mosaic_tiles.csv")}}}
json.dump(man, open(os.path.join(run, "arf_reproducibility_manifest.json"), "w"))
pd.DataFrame([{"run_id": "w1", "status": "done", "release": "1.5.3", "out_dir": run}]).to_csv(
    os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)

spec = importlib.util.spec_from_file_location("ex", S)
ex = importlib.util.module_from_spec(spec); spec.loader.exec_module(ex)
df, summ = ex.main(["--code", code, "--batch-dir", batch])
check(summ["windows"][0]["replay_identical"], "replay reproduces the folder's arf_onset_robustness.csv")
on = pd.read_csv(ONSET).set_index("onset_shift_min")
d = df.set_index("onset_shift_min")
d = d[d.reachable_tiles > 0]
eq = lambda a, b: np.allclose(pd.to_numeric(a, errors="coerce").values, pd.to_numeric(b, errors="coerce").values,
                              equal_nan=True)
check(len(d) == int((on.reachable_clear_tiles > 0).sum()), f"{len(d)} onsets with reachable tiles")
check(eq(d.ecc_images, on.loc[d.index, "earliest_completion_clear_images95r"]), "ECC images equal the onset CSV")
check(eq(d.fcfs_images, on.loc[d.index, "fcfs_unbudgeted_clear_images95r"]), "FCFS images equal the onset CSV")
check(eq(d.floor_min, on.loc[d.index, "best_possible_clear_t95r_min"]), "floor equals the onset CSV")
check(eq(d.ecc_min, on.loc[d.index, "earliest_completion_clear_t95r_min"]), "ECC time equals the onset CSV")
s = d[d.exact_status == "optimal"]
exp = d[(d.index.values % 80 == 0) | (d.ecc_images >= d.fcfs_images)]
check(set(s.index) == set(exp.index) and (d.exact_status.notna() == d.index.isin(exp.index)).all(),
      f"exact solves exactly at the ties/losses and the 80-min sample, all optimal ({len(s)} onsets)")
check(summ["windows"][0]["sample_solved"] == int((d.index.values % 80 == 0).sum()), "sample count recorded")
check((s.exact_images <= s.ecc_images).all(), "exact minimum <= ECC images at every onset")
check((s.exact_images <= s.fcfs_images).all(), "exact minimum <= FCFS images at every onset")
print(summ["windows"][0])
print(d[["floor_min", "ecc_min", "ecc_images", "fcfs_images", "exact_images", "exact_pairs", "exact_seconds"]].describe().round(2))

# tampered input refused
open(os.path.join(run, "arf_refined_events.csv"), "a").write("\n")
try:
    ex.main(["--code", code, "--batch-dir", batch]); refused = False
except AssertionError as e:
    refused = "differs from what Cell 9 recorded" in str(e)
check(refused, "a changed arf_refined_events.csv is refused")
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
