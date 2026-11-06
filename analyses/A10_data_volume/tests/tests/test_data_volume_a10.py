"""Test of arf_data_volume_a6.py on the REAL saved inputs of the 12 October 2026 CKGR+Okavango run (97
satellites), presented as one primary run made with Cell 2's exact altitude and period (as in
test_exact_economy_v2.py), with the REAL arf_exact_economy_a4.py and arf_colab_parity.py.
  - the script's own checks pass (replay acceptable; image counts equal the run's onset CSV)
  - for three onsets, the images and volumes are recomputed independently from the replay: ECC's images
    by earliest_completion, FCFS's by select_fcfs in delivery order up to the level, each usable tile's own h3.cell_area
  - every satellite maps to a sensor rate; low <= high; saving = FCFS - ECC
  - a10_table.tex compiles at IEEE width without overfull boxes
Usage: test_data_volume_a6.py SCRIPT EXACT_A4 PIPELINE_V153 PARITY_V153 PROTOCOL REFINED TILES ORBIT_MIN"""
import sys, os, json, shutil, tempfile, hashlib, importlib.util, math, subprocess, io, contextlib, runpy
import numpy as np, pandas as pd, h3

S, EX, P, PAR, PROT, REF, TILES, ORB = sys.argv[1:9]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
ex = runpy.run_path(EX, run_name="x")
par = ex["load"](PAR, "arf_colab_parity")
mu, Re = 398600.4418, 6378.137
ALT0 = (mu * (float(ORB) * 60.0 / (2 * math.pi)) ** 2) ** (1.0 / 3.0) - Re + 0.123456789
T0 = float(2 * np.pi * np.sqrt((Re + ALT0) ** 3 / mu) / 60.0)
tiles = list(pd.read_csv(TILES, usecols=["tile"])["tile"])
args = (P, REF, tiles, T0, "ckgr_okavango", "earth-resources+planet", "optimistic", 97)
ns = ex["replay_exact"](par, ALT0, T0, *args)

root = tempfile.mkdtemp(prefix="a6_")
code, batch, run = (os.path.join(root, x) for x in ("ARF_code", "ARF_batch", "run_w1"))
for x in (code, batch, run):
    os.makedirs(x)
for src, n in ((P, "algorithmic_reconstruction_framework.py"), (PAR, "arf_colab_parity.py"),
               (PROT, "arf_confirmatory_protocol.json"), (EX, "arf_exact_economy_a4.py")):
    shutil.copy(src, os.path.join(code, n))
proto = json.load(open(PROT))
shutil.copy(REF, os.path.join(run, "arf_refined_events.csv"))
shutil.copy(TILES, os.path.join(run, "arf_mosaic_tiles.csv"))
shutil.copy(os.path.join(ns["_out_dir"], "arf_onset_robustness.csv"), run)
json.dump({"median_alt_km": ALT0, "t_max_wait_min": T0}, open(os.path.join(run, "arf_cell2_summary.json"), "w"))
json.dump({"framework_version": "1.5.3",
           "config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                                "planet_swath_case": "optimistic", "analysis_epoch": proto["windows"][0]},
           "inputs": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "n_satellites": 97,
                      "config": {"planet_swath_case": "optimistic"}},
           "scenario_results": {"k": {"orbital_period_min": T0}},
           "output_validation": {"files": {f: {"sha256": sha(os.path.join(run, f))}
                                           for f in ("arf_refined_events.csv", "arf_mosaic_tiles.csv")}}},
          open(os.path.join(run, "arf_reproducibility_manifest.json"), "w"))
pd.DataFrame([{"run_id": "w1", "status": "done", "release": "1.5.3", "out_dir": run}]).to_csv(
    os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
spec = importlib.util.spec_from_file_location("v", S); v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
with contextlib.redirect_stdout(io.StringIO()):
    df, summ = v.main(["--code", code, "--batch-dir", batch])
check(len(df) == 73, f"73 onsets evaluated, the script's own replay and count checks passed")
# independent recomputation for three onsets
H = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
ok = True
for d in (0, 560, 1440):
    sub, _ = ns["subset_by_time"](ns["event_data"], None, d, d + H)
    R = ns["reachable_tile_count"](sub, None); need = int(np.ceil(0.95 * R - 1e-9))
    _, ke, _ = ns["earliest_completion"](sub, None, 30, 95.0, denominator=R)
    order = sorted(ns["select_fcfs"](sub, None, len(sub))[0], key=ns["delivery_sort_key"](sub))
    util = ns["_effective_util"](sub, None)
    cov, kf = set(), []
    for key in order:
        rows = sub.rows(sub.index[key]); t = set(sub.tile_idx[rows][util[rows] > 0].tolist())
        kf.append(key); cov |= t
        if len(cov) >= need:
            break
    def vol(keys, j):
        tot = 0.0
        for key in keys:
            rows = sub.rows(sub.index[key])
            use = util[rows] > 0
            tot += sum(h3.cell_area(str(ns["_tile_ids"][i]), "km^2") for i in sub.tile_idx[rows][use]) * v.RATES[v.sensor(key[0])][j]
        return tot
    r = df.set_index("onset_shift_min").loc[d]
    ok &= len(kf) == r.fcfs_images and len(ke) == r.ecc_images
    ok &= abs(vol(ke, 0) - r.ecc_mb_low) < 1e-6 and abs(vol(ke, 1) - r.ecc_mb_high) < 1e-6
    ok &= abs(vol(kf, 0) - r.fcfs_mb_low) < 1e-6 and abs(vol(kf, 1) - r.fcfs_mb_high) < 1e-6
check(ok, "images and volumes at onsets 0, 560, 1440 recomputed independently (FCFS by tile union)")
names = pd.read_csv(REF, usecols=["sat"]).sat.unique()
check(all(v.sensor(n) in v.RATES for n in names), f"all {len(names)} satellites map to a sensor rate")
check((df.ecc_mb_low <= df.ecc_mb_high).all() and np.allclose(df.saving_mb_low, df.fcfs_mb_low - df.ecc_mb_low),
      "low <= high; saving = FCFS - ECC")
print(df[["onset_shift_min", "ecc_images", "fcfs_images", "ecc_mb_low", "ecc_mb_high", "fcfs_mb_low", "fcfs_mb_high"]].describe().round(1))
tex = os.path.join(root, "t.tex")
open(tex, "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,"
                     "textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
                     + open(os.path.join(batch, "data_volume_a10", "a10_table.tex")).read() + "\\end{document}\n")
r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=root, capture_output=True)
log = open(os.path.join(root, "t.log"), errors="ignore").read()
check(r.returncode == 0 and "Overfull \\hbox" not in log, "a10_table.tex compiles at IEEE width without overfull boxes")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
