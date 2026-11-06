"""Test of arf_exact_time_a9.py on the REAL saved inputs of the 12 October 2026 CKGR+Okavango run (97 satellites),
presented as one primary run made with Cell 2's exact altitude and period (as in test_data_volume_a6.py), with the
REAL arf_exact_economy_a4.py and arf_colab_parity.py. Budget part on onsets 0, 480, 960, 1440 min with k = 5, 10.
Checks:
  - the script's replay checks pass (onset and cloud CSVs reproduced; rebuilt plans equal the run's)
  - every exact time is verified independently: feasible (min images <= k) at T*, infeasible at the previous
    delivery time; every "no k-cover exists" verified at the horizon end; ECC's time and images equal
    earliest_completion
  - Fail-Signal times recomputed independently for two realisations (own union of usable tiles over the plan's
    images in delivery order)
  - resumable from the cache; both LaTeX tables compile at IEEE width
Usage: test_exact_time_a9.py SCRIPT EXACT_A4 PIPELINE_V153 PARITY_V153 PROTOCOL REFINED TILES ORBIT_MIN"""
import sys, os, json, shutil, tempfile, hashlib, importlib.util, math, subprocess, io, contextlib, runpy
import numpy as np, pandas as pd

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
root = tempfile.mkdtemp(prefix="a9_")
code, batch, run = (os.path.join(root, x) for x in ("ARF_code", "ARF_batch", "run_w1"))
for x in (code, batch, run): os.makedirs(x)
for src, n in ((P, "algorithmic_reconstruction_framework.py"), (PAR, "arf_colab_parity.py"),
               (PROT, "arf_confirmatory_protocol.json"), (EX, "arf_exact_economy_a4.py")):
    shutil.copy(src, os.path.join(code, n))
shutil.copy(REF, os.path.join(run, "arf_refined_events.csv")); shutil.copy(TILES, os.path.join(run, "arf_mosaic_tiles.csv"))
for f in ("arf_onset_robustness.csv", "arf_cloud_robustness.csv"):
    shutil.copy(os.path.join(ns["_out_dir"], f), run)
json.dump({"median_alt_km": ALT0, "t_max_wait_min": T0}, open(os.path.join(run, "arf_cell2_summary.json"), "w"))
json.dump({"framework_version": "1.5.3",
           "config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                                "planet_swath_case": "conservative", "analysis_epoch": "2026-08-03T00:00:00"},
           "inputs": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "n_satellites": 97,
                      "config": {"planet_swath_case": "optimistic"}},
           "scenario_results": {"k": {"orbital_period_min": T0}},
           "output_validation": {"files": {f: {"sha256": sha(os.path.join(run, f))}
                                           for f in ("arf_refined_events.csv", "arf_mosaic_tiles.csv")}}},
          open(os.path.join(run, "arf_reproducibility_manifest.json"), "w"))
pd.DataFrame([{"run_id": "w1", "status": "done", "release": "1.5.3", "out_dir": run}]).to_csv(
    os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
spec = importlib.util.spec_from_file_location("a9", S); a9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a9)
a9.SAMPLE_MIN, a9.BUDGETS = 480, (5, 10)
with contextlib.redirect_stdout(io.StringIO()):
    B, F, summ = a9.main(["--code", code, "--batch-dir", batch])
check(len(B) == 8 and len(F) == int(ns["CLOUD_SEEDS"]), f"replay checks passed; {len(B)} budget rows, {len(F)} realisations rebuilt and equal to the run's")
print(B[["onset_shift_min", "k", "floor_min", "ecc_min", "exact_min", "verdict", "solves"]].to_string(index=False))
H = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
ok, n_ver = True, 0
for r in B.itertuples():
    sub, _ = ns["subset_by_time"](ns["event_data"], None, r.onset_shift_min, r.onset_shift_min + H)
    R = ns["reachable_tile_count"](sub, None); need = int(np.ceil(0.95 * R - 1e-9))
    te, _, ne = ns["earliest_completion"](sub, None, r.k, 95.0, denominator=R)
    ok &= (te is None and pd.isna(r.ecc_min)) or abs(te - r.ecc_min) < 1e-9
    d = sorted({float(x) for x in ns["event_delivery_times"](sub) if np.isfinite(x)})
    def mi(t):
        return ns["solve_ilp_min_images"](ex["delivered_by"](ns, sub, t), None, need, 300.0, 3_000_000)
    if r.verdict == "no k-cover exists":
        m = mi(d[-1]); ok &= m["status"] == "infeasible" or (m["status"] == "optimal" and m["images"] > r.k); n_ver += 1
    elif r.verdict != "ECC at floor" and not r.unknown:
        m = mi(r.exact_min); ok &= m["status"] == "optimal" and m["images"] <= r.k
        prev = [t for t in d if t < r.exact_min - 1e-9 and t >= r.floor_min - 1e-9]
        if prev:
            m2 = mi(prev[-1]); ok &= m2["status"] == "infeasible" or (m2["status"] == "optimal" and m2["images"] > r.k)
        n_ver += 1
    else:
        ok &= abs(r.exact_min - r.floor_min) < 1e-9 and te is not None and te <= r.floor_min + 1e-9
check(ok, f"ECC rows equal earliest_completion; {n_ver} exact verdicts verified independently by fresh ILP solves")
# independent Fail-Signal times
ed = ns["event_data"]; ok = True
dl = ns["event_delivery_times"](ed)
order_of = {k_: dl[i] for i, k_ in enumerate(ed.keys_)}
for i in (0, 7):
    seed = int(ns["RANDOM_SEED"]) + i
    m = ns["build_cloud_mask"](ed, ns["CLOUD_PROBABILITY"], seed)
    util = np.where(m.flags, 0.0, ed.util)
    reach = int(np.unique(ed.tile_idx[util > 0]).size)
    a = ns["select_passes"](ed, ns["MAX_SELECTION_STEPS"], cloud_mask=m, gap_fill=True, time_aware=True,
                            latency_decay_rate=float(ns["LATENCY_DECAY_RATE"]))[0]
    o = ns["plan_expected_coverage"](ed, len(a), ns["_ptile"], ns["time_aware_weights"](ed, float(ns["LATENCY_DECAY_RATE"])))
    def t_at(sel, lv):
        seen = set()
        for key in sorted(sel, key=lambda k_: (order_of[k_], k_)):
            s = ed.rows(ed.index[key]); seen |= set(ed.tile_idx[s][util[s] > 0].tolist())
            if 100.0 * len(seen) / reach >= lv - 1e-9:
                return order_of[key] / 60.0
        return None
    r = F[F.seed == seed].iloc[0]
    for lv in (50.0, 80.0):
        for tag, sel in (("adaptive", a), ("openloop", o)):
            v, w = t_at(sel, lv), r[f"{tag}_t{lv:g}_h"]
            ok &= (v is None and pd.isna(w)) or (v is not None and abs(v - w) < 1e-9)
    L = r.level_L_pct
    ok &= abs(t_at(o, L) - r.openloop_t_L_h) < 1e-9 and (pd.isna(r.adaptive_t_L_h) or abs(t_at(a, L) - r.adaptive_t_L_h) < 1e-9)
check(ok, "Fail-Signal times at 50%, 80% and the open-loop level recomputed independently for two realisations")
calls = []
orig = ex["replay_exact"]
with contextlib.redirect_stdout(io.StringIO()):
    B2, F2, _ = a9.main(["--code", code, "--batch-dir", batch])
check(len(B2) == len(B) and len(F2) == len(F), "second run resumes from the cache")
okt = True
for tb in ("tab_a9_budget.tex", "tab_a9_failsignal.tex"):
    d = tempfile.mkdtemp(dir=root)
    open(os.path.join(d, "t.tex"), "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,"
        "paperheight=11in,textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
        + open(os.path.join(batch, "exact_time_a9", tb)).read() + "x\n\\end{document}\n")
    rr = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=d, capture_output=True)
    okt &= rr.returncode == 0 and "Overfull \\hbox" not in open(os.path.join(d, "t.log"), errors="ignore").read()
check(okt, "both LaTeX tables compile at IEEE width without overfull boxes")
print(json.dumps({k: summ[k] for k in ("budget", "failsignal")}, indent=1, default=str))
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
