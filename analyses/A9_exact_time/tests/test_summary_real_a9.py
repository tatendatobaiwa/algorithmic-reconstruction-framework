"""Summary of arf_exact_time_a9.py on the REAL A9 results of 30 October 2026 (a9_budget.csv, a9_failsignal.csv from
Drive), fed through the script's own cache: no replay, no solve. Checks the corrected Fail-Signal summary against an
independent computation in which a realisation whose adaptive plan never reaches the open-loop level ranks as worst,
and that both tables compile. Usage: test_summary_real_a9.py SCRIPT EXACT_A4 PIPELINE PARITY PROTOCOL BUDGET FAIL"""
import sys, os, json, shutil, tempfile, importlib.util, io, contextlib, subprocess
import numpy as np, pandas as pd
S, EX, P, PAR, PROT, BUD, FAIL = sys.argv[1:8]
RES = []
def check(c, m): RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m)
B, F = pd.read_csv(BUD), pd.read_csv(FAIL)
root = tempfile.mkdtemp(); code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); cache = os.path.join(batch, "exact_time_a9", "cache"); os.makedirs(cache)
for src, n in ((P, "algorithmic_reconstruction_framework.py"), (PAR, "arf_colab_parity.py"), (EX, "arf_exact_economy_a4.py"),
               (PROT, "arf_confirmatory_protocol.json")):
    shutil.copy(src, os.path.join(code, n))
rows = {"dry": [], "wet": []}
nn = lambda r: {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in r.items()}
for (s, w), g in F.groupby(["season", "window"]):
    d = os.path.join(root, w); os.makedirs(d)
    json.dump({"config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
               "planet_swath_case": "conservative", "analysis_epoch": w + "T00:00:00"}},
              open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    rows[s].append({"run_id": w, "status": "done", "release": "1.5.3", "out_dir": d})
    json.dump({"budget": [nn(r) for r in B[B.window == w].to_dict("records")],
               "fail": [nn(r) for r in g.to_dict("records")],
               "info": {"season": s, "window": w}}, open(os.path.join(cache, f"{s}_{w}.json"), "w"))
pd.DataFrame(rows["dry"]).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
pd.DataFrame(rows["wet"]).to_csv(os.path.join(batch, "wetseason_v153_a5_index.csv"), index=False)
spec = importlib.util.spec_from_file_location("a9", S); a9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a9)
with contextlib.redirect_stdout(io.StringIO()):
    _, _, summ = a9.main(["--code", code, "--batch-dir", batch])
ok = True
for s, g in F.groupby("season"):
    x = np.where(g.dT_L_h.notna(), g.dT_L_h, -np.inf)
    v = summ["failsignal"][s]
    ok &= v["adaptive_never_reaches_openloop_level"] == int(g.dT_L_h.isna().sum())
    ok &= abs(v["dT_L_h_median"] - np.median(x)) < 1e-12
    ok &= abs(v["pct_adaptive_later_or_never"] - 100 * np.mean(x < -1e-9)) < 1e-9
    ok &= abs(v["pct_adaptive_sooner"] - 100 * np.mean(x > 1e-9)) < 1e-9
check(ok, "corrected Fail-Signal summary equals an independent computation (never-reaching ranked worst)")
check(summ["budget"]["wet_k15"]["verdicts"].get("ECC missed a cover") == 7 and sum(s["unknown"] for s in summ["budget"].values()) == 0,
      "budget summary unchanged: 7 greedy misses (wet, k = 15), no unknown")
okt = True
for tb in ("tab_a9_budget.tex", "tab_a9_failsignal.tex"):
    d = tempfile.mkdtemp(dir=root)
    open(os.path.join(d, "t.tex"), "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,"
        "paperheight=11in,textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
        + open(os.path.join(batch, "exact_time_a9", tb)).read() + "x\n\\end{document}\n")
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=d, capture_output=True)
    okt &= r.returncode == 0 and "Overfull \\hbox" not in open(os.path.join(d, "t.log"), errors="ignore").read()
check(okt, "both tables compile at IEEE width without overfull boxes")
shutil.copy(os.path.join(batch, "exact_time_a9", "tab_a9_failsignal.tex"), "/home/claude/arf/a9results/tab_a9_failsignal.tex")
shutil.copy(os.path.join(batch, "exact_time_a9", "tab_a9_budget.tex"), "/home/claude/arf/a9results/tab_a9_budget.tex")
print(json.dumps(summ["failsignal"], indent=1))
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
