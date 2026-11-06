"""Test of arf_wetseason_summary_a5.py: 16 run folders with REAL Cell 6-8 outputs (release 1.5.3 replayed on
the saved 12 October inputs) and the REAL Cell 7D files of window 2026-08-03, given A5 manifests; in the
earthenv folders the cloud record is set to a climatology (stand-in values) and H2 is changed per window so
the statistics are not trivial. Checks: 16 runs read; H2 per window equals the run's own summary; mean, SD,
t lower bounds (alpha, alpha/3) and exact tests equal scipy; observed scene cloud recomputed from the CSV;
the dry-season comparison read when present; the LaTeX table compiles at IEEE width without overfull boxes.
Usage: test_wet_summary_a5.py SUMMARY ANALYSIS_V153 PROTOCOL FIXTURE_7D BASE_RUN"""
import sys, os, json, shutil, tempfile, subprocess, importlib.util, io, contextlib, math
import numpy as np, pandas as pd
from scipy import stats

S, AN, PROT, FIX, BASE = sys.argv[1:6]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


root = tempfile.mkdtemp(prefix="a5sum_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); os.makedirs(batch)
shutil.copy(AN, os.path.join(code, "arf_confirmatory_analysis.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
W = ["2026-01-05", "2026-01-12", "2026-01-19", "2026-01-26", "2026-02-02", "2026-02-09", "2026-02-16", "2026-02-23"]
h2e = [0.4, 1.9, 3.2, -0.3, 2.2, 0.0, 1.1, 2.6]
rows = []
for v in ("p035", "earthenv"):
    for i, w in enumerate(W):
        d = os.path.join(root, v, w)
        shutil.copytree(BASE, d)
        for f in os.listdir(FIX):
            if f.startswith("arf_acquisition_validation"):
                shutil.copy(os.path.join(FIX, f), d)
        cfg = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative",
               "analysis_epoch": w + "T00:00:00", **({"cloud_probability_source": "earthenv"} if v == "earthenv" else {})}
        json.dump({"framework_version": "1.5.3", "cell_versions": {"consistent": True}, "config_overrides": cfg},
                  open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
        if v == "earthenv":
            p = os.path.join(d, "arf_cloud_robustness_summary.json")
            c = json.load(open(p))
            c["cloud_model"]["probability"] = {"source": "earthenv", "mean": 0.30 + 0.01 * i}
            c["h2"]["time_aware"]["difference_vs_openloop_expected"]["mean"] = h2e[i]
            json.dump(c, open(p, "w"))
        rows.append({"run_id": f"{v}_{w}", "status": "done", "release": "1.5.3", "out_dir": d, "minutes": "4",
                     "config": json.dumps(cfg)})
pd.DataFrame(rows).to_csv(os.path.join(batch, "wetseason_v153_a5_index.csv"), index=False)
pd.DataFrame([{"roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative",
               "exclude_satellites": "", "h1a_pct": 100.0, "h1b_pct": x, "h2_diff": y}
              for x, y in zip([100] * 7 + [64.38], [0.85, 2.535, 4.669, 5.041, 0.867, 1.047, 0.206, 3.018])]
             ).to_csv(os.path.join(batch, "arf_confirmatory_results_windows.csv"), index=False)
spec = importlib.util.spec_from_file_location("s", S); s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
with contextlib.redirect_stdout(io.StringIO()):
    df, summ = s.main(["--code", code, "--batch-dir", batch])
check(len(df) == 16 and set(df.variant) == {"p035", "earthenv"}, "16 runs read, two cloud settings")
e = df[df.variant == "earthenv"].sort_values("window")
check(np.allclose(e.h2_diff.values, h2e), "earthenv H2 per window equals the runs' summaries")
x = np.array(h2e); n = len(x); m = x.mean(); sd = x.std(ddof=1)
p = summ["by_variant"]["earthenv"]["h2_diff"]
ok = (abs(p["mean"] - m) < 1e-12 and abs(p["sd"] - sd) < 1e-12 and
      abs(p["lcb_alpha"] - (m - stats.t.ppf(0.95, n - 1) * sd / math.sqrt(n))) < 1e-12 and
      abs(p["lcb_alpha_over_3"] - (m - stats.t.ppf(1 - 0.05 / 3, n - 1) * sd / math.sqrt(n))) < 1e-12)
check(ok, "mean, SD and both Student-t lower bounds equal scipy")
nz = x[x != 0]
check(p["above"] == 6 and p["equal"] == 1 and p["below"] == 1 and
      abs(p["sign_p"] - stats.binomtest(6, 7, 0.5, alternative="greater").pvalue) < 1e-15 and
      abs(p["wilcoxon_p"] - stats.wilcoxon(nz, alternative="greater", method="exact").pvalue) < 1e-15,
      "exact tests: 6 above, 1 equal (dropped), 1 below; p values equal scipy")
sc = pd.read_csv(os.path.join(FIX, "arf_acquisition_validation.csv"))
check(abs(e.observed_cloud_mean_pct.iloc[0] - sc.cloud_pct.mean()) < 1e-12, "observed scene cloud recomputed")
check(abs(summ["dry_season"]["p035"]["h1b_pct"]["mean"] - (700 + 64.38) / 8) < 1e-9, "dry-season comparison read")
tex = os.path.join(root, "t.tex")
open(tex, "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,"
                     "textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
                     + open(os.path.join(batch, "wetseason_a5", "a5_table.tex")).read() + "\\end{document}\n")
r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=root, capture_output=True)
log = open(os.path.join(root, "t.log"), errors="ignore").read()
check(r.returncode == 0 and "Overfull \\hbox" not in log, "a5_table.tex compiles at IEEE width without overfull boxes")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
