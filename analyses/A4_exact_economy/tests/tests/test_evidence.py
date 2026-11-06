"""Tests of arf_paper1_evidence.py on run folders with REAL Cell 6-8 outputs (Cells 6-8 of release 1.5.3
replayed on the saved inputs of 12 October 2026) and the REAL Cell 7D files of window 2026-08-03 (downloaded
from the primary run folder), each given a manifest for one protocol window.
Checks: every window read; per-window values equal the source files; the protocol reader's H2 matches;
exact tests equal an independent computation; LaTeX tables compile (pdflatex).
Usage: test_evidence.py EVIDENCE ANALYSIS_V153 PROTOCOL_A3 FIXTURE_7D RUN_DIR x8"""
import sys, os, json, shutil, tempfile, subprocess, importlib.util
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon

EV, AN, PROT, FIX7D, *RUNS = sys.argv[1:]
assert len(RUNS) == 8
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


root = tempfile.mkdtemp(prefix="p1test_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); os.makedirs(batch)
shutil.copy(AN, os.path.join(code, "arf_confirmatory_analysis.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
proto = json.load(open(PROT))
rows = []
for w, src in zip(proto["windows"], RUNS):
    d = os.path.join(root, "run_" + w[:10])
    shutil.copytree(src, d)
    for f in os.listdir(FIX7D):
        shutil.copy(os.path.join(FIX7D, f), d)
    json.dump({"framework_version": "1.5.3", "cell_versions": {"consistent": True},
               "config_overrides": {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                                    "planet_swath_case": "conservative", "analysis_epoch": w}},
              open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    rows.append({"run_id": w, "status": "done", "release": "1.5.3", "out_dir": d, "minutes": "3.5"})
# one run of another release and one failed run: both must be ignored
rows.append({"run_id": "other", "status": "done", "release": "1.5.2", "out_dir": d, "minutes": "1"})
rows.append({"run_id": "err", "status": "error", "release": "1.5.3", "out_dir": "", "minutes": "1"})
pd.DataFrame(rows).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)

spec = importlib.util.spec_from_file_location("ev", EV)
ev = importlib.util.module_from_spec(spec); spec.loader.exec_module(ev)
df, summ = ev.main(["--code", code, "--batch-dir", batch])
out = os.path.join(batch, "paper1_evidence")
check(len(df) == 8 and sorted(df.window) == [w[:10] for w in proto["windows"]], "8 windows read, others ignored")
ok = True
for r in df.itertuples():
    cr = json.load(open(os.path.join(r.run_dir, "arf_cloud_robustness_summary.json")))["h2"]["time_aware"]
    ok &= abs(r.diff_vs_expected_mean - cr["difference_vs_openloop_expected"]["mean"]) < 1e-12
    ok &= abs(r.cov_oracle - cr["coverage_pct_mean"]["oracle"]) < 1e-12
    g = pd.read_csv(os.path.join(r.run_dir, "arf_ilp_gap.csv"))
    e = g[g.ecc_exact_status == "optimal"]
    ok &= r.ecc_matches_exact == int((e.ecc_greedy_images == e.ecc_exact_images).sum())
    on = pd.read_csv(os.path.join(r.run_dir, "arf_onset_robustness.csv"))
    on = on[on.reachable_clear_tiles > 0]
    ok &= abs(r.earliest_completion_t_median_min - on.earliest_completion_clear_t95r_min.dropna().median()) < 1e-9
    ok &= abs(r.fcfs_unbudgeted_images_median - on.fcfs_unbudgeted_clear_images95r.dropna().median()) < 1e-9
check(ok, "per-window values equal the source files (H2, oracle, ILP matches, ECC median, FCFS images)")
x = df.h2_diff.astype(float).values
nz = x[x != 0]
ref_sign = binomtest(int((nz > 0).sum()), nz.size, 0.5, alternative="greater").pvalue
ref_w = wilcoxon(nz, alternative="greater", method="exact").pvalue
t = summ["tests"]["H2_vs_zero"]
check(abs(t["sign_p"] - ref_sign) < 1e-15 and abs(t["wilcoxon_p"] - ref_w) < 1e-15, "H2 exact tests equal scipy directly")
if (nz > 0).all() and nz.size == 8:
    check(abs(t["sign_p"] - 1 / 256) < 1e-15 and abs(t["wilcoxon_p"] - 1 / 256) < 1e-15,
          "all 8 windows positive: p = 1/256 for both tests")
h1a = summ["tests"]["H1a_vs_threshold"]
check(h1a["above"] + h1a["equal"] + h1a["below"] == 8, "H1a counts add to 8")
sc = pd.read_csv(os.path.join(out, "p1_archive_scenes.csv"))
check(len(sc) == 24 and abs(summ["archive_cloud"]["mean_pct"] - sc.cloud_pct.mean()) < 1e-12,
      "archive scenes pooled (3 per window) and observed cloud mean recomputed")
pp = pd.read_csv(os.path.join(out, "p1_archive_pooled.csv")).set_index("platform")
check(pp.loc["sentinel-2a", "matched"] == 8 and pp.loc["landsat-8", "predicted"] == 0,
      "archive pooled per platform (S2A matched 8 of 8; Landsat-8 none predicted)")
# LaTeX
tex = os.path.join(root, "t.tex")
open(tex, "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\usepackage{amsmath}\n\\begin{document}\n"
                     + open(os.path.join(out, "p1_tables.tex")).read() + "\n\\end{document}\n")
r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=root,
                   capture_output=True, text=True)
log = open(os.path.join(root, "t.log"), errors="ignore").read()
check(r.returncode == 0, "p1_tables.tex compiles (two columns of 3.5 in, IEEE text width 7.16 in)")
check("Overfull \\hbox" not in log, "no overfull boxes")
print(open(os.path.join(out, "p1_tables.tex")).read()[:1800])
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
