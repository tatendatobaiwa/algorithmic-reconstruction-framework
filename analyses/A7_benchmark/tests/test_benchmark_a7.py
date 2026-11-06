"""Test of arf_benchmark_a7.py on run folders holding REAL Cell 6-8 outputs (release 1.5.3 replayed on the saved
12 October inputs, 8 folders) and the REAL arf_mosaic_tiles.csv of that run, given manifests for three arms
(primary, Lower Zambezi, wet season A5) plus one run of another release and one arm not in the list.
Checks: arms and windows read, the others ignored; per-window values recomputed independently from the onset CSV;
bootstrap interval reproducible and containing the mean; cloud values from arf_cloud_robustness.csv; both LaTeX
tables compile at IEEE width without overfull boxes.
Usage: test_benchmark_a7.py SCRIPT ANALYSIS_V153 PROTOCOL TILES RUN_DIR x8"""
import sys, os, json, shutil, tempfile, subprocess, importlib.util, io, contextlib
import numpy as np, pandas as pd

S, AN, PROT, TILES, *RUNS = sys.argv[1:]
assert len(RUNS) == 8
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


root = tempfile.mkdtemp(prefix="a7_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); os.makedirs(batch)
shutil.copy(AN, os.path.join(code, "arf_confirmatory_analysis.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
proto = json.load(open(PROT))
WET = ["2026-01-05T00:00:00", "2026-01-12T00:00:00"]


def mk(name, src, cfg, release="1.5.3"):
    d = os.path.join(root, name)
    shutil.copytree(src, d)
    shutil.copy(TILES, os.path.join(d, "arf_mosaic_tiles.csv"))
    json.dump({"framework_version": release, "cell_versions": {"consistent": True}, "config_overrides": cfg},
              open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    return d


base = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative"}
prim, sec, wet = [], [], []
for i, (w, src) in enumerate(zip(proto["windows"], RUNS)):
    prim.append({"run_id": f"p{i}", "status": "done", "release": "1.5.3",
                 "out_dir": mk(f"p_{w[:10]}", src, {**base, "analysis_epoch": w})})
    sec.append({"run_id": f"lz{i}", "status": "done", "release": "1.5.3",
                "out_dir": mk(f"lz_{w[:10]}", RUNS[(i + 3) % 8], {**base, "active_roi": "lower_zambezi", "analysis_epoch": w})})
sec.append({"run_id": "old", "status": "done", "release": "1.5.3",
            "out_dir": mk("old", RUNS[0], {**base, "analysis_epoch": proto["windows"][0]}, release="1.5.2")})
sec.append({"run_id": "other", "status": "done", "release": "1.5.3",
            "out_dir": mk("other", RUNS[1], {**base, "fleet": "planet", "analysis_epoch": proto["windows"][0]})})
for i, w in enumerate(WET):
    wet.append({"run_id": f"w{i}", "status": "done", "release": "1.5.3",
                "out_dir": mk(f"wet_{w[:10]}", RUNS[i], {**base, "analysis_epoch": w})})
for i, w in enumerate(WET):     # the A5 earthenv runs of the same windows: must not count twice
    wet.append({"run_id": f"we{i}", "status": "done", "release": "1.5.3",
                "out_dir": mk(f"wete_{w[:10]}", RUNS[i], {**base, "analysis_epoch": w, "cloud_probability_source": "earthenv"})})
pd.DataFrame(prim).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
pd.DataFrame(sec).to_csv(os.path.join(batch, "confirmatory_v153_secondary_index.csv"), index=False)
pd.DataFrame(wet).to_csv(os.path.join(batch, "wetseason_v153_a5_index.csv"), index=False)

spec = importlib.util.spec_from_file_location("a7", S); a7 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a7)
with contextlib.redirect_stdout(io.StringIO()):
    W, summ = a7.main(["--code", code, "--batch-dir", batch])
labs = {l for l, *_ in a7.ARM_ORDER}
P, LZ, WT = a7.ARM_ORDER[0][0], "Lower Zambezi", a7.ARM_ORDER[-1][0]
check(set(W.arm) == {P, LZ, WT} and (W.arm == P).sum() == 8 and (W.arm == LZ).sum() == 8 and (W.arm == WT).sum() == 2,
      "3 arms read (8 + 8 + 2 windows); other release, unlisted fleet and earthenv duplicates ignored")
n = 12231 if TILES else None
n = len(pd.read_csv(TILES, usecols=["tile"]))
ok = True
for r in W.itertuples():
    on = pd.read_csv(os.path.join(r.run_dir, "arf_onset_robustness.csv"))
    on = on[on.reachable_clear_tiles > 0]
    ok &= abs(r.floor_h_median - (on.best_possible_clear_t95r_min / 60).median()) < 1e-12
    ok &= abs(r.reach_pct_median - (100 * on.reachable_clear_tiles / n).median()) < 1e-9
    ok &= r.ecc_images_median == on.earliest_completion_clear_images95r.median()
    ok &= r.coverage_first_images_median == on.coverage_first_clear_images95r.median()
    ok &= r.fcfs_unbudgeted_images_median == on.fcfs_unbudgeted_clear_images95r.median()
    ok &= abs(r.ecc_saving_vs_coverage_first_mean - (on.coverage_first_clear_images95r - on.earliest_completion_clear_images95r).mean()) < 1e-12
    ok &= abs(r.ecc_gap_min_max - (on.earliest_completion_clear_t95r_min - on.best_possible_clear_t95r_min).max()) < 1e-12
    c = pd.read_csv(os.path.join(r.run_dir, "arf_cloud_robustness.csv"))
    ok &= abs(r.cloud_floor_h_median - (c.floor_t95r_min / 60).median()) < 1e-12
check(ok, "per-window floor, reachable share, images, savings, gap and cloud floor recomputed independently")
s = summ["arms"][P]["floor_h_median"]
x = W[W.arm == P].floor_h_median.values
lo, hi = a7.boot_ci(x)
check(abs(s["mean_over_windows"] - x.mean()) < 1e-12 and s["ci95_window_bootstrap"] == [lo, hi] and lo <= x.mean() <= hi,
      "window-block bootstrap: reproducible with the fixed seed and contains the mean")
check(summ["cloud_source_variants_skipped"] == 2 and summ["arms"][WT]["windows"] == 2 and summ["arms"][P]["onsets"] == int(W[W.arm == P].onsets.sum()),
      "wet-season arm and onset totals reported")
ok = True
for t in ("tab_a7_floor.tex", "tab_a7_images.tex"):
    d = tempfile.mkdtemp(dir=root)
    open(os.path.join(d, "t.tex"), "w").write(
        "\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,textwidth=7.16in,"
        "columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
        + open(os.path.join(batch, "benchmark_a7", t)).read() + "x\n\\end{document}\n")
    rr = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=d, capture_output=True)
    ok &= rr.returncode == 0 and "Overfull \\hbox" not in open(os.path.join(d, "t.log"), errors="ignore").read()
check(ok, "both LaTeX tables compile at IEEE width without overfull boxes")
print(W[W.arm == P][["window", "floor_h_median", "reach_pct_median", "ecc_images_median",
                     "coverage_first_images_median", "fcfs_unbudgeted_images_median"]].round(2).to_string(index=False))
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
