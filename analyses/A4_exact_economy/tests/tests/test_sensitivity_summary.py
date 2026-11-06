"""Test of arf_sensitivity_summary_a4.py with REAL replayed outputs (release 1.5.3 Cells 6-8 on the saved
12 October 2026 CKGR+Okavango inputs): protocol settings (BASE), k = 10 (K10) and the full-window exact
optimum (ILP48); the earthenv variant is a labelled stand-in (BASE with its cloud probability record
changed), used only to exercise the table code. Each is presented as all eight protocol windows.
Usage: test_sensitivity_summary.py SUMMARY ANALYSIS_V153 PROTOCOL_A3 FIXTURE_7D BASE K10 ILP48"""
import sys, os, json, shutil, tempfile, subprocess, importlib.util
import pandas as pd

S, AN, PROT, FIX, BASE, K10, ILP = sys.argv[1:8]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


root = tempfile.mkdtemp(prefix="a4sum_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); os.makedirs(batch)
shutil.copy(AN, os.path.join(code, "arf_confirmatory_analysis.py"))
shutil.copy(PROT, os.path.join(code, "arf_confirmatory_protocol.json"))
proto = json.load(open(PROT))
design = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative"}


def make(src, name, w, extra):
    d = os.path.join(root, name, w[:10])
    shutil.copytree(src, d)
    for f in os.listdir(FIX):
        if f.startswith("arf_acquisition_validation"):    # Cell 7D only; each run keeps its own ILP file
            shutil.copy(os.path.join(FIX, f), d)
    cfg = {**design, "analysis_epoch": w, **extra}
    json.dump({"framework_version": "1.5.3", "cell_versions": {"consistent": True}, "config_overrides": cfg},
              open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    return {"run_id": f"{name}_{w}", "status": "done", "release": "1.5.3", "out_dir": d, "minutes": "3",
            "config": json.dumps(cfg)}


conf, sens = [], []
for w in proto["windows"]:
    conf.append(make(BASE, "conf", w, {}))
    sens.append(make(K10, "k10", w, {"max_selection_steps": 10}))
    sens.append(make(ILP, "ilp48", w, {"ilp_horizon_hours": 48, "ilp_instances": 1, "ilp_max_pairs": 3000000,
                                       "ilp_time_limit_s": 900.0}))
    r = make(BASE, "earthenv", w, {"cloud_probability_source": "earthenv"})
    p = os.path.join(r["out_dir"], "arf_cloud_robustness_summary.json")
    c = json.load(open(p)); c["cloud_model"]["probability"] = {"source": "earthenv", "mean": 0.08}
    json.dump(c, open(p, "w")); sens.append(r)
pd.DataFrame(conf).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
pd.DataFrame(sens).to_csv(os.path.join(batch, "sensitivity_v153_a4_index.csv"), index=False)

spec = importlib.util.spec_from_file_location("s", S)
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
df, summ = s.main(["--code", code, "--batch-dir", batch])
check(summ["variants"] == {"protocol": 8, "k10": 8, "ilp48": 8, "earthenv": 8}, f"variants read: {summ['variants']}")
an = __import__("runpy").run_path(os.path.join(code, "arf_confirmatory_analysis.py"))
w = an["window_outcomes"](sens[0]["out_dir"], int(proto["min_evaluable_onsets"]))
k = [b for b in summ["budget"] if b["k"] == 10][0]
check(abs(k["h1a_mean"] - w["h1a_pct"]) < 1e-12 and abs(k["h1b_mean"] - w["h1b_pct"]) < 1e-12,
      "k = 10: H1a and H1b equal the protocol reader's window outcomes")
on = pd.read_csv(os.path.join(K10, "arf_onset_robustness.csv"))
on = on[on.reachable_clear_tiles > 0]
gap = (on.earliest_completion_clear_t95r_min - on.best_possible_clear_t95r_min).dropna()
check(abs(k["gap_median_min"] - gap.median()) < 1e-9 and abs(k["gap_max_min"] - gap.max()) < 1e-9,
      f"k = 10: ECC gap above the floor recomputed (median {gap.median():g}, max {gap.max():g} min)")
base = [b for b in summ["budget"] if b["k"] == 30][0]
print("budget rows:", json.dumps(summ["budget"], default=str)[:900])
il = summ["ilp48"][0]
check(il["files_compared"] == 7 and il["files_identical"] == 7, "ilp48 replay: the 7 result files equal the protocol run's")
check(str(il["horizon_min"]) == "2880" and il["ecc_exact_status"] == "optimal", "ilp48: one 48-h instance, solved to optimality")
check(len(summ["cloud"]) == 8 and summ["cloud"][0]["p_earthenv_mean"] == 0.08, "cloud table built from the variant's record")
tex = os.path.join(root, "t.tex")
open(tex, "w").write("\\documentclass[10pt,twocolumn]{article}\n\\usepackage[paperwidth=8.5in,paperheight=11in,"
                     "textwidth=7.16in,columnsep=0.16in,textheight=9.5in]{geometry}\n\\begin{document}\n"
                     + open(os.path.join(batch, "sensitivity_a4", "a4_tables.tex")).read() + "\n\\end{document}\n")
r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "t.tex"], cwd=root, capture_output=True, text=True)
log = open(os.path.join(root, "t.log"), errors="ignore").read()
check(r.returncode == 0 and "Overfull \\hbox" not in log, "a4_tables.tex compiles at IEEE width without overfull boxes")
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
