"""Release 1.5.3 presentation changes, on one end-to-end run (e2e15 harness:
synthetic Earth-observation fleet over the real Lower Zambezi tiles; Cells 3-9
from the pipeline file) plus direct calls. Usage: test_presentation_v153.py PIPELINE

Figure 3   y-axis spans the plotted values (not from 0 h) for curves at 24-26 h;
           legend handles returned, no legend inside the axes; the floor drawn
           over the curves
Figure 5   at the H1a level of reachable tiles, hours on a linear axis
Figure 9   minutes on a linear axis, gap per selector written on the panel
Figure 7   at most 7.16 in wide; short panel (b) title, cloud settings in the footer
Cell 8     superseded hypotheses labelled (console, report, JSON); this run's
           statuses first (console, report, JSON); with two scenarios in the
           ledger the overall statuses say they pool them
Cell 9     pipeline log labels the superseded hypotheses"""
import sys, os, io, re, json, contextlib, logging
import numpy as np
sys.path.insert(0, "/home/claude/validate")
import e2e15
import pymupdf
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = sys.argv[1]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


src = open(P, encoding="utf-8").read()
C = e2e15.cells_of(src)

# --- Figure 3: _draw_latency on its own (Cell 7's definitions, no run) ----------------------------
import ast
tree = ast.parse(C["7"])
keep = {"_grid", "_draw_latency"}
tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in keep]
n7 = {"np": np, "COLORS": {"clear_sky": "g", "level": "m", "h1_target": "y", "grid": "0.8"}}
exec(compile(tree, "cell7_fns", "exec"), n7)
fig, ax = plt.subplots()
sweep = [{"k": k, "latency_clear_min": (24.2 if k < 10 or k >= 20 else 25.1) * 60} for k in range(5, 55, 5)]
lvl = {k: 24.17 * 60 for k in range(10, 55, 5)}
h, l = n7["_draw_latency"](ax, sweep, lvl, 24.17 * 60, 94.0, 95.0)
lo, hi = ax.get_ylim()
z_floor = [ln.get_zorder() for ln in ax.get_lines() if ln.get_label().startswith("Physical floor")]
z_curve = [ln.get_zorder() for ln in ax.get_lines() if ln.get_label().startswith("95% of reachable")]
check(23.9 < lo < 24.17 and 25.73 < hi < 26.1 and ax.get_legend() is None and len(h) == 4
      and z_floor and z_curve and z_floor[0] > z_curve[0],
      f"Figure 3: y-axis {lo:.2f}-{hi:.2f} h for values 24.17-25.74 h (not from 0); 4 legend entries returned, "
      f"none inside the axes; floor drawn over the level curve (zorder {z_floor} > {z_curve})")
sweep2 = sweep[:2] + [{"k": 15}] + sweep[3:]
fig2, ax2 = plt.subplots()
n7["_draw_latency"](ax2, sweep2, lvl, 24.17 * 60, 94.0, 95.0)
check(ax2.get_ylim()[0] > 20, "Figure 3: a budget without a completion time is left out, not drawn at 0 h")
plt.close("all")

# --- end-to-end run -----------------------------------------------------------------------------
tiles = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
ns = e2e15.run(P, tiles, overrides={"window_hours": 48, "onset_shift_hours": 24, "onset_shift_step_min": 120},
               cells=("3", "4", "5", "6", "6B", "7", "7B", "7C", "8", "9"))
out = ns["_out_dir"]


def pdf(name):
    d = pymupdf.open(os.path.join(out, name))
    pg = d[0]
    w_in, h_in, txt = pg.rect.width / 72.0, pg.rect.height / 72.0, pg.get_text()
    d.close()
    return w_in, h_in, txt


w7, h7, t7 = pdf("arf_fig7_mosaic_map.pdf")
from PIL import Image
w7png = Image.open(os.path.join(out, "arf_fig7_mosaic_map.png")).size[0] / 300.0
check(w7 <= 7.16 and w7png <= 7.16 and "(b) Cloud: Fail-Signal re-planning" in t7 and "Panel (b) cloud:" in t7
      and "p = " in t7.split("Panel (b) cloud:")[1][:80],
      f"Figure 7: {w7:.2f} in (PDF), {w7png:.2f} in (PNG) wide, at most 7.16; short panel (b) title; cloud settings "
      f"in the footer")
_, _, t5 = pdf("arf_fig5_onset_robustness.pdf")
_, _, t9 = pdf("arf_fig9_h1a_gap.pdf")
_, _, t3 = pdf("arf_fig3_latency_analysis.pdf")
check("of reachable tiles (h)" in t5 and "symlog" not in t5 and "ROI tiles" not in t5 and "Physical floor" in t5,
      "Figure 5: time from onset to the H1a level of reachable tiles, hours, linear axis, floor drawn")
check("minus the floor (min)" in t9 and "Gap to the floor, per selector" in t9 and "symlog" not in t9
      and "H1a margin" in t9,
      "Figure 9: gap to the floor in minutes, linear axis, gap per selector on the panel, H1a margin drawn")
check("H1a margin" in t3 and "Physical floor" in t3 and "10 min" not in t3, "Figure 3: floor and H1a margin, no 10-min target")

hv = json.load(open(os.path.join(out, "arf_hypothesis_validation.json")))
lab = ns["_scenario_label"](ns["_scenario_key"]())
tr = hv.get("this_run") or {}
exp = {"H1a": hv["H1_revised"]["assessments"][lab]["status"], "H1b": hv["H1b_image_economy"]["assessments"][lab]["status"],
       "H1 (H1a and H1b)": hv["H1_reframed"]["assessments"][lab]["status"], "H2": hv["H2"]["assessments"][lab]["status"],
       "H3": hv["H3"]["assessments"][lab]["status"]}
check(tr.get("scenario") == lab and tr.get("statuses") == exp and tr.get("ledger_scenarios") == [lab],
      f"JSON this_run: {lab}: {tr.get('statuses')}")
check(all(hv[k].get("superseded") is True and hv[k].get("superseded_note") for k in ("H1", "H1_time_aware", "H2_original_ratio"))
      and not any(hv[k].get("superseded") for k in ("H1_revised", "H1b_image_economy", "H1_reframed", "H2", "H3")),
      "JSON: H1, H1_time_aware and H2_original_ratio marked superseded, the protocol hypotheses not")
rep = open(os.path.join(out, "arf_final_report.txt")).read()
check("THIS RUN" in rep and f"Scenario: {lab}" in rep and rep.count("SUPERSEDED ") >= 3
      and "ORIGINAL H1 (kept for the record):" not in rep and "under H1 in" not in rep and "pools the" not in rep,
      "report: THIS RUN block first; the three superseded hypotheses labelled; one scenario, so no pooled label")
so = ns["_stdout"]
check(f"This run ({lab}):" in so and "SUPERSEDED (kept for the record): H2 (Reliability via Fail-Signal)" in so
      and "SUPERSEDED (kept for the record): H1 (Latency Reduction)" in so
      and "superseded H1 target" in so,
      "console: this run's statuses first; superseded hypotheses labelled (Cells 7C and 8)")
log = open(os.path.join(out, "arf_pipeline_log.txt")).read()
check("SUPERSEDED (kept for the record): H1 (Latency Reduction):" in log
      and "SUPERSEDED (kept for the record): H2 (Reliability via Fail-Signal):" in log
      and not re.search(r"^H2 \(Reliability via Fail-Signal\):", log, re.M),
      "pipeline log (Cell 9): superseded hypotheses labelled")

# --- two scenarios in the ledger: overall statuses say they pool --------------------------------
key = ns["_scenario_key"]()
ns["scenario_results"][(key[0] + "_copy", key[1])] = dict(ns["scenario_results"][key])
sink = io.StringIO()
old = os.getcwd()
os.chdir(ns["_work_dir"])
try:
    with contextlib.redirect_stdout(sink):
        exec(compile(e2e15.strip_org(C["8"]), "cell8", "exec"), ns)
    rep2 = open(os.path.join(out, "arf_final_report.txt")).read()
    hv2 = json.load(open(os.path.join(out, "arf_hypothesis_validation.json")))
    ok = ("pools the 2 scenarios in this kernel" in rep2 and "NOTE: the ledger holds 2 scenarios" in rep2
          and hv2["this_run"]["scenario"] == lab and len(hv2["this_run"]["ledger_scenarios"]) == 2
          and "pools them all" in sink.getvalue())
except Exception as e:
    ok, rep2 = False, repr(e)
finally:
    os.chdir(old)
    del ns["scenario_results"][(key[0] + "_copy", key[1])]
check(ok, "two scenarios in the ledger: overall statuses labelled as pooled (report, console); this_run names the "
          "run's own scenario")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
sys.exit(1 if n_fail else 0)
