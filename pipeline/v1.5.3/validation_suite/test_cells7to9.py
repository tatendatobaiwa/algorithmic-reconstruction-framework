"""Cells 7, 7B, 7C (figures), 8 and 9 of v1.5.2 on one end-to-end run
(e2e15 harness: synthetic Earth-observation fleet over the real Lower Zambezi
tiles; Cells 4-9 from the pipeline file). Usage: test_cells7to9.py PIPELINE
Cell 7  figures drawn at IEEE widths (PDF page sizes); Figure 3 shows the
        physical floor and the H1a margin, not the falsified 10-min target
        (PDF text); its floor and per-k times to the H1a level equal an
        independent recomputation.
Cell 7B arf_coverage_vs_time.csv: the time each curve reaches its final
        coverage equals Cell 6's completion latency; an image with no downlink
        is left out instead of crashing.
Cell 7C Figures 5 and 9 exist at two-column width; no 10-min line.
Cell 8  H3 overall status: SUPPORTED only if every evaluable scenario is; a
        ledger of two scenarios produces the full report.
Cell 9  manifest: JPype1 and Java recorded, cell versions consistent, every
        expected file (with the new 7B/7C outputs) found; scenario key by fleet tag."""
import sys, os, io, re, ast, json, math, contextlib, logging
import numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/validate")
import e2e15
import pymupdf

P = sys.argv[1]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


tiles = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
# Cell 3 runs too (its map helpers draw Cell 7's mosaic map), on the harness's stand-ins
ns = e2e15.run(P, tiles, overrides={"window_hours": 48, "onset_shift_hours": 24, "onset_shift_step_min": 120},
               cells=("3", "4", "5", "6", "6B", "7", "7B", "7C", "8", "9"))
out = ns["_out_dir"]


def pdf(name):
    d = pymupdf.open(os.path.join(out, name))
    pg = d[0]
    w_in, h_in, txt = pg.rect.width / 72.0, pg.rect.height / 72.0, pg.get_text()
    d.close()
    return w_in, h_in, txt


# Cell 7: sizes and Figure 3 content
sizes = {n: pdf(n)[:2] for n in ("arf_fig1_coverage_progression.pdf", "arf_fig2_quality_progression.pdf",
                                 "arf_fig3_latency_analysis.pdf", "arf_fig4_coverage_budget.pdf",
                                 "arf_manuscript_figures.pdf", "arf_fig5_onset_robustness.pdf", "arf_fig9_h1a_gap.pdf",
                                 "arf_fig6_cost_latency.pdf", "arf_fig8_cloud_robustness.pdf")}
single = all(2.9 <= sizes[n][0] <= 3.75 for n in list(sizes)[:4])
double = all(6.6 <= sizes[n][0] <= 7.45 for n in list(sizes)[4:])
check(single and double, "figure widths (in): " + ", ".join(f"{k.split('_')[1]} {v[0]:.2f}" for k, v in sizes.items())
      + " (one column 3.5, two columns 7.16; bbox_inches='tight' trims or pads a little)")
_, _, t3 = pdf("arf_fig3_latency_analysis.pdf")
_, _, t5 = pdf("arf_fig5_onset_robustness.pdf")
_, _, t9 = pdf("arf_fig9_h1a_gap.pdf")
check("H1a margin" in t3 and "Physical floor" in t3 and "10 min" not in t3 and "H1 target" not in t5.replace("H1a", "")
      and "H1a margin" in t9 and "(h)" in t3,
      "Figure 3: physical floor and H1a margin in hours, no 10-min target; Figure 5 without the 10-min line; "
      "Figure 9 shows the H1a margin")
meta = json.load(open(os.path.join(out, "arf_cell7_metadata.json")))
f3 = meta["figure3"]
st = ns["event_data"]
lvl = f3["h1_level_pct_of_reachable"]
reach = int(np.unique(st.tile_idx[st.util > 0]).size)
need = int(np.ceil(lvl / 100.0 * reach - 1e-9))


def t_level(keys):
    """Independent: deliveries sorted by (time, original position), count usable tiles."""
    order = sorted(range(len(keys)), key=lambda i: (keys[i][1] + st.wait_lookup[keys[i]], i))
    seen = set()
    for i in order:
        k = keys[i]
        s = st.rows(st.index[k])
        for t, u in zip(st.tile_idx[s], st.util[s]):
            if u > 0:
                seen.add(int(t))
        if len(seen) >= need:
            return float(k[1] + st.wait_lookup[k])
    return None


floor_ind = t_level(list(st.keys_))
hist = ns["run_submodular_greedy_selection"](st, max(int(k) for k in f3["time_to_level_min_by_k"]), gap_fill=False)[3]
sel = [h["selected"] for h in hist]
per_k_ok = all((v is None and t_level(sel[:int(k)]) is None) or (v is not None and abs(v - t_level(sel[:int(k)])) < 1e-9)
               for k, v in f3["time_to_level_min_by_k"].items())
check(f3["reachable_tiles"] == reach and abs(f3["physical_floor_min"] - floor_ind) < 1e-9 and per_k_ok
      and abs(f3["h1a_margin_min"] - ns["H1_FLOOR_MARGIN_MIN"]) < 1e-9,
      f"Figure 3 data: floor {f3['physical_floor_min']:.0f} min and the time to {lvl:g}% of {reach} reachable tiles at "
      f"each of {len(f3['time_to_level_min_by_k'])} budgets equal an independent recomputation; margin = Cell 7C's "
      f"({f3['h1a_margin_min']:.1f} min)")

# Cell 7B
cvt = pd.read_csv(os.path.join(out, "arf_coverage_vs_time.csv"))
d = ns["scenario_results"][ns["_scenario_key"]()]
ok = True
for sky, lat in (("clear", d["latency_clear_min"]), ("cloud", d["latency_cloud_min"])):
    c = cvt[cvt.sky == sky]
    final = c.coverage_pct.max()
    t_final = c[c.coverage_pct >= final - 1e-9].delivery_time_min.min()
    ok &= abs(t_final - lat) < 1e-9
check(ok and set(cvt.columns) == {"sky", "delivery_time_min", "coverage_pct", "sat", "t_min"},
      "arf_coverage_vs_time.csv: each curve reaches its final coverage at Cell 6's completion latency (clear, cloud)")
k0 = st.keys_[0]
st.wait_lookup[k0] = None
try:
    curve = ns["compute_coverage_vs_time"]([k0, st.keys_[1]], st, ns["refined_events_df"], ns["total_target_tiles"])
    ok = [p["event"] for p in curve] == [st.keys_[1]]
except Exception as e:
    ok = False
check(ok, "Cell 7B: an image with no downlink (wait None) is left out of the curve, no crash")

# Cell 8: H3 overall rule and a two-scenario ledger
C = e2e15.cells_of(open(P, encoding="utf-8").read())
tree = ast.parse(C["8"])
tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ("_scenario_label",
                                                                                    "validate_h3_redundancy_minimisation")]
n8 = {"logger": logging.getLogger("t8"), "H3_COVERAGE_THRESHOLD_PCT": 95.0, "H3_REDUNDANCY_THRESHOLD": 1.5}
exec(compile(tree, "c8", "exec"), n8)
led = {("a", "f1"): {"coverage_clear_pct": 99.0, "h3_coverage_pct": 95.0, "reachable_usable_clear_pct": 99.0,
                     "h3_redundancy_clear": 1.2},
       ("b", "f1"): {"coverage_clear_pct": 99.0, "h3_coverage_pct": 95.0, "reachable_usable_clear_pct": 99.0,
                     "h3_redundancy_clear": 2.0},
       ("c", "f1"): {"coverage_clear_pct": 60.0, "h3_coverage_pct": 95.0, "reachable_usable_clear_pct": 60.0,
                     "h3_redundancy_clear": None}}
r_mixed = n8["validate_h3_redundancy_minimisation"](led)["overall_status"]
r_ok = n8["validate_h3_redundancy_minimisation"]({k: v for k, v in led.items() if k != ("b", "f1")})["overall_status"]
r_none = n8["validate_h3_redundancy_minimisation"]({("c", "f1"): led[("c", "f1")]})["overall_status"]
check(r_mixed.startswith("NOT SUPPORTED") and r_ok.startswith("SUPPORTED") and r_none == "NOT EVALUABLE",
      f"H3 overall: one supported + one not -> {r_mixed}; supported + not evaluable -> {r_ok}; only not evaluable -> {r_none}")
key = ns["_scenario_key"]()
ns["scenario_results"][(key[0] + "_copy", key[1])] = dict(ns["scenario_results"][key])
sink = io.StringIO()
old = os.getcwd()
os.chdir(ns["_work_dir"])
try:
    with contextlib.redirect_stdout(sink):
        exec(compile(e2e15.strip_org(C["8"]), "cell8", "exec"), ns)
    rep = open(os.path.join(out, "arf_final_report.txt")).read()
    ok = rep.count("Policy comparison, ") == 2 and "lower_zambezi_copy" in rep
except Exception as e:
    ok, rep = False, repr(e)
finally:
    os.chdir(old)
check(ok, "Cell 8 with two scenarios in the ledger: full report, both scenarios listed")
del ns["scenario_results"][(key[0] + "_copy", key[1])]

# Cell 9 (from the e2e run)
man = json.load(open(os.path.join(out, "arf_reproducibility_manifest.json")))
env = man["environment"]
val = man["output_validation"]
check("JPype1" in env["packages"] and "java" in env and env.get("validated_versions") is None or "JPype1" in env["packages"],
      f"manifest records JPype1 ({env['packages'].get('JPype1')}) and Java ({env.get('java')}; no JVM in this harness)")
check(man["cell_versions"]["consistent"] and val["n_missing"] == 0 and "arf_coverage_vs_time.csv" in val["files"]
      and "arf_fig9_h1a_gap.pdf" in val["files"],
      f"cell versions consistent; {val['n_found']}/{val['n_expected']} expected files found, including "
      f"arf_coverage_vs_time.csv and arf_fig9_h1a_gap")
check(list(man["scenario_results"]) == [f"{key[0]}__N{key[1]}"] and key[1] == ns["FLEET_TAG"]
      and f"fleet {ns['FLEET_TAG']}" in ns["_stdout"],
      f"scenario key by fleet tag: {list(man['scenario_results'])}; Cell 9 prints the fleet tag")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cells7to9.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
