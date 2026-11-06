"""Release 1.5.3 changes no computed value: Cells 6-8 of 1.5.2 and of 1.5.3 are
replayed on the same saved inputs and every output file is compared.

  CSV            byte for byte
  JSON           value for value (exact); the only differences allowed are the
                 keys 1.5.3 adds: hypothesis_validation["this_run"] and
                 "superseded" / "superseded_note" in H1, H1_time_aware and
                 H2_original_ratio
  report (.txt)  every 1.5.2 line that is not in the 1.5.3 report is listed;
                 they must all be label lines (checked by a pattern)
  PNG            byte for byte for the figures 1.5.3 does not touch; the
                 changed figures are listed

Inputs: the 15 runs of 12 October 2026 (exact Cell 5 output of 1.5.2, as in the
1.5.2 audit; ckgr_open from Colab) through replay150.run, and the five complete
Colab runs of 15 October 2026 (their saved arf_refined_events.csv and tile list)
through arf_colab_parity.replay.

Usage: compare_v152_v153.py P152 P153 PARITY_PY OUT_JSONL (oct12 RUN... | oct15 RUN_DIR...)
"""
import sys, os, re, json, math, hashlib, importlib.util, time, traceback
sys.path.insert(0, "/home/claude/validate")
import pandas as pd

if __name__ == "__main__" and sys.argv[1] == "--recheck":
    P2 = P3 = PAR = OUT = mode = None
    items = []
else:
    P2, P3, PAR, OUT, mode, *items = sys.argv[1:]
    spec = importlib.util.spec_from_file_location("arf_colab_parity", PAR)
    par = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(par)

FIGS_CHANGED = {"arf_fig1_coverage_progression.png", "arf_fig2_quality_progression.png",   # whole-number k ticks
                "arf_fig3_latency_analysis.png", "arf_fig4_coverage_budget.png", "arf_fig5_onset_robustness.png",
                "arf_fig7_mosaic_map.png", "arf_fig8_cloud_robustness.png", "arf_fig9_h1a_gap.png",
                "arf_manuscript_figures.png"}
# Wall-clock measurements differ from run to run whatever the release (e.g. arf_ilp_gap.csv, ilp_seconds):
# compared everywhere except these columns, and listed.
TIMING = re.compile(r"seconds|_time_s$|^time_s$|runtime|elapsed")
ADDED_KEYS = {("this_run",), ("H1", "superseded"), ("H1", "superseded_note"), ("H1_time_aware", "superseded"),
              ("H1_time_aware", "superseded_note"), ("H2_original_ratio", "superseded"),
              ("H2_original_ratio", "superseded_note")}
LABEL_LINE = re.compile(r"ORIGINAL|kept for the record|under H1 in|^H1 \(Latency Reduction\)|"
                        r"^H\d[ab]?[ :(]|^H1 \(v1\.4\.0\)|^H2 \(v1\.5\.0\)|^H3 \(Redundancy")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def flat(d, pre=()):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, pre + (str(k),)))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(flat(v, pre + (f"[{i}]",)))
    else:
        out[pre] = d
    return out


def same(a, b):
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def compare_dirs(d2, d3):
    f2, f3 = set(os.listdir(d2)), set(os.listdir(d3))
    res = {"only_152": sorted(f2 - f3), "only_153": sorted(f3 - f2), "csv": {}, "json": {}, "png_same": [],
           "png_changed": [], "report_lines_removed": [], "report_lines_unexplained": []}
    for f in sorted(f2 & f3):
        p2, p3 = os.path.join(d2, f), os.path.join(d3, f)
        if f.endswith(".csv"):
            if sha(p2) == sha(p3):
                res["csv"][f] = True
            else:
                a, b = pd.read_csv(p2), pd.read_csv(p3)
                tcols = [c for c in a.columns if TIMING.search(str(c))]
                res["csv"][f] = (list(a.columns) == list(b.columns) and a.drop(columns=tcols).equals(b.drop(columns=tcols)))
                res.setdefault("csv_timing_columns_differ", {})[f] = [c for c in tcols if not a[c].equals(b[c])]
        elif f.endswith(".json"):
            a, b = flat(json.load(open(p2))), flat(json.load(open(p3)))
            diff = [k for k in set(a) | set(b)
                    if not (k in a and k in b and same(a[k], b[k]))]
            bad = [k for k in diff if not (f == "arf_hypothesis_validation.json" and
                                           any(k[:len(x)] == x for x in ADDED_KEYS) and k not in a)]
            res["json"][f] = {"values": len(a), "differences": len(diff), "unexpected": ["/".join(k) for k in bad][:20]}
        elif f.endswith(".png"):
            (res["png_same"] if sha(p2) == sha(p3) else res["png_changed"]).append(f)
        elif f == "arf_final_report.txt":
            l2 = open(p2).read().splitlines()
            l3 = set(open(p3).read().splitlines())
            gone = [l for l in l2 if l not in l3 and not l.startswith("Date (UTC):") and not l.startswith("Framework Version")]
            res["report_lines_removed"] = gone
            res["report_lines_unexplained"] = [l for l in gone if not LABEL_LINE.search(l.strip())]
    res["ok"] = (all(res["csv"].values()) and all(not v["unexpected"] for v in res["json"].values())
                 and set(res["png_changed"]) <= FIGS_CHANGED and not res["report_lines_unexplained"]
                 and not res["only_152"] and not res["only_153"])
    return res


def oct12(name):
    import replay150 as R
    R.RD = "/home/claude/realdata/" if name == "ckgr_open" else "/home/claude/runs/exact152/"
    ns2 = R.run(P2, name)
    ns3 = R.run(P3, name)
    return ns2, ns3


def oct15(run_dir):
    man = json.load(open(os.path.join(run_dir, "arf_reproducibility_manifest.json")))
    inp = man["inputs"]
    fleet, case, roi = inp["fleet"], inp["config"]["planet_swath_case"], inp["active_roi"]
    key = f"{roi}__N{par.fleet_tag(fleet, case)}"
    refined = os.path.join(run_dir, "arf_refined_events.csv")
    tiles = list(pd.read_csv(os.path.join(run_dir, "arf_mosaic_tiles.csv"), usecols=["tile"])["tile"])
    orbit = man["scenario_results"][key]["orbital_period_min"]
    ns2 = par.replay(P2, refined, tiles, orbit, roi, fleet, case, inp["n_satellites"])
    ns3 = par.replay(P3, refined, tiles, orbit, roi, fleet, case, inp["n_satellites"])
    return ns2, ns3


if __name__ == "__main__" and sys.argv[1] == "--recheck":
    # Re-apply the file comparison to the output folders recorded in a JSONL (no replay).
    src_jsonl, dst_jsonl = sys.argv[2], sys.argv[3]
    with open(dst_jsonl, "w") as fo:
        for line in open(src_jsonl):
            rec = json.loads(line)
            if "dirs" in rec:
                sr_ok = not rec.get("scenario_results_differences") and rec.get("releases") == ["1.5.2", "1.5.3"]
                for k in ("csv", "json", "png_same", "png_changed", "report_lines_removed",
                          "report_lines_unexplained", "only_152", "only_153", "csv_timing_columns_differ"):
                    rec.pop(k, None)
                rec.update(compare_dirs(*rec["dirs"]))
                rec["ok"] = rec["ok"] and sr_ok
            fo.write(json.dumps(rec, default=str) + "\n")
            print(rec["input"], "IDENTICAL" if rec.get("ok") else "CHECK", rec.get("csv_timing_columns_differ", {}))
    sys.exit(0)

for it in items:
    t0 = time.time()
    rec = {"input": it, "mode": mode}
    try:
        ns2, ns3 = (oct12 if mode == "oct12" else oct15)(it)
        rec["releases"] = [ns2.get("__version__"), ns3.get("__version__")]
        rec.update(compare_dirs(ns2["_out_dir"], ns3["_out_dir"]))
        s2 = ns2["scenario_results"][ns2["_scenario_key"]()]
        s3 = ns3["scenario_results"][ns3["_scenario_key"]()]
        a, b = flat(json.loads(json.dumps(s2, default=str))), flat(json.loads(json.dumps(s3, default=str)))
        rec["scenario_results_fields"] = len(a)
        rec["scenario_results_differences"] = ["/".join(k) for k in set(a) | set(b)
                                               if not (k in a and k in b and same(a[k], b[k]))][:20]
        rec["ok"] = rec["ok"] and not rec["scenario_results_differences"] and rec["releases"] == ["1.5.2", "1.5.3"]
        rec["dirs"] = [ns2["_out_dir"], ns3["_out_dir"]]
    except Exception:
        rec["ok"], rec["error"] = False, traceback.format_exc()[-2500:]
    rec["seconds"] = round(time.time() - t0, 1)
    with open(OUT, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")
    print(f"{it}: {'IDENTICAL' if rec.get('ok') else 'CHECK'} | csv {sum(rec.get('csv', {}).values())}/"
          f"{len(rec.get('csv', {}))} same | json {len(rec.get('json', {}))} files, unexpected "
          f"{sum(len(v['unexpected']) for v in rec.get('json', {}).values())} | png same {len(rec.get('png_same', []))}, "
          f"changed {rec.get('png_changed')} | report lines changed {len(rec.get('report_lines_removed', []))}, "
          f"unexplained {len(rec.get('report_lines_unexplained', []))} | {rec['seconds']} s", flush=True)
