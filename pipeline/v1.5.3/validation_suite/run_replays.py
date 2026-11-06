"""Run replay150 for a list of runs with a pipeline file and store a compact
summary per run (JSON lines). Usage: run_replays.py PIPE OUT.jsonl run1 run2 ...
Options via env: ARF_DEDUPE=1 (end-cap correction), ARF_FAMILY=open|planet,
ARF_QUICK=1 (2 cloud seeds, no sweep, 4-hourly onsets)."""
import sys, os, json, time, traceback
sys.path.insert(0, "/home/claude/validate")
import numpy as np
import replay150 as R

pipe, out_path, names = sys.argv[1], sys.argv[2], sys.argv[3:]
quick = os.environ.get("ARF_QUICK") == "1"
# v1.5.2 audit: ARF_REFINED_DIR replaces the refined-events folder (exact Cell 5
# outputs of a release, as <run>_refined_events.csv); ARF_DROP_SATS drops the
# rows of the named satellites (comma-separated), e.g. SENTINEL-2A.
if os.environ.get("ARF_REFINED_DIR"):
    R.RD = os.environ["ARF_REFINED_DIR"].rstrip("/") + "/"
_drop = [s.strip() for s in os.environ.get("ARF_DROP_SATS", "").split(",") if s.strip()]
if _drop:
    _read_csv = R.pd.read_csv

    def _read_without(path, *a, **k):
        df = _read_csv(path, *a, **k)
        if str(path).endswith("_refined_events.csv") and "sat" in df.columns:
            df = df[~df["sat"].astype(str).isin(_drop)].reset_index(drop=True)
        return df
    R.pd.read_csv = _read_without
ov = dict(cloud_seeds=2, cloud_sweep_lengths_km=(), cloud_sweep_times_h=(), onset_shift_step_min=240,
          ilp_instances=1, random_baseline_draws=2) if quick else {}
extra_ov = json.loads(os.environ.get("ARF_OVERRIDES", "{}"))
ov.update(extra_ov)


def g(d, *ks):
    for k in ks:
        if d is None:
            return None
        d = d.get(k) if isinstance(d, dict) else None
    return d


for name in names:
    t0 = time.time()
    rec = {"run": name, "pipe": os.path.basename(pipe), "dedupe": os.environ.get("ARF_DEDUPE") == "1",
           "family": os.environ.get("ARF_FAMILY"), "overrides": ov, "refined_dir": R.RD, "dropped_sats": _drop}
    try:
        ns = R.run(pipe, name, overrides=ov, dedupe=rec["dedupe"], family=rec["family"])
        rec["out_dir"] = ns["_out_dir"]
        d = R.result(ns)
        hv = json.load(open(os.path.join(ns["_out_dir"], "arf_hypothesis_validation.json")))
        lab = list(hv["H1_revised"]["assessments"])[0]
        ob = d.get("onset_robustness") or {}
        cr = d.get("cloud_robustness") or {}
        il = d.get("ilp_gap") or {}
        rec.update({
            "n_events": int(len(ns["event_data"])), "n_rows": int(ns["event_data"].n_rows),
            "reach_pct": d.get("reachable_usable_clear_pct"),
            "cov_clear": d.get("coverage_clear_pct"), "ta_cov_clear": d.get("ta_coverage_clear_pct"),
            "cov_cloud": d.get("coverage_cloud_pct"), "ta_cov_cloud": d.get("ta_coverage_cloud_pct"),
            "lat_clear": d.get("latency_clear_min"), "ta_lat_clear": d.get("ta_latency_clear_min"),
            "tau_min": d.get("latency_decay_tau_min"), "h3_red_clear": d.get("h3_redundancy_clear"),
            "h3_passes_clear": d.get("h3_passes_clear"), "sel_clear": d.get("selected_events_clear"),
            "status": {k: hv[k]["assessments"][lab]["status"] for k in ("H1_reframed", "H1_revised", "H1b_image_economy", "H2", "H3")},
            "h1a": g(ob, "h1_revised", "series", "earliest_completion_clear"),
            "h1a_ta": g(ob, "h1_revised", "series", "time_aware_clear"),
            "h1b": g(ob, "h1b_image_economy", "series", "clear"),
            "h2_ta": {k: g(cr, "h2", "time_aware", k) for k in ("status", "coverage_pct_mean", "images_mean",
                                                                  "difference_vs_openloop_expected", "not_evaluable_reasons",
                                                                  "replans_total", "compensatory_total")},
            "h2_cf": {k: g(cr, "h2", "coverage_first", k) for k in ("status", "coverage_pct_mean", "images_mean",
                                                                      "difference_vs_openloop_expected")},
            "sensitivity": cr.get("sensitivity"),
            "ilp": {k: il.get(k) for k in il if "ratio" in k or "solved" in k or "ecc" in k or "instances" in k},
            "floor": {k: d.get(k) for k in d if k.startswith("floor")},
            "seconds": round(time.time() - t0, 1),
        })
        stray = os.listdir(ns["_work_dir"])
        rec["stray_files"] = stray
    except Exception:
        rec["error"] = traceback.format_exc()[-3000:]
    with open(out_path, "a") as f:
        f.write(json.dumps(rec, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    print(name, "ok" if "error" not in rec else "ERROR", rec.get("seconds"), flush=True)
