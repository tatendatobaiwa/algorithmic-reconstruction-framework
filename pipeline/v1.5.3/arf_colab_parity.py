# =============================================================================
# ARF COLAB PARITY CHECK (v1.5.3)
# =============================================================================
# Does the release compute the same numbers in Colab as in the environment it was
# tested in? Re-runs Cells 6-8 (6, 6B, 7, 7B, 7C, 8) of the release on the SAVED
# inputs of earlier runs (their arf_refined_events.csv and tile list on Drive)
# and compares the results with the reference values shipped with the release,
# computed from the same inputs outside Colab. Cells 0-5 (orbital elements,
# propagation) are not re-run: their saved output is the input here, checked by
# SHA-256 so a changed file is reported, not used.
#
# Reference values (v1.5.3):
#   arf_parity_expected_v153.json        the five complete runs of 15 October 2026,
#       in their own Drive folders (MyDrive/ARF_Outputs_ckgr_okavango_...): nothing
#       to restore. Cells 0-5 of these runs were reproduced outside Colab with real
#       Orekit row for row, and Cells 6-8 field for field.
#   arf_parity_expected_v153_12oct.json  the 15 runs of 12 October 2026 (three
#       regions). Their folders were deleted from Drive; unzip
#       ARF_12Oct2026_saved_inputs.zip into a NEW folder, e.g.
#       MyDrive/ARF_parity_12Oct, and pass --drive pointing at it. Never unzip it
#       into MyDrive itself: two of its folder names are now used by the runs of
#       15 October, whose inputs would be overwritten.
#
# Use (Colab, after mounting Drive; no Orekit or Space-Track needed):
#   !python /content/drive/MyDrive/ARF_code/arf_colab_parity.py \
#       --pipeline /content/drive/MyDrive/ARF_code/algorithmic_reconstruction_framework.py \
#       --expected /content/drive/MyDrive/ARF_code/arf_parity_expected_v153.json
#   and, for the 12 October set:
#   !python /content/drive/MyDrive/ARF_code/arf_colab_parity.py \
#       --pipeline /content/drive/MyDrive/ARF_code/algorithmic_reconstruction_framework.py \
#       --expected /content/drive/MyDrive/ARF_code/arf_parity_expected_v153_12oct.json \
#       --drive /content/drive/MyDrive/ARF_parity_12Oct
# Writes arf_parity_report.json next to --expected and prints one line per run:
#   MATCH        every compared value equal (floats to a relative 1e-9)
#   DIFFERS      some values differ -- listed; ILP fields depend on the solver
#                version (HiGHS via SciPy) and are listed separately
#   SKIPPED      the run's saved input is not the one the reference used
# v1.5.2: the fleet tag of an open-data fleet carries no swath case (as in
# Cell 0), and the scenario key is the fleet tag (Cell 0, v1.5.2).
# v1.5.3: the reference values of the 15 October runs added; the procedure is
# unchanged.
# =============================================================================
import os
import io
import re
import sys
import json
import math
import time
import typing
import hashlib
import logging
import argparse
import tempfile
import contextlib

import numpy as np
import pandas as pd

Re_KM, MU = 6378.137, 398600.4418

# Scalar results compared, as paths into scenario_results[<scenario>]
METRICS = (
    "candidate_events", "coverage_clear_pct", "coverage_cloud_pct", "ta_coverage_clear_pct", "ta_coverage_cloud_pct",
    "latency_clear_min", "latency_cloud_min", "ta_latency_clear_min", "ta_latency_cloud_min", "mean_quality_clear",
    "mean_quality_cloud", "ta_mean_quality_clear", "ta_mean_quality_cloud", "h3_redundancy_clear", "h3_passes_clear",
    "reachable_usable_clear_pct", "latency_decay_tau_min", "cloud_replans", "cloud_compensatory_passes",
    "selected_events_clear", "selected_events_cloud",
    "onset_robustness/h1_revised/series/earliest_completion_clear/n_onsets_evaluated",
    "onset_robustness/h1_revised/series/earliest_completion_clear/pct_within_margin",
    "onset_robustness/h1_revised/series/time_aware_clear/pct_within_margin",
    "onset_robustness/h1b_image_economy/series/clear/n_onsets_evaluated",
    "onset_robustness/h1b_image_economy/series/clear/pct_success",
    "onset_robustness/h1b_image_economy/series/clear/median_images_ecc",
    "onset_robustness/h1b_image_economy/series/clear/median_images_fcfs",
    "cloud_robustness/h2/time_aware/difference_vs_openloop_expected/mean",
    "cloud_robustness/h2/time_aware/difference_vs_openloop_expected/lcb",
    "cloud_robustness/h2/time_aware/images_mean",
    "cloud_robustness/h2/time_aware/coverage_pct_mean/adaptive",
    "cloud_robustness/h2/time_aware/coverage_pct_mean/openloop_expected",
    "cloud_robustness/h2/coverage_first/difference_vs_openloop_expected/mean",
    "ilp_gap/time_aware_tw_ratio_min", "ilp_gap/coverage_first_ratio_min", "ilp_gap/time_aware_ratio_min",
)
STATUSES = ("H1_reframed", "H1_revised", "H1b_image_economy", "H2", "H3")
SOLVER_DEPENDENT = ("ilp_gap/",)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tile_set_sha(tiles):
    return hashlib.sha256("\n".join(sorted(tiles)).encode()).hexdigest()


def _get(d, path):
    for k in path.split("/"):
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def fleet_tag(fleet, case):
    """Cell 0's folder tag of a real fleet (48 h live window): the swath case only with PlanetScope."""
    return fleet.replace("+", "_") + (f"-{case}" if "planet" in fleet else "")


def replay(pipeline_path, refined_csv, tiles, orbit_min, roi, fleet, case, n_sats, quiet=True):
    """Run Cells 6-8 of the pipeline on saved inputs; returns the namespace.
    The same procedure produced the reference values."""
    import matplotlib
    matplotlib.use("Agg")
    logging.disable(logging.INFO)
    src = open(pipeline_path, encoding="utf-8").read()
    m = re.search(r'^__version__ = "([^"]+)"', src, flags=re.M)
    release = m.group(1) if m else "unknown"

    def section(start, end):
        i = src.index(start)
        return src[i:src.index(end, i)]
    body = section("# CELL 6:", "# CELL 7D:") + section("# CELL 8:", "# CELL 9:")
    alt_km = (MU * (orbit_min * 60.0 / (2 * math.pi)) ** 2) ** (1.0 / 3.0) - Re_KM
    t_orbit = 2 * np.pi * np.sqrt((Re_KM + alt_km) ** 3 / MU) / 60.0
    refined = pd.read_csv(refined_csv)
    out_dir, work_dir = tempfile.mkdtemp(prefix="ARF_PARITY_OUT_"), tempfile.mkdtemp(prefix="ARF_PARITY_CWD_")
    cwd = os.getcwd()
    os.chdir(work_dir)
    tag = fleet_tag(fleet, case)
    ns = dict(np=np, pd=pd, logging=logging, os=os, json=json, math=math, display=lambda *a, **k: None,
              Re_km=Re_KM, mu_earth=MU, OUT_DIR=out_dir, CONSTELLATION_ALTITUDE_KM=alt_km, T_MAX_WAIT_MINUTES=t_orbit,
              SWATH_HALF_ANGLE_DEG=9.22, WINDOW_HOURS=48, H3_RESOLUTION=7, target_tiles=set(tiles),
              _tile_ids=np.array(sorted(tiles)), refined_events_df=refined, ACTIVE_N=1000, ACTIVE_ROI=roi,
              CLOUD_PROBABILITY=0.35, MAX_SELECTION_STEPS=30, RANDOM_SEED=42, LATENCY_DECAY_RATE=3.0,
              DECAY_SWEEP=(0.5, 1.0, 2.0, 3.0, 5.0), ONSET_ROBUSTNESS_ENABLED=True, ONSET_SHIFT_HOURS=24,
              ONSET_SHIFT_STEP_MIN=20, RANDOM_BASELINE_DRAWS=10, ILP_ENABLED=True, ILP_INSTANCES=5,
              ILP_HORIZON_HOURS=6, ILP_BUDGET_K=10, ILP_TIME_LIMIT_S=120.0, ILP_MAX_PAIRS=100000,
              COVERAGE_TARGETS_PCT=(95.0, 99.0), H1_REACHABLE_COVERAGE_PCT=95.0, H1_FLOOR_MARGIN_MIN_CONFIG=None,
              H1_MIN_ONSET_PCT=90.0, CLOUD_KNOWLEDGE="realistic", FLEET=fleet, SCENARIO_N=tag,
              PLANET_SWATH_CASE=case, FLEET_TAG=tag,
              TIME_AWARE_OBJECTIVE="facility", CLOUD_MODEL="correlated", CLOUD_CORR_LENGTH_KM=50.0,
              CLOUD_CORR_TIME_H=3.0, CLOUD_PROBABILITY_SOURCE="constant", CLOUD_CLIMATOLOGY_MONTH=None,
              CLOUD_PLANNING="expected", CLOUD_SEEDS=20, CLOUD_SWEEP_LENGTHS_KM=(0.0, 10.0, 25.0, 50.0, 100.0, 200.0),
              CLOUD_SWEEP_TIMES_H=(1.0, 3.0, 6.0, 24.0), CLOUD_SWEEP_SEEDS=10, LATENCY_DECAY_SCALE="fleet",
              constellation=[{"id": f"S{i}"} for i in range(int(n_sats))], scenario_results={},
              __name__="__main__", __version__=release,
              _ARF_CELLS_RUN={c: release for c in ("0", "1", "2", "3")},
              **{n: getattr(typing, n) for n in ("List", "Dict", "Tuple", "Any", "Optional", "Set")})
    import h3
    ns["h3"] = h3
    sink = io.StringIO()
    try:
        with contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext():
            exec(compile(body, f"cells6to8:{os.path.basename(pipeline_path)}", "exec"), ns)
    finally:
        os.chdir(cwd)
    ns["_stdout"], ns["_out_dir"] = sink.getvalue(), out_dir
    return ns


def collect(ns):
    """Compared values of one replay: METRICS and the hypothesis statuses."""
    d = ns["scenario_results"][ns["_scenario_key"]()]
    out = {p: _get(d, p) for p in METRICS}
    hv = ns.get("hypothesis_validation") or {}
    for h in STATUSES:
        a = (hv.get(h) or {}).get("assessments") or {}
        out[f"status/{h}"] = list(a.values())[0].get("status") if a else None
    return json.loads(json.dumps(out, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def compare(expected, got, rtol=1e-9):
    diffs = []
    for k, e in expected.items():
        g = got.get(k)
        same = (e == g)
        if not same and isinstance(e, (int, float)) and isinstance(g, (int, float)):
            same = math.isclose(float(e), float(g), rel_tol=rtol, abs_tol=1e-12)
        if not same:
            diffs.append({"field": k, "expected": e, "got": g,
                          "solver_dependent": any(k.startswith(p) for p in SOLVER_DEPENDENT)})
    return diffs


def main(argv=None):
    ap = argparse.ArgumentParser(description="ARF Colab parity check (Cells 6-8 on saved inputs).")
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--drive", default="/content/drive/MyDrive")
    ap.add_argument("--runs", nargs="*", default=None, help="run folder names (default: all in --expected)")
    a = ap.parse_args(argv)
    exp = json.load(open(a.expected))
    src_sha = sha256_file(a.pipeline)
    report = {"pipeline": a.pipeline, "pipeline_sha256": src_sha, "reference_pipeline_sha256": exp.get("pipeline_sha256"),
              "python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__, "runs": {}}
    try:
        import scipy
        report["scipy"] = scipy.__version__
    except Exception:
        pass
    if src_sha != exp.get("pipeline_sha256"):
        print(f"NOTE: the pipeline file differs from the one the reference used "
              f"({src_sha[:16]} vs {str(exp.get('pipeline_sha256'))[:16]}); differences may come from the code.")
    names = a.runs or list(exp["runs"])
    for name in names:
        r = exp["runs"][name]
        folder = os.path.join(a.drive, name)
        rec = {"folder": folder}
        refined = os.path.join(folder, "arf_refined_events.csv")
        tiles_csv = os.path.join(folder, "arf_mosaic_tiles.csv")
        if not (os.path.exists(refined) and os.path.exists(tiles_csv)):
            rec["status"] = "SKIPPED (saved inputs not found)"
        elif sha256_file(refined) != r["refined_events_sha256"]:
            rec["status"] = "SKIPPED (arf_refined_events.csv differs from the reference input: re-run since?)"
        else:
            tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
            if tile_set_sha(tiles) != r["tile_set_sha256"]:
                rec["status"] = "SKIPPED (tile list differs from the reference input)"
            else:
                t0 = time.time()
                ns = replay(a.pipeline, refined, tiles, r["orbital_period_min"], r["roi"], r["fleet"],
                            r["planet_swath_case"], r["n_satellites"])
                got = collect(ns)
                diffs = compare(r["values"], got)
                core = [x for x in diffs if not x["solver_dependent"]]
                rec.update(status="MATCH" if not diffs else ("DIFFERS" if core else "DIFFERS (solver fields only)"),
                           seconds=round(time.time() - t0, 1), differences=diffs)
        report["runs"][name] = rec
        print(f"{name}: {rec['status']}" + "".join(f"\n    {x['field']}: expected {x['expected']}, got {x['got']}"
                                                   + (" [solver]" if x["solver_dependent"] else "")
                                                   for x in rec.get("differences", [])[:12]), flush=True)
    out = os.path.join(os.path.dirname(os.path.abspath(a.expected)), "arf_parity_report.json")
    json.dump(report, open(out, "w"), indent=2, default=str)
    n = sum(v["status"] == "MATCH" for v in report["runs"].values())
    print(f"{n}/{len(report['runs'])} runs MATCH. Report: {out}")
    return report


if __name__ == "__main__":
    main()
