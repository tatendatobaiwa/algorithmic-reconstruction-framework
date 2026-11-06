# =============================================================================
# ARF DENSE-CONSTELLATION SUMMARY (v1.5.4, protocol amendment A3): descriptive
# =============================================================================
# Reads the run folders of the batch exploratory_v154_dense (index in ARF_batch)
# and tabulates, per window and fleet size N, the same window outcomes as the
# confirmatory analysis (window_outcomes of arf_confirmatory_analysis.py: H1a and
# H1b shares over evaluable onsets, H2 coverage difference vs the open-loop plan),
# with coverage, latency, images, redundancy and run time from the run's manifest
# and the batch index. Two summaries, both DESCRIPTIVE (A3: no hypothesis
# decision): the scaling over N in window 1 (log-log slope of run time and of
# candidate events against N), and the spread across the eight windows
# for N <= 1000 (mean, SD, min, max). Only runs whose manifest records the arm's
# release and settings are used. Writes CSV files and a LaTeX table (fully
# enclosed grid) to <BATCH_DIR>/dense_summary/.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_dense_summary.py
# =============================================================================
import os
import sys
import json
import math
import runpy
import argparse

import numpy as np
import pandas as pd

# latency_clear_min of a run = completion time of the COVERAGE-FIRST selector (no preference for early passes):
# not a scaling measure; kept in the CSV as coverage_first_latency_clear_min. Latency across N is read per onset
# at the H1a level (95% of reachable tiles) from arf_onset_robustness.csv: the physical floor and ECC (v1.5.4 fix).
SR_FIELDS = ("coverage_clear_pct", "coverage_cloud_pct", "latency_clear_min", "ta_latency_clear_min",
             "latency_cloud_min", "selected_events_clear", "h3_redundancy_clear", "reachable_tiles_pct",
             "candidate_events", "orbital_period_min", "real_elbow_k")


def run_rows(index_csv, analysis, a3, min_onsets):
    idx = pd.read_csv(index_csv, dtype=str).fillna("")
    rows, skipped = [], []
    for _, r in idx.iterrows():
        d = r["out_dir"]
        if r["status"] != "done" or not os.path.exists(os.path.join(d, "arf_reproducibility_manifest.json")):
            skipped.append((r["run_id"], f"status {r['status']}"))
            continue
        m = json.load(open(os.path.join(d, "arf_reproducibility_manifest.json")))
        cfg = m.get("config_overrides") or {}
        if (m.get("framework_version") != a3["release_of_arm"] or cfg.get("fleet") != "starlink"
                or r["pipeline_sha256"] != a3["pipeline_sha256_of_arm"]):
            skipped.append((r["run_id"], "not a run of the arm's release"))
            continue
        w = analysis["window_outcomes"](d, min_onsets)
        sr = list((m.get("scenario_results") or {}).values())
        sr = sr[0] if sr else {}
        row = {"window": str(cfg.get("analysis_epoch"))[:10], "N": int(cfg.get("active_n")),
               "n_satellites": (m.get("inputs") or {}).get("n_satellites"), "run_minutes": float(r["minutes"]),
               **{k: w.get(k) for k in ("h1a_pct", "h1a_onsets", "h1b_pct", "h1b_onsets", "h1b_median_images_ecc",
                                       "h1b_median_images_fcfs", "h2_diff", "h2_evaluable")},
               **{k: sr.get(k) for k in SR_FIELDS}, "run_dir": d}
        row["coverage_first_latency_clear_min"] = row.pop("latency_clear_min")
        on = pd.read_csv(os.path.join(d, "arf_onset_robustness.csv"))
        for col, name in (("best_possible_clear_t95r_min", "floor_t95r"), ("earliest_completion_clear_t95r_min", "ecc_t95r")):
            v = pd.to_numeric(on[col], errors="coerce").dropna() if col in on.columns else pd.Series(dtype=float)
            row[f"{name}_median_min"] = float(v.median()) if v.size else None
            row[f"{name}_p90_min"] = float(v.quantile(0.9)) if v.size else None
        row["onsets"] = int(len(on))
        rows.append(row)
    return pd.DataFrame(rows), skipped


def loglog_slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = (x > 0) & (y > 0) & np.isfinite(y)
    if ok.sum() < 2:
        return None
    return float(np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)[0])


def fmt(v, nd=1):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    return f"{v:,.{nd}f}" if isinstance(v, float) else f"{v:,}"


def latex_table(df):
    cols = (("N", "$N$", 0), ("coverage_cloud_pct", "Cov. cloud (\\%)", 1),
            ("floor_t95r_median_min", "Floor (min)", 0), ("ecc_t95r_median_min", "ECC (min)", 0),
            ("h1a_pct", "H1a (\\%)", 1), ("h1b_pct", "H1b (\\%)", 1), ("h2_diff", "H2 (pts)", 2),
            ("selected_events_clear", "Images", 0), ("run_minutes", "Run (min)", 0))
    lines = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
             "\\caption{Dense-constellation projection (exploratory, amendment A3): window of 3 August 2026, "
             "Starlink geometry as a hypothetical Earth-observation constellation (imaging parameters assumed). "
             "Floor and ECC: median over onsets of the time to 95\\% of reachable tiles, clear sky; "
             "the physical floor uses every candidate pass.}",
             "\\label{tab:dense_window1}", "\\begin{tabular}{|" + "|".join("r" for _ in cols) + "|}", "\\hline",
             " & ".join(h for _, h, _ in cols) + " \\\\", "\\hline"]
    for _, r in df.sort_values("N").iterrows():
        def cell(v, nd):
            if v is None or not np.isfinite(float(v)):
                return "--"
            return f"{int(round(float(v))):,}" if nd == 0 else f"{float(v):,.{nd}f}"
        lines.append(" & ".join(cell(r[c], nd) for c, _, nd in cols) + " \\\\")
        lines.append("\\hline")
    lines += ["\\end{tabular}", "\\end{table*}"]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Descriptive summary of the dense-constellation arm (A3).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    a, _ = ap.parse_known_args(argv)
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    a3 = [x for x in proto["amendments"] if x["id"] == "A3"][0]
    analysis = runpy.run_path(os.path.join(a.code, "arf_confirmatory_analysis.py"), run_name="arf_confirmatory_analysis")
    df, skipped = run_rows(os.path.join(a.batch_dir, "exploratory_v154_dense_index.csv"), analysis, a3,
                           int(proto["min_evaluable_onsets"]))
    out = os.path.join(a.batch_dir, "dense_summary")
    os.makedirs(out, exist_ok=True)
    df.sort_values(["window", "N"]).to_csv(os.path.join(out, "dense_runs.csv"), index=False)
    first = str(proto["windows"][0])[:10]
    w1 = df[df.window == first].sort_values("N")
    scaling = {"window": first, "N": w1.N.tolist(),
               "slope_run_minutes_vs_N": loglog_slope(w1.N, w1.run_minutes),
               "slope_candidate_events_vs_N": loglog_slope(w1.N, w1.candidate_events),
               "slope_floor_t95r_median_vs_N": loglog_slope(w1.N, w1.floor_t95r_median_min),
               "slope_ecc_t95r_median_vs_N": loglog_slope(w1.N, w1.ecc_t95r_median_min),
               "note": "descriptive log-log least-squares slope over the fleet sizes run in window 1"}
    spread = []
    for n, g in df[df.N <= 1000].groupby("N"):
        for col in ("h1a_pct", "h1b_pct", "h2_diff", "coverage_cloud_pct", "floor_t95r_median_min", "ecc_t95r_median_min",
                    "selected_events_clear"):
            v = pd.to_numeric(g[col], errors="coerce").dropna()
            spread.append({"N": int(n), "quantity": col, "windows": int(v.size),
                           "mean": float(v.mean()) if v.size else None,
                           "sd": float(v.std(ddof=1)) if v.size > 1 else None,
                           "min": float(v.min()) if v.size else None, "max": float(v.max()) if v.size else None})
    pd.DataFrame(spread).to_csv(os.path.join(out, "dense_spread_windows.csv"), index=False)
    json.dump(scaling, open(os.path.join(out, "dense_scaling_window1.json"), "w"), indent=1)
    open(os.path.join(out, "dense_window1_table.tex"), "w").write(latex_table(w1))
    n_exp = a3["arm"]["n_runs"]
    print(f"Dense-constellation arm (A3, release {a3['release_of_arm']}): {len(df)}/{n_exp} runs summarised"
          + (f"; skipped: {skipped}" if skipped else ""))
    if len(w1):
        print(w1[["N", "coverage_cloud_pct", "floor_t95r_median_min", "ecc_t95r_median_min", "h1a_pct", "h1b_pct", "h2_diff",
                  "selected_events_clear", "run_minutes"]].to_string(index=False))
        print(f"log-log slope over N (window {first}): run time {fmt(scaling['slope_run_minutes_vs_N'], 2)}, "
              f"candidate events {fmt(scaling['slope_candidate_events_vs_N'], 2)}, floor (median, 95% of reachable) "
              f"{fmt(scaling['slope_floor_t95r_median_vs_N'], 2)}, ECC (median, 95% of reachable) {fmt(scaling['slope_ecc_t95r_median_vs_N'], 2)}")
    print(f"Files: {out}/dense_runs.csv, dense_spread_windows.csv, dense_scaling_window1.json, dense_window1_table.tex")
    return df


if __name__ == "__main__":
    main()
