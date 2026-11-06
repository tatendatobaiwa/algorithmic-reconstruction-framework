# =============================================================================
# ARF SENSITIVITY SUMMARY (protocol amendment A4, POST HOC; read-only)
# =============================================================================
# Reads the runs of batch sensitivity_v153_a4 (arf_sensitivity_runner_a4.py)
# and the confirmatory primary runs of the same windows, and tabulates:
#   S1 budget k = 5, 10, 15, 20 against the protocol's 30: per window the H1a
#      and H1b shares (the protocol's own rules, window_outcomes), the ECC gap
#      above the physical floor at the H1a level (median and maximum over
#      onsets, minutes), ECC's median images, and the share of onsets at which
#      ECC reaches the level at all
#   S2 cloud probability: constant p = 0.35 (protocol) against the EarthEnv
#      monthly climatology (per tile): H2 difference and the four-plan bracket
#   S3 exact optimum on the full 48-h instance: ECC's greedy cover against the
#      exact minimum number of images (and solver status and time); and a
#      determinism check, every Cell 7C result file of the ilp48 run except
#      arf_ilp_gap.csv equal to the confirmatory run's
# Writes CSV files, a JSON summary and LaTeX tables (fully enclosed grid) to
# <BATCH_DIR>/sensitivity_a4/.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_sensitivity_summary_a4.py
# =============================================================================
import os
import json
import runpy
import hashlib
import argparse

import numpy as np
import pandas as pd

DETERMINISM_FILES = ("arf_onset_robustness.csv", "arf_onset_robustness_summary.json", "arf_cloud_robustness.csv",
                     "arf_coverage_vs_time.csv", "arf_selection_history_clear.csv", "arf_selection_history_cloud.csv",
                     "arf_policy_comparison.csv")


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def variant_of(cfg):
    if "max_selection_steps" in cfg:
        return f"k{int(cfg['max_selection_steps']):02d}"
    if cfg.get("cloud_probability_source") == "earthenv":
        return "earthenv"
    if "ilp_horizon_hours" in cfg:
        return "ilp48"
    return "protocol"


def run_row(d, an, min_onsets, level="95"):
    w = an["window_outcomes"](d, min_onsets)
    on = pd.read_csv(os.path.join(d, "arf_onset_robustness.csv"))
    on = on[_num(on.get("reachable_clear_tiles")).fillna(0) > 0]
    fl = _num(on.get(f"best_possible_clear_t{level}r_min"))
    te = _num(on.get(f"earliest_completion_clear_t{level}r_min"))
    gap = (te - fl).dropna()
    ons = json.load(open(os.path.join(d, "arf_onset_robustness_summary.json")))
    cr = json.load(open(os.path.join(d, "arf_cloud_robustness_summary.json")))
    h = cr["h2"]["time_aware"]
    cov = h.get("coverage_pct_mean") or {}
    prob = (cr.get("cloud_model") or {}).get("probability") or {}
    row = {"window": str(w["analysis_epoch"])[:10], "h1a_pct": w["h1a_pct"], "h1b_pct": w["h1b_pct"], "h2_diff": w["h2_diff"],
           "margin_min": ons["h1_revised"].get("floor_margin_min"),
           "ecc_reached_pct": float(100 * te.notna().mean()) if len(on) else None,
           "gap_median_min": float(gap.median()) if gap.size else None,
           "gap_max_min": float(gap.max()) if gap.size else None,
           "floor_median_min": float(fl.median()) if fl.notna().any() else None,
           "ecc_median_min": float(te.median()) if te.notna().any() else None,
           "ecc_images_median": float(_num(on.get(f"earliest_completion_clear_images{level}r")).median()),
           "fcfs_images_median": float(_num(on.get(f"fcfs_unbudgeted_clear_images{level}r")).median()),
           **{f"cov_{k}": cov.get(k) for k in ("adaptive", "openloop_expected", "openloop_clear", "oracle")},
           "cloud_p_source": prob.get("source"), "cloud_p_mean": prob.get("mean", prob.get("value")),
           "run_dir": d}
    g = os.path.join(d, "arf_ilp_gap.csv")
    if os.path.exists(g):
        ig = pd.read_csv(g)
        row.update({"ilp_rows": int(len(ig)), "ilp_horizons_min": ";".join(str(int(x)) for x in ig["horizon_min"]),
                    "ilp_pairs_max": int(_num(ig["n_pairs"]).max()) if len(ig) else None,
                    "ilp_seconds_max": float(_num(ig["ilp_seconds"]).max()) if len(ig) else None,
                    "ecc_greedy_images": ";".join(str(x) for x in ig.get("ecc_greedy_images", [])),
                    "ecc_exact_images": ";".join(str(x) for x in ig.get("ecc_exact_images", [])),
                    "ecc_exact_status": ";".join(str(x) for x in ig.get("ecc_exact_status", [])),
                    "maxcov_ilp_status": ";".join(str(x) for x in ig["ilp_status"]),
                    "time_aware_ratio": ";".join(str(x) for x in ig.get("time_aware_ratio", []))})
    return row


def table(caption, label, header, rows, align):
    L = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         f"\\caption{{{caption}}}", f"\\label{{{label}}}", "\\begin{tabular}{|" + "|".join(align) + "|}", "\\hline",
         " & ".join(header) + " \\\\", "\\hline"]
    for r in rows:
        L += [" & ".join(r) + " \\\\", "\\hline"]
    return "\n".join(L + ["\\end{tabular}", "\\end{table*}"])


def f(v, nd=1):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "--"
    return "--" if not np.isfinite(v) else f"{v:,.{nd}f}"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Post-hoc sensitivity summary (A4).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    a, _ = ap.parse_known_args(argv)
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    an = runpy.run_path(os.path.join(a.code, "arf_confirmatory_analysis.py"), run_name="arf_confirmatory_analysis")
    mo = int(proto["min_evaluable_onsets"])
    rows = []
    for name, var in (("confirmatory_v153_primary_index.csv", None), ("sensitivity_v153_a4_index.csv", "x")):
        p = os.path.join(a.batch_dir, name)
        if not os.path.exists(p):
            raise SystemExit(f"missing {p}")
        idx = pd.read_csv(p, dtype=str).fillna("")
        for _, r in idx[(idx.status == "done") & (idx.release == proto["release"])].iterrows():
            cfg = json.loads(r["config"]) if r.get("config") else {}
            if cfg.get("analysis_epoch") not in proto["windows"]:
                continue
            rows.append({"variant": "protocol" if var is None else variant_of(cfg), **run_row(r["out_dir"], an, mo)})
    df = pd.DataFrame(rows).sort_values(["variant", "window"])
    out = os.path.join(a.batch_dir, "sensitivity_a4")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "a4_runs.csv"), index=False)
    base = df[df.variant == "protocol"].set_index("window")

    # S1 budget
    s1 = []
    for v in ("k05", "k10", "k15", "k20", "protocol"):
        g = df[df.variant == v]
        if not len(g):
            continue
        s1.append({"k": 30 if v == "protocol" else int(v[1:]), "windows": int(len(g)),
                   "h1a_mean": _num(g.h1a_pct).mean(), "h1a_min": _num(g.h1a_pct).min(),
                   "h1b_mean": _num(g.h1b_pct).mean(), "h1b_min": _num(g.h1b_pct).min(),
                   "gap_median_min": _num(g.gap_median_min).median(), "gap_max_min": _num(g.gap_max_min).max(),
                   "ecc_reached_pct_min": _num(g.ecc_reached_pct).min(),
                   "ecc_images_median": _num(g.ecc_images_median).median(), "h2_mean": _num(g.h2_diff).mean()})
    s1 = pd.DataFrame(s1)
    s1.to_csv(os.path.join(out, "a4_budget.csv"), index=False)
    # S2 cloud
    s2 = []
    for w, g in df[df.variant.isin(["protocol", "earthenv"])].groupby("window"):
        e = g[g.variant == "earthenv"]
        b = g[g.variant == "protocol"]
        if len(e) and len(b):
            e, b = e.iloc[0], b.iloc[0]
            s2.append({"window": w, "p_const": b.cloud_p_mean, "p_earthenv_mean": e.cloud_p_mean,
                       "h2_const": b.h2_diff, "h2_earthenv": e.h2_diff,
                       **{f"earthenv_cov_{k}": e[f"cov_{k}"] for k in ("adaptive", "openloop_expected", "openloop_clear", "oracle")}})
    s2 = pd.DataFrame(s2)
    s2.to_csv(os.path.join(out, "a4_cloud.csv"), index=False)
    # S3 exact optimum + determinism
    s3 = []
    for _, r in df[df.variant == "ilp48"].iterrows():
        same = {}
        if r.window in base.index:
            bd = base.loc[r.window, "run_dir"]
            same = {fn: _sha(os.path.join(bd, fn)) == _sha(os.path.join(r.run_dir, fn)) for fn in DETERMINISM_FILES
                    if os.path.exists(os.path.join(bd, fn))}
        s3.append({"window": r.window, "horizon_min": r.get("ilp_horizons_min"), "pairs": r.get("ilp_pairs_max"),
                   "ecc_greedy_images": r.get("ecc_greedy_images"), "ecc_exact_images": r.get("ecc_exact_images"),
                   "ecc_exact_status": r.get("ecc_exact_status"), "maxcov_status": r.get("maxcov_ilp_status"),
                   "time_aware_ratio": r.get("time_aware_ratio"), "seconds_max": r.get("ilp_seconds_max"),
                   "files_compared": len(same), "files_identical": int(sum(same.values())),
                   "files_differing": ";".join(k for k, v in same.items() if not v)})
    s3 = pd.DataFrame(s3)
    s3.to_csv(os.path.join(out, "a4_ilp48.csv"), index=False)
    summary = {"label": "POST HOC sensitivity (protocol amendment A4); not part of the confirmatory decisions",
               "runs": int(len(df)), "variants": df.variant.value_counts().to_dict(),
               "budget": s1.to_dict("records"), "cloud": s2.to_dict("records"), "ilp48": s3.to_dict("records")}
    json.dump(summary, open(os.path.join(out, "a4_summary.json"), "w"), indent=1, default=str)

    T = []
    if len(s1):
        T.append(table("Post hoc budget sensitivity (amendment A4), primary design, eight windows: ECC with an image "
                       "budget $k$ (protocol: 30). Gap: ECC minus the physical floor at the H1a level (95\\% of "
                       "reachable tiles, clear sky), median and maximum over onsets, in minutes; H1a margin one "
                       "orbital period (94 min). H1a, H1b: mean and minimum over windows, \\%. Reach: smallest share "
                       "of onsets (\\%) at which ECC reaches the level at all; images: median.", "tab:a4_budget",
                       ["$k$", "H1a mean", "H1a min", "H1b mean", "Gap median", "Gap max", "Reach", "Images"],
                       [[str(int(r.k)), f(r.h1a_mean), f(r.h1a_min), f(r.h1b_mean), f(r.gap_median_min, 0),
                         f(r.gap_max_min, 0), f(r.ecc_reached_pct_min), f(r.ecc_images_median, 0)] for r in s1.itertuples()],
                       ["r"] * 8))
    if len(s2):
        T.append(table("Post hoc cloud-probability sensitivity (amendment A4): H2 difference (adaptive minus "
                       "open-loop on expected coverage, same images, coverage points) with the protocol's constant "
                       "$p = 0.35$ and with the EarthEnv monthly climatology per tile (Wilson and Jetz, 2016).",
                       "tab:a4_cloud", ["Window", "$p$ climatology (mean)", "H2, $p=0.35$", "H2, climatology",
                                        "Adaptive (\\%)", "Open-loop exp. (\\%)", "Oracle (\\%)"],
                       [[r.window, f(r.p_earthenv_mean, 3), f(r.h2_const, 2), f(r.h2_earthenv, 2),
                         f(r.earthenv_cov_adaptive), f(r.earthenv_cov_openloop_expected), f(r.earthenv_cov_oracle)]
                        for r in s2.itertuples()], ["l"] + ["r"] * 6))
    if len(s3):
        T.append(table("Post hoc exact optimum on the full 48-h instance (amendment A4): images ECC's greedy cover "
                       "needs to reach 95\\% of reachable tiles against the exact minimum (integer programme, "
                       "HiGHS via SciPy).", "tab:a4_ilp48",
                       ["Window", "Pairs", "ECC greedy", "Exact minimum", "Status", "Solve (s)"],
                       [[r.window, f(r.pairs, 0), str(r.ecc_greedy_images), str(r.ecc_exact_images),
                         str(r.ecc_exact_status), f(r.seconds_max, 1)] for r in s3.itertuples()], ["l"] + ["r"] * 5))
    open(os.path.join(out, "a4_tables.tex"), "w").write("\n\n".join(T) + "\n")
    pd.set_option("display.width", 220)
    print(f"A4 sensitivity (POST HOC): {len(df)} runs; variants {summary['variants']}")
    print("\nS1 budget\n" + (s1.round(2).to_string(index=False) if len(s1) else "  none"))
    print("\nS2 cloud probability\n" + (s2.round(3).to_string(index=False) if len(s2) else "  none"))
    print("\nS3 exact optimum, full window, and determinism\n" + (s3.to_string(index=False) if len(s3) else "  none"))
    print(f"\nFiles: {out}/a4_runs.csv, a4_budget.csv, a4_cloud.csv, a4_ilp48.csv, a4_summary.json, a4_tables.tex")
    return df, summary


if __name__ == "__main__":
    main()
