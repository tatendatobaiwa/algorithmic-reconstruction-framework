# =============================================================================
# ARF PAPER 1 EVIDENCE (read-only; for the manuscript, labelled POST HOC where
# it goes beyond the frozen protocol's decision rules)
# =============================================================================
# Reads the run folders of the primary confirmatory batch (release 1.5.3; the
# index confirmatory_v153_primary_index.csv) and gathers, per window, what a
# referee will ask for beyond the three decisions, from the files each run
# already wrote. Nothing is recomputed from the orbits; nothing is written into
# the run folders.
#   E1 Baselines at the H1a level (95% of reachable tiles, clear sky): median
#      over onsets of the time to the level and of the images used, for the
#      physical floor, ECC, FCFS (unbudgeted), coverage-first and time-aware
#      greedy, with the share of onsets at which each reaches the level
#      (arf_onset_robustness.csv)
#   E2 Exact optimum (arf_ilp_gap.csv, Cell 7C sub-problems): ECC's greedy
#      cover against the exact minimum number of images; greedy coverage and
#      FCFS against the exact maximum coverage (budget ILP_BUDGET_K); solver time
#   E3 H2 bracket (arf_cloud_robustness_summary.json): coverage under cloud of
#      the adaptive (Fail-Signal) selector, the open-loop plan on expected
#      coverage (cloud-aware, no observations), the open-loop clear-sky plan
#      and the cloud oracle, same images; the share of the oracle's gain
#      captured; and the H2 difference across the cloud-scale sweep
#   E4 Archive validation (arf_acquisition_validation_summary.json and .csv,
#      Cell 7D): precision, recall and time offset per platform; observed scene
#      cloud cover of the real acquisitions against the model's probability
#   E5 Window-level exact tests (POST HOC complement to the protocol's
#      Student-t bounds): one-sided sign test and Wilcoxon signed-rank test over
#      the windows, window = unit of inference; with 8 windows the smallest
#      attainable one-sided p of either is 1/256
# Every window outcome is taken from arf_confirmatory_analysis.window_outcomes,
# the protocol's own reader, and the script stops if its H1a/H1b/H2 differ
# from the run's own summaries.
# Writes CSV files, a JSON summary and LaTeX tables (fully enclosed grid) to
# <BATCH_DIR>/paper1_evidence/.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_paper1_evidence.py
# =============================================================================
import os
import json
import math
import runpy
import argparse

import numpy as np
import pandas as pd

LEVEL = "95"   # the H1a level of the protocol (h1_reachable_coverage_pct = 95)
POLICIES = (("best_possible", "Physical floor (all passes)", False),
            ("earliest_completion", "ECC ($k=30$)", True),
            ("fcfs_unbudgeted", "FCFS", True),
            ("coverage_first", "Coverage-first greedy ($k=30$)", True),
            ("time_aware", "Time-aware greedy ($k=30$)", True))


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _med(s):
    s = _num(s).dropna()
    return float(s.median()) if s.size else None


def baselines(run_dir, sky="clear"):
    on = pd.read_csv(os.path.join(run_dir, "arf_onset_robustness.csv"))
    ok = _num(on.get(f"reachable_{sky}_tiles")).fillna(0) > 0
    on = on[ok]
    out = {"onsets": int(len(on))}
    for p, _, has_images in POLICIES:
        t = _num(on.get(f"{p}_{sky}_t{LEVEL}r_min"))
        out[f"{p}_t_median_min"] = _med(t)
        out[f"{p}_reached_pct"] = float(100.0 * t.notna().mean()) if len(on) else None
        if has_images:
            out[f"{p}_images_median"] = _med(on.get(f"{p}_{sky}_images{LEVEL}r"))
    out["random_cov_pct_median"] = _med(on.get(f"random_{sky}_cov_pct"))
    return out


def optimality(run_dir):
    p = os.path.join(run_dir, "arf_ilp_gap.csv")
    if not os.path.exists(p):
        return {"ilp_instances": 0}
    g = pd.read_csv(p)
    opt = g[g["ilp_status"] == "optimal"]
    e = g[(g.get("ecc_exact_status") == "optimal") & g["ecc_greedy_images"].notna()] if "ecc_exact_status" in g else g.iloc[:0]
    extra = (_num(e["ecc_greedy_images"]) - _num(e["ecc_exact_images"])) if len(e) else pd.Series(dtype=float)
    out = {"ilp_instances": int(len(g)), "ilp_optimal": int(len(opt)),
           "ecc_exact_instances": int(len(e)), "ecc_matches_exact": int((extra == 0).sum()),
           "ecc_extra_images_max": float(extra.max()) if extra.size else None,
           "ecc_exact_images_sum": float(_num(e["ecc_exact_images"]).sum()) if len(e) else None,
           "ecc_greedy_images_sum": float(_num(e["ecc_greedy_images"]).sum()) if len(e) else None,
           "ilp_seconds_median": _med(g["ilp_seconds"]), "ilp_seconds_max": float(_num(g["ilp_seconds"]).max()),
           "horizons_min": sorted(int(h) for h in g["horizon_min"].unique())}
    for n in ("coverage_first", "time_aware", "fcfs"):
        r = _num(opt[f"{n}_ratio"]).dropna()
        out[f"{n}_ratio_min"] = float(r.min()) if r.size else None
        out[f"{n}_ratio_median"] = float(r.median()) if r.size else None
    r = _num(g.get("time_aware_tw_ratio")).dropna() if "time_aware_tw_ratio" in g else pd.Series(dtype=float)
    out["time_aware_tw_ratio_min"] = float(r.min()) if r.size else None
    return out


def h2_bracket(run_dir):
    cr = json.load(open(os.path.join(run_dir, "arf_cloud_robustness_summary.json")))
    h = cr["h2"]["time_aware"]
    cov = h.get("coverage_pct_mean") or {}
    dx, dc = h.get("difference_vs_openloop_expected") or {}, h.get("difference_vs_openloop_clear") or {}
    out = {"h2_evaluable": bool(h.get("evaluable")), "realisations": h.get("n_realisations"),
           "images_mean": h.get("images_mean"),
           **{f"cov_{k}": cov.get(k) for k in ("adaptive", "openloop_expected", "openloop_clear", "oracle")},
           "diff_vs_expected_mean": dx.get("mean"), "diff_vs_expected_lcb95": dx.get("lcb"),
           "diff_vs_expected_pct_negative": dx.get("pct_negative"),
           "diff_vs_clear_mean": dc.get("mean"), "diff_vs_clear_lcb95": dc.get("lcb"),
           "share_oracle_gain_median": h.get("share_of_oracle_gain_captured_median"),
           "replans_total": h.get("replans_total"), "compensatory_total": h.get("compensatory_total"),
           "cloud_probability": ((cr.get("cloud_model") or {}).get("probability") or {}).get("value"),
           "cloud_probability_source": ((cr.get("cloud_model") or {}).get("probability") or {}).get("source")}
    sweep = []
    for s in cr.get("sensitivity") or []:
        sweep.append({"model": s.get("model"), "length_km": s.get("length_km"), "time_h": s.get("time_h"),
                      "primary": s.get("primary"), "diff_mean": s.get("diff_mean"), "diff_lcb95": s.get("diff_lcb95"),
                      "realisations": s.get("n_realisations")})
    return out, sweep


def archive(run_dir):
    p = os.path.join(run_dir, "arf_acquisition_validation_summary.json")
    if not os.path.exists(p):
        return {"archive_status": "missing"}, [], []
    v = json.load(open(p))
    out = {"archive_status": v.get("status"), **{f"archive_{k}": v.get(k) for k in (
        "predicted_overpasses", "actual_overpasses", "matched", "precision", "recall", "median_time_offset_min",
        "predicted_tiles_confirmed_pct", "actual_tiles_predicted_pct", "actual_cloud_cover_mean_pct",
        "actual_cloud_cover_median_pct", "model_cloud_probability_pct")}}
    plats = [{"platform": k, **b} for k, b in (v.get("by_platform") or {}).items()]
    rows = []
    c = os.path.join(run_dir, "arf_acquisition_validation.csv")
    if os.path.exists(c):
        rows = pd.read_csv(c).to_dict("records")
    return out, plats, rows


def exact_tests(x, ref=0.0):
    """One-sided tests that the window values exceed `ref`: sign test (ties dropped) and Wilcoxon
    signed-rank (exact; zero differences dropped). With all differences tied at zero, no test."""
    from scipy.stats import binomtest, wilcoxon
    d = np.asarray([v - ref for v in x if v is not None and np.isfinite(v)], float)
    nz = d[d != 0]
    out = {"windows": int(d.size), "above": int((d > 0).sum()), "equal": int((d == 0).sum()), "below": int((d < 0).sum())}
    if nz.size == 0:
        return {**out, "sign_p": None, "wilcoxon_p": None, "min_attainable_p": None}
    out["sign_p"] = float(binomtest(int((nz > 0).sum()), int(nz.size), 0.5, alternative="greater").pvalue)
    out["wilcoxon_p"] = float(wilcoxon(nz, alternative="greater", method="exact").pvalue) if nz.size >= 1 else None
    out["min_attainable_p"] = 0.5 ** int(nz.size)
    return out


def fmt(v, nd=1, pct=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    if isinstance(v, (int, np.integer)) and nd == 0:
        return f"{int(v):,}"
    return f"{float(v):,.{nd}f}"


def grid_table(caption, label, header, rows, align=None, star=True):
    env = "table*" if star else "table"
    align = align or ["l"] + ["r"] * (len(header) - 1)
    L = [f"\\begin{{{env}}}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         f"\\caption{{{caption}}}", f"\\label{{{label}}}", "\\begin{tabular}{|" + "|".join(align) + "|}", "\\hline",
         " & ".join(header) + " \\\\", "\\hline"]
    for r in rows:
        L += [" & ".join(r) + " \\\\", "\\hline"]
    L += ["\\end{tabular}", f"\\end{{{env}}}"]
    return "\n".join(L)


def rng_str(vals, nd=0):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    if not v:
        return "--"
    med = float(np.median(v))
    return f"{fmt(med, nd)} [{fmt(min(v), nd)}--{fmt(max(v), nd)}]"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Read-only manuscript evidence from the primary confirmatory runs.")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    ap.add_argument("--index", default="confirmatory_v153_primary_index.csv")
    a, _ = ap.parse_known_args(argv)
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    an = runpy.run_path(os.path.join(a.code, "arf_confirmatory_analysis.py"), run_name="arf_confirmatory_analysis")
    idx = pd.read_csv(os.path.join(a.batch_dir, a.index), dtype=str).fillna("")
    idx = idx[(idx.status == "done") & (idx.release == proto["release"])]
    W, SW, PL, SC = [], [], [], []
    for _, r in idx.iterrows():
        d = r["out_dir"]
        w = an["window_outcomes"](d, int(proto["min_evaluable_onsets"]))
        if str(w.get("analysis_epoch")) not in proto["windows"]:
            continue
        b = baselines(d)
        o = optimality(d)
        h, sweep = h2_bracket(d)
        av, plats, scenes = archive(d)
        # consistency with the protocol reader
        assert (w["h2_diff"] is None) == (not h["h2_evaluable"]) and (
            w["h2_diff"] is None or abs(w["h2_diff"] - h["diff_vs_expected_mean"]) < 1e-12), f"{d}: H2 differs"
        win = str(w["analysis_epoch"])[:10]
        W.append({"window": win, "h1a_pct": w["h1a_pct"], "h1b_pct": w["h1b_pct"], "h2_diff": w["h2_diff"],
                  "run_minutes": float(r["minutes"]) if r.get("minutes") else None, **b, **o, **h, **av, "run_dir": d})
        SW += [{"window": win, **s} for s in sweep]
        PL += [{"window": win, **p} for p in plats]
        SC += [{"window": win, **s} for s in scenes if s.get("kind") in ("matched", "actual_only")]
    if not W:
        raise SystemExit(f"no completed run of release {proto['release']} on a protocol window in {a.index}")
    df = pd.DataFrame(W).sort_values("window")
    out = os.path.join(a.batch_dir, "paper1_evidence")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "p1_windows.csv"), index=False)
    sw = pd.DataFrame(SW)
    sw.to_csv(os.path.join(out, "p1_cloud_scale_sweep.csv"), index=False)
    pl = pd.DataFrame(PL)
    pl.to_csv(os.path.join(out, "p1_archive_by_platform.csv"), index=False)
    sc = pd.DataFrame(SC)
    sc.to_csv(os.path.join(out, "p1_archive_scenes.csv"), index=False)

    th = float(proto["threshold_pct"])
    tests = {"note": "POST HOC: exact one-sided tests over windows (window = unit of inference); the protocol's "
                     "decisions are the Student-t lower bounds of arf_confirmatory_analysis.py",
             "H1a_vs_threshold": exact_tests(df.h1a_pct.tolist(), th),
             "H1b_vs_threshold": exact_tests(df.h1b_pct.tolist(), th),
             "H2_vs_zero": exact_tests(df.h2_diff.tolist(), 0.0)}
    # pooled archive: per platform over windows
    pp = []
    if len(pl):
        for plat, g in pl.groupby("platform"):
            n_p, n_r, n_m = (int(_num(g[c]).fillna(0).sum()) for c in ("predicted_overpasses", "actual_overpasses", "matched"))
            pp.append({"platform": plat, "windows_evaluated": int(g["evaluated"].astype(str).str.lower().eq("true").sum()),
                       "predicted": n_p, "actual": n_r, "matched": n_m,
                       "precision": n_m / n_p if n_p else None, "recall": n_m / n_r if n_r else None,
                       "archive_lag_h": _med(g.get("archive_lag_h"))})
    pp = pd.DataFrame(pp)
    pp.to_csv(os.path.join(out, "p1_archive_pooled.csv"), index=False)
    cloud = {}
    if len(sc) and "cloud_pct" in sc:
        c = _num(sc["cloud_pct"]).dropna()
        cloud = {"scenes_with_cloud": int(c.size), "mean_pct": float(c.mean()) if c.size else None,
                 "median_pct": float(c.median()) if c.size else None,
                 "share_below_10pct": float(100 * (c < 10).mean()) if c.size else None,
                 "share_above_50pct": float(100 * (c > 50).mean()) if c.size else None,
                 "model_probability_pct": _med(df.get("archive_model_cloud_probability_pct")),
                 "note": "observed: real acquisitions' scene cloud cover from the archive catalogue (Cell 7D); "
                         "model: the protocol's constant per-tile cloud probability"}
    summary = {"release": proto["release"], "windows": df.window.tolist(), "tests": tests,
               "archive_cloud": cloud,
               "optimality_pooled": {
                   "ecc_instances": int(df.ecc_exact_instances.sum()), "ecc_matches": int(df.ecc_matches_exact.sum()),
                   "ecc_extra_images_max": float(np.nanmax(df.ecc_extra_images_max.astype(float))) if df.ecc_extra_images_max.notna().any() else None,
                   **{f"{n}_ratio_min": (float(np.nanmin(df[f"{n}_ratio_min"].astype(float))) if df[f"{n}_ratio_min"].notna().any() else None)
                      for n in ("coverage_first", "time_aware", "fcfs")},
                   "ilp_seconds_max": float(df.ilp_seconds_max.max())},
               "h2_pooled": {k: rng_str(df[k].astype(float).tolist(), 2) for k in (
                   "cov_adaptive", "cov_openloop_expected", "cov_openloop_clear", "cov_oracle", "diff_vs_expected_mean",
                   "diff_vs_clear_mean", "share_oracle_gain_median")}}
    json.dump(summary, open(os.path.join(out, "p1_summary.json"), "w"), indent=1, default=str)

    # ---- LaTeX tables -------------------------------------------------------
    T = []
    rows = []
    for p, lab, has_images in POLICIES:
        rows.append([lab, rng_str(df[f"{p}_t_median_min"].astype(float).tolist(), 0),
                     rng_str(df[f"{p}_images_median"].astype(float).tolist(), 0) if has_images else "all",
                     rng_str(df[f"{p}_reached_pct"].astype(float).tolist(), 1)])
    T.append(grid_table(
        "Selectors at the H1a level (95\\% of reachable tiles, clear sky), primary design, eight windows. Each cell: "
        "median over windows of the per-window median over onsets [range over windows]. FCFS: first-come-first-served, "
        "unbudgeted.", "tab:p1_baselines",
        ["Selector", "Time to level (min)", "Images", "Onsets reaching level (\\%)"], rows, ["l", "r", "r", "r"]))
    po = summary["optimality_pooled"]
    rows = [["Sub-problems where ECC's greedy cover uses the exact minimum of images", f"{po['ecc_matches']}/{po['ecc_instances']}"],
            ["Largest excess of ECC over the exact minimum", fmt(po["ecc_extra_images_max"], 0) + " images"],
            ["Coverage-first greedy / exact maximum coverage (worst)", fmt(po["coverage_first_ratio_min"], 3)],
            ["Time-aware greedy / exact maximum coverage (worst)", fmt(po["time_aware_ratio_min"], 3)],
            ["FCFS / exact maximum coverage (worst)", fmt(po["fcfs_ratio_min"], 3)],
            ["Longest exact solve (s)", fmt(po["ilp_seconds_max"], 1)]]
    T.append(grid_table("Greedy selection against the exact integer optimum on the Cell 7C sub-problems "
                        "(6-h horizons, clear sky), pooled over the eight primary windows.", "tab:p1_ilp",
                        ["Quantity", "Value"], rows, ["p{2.45in}", "r"], star=False))
    rows = []
    for _, r in df.iterrows():
        rows.append([r.window, fmt(r.cov_adaptive, 1), fmt(r.cov_openloop_expected, 1), fmt(r.cov_openloop_clear, 1),
                     fmt(r.cov_oracle, 1), fmt(r.diff_vs_expected_mean, 2), fmt(r.diff_vs_expected_lcb95, 2),
                     fmt(None if r.share_oracle_gain_median is None else 100 * float(r.share_oracle_gain_median), 0)])
    T.append(grid_table(
        "Coverage under cloud (\\%), same number of images, time-aware objective, 20 realisations per window "
        "(correlated cloud, $L=50$ km, $T=3$ h, $p=0.35$). Open-loop (exp.): the same planner on expected coverage "
        "without observations; oracle: the greedy told every image's cloud in advance. $\\Delta$: adaptive minus "
        "open-loop (exp.), with its one-sided 95\\% lower bound over realisations.", "tab:p1_h2",
        ["Window", "Adaptive", "Open-loop (exp.)", "Open-loop (clear)", "Oracle", "$\\Delta$ (pts)", "LB (pts)",
         "Oracle gain (\\%)"], rows))
    if len(pp):
        rows = [[str(r.platform), fmt(r.predicted, 0), fmt(r.actual, 0), fmt(r.matched, 0), fmt(r.precision, 2),
                 fmt(r.recall, 2)] for r in pp.itertuples()]
        T.append(grid_table("Predicted against archived acquisitions, pooled over the eight primary windows "
                            "(Cell 7D; PlanetScope not validated, no open archive). Pred.: predicted overpasses; "
                            "Arch.: archived; Match.: matched within 15 min; Prec.: precision.", "tab:p1_archive",
                            ["Platform", "Pred.", "Arch.", "Match.", "Prec.", "Recall"], rows,
                            ["l", "r", "r", "r", "r", "r"], star=False))
    rows = []
    for h, k in (("H1a", "H1a_vs_threshold"), ("H1b", "H1b_vs_threshold"), ("H2", "H2_vs_zero")):
        t = tests[k]
        rows.append([h, f"{t['above']}/{t['equal']}/{t['below']}", fmt(t["sign_p"], 4), fmt(t["wilcoxon_p"], 4)])
    T.append(grid_table(f"Post hoc exact one-sided tests over the eight windows (window = unit of inference): H1a "
                        f"and H1b against {th:g}\\%, H2 against 0. Counts: windows above / equal / below; equal windows "
                        "are dropped, so all-equal gives no test.", "tab:p1_exact",
                        ["Hyp.", "Above/eq./below", "Sign $p$", "Wilcoxon $p$"], rows,
                        ["l", "r", "r", "r"], star=False))
    open(os.path.join(out, "p1_tables.tex"), "w").write("\n\n".join(T) + "\n")

    pd.set_option("display.width", 220)
    print(f"Paper 1 evidence, release {proto['release']}, {len(df)} windows ({', '.join(df.window)})\n")
    print("E1 selectors at the H1a level (median over onsets per window)")
    print(df[["window"] + [f"{p}_t_median_min" for p, _, _ in POLICIES] +
             [f"{p}_images_median" for p, _, h in POLICIES if h]].round(1).to_string(index=False))
    print("\nE2 exact optimum:", json.dumps(po))
    print("\nE3 H2 bracket")
    print(df[["window", "cov_adaptive", "cov_openloop_expected", "cov_openloop_clear", "cov_oracle",
              "diff_vs_expected_mean", "diff_vs_expected_lcb95", "share_oracle_gain_median"]].round(2).to_string(index=False))
    if len(sw):
        print("\n   H2 difference across the cloud-scale sweep (median over windows)")
        print(sw.groupby(["model", "length_km"], dropna=False)["diff_mean"].median().round(2).to_string())
    print("\nE4 archive:")
    print(pp.round(3).to_string(index=False) if len(pp) else "   no archive validation found")
    print("   observed scene cloud vs model:", json.dumps(cloud))
    print("\nE5 exact tests (post hoc):", json.dumps(tests, indent=1))
    print(f"\nFiles: {out}/p1_windows.csv, p1_cloud_scale_sweep.csv, p1_archive_by_platform.csv, p1_archive_pooled.csv, "
          "p1_archive_scenes.csv, p1_summary.json, p1_tables.tex")
    return df, summary


if __name__ == "__main__":
    main()
