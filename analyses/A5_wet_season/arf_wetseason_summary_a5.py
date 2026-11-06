# =============================================================================
# ARF WET-SEASON SUMMARY (protocol amendment A5, EXPLORATORY; read-only)
# =============================================================================
# Reads the runs of batch wetseason_v153_a5 and tabulates, per cloud setting
# (p035: the protocol's constant p = 0.35; earthenv: monthly climatology) and
# window, the protocol's own window outcomes (H1a, H1b, H2 via
# arf_confirmatory_analysis.window_outcomes), the cloud probability used, the
# four-plan coverage bracket under cloud, and the archive check of Cell 7D
# (matched overpasses and the OBSERVED scene cloud cover of the real
# acquisitions). Over the eight windows: mean, SD, the one-sided Student-t lower
# bounds at alpha = 0.05 and alpha/3 (the protocol's rule, applied
# descriptively), and exact one-sided sign and Wilcoxon tests. The dry-season
# values (primary confirmatory windows; A4 climatology runs) are placed beside
# them when present. No hypothesis decision is made (exploratory arm).
# Writes ARF_batch/wetseason_a5/ (CSV, JSON, LaTeX table with a fully enclosed grid).
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_wetseason_summary_a5.py
# =============================================================================
import os
import json
import math
import runpy
import argparse

import numpy as np
import pandas as pd


def lcb(x, alpha):
    from scipy.stats import t
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    if x.size < 2:
        return None
    return float(x.mean() - t.ppf(1 - alpha, x.size - 1) * x.std(ddof=1) / math.sqrt(x.size))


def exact(x, ref):
    from scipy.stats import binomtest, wilcoxon
    d = np.asarray([v - ref for v in x if v is not None and np.isfinite(v)], float)
    nz = d[d != 0]
    if nz.size == 0:
        return {"above": int((d > 0).sum()), "equal": int((d == 0).sum()), "below": int((d < 0).sum()),
                "sign_p": None, "wilcoxon_p": None}
    return {"above": int((d > 0).sum()), "equal": int((d == 0).sum()), "below": int((d < 0).sum()),
            "sign_p": float(binomtest(int((nz > 0).sum()), int(nz.size), 0.5, alternative="greater").pvalue),
            "wilcoxon_p": float(wilcoxon(nz, alternative="greater", method="exact").pvalue)}


def run_row(d, an, min_onsets):
    w = an["window_outcomes"](d, min_onsets)
    cr = json.load(open(os.path.join(d, "arf_cloud_robustness_summary.json")))
    h = cr["h2"]["time_aware"]
    cov = h.get("coverage_pct_mean") or {}
    prob = (cr.get("cloud_model") or {}).get("probability") or {}
    row = {"window": str(w["analysis_epoch"])[:10], "h1a_pct": w["h1a_pct"], "h1b_pct": w["h1b_pct"],
           "h2_diff": w["h2_diff"], "h2_evaluable": w["h2_evaluable"],
           "p_source": prob.get("source"), "p_mean": prob.get("mean", prob.get("value")),
           **{f"cov_{k}": cov.get(k) for k in ("adaptive", "openloop_expected", "openloop_clear", "oracle")},
           "archive_status": w.get("archive_status"), "archive_matched": w.get("archive_matched"),
           "archive_precision": w.get("archive_precision"), "archive_recall": w.get("archive_recall"),
           "run_dir": d}
    c = os.path.join(d, "arf_acquisition_validation.csv")
    if os.path.exists(c):
        s = pd.read_csv(c)
        s = s[s.kind.isin(["matched", "actual_only"])]
        cc = pd.to_numeric(s.cloud_pct, errors="coerce").dropna()
        row.update(observed_scenes=int(cc.size), observed_cloud_mean_pct=float(cc.mean()) if cc.size else None,
                   observed_cloud_median_pct=float(cc.median()) if cc.size else None)
    return row


def pooled(df, alpha, th):
    out = {}
    for col, ref in (("h1a_pct", th), ("h1b_pct", th), ("h2_diff", 0.0)):
        x = pd.to_numeric(df[col], errors="coerce").dropna().tolist()
        out[col] = {"windows": len(x), "mean": float(np.mean(x)) if x else None,
                    "sd": float(np.std(x, ddof=1)) if len(x) > 1 else None,
                    "lcb_alpha": lcb(x, alpha), "lcb_alpha_over_3": lcb(x, alpha / 3), **exact(x, ref)}
    return out


def f(v, nd=1):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "--"
    return "--" if not np.isfinite(v) else f"{v:,.{nd}f}"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Wet-season arm summary (A5, exploratory).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    a, _ = ap.parse_known_args(argv)
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    an = runpy.run_path(os.path.join(a.code, "arf_confirmatory_analysis.py"), run_name="arf_confirmatory_analysis")
    idx = pd.read_csv(os.path.join(a.batch_dir, "wetseason_v153_a5_index.csv"), dtype=str).fillna("")
    rows = []
    for _, r in idx[(idx.status == "done") & (idx.release == proto["release"])].iterrows():
        cfg = json.loads(r["config"]) if r.get("config") else {}
        rows.append({"variant": "earthenv" if cfg.get("cloud_probability_source") == "earthenv" else "p035",
                     **run_row(r["out_dir"], an, int(proto["min_evaluable_onsets"]))})
    df = pd.DataFrame(rows).sort_values(["variant", "window"])
    out = os.path.join(a.batch_dir, "wetseason_a5")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "a5_runs.csv"), index=False)
    alpha, th = float(proto["alpha"]), float(proto["threshold_pct"])
    summ = {"label": "EXPLORATORY (amendment A5): no hypothesis decision", "runs": int(len(df)),
            "by_variant": {v: pooled(g, alpha, th) for v, g in df.groupby("variant")}}
    obs = pd.to_numeric(df.get("observed_cloud_mean_pct"), errors="coerce")
    summ["observed_scene_cloud_mean_pct_over_windows"] = float(obs.dropna().mean()) if obs.notna().any() else None
    # dry season beside it, when the records exist
    dry = {}
    p = os.path.join(a.batch_dir, "arf_confirmatory_results_windows.csv")
    if os.path.exists(p):
        c = pd.read_csv(p)
        c = c[(c.roi == "ckgr_okavango") & (c.fleet == "earth-resources+planet") &
              (c.planet_swath_case.fillna("conservative") == "conservative") & (c.exclude_satellites.fillna("") == "")]
        dry["p035"] = pooled(c, alpha, th)
    p = os.path.join(a.batch_dir, "sensitivity_a4", "a4_cloud.csv")
    if os.path.exists(p):
        c = pd.read_csv(p).rename(columns={"h2_earthenv": "h2_diff"})
        dry["earthenv_h2"] = pooled(c.assign(h1a_pct=np.nan, h1b_pct=np.nan), alpha, th)["h2_diff"]
    summ["dry_season"] = dry
    json.dump(summ, open(os.path.join(out, "a5_summary.json"), "w"), indent=1, default=str)

    # LaTeX: per window, both cloud settings
    L = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         "\\caption{Exploratory wet-season arm (amendment A5): primary design, eight 48-h windows in January and "
         "February 2026, run retrospectively with release 1.5.3. H1a, H1b: share of onsets (\\%). H2: adaptive "
         "minus open-loop coverage under cloud, same images (points), with the protocol's constant $p = 0.35$ and "
         "with the EarthEnv monthly climatology ($\\bar p$: its mean over the tiles). Observed cloud: mean scene "
         "cloud cover of the real Sentinel-2 and Landsat acquisitions in the window (Cell 7D).}",
         "\\label{tab:a5_wet}", "\\begin{tabular}{|l|r|r|r|r|r|r|}", "\\hline",
         "Window & H1a & H1b & H2, $p=0.35$ & $\\bar p$ climatology & H2, climatology & Observed cloud (\\%) \\\\",
         "\\hline"]
    wide = df.pivot_table(index="window", columns="variant",
                          values=["h1a_pct", "h1b_pct", "h2_diff", "p_mean", "observed_cloud_mean_pct"], aggfunc="first")
    g = lambda w, col, v: wide.loc[w, (col, v)] if (col, v) in wide.columns else None
    for w in wide.index:
        L += [" & ".join([w, f(g(w, "h1a_pct", "p035")), f(g(w, "h1b_pct", "p035")), f(g(w, "h2_diff", "p035"), 2),
                          f(g(w, "p_mean", "earthenv"), 3), f(g(w, "h2_diff", "earthenv"), 2),
                          f(g(w, "observed_cloud_mean_pct", "p035"), 1)]) + " \\\\", "\\hline"]
    bv = summ["by_variant"]
    if "p035" in bv:
        L += [" & ".join(["Mean", f(bv["p035"]["h1a_pct"]["mean"]), f(bv["p035"]["h1b_pct"]["mean"]),
                          f(bv["p035"]["h2_diff"]["mean"], 2), "--",
                          f(bv.get("earthenv", {}).get("h2_diff", {}).get("mean"), 2),
                          f(summ["observed_scene_cloud_mean_pct_over_windows"], 1)]) + " \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "a5_table.tex"), "w").write("\n".join(L) + "\n")

    pd.set_option("display.width", 220)
    print(f"A5 wet-season arm (EXPLORATORY): {len(df)} runs")
    print(df[["variant", "window", "h1a_pct", "h1b_pct", "h2_diff", "p_mean", "cov_adaptive", "cov_openloop_expected",
              "cov_oracle", "archive_matched", "observed_scenes", "observed_cloud_mean_pct"]].round(3).to_string(index=False))
    print("\nOver the eight windows (descriptive; alpha = %.2f, alpha/3 = %.4f):" % (alpha, alpha / 3))
    print(json.dumps(summ["by_variant"], indent=1, default=str))
    print("\nDry season beside it:", json.dumps(dry, indent=1, default=str)[:3000])
    print(f"\nFiles: {out}/a5_runs.csv, a5_summary.json, a5_table.tex")
    return df, summ


if __name__ == "__main__":
    main()
