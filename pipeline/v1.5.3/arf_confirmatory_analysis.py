# =============================================================================
# ARF CONFIRMATORY ANALYSIS (v1.5.2; unchanged in 1.5.3): the frozen protocol's decision rules
# =============================================================================
# Reads the run folders of a batch (arf_batch_runner.py) and applies the
# decision rules of arf_confirmatory_protocol.json, which were fixed before any
# of these windows was run. Nothing here is tuned to the data.
#
# Unit of analysis: the analysis WINDOW. Onsets within a window share most of
# their passes and cloud realisations share the window's passes, so neither is
# an independent replicate; each window contributes one number per hypothesis:
#   H1a_w  % of evaluable onsets where ECC (k = 30) reaches 95% of reachable
#          tiles within one orbital period of the physical floor (clear sky)
#   H1b_w  % of eligible onsets (FCFS needs >= 2 images) where ECC does so with
#          fewer images than first-come-first-served
#   H2_w   mean over cloud realisations of (adaptive - open-loop) coverage
#          points with the same images (time-aware objective, primary cloud model)
# A window counts for H1a/H1b only with >= min_evaluable_onsets evaluated
# onsets, for H2 only when H2 was evaluable in it. A region-fleet design is NOT
# EVALUABLE with fewer than min_windows evaluable windows. Otherwise:
#   H1a, H1b SUPPORTED if the one-sided lower confidence bound (Student t over
#            windows) of the mean of H1x_w is >= threshold_pct
#   H2       SUPPORTED if the one-sided lower confidence bound of the mean of
#            H2_w is > 0
# Primary design: the protocol's three primary hypotheses at a Bonferroni-
# adjusted level (alpha / 3). Secondary and exploratory designs: unadjusted,
# labelled exploratory.
# v1.5.2: designs are also matched on exclude_satellites (amendment 1 of the
# protocol adds exploratory designs without Sentinel-2A), so a run without a
# satellite is never pooled with the full-fleet runs of the same region and
# fleet; per-platform archive validation (Cell 7D) is carried into the window
# table.
# v1.5.3: no change to the analysis. The release of record is read from the
# protocol file (amendment A2: 1.5.3), and runs of any other release are
# excluded and listed.
# =============================================================================
import os
import json
import math
import argparse

import numpy as np
import pandas as pd


def _t_lcb(x, alpha):
    """Mean, SD and one-sided (1 - alpha) lower confidence bound of the mean (Student t)."""
    from scipy.stats import t as _t
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], dtype=float)
    n = int(x.size)
    if n == 0:
        return n, None, None, None
    m = float(x.mean())
    if n == 1:
        return n, m, None, None
    sd = float(x.std(ddof=1))
    return n, m, sd, m - float(_t.ppf(1.0 - alpha, n - 1)) * sd / math.sqrt(n)


def _excl_key(names):
    """Canonical key of an exclude_satellites list ('' = full fleet)."""
    return ",".join(sorted(str(n).strip().upper() for n in (names or [])))


def window_outcomes(run_dir, min_evaluable_onsets):
    """The three window-level outcomes of one run folder (None when not evaluable)."""
    out = {"run_dir": run_dir}
    man = json.load(open(os.path.join(run_dir, "arf_reproducibility_manifest.json")))
    cfg = man.get("config_overrides") or {}
    out.update(release=man.get("framework_version"), roi=cfg.get("active_roi"), fleet=cfg.get("fleet"),
               planet_swath_case=cfg.get("planet_swath_case"), analysis_epoch=cfg.get("analysis_epoch"),
               exclude_satellites=_excl_key(cfg.get("exclude_satellites")),
               versions_consistent=(man.get("cell_versions") or {}).get("consistent"))
    ons = json.load(open(os.path.join(run_dir, "arf_onset_robustness_summary.json")))
    a = ons["h1_revised"]["series"]["earliest_completion_clear"]
    b = ons["h1b_image_economy"]["series"]["clear"]
    out["h1a_onsets"] = a["n_onsets_evaluated"]
    out["h1a_pct"] = a["pct_within_margin"] if (a["n_onsets_evaluated"] or 0) >= min_evaluable_onsets else None
    out["h1b_onsets"] = b["n_onsets_evaluated"]
    out["h1b_pct"] = b["pct_success"] if (b["n_onsets_evaluated"] or 0) >= min_evaluable_onsets else None
    out["h1b_median_images_ecc"], out["h1b_median_images_fcfs"] = b.get("median_images_ecc"), b.get("median_images_fcfs")
    cr = json.load(open(os.path.join(run_dir, "arf_cloud_robustness_summary.json")))
    h2 = cr["h2"]["time_aware"]
    out["h2_evaluable"] = bool(h2.get("evaluable"))
    out["h2_diff"] = h2["difference_vs_openloop_expected"]["mean"] if h2.get("evaluable") else None
    out["h2_realisations"] = h2.get("n_realisations")
    out["h2_cloud_model"] = (cr.get("cloud_model") or {}).get("model")
    av = os.path.join(run_dir, "arf_acquisition_validation_summary.json")
    if os.path.exists(av):
        v = json.load(open(av))
        out.update({f"archive_{k}": v.get(k) for k in ("status", "predicted_overpasses", "actual_overpasses",
                                                       "matched", "precision", "recall", "median_time_offset_min",
                                                       "predicted_tiles_confirmed_pct", "actual_tiles_predicted_pct",
                                                       "actual_cloud_cover_mean_pct")})
        for plat, bp in (v.get("by_platform") or {}).items():          # v1.5.2 (Cell 7D)
            out[f"archive_{plat}_precision"] = bp.get("precision")
            out[f"archive_{plat}_recall"] = bp.get("recall")
    return out


def decide(rows, design, protocol, alpha):
    th = float(protocol["threshold_pct"])
    mw = int(protocol["min_windows"])
    res = {"design": design, "alpha": alpha, "n_windows": len(rows)}
    for h, col, rule in (("H1a", "h1a_pct", "ge_threshold"), ("H1b", "h1b_pct", "ge_threshold"),
                         ("H2", "h2_diff", "gt_zero")):
        vals = [r[col] for r in rows if r.get(col) is not None]
        n, m, sd, lcb = _t_lcb(vals, alpha)
        if n < mw or lcb is None:          # a lower bound needs at least two windows
            st = "NOT EVALUABLE"
        elif rule == "ge_threshold":
            st = "SUPPORTED" if lcb is not None and lcb >= th else "NOT SUPPORTED"
        else:
            st = "SUPPORTED" if lcb is not None and lcb > 0 else "NOT SUPPORTED"
        res[h] = {"status": st, "evaluable_windows": n, "mean": m, "sd": sd, f"lcb_{1 - alpha:.4f}": lcb,
                  "rule": (f"one-sided {100 * (1 - alpha):.2f}% lower bound of the mean over windows >= {th:g}%"
                           if rule == "ge_threshold" else
                           f"one-sided {100 * (1 - alpha):.2f}% lower bound of the mean over windows > 0 points")}
    return res


def analyse(run_dirs, protocol):
    rows = [window_outcomes(d, int(protocol["min_evaluable_onsets"])) for d in run_dirs]
    df = pd.DataFrame(rows)
    results = {"protocol": protocol.get("id"), "frozen_utc": protocol.get("frozen_utc"),
               "amendments": [a.get("id") for a in protocol.get("amendments", [])],
               "release_expected": protocol.get("release"), "designs": []}
    if df.empty:
        return results, df
    bad = df[(df["release"] != protocol.get("release")) | (df["versions_consistent"] != True)]  # noqa: E712
    results["runs_excluded_release_or_mixed"] = bad["run_dir"].tolist()
    df = df.drop(bad.index)
    for d in protocol["designs"]:
        sel = df[(df.roi == d["active_roi"]) & (df.fleet == d["fleet"]) &
                 (df.planet_swath_case.fillna("conservative") == d.get("planet_swath_case", "conservative")) &
                 (df.exclude_satellites == _excl_key(d.get("exclude_satellites")))]
        sel = sel[sel.analysis_epoch.isin(protocol["windows"])]
        alpha = (float(protocol["alpha"]) / len(protocol["primary_hypotheses"]) if d["role"] == "primary"
                 else float(protocol["alpha"]))
        r = decide(sel.to_dict("records"), d, protocol, alpha)
        r["role"] = d["role"]
        r["windows_found"] = sorted(sel.analysis_epoch.tolist())
        r["windows_missing"] = sorted(set(protocol["windows"]) - set(sel.analysis_epoch))
        results["designs"].append(r)
    return results, df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protocol", required=True, help="arf_confirmatory_protocol.json")
    ap.add_argument("--index", nargs="+", help="batch index CSV(s) from arf_batch_runner.py")
    ap.add_argument("--run-dirs", nargs="*", default=[], help="run folders (instead of / in addition to --index)")
    ap.add_argument("--out", default="arf_confirmatory_results.json")
    a = ap.parse_args(argv)
    protocol = json.load(open(a.protocol))
    dirs = list(a.run_dirs)
    for p in a.index or []:
        idx = pd.read_csv(p, dtype=str).fillna("")
        dirs += idx.loc[idx.status == "done", "out_dir"].tolist()
    dirs = [d for d in dict.fromkeys(dirs) if os.path.exists(os.path.join(d, "arf_reproducibility_manifest.json"))]
    results, df = analyse(dirs, protocol)
    with open(a.out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    df.to_csv(os.path.splitext(a.out)[0] + "_windows.csv", index=False)
    for r in results["designs"]:
        d = r["design"]
        print(f"[{r['role']}] {d['active_roi']} / {d['fleet']} / {d.get('planet_swath_case')}"
              + (f" / without {', '.join(d['exclude_satellites'])}" if d.get("exclude_satellites") else "")
              + f": {len(r['windows_found'])} windows (missing {len(r['windows_missing'])}), alpha {r['alpha']:.4f}")
        for h in ("H1a", "H1b", "H2"):
            x = r[h]
            lcb = [v for k, v in x.items() if k.startswith("lcb_")][0]
            fmt = lambda v: "n/a" if v is None else f"{v:.2f}"
            print(f"    {h}: {x['status']} (evaluable windows {x['evaluable_windows']}, mean {fmt(x['mean'])}, "
                  f"lower bound {fmt(lcb)})")
    if results.get("runs_excluded_release_or_mixed"):
        print("Excluded (other release or mixed cells):", results["runs_excluded_release_or_mixed"])
    print("Written:", a.out)
    return results


if __name__ == "__main__":
    main()
