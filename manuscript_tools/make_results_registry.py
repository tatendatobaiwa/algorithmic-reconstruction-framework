"""Builds facts_results.json: every value the Results may quote, formatted, with its source file.
All inputs are result files of the runs (copied from Drive or produced by the post hoc scripts)."""
import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..")) if _o.path.basename(_o.path.dirname(_o.path.abspath(__file__))) in ("figscripts", "checks") else None; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
from paths import D as _D, A as _A, FACTS as _FACTS
import json, numpy as np, pandas as pd
D, A = _D, _A
R = {}
def put(val, src, *fmts):
    for f in fmts or ("{:.1f}",):
        R.setdefault(f.format(val), []).append(src)
conf = json.load(open(D + "arf_confirmatory_results.json"))
for d in conf["designs"]:
    tag = f"confirmatory {d['design'].get('active_roi')} {d['design'].get('fleet')} {d['design'].get('planet_swath_case')} {d['design'].get('exclude_satellites')}"
    for h in ("H1a", "H1b", "H2"):
        x = d[h]; put(x["evaluable_windows"], tag + f" {h} windows", "{:d}")
        put(x["mean"], tag + f" {h} mean", "{:.1f}", "{:.2f}", "{:.0f}")
        for k, v in x.items():
            if k.startswith("lcb") and v is not None:
                put(v, tag + f" {h} {k}", "{:.1f}", "{:.2f}", "{:.0f}")
        if x.get("sd") is not None: put(x["sd"], tag + f" {h} sd", "{:.1f}", "{:.2f}")
f = json.load(open(_FACTS))
for i, v in enumerate(f["h1b_by_window"]): put(v, f"h1b window {i}", "{:.1f}")
for i, v in enumerate(f["h2_by_window"]): put(v, f"h2 window {i}", "{:.2f}", "{:.1f}")
put(min(f["h2_by_window"]), "h2 min", "{:.2f}"); put(max(f["h2_by_window"]), "h2 max", "{:.2f}")
p1 = json.load(open(D + "p1_summary.json"))
for k, v in p1["tests"].items():
    if isinstance(v, dict):
        for kk in ("sign_p", "wilcoxon_p"): put(v[kk], f"p1 {k} {kk}", "{:.3f}", "{:.4f}")
for k in ("pred_overpasses", "matched"): put(f[f"archive|{k}"], f"archive {k}", "{:d}")
for v in f["archive|offset_range"]: put(v, "archive offset", "{:.2f}")
for v in f["archive|actual_tiles_predicted_range"]: put(v, "archive tiles predicted", "{:.1f}")
a7 = json.load(open(A + "a7results/a7_summary.json"))
for arm, v in a7["arms"].items():
    for k, x in v.items():
        if isinstance(x, dict) and x.get("mean_over_windows") is not None:
            put(x["mean_over_windows"], f"a7 {arm} {k} mean", "{:.1f}", "{:.0f}", "{:.2f}")
            for c in (x.get("ci95_window_bootstrap") or []):
                if c is not None: put(c, f"a7 {arm} {k} ci", "{:.1f}", "{:.0f}")
a4e = json.load(open(D + "a4_exact_economy_summary.json"))["windows"]
put(sum(w["onsets_solved"] for w in a4e), "a4 solved", "{:d}"); put(sum(w["ecc_optimal"] for w in a4e), "a4 optimal", "{:d}")
put(sum(w["onsets_tie_or_loss"] for w in a4e), "a4 ties", "{:d}"); put(int(max(w["ecc_extra_images_max"] for w in a4e)), "a4 excess", "{:d}")
for w in a4e: put(w["onsets_solved"], f"a4 {w['window']} solved", "{:d}"); put(w["ecc_optimal"], f"a4 {w['window']} optimal", "{:d}")
a4 = json.load(open(D + "a4_summary.json"))
for b in a4["budget"]:
    put(b["k"], "a4 k", "{:d}")
    for k in ("h1a_mean", "h1b_mean", "h1a_min", "ecc_images_median"): put(b[k], f"a4 k{b['k']} {k}", "{:.1f}", "{:.0f}")
    put(b["gap_max_min"], f"a4 k{b['k']} gap max", "{:.0f}")
for c in a4["cloud"]: put(100 * c["p_earthenv_mean"], f"a4 earthenv p {c['window']}", "{:.1f}")
for r in a4["ilp48"]:
    put(int(r["ecc_greedy_images"]), "a4 ilp48 ecc images", "{:d}"); put(int(r["ecc_exact_images"]), "a4 ilp48 exact", "{:d}")
a5 = json.load(open(D + "a5_summary.json"))
for var, v in a5["by_variant"].items():
    for h, x in v.items():
        for k in ("mean", "lcb_alpha", "lcb_alpha_over_3"): put(x[k], f"a5 {var} {h} {k}", "{:.1f}", "{:.2f}", "{:.3f}")
        put(x["above"], f"a5 {var} {h} above", "{:d}"); put(x["wilcoxon_p"], f"a5 {var} {h} wilcoxon", "{:.3f}", "{:.4f}")
x = a5["dry_season"]["earthenv_h2"]
for k in ("mean", "lcb_alpha", "lcb_alpha_over_3"): put(x[k], f"a5 dry earthenv h2 {k}", "{:.1f}", "{:.2f}", "{:.3f}")
for s in ("dry", "wet"):
    for k in ("overpasses", "full", "fail", "le5", "ge50", "scenes", "t95reached"): put(f[f"a8|{s}|{k}"], f"a8 {s} {k}", "{:d}")
    for k in ("cloudy_pooled", "obs48", "clear48"): put(f[f"a8|{s}|{k}"], f"a8 {s} {k}", "{:.1f}")
    put(f[f"a8|{s}|corr"], f"a8 {s} corr", "{:.2f}", "{:.3f}")
for k in ("dry|08", "dry|09", "wet|01", "wet|02"): put(f[f"a8|{k}|cloudy"], f"a8 {k}", "{:.1f}")
put(f["a9|pairs"], "a9 pairs", "{:d}"); put(f["a9|unknown"], "a9 unknown", "{:d}")
for s in ("dry", "wet"):
    for k in (5, 10, 15, 20):
        for verdict, n in f[f"a9|{s}|k{k}"].items(): put(n, f"a9 {s} k{k} {verdict}", "{:d}")
    fs = f[f"a9|{s}|fs"]
    for k in ("n", "never"): put(fs[k], f"a9 {s} fs {k}", "{:d}")
    for k in ("sooner", "same", "later"): put(fs[k], f"a9 {s} fs {k}", "{:.1f}", "{:.0f}")
    put(fs["gain_mean"], f"a9 {s} fs gain", "{:.1f}", "{:.2f}")
    for q in fs["q"]: put(q, f"a9 {s} fs quartile", "{:.2f}")
fsd = pd.read_csv(A + "a9results/a9_failsignal.csv")
for s, g in fsd.groupby("season"):
    worse = int(((g.adaptive_reach_pct < g.openloop_reach_pct - 1e-9)).sum()); put(worse, f"a9 {s} adaptive ends lower", "{:d}")
    put(100 * worse / len(g), f"a9 {s} adaptive ends lower pct", "{:.0f}", "{:.1f}")
    put(float((g.openloop_t50_h - g.adaptive_t50_h).median()), f"a9 {s} t50 diff median", "{:.0f}", "{:.1f}")
    put(float((g.openloop_t80_h - g.adaptive_t80_h).median()), f"a9 {s} t80 diff median", "{:.0f}", "{:.1f}")
s11 = json.load(open(A + "a11results/a11_manuscript_stats.json"))
for v in s11["saving_gb_median_per_onset"] + s11["reduction_pct_median_per_onset"] + s11["reduction_pct_iqr_low_rate"] + s11["reduction_pct_range"]:
    put(v, "a11", "{:.1f}", "{:.0f}")
for k in ("onsets_saving_positive", "onsets_saving_zero", "onsets_ecc_more_data", "onsets"): put(s11[k], f"a11 {k}", "{:d}")
put(100 * s11["planetscope_share_of_saving_high_pooled"], "a11 ps share pooled", "{:.1f}")
put(100 * s11["planetscope_share_min_window"], "a11 ps share min", "{:.1f}")
for k in ("window_mean_of_median_saving_gb", "window_mean_of_median_reduction_pct"):
    for lab in ("low", "high"):
        put(s11[k][lab]["mean"], f"a11 {k} {lab}", "{:.1f}")
        for c in s11[k][lab]["ci95"]: put(c, f"a11 {k} {lab} ci", "{:.1f}")
W = pd.read_csv(A + "a11results/a11_windows_recomputed.csv")
for r in W.itertuples():
    for c in ("ecc_img", "fcfs_img"): put(getattr(r, c), f"a11 {r.window} {c}", "{:.0f}")
    for c in ("saving_low", "saving_high", "red_low", "red_high", "ecc_low", "ecc_high", "fcfs_low", "fcfs_high"): put(getattr(r, c), f"a11 {r.window} {c}", "{:.1f}")
json.dump(R, open("facts_results.json", "w"), indent=0)
print(len(R), "formatted values registered")
