"""Result figures generated from the study's result files (no number typed by hand).

  fig_decisions.png   per-window values of H1a, H1b and H2 (primary design), mean, one-sided lower bound at
                      1 - alpha/3 (Student t over windows) and the threshold; the bound is recomputed here and
                      asserted equal to arf_confirmatory_results.json
  fig_images.png      (a) images per onset at the H1a level by selector, window medians (A7);
                      (b) data volume ECC vs FCFS per onset at the upper rates (A11)
  fig_budget.png      A9 verdicts by budget k and season (share of sampled onsets)
  fig_failsignal.png  A9: (a) ECDF of Delta T by season; (b) adaptive vs open-loop final coverage per realisation
  fig_cloud.png       A8: (a) cloudy share per Sentinel-2 overpass by season against p = 0.35 and EarthEnv;
                      (b) per-overpass cloudy share from the SCL reader against the catalogue's eo:cloud_cover

Usage: fig_results.py OUT_DIR   (input paths below; each figure checks its inputs against stored summaries)
"""
import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..")) if _o.path.basename(_o.path.dirname(_o.path.abspath(__file__))) in ("figscripts", "checks") else None; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
from paths import D as _D, A as _A, FACTS as _FACTS
import sys, json, os
import numpy as np, pandas as pd
from scipy import stats
import figstyle
import matplotlib.pyplot as plt

OUT = sys.argv[1]
D, A = _D, _A
facts = json.load(open(_FACTS))
figstyle.setup()
C = figstyle.OKABE
WLAB = ["3 Aug", "10 Aug", "17 Aug", "24 Aug", "31 Aug", "7 Sep", "14 Sep", "21 Sep"]


def lcb(x, alpha):
    x = np.asarray(x, float); n = len(x)
    return x.mean() - stats.t.ppf(1 - alpha, n - 1) * x.std(ddof=1) / np.sqrt(n)


# ---------------------------------------------------------------- decisions
conf = json.load(open(D + "arf_confirmatory_results.json"))
prim = [d for d in conf["designs"] if d["design"]["role"] == "primary"][0]
a4s = json.load(open(D + "a4_summary.json"))
h1b_min = [r for r in a4s["budget"] if r["k"] == 30][0]["h1b_min"]          # the one window below 100 (exact value)
h1b = [100.0 if x == 100.0 else h1b_min for x in facts["h1b_by_window"]]
assert all(abs(a - b) < 0.01 for a, b in zip(h1b, facts["h1b_by_window"]))
h2 = [c["h2_const"] for c in sorted(a4s["cloud"], key=lambda c: c["window"])]   # exact per-window H2 (A4 file)
assert all(abs(a - b) < 0.001 for a, b in zip(h2, facts["h2_by_window"]))
vals = {"H1a": [100.0] * 8, "H1b": h1b, "H2": h2}
assert prim["H1a"]["mean"] == 100.0 and prim["H1a"]["sd"] == 0.0     # every window 100%
fig, axs = plt.subplots(1, 3, figsize=(figstyle.FULL_W, 2.5))
spec = [("H1a", "(a) H1a: latency at the floor", "Onsets within margin (%)", 90, (55, 104)),
        ("H1b", "(b) H1b: fewer images than FCFS", "Eligible onsets (%)", 90, (55, 104)),
        ("H2", "(c) H2: value of the Fail-Signal", "Coverage gain (points)", 0, (-0.6, 5.6))]
for ax, (h, title, ylab, thr, ylim) in zip(axs, spec):
    x = np.asarray(vals[h], float)
    lb = lcb(x, 0.05 / 3) if x.std() > 0 else x.mean()
    lbj = [v for k, v in prim[h].items() if k.startswith("lcb")][0]
    assert abs(lb - lbj) < 1e-6 and abs(x.mean() - prim[h]["mean"]) < 1e-6, (h, lb, lbj)
    ax.scatter(range(8), x, s=16, color=C["blue"], zorder=3, label="Window")
    ax.axhline(x.mean(), color=C["black"], lw=0.9, label="Mean over windows")
    ax.axhline(lb, color=C["vermillion"], lw=0.9, ls="--", label=r"Lower bound ($1-\alpha/3$)")
    ax.axhline(thr, color=C["grey"], lw=0.9, ls=":", label="Threshold")
    ax.set_xticks(range(8)); ax.set_xticklabels(WLAB, rotation=60, fontsize=6.5)
    ax.set_ylim(*ylim); ax.set_title(title, fontsize=8); ax.set_ylabel(ylab)
    verdict = "supported" if lb >= thr and (h != "H2" or lb > thr) else "not supported"
    ax.text(0.98, 0.04, f"bound {lb:.1f} {'%' if h != 'H2' else 'pts'}: {verdict}" if h != "H2" else
            f"bound {lb:.2f} pts: {verdict}", transform=ax.transAxes, ha="right", fontsize=6.5)
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=4, fontsize=6.8, frameon=False, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=(0, 0.07, 1, 1))
figstyle.save(fig, f"{OUT}/fig_decisions.png"); print("fig_decisions ok")

# ---------------------------------------------------------------- images and volume
w = pd.read_csv(A + "a7results/a7_windows.csv")
w = w[w.arm == "Primary: CKGR+Okavango, all 96, conservative"].sort_values("window")
assert len(w) == 8
a7 = json.load(open(A + "a7results/a7_summary.json"))["arms"]["Primary: CKGR+Okavango, all 96, conservative"]
assert abs(w.ecc_images_median.mean() - a7["ecc_images_median"]["mean_over_windows"]) < 1e-9
v = pd.read_csv(A + "a11results/a11_onsets.csv")
assert len(v) == 584
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 2.7), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]
cols = [("ecc_images_median", "ECC"), ("coverage_first_images_median", "Coverage-first"),
        ("time_aware_images_median", "Time-aware"), ("fcfs_unbudgeted_images_median", "FCFS")]
for j, (c, lab) in enumerate(cols):
    ax.plot(np.arange(8) + (j - 1.5) * 0.12, w[c], "o", ms=3.5, color=figstyle.SELECTOR[lab], label=lab)
ax.set_xticks(range(8)); ax.set_xticklabels(WLAB, rotation=60, fontsize=6.5)
ax.set_ylabel("Images to 95% of reachable cells\n(median over onsets)")
ax.set_title("(a) Images per selector", fontsize=8); ax.legend(fontsize=6.3, frameon=False, ncol=2, loc="upper left")
ax.set_ylim(0, max(w[c].max() for c, _ in cols) * 1.35)
ax = axs[1]
for s, col in (("2026-09-14", C["purple"]), ):
    pass
ax.scatter(v.fcfs_mb_high / 1000, v.ecc_mb_high / 1000, s=4, alpha=0.5, color=C["blue"], edgecolors="none")
m = max(v.fcfs_mb_high.max(), v.ecc_mb_high.max()) / 1000 * 1.05
ax.plot([0, m], [0, m], color=C["grey"], lw=0.8, ls=":")
ax.text(m * 0.62, m * 0.68, "equal volume", rotation=38, fontsize=6.5, color="0.4")
ax.set_xlim(0, m); ax.set_ylim(0, m)
ax.set_xlabel("FCFS volume (GB)"); ax.set_ylabel("ECC volume (GB)")
ax.set_title("(b) Data volume per onset (upper rates)", fontsize=8)
n_below = int((v.saving_mb_high > 0).sum())
ax.text(0.97, 0.04, f"{n_below} of {len(v)} onsets below the line", transform=ax.transAxes, ha="right", fontsize=6.5)
fig.tight_layout()
figstyle.save(fig, f"{OUT}/fig_images.png"); print("fig_images ok")

# ---------------------------------------------------------------- budget verdicts (A9)
b = pd.read_csv(A + "a9results/a9_budget.csv")
assert len(b) == facts["a9|pairs"] and int(b.unknown.sum()) == facts["a9|unknown"]
order = [("ECC at floor", C["green"]), ("ECC at exact time", C["sky"]), ("ECC missed a cover", C["vermillion"]),
         ("no k-cover exists", "0.8")]
for s in ("dry", "wet"):
    for k in (5, 10, 15, 20):
        got = b[(b.season == s) & (b.k == k)].verdict.value_counts().to_dict()
        assert got == facts[f"a9|{s}|k{k}"], (s, k, got)
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 2.4), sharey=True)
for ax, s in zip(axs, ("dry", "wet")):
    bottom = np.zeros(4)
    for verdict, col in order:
        h = np.array([100 * ((b.season == s) & (b.k == k) & (b.verdict == verdict)).sum() / ((b.season == s) & (b.k == k)).sum()
                      for k in (5, 10, 15, 20)])
        ax.bar(range(4), h, bottom=bottom, color=col, width=0.65, label=verdict, edgecolor="white", lw=0.4)
        bottom += h
    assert np.allclose(bottom, 100)
    ax.set_xticks(range(4)); ax.set_xticklabels([f"k = {k}" for k in (5, 10, 15, 20)])
    ax.set_title(f"({'a' if s == 'dry' else 'b'}) {s.capitalize()} season", fontsize=8)
axs[0].set_ylabel("Sampled onsets (%)")
axs[1].legend(fontsize=6.3, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout()
figstyle.save(fig, f"{OUT}/fig_budget.png"); print("fig_budget ok")

# ---------------------------------------------------------------- Fail-Signal in hours (A9)
f = pd.read_csv(A + "a9results/a9_failsignal.csv")
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 2.6))
ax = axs[0]
for s in ("dry", "wet"):
    g = f[f.season == s]
    assert len(g) == facts[f"a9|{s}|fs"]["n"]
    worse = g.adaptive_reach_pct < g.openloop_reach_pct - 1e-9
    assert int(worse.sum()) == facts[f"a9|{s}|fs"]["never"]
    dt = np.where(worse, -np.inf, g.dT_L_h.to_numpy(float))
    x = np.sort(np.clip(dt, -3, 3))
    ax.step(x, np.arange(1, len(x) + 1) / len(x) * 100, where="post", color=figstyle.SEASON[s], lw=1.1,
            label=f"{s} ({int(worse.sum())} never reach it, at $-3$ h)")
ax.axvline(0, color=C["grey"], lw=0.7, ls=":")
ax.set_xlim(-3.1, 3.1); ax.set_xlabel(r"$\Delta T$: open-loop minus adaptive time (h)")
ax.set_ylabel("Realisations (%)"); ax.set_title("(a) Time gained by the Fail-Signal", fontsize=8)
ax.legend(fontsize=6.3, frameon=False, loc="upper left")
ax = axs[1]
for s in ("dry", "wet"):
    g = f[f.season == s]
    ax.scatter(g.openloop_reach_pct, g.adaptive_reach_pct, s=6, alpha=0.6, color=figstyle.SEASON[s], edgecolors="none", label=s)
lo = min(f.openloop_reach_pct.min(), f.adaptive_reach_pct.min()) - 1
ax.plot([lo, 100.5], [lo, 100.5], color=C["grey"], lw=0.8, ls=":")
ax.set_xlim(lo, 100.5); ax.set_ylim(lo, 100.5)
ax.set_xlabel("Open-loop final coverage (% reachable)"); ax.set_ylabel("Adaptive final coverage (% reachable)")
ax.set_title("(b) Coverage recovered", fontsize=8); ax.legend(fontsize=6.5, frameon=False, loc="upper left")
fig.tight_layout()
figstyle.save(fig, f"{OUT}/fig_failsignal.png"); print("fig_failsignal ok")

# ---------------------------------------------------------------- observed cloud (A8)
o = pd.read_csv(A + "a8results/a8_overpasses.csv")
for s in ("dry", "wet"):
    assert (o.season == s).sum() == facts[f"a8|{s}|overpasses"]
a4 = json.load(open(D + "a4_summary.json"))
ee = {c["window"][:7]: 100 * c["p_earthenv_mean"] for c in a4["cloud"]}
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 2.6))
ax = axs[0]
rng = np.random.default_rng(1)
for i, s in enumerate(("dry", "wet")):
    g = o[o.season == s]
    ax.scatter(i + rng.uniform(-0.15, 0.15, len(g)), g.cloudy_pct_of_observed, s=8 + 40 * g.tiles_observed / o.tiles_observed.max(),
               color=figstyle.SEASON[s], alpha=0.75, edgecolors="k", linewidths=0.3)
    pooled = 100 * g.tiles_cloudy.sum() / g.tiles_observed.sum()
    assert abs(pooled - facts[f"a8|{s}|cloudy_pooled"]) < 1e-6
    ax.plot([i - 0.28, i + 0.28], [pooled, pooled], color="k", lw=1.2)
ax.axhline(35, color=C["vermillion"], lw=0.9, ls="--", label="Model $p$ = 0.35")
ax.fill_between([-0.4, 0.4], min(ee.values()), max(ee.values()), color=C["green"], alpha=0.3, lw=0, label="EarthEnv, Aug-Sep (dry)")
ax.set_xticks([0, 1]); ax.set_xlim(-0.5, 1.5); ax.set_xticklabels(["Dry (Aug-Sep 2026)", "Wet (Jan-Feb 2026)"])
ax.set_ylabel("Cloudy cells per overpass (%)"); ax.set_ylim(-3, 103)
ax.set_title("(a) Observed cloud per overpass", fontsize=8); ax.legend(fontsize=6.3, frameon=False, loc="upper left")
ax = axs[1]
g = o.dropna(subset=["eo_cloud_cover_mean", "cloudy_pct_of_observed"])
for s in ("dry", "wet"):
    gg = g[g.season == s]
    ax.scatter(gg.eo_cloud_cover_mean, gg.cloudy_pct_of_observed, s=10, color=figstyle.SEASON[s], edgecolors="k", linewidths=0.3, label=s)
    r = np.corrcoef(gg.eo_cloud_cover_mean, gg.cloudy_pct_of_observed)[0, 1]
    assert abs(r - facts[f"a8|{s}|corr"]) < 1e-3, (s, r)
ax.plot([0, 100], [0, 100], color=C["grey"], lw=0.8, ls=":")
ax.set_xlim(-3, 103); ax.set_ylim(-3, 103)
ax.set_xlabel("Catalogue scene cloud cover (%)"); ax.set_ylabel("Cloudy cells, SCL reader (%)")
ax.set_title("(b) Reader against the catalogue", fontsize=8); ax.legend(fontsize=6.5, frameon=False, loc="upper left")
fig.tight_layout()
figstyle.save(fig, f"{OUT}/fig_cloud.png"); print("fig_cloud ok")
