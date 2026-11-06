"""Figure for post hoc amendments A12 (what sets the floor) and A13 (PlanetScope
acquisition fraction), and Table A6. Inputs: the summary JSON and the two CSVs
written by arf/a12a13/summarise_a12_a13.py. Every plotted mean is re-derived from
the CSVs and asserted equal to the summary before anything is drawn.
Usage: fig_a12_a13.py SUMMARY_JSON FIG_OUT TABLE_OUT
"""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from figstyle import setup, save, OKABE, SEASON, FULL_W
import matplotlib.pyplot as plt

SUMM, FIG, TAB = sys.argv[1:4]
s = json.load(open(SUMM))
on = pd.read_csv(SUMM.replace(".json", "_a12_onsets.csv"))
win = pd.read_csv(SUMM.replace(".json", "_a13_windows.csv"))
FR = [1.0, 0.75, 0.5, 0.25]
# cross-checks against the summary
for season in ("dry", "wet"):
    g = on[on.season == season].groupby("window").floor_min.median() / 60
    assert abs(g.mean() - s["A12"][season]["floor_h_median"]["mean"]) < 1e-9
    for f in FR:
        w = win[(win.season == season) & (win.f == f)]
        for k in ("reach", "saving", "h1b_pct"):
            assert abs(w[k].mean() - s["A13"][season][f"{f:g}"][k]["mean"]) < 1e-9, (season, f, k)
OPEN_DATA_REACH_DRY = 53     # A7, open-data fleet only (dry season), as reported in Table 3 of the manuscript

setup()
fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.55), gridspec_kw={"width_ratios": [1.35, 1, 1]})
# (a) floor against onset local mean solar time
lo = min(s["A12"][k]["completion_lst_h_range"][0] for k in ("dry", "wet"))
hi = max(s["A12"][k]["completion_lst_h_range"][1] for k in ("dry", "wet"))
ax[0].axvspan(lo, hi, color=OKABE["grey"], alpha=0.18, lw=0)
for season, mk in (("dry", "o"), ("wet", "s")):
    g = on[on.season == season]
    ax[0].scatter(g.local_mean_solar_time_h, g.floor_min / 60, s=3, marker=mk, color=SEASON[season],
                  alpha=0.45, lw=0, label=f"{season.capitalize()} season")
ax[0].set_xlim(0, 24); ax[0].set_xticks(range(0, 25, 6))
ax[0].set_xlabel("Onset, local mean solar time (h)"); ax[0].set_ylabel("Floor (h)")
ax[0].text((lo + hi) / 2, 25.3, "mosaics\ncompleted", ha="center", va="top", fontsize=6.5)
ax[0].set_ylim(0, 25.5)
ax[0].set_title("(a) What sets the floor (A12)")
ax[0].legend(loc="lower left", frameon=False, markerscale=2.5, handletextpad=0.2)
# (b) reachable share and (c) image saving against the PlanetScope fraction
for season, mk, dx in (("dry", "o", -0.012), ("wet", "s", 0.012)):
    for i, key in ((1, "reach"), (2, "saving")):
        m = [s["A13"][season][f"{f:g}"][key]["mean"] for f in FR]
        lo_ = [s["A13"][season][f"{f:g}"][key]["ci95"][0] for f in FR]
        hi_ = [s["A13"][season][f"{f:g}"][key]["ci95"][1] for f in FR]
        x = np.array(FR) + dx
        ax[i].errorbar(x, m, yerr=[np.array(m) - lo_, np.array(hi_) - m], marker=mk, ms=3.5, lw=1, capsize=2,
                       color=SEASON[season], label=f"{season.capitalize()} season")
ax[1].axhline(OPEN_DATA_REACH_DRY, color=OKABE["black"], ls=":", lw=0.9)
ax[1].text(1.05, OPEN_DATA_REACH_DRY + 1.2, "open data only (dry)", fontsize=6.5, ha="left")
ax[1].set_ylim(30, 100); ax[1].set_ylabel("Region reachable in 24 h (%)")
ax[1].set_title("(b) Reach (A13)")
ax[2].axhline(0, color=OKABE["black"], lw=0.6)
ax[2].set_ylim(0, 10.5); ax[2].set_ylabel("Images saved by ECC vs FCFS")
ax[2].set_title("(c) Image saving (A13)")
for i in (1, 2):
    ax[i].set_xlim(1.07, 0.18); ax[i].set_xticks(FR); ax[i].set_xticklabels(["1", "0.75", "0.5", "0.25"])
    ax[i].set_xlabel("PlanetScope passes acquiring")
ax[2].legend(loc="upper right", frameon=False)
fig.tight_layout(w_pad=1.0)
print(FIG, save(fig, FIG))

# Table A6
def c(x, d=1):
    return f"{x['mean']:.{d}f} [{x['ci95'][0]:.{d}f}, {x['ci95'][1]:.{d}f}]"
L = [r"\begin{table}[H]", r"\centering", r"\footnotesize", r"\renewcommand{\arraystretch}{1.2}", r"\setlength{\tabcolsep}{3pt}",
     r"\caption{PlanetScope acquisition fraction (post hoc, A13): each PlanetScope pass acquires with probability $f$ "
     r"(20 seeds, nested fractions); Sentinel-2 and Landsat are unchanged. Means over the eight windows of each season "
     r"with window-block bootstrap 95\% intervals. Reach: median share of the region reachable within 24~h. "
     r"Region: share of onsets at which 95\% of the whole region is reachable. H1a and H1b: shares of onsets by the "
     r"protocol's rules. Saving: FCFS minus ECC images per onset. At $f = 1$ every value equals the run's own.}",
     r"\label{tab:a13}",
     r"\begin{tabular}{|c|c|c|c|c|c|c|}", r"\hline",
     r"Season & $f$ & Reach (\%) & Region (\%) & H1a (\%) & H1b (\%) & Saving (images) \\", r"\hline"]
for season in ("dry", "wet"):
    for f in FR:
        x = s["A13"][season][f"{f:g}"]
        L += [f"{season.capitalize()} & {f:g} & {c(x['reach'], 0)} & {c(x['region95'], 1)} & "
              f"{x['h1a_pct']['mean']:.0f} & {c(x['h1b_pct'], 1)} & {c(x['saving'], 1)} \\\\", r"\hline"]
L += [r"\end{tabular}", r"\end{table}"]
open(TAB, "w").write("\n".join(L) + "\n")
print(TAB)
