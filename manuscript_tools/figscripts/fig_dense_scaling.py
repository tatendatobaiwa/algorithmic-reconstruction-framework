"""Paper 2 preliminary figure: physical floor and images versus fleet size (A3 dense-constellation arm, release 1.5.4),
against the real 96-satellite fleet (release 1.5.3). Inputs: dense_runs.csv, eo_vs_dense_floor.csv (dense_summary on Drive).
Checks: window means equal dense_spread_windows.csv; real floors equal the confirmatory runs' floors in eo_vs_dense_floor.csv.
Usage: fig_dense_scaling.py DATA_DIR OUT_PNG"""
import sys, numpy as np, pandas as pd
import figstyle
import matplotlib.pyplot as plt
D, OUT = sys.argv[1], sys.argv[2]
r = pd.read_csv(f"{D}/dense_runs.csv"); s = pd.read_csv(f"{D}/dense_spread_windows.csv"); e = pd.read_csv(f"{D}/eo_vs_dense_floor.csv")
for N in (100, 500, 1000):
    assert abs(r[r.N == N].floor_t95r_median_min.mean() - s[(s.N == N) & (s.quantity == "floor_t95r_median_min")]["mean"].iloc[0]) < 1e-9
figstyle.setup(); C = figstyle.OKABE
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 2.7))
ax = axs[0]
for N, g in r.groupby("N"):
    ax.scatter(np.full(len(g), N) * np.exp(np.random.default_rng(N).uniform(-0.06, 0.06, len(g))), g.floor_t95r_median_min / 60,
               s=10, color=C["blue"], alpha=0.7, edgecolors="none")
m = r.groupby("N").floor_t95r_median_min.median() / 60
ax.plot(m.index, m.values, color=C["blue"], lw=1, label="Starlink geometry (median of windows)")
ax.scatter(np.full(len(e), 96), e.real_floor / 60, s=12, marker="s", color=C["vermillion"], label="Real EO fleet, 96 satellites", zorder=3)
ax.set_xscale("log"); ax.set_xlabel("Satellites in the fleet"); ax.set_ylabel("Floor to 95% of reachable cells (h)")
ax.set_title("(a) Physical floor against fleet size", fontsize=8); ax.legend(fontsize=6.3, frameon=False)
ax = axs[1]
mi = r.groupby("N")[["h1b_median_images_ecc", "h1b_median_images_fcfs"]].median()
ax.plot(mi.index, mi.h1b_median_images_ecc, "o-", color=figstyle.SELECTOR["ECC"], ms=3.5, label="ECC")
ax.plot(mi.index, mi.h1b_median_images_fcfs, "s-", color=figstyle.SELECTOR["FCFS"], ms=3.5, label="FCFS")
ax.set_xscale("log"); ax.set_xlabel("Satellites in the fleet"); ax.set_ylabel("Images to 95% of reachable cells\n(median over windows)")
ax.set_title("(b) Images against fleet size", fontsize=8); ax.legend(fontsize=6.5, frameon=False); ax.set_ylim(0, None)
fig.text(0.01, -0.02, "Exploratory (amendment A3): imaging parameters of the hypothetical constellation are assumptions. N = 2000 and 4000: window of 3 August 2026 only.",
         fontsize=6.3, color="0.35")
fig.tight_layout(); figstyle.save(fig, OUT)
print(m.round(2).to_dict()); print(mi.to_dict())
