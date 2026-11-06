"""Figure: where the fleet can and cannot see within 24 h of an onset (one primary window).

Inputs (one run folder of release 1.5.3): arf_refined_events.csv (Cell 5) and arf_onset_robustness.csv (Cell 7C).
For each of the 73 onsets (0..1440 min, 20-min step) the events with onset <= t_min < onset + 1440 are kept
(the rule of subset_by_time, Cell 7C); a cell is reachable when some kept event sees it with utility
u = max(0, 1 - A/theta_max) * D > 0 (Cell 6). Panel (a): share of onsets at which each cell is reachable.
Panel (b): median over onsets of the time from onset to the cell's first usable view (h), where reachable.
Check: the number of reachable cells per onset must equal the run's reachable_clear_tiles column exactly.

Usage: fig_reachability_map.py RUN_DIR OUT_PNG [--window "3 August 2026"]
"""
import sys, argparse
import numpy as np, pandas as pd, h3
import figstyle
from matplotlib.collections import PolyCollection
from matplotlib import colors, cm
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("run_dir"); ap.add_argument("out_png"); ap.add_argument("--window", default="3 August 2026")
a = ap.parse_args()
ev = pd.read_csv(f"{a.run_dir}/arf_refined_events.csv")
on = pd.read_csv(f"{a.run_dir}/arf_onset_robustness.csv")
ev["u"] = (1.0 - ev.A_nadir_deg / ev.theta_max_deg).clip(lower=0) * ev.D_downlink
ev = ev[ev.u > 0]
ev["deliv"] = ev.t_min + ev.T_wait_min.fillna(np.inf)
H, STEP, SHIFT = 1440, 20, 1440
onsets = np.arange(0, SHIFT + 1, STEP)
cells = sorted(set(pd.read_csv(f"{a.run_dir}/arf_mosaic_tiles.csv", usecols=["tile"]).tile)) \
    if __import__("os").path.exists(f"{a.run_dir}/arf_mosaic_tiles.csv") else None
reach = {}; first = {}
counts = []
for d in onsets:
    s = ev[(ev.t_min >= d) & (ev.t_min < d + H)]
    f = (s.deliv - d).groupby(s.tile).min()          # time from onset to first usable delivery (min)
    counts.append(len(f))
    for t, v in f.items():
        reach[t] = reach.get(t, 0) + 1
        first.setdefault(t, []).append(v / 60.0)
run_counts = on.set_index("onset_shift_min").loc[onsets, "reachable_clear_tiles"].to_numpy()
assert np.array_equal(np.array(counts), run_counts), "reachable counts differ from the run's onset results"
print(f"check passed: reachable cells per onset equal the run's ({min(counts)}-{max(counts)} cells)")

if cells is None:
    raise SystemExit("arf_mosaic_tiles.csv (the region's cell list) is needed to draw unreachable cells")
share = np.array([reach.get(t, 0) / len(onsets) for t in cells])
med = np.array([np.median(first[t]) if t in first else np.nan for t in cells])
polys = []
lat0 = np.mean([h3.cell_to_latlng(t)[0] for t in cells[::50]])
for t in cells:
    b = np.array(h3.cell_to_boundary(t))          # (lat, lng)
    polys.append(b[:, ::-1])
figstyle.setup()
fig, axs = plt.subplots(1, 2, figsize=(figstyle.FULL_W, 4.4), sharey=True)
panels = [(share * 100, "(a) Reachable within 24 h", cm.viridis, colors.Normalize(0, 100), "Share of the 73 onsets (%)"),
          (med, "(b) Time to first usable view", cm.magma_r, colors.Normalize(0, 24), "Median over onsets (h from onset)")]
for ax, (vals, title, cmap, norm, clab) in zip(axs, panels):
    fc = [cmap(norm(v)) if np.isfinite(v) and (cmap is not cm.viridis or v > 0) else (0.85, 0.85, 0.85, 1) for v in vals]
    ax.add_collection(PolyCollection(polys, facecolors=fc, edgecolors="none", linewidths=0, rasterized=True))
    ax.set_xlim(21.6, 25.6); ax.set_ylim(-23.45, -18.1)
    ax.set_aspect(1 / np.cos(np.radians(lat0)))
    ax.set_title(title, fontsize=8.5)
    ax.set_xticks([22, 23, 24, 25]); ax.set_xticklabels([f"{x}°E" for x in (22, 23, 24, 25)])
    ax.set_yticks([-23, -22, -21, -20, -19]); ax.set_yticklabels([f"{-y}°S" for y in (-23, -22, -21, -20, -19)])
    ax.grid(lw=0.3, ls=":", color="0.6")
    for sp in ("top", "right"): ax.spines[sp].set_visible(True)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal", pad=0.07, fraction=0.045)
    cb.ax.tick_params(labelsize=7); cb.set_label(clab, fontsize=7.5)
# scale bar (100 km) and north arrow on panel (a)
ax = axs[0]
km_per_deg_lon = 111.32 * np.cos(np.radians(-22.8))
x0, y0 = 23.85, -23.25
ax.plot([x0, x0 + 100 / km_per_deg_lon], [y0, y0], color="k", lw=2, solid_capstyle="butt")
ax.text(x0 + 50 / km_per_deg_lon, y0 + 0.09, "100 km", ha="center", fontsize=7)
ax.annotate("N", xy=(21.95, -18.45), xytext=(21.95, -18.95), ha="center", fontsize=8, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color="k", lw=0.8))
fig.text(0.01, -0.01, f"Window starting {a.window}; 73 onsets, 24-h horizon each; grey: never reachable (a) or not reachable (b). "
         "Cells: H3 resolution 7. Coordinates: WGS 84.", fontsize=6.5, color="0.35")
print("size", figstyle.save(fig, a.out_png))
print(f"cells never reachable: {(share == 0).sum()} of {len(cells)}; reachable at every onset: {(share == 1).sum()}")
print(f"median time to first view over cells (h): {np.nanmedian(med):.2f}")
