# =============================================================================
# CELL 7: MANUSCRIPT GRAPH GENERATOR
# =============================================================================
#
# This module generates publication-quality figures for the LaTeX manuscript
# from the simulation results produced by Cell 6. It produces four key
# visualisations, drawn at IEEE column widths (v1.5.2):
#   1. Submodular Coverage Progression (Scenario A vs. B)
#   2. Mean Observation Quality Progression
#   3. Latency vs. Budget, in hours, against the physical floor and the H1a
#      margin (v1.5.2; was the falsified 10-minute H1 target, in minutes)
#   4. Coverage vs. Budget (k) Trade-off
#   7. Mosaic map: which image supplied each tile, and when
#
# Dependencies:
#   - Cell 0: ARFConfig (provides ACTIVE_N, etc.)
#   - Cell 1: target_tiles
#   - Cell 5: refined_events_df
#   - Cell 6: hist_clear, hist_cloud, sweep_results, latency metrics
#
# Outputs:
#   - arf_manuscript_figures.png (Composite raster preview)
#   - arf_manuscript_figures.pdf (Vector output for LaTeX)
#   - Individual figure files (PNG/PDF per subplot)
# =============================================================================
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
#
# References:
#   [9] Hunter, J.D. (2007). "Matplotlib: A 2D Graphics Environment."
#       Computing in Science & Engineering, 9(3), 90-95.
#   [21] Nemhauser, G.L., Wolsey, L.A. & Fisher, M.L. (1978). "An analysis
#        of approximations for maximizing submodular set functions—I."
#        Mathematical Programming, 14, 265–294.
# =============================================================================

"""
ARF Manuscript Graph Generator
================================

This module generates publication-quality figures for the research
manuscript: Times New Roman (or a metric-compatible clone), embedded
TrueType fonts, 300 dpi PNG plus vector PDF, drawn at IEEE column widths
(3.5 in single, 7.16 in double) with 8-9 pt text (v1.5.2).

Figure Descriptions:
    1. **Coverage Progression**: Shows how cumulative spatial coverage
       increases with each greedy selection step, comparing clear-sky
       (Scenario A) and cloud-obscured (Scenario B) conditions. The
       geometric ceiling is annotated.

    2. **Quality Progression**: Shows how the mean observation quality
       of the covered tiles evolves with each selection step.

    3. **Latency vs. Budget**: clear-sky completion time and time to the
       H1a level (95% of reachable tiles) for each budget k, in hours,
       against the physical floor and the floor plus the H1a margin.

    4. **Coverage vs. Budget Trade-off**: Shows the coverage achieved
       for different cardinality constraints (k), justifying the
       choice of MAX_SELECTION_STEPS.

Part of ARF (the cell's release is _ARF_CELL_VERSION below).
"""

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 7 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["7"] = _ARF_CELL_VERSION

import os
import sys
import logging
from typing import Dict, List, Tuple, Optional, Any

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logger = logging.getLogger("arf.cell7")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# =============================================================================
# CUSTOM EXCEPTIONS
# =============================================================================

class FigureGenerationError(Exception):
    """Raised when figure generation fails."""
    pass


# =============================================================================
# MATPLOTLIB CONFIGURATION
# =============================================================================
# Publication-quality settings for LaTeX manuscript integration.

plt.rcParams.update({
    # Font configuration (same fallback list as Cell 3: Times New Roman, else a
    # metric-compatible Times clone such as Liberation Serif on Colab)
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Liberation Serif', 'Nimbus Roman', 'FreeSerif', 'DejaVu Serif'],
    'font.size': 10,
    'mathtext.fontset': 'stix',

    # Axis configuration
    'axes.linewidth': 1.0,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',

    # Tick configuration
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,

    # Legend configuration
    'legend.fontsize': 8,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',

    # Figure configuration
    'figure.dpi': 100,
    'figure.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    # Embed fonts in PDFs as TrueType (Type 42) rather than matplotlib's default Type 3
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# =============================================================================
# COLOR PALETTE
# =============================================================================
# Consistent color scheme across all figures.

COLORS = {
    'clear_sky': '#1b9e77',       # Teal green
    'cloud_cover': '#d95f02',     # Orange
    'ceiling': '#7570b3',         # Purple
    'guarantee': '#e7298a',       # Pink (unused since the (1-1/e) line was removed)
    'level': '#e7298a',           # Pink: time to the H1a level (v1.5.2)
    'h1_target': '#66a61e',       # Green: floor + H1a margin (v1.5.2; was the 10-min target)
    'grid': '#cccccc',            # Light gray
}


# =============================================================================
# DATA EXTRACTION
# =============================================================================

# Verify required globals exist
_required_globals = {
    'hist_clear': 'Cell 6',
    'hist_cloud': 'Cell 6',
    'refined_events_df': 'Cell 5',
    'target_tiles': 'Cell 1',
}

_missing = []
for var, cell in _required_globals.items():
    if var not in globals():
        _missing.append(f"{var} (from {cell})")

if _missing:
    raise FigureGenerationError(
        f"Missing required globals: {', '.join(_missing)}. "
        f"Please run the required cells first."
    )

# --- Extract Real Ceiling Line ---
# Derived from actual access-event geometry, not the greedy algorithm: tiles
# with at least one USABLE clear-sky observation (utility > 0), the same rule
# every coverage figure uses (v1.5.2; was any refined row, which also counted
# rows of zero utility -- a pass with no downlink -- identical in all runs so far).
max_achievable_tiles = (int(np.unique(event_data.tile_idx[event_data.util > 0.0]).size)
                        if 'event_data' in globals() else refined_events_df['tile'].nunique())
total_target_tiles = len(target_tiles)
max_achievable_pct = 100 * max_achievable_tiles / total_target_tiles

logger.info(f"Max achievable coverage (geometric ceiling): "
            f"{max_achievable_tiles} / {total_target_tiles} tiles = "
            f"{max_achievable_pct:.2f}%")

# --- Extract Selection History Metrics ---
steps_clear = [h['step'] for h in hist_clear]
cov_clear_pct = [h['coverage_pct'] for h in hist_clear]
qual_clear = [h['mean_quality'] for h in hist_clear]

steps_cloud = [h['step'] for h in hist_cloud]
cov_cloud_pct = [h['coverage_pct'] for h in hist_cloud]
qual_cloud = [h['mean_quality'] for h in hist_cloud]

# --- Extract Sweep Results (if available) ---
sweep_data = None
if 'sweep_results' in globals() and len(sweep_results) > 0:
    sweep_data = sweep_results
    logger.info(f"Coverage sweep data available: {len(sweep_data)} k-values")

# --- Extract Latency Data ---
latency_clear = scenario_results.get(_scenario_key(), {}).get('latency_clear_min')
latency_cloud = scenario_results.get(_scenario_key(), {}).get('latency_cloud_min')

# --- H1a reference for the latency figure (v1.5.2) ---------------------------
# The original H1 target (mosaic complete within 10 min) was falsified by the
# real-fleet runs and is kept only in the reports. The latency figure instead
# shows the quantities H1a tests: the time at which the k-image selection's
# delivered mosaic reaches H1_REACHABLE_COVERAGE_PCT of the reachable tiles,
# the physical floor for that level (every pass in the window, floor_time of
# Cell 6) and the floor plus the H1a margin (one orbital period at the fleet's
# median altitude unless set in Cell 0 -- the same rule as Cell 7C). Full
# window, clear sky, coverage-first greedy (Cell 6's sweep); Cell 7C tests H1a
# over shifted onsets with the earliest-completion cover.
_H1L7 = float(globals().get("H1_REACHABLE_COVERAGE_PCT", 95.0))
_reach7 = int(np.unique(event_data.tile_idx[event_data.util > 0.0]).size) if len(event_data) else 0
_need7 = int(np.ceil(_H1L7 / 100.0 * _reach7 - 1e-9)) if _reach7 else 0
_floor7 = floor_time(event_data, _H1L7) if ('floor_time' in globals() and _reach7) else None
_margin7 = (float(H1_FLOOR_MARGIN_MIN_CONFIG) if globals().get("H1_FLOOR_MARGIN_MIN_CONFIG") is not None else
            float(2.0 * np.pi * np.sqrt((Re_km + CONSTELLATION_ALTITUDE_KM) ** 3 / mu_earth) / 60.0))


def time_to_level(selected, store, need):
    """Delivery time (min) at which the images in `selected`, in delivery order
    (ties in selection order), first cover `need` tiles with a usable clear-sky
    observation (utility > 0); None if they never do."""
    if need <= 0 or not selected:
        return None
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n = 0
    for key in sorted(selected, key=delivery_sort_key(store)):
        w = store.wait_lookup.get(key)
        if w is None or np.isnan(w):
            break                                          # undelivered images sort last
        s = store.rows(store.index[key])
        t = store.tile_idx[s][store.util[s] > 0.0]
        new = t[~seen[t]]
        if new.size:
            seen[new] = True
            n += new.size
            if n >= need:
                return float(key[1] + w)
    return None


_t_level_by_k = {}
if sweep_data:
    _k_all = [d['k'] for d in sweep_data]
    _hist7 = run_submodular_greedy_selection(event_data, max(_k_all), cloud_mask=None, gap_fill=False)[3]
    _sel7 = [h['selected'] for h in _hist7]          # greedy is incremental: budget k = the first k
    _t_level_by_k = {k: time_to_level(_sel7[:k], event_data, _need7) for k in _k_all}


# =============================================================================
# FIGURE STYLE (v1.5.2): drawn at the printed size
# =============================================================================
# IEEE Author Center, "Resolution and Size": one column 3.5 in (88.9 mm), two
# columns 7.16 in (182 mm), colour and greyscale above 300 dpi; "File
# Formatting": type "approximately 9-10 point when viewed at full size", fonts
# embedded. Earlier versions drew 7 x 5 in single figures (and a 14 x 10 in
# composite) with 10 pt text, which prints at about 5 pt once fitted to a
# column. Each figure is now drawn at its column width with 8-9 pt text.
IEEE_ONE_COL_IN, IEEE_TWO_COL_IN = 3.5, 7.16
_IEEE_RC = {"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "axes.titleweight": "bold",
            "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "legend.framealpha": 0.9,
            "lines.linewidth": 1.2, "lines.markersize": 3.5, "axes.linewidth": 0.6,
            "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42}


def _save_figure(fig, prefix):
    for fmt in ('png', 'pdf'):
        filepath = os.path.join(OUT_DIR if 'OUT_DIR' in globals() else os.getcwd(), f"{prefix}.{fmt}")
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        logger.info(f"Saved: {filepath}")


def _grid(ax):
    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5, color=COLORS['grid'])
    # v1.5.3: every figure that calls this counts images (k) along x: whole-number ticks
    # (a run with three candidate images showed k = 1.0, 1.5, 2.0, ...)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))


def _no_sweep_note(ax):
    """v1.5.3: say why a budget panel (Figures 3 and 4) is empty: the sweep
    starts at k = 5, and a window with fewer candidate images has no budget to
    sweep (open data alone, CKGR + Okavango, 15 October 2026: 3 images)."""
    _ed = globals().get("event_data")
    n = len(_ed) if _ed is not None else 0
    ax.text(0.5, 0.5, f"No budget sweep: the window offers {n} candidate image{'' if n == 1 else 's'},\n"
                      f"fewer than the smallest budget swept (k = 5)", transform=ax.transAxes,
            ha="center", va="center", fontsize=7.5, color="0.3")


# =============================================================================
# FIGURE 1: COVERAGE PROGRESSION
# =============================================================================
# The (1-1/e) guarantee is NOT drawn: the NWF (1978) bound compares the greedy
# with OPT_k, the best coverage achievable with exactly k selections, not with
# the geometric ceiling drawn here; OPT_k is not efficiently computable, so no
# fixed line on this axis represents the guarantee (Cell 7C measures the gap
# exactly with an ILP on sub-problems).

def _draw_coverage_progression(ax, steps_clear, cov_clear_pct, steps_cloud, cov_cloud_pct, max_achievable_pct,
                               cloud_label='Scenario B (Cloud / Fail-Signal)'):
    ax.plot(steps_clear, cov_clear_pct, color=COLORS['clear_sky'], linestyle='-', marker='o',
            label='Scenario A (Clear Sky)', zorder=3)
    ax.plot(steps_cloud, cov_cloud_pct, color=COLORS['cloud_cover'], linestyle='--', marker='s',
            label=cloud_label, zorder=3)
    ax.axhline(y=max_achievable_pct, color=COLORS['ceiling'], linestyle=':', linewidth=1.0, alpha=0.8,
               label=f'Geometric Ceiling ({max_achievable_pct:.1f}%)', zorder=2)
    ax.set_xlabel('Greedy Selection Step (k)')
    ax.set_ylabel('Cumulative Spatial Coverage (%)')
    ax.legend(loc='lower right')
    _grid(ax)
    ax.set_ylim(bottom=0, top=min(100, max_achievable_pct * 1.1))


def plot_coverage_progression(steps_clear, cov_clear_pct, steps_cloud, cov_cloud_pct, max_achievable_pct,
                              save_prefix: str = "arf_fig1"):
    """Cumulative coverage after each greedy step, clear sky and cloud, with
    the geometric ceiling. One IEEE column wide."""
    with plt.rc_context(_IEEE_RC):
        fig, ax = plt.subplots(figsize=(IEEE_ONE_COL_IN, 2.7))
        _draw_coverage_progression(ax, steps_clear, cov_clear_pct, steps_cloud, cov_cloud_pct, max_achievable_pct)
        ax.set_title('Submodular Coverage Progression')
        fig.tight_layout()
        _save_figure(fig, save_prefix)
    return fig, ax


# =============================================================================
# FIGURE 2: QUALITY PROGRESSION
# =============================================================================

def _draw_quality_progression(ax, steps_clear, qual_clear, steps_cloud, qual_cloud,
                              cloud_label='Scenario B (Cloud / Fail-Signal)'):
    ax.plot(steps_clear, qual_clear, color=COLORS['clear_sky'], linestyle='-', marker='o',
            label='Scenario A (Clear Sky)', zorder=3)
    ax.plot(steps_cloud, qual_cloud, color=COLORS['cloud_cover'], linestyle='--', marker='s',
            label=cloud_label, zorder=3)
    ax.set_xlabel('Greedy Selection Step (k)')
    ax.set_ylabel('Mean Observation Quality (Q)')
    ax.legend(loc='lower right')
    _grid(ax)
    ax.set_ylim(bottom=0, top=1.05)


def plot_quality_progression(steps_clear, qual_clear, steps_cloud, qual_cloud, save_prefix: str = "arf_fig2"):
    """Mean utility of the covered tiles after each greedy step. One column wide."""
    with plt.rc_context(_IEEE_RC):
        fig, ax = plt.subplots(figsize=(IEEE_ONE_COL_IN, 2.7))
        _draw_quality_progression(ax, steps_clear, qual_clear, steps_cloud, qual_cloud)
        ax.set_title('Observation Quality Progression')
        fig.tight_layout()
        _save_figure(fig, save_prefix)
    return fig, ax


# =============================================================================
# FIGURE 3: LATENCY AGAINST THE BUDGET (H1a reference, hours)
# =============================================================================

def _draw_latency(ax, sweep_data, t_level_by_k, floor_min, margin_min, level_pct):
    """Curves and floor lines of Figure 3 (and composite panel (c)). Returns
    the legend handles and labels; the caller puts the legend under the axes
    (v1.5.3: inside them it covered the curves).
    v1.5.3: the y-axis spans the plotted values, the floor and the floor plus
    the margin, instead of starting at 0 h: in a window whose first passes
    come a day in (15 October 2026: every value between 24 and 26 h) the
    curves lay pressed against the top of the axes. A budget with no
    completion time is left out instead of being drawn at 0 h."""
    ys = []
    if sweep_data:
        k_vals = [d['k'] for d in sweep_data]
        y_done = [np.nan if d.get('latency_clear_min') is None else d['latency_clear_min'] / 60.0
                  for d in sweep_data]
        ax.plot(k_vals, y_done, color=COLORS['clear_sky'], linestyle='-', marker='o',
                label='Mosaic complete (own final coverage)', zorder=3)
        ys += y_done
        kk = [k for k in k_vals if t_level_by_k.get(k) is not None]
        if kk:
            y_lvl = [t_level_by_k[k] / 60.0 for k in kk]
            ax.plot(kk, y_lvl, color=COLORS['level'], linestyle='-', marker='^',
                    label=f'{level_pct:g}% of reachable tiles', zorder=3)
            ys += y_lvl
    if floor_min is not None:
        # v1.5.3: the floor and the margin are drawn over the curves (zorder 4,
        # was 2): the curve of the H1a level often lies on the floor.
        ax.axhline(floor_min / 60.0, color='k', linestyle=':', linewidth=1.0,
                   label=f'Physical floor, {level_pct:g}% (all passes)', zorder=4)
        ax.axhline((floor_min + margin_min) / 60.0, color=COLORS['h1_target'], linestyle='--', linewidth=1.0,
                   label=f'Floor + {margin_min:.0f} min (H1a margin)', zorder=4)
        ys += [floor_min / 60.0, (floor_min + margin_min) / 60.0]
    ax.set_xlabel('Selection Budget (k)')
    ax.set_ylabel('Time from analysis start (h)')
    _grid(ax)
    if not sweep_data:
        _no_sweep_note(ax)
    ys = [float(y) for y in ys if y is not None and np.isfinite(y)]
    if ys:
        lo, hi = min(ys), max(ys)
        pad = max(0.08 * (hi - lo), 0.1)                # at least 6 min above and below the values
        ax.set_ylim(max(lo - pad, 0.0), hi + pad)
    return ax.get_legend_handles_labels()


def plot_latency_analysis(sweep_data, t_level_by_k, floor_min, margin_min, level_pct, save_prefix: str = "arf_fig3"):
    """
    Clear-sky latency of the coverage-first greedy against the budget k, in
    hours: completion time (each budget's own final coverage) and the time to
    the H1a level of the reachable tiles, with the physical floor for that
    level and the floor plus the H1a margin. (Until v1.5.1 this figure drew
    the falsified 10-minute target, in minutes.) One column wide; from v1.5.3
    the y-axis spans the data and the legend sits under the axes.
    """
    with plt.rc_context(_IEEE_RC):
        fig, ax = plt.subplots(figsize=(IEEE_ONE_COL_IN, 3.45))
        handles, labels = _draw_latency(ax, sweep_data, t_level_by_k, floor_min, margin_min, level_pct)
        ax.set_title('Latency vs. Budget (clear sky)')
        fig.tight_layout(rect=(0, 0.22 if handles else 0, 1, 1))
        if handles:
            fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=1, frameon=False,
                       fontsize=7.5)
        _save_figure(fig, save_prefix)
    return fig, ax


# =============================================================================
# FIGURE 4: COVERAGE VS. BUDGET TRADE-OFF
# =============================================================================

def _draw_budget(ax, sweep_data, max_achievable_pct, cloud_label='Scenario B (Cloud / Fail-Signal)'):
    if sweep_data:
        k_vals = [d['k'] for d in sweep_data]
        ax.plot(k_vals, [d['coverage_clear_pct'] for d in sweep_data], color=COLORS['clear_sky'], linestyle='-',
                marker='o', label='Scenario A (Clear Sky)', zorder=3)
        ax.plot(k_vals, [d['coverage_cloud_pct'] for d in sweep_data], color=COLORS['cloud_cover'], linestyle='--',
                marker='s', label=cloud_label, zorder=3)
    ax.axhline(y=max_achievable_pct, color=COLORS['ceiling'], linestyle=':', linewidth=1.0, alpha=0.8,
               label=f'Geometric Ceiling ({max_achievable_pct:.1f}%)', zorder=2)
    if not sweep_data:
        _no_sweep_note(ax)
    ax.set_xlabel('Selection Budget (k)')
    ax.set_ylabel('Coverage (%)')
    ax.legend(loc='lower right')
    _grid(ax)
    ax.set_ylim(bottom=0, top=min(100, max_achievable_pct * 1.1))


def plot_coverage_budget_tradeoff(sweep_data, max_achievable_pct, save_prefix: str = "arf_fig4"):
    """Coverage at each budget k (no gap-fill), clear sky and cloud, with the
    geometric ceiling (no (1-1/e) line: see Figure 1). One column wide."""
    with plt.rc_context(_IEEE_RC):
        fig, ax = plt.subplots(figsize=(IEEE_ONE_COL_IN, 2.7))
        _draw_budget(ax, sweep_data, max_achievable_pct)
        ax.set_title('Coverage vs. Budget Trade-off')
        fig.tight_layout()
        _save_figure(fig, save_prefix)
    return fig, ax


# =============================================================================
# GENERATE ALL FIGURES
# =============================================================================

logger.info("Generating manuscript figures...")
fig1, ax1 = plot_coverage_progression(steps_clear, cov_clear_pct, steps_cloud, cov_cloud_pct, max_achievable_pct,
                                      save_prefix="arf_fig1_coverage_progression")
fig2, ax2 = plot_quality_progression(steps_clear, qual_clear, steps_cloud, qual_cloud,
                                     save_prefix="arf_fig2_quality_progression")
fig3, ax3 = plot_latency_analysis(sweep_data, _t_level_by_k, _floor7, _margin7, _H1L7,
                                  save_prefix="arf_fig3_latency_analysis")
fig4, ax4 = plot_coverage_budget_tradeoff(sweep_data, max_achievable_pct, save_prefix="arf_fig4_coverage_budget")

plt.show()
# v1.5.1: release the four single figures once shown (batches of runs otherwise
# accumulate open figures: matplotlib warns beyond 20)
for _f in [globals().get(n) for n in ("fig1", "fig2", "fig3", "fig4")]:
    if _f is not None:
        plt.close(_f)


# =============================================================================
# COMPOSITE FIGURE (for manuscript): two IEEE columns wide
# =============================================================================

logger.info("Generating composite manuscript figure...")
with plt.rc_context(_IEEE_RC):
    fig_composite, axes = plt.subplots(2, 2, figsize=(IEEE_TWO_COL_IN, 6.0))
    _draw_coverage_progression(axes[0, 0], steps_clear, cov_clear_pct, steps_cloud, cov_cloud_pct,
                               max_achievable_pct, cloud_label='Scenario B (Cloud)')
    axes[0, 0].set_title('(a) Submodular Coverage Progression', loc='left')
    _draw_quality_progression(axes[0, 1], steps_clear, qual_clear, steps_cloud, qual_cloud,
                              cloud_label='Scenario B (Cloud)')
    axes[0, 1].set_title('(b) Observation Quality Progression', loc='left')
    _h3c, _l3c = _draw_latency(axes[1, 0], sweep_data, _t_level_by_k, _floor7, _margin7, _H1L7)
    axes[1, 0].set_title('(c) Latency vs. Budget (clear sky)', loc='left')
    _draw_budget(axes[1, 1], sweep_data, max_achievable_pct, cloud_label='Scenario B (Cloud)')
    axes[1, 1].set_title('(d) Coverage vs. Budget Trade-off', loc='left')
    # v1.5.3: panel (c)'s legend under the panels, from panel (c)'s left edge
    fig_composite.tight_layout(rect=(0, 0.07 if _h3c else 0, 1, 1))
    if _h3c:
        fig_composite.legend(_h3c, _l3c, loc='lower left', ncol=2, frameon=False, fontsize=7.5,
                             bbox_to_anchor=(axes[1, 0].get_position().x0, 0.0))
    _save_figure(fig_composite, "arf_manuscript_figures")

plt.show()
plt.close(fig_composite)


# =============================================================================
# FIGURE 7: MOSAIC MAP (which image supplied each tile, and when)
# =============================================================================
# The mosaic the selector builds from the Cell 6 selections, drawn on the
# ROI's H3 cells. Each covered tile is coloured by the delivery time (pass
# time + downlink wait, hours from the analysis start) of the FIRST selected
# image that covers it under the selection's own rule (utility > 0, not
# cloud-obscured): the mosaic as it stands when it is complete -- Cell 6's
# H1 completion time is the latest of these times; later images can only
# improve quality. Thin dark lines separate tiles supplied by different images
# (seamlines). Grey: no image in the window could cover the tile. Hatched: an
# image in the window could have covered it, but none of the selected images
# did. (a) clear sky; (b) cloud, with the configured cloud knowledge.
# Uses the map helpers of Cell 3. Also writes arf_mosaic_tiles.csv (one row per
# tile: delivery time and source image in each panel), which loads into any GIS
# through the H3 cell IDs.

MOSAIC_MAP_POLICY = "time_aware"     # "time_aware" (the selector H1 tests) or "coverage_first"
MOSAIC_COLORMAP = "viridis"          # perceptually uniform and colour-vision-deficiency safe
MOSAIC_NONE_COLOR = "#bdbdbd"        # tiles with no usable image in the window
MOSAIC_MISS_EDGE = "#6e6e6e"         # hatching of reachable tiles left out of the mosaic
MOSAIC_MAX_SEAM_RINGS = 400          # seamlines are drawn only up to this many image fragments per
                                     # panel; beyond it (e.g. a cloud-fragmented mosaic) colour alone is legible


def tile_first_delivery(selected, store, mask=None):
    """Per ROI tile (index into _tile_ids): delivery time (minutes) of the first
    selected image that covers it (utility > 0, not cloud-obscured) and that
    image's position in delivery order; NaN / -1 where no selected image does.
    Images are ordered exactly as in compute_completion_latency (delivery time,
    ties keep selection order). Returns (first, source, ordered keys)."""
    util = np.where(mask.flags, 0.0, store.util) if mask is not None else store.util

    def _delivery(key):
        w = store.wait_lookup.get(key)
        return key[1] + (np.nan if w is None else w)

    first = np.full(len(_tile_ids), np.nan)
    source = np.full(len(_tile_ids), -1, dtype=np.int64)
    order = sorted(selected, key=delivery_sort_key(store)) if 'delivery_sort_key' in globals() else \
        sorted(selected, key=_delivery)
    for j, key in enumerate(order):
        d = _delivery(key)
        if np.isnan(d):
            continue
        s = store.rows(store.index[key])
        t = store.tile_idx[s][util[s] > 0.0]
        new = t[np.isnan(first[t])]
        first[new] = d
        source[new] = j
    return first, source, order


def reachable_tiles(store, mask=None):
    """Tiles with at least one usable observation among ALL images (as Cell 7C's reachable_tile_count)."""
    util = np.where(mask.flags, 0.0, store.util) if mask is not None else store.util
    out = np.zeros(len(_tile_ids), dtype=bool)
    out[store.tile_idx[util > 0.0]] = True
    return out


def _hour_bounds(t_max_h: float, max_bins: int = 8):
    """Colour-bin edges in hours from 0 to a 'nice' value >= t_max_h."""
    for step in (0.25, 0.5, 1, 2, 3, 4, 6, 12, 24, 48, 96, 168):
        if t_max_h / step <= max_bins:
            top = step * max(1, math.ceil(t_max_h / step - 1e-9))
            return np.arange(0.0, top + 1e-9, step)
    return np.linspace(0.0, max(t_max_h, 1.0), max_bins + 1)


def plot_mosaic_map(store, mask, sel_clear, sel_cloud, policy_label: str, out_dir: str,
                    stem: str = "arf_fig7_mosaic_map"):
    """Two-panel mosaic map (clear sky, cloud) and the per-tile CSV. Returns a summary dict."""
    from matplotlib.colors import BoundaryNorm
    tiles = list(_tile_ids)
    first_c, src_c, order_c = tile_first_delivery(sel_clear, store)
    first_k, src_k, order_k = tile_first_delivery(sel_cloud, store, mask)
    reach_c, reach_k = reachable_tiles(store), reachable_tiles(store, mask)

    outline_geoms, outline_label, outline_source = roi_outline_for_maps()
    ext = roi_map_extent(target_tiles, outline_geoms)
    items = [("scale", _SCALE_BAR_BOX, ["lower right", "lower left", "upper left", "upper right"]),
             ("north", _NORTH_BOX, ["upper right", "upper left", "lower right", "lower left"])]
    if "plan_furniture" in globals():              # v1.4.0: placed off the ROI, frame grown if needed
        ext, corners, furniture = plan_furniture(ext, map_obstacle(target_tiles, outline_geoms), items)
    else:                                          # Cell 3 of v1.3.1
        corners = choose_corners(ext, np.array([h3.cell_to_latlng(t)[::-1] for t in tiles]), items)
        furniture = None
    aspect = panel_aspect(ext)
    polys = tile_polygons(tiles)
    done = np.concatenate([first_c[~np.isnan(first_c)], first_k[~np.isnan(first_k)]]) / 60.0
    bounds = _hour_bounds(float(done.max()) if done.size else float(globals().get("WINDOW_HOURS", 48)))
    cmap = plt.get_cmap(MOSAIC_COLORMAP)
    norm = BoundaryNorm(bounds, cmap.N)
    big = len(tiles) > MAP_RASTER_TILES
    knowledge = globals().get("CLOUD_KNOWLEDGE", "realistic")
    _cm = globals().get("cloud_model_info") or {}
    _cm_txt = (f"correlated, L = {_cm.get('length_km'):g} km, T = {_cm.get('time_h'):g} h, "
               if _cm.get("model") == "correlated" else "independent, " if _cm.get("model") else "")
    # v1.5.3: the cloud settings moved from panel (b)'s title, which made the map
    # 7.75 in wide (two IEEE columns: 7.16 in), to the footer.
    cloud_title = "(b) Cloud: " + ("Fail-Signal re-planning" if knowledge != "oracle" else "cloud known in advance")
    cloud_note = (f"Panel (b) cloud: {_cm_txt}p = {CLOUD_PROBABILITY:g}; "
                  + ("each image's cloud is seen on arrival and the plan is updated (Fail-Signal)"
                     if knowledge != "oracle" else "cloud known to the selector in advance"))
    panels = [("(a) Clear sky", first_c, src_c, order_c, reach_c, sel_clear),
              (cloud_title, first_k, src_k, order_k, reach_k, sel_cloud)]

    summary = {}
    with plt.rc_context(_MAP_RC):
        fig, axes, layout = map_figure([aspect, aspect], top=0.52, bottom=1.12)
        for i, (ax, (title, first, src, order, reach, sel)) in enumerate(zip(axes, panels)):
            setup_map_axes(ax, ext)
            ax.set_facecolor("white")
            cov = ~np.isnan(first)
            miss = reach & ~cov
            none = ~reach & ~cov
            if cov.any():
                pc = PolyCollection([polys[j] for j in np.flatnonzero(cov)], array=first[cov] / 60.0,
                                    cmap=cmap, norm=norm, edgecolors="face", linewidths=0.2, zorder=2)
                pc.set_rasterized(big)
                ax.add_collection(pc)
            if miss.any():
                pc = PolyCollection([polys[j] for j in np.flatnonzero(miss)], facecolors="white",
                                    edgecolors=MOSAIC_MISS_EDGE, linewidths=0.0, hatch="//////", zorder=2)
                pc.set_rasterized(big)
                ax.add_collection(pc)
            if none.any():
                pc = PolyCollection([polys[j] for j in np.flatnonzero(none)], facecolors=MOSAIC_NONE_COLOR,
                                    edgecolors="face", linewidths=0.2, zorder=2)
                pc.set_rasterized(big)
                ax.add_collection(pc)
            seams = []
            for j in np.unique(src[src >= 0]):
                seams += geo_lines([tiles_outline(_tile_ids[src == j])])
            seams_drawn = 0 < len(seams) <= MOSAIC_MAX_SEAM_RINGS
            if seams_drawn:
                lc = LineCollection(seams, colors="0.12", linewidths=0.45 if len(tiles) <= MAP_EDGE_TILES else 0.35,
                                    zorder=3)
                ax.add_collection(lc)
            draw_outline(ax, outline_geoms, lw=0.8, zorder=4)
            if i == 0:
                add_scale_bar(ax, ext, corner=corners["scale"])
                add_north_arrow(ax, corner=corners["north"])
            n_src = int(np.unique(src[src >= 0]).size)
            pct = 100.0 * float(cov.sum()) / len(tiles)
            t_done = float(np.nanmax(first)) / 60.0 if cov.any() else None
            reach_pct = 100.0 * float(reach.sum()) / len(tiles)
            stats = (f"{n_src} of {len(sel)} selected images supply tiles"
                     + ("" if seams_drawn or not seams else f" (seamlines omitted: {len(seams):,} fragments)") + "\n"
                     + (f"{pct:.1f}% of tiles covered ({reach_pct:.1f}% reachable); complete at {t_done:.1f} h"
                        if t_done is not None else f"No tile covered ({reach_pct:.1f}% reachable)"))
            ax.text(0.0, 1.012, stats, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.3,
                    linespacing=1.3, color="0.15")
            ax.set_title(title, loc="left", pad=21)
            summary["clear" if i == 0 else "cloud"] = {
                "n_selected": len(sel), "n_source_images": n_src, "tiles_covered_pct": pct,
                "tiles_reachable_pct": reach_pct, "mosaic_complete_h": t_done,
                "image_fragments": len(seams), "seamlines_drawn": seams_drawn}

        # Colour bar (left) and legend (right) under the panels.
        fw, fh, bot = layout["fig_w"], layout["fig_h"], layout["bottom"]
        cax = fig.add_axes([0.52 / fw, (bot - 0.52) / fh, 2.7 / fw, 0.085 / fh])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                          ticks=bounds, spacing="proportional")
        cb.ax.tick_params(labelsize=6.5, width=0.5, length=2.5)
        cb.outline.set_linewidth(0.5)
        cb.set_ticklabels([f"{b:g}" for b in bounds])
        # v1.5.0: with no downlink wait anywhere, delivery time IS acquisition time; say so.
        _acq_only = bool(np.all(np.nan_to_num(np.array(list(store.wait_lookup.values()), dtype=np.float64),
                                              nan=1.0) == 0.0)) if len(store) else False
        cb.set_label(("Acquisition" if _acq_only else "Delivery") + " time of the first usable image "
                     "(h from analysis start)", fontsize=7, labelpad=2)
        handles = [mpatches.Patch(facecolor="white", edgecolor=MOSAIC_MISS_EDGE, hatch="//////", lw=0.4,
                                  label="Reachable, not in mosaic"),
                   mpatches.Patch(facecolor=MOSAIC_NONE_COLOR, edgecolor="none", label="No usable image in window"),
                   mlines.Line2D([], [], color="0.12", lw=0.6, label="Image boundary (seamline)"),
                   mlines.Line2D([], [], color=MAP_COLORS["boundary"], lw=0.8, label=outline_label)]
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.50, (bot - 0.36) / fh), ncol=2,
                   frameon=False, handlelength=1.8, columnspacing=1.2, handletextpad=0.5, fontsize=7)
        cite, short = _wdpa_citation(globals().get("roi_boundary_report")) if '_wdpa_citation' in globals() else (None, None)
        credits = [policy_label, f"Boundary: {short} (UNEP-WCMC & IUCN)" if short else None, "Coordinates: WGS 84"]
        fig.text(0.52 / fw, 0.04 / fh, cloud_note + "\n" + " · ".join(c for c in credits if c), fontsize=6,
                 color="0.4", ha="left", va="bottom", linespacing=1.45)
        for fmt in ("png", "pdf"):
            fig.savefig(os.path.join(out_dir, f"{stem}.{fmt}"), dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig)

    # Per-tile table (H3 IDs are derived grid cells, not WDPA data).
    def _src_cols(src, order):
        sat = np.array([order[j][0] if j >= 0 else "" for j in src], dtype=object)
        tmin = np.array([order[j][1] if j >= 0 else -1 for j in src], dtype=np.int64)
        return sat, tmin
    sat_c, t_c = _src_cols(src_c, order_c)
    sat_k, t_k = _src_cols(src_k, order_k)
    pd.DataFrame({"tile": _tile_ids, "reachable_clear": reach_c, "first_delivery_min_clear": first_c,
                  "source_sat_clear": sat_c, "source_t_min_clear": t_c,
                  "reachable_cloud": reach_k, "first_delivery_min_cloud": first_k,
                  "source_sat_cloud": sat_k, "source_t_min_cloud": t_k}).to_csv(
        os.path.join(out_dir, "arf_mosaic_tiles.csv"), index=False)
    summary.update(policy=policy_label, colour_bins_h=[float(b) for b in bounds], outline_source=outline_source,
                   extent_lon_lat=[round(v, 5) for v in ext],
                   furniture=None if furniture is None else dict(furniture, anchors=corners))
    return summary


_mosaic_summary = None
_map_helpers = ("roi_map_extent", "map_figure", "setup_map_axes", "tile_polygons", "roi_outline_for_maps",
                "add_scale_bar", "add_north_arrow", "draw_outline", "tiles_outline", "geo_lines")
if all(h in globals() for h in _map_helpers):
    if MOSAIC_MAP_POLICY == "time_aware":
        _mosaic_sel = (sel_ta_clear, sel_ta_cloud)
        _tau = globals().get("LATENCY_DECAY_TAU_MIN")
        _policy_label = (f"Selector: time-aware greedy (k = {MAX_SELECTION_STEPS}, decay {LATENCY_DECAY_RATE:g}"
                         + ("" if _tau is None else f", time scale {_tau / 60:.1f} h") + ") with gap-filling")
    else:
        _mosaic_sel = (sel_events_clear, sel_events_cloud)
        _policy_label = f"Selector: coverage-first greedy (k = {MAX_SELECTION_STEPS}) with gap-filling"
    try:
        _mosaic_summary = plot_mosaic_map(event_data, cloud_mask, *_mosaic_sel, _policy_label,
                                          OUT_DIR if 'OUT_DIR' in globals() else os.getcwd())
        logger.info("Saved: arf_fig7_mosaic_map.png/.pdf and arf_mosaic_tiles.csv")
    except Exception as e:
        logger.warning(f"Mosaic map not generated: {e}")
else:
    logger.warning("Mosaic map skipped: run Cell 3 first (it defines the map helpers this figure uses).")


# =============================================================================
# FIGURE METADATA
# =============================================================================

_metadata = {
    "figures_generated": [
        "arf_fig1_coverage_progression",
        "arf_fig2_quality_progression",
        "arf_fig3_latency_analysis",
        "arf_fig4_coverage_budget",
        "arf_manuscript_figures (composite)",
    ] + (["arf_fig7_mosaic_map"] if _mosaic_summary else []),
    "mosaic_map": _mosaic_summary,
    # Satellites actually simulated (Cell 2). Until v1.3.1 this was ACTIVE_N,
    # which real-fleet runs leave at its default (1000) while the fleet is all
    # operational satellites of the chosen fleet (e.g. 97).
    "n_constellation": len(constellation) if 'constellation' in globals() else 'unknown',
    "fleet": globals().get("FLEET", "starlink"),
    "scenario": str(globals().get("SCENARIO_N", globals().get("ACTIVE_N", "unknown"))),
    "max_achievable_pct": max_achievable_pct,
    # "guarantee_pct" removed: max_achievable_pct * (1 - 1/e) is not a
    # valid guarantee figure -- see plot_coverage_progression's
    # docstring. The qualitative (1-1/e) guarantee claim (without a
    # fabricated percentage) is still recorded correctly under
    # "approximation_guarantee" in the hypothesis_validation export.
    "final_coverage_clear": cov_clear_pct[-1] if cov_clear_pct else None,
    "final_coverage_cloud": cov_cloud_pct[-1] if cov_cloud_pct else None,
    "latency_clear_min": latency_clear,
    "latency_cloud_min": latency_cloud,
    # v1.5.2: what Figure 3 draws (full window, clear sky, coverage-first)
    "figure3": {"h1_level_pct_of_reachable": _H1L7, "reachable_tiles": _reach7,
                "physical_floor_min": _floor7, "h1a_margin_min": _margin7,
                "time_to_level_min_by_k": {str(k): v for k, v in _t_level_by_k.items()},
                "figure_widths_in": {"single": IEEE_ONE_COL_IN, "composite": IEEE_TWO_COL_IN}},
    "h1_target_min": 10.0,
    # 10.0 matches H1_LATENCY_THRESHOLD_MIN (Cell 8) -- duplicated here
    # rather than referenced, since Cell 8 runs after this cell and
    # that constant doesn't exist yet at this point. If the H1
    # threshold ever changes, both this line and Cell 8's constant
    # need updating together.
    "h1_supported": latency_clear is not None and latency_clear < 10.0,
    "h1_supported_note": "original H1 (< 10 min) only; H1a/H1b are in arf_hypothesis_validation.json (Cell 8)",
}

# Save metadata for audit trail
try:
    import json
    with open(os.path.join(OUT_DIR if 'OUT_DIR' in globals() else os.getcwd(), 'arf_cell7_metadata.json'), 'w') as f:
        json.dump(_metadata, f, indent=2, default=str)
    logger.info("Figure metadata saved to arf_cell7_metadata.json")
except Exception as e:
    logger.warning(f"Could not save metadata: {e}")


# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("  CELL 7: Manuscript Graph Generator")
print("=" * 70)
print(f"  Figures Generated:")
print(f"    1. Coverage Progression (Scenario A vs. B)")
print(f"    2. Quality Progression (Scenario A vs. B)")
print(f"    3. Latency vs. Budget (hours; physical floor and H1a margin)")
print(f"    4. Coverage vs. Budget Trade-off")
print(f"    5. Composite Manuscript Figure (2x2, {IEEE_TWO_COL_IN} in wide)")
print(f"    7. Mosaic Map ({MOSAIC_MAP_POLICY.replace('_', '-')} selector)" if _mosaic_summary else "    7. Mosaic Map: not generated (see log)")
print(f"  ─────────────────────────────────────────────────────────")
print(f"  Geometric Ceiling:         {max_achievable_pct:.2f}%")
# v1.5.2: "is not None" (a latency of 0.0 min is a valid value and was not printed)
if cov_clear_pct:
    print(f"  Final Coverage (Clear):    {cov_clear_pct[-1]:.2f}%")
if cov_cloud_pct:
    print(f"  Final Coverage (Cloud):    {cov_cloud_pct[-1]:.2f}%")
if latency_clear is not None:
    print(f"  Latency (Clear):           {latency_clear:.1f} min ({latency_clear / 60:.1f} h)")
if latency_cloud is not None:
    print(f"  Latency (Cloud):           {latency_cloud:.1f} min ({latency_cloud / 60:.1f} h)")
if _floor7 is not None:
    print(f"  Floor to {_H1L7:g}% reachable: {_floor7:.1f} min; H1a margin {_margin7:.1f} min")
print(f"  Original H1 (< 10 min):    {'SUPPORTED' if _metadata.get('h1_supported') else 'NOT SUPPORTED'}"
      f"  (superseded, kept for the record; H1a/H1b: Cell 8)")
print(f"  ─────────────────────────────────────────────────────────")
print(f"  Output Formats:            PNG (300 DPI), PDF (vector)")
print(f"  Composite Figure:          arf_manuscript_figures.png/pdf")
print("=" * 70)

print(f"\n=== CELL 7 COMPLETE ===")

