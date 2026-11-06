"""Build release 1.5.3 from release 1.5.2 by exact, counted text replacements.

Presentation only: figures (Cells 7, 7C), labels of superseded hypotheses and
the per-run headline (Cells 7, 7C, 8, 9), the version strings of all cells and
the Cell 0 workflow note. No line that computes a saved number is touched; the
replays in compare_v152_v153.py check that every CSV and every JSON value is
unchanged.

Usage: make_v153.py IN_1.5.2 OUT_1.5.3
"""
import re
import sys

SRC, DST = sys.argv[1], sys.argv[2]
src = open(SRC, encoding="utf-8").read()
EXPECTED_SHA_IN = "22896e4f1e834e355fa97bcd25cb6f1ccef48a1491ab7f7189f47774d882c313"
import hashlib
assert hashlib.sha256(src.encode("utf-8")).hexdigest() == EXPECTED_SHA_IN, "input is not the shipped 1.5.2"


def rep(old, new, count=1):
    global src
    n = src.count(old)
    assert n == count, f"expected {count} occurrence(s), found {n}: {old[:90]!r}"
    src = src.replace(old, new)


def block(start, end, new):
    """Replace the text from `start` (included) to `end` (excluded); both unique."""
    global src
    assert src.count(start) == 1, f"start marker not unique: {start[:80]!r}"
    i = src.index(start)
    assert src.count(end, i) >= 1, f"end marker missing: {end[:80]!r}"
    j = src.index(end, i)
    src = src[:i] + new + src[j:]


# ----------------------------------------------------------------------------- version strings
rep('__version__ = "1.5.2"   # 1.5.2: per-cell audit --',
    '__version__ = "1.5.3"   # 1.5.3: presentation only, every computed value as in 1.5.2 -- Figure 3 and '
    'composite panel (c) zoomed to their data with the legend under the axes, Figure 7 panel (b) title '
    'shortened so the map fits two IEEE columns, its cloud settings moved to the footer (Cell 7); Figures 5 '
    'and 9 at the H1a level of reachable tiles, Figure 5 in hours and Figure 9 in minutes on linear axes, '
    'legends under the panels (Cell 7C); superseded hypotheses labelled as such in the console, the report, '
    'the JSON and the pipeline log, and each run\'s own statuses printed first with ledger-wide statuses '
    'labelled as pooled (Cells 7, 7C, 8, 9); 1.5.2: per-cell audit --')
rep('_ARF_CELL_VERSION = "1.5.2"', '_ARF_CELL_VERSION = "1.5.3"', count=13)

# ----------------------------------------------------------------------------- Cell 0: workflow note
rep("#   3. Results accumulate in the scenario_results ledger until the kernel\n"
    "#      restarts; each scenario also writes its own output folder.\n",
    "#   3. Results accumulate in the scenario_results ledger until the kernel\n"
    "#      restarts; each scenario also writes its own output folder. From v1.5.3\n"
    "#      each run's report leads with that run's own statuses; an \"overall\"\n"
    "#      status pools every scenario in the ledger and says so.\n")

# ----------------------------------------------------------------------------- Cell 7: Figure 3
block("def _draw_latency(ax, sweep_data, t_level_by_k, floor_min, margin_min, level_pct):\n",
      "# =============================================================================\n"
      "# FIGURE 4: COVERAGE VS. BUDGET TRADE-OFF",
      '''def _draw_latency(ax, sweep_data, t_level_by_k, floor_min, margin_min, level_pct):
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


''')

# k axes: whole numbers; an empty budget sweep says why (v1.5.3)
rep("def _grid(ax):\n"
    "    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5, color=COLORS['grid'])\n",
    "def _grid(ax):\n"
    "    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5, color=COLORS['grid'])\n"
    "    # v1.5.3: every figure that calls this counts images (k) along x: whole-number ticks\n"
    "    # (a run with three candidate images showed k = 1.0, 1.5, 2.0, ...)\n"
    "    from matplotlib.ticker import MaxNLocator\n"
    "    ax.xaxis.set_major_locator(MaxNLocator(integer=True))\n"
    "\n"
    "\n"
    "def _no_sweep_note(ax):\n"
    "    \"\"\"v1.5.3: say why a budget panel (Figures 3 and 4) is empty: the sweep\n"
    "    starts at k = 5, and a window with fewer candidate images has no budget to\n"
    "    sweep (open data alone, CKGR + Okavango, 15 October 2026: 3 images).\"\"\"\n"
    "    _ed = globals().get(\"event_data\")\n"
    "    n = len(_ed) if _ed is not None else 0\n"
    "    ax.text(0.5, 0.5, f\"No budget sweep: the window offers {n} candidate image{'' if n == 1 else 's'},\\n\"\n"
    "                      f\"fewer than the smallest budget swept (k = 5)\", transform=ax.transAxes,\n"
    "            ha=\"center\", va=\"center\", fontsize=7.5, color=\"0.3\")\n")
rep("    ax.set_xlabel('Selection Budget (k)')\n"
    "    ax.set_ylabel('Coverage (%)')\n",
    "    if not sweep_data:\n"
    "        _no_sweep_note(ax)\n"
    "    ax.set_xlabel('Selection Budget (k)')\n"
    "    ax.set_ylabel('Coverage (%)')\n")
# Figure 8: 7.0 in, so labels that grow (a flat H2 difference prints -0.04 ... 0.04) stay within 7.16 in
rep("        fig8, axs8 = plt.subplots(1, _n_pan, figsize=(7.16, 2.6),\n",
    "        # v1.5.3: 7.0 in (was 7.16): with the tight bounding box the saved figure grew to 7.18 in when\n"
    "        # panels (b) and (c) carry longer tick labels (open data alone, 15 October 2026)\n"
    "        fig8, axs8 = plt.subplots(1, _n_pan, figsize=(7.0, 2.6),\n")

# composite figure: panel (c) legend under the panels
rep("    fig_composite, axes = plt.subplots(2, 2, figsize=(IEEE_TWO_COL_IN, 5.6))\n",
    "    fig_composite, axes = plt.subplots(2, 2, figsize=(IEEE_TWO_COL_IN, 6.0))\n")
rep("    _draw_latency(axes[1, 0], sweep_data, _t_level_by_k, _floor7, _margin7, _H1L7)\n"
    "    axes[1, 0].set_title('(c) Latency vs. Budget (clear sky)', loc='left')\n"
    "    _draw_budget(axes[1, 1], sweep_data, max_achievable_pct, cloud_label='Scenario B (Cloud)')\n"
    "    axes[1, 1].set_title('(d) Coverage vs. Budget Trade-off', loc='left')\n"
    "    fig_composite.tight_layout()\n",
    "    _h3c, _l3c = _draw_latency(axes[1, 0], sweep_data, _t_level_by_k, _floor7, _margin7, _H1L7)\n"
    "    axes[1, 0].set_title('(c) Latency vs. Budget (clear sky)', loc='left')\n"
    "    _draw_budget(axes[1, 1], sweep_data, max_achievable_pct, cloud_label='Scenario B (Cloud)')\n"
    "    axes[1, 1].set_title('(d) Coverage vs. Budget Trade-off', loc='left')\n"
    "    # v1.5.3: panel (c)'s legend under the panels, from panel (c)'s left edge\n"
    "    fig_composite.tight_layout(rect=(0, 0.07 if _h3c else 0, 1, 1))\n"
    "    if _h3c:\n"
    "        fig_composite.legend(_h3c, _l3c, loc='lower left', ncol=2, frameon=False, fontsize=7.5,\n"
    "                             bbox_to_anchor=(axes[1, 0].get_position().x0, 0.0))\n")

# Cell 7 console summary: the original H1 is superseded
rep("print(f\"  Original H1 (< 10 min):    {'SUPPORTED ✅' if _metadata.get('h1_supported') else 'NOT SUPPORTED ❌'}\"\n"
    "      f\"  (H1a/H1b: Cell 8)\")\n",
    "print(f\"  Original H1 (< 10 min):    {'SUPPORTED' if _metadata.get('h1_supported') else 'NOT SUPPORTED'}\"\n"
    "      f\"  (superseded, kept for the record; H1a/H1b: Cell 8)\")\n")

# ----------------------------------------------------------------------------- Cell 7: Figure 7 title
rep('    cloud_title = (f"(b) Cloud ({_cm_txt}p = {CLOUD_PROBABILITY:g}): "\n'
    '                   + ("Fail-Signal re-planning" if knowledge != "oracle" else "cloud known in advance"))\n',
    '    # v1.5.3: the cloud settings moved from panel (b)\'s title, which made the map\n'
    '    # 7.75 in wide (two IEEE columns: 7.16 in), to the footer.\n'
    '    cloud_title = "(b) Cloud: " + ("Fail-Signal re-planning" if knowledge != "oracle" else "cloud known in advance")\n'
    '    cloud_note = (f"Panel (b) cloud: {_cm_txt}p = {CLOUD_PROBABILITY:g}; "\n'
    '                  + ("each image\'s cloud is seen on arrival and the plan is updated (Fail-Signal)"\n'
    '                     if knowledge != "oracle" else "cloud known to the selector in advance"))\n')
rep('        fig.text(0.52 / fw, 0.04 / fh, " · ".join(c for c in credits if c), fontsize=6, color="0.4",\n'
    '                 ha="left", va="bottom")\n',
    '        fig.text(0.52 / fw, 0.04 / fh, cloud_note + "\\n" + " · ".join(c for c in credits if c), fontsize=6,\n'
    '                 color="0.4", ha="left", va="bottom", linespacing=1.45)\n')

# ----------------------------------------------------------------------------- Cell 7C: Figures 5 and 9
block("    # Figure 5: empirical CDF, across onsets, of the time to the highest common\n",
      "else:\n    logger.info(\"Onset robustness analysis disabled (onset_robustness_enabled=False in Cell 0).\")\n",
      '''    # Figure 5: empirical CDF, across the evaluated onsets, of the time from the
    # onset to the H1a level (H1_REACHABLE_COVERAGE_PCT of the tiles reachable in
    # that onset's horizon): the selectors against the physical floor, like for
    # like. An onset where a selector never reaches the level leaves its curve
    # below 100%.
    # v1.5.3: until v1.5.2 the level was the highest COVERAGE_TARGETS_PCT of ALL
    # ROI tiles (99%), which no curve reaches when the window's passes cannot
    # (15 October 2026: 65-89% of the ROI reachable), so both panels were empty;
    # hours now on a linear axis over the onset horizon (was symlog), one legend
    # under the panels. The random baseline has no time to the reachable level
    # (Cell 7C records it at the ROI targets only), so it is not drawn here.
    _rc7c = {"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
             "legend.fontsize": 7.5, "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42}
    from matplotlib.ticker import MaxNLocator as _MaxNLocator
    _sel_styles = {"earliest_completion": ("Earliest-completion cover (H1a)", "#e7298a", "-"),
                   "time_aware": ("Time-aware greedy", "#1b9e77", "-"),
                   "coverage_first": ("Coverage-first greedy", "#d95f02", "--"),
                   "fcfs": ("First-come-first-served (budget k)", "#7570b3", "-.")}

    def _legend_under(fig, axes, ncol, order):
        """One legend for both panels, under them, in the given order of labels
        (a series may appear in one panel only)."""
        hl = {}
        for _ax in axes:
            for _h, _l in zip(*_ax.get_legend_handles_labels()):
                hl.setdefault(_l, _h)
        labs = [l for l in order if l in hl] + [l for l in hl if l not in order]
        if labs:
            fig.legend([hl[l] for l in labs], labs, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=ncol,
                       frameon=False)

    _floor_lab = "Physical floor (all passes)"
    _hz_h = float(WINDOW_HOURS - ONSET_SHIFT_HOURS)
    with plt.rc_context(_rc7c):
        fig5, axes5 = plt.subplots(1, 2, figsize=(7.16, 3.35), sharey=True)
        for ax, sky in zip(axes5, ("clear", "cloud")):
            fl = pd.to_numeric(onset_df[f"best_possible_{sky}_t{_H1R:g}r_min"], errors="coerce")
            n_ev = int(fl.notna().sum())
            at_floor, never = [], []
            # The floor as a wide grey band under the selectors (a selector at the
            # floor runs inside it); the H1a selector drawn last, on top.
            for name, (lab, col, ls, lw) in [("best_possible", (_floor_lab, "0.72", "-", 4.0))] + \\
                    [(n, (*_sel_styles[n], 1.3)) for n in reversed(list(_sel_styles))]:
                v = pd.to_numeric(onset_df[f"{name}_{sky}_t{_H1R:g}r_min"], errors="coerce")
                x = np.sort(v[fl.notna() & v.notna()].to_numpy(dtype=float)) / 60.0
                if x.size and n_ev:
                    ax.step(np.r_[0.0, x], 100.0 * np.arange(0, x.size + 1) / n_ev, where="post", color=col,
                            ls=ls, lw=lw, label=lab, solid_capstyle="butt")
                    _v, _f = v[fl.notna()], fl[fl.notna()]
                    if name != "best_possible" and _v.notna().all() and bool(np.all(np.abs(_v - _f) <= 1e-9)):
                        at_floor.insert(0, lab.split(" (")[0])
                elif n_ev:
                    never.insert(0, lab.split(" (")[0])
            _notes = (["At the floor at every onset:"] + at_floor if at_floor else []) + \\
                     (["Never at the level:"] + never if never else [])
            if _notes:
                ax.text(0.98, 0.04, "\\n".join(_notes), transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
                        color="0.25")
            ax.set_xlim(0, _hz_h)
            ax.xaxis.set_major_locator(_MaxNLocator(nbins=7, steps=[1, 2, 3, 4, 6, 10]))
            ax.set_ylim(0, 102)
            ax.set_xlabel(f"Time from onset to {_H1R:g}% of reachable tiles (h)")
            ax.set_title(f"({'a' if sky == 'clear' else 'b'}) {'Clear sky' if sky == 'clear' else 'Cloud / Fail-Signal'}"
                         f" ({n_ev} onsets)", loc="left")
            ax.grid(alpha=0.3, lw=0.5)
        axes5[0].set_ylabel("Share of onsets (%)")
        fig5.tight_layout(rect=(0, 0.15, 1, 1))
        _legend_under(fig5, axes5, ncol=3, order=[_floor_lab] + [s[0] for s in _sel_styles.values()])
        for fmt in ("png", "pdf"):
            fig5.savefig(os.path.join(_out_dir, f"arf_fig5_onset_robustness.{fmt}"), dpi=300, bbox_inches="tight")
        plt.show()

        # Figure 9 (v1.5.2): the H1a quantity itself -- per onset, the time to
        # H1_REACHABLE_COVERAGE_PCT of the reachable tiles minus the physical
        # floor's time to it, against the H1a margin. Empirical CDF over the
        # evaluated onsets; an onset where a selector never reaches the level
        # leaves its curve below 100%.
        # v1.5.3: minutes on a linear axis that ends a quarter past the margin
        # (was hours on a symlog axis, where gaps of a few minutes against a
        # 94-min margin lay on the y-axis); the axis starts just left of 0 so a
        # selector at the floor (gap 0) shows off the spine; the largest gap of
        # each selector (or how many onsets fall beyond the axis or never reach
        # the level) is written on the panel; one legend under the panels.
        _xmax9 = max(1.25 * H1_FLOOR_MARGIN_MIN, 10.0)
        _margin_lab = f"H1a margin ({H1_FLOOR_MARGIN_MIN:.0f} min)"
        fig9, axes9 = plt.subplots(1, 2, figsize=(7.16, 3.15), sharey=True)
        for ax, sky in zip(axes9, ("clear", "cloud")):
            fl = pd.to_numeric(onset_df[f"best_possible_{sky}_t{_H1R:g}r_min"], errors="coerce")
            n_ev = int(fl.notna().sum())
            notes = []
            for name in reversed(list(_sel_styles)):               # the H1a selector last, on top
                lab, col, ls = _sel_styles[name]
                v = pd.to_numeric(onset_df[f"{name}_{sky}_t{_H1R:g}r_min"], errors="coerce")
                gap = np.sort(((v - fl)[fl.notna() & v.notna()]).to_numpy(dtype=float))
                if gap.size and n_ev:
                    ax.step(np.r_[0.0, gap], 100.0 * np.arange(0, gap.size + 1) / n_ev, where="post", color=col,
                            ls=ls, lw=1.3, label=lab)
                n_out = n_ev - int((gap <= _xmax9).sum())
                if n_ev:
                    notes.insert(0, f"{lab.split(' (')[0]}: " + (f"largest {round(float(gap.max()), 1):g} min"
                                                                 if not n_out else
                                                                 f"{n_out} beyond {_xmax9:.0f} min or never"))
            ax.axvline(H1_FLOOR_MARGIN_MIN, color="k", ls=":", lw=1.0, label=_margin_lab)
            if notes:
                ax.text(0.97 * H1_FLOOR_MARGIN_MIN, 0.04, "Gap to the floor, per selector\\n" + "\\n".join(notes),
                        transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=7, color="0.25")
            ax.set_xlim(-0.03 * _xmax9, _xmax9)
            ax.set_ylim(0, 102)
            ax.set_xlabel(f"Time to {_H1R:g}% of reachable tiles minus the floor (min)")
            ax.set_title(f"({'a' if sky == 'clear' else 'b'}) {'Clear sky' if sky == 'clear' else 'Cloud'}"
                         f" ({n_ev} onsets)", loc="left")
            ax.grid(alpha=0.3, lw=0.5)
        axes9[0].set_ylabel("Share of onsets (%)")
        fig9.tight_layout(rect=(0, 0.15, 1, 1))
        _legend_under(fig9, axes9, ncol=3, order=[s[0] for s in _sel_styles.values()] + [_margin_lab])
        for fmt in ("png", "pdf"):
            fig9.savefig(os.path.join(_out_dir, f"arf_fig9_h1a_gap.{fmt}"), dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig5)
        plt.close(fig9)
''')

# Cell 7C console: the 10-min target is the superseded H1
rep('    print(f"\\n  {len(onset_df)} onsets | H1 target {_H1_TARGET_MIN:g} min")\n'
    '    print(f"  {\'series\':26s} {\'median\':>8s} {\'p90\':>8s} {\'max\':>8s} {\'% < H1\':>8s}")\n',
    '    print(f"\\n  {len(onset_df)} onsets | superseded H1 target {_H1_TARGET_MIN:g} min (kept for the record)")\n'
    '    print(f"  {\'series\':26s} {\'median\':>8s} {\'p90\':>8s} {\'max\':>8s} {\'% <10min\':>8s}")\n')

# ----------------------------------------------------------------------------- Cell 8: superseded labels
rep("h3_result = validate_h3_redundancy_minimisation(scenario_results)\n",
    "h3_result = validate_h3_redundancy_minimisation(scenario_results)\n"
    "# v1.5.3: the superseded hypotheses carry the label wherever they are printed or\n"
    "# saved. Until v1.5.2 the console and the pipeline log printed, e.g., \"H2\n"
    "# (Reliability via Fail-Signal): SUPPORTED\" with nothing to say that this test\n"
    "# was replaced (and passes whenever p < 0.5).\n"
    "for _r, _note in ((h1_result, \"superseded by H1a and H1b: the original < 10 min completion target, falsified \"\n"
    "                              \"by the real-fleet runs of 2026-10-09\"),\n"
    "                  (h1_ta_result, \"superseded by H1a and H1b: the < 10 min completion target for the \"\n"
    "                                 \"time-aware selector\"),\n"
    "                  (h2_result, \"superseded by H2 (v1.5.0, comparative): a plan with no Fail-Signal passes \"\n"
    "                              \"this ratio test whenever p < 0.5\")):\n"
    "    _r[\"superseded\"], _r[\"superseded_note\"] = True, _note\n")

rep('print("=" * 70)\nprint("HYPOTHESIS VALIDATION")\nprint("=" * 70)\n',
    '# v1.5.3: this run\'s own statuses first. An "overall" status pools every\n'
    '# scenario in the kernel\'s ledger (SUPPORTED only if every evaluable one is),\n'
    '# so a run made after others in one kernel could read "H2 NOT SUPPORTED"\n'
    '# overall while its own scenario was supported (15 October 2026).\n'
    '_this_label = _scenario_label(_scenario_key()) if "_scenario_key" in globals() else None\n'
    '_this_run = {"scenario": _this_label, "ledger_scenarios": [_scenario_label(k) for k in scenario_results],\n'
    '             "statuses": {_hk: ((_hr.get("assessments") or {}).get(_this_label) or {}).get("status", "NOT RUN")\n'
    '                          for _hk, _hr in (("H1a", h1_revised_result), ("H1b", h1b_result),\n'
    '                                           ("H1 (H1a and H1b)", h1_reframed_result), ("H2", h2c_result),\n'
    '                                           ("H3", h3_result))},\n'
    '             "note": "statuses of this run\'s own scenario; the overall_status of each hypothesis pools "\n'
    '                     "every scenario in ledger_scenarios"}\n'
    'hypothesis_validation["this_run"] = _this_run\n'
    '_ov = (f" (overall: pools the {len(scenario_results)} scenarios in this kernel)"\n'
    '       if len(scenario_results) > 1 else "")\n'
    'print("=" * 70)\nprint("HYPOTHESIS VALIDATION")\nprint("=" * 70)\n'
    'print(f"\\nThis run ({_this_label}): " + ", ".join(f"{k} {v}" for k, v in _this_run["statuses"].items()))\n'
    'if len(scenario_results) > 1:\n'
    '    print(f"NOTE: this kernel\'s ledger holds {len(scenario_results)} scenarios; each \'overall\' status below "\n'
    '          f"pools them all. Quote this run\'s line above or the per-scenario lines.")\n')
rep('print(f"\\n{h1_reframed_result[\'hypothesis\']}: {h1_reframed_result[\'overall_status\']}")\n',
    'print(f"\\n{h1_reframed_result[\'hypothesis\']}: {h1_reframed_result[\'overall_status\']}{_ov}")\n')
rep('print(f"\\n{h1_revised_result[\'hypothesis\']}: {h1_revised_result[\'overall_status\']}")\n',
    'print(f"\\n{h1_revised_result[\'hypothesis\']}: {h1_revised_result[\'overall_status\']}{_ov}")\n')
rep('print(f"\\n{h2c_result[\'hypothesis\']}: {h2c_result[\'overall_status\']}")\n',
    'print(f"\\n{h2c_result[\'hypothesis\']}: {h2c_result[\'overall_status\']}{_ov}")\n')
rep('print(f"\\n{h1b_result[\'hypothesis\']}: {h1b_result[\'overall_status\']}")\n',
    'print(f"\\n{h1b_result[\'hypothesis\']}: {h1b_result[\'overall_status\']}{_ov}")\n')
rep("    print(f\"\\n{h_data['hypothesis']}: {h_data['overall_status']}\")\n"
    "    for n_label, a in h_data['assessments'].items():\n"
    "        print(f\"  {n_label}: {a}\")\n",
    "    _sup = \"SUPERSEDED (kept for the record): \" if h_data.get(\"superseded\") else \"\"\n"
    "    print(f\"\\n{_sup}{h_data['hypothesis']}: {h_data['overall_status']}{_ov}\")\n"
    "    if h_data.get(\"superseded\"):\n"
    "        print(f\"  ({h_data['superseded_note']})\")\n"
    "    for n_label, a in h_data['assessments'].items():\n"
    "        print(f\"  {n_label}: {a}\")\n")

# report: this run first, pooled statuses labelled, superseded hypotheses labelled
rep("          f\"Framework Version: {globals().get('__version__', 'unknown')}\", \"\",\n",
    "          f\"Framework Version: {globals().get('__version__', 'unknown')}\", \"\",\n"
    "          _sub, \"THIS RUN\", _sub, \"\", f\"Scenario: {_this_label}\",\n"
    "          *[f\"  {k}: {v}\" for k, v in _this_run[\"statuses\"].items()],\n"
    "          *([f\"NOTE: the ledger holds {len(scenario_results)} scenarios (this kernel ran more than one); each \"\n"
    "             f\"'overall' status under HYPOTHESIS VALIDATION pools them all. This folder's run is the one above.\"]\n"
    "            if len(scenario_results) > 1 else []), \"\",\n")
rep("_lines.append(f\"{h1_reframed_result['hypothesis']}: {h1_reframed_result['overall_status']}\")\n",
    "_lines.append(f\"{h1_reframed_result['hypothesis']}: {h1_reframed_result['overall_status']}{_ov}\")\n")
rep("_lines.append(f\"{h1_revised_result['hypothesis']}: {h1_revised_result['overall_status']}\")\n",
    "_lines.append(f\"{h1_revised_result['hypothesis']}: {h1_revised_result['overall_status']}{_ov}\")\n")
rep("_lines.append(f\"{h1b_result['hypothesis']}: {h1b_result['overall_status']}\")\n",
    "_lines.append(f\"{h1b_result['hypothesis']}: {h1b_result['overall_status']}{_ov}\")\n")
rep("_lines.append(\"ORIGINAL H1 (kept for the record):\")\n"
    "_lines.append(f\"{h1_result['hypothesis']}: {h1_result['overall_status']}\")\n",
    "_lines.append(\"SUPERSEDED HYPOTHESES (kept for the record; not tested by the protocol, not results to quote):\")\n"
    "_lines.append(f\"SUPERSEDED {h1_result['hypothesis']}: {h1_result['overall_status']}{_ov} \"\n"
    "              f\"({h1_result['superseded_note']})\")\n")
rep("_lines.append(f\"{h1_ta_result['hypothesis']}: {h1_ta_result['overall_status']}\")\n",
    "_lines.append(f\"SUPERSEDED {h1_ta_result['hypothesis']}: {h1_ta_result['overall_status']}{_ov} \"\n"
    "              f\"({h1_ta_result['superseded_note']})\")\n")
rep("_lines.append(f\"{h2c_result['hypothesis']}: {h2c_result['overall_status']}\")\n",
    "_lines.append(f\"{h2c_result['hypothesis']}: {h2c_result['overall_status']}{_ov}\")\n")
rep("_lines.append(f\"ORIGINAL {h2_result['hypothesis']} (kept for the record): {h2_result['overall_status']} \"\n"
    "              f\"(requires cloud coverage to maintain > {h2_result['maintenance_threshold_pct']:g}% of clear-sky coverage)\")\n",
    "_lines.append(f\"SUPERSEDED original {h2_result['hypothesis']} (kept for the record): \"\n"
    "              f\"{h2_result['overall_status']}{_ov} (requires cloud coverage to maintain > \"\n"
    "              f\"{h2_result['maintenance_threshold_pct']:g}% of clear-sky coverage; {h2_result['superseded_note']})\")\n")
rep("_lines.append(f\"{h3_result['hypothesis']}: {h3_result['overall_status']} (requires coverage \"\n",
    "_lines.append(f\"{h3_result['hypothesis']}: {h3_result['overall_status']}{_ov} (requires coverage \"\n")
rep("| under H1 in {_s['pct_onsets_under_h1']:.1f}% of onsets\")",
    "| under the superseded 10-min H1 in {_s['pct_onsets_under_h1']:.1f}% of onsets\")", count=2)

# ----------------------------------------------------------------------------- Cell 9: pipeline log
rep("                    if isinstance(h_data, dict) and \"overall_status\" in h_data:\n"
    "                        f.write(f\"{h_data.get('hypothesis', h_key)}: {h_data['overall_status']}\\n\")\n",
    "                    if isinstance(h_data, dict) and \"overall_status\" in h_data:\n"
    "                        _sup = \"SUPERSEDED (kept for the record): \" if h_data.get(\"superseded\") else \"\"\n"
    "                        f.write(f\"{_sup}{h_data.get('hypothesis', h_key)}: {h_data['overall_status']}\\n\")\n")

# ----------------------------------------------------------------------------- checks
assert src.count('_ARF_CELL_VERSION = "1.5.3"') == 13 and '_ARF_CELL_VERSION = "1.5.2"' not in src
assert re.search(r'^__version__ = "1\.5\.3"', src, re.M)
compile(src, DST, "exec")
open(DST, "w", encoding="utf-8").write(src)
print(f"{DST}: written, sha256 {hashlib.sha256(src.encode('utf-8')).hexdigest()}")
