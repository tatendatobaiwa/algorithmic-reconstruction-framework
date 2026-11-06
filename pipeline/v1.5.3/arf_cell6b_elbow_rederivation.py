# =============================================================================
# CELL 6B: RE-DERIVE SELECTION BUDGET (k) AGAINST THE CURRENT SCENARIO
# =============================================================================
# Cell 6 runs with a fixed selection budget k = MAX_SELECTION_STEPS
# (Cell 0). This cell independently re-derives a data-driven k for
# whichever scenario (ACTIVE_ROI, ACTIVE_N) is currently loaded, using
# the Kneedle elbow method (Satopaa, Albrecht, Irwin & Raghavan 2011, ICDCS
# Workshops, 166-171, doi:10.1109/ICDCSW.2011.20; here its chord-distance
# form) on a real coverage-vs-k sweep -- the point of diminishing returns for
# this specific candidate-event set, not an assumption carried over from a
# different scenario's derivation.
# This matters because candidate-event counts vary substantially
# between scenarios (region size, fleet density, and -- for
# spacetrack_history ROIs -- how many satellites existed at the
# historical epoch all affect this), so a k tuned for one scenario is
# not automatically appropriate for another.
#
# Run this after Cell 6 (needs event_data, run_submodular_greedy_selection,
# total_target_tiles already in the kernel namespace). Results are
# recorded in scenario_results under this run's scenario key (ROI and
# N or fleet tag, Cell 6), alongside Cell 6's own entry, so re-running
# this cell for several scenarios builds up a comparable ledger.
# The swept budgets stop at 100 (or the number of candidate events): a
# fleet that saturates beyond k = 100 shows no knee here (v1.5.2 note).
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 6B is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["6B"] = _ARF_CELL_VERSION

import numpy as np


def compute_elbow_point(k_values, coverage_values):
    """Kneedle chord-distance method -- point of maximum perpendicular
    distance from the line connecting the first and last (k, coverage)
    points.

    Degenerate cases -- fewer than 3 points, or a flat curve (coverage
    already saturated at the smallest budget swept, as for a small ROI a
    single strip can cover) -- have no knee. Normalising a flat curve
    divides 0 by 0, and argmax over the resulting NaNs silently returns
    index 0, reporting an arbitrary "elbow". Instead, return the smallest k
    that reaches the maximum coverage (the saturation budget).
    """
    k_arr = np.array(k_values, dtype=float)
    cov_arr = np.array(coverage_values, dtype=float)
    if k_arr.size < 3 or np.ptp(k_arr) == 0 or np.ptp(cov_arr) == 0:
        i = int(np.argmax(cov_arr >= cov_arr.max() - 1e-9))
        return int(k_arr[i]), float(cov_arr[i]), np.zeros_like(cov_arr)

    k_norm = (k_arr - k_arr.min()) / (k_arr.max() - k_arr.min())
    cov_norm = (cov_arr - cov_arr.min()) / (cov_arr.max() - cov_arr.min())

    x1, y1 = k_norm[0], cov_norm[0]
    x2, y2 = k_norm[-1], cov_norm[-1]

    distances = np.abs((y2 - y1) * k_norm - (x2 - x1) * cov_norm + x2 * y1 - y2 * x1) / \
                np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)

    elbow_idx = int(np.argmax(distances))
    return int(k_arr[elbow_idx]), float(cov_arr[elbow_idx]), distances


print(f"=== RE-DERIVING SELECTION BUDGET FOR ROI={ACTIVE_ROI!r}, "
      f"{'N=' + str(ACTIVE_N) if globals().get('FLEET', 'starlink') == 'starlink' else 'fleet=' + str(FLEET)} ===")
print(f"Candidate events available: {len(event_data)}")

# Sweep k from small values up to full exhaustion (capped for speed --
# same size-aware cap used elsewhere in this pipeline).
_sweep_ceiling = min(len(event_data), 100)
k_values = sorted(set(
    list(range(2, min(20, _sweep_ceiling) + 1, 2)) +
    list(range(20, _sweep_ceiling + 1, 5)) +
    [_sweep_ceiling]
))
k_values = [k for k in k_values if k <= len(event_data)]

print(f"Sweeping k = {k_values}")

# One run at the largest budget: greedy is incremental, so budget k's
# coverage is the k-th entry of this history (see build_sweep_results).
_, _, _, _hist_full = (run_submodular_greedy_selection(event_data, max_steps=max(k_values))
                       if k_values else (None, None, None, []))
coverage_at_k = [(_hist_full[:k][-1]['coverage_pct'] if _hist_full[:k] else 0.0) for k in k_values]
# Exact saturation budget: first greedy step reaching the run's final coverage.
saturation_k = (next(h['step'] for h in _hist_full
                     if h['coverage_pct'] >= _hist_full[-1]['coverage_pct'] - 1e-9)
                if _hist_full else None)

print("\nk     coverage%")
for k, cov in zip(k_values, coverage_at_k):
    print(f"{k:<5} {cov:.2f}")

if k_values:
    elbow_k, elbow_cov, distances = compute_elbow_point(k_values, coverage_at_k)
else:
    elbow_k, elbow_cov, distances = None, 0.0, np.zeros(0)

print(f"\n=== ELBOW POINT: ROI={ACTIVE_ROI!r} ===")
print(f"k = {elbow_k}, coverage = {elbow_cov:.2f}%")
print(f"Saturation budget (first k reaching final coverage "
      f"{_hist_full[-1]['coverage_pct'] if _hist_full else 0.0:.2f}%): k = {saturation_k}")
if k_values and not np.any(distances):
    print("(No diminishing-returns knee: coverage is flat across the swept budgets, "
          "so the elbow reported is the saturation budget, not a Kneedle knee.)")
print(f"(Configured MAX_SELECTION_STEPS for this run: {MAX_SELECTION_STEPS} -- "
      f"compare against the elbow above to see whether the fixed budget "
      f"Cell 6 actually used is under-, over-, or well-matched to this "
      f"scenario's real diminishing-returns point.)")

# Record in the ledger, keyed the same way Cell 6 records its own
# results, so both entries for this (ROI, N) merge into one record
# rather than each cell silently overwriting the other's fields.
if 'scenario_results' not in globals():
    scenario_results = {}
scenario_results.setdefault(_scenario_key(), {})
scenario_results[_scenario_key()]['real_elbow_k'] = elbow_k
scenario_results[_scenario_key()]['real_elbow_coverage_pct'] = elbow_cov
scenario_results[_scenario_key()]['real_k_sweep'] = dict(zip(k_values, coverage_at_k))
scenario_results[_scenario_key()]['saturation_k'] = saturation_k

print(f"\n=== CELL 6B COMPLETE ===")

