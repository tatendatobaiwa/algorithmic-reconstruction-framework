# =============================================================================
# CELL 7C: OBJECTIVE COMPARISON & ONSET-ROBUSTNESS ANALYSIS (H1)
# =============================================================================
# Two analyses of the H1 latency (mosaic completion time, Cell 6):
#
# 1. DECAY SWEEP: the time-aware objective (Cell 6) at each rate in
#    DECAY_SWEEP, plus rate 0 -- which weights every event equally and is
#    therefore identical to the coverage-first objective. Reports latency
#    against final coverage and mean quality, so LATENCY_DECAY_RATE is one
#    reported point on a trade-off curve, not an unexplained tuned value.
#
# 2. ONSET ROBUSTNESS: one onset time is one sample of orbital geometry; the
#    first pass may arrive early or late by chance. The request time is
#    shifted from 0 to ONSET_SHIFT_HOURS in ONSET_SHIFT_STEP_MIN steps, and for
#    each shift only the passes in the following
#    (WINDOW_HOURS - ONSET_SHIFT_HOURS) hours are used, times re-expressed
#    from the shifted onset. For each shift it records:
#      - best possible: earliest delivery time at which ALL passes so far
#        together reach full reachable coverage (a physical lower bound,
#        independent of budget k and of any algorithm)
#      - coverage-first and time-aware selections (k = MAX_SELECTION_STEPS)
#    under clear sky and the Fail-Signal cloud mask. The cloud mask is the
#    same per-(event, tile) mask as Cell 6, not redrawn per shift.
#
# 3. BASELINES: first-come-first-served (passes in delivery order, kept if
#    they add a new tile) and random selection (median of
#    RANDOM_BASELINE_DRAWS seeded draws), compared with both greedy objectives
#    on the full window and in the onset-robustness analysis.
#
# 4. OPTIMALITY GAP: on ILP_INSTANCES sub-problems (ILP_HORIZON_HOURS of
#    passes, budget ILP_BUDGET_K, clear sky), the exact maximum-utility
#    selection from a mixed-integer program (SciPy/HiGHS) -- the objective the
#    greedy approximates -- reported as each policy's utility / optimum.
#
# EQUAL-COVERAGE LATENCY: completion time is measured to each policy's OWN
# final coverage, so a policy that covers less can look faster (under cloud,
# a budget-limited policy may stop well short of full coverage). Every policy
# and the physical floor are therefore also timed to common coverage levels,
# COVERAGE_TARGETS_PCT (default 95% and 99% of ROI tiles); a policy that never
# reaches a level is recorded as not reaching it. Figure 5 plots the highest
# level. At equal coverage no policy can beat the physical floor.
#
# All policy comparisons in this cell use the budget k WITHOUT the
# Fail-Signal gap-fill: gap-fill adds compensatory passes beyond k, which the
# baselines have no counterpart for, so including it would favour the greedy
# variants. Cell 6 reports H1 with gap-fill; where gap-fill never triggers
# (every tile covered within k) the two agree.
#
# H1 FROM v1.4.0 (tested in Cell 8), per onset, for the level L =
# H1_REACHABLE_COVERAGE_PCT of the tiles reachable in that onset's horizon:
#   H1a, latency non-inferiority: time-aware time to L minus the floor's time
#       to L <= H1_FLOOR_MARGIN_MIN (the pre-specified revised H1).
#   H1b, image economy: images first-come-first-served (FCFS, no budget) has
#       taken delivery of when its mosaic reaches L, against the images of the
#       earliest-completion cover (ECC, budget k), which must also reach L
#       within the margin of the floor. FCFS keeps every image that adds a new
#       tile, so its mosaic at any time is the union of all delivered images:
#       it reaches L exactly at the floor. Onsets where FCFS needs one image
#       cannot show economy and are counted separately.
#   Quality at L (reported, not tested): mean utility of the tiles covered
#       when each policy's delivered mosaic first reaches L.
# A share of onsets is NOT EVALUABLE below H1_MIN_EVALUABLE_ONSETS (Cell 0),
# and is reported with its exact one-sided 95% lower confidence bound
# (Clopper & Pearson 1934). Consecutive onsets share most of their passes, so
# they are not independent trials: the bound is optimistic and the minimum
# count is necessary, not sufficient.
#
# Cost: 2 greedy selections + baselines per shift and sky. Trivial for small
# ROIs; for large ROIs at N=1000, increase ONSET_SHIFT_STEP_MIN or set
# onset_robustness_enabled=False in Cell 0 to trade resolution for time.
#
# Outputs (OUT_DIR): arf_decay_sweep.csv, arf_policy_comparison.csv,
#   arf_onset_robustness.csv, arf_onset_robustness_summary.json,
#   arf_fig5_onset_robustness.png/.pdf, arf_fig9_h1a_gap.png/.pdf (v1.5.2),
#   arf_ilp_gap.csv, arf_cost_latency_frontier.csv, arf_fig6_cost_latency.png/.pdf,
#   arf_cloud_robustness.csv, arf_cloud_robustness_summary.json,
#   arf_cloud_sensitivity.csv, arf_fig8_cloud_robustness.png/.pdf
# =============================================================================
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 7C is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["7C"] = _ARF_CELL_VERSION

import os
import json
import math
import logging
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger("arf.cell7c")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()
# Same H1 threshold as Cell 8's H1_LATENCY_THRESHOLD_MIN (Cell 8 runs later,
# so it is restated here; change both together).
_H1_TARGET_MIN = 10.0

# Revised H1 margin: one orbital period at the fleet's median altitude, unless
# set explicitly in Cell 0. T = 2*pi*sqrt(a^3/mu) (Kepler's third law).
H1_ORBITAL_PERIOD_MIN = float(2.0 * np.pi * np.sqrt((Re_km + CONSTELLATION_ALTITUDE_KM) ** 3 / mu_earth) / 60.0)
H1_FLOOR_MARGIN_MIN = (H1_ORBITAL_PERIOD_MIN if globals().get("H1_FLOOR_MARGIN_MIN_CONFIG") is None
                       else float(H1_FLOOR_MARGIN_MIN_CONFIG))
H1_MARGIN_SOURCE = ("derived: one orbital period at the fleet median altitude "
                    f"{CONSTELLATION_ALTITUDE_KM:.1f} km" if globals().get("H1_FLOOR_MARGIN_MIN_CONFIG") is None
                    else "set in Cell 0")
_H1R = float(globals().get("H1_REACHABLE_COVERAGE_PCT", 95.0))
_H1MIN = float(globals().get("H1_MIN_ONSET_PCT", 90.0))

# Evaluability (v1.4.0). Clopper, C. J. & Pearson, E. S. (1934). The use of
# confidence or fiducial limits illustrated in the case of the binomial.
# Biometrika 26(4), 404-413, doi:10.1093/biomet/26.4.404. The exact one-sided
# lower bound for x successes in n trials at confidence 1 - alpha is the alpha
# quantile of Beta(x, n - x + 1); with x = n it is alpha**(1/n).
H1_CONFIDENCE_ALPHA = 0.05


def cp_lower_bound(x, n, alpha=H1_CONFIDENCE_ALPHA):
    """Exact (Clopper-Pearson) one-sided lower confidence bound, as a fraction,
    for x successes in n trials; 0 if x == 0, None if n == 0 or SciPy is missing."""
    if not n:
        return None
    if x <= 0:
        return 0.0
    try:
        from scipy.stats import beta
    except ImportError:
        return None
    return float(beta.ppf(alpha, x, n - x + 1))


def min_onsets_for(p0_pct, alpha=H1_CONFIDENCE_ALPHA):
    """Smallest n for which n successes in n trials give a one-sided lower bound
    >= p0_pct: alpha**(1/n) >= p0, i.e. n >= ln(alpha) / ln(p0). None for p0 = 100%
    (no finite sample reaches it)."""
    p0 = float(p0_pct) / 100.0
    if p0 >= 1.0:
        return None
    return int(math.ceil(math.log(alpha) / math.log(p0) - 1e-12))


H1_MIN_EVALUABLE_ONSETS = (int(H1_MIN_EVALUABLE_ONSETS_CONFIG)
                           if globals().get("H1_MIN_EVALUABLE_ONSETS_CONFIG") is not None
                           else (min_onsets_for(_H1MIN) or 1))
H1_EVALUABILITY_SOURCE = ("set in Cell 0" if globals().get("H1_MIN_EVALUABLE_ONSETS_CONFIG") is not None else
                          f"derived: smallest n with {H1_CONFIDENCE_ALPHA:g}**(1/n) >= {_H1MIN / 100:g} "
                          f"(all-pass Clopper-Pearson one-sided {100 * (1 - H1_CONFIDENCE_ALPHA):g}% lower bound)")


def subset_by_time(store, mask, t0, t1):
    """
    EventStore (and aligned CloudMask) of the events with t0 <= t_min < t1,
    with times re-expressed relative to t0. Rows and cloud flags are copied
    from the full store, so utilities and cloud cover are unchanged.
    """
    keep = [e for e, k in enumerate(store.keys_) if t0 <= k[1] < t1]
    if not keep:
        empty = EventStore([], np.zeros(1, np.int64), np.zeros(0, np.int32), np.zeros(0), [])
        return empty, (CloudMask(empty, np.zeros(0, bool)) if mask is not None else None)
    rows = np.concatenate([np.arange(store.offsets[e], store.offsets[e + 1]) for e in keep])
    counts = store.counts[keep]
    sub = EventStore([(store.keys_[e][0], store.keys_[e][1] - t0) for e in keep],
                     np.r_[0, np.cumsum(counts)], store.tile_idx[rows], store.util[rows],
                     [store.wait_lookup[store.keys_[e]] for e in keep])
    return sub, (CloudMask(sub, mask.flags[rows]) if mask is not None else None)


def best_possible_completion(store, mask=None):
    """
    Physical lower bound on mosaic completion: with ALL events ordered by
    delivery time, the delivery time at which their union first reaches the
    full coverage reachable from this store (utility > 0, not cloud-obscured).
    Returns (minutes, coverage_pct), or (None, 0.0) if nothing is reachable.
    """
    if len(store) == 0:
        return None, 0.0
    util = np.where(mask.flags, 0.0, store.util) if mask is not None else store.util
    reachable = np.unique(store.tile_idx[util > 0.0]).size
    if reachable == 0:
        return None, 0.0
    deliveries = []
    for e, (sat, t) in enumerate(store.keys_):
        w = store.wait_lookup.get((sat, t))
        if w is not None and not np.isnan(w):
            deliveries.append((t + w, e))
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n_seen = 0
    for d, e in sorted(deliveries):
        s = store.rows(e)
        t = store.tile_idx[s][util[s] > 0.0]
        new = t[~seen[t]]
        if new.size:
            seen[new] = True
            n_seen += new.size
            if n_seen == reachable:
                return float(d), 100.0 * reachable / total_target_tiles
    return None, 100.0 * n_seen / total_target_tiles


def _run_objective(store, mask, decay):
    """Selection at budget MAX_SELECTION_STEPS (no gap-fill -- see header) with
    the time-aware objective at `decay` (decay 0 == coverage-first). Returns
    completion latency, final coverage % and final mean quality of the
    covered tiles."""
    sel, cov, _, _ = select_passes(
        store, MAX_SELECTION_STEPS, cloud_mask=mask, gap_fill=False,
        time_aware=decay > 0, latency_decay_rate=decay)
    lat = compute_completion_latency(sel, store, mask)
    q = float(np.mean(list(cov.values()))) if cov else 0.0
    return lat, 100.0 * len(cov) / total_target_tiles, q


def _effective_util(store, mask):
    """Per-row utility with cloud-obscured rows zeroed (same rule as Cell 6)."""
    return np.where(mask.flags, 0.0, store.util) if mask is not None else store.util


def _delivery_time(store, key):
    w = store.wait_lookup.get(key)
    return key[1] + (np.nan if w is None else w)


def select_fcfs(store, mask, k):
    """
    First-come-first-served baseline: take passes in order of delivery time
    (t_min + downlink wait), keeping any pass that adds at least one tile not
    yet covered, until k passes are kept. The obvious operational rule: no
    optimisation, just "use the next pass that shows something new".
    Returns (selected_keys, covered utility per tile index).
    """
    util = _effective_util(store, mask)
    best = np.zeros(len(_tile_ids))
    selected = []
    for key in sorted(store.keys_, key=delivery_sort_key(store)):
        if len(selected) >= k:
            break
        s = store.rows(store.index[key])
        t, u = store.tile_idx[s], util[s]
        if np.any((u > 0.0) & (best[t] == 0.0)):
            selected.append(key)
            np.maximum.at(best, t, u)
    return selected, best


def select_random(store, mask, k, rng):
    """Random baseline: k passes drawn uniformly (without replacement) from those
    that could cover at least one tile. With realistic cloud knowledge the draw
    uses clear-sky geometry (which passes turn out cloudy is not known in
    advance); coverage is then scored with the cloud mask."""
    util = _effective_util(store, mask)
    draw_util = (store.util if mask is not None and globals().get("CLOUD_KNOWLEDGE", "realistic") != "oracle"
                 else util)
    useful = [e for e in range(len(store)) if np.any(draw_util[store.rows(e)] > 0.0)]
    pick = rng.choice(len(useful), size=min(k, len(useful)), replace=False) if useful else []
    best = np.zeros(len(_tile_ids))
    selected = []
    for i in pick:
        e = useful[int(i)]
        s = store.rows(e)
        np.maximum.at(best, store.tile_idx[s], util[s])
        selected.append(store.keys_[e])
    return selected, best


def solve_ilp_max_utility(store, mask, k, time_limit_s=60.0, max_pairs=300_000, weight=None):
    """
    Exact budgeted maximum-utility selection (the objective Cell 6's greedy
    approximates), solved as a mixed-integer program with SciPy's HiGHS.
    weight (per event, v1.5.0): solve the time-weighted coverage instead,
    u[e,t] replaced by weight[e] * u[e,t] -- the exact optimum of the
    facility objective of the time-aware selector (Cell 0):
        maximise   sum_{e,t} u[e,t] * z[e,t]
        subject to sum_e z[e,t] <= 1          (a tile takes its utility once)
                   z[e,t] <= x[e]             (only from a selected pass)
                   sum_e x[e] <= k            (budget)
                   x binary, 0 <= z <= 1
    Given binary x, the optimal z simply assigns each tile to its best selected
    pass, so the objective equals the greedy's sum of best utilities per tile.
    Returns dict(status, objective, bound, selected, seconds, n_pairs).
    """
    import time as _time
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix, vstack
    util = _effective_util(store, mask)
    scale = 1.0
    if weight is not None:
        w = np.asarray(weight, dtype=np.float64)
        has = np.add.reduceat((util > 0.0).astype(np.int64), store.offsets[:-1]) > 0 if len(store) else np.zeros(0, bool)
        # v1.5.1: weights normalised to a maximum of 1 over the passes that can
        # contribute, so the solver's absolute tolerances (1e-6 to 1e-7) do not
        # swallow a small objective; the objective is scaled back on return.
        scale = float(w[has].max()) if has.any() and w[has].max() > 0 else 1.0
        util = util * np.repeat(w / scale, store.counts)
    keep = util > 0.0
    pair_event = np.repeat(np.arange(len(store)), store.counts)[keep]
    pair_tile = store.tile_idx[keep]
    pair_u = util[keep]
    n_e, n_p = len(store), int(keep.sum())
    if n_p == 0:
        return {"status": "empty", "objective": 0.0, "bound": 0.0, "selected": [], "seconds": 0.0, "n_pairs": 0}
    if n_p > max_pairs:
        return {"status": f"skipped (>{max_pairs} pairs)", "objective": None, "bound": None,
                "selected": [], "seconds": 0.0, "n_pairs": n_p}
    _, tile_row = np.unique(pair_tile, return_inverse=True)
    n_t = int(tile_row.max()) + 1
    c = np.concatenate([np.zeros(n_e), -pair_u])
    budget = coo_matrix((np.ones(n_e), (np.zeros(n_e), np.arange(n_e))), shape=(1, n_e + n_p))
    per_tile = coo_matrix((np.ones(n_p), (tile_row, n_e + np.arange(n_p))), shape=(n_t, n_e + n_p))
    link = coo_matrix((np.r_[np.ones(n_p), -np.ones(n_p)],
                       (np.r_[np.arange(n_p), np.arange(n_p)], np.r_[n_e + np.arange(n_p), pair_event])),
                      shape=(n_p, n_e + n_p))
    A = vstack([budget, per_tile, link]).tocsr()
    ub = np.r_[k, np.ones(n_t), np.zeros(n_p)]
    t0 = _time.time()
    # v1.5.1: relative MIP gap 1e-9 (HiGHS default 1e-4: "optimal" then meant
    # within 0.01% of the optimum)
    res = milp(c, constraints=LinearConstraint(A, -np.inf, ub), bounds=Bounds(0, 1),
               integrality=np.r_[np.ones(n_e), np.zeros(n_p)],
               options={"time_limit": float(time_limit_s), "mip_rel_gap": 1e-9})
    secs = _time.time() - t0
    if res.x is None:
        return {"status": f"no solution ({res.message})", "objective": None, "bound": None,
                "selected": [], "seconds": secs, "n_pairs": n_p}
    chosen = [store.keys_[e] for e in np.flatnonzero(res.x[:n_e] > 0.5)]
    bound = getattr(res, "mip_dual_bound", None)
    return {"status": "optimal" if res.status == 0 else "time limit (best found)",
            "objective": float(-res.fun) * scale, "bound": None if bound is None else float(-bound) * scale,
            "selected": chosen, "seconds": secs, "n_pairs": n_p}


def solve_ilp_min_images(store, mask, need, time_limit_s=60.0, max_pairs=300_000):
    """
    Exact minimum number of images covering at least `need` tiles (a tile is
    covered by any usable observation), as a mixed-integer program (HiGHS):
        minimise sum_e x[e]
        subject to y[t] <= sum_{e covers t} x[e],  sum_t y[t] >= need,
                   x, y binary.
    The check for ECC's greedy cover step. Returns dict(status, images, seconds, n_pairs).
    """
    import time as _time
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix, vstack
    util = _effective_util(store, mask)
    keep = util > 0.0
    pair_event = np.repeat(np.arange(len(store)), store.counts)[keep]
    pair_tile = store.tile_idx[keep]
    n_e, n_p = len(store), int(keep.sum())
    if n_p == 0:
        return {"status": "empty", "images": None, "seconds": 0.0, "n_pairs": 0}
    if n_p > max_pairs:
        return {"status": f"skipped (>{max_pairs} pairs)", "images": None, "seconds": 0.0, "n_pairs": n_p}
    _, tile_row = np.unique(pair_tile, return_inverse=True)
    n_t = int(tile_row.max()) + 1
    if n_t < need:
        return {"status": "infeasible", "images": None, "seconds": 0.0, "n_pairs": n_p}
    c = np.r_[np.ones(n_e), np.zeros(n_t)]
    cover = coo_matrix((np.r_[-np.ones(n_p), np.ones(n_t)],
                        (np.r_[tile_row, np.arange(n_t)], np.r_[pair_event, n_e + np.arange(n_t)])),
                       shape=(n_t, n_e + n_t))
    total = coo_matrix((np.ones(n_t), (np.zeros(n_t), n_e + np.arange(n_t))), shape=(1, n_e + n_t))
    t0 = _time.time()
    res = milp(c, constraints=[LinearConstraint(cover.tocsr(), -np.inf, 0.0), LinearConstraint(total.tocsr(), need, np.inf)],
               bounds=Bounds(0, 1), integrality=np.ones(n_e + n_t), options={"time_limit": float(time_limit_s)})
    secs = _time.time() - t0
    if res.x is None or res.status != 0:
        return {"status": "time limit or no solution", "images": None, "seconds": secs, "n_pairs": n_p}
    return {"status": "optimal", "images": int(round(res.fun)), "seconds": secs, "n_pairs": n_p}


def reachable_tile_count(store, mask):
    """Tiles with at least one usable observation (utility > 0, not
    cloud-obscured) among ALL events in this store: the most any policy
    could cover from these passes."""
    util = _effective_util(store, mask)
    return int(np.unique(store.tile_idx[util > 0.0]).size)


def coverage_times(store, mask, targets, selected=None, denominator=None):
    """
    Delivery time at which delivered coverage first reaches each target (% of
    ROI tiles), with events ordered by delivery time (t_min + downlink wait).
    selected=None uses ALL events -- the physical floor for that coverage
    level; otherwise only the selected events. A tile counts under the same
    rule as everywhere else (utility > 0, not cloud-obscured). Returns
    {target: minutes or None if never reached}. denominator: tile count the
    percentages refer to (default: all ROI tiles; the revised H1 passes the
    reachable-tile count).
    """
    util = _effective_util(store, mask)
    denom = total_target_tiles if denominator is None else denominator
    if denom <= 0:
        return {tg: None for tg in targets}
    keys = store.keys_ if selected is None else selected
    order = sorted(keys, key=delivery_sort_key(store))
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n_seen, out, pending = 0, {}, sorted(targets)
    for key in order:
        if not pending:
            break
        d = _delivery_time(store, key)
        if np.isnan(d):
            continue
        s = store.rows(store.index[key])
        t = store.tile_idx[s][util[s] > 0.0]
        new = t[~seen[t]]
        if new.size:
            seen[new] = True
            n_seen += new.size
            pct = 100.0 * n_seen / denom
            while pending and pct >= pending[0] - 1e-9:
                out[pending.pop(0)] = float(d)
    out.update({tg: None for tg in pending})
    return out


def _count_cover(store, util_eff, events, k, need):
    """Greedy maximum coverage, UNWEIGHTED (a tile counts once it has any usable
    observation: utility > 0, not cloud-obscured), using at most k of
    `events` (event indices, in delivery order -- ties go to the earlier
    image). Stops as soon as `need` tiles are covered.
    Returns (tiles covered, chosen event indices in greedy order)."""
    events = np.asarray(events, dtype=np.int64)
    if events.size == 0 or k <= 0:
        return 0, []
    counts = store.counts[events]
    rows = np.concatenate([np.arange(store.offsets[e], store.offsets[e + 1]) for e in events])
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    t = store.tile_idx[rows]
    ok = util_eff[rows] > 0.0
    covered = np.zeros(len(_tile_ids), dtype=bool)
    avail = np.ones(events.size, dtype=bool)
    chosen, n_cov = [], 0
    for _ in range(k):
        gain = np.add.reduceat((ok & ~covered[t]).astype(np.int64), starts)
        gain[~avail] = -1
        j = int(np.argmax(gain))
        if gain[j] <= 0:
            break
        avail[j] = False
        chosen.append(int(events[j]))
        sl = slice(int(starts[j]), int(starts[j] + counts[j]))
        tt = t[sl][ok[sl]]
        n_cov += int((~covered[tt]).sum())
        covered[tt] = True
        if n_cov >= need:
            break
    return n_cov, chosen


def earliest_completion(store, mask, k, target_pct, denominator=None):
    """
    Earliest-completion cover (ECC). The earliest delivery time t at which at
    most k images DELIVERED BY t cover target_pct of the tiles (of the ROI, or
    of `denominator` tiles), found by binary search over delivery times with
    a greedy set cover at each probe; the floor (all images) bounds the search
    from below. Uses only images already delivered at t, so under cloud it is
    causal: their cloud is known by then (from the image itself).
    Returns (completion time = latest delivery among the chosen images or
    None if never reached, chosen keys, images used). The greedy cover step
    is approximate (Chvatal 1979: at most H(n) times the minimum number of
    images), so the time is the earliest at which the GREEDY finds a cover.
    Bisection assumes that once the greedy succeeds with the images delivered
    by t it also succeeds with those delivered by any later t. Greedy maximum
    coverage is not monotone in its ground set in general, so this is an
    assumption (v1.5.2 note): it held at every onset of the October 2026 runs
    (bisection equal to a linear scan, test_ecc_monotone.py), and a violation
    could only make the reported time LATER than the greedy's earliest success
    -- every reported time is achieved by its returned images -- so H1a can
    only be understated by it, never overstated.
    """
    util_eff = _effective_util(store, mask)
    denom = total_target_tiles if denominator is None else denominator
    if denom <= 0 or len(store) == 0:
        return None, [], 0
    need = int(np.ceil(target_pct / 100.0 * denom - 1e-9))
    d = event_delivery_times(store)
    order = np.argsort(d, kind="stable")
    order = order[np.isfinite(d[order])]
    times = d[order]
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n_seen, lo = 0, None
    for pos, e in enumerate(order):                       # floor: all delivered images
        s = store.rows(e)
        tt = store.tile_idx[s][util_eff[s] > 0.0]
        new = tt[~seen[tt]]
        seen[new] = True
        n_seen += new.size
        if n_seen >= need:
            lo = pos
            break
    if lo is None:
        return None, [], 0

    def probe(pos):
        end = int(np.searchsorted(times, times[pos], side="right"))   # all images delivered by times[pos]
        n_cov, chosen = _count_cover(store, util_eff, order[:end], k, need)
        return n_cov >= need, chosen

    hi = len(order) - 1
    ok_hi, chosen_hi = probe(hi)
    if not ok_hi:
        return None, [], 0
    best = chosen_hi
    a, b = lo, hi
    while a < b:
        mid = (a + b) // 2
        ok_mid, ch = probe(mid)
        if ok_mid:
            b, best = mid, ch
        else:
            a = mid + 1
    # Invariant: probe(b) succeeded and `best` is its cover, so at exit (a == b)
    # `best` is the cover at the earliest successful delivery time found.
    return float(max(d[e] for e in best)), [store.keys_[e] for e in best], len(best)


def _ecc_best_effort(store, mask, k):
    """ECC when no target is reachable: the greedy k-image cover of every
    delivered image (maximum coverage by the end of the window)."""
    util_eff = _effective_util(store, mask)
    d = event_delivery_times(store)
    order = np.argsort(d, kind="stable")
    order = order[np.isfinite(d[order])]
    _, chosen = _count_cover(store, util_eff, order, k, need=10 ** 12)
    return [store.keys_[e] for e in chosen]


def level_stats(store, mask, selected, need):
    """The delivered mosaic of `selected` (images in delivery order, ties in
    selection order, as compute_completion_latency) when it first covers `need`
    tiles (usable: utility > 0, not cloud-obscured): (delivery time of the image
    at which that happens, images delivered up to and including it, mean
    utility of the tiles covered then). (None, None, None) if never."""
    if need <= 0 or not selected:
        return None, None, None
    util = _effective_util(store, mask)
    best = np.zeros(len(_tile_ids))
    n_cov = 0
    for i, key in enumerate(sorted(selected, key=delivery_sort_key(store)), 1):
        d = _delivery_time(store, key)
        if np.isnan(d):
            break                                   # undelivered images sort last
        s = store.rows(store.index[key])
        t, u = store.tile_idx[s], util[s]
        n_cov += int(np.unique(t[(u > 0.0) & (best[t] == 0.0)]).size)
        np.maximum.at(best, t, u)
        if n_cov >= need:
            cov = best > 0.0
            return float(d), i, float(best[cov].mean())
    return None, None, None


def mosaic_quality(store, mask, selected):
    """Mean utility of the tiles covered by `selected` (best usable observation per tile)."""
    util = _effective_util(store, mask)
    best = np.zeros(len(_tile_ids))
    for key in selected:
        s = store.rows(store.index[key])
        np.maximum.at(best, store.tile_idx[s], util[s])
    cov = best > 0.0
    return float(best[cov].mean()) if cov.any() else None


def facility_value(store, mask, selected, weight):
    """Time-weighted coverage f(S) = sum_t max_{e in S} w_e * u_e,t (usable observations)."""
    util = _effective_util(store, mask)
    best = np.zeros(len(_tile_ids))
    for key in selected:
        e = store.index[key]
        s = store.rows(e)
        np.maximum.at(best, store.tile_idx[s], weight[e] * util[s])
    return float(best.sum())


def realized_coverage(store, mask, selected):
    """Tiles with a usable observation (utility > 0, not cloud-obscured) among `selected`."""
    util = _effective_util(store, mask)
    seen = np.zeros(len(_tile_ids), dtype=bool)
    for key in selected:
        s = store.rows(store.index[key])
        seen[store.tile_idx[s][util[s] > 0.0]] = True
    return int(seen.sum())


def plan_expected_coverage(store, k, p_tile, weight=None):
    """
    Open-loop plan (v1.5.0): k passes chosen IN ADVANCE, with no cloud
    observation, by the same planner the adaptive selector uses
    (_expected_plan, Cell 6): expected coverage under cloud from pass geometry
    and each tile's cloud probability, passes treated as independent,
    time-weighted when `weight` is given. The adaptive selector differs only
    in re-planning from what it has seen, so the comparison isolates the value
    of acting on cloud observations.
    """
    if len(store) == 0 or k <= 0:
        return []
    d = event_delivery_times(store)
    plan = _expected_plan(store, np.ones(len(_tile_ids), dtype=bool), np.isfinite(d), int(k), p_tile, weight)
    return [store.keys_[e] for e in plan]


def paired_summary(d, alpha=0.05):
    """Mean of paired differences with its one-sided (1 - alpha) lower confidence bound
    (Student t, n - 1 degrees of freedom), SD, and the share of positive / negative pairs."""
    d = np.asarray([x for x in d if x is not None and np.isfinite(x)], dtype=np.float64)
    n = int(d.size)
    if n == 0:
        return {"n": 0, "mean": None, "sd": None, "lcb": None, "pct_positive": None, "pct_negative": None}
    m = float(d.mean())
    sd = float(d.std(ddof=1)) if n > 1 else 0.0
    if n > 1:
        from scipy.stats import t as _t
        lcb = m - float(_t.ppf(1.0 - alpha, n - 1)) * sd / math.sqrt(n)
    else:
        lcb = None
    return {"n": n, "mean": m, "sd": sd, "lcb": lcb, "alpha": alpha,
            "pct_positive": float(100.0 * (d > 1e-9).mean()), "pct_negative": float(100.0 * (d < -1e-9).mean())}


def _metrics(selected, store, mask):
    """Completion latency, coverage %, mean quality and per-target coverage
    times of one selection (budget-limited, no gap-fill)."""
    util = _effective_util(store, mask)
    best = np.zeros(len(_tile_ids))
    for key in selected:
        s = store.rows(store.index[key])
        np.maximum.at(best, store.tile_idx[s], util[s])
    covered = best > 0.0
    return {"lat": compute_completion_latency(selected, store, mask),
            "cov": 100.0 * int(covered.sum()) / total_target_tiles,
            "q": float(best[covered].mean()) if covered.any() else 0.0,
            "t": coverage_times(store, mask, COVERAGE_TARGETS_PCT, selected)}


POLICIES = ("coverage_first", "time_aware", "fcfs", "random", "earliest_completion")


def _policy_selections(store, mask, policy, seed=None):
    """One selection (deterministic policies) or RANDOM_BASELINE_DRAWS seeded
    selections (random), each at budget MAX_SELECTION_STEPS, no gap-fill."""
    k = MAX_SELECTION_STEPS
    if policy in ("coverage_first", "time_aware"):
        rate = 0.0 if policy == "coverage_first" else LATENCY_DECAY_RATE
        sel, _, _, _ = select_passes(store, k, cloud_mask=mask, gap_fill=False,
                                                       time_aware=rate > 0, latency_decay_rate=rate)
        return [sel]
    if policy == "fcfs":
        return [select_fcfs(store, mask, k)[0]]
    rng = np.random.default_rng(RANDOM_SEED if seed is None else seed)
    return [select_random(store, mask, k, rng)[0] for _ in range(RANDOM_BASELINE_DRAWS)]


def _median_or_none(values):
    """Median treating None as 'never' (+inf); None if the median is never."""
    a = np.array([np.inf if v is None else v for v in values], dtype=float)
    m = float(np.median(a)) if a.size else np.inf
    return None if np.isinf(m) else m


def _policy_outcome(store, mask, policy, seed=None):
    """Metrics of one policy: completion latency, coverage %, mean quality and
    time to each coverage target. Random: median over its seeded draws."""
    if len(store) == 0:
        return {"lat": None, "cov": 0.0, "q": 0.0, "t": {tg: None for tg in COVERAGE_TARGETS_PCT}}
    if policy == "earliest_completion":
        # One ECC run per coverage target (it stops once its target is met);
        # coverage/quality/latency are those of the highest target reached,
        # or of the best k-image cover if none is reached.
        t, best_sel = {}, None
        for tg in sorted(COVERAGE_TARGETS_PCT):
            tc, sel, _ = earliest_completion(store, mask, MAX_SELECTION_STEPS, tg)
            t[tg] = tc
            if sel:
                best_sel = sel
        m = _metrics(best_sel if best_sel is not None else _ecc_best_effort(store, mask, MAX_SELECTION_STEPS),
                     store, mask)
        m["t"] = t
        return m
    ms = [_metrics(sel, store, mask) for sel in _policy_selections(store, mask, policy, seed)]
    if len(ms) == 1:
        return ms[0]
    return {"lat": _median_or_none([m["lat"] for m in ms]),
            "cov": float(np.median([m["cov"] for m in ms])),
            "q": float(np.median([m["q"] for m in ms])),
            "t": {tg: _median_or_none([m["t"][tg] for m in ms]) for tg in COVERAGE_TARGETS_PCT}}


# =============================================================================
# 1. DECAY SWEEP
# =============================================================================
_sweep_rates = [0.0] + sorted(set(float(x) for x in DECAY_SWEEP) - {0.0})
decay_sweep_rows = []
for rate in _sweep_rates:
    lc, cc, qc = _run_objective(event_data, None, rate)
    ll, cl, ql = _run_objective(event_data, cloud_mask, rate)
    decay_sweep_rows.append({
        "decay_rate": rate, "objective": "coverage-first" if rate == 0 else "time-aware",
        "h1_latency_clear_min": lc, "coverage_clear_pct": cc, "mean_quality_clear": qc,
        "h1_latency_cloud_min": ll, "coverage_cloud_pct": cl, "mean_quality_cloud": ql,
    })
decay_sweep_df = pd.DataFrame(decay_sweep_rows)
decay_sweep_df.to_csv(os.path.join(_out_dir, "arf_decay_sweep.csv"), index=False)
logger.info(f"Decay sweep over {len(_sweep_rates)} rates saved to arf_decay_sweep.csv")

# =============================================================================
# 1B. POLICY COMPARISON (full window, budget k, no gap-fill)
# =============================================================================
policy_rows = []
for sky, m in (("clear", None), ("cloud", cloud_mask)):
    floor = coverage_times(event_data, m, COVERAGE_TARGETS_PCT)
    policy_rows.append({"sky": sky, "policy": "best_possible", "h1_latency_min": None, "coverage_pct": None,
                        "mean_quality": None, **{f"t{tg:g}_min": floor[tg] for tg in COVERAGE_TARGETS_PCT}})
    for name in POLICIES:
        o = _policy_outcome(event_data, m, name)
        policy_rows.append({"sky": sky, "policy": name, "h1_latency_min": o["lat"],
                            "coverage_pct": o["cov"], "mean_quality": o["q"],
                            **{f"t{tg:g}_min": o["t"][tg] for tg in COVERAGE_TARGETS_PCT}})
policy_df = pd.DataFrame(policy_rows)
policy_df.to_csv(os.path.join(_out_dir, "arf_policy_comparison.csv"), index=False)

# =============================================================================
# 1C. COST-LATENCY FRONTIER (earliest-completion cover at each image budget)
# =============================================================================
# For a budget of k images: the earliest time a mosaic of at most k delivered
# images reaches the target. k -> infinity is the floor. Targets: the revised-
# H1 level of REACHABLE tiles, and each COVERAGE_TARGETS_PCT level of the ROI.
# Grid up to the number of candidate images, plus the run's own budget k always
# (a budget above the image count behaves like the image count).
_frontier_k = sorted({k for k in (1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
                      if k <= len(event_data)} | {MAX_SELECTION_STEPS})
frontier_rows = []
_h1r_level = float(globals().get("H1_REACHABLE_COVERAGE_PCT", 95.0))
for sky, m in (("clear", None), ("cloud", cloud_mask)):
    _reach = reachable_tile_count(event_data, m)
    targets = [("reachable", _h1r_level, _reach)] + [("roi", tg, total_target_tiles) for tg in COVERAGE_TARGETS_PCT]
    for kind, tg, denom in targets:
        floor_t = coverage_times(event_data, m, (tg,), denominator=denom)[tg] if denom else None
        for k in _frontier_k:
            tc, _, n_used = earliest_completion(event_data, m, k, tg, denominator=denom) if denom else (None, [], 0)
            frontier_rows.append({"sky": sky, "target": kind, "target_pct": tg, "denominator_tiles": denom,
                                  "k": k, "t_complete_min": tc, "images_used": n_used if tc is not None else None,
                                  "floor_min": floor_t})
frontier_df = pd.DataFrame(frontier_rows)
frontier_df.to_csv(os.path.join(_out_dir, "arf_cost_latency_frontier.csv"), index=False)
logger.info(f"Cost-latency frontier over k={_frontier_k} saved to arf_cost_latency_frontier.csv")

# v1.5.2: two IEEE columns wide (7.16 in) with 8-9 pt text (was 7.5 x 4.6 in).
with plt.rc_context({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "ps.fonttype": 42}):
    fig6, ax6 = plt.subplots(figsize=(7.16, 3.0))
    for sky, col in (("clear", "#1b9e77"), ("cloud", "#7570b3")):
        f = frontier_df[(frontier_df.sky == sky) & (frontier_df.target == "reachable")]
        ok = f.t_complete_min.notna()
        if ok.any():
            ax6.plot(f.k[ok], f.t_complete_min[ok] / 60.0, marker="o", ms=3.5, lw=1.2, color=col,
                     label=f"{'Clear sky' if sky == 'clear' else 'Cloud (realised)'}: ECC")
            fl = f.floor_min.dropna()
            if len(fl):
                ax6.axhline(fl.iloc[0] / 60.0, color=col, ls=":", lw=1.0,
                            label=f"{'Clear' if sky == 'clear' else 'Cloud'} floor (all images)")
    ax6.set_xscale("log")
    ax6.set_xlabel("Image budget k (images downloaded / processed / purchased)")
    ax6.set_ylabel(f"Time to {_h1r_level:g}% of reachable tiles (h)")
    _scen6 = globals().get("SCENARIO_N", ACTIVE_N)
    ax6.set_title(f"Cost-latency frontier: {roi_short_name() if 'roi_short_name' in globals() else ACTIVE_ROI}, "
                  f"{'N=' + str(_scen6) if isinstance(_scen6, int) else _scen6}", loc="left")
    ax6.grid(alpha=0.3, which="both", lw=0.5)
    ax6.legend()
    fig6.tight_layout()
    for fmt in ("png", "pdf"):
        fig6.savefig(os.path.join(_out_dir, f"arf_fig6_cost_latency.{fmt}"), dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig6)

# =============================================================================
# 2. ONSET ROBUSTNESS
# =============================================================================
onset_robustness_summary = {}
if ONSET_ROBUSTNESS_ENABLED:
    _horizon = (WINDOW_HOURS - ONSET_SHIFT_HOURS) * 60
    _shifts = list(range(0, ONSET_SHIFT_HOURS * 60 + 1, ONSET_SHIFT_STEP_MIN))
    logger.info(f"Onset robustness: {len(_shifts)} onsets over {ONSET_SHIFT_HOURS} h, "
                f"{_horizon / 60:.0f} h of passes each, decay={LATENCY_DECAY_RATE}")
    onset_rows = []
    for d in _shifts:
        sub, sub_mask = subset_by_time(event_data, cloud_mask, d, d + _horizon)
        row = {"onset_shift_min": d, "n_events": len(sub),
               "decay_tau_min": decay_time_scale(sub) if len(sub) else None}   # time-aware weight scale here
        for sky, m in (("clear", None), ("cloud", sub_mask)):
            row[f"best_possible_{sky}_min"], row[f"best_possible_{sky}_cov_pct"] = best_possible_completion(sub, m)
            for tg, v in coverage_times(sub, m, COVERAGE_TARGETS_PCT).items():
                row[f"best_possible_{sky}_t{tg:g}_min"] = v
            for name in POLICIES:
                o = _policy_outcome(sub, m, name, seed=RANDOM_SEED + d)
                row[f"{name}_{sky}_min"], row[f"{name}_{sky}_cov_pct"] = o["lat"], o["cov"]
                for tg, v in o["t"].items():
                    row[f"{name}_{sky}_t{tg:g}_min"] = v
            # Revised H1: time to H1_REACHABLE_COVERAGE_PCT of the tiles reachable
            # in this onset's horizon, floor vs each deterministic policy.
            n_reach = reachable_tile_count(sub, m)
            row[f"reachable_{sky}_tiles"] = n_reach
            row[f"best_possible_{sky}_t{_H1R:g}r_min"] = coverage_times(sub, m, (_H1R,), denominator=n_reach)[_H1R]
            # H1a/H1b (v1.4.0): time, images and mosaic quality when each
            # policy's delivered mosaic first covers `need` reachable tiles.
            need = int(np.ceil(_H1R / 100.0 * n_reach - 1e-9)) if n_reach else 0
            for name in ("coverage_first", "time_aware", "fcfs"):
                sel = _policy_selections(sub, m, name)[0] if len(sub) else []
                row[f"{name}_{sky}_t{_H1R:g}r_min"] = (coverage_times(sub, m, (_H1R,), sel, denominator=n_reach)[_H1R]
                                                       if n_reach else None)
                if name != "fcfs":
                    _, n_img, q = level_stats(sub, m, sel, need) if n_reach else (None, None, None)
                    row[f"{name}_{sky}_images{_H1R:g}r"], row[f"{name}_{sky}_quality{_H1R:g}r"] = n_img, q
            # FCFS without a budget: the images it takes to reach the level (at the floor).
            t_f, n_f, q_f = (level_stats(sub, m, select_fcfs(sub, m, len(sub))[0], need) if n_reach
                             else (None, None, None))
            row[f"fcfs_unbudgeted_{sky}_t{_H1R:g}r_min"] = t_f
            row[f"fcfs_unbudgeted_{sky}_images{_H1R:g}r"], row[f"fcfs_unbudgeted_{sky}_quality{_H1R:g}r"] = n_f, q_f
            t_e, sel_e, n_e = (earliest_completion(sub, m, MAX_SELECTION_STEPS, _H1R, denominator=n_reach) if n_reach
                               else (None, [], 0))
            row[f"earliest_completion_{sky}_t{_H1R:g}r_min"] = t_e
            row[f"earliest_completion_{sky}_images{_H1R:g}r"] = n_e if t_e is not None else None
            row[f"earliest_completion_{sky}_quality{_H1R:g}r"] = mosaic_quality(sub, m, sel_e) if t_e is not None else None
        onset_rows.append(row)
    onset_df = pd.DataFrame(onset_rows)
    onset_df.to_csv(os.path.join(_out_dir, "arf_onset_robustness.csv"), index=False)

    def _stats(col):
        x = pd.to_numeric(onset_df[col], errors="coerce").to_numpy(dtype=float)
        ok = ~np.isnan(x)
        return {
            "median_min": float(np.median(x[ok])) if ok.any() else None,
            "p90_min": float(np.percentile(x[ok], 90)) if ok.any() else None,
            "max_min": float(np.max(x[ok])) if ok.any() else None,
            # Onsets where the series never completes count as NOT meeting H1.
            "pct_onsets_under_h1": float(100.0 * np.mean(ok & (np.nan_to_num(x, nan=np.inf) < _H1_TARGET_MIN))),
            "n_onsets": int(len(x)), "n_completed": int(ok.sum()),
            "pct_onsets_reached": float(100.0 * ok.mean()) if len(x) else 0.0,
        }

    _taus = pd.to_numeric(onset_df["decay_tau_min"], errors="coerce").dropna()
    onset_robustness_summary = {
        "h1_target_min": _H1_TARGET_MIN, "latency_decay_rate": LATENCY_DECAY_RATE,
        "latency_decay_scale": globals().get("LATENCY_DECAY_SCALE", "fleet"),
        "decay_tau_min_over_onsets": ({"median": float(_taus.median()), "min": float(_taus.min()),
                                       "max": float(_taus.max())} if len(_taus) else None),
        "onset_shift_hours": ONSET_SHIFT_HOURS, "onset_shift_step_min": ONSET_SHIFT_STEP_MIN,
        "horizon_hours": _horizon / 60,
        "series": {f"{name}_{sky}": _stats(f"{name}_{sky}_min")
                   for sky in ("clear", "cloud") for name in ("best_possible",) + POLICIES},
        # Equal-coverage series: time to each common coverage level.
        "coverage_targets_pct": list(COVERAGE_TARGETS_PCT),
        "equal_coverage_series": {f"{name}_{sky}_t{tg:g}": _stats(f"{name}_{sky}_t{tg:g}_min")
                                  for tg in COVERAGE_TARGETS_PCT for sky in ("clear", "cloud")
                                  for name in ("best_possible",) + POLICIES},
    }
    # Revised H1 per policy and sky: share of onsets (with any reachable tile)
    # where the policy reaches _H1R% of reachable tiles within the margin of
    # the floor. An onset where the policy never reaches it counts as a miss.
    def _h1_revised(name, sky):
        n_ok = n_eval = 0
        gaps = []
        for _, r in onset_df.iterrows():
            if not r.get(f"reachable_{sky}_tiles"):
                continue
            fl = r.get(f"best_possible_{sky}_t{_H1R:g}r_min")
            if fl is None or pd.isna(fl):
                continue
            n_eval += 1
            v = r.get(f"{name}_{sky}_t{_H1R:g}r_min")
            if v is not None and not pd.isna(v):
                gaps.append(float(v) - float(fl))
                n_ok += int(float(v) - float(fl) <= H1_FLOOR_MARGIN_MIN + 1e-9)
        lcb = cp_lower_bound(n_ok, n_eval)
        return {"n_onsets_evaluated": n_eval, "n_within_margin": n_ok,
                "pct_within_margin": (100.0 * n_ok / n_eval) if n_eval else None,
                "pct_within_margin_lcb95": None if lcb is None else 100.0 * lcb,
                "evaluable": n_eval >= H1_MIN_EVALUABLE_ONSETS,
                "median_gap_min": float(np.median(gaps)) if gaps else None,
                "p90_gap_min": float(np.percentile(gaps, 90)) if gaps else None,
                "pct_reached": (100.0 * len(gaps) / n_eval) if n_eval else None}

    _evaluability = {"min_evaluable_onsets": H1_MIN_EVALUABLE_ONSETS, "source": H1_EVALUABILITY_SOURCE,
                     "confidence_alpha": H1_CONFIDENCE_ALPHA,
                     "reference": "Clopper, C. J. & Pearson, E. S. (1934). Biometrika 26(4), 404-413, "
                                  "doi:10.1093/biomet/26.4.404",
                     "caveat": "consecutive onsets share most of their passes, so they are not independent "
                               "trials: the lower bound is optimistic and the minimum count necessary, not sufficient"}
    onset_robustness_summary["h1_revised"] = {
        "framing": ("H1a: latency non-inferiority to the physical floor. Pre-specified 2026-10-09 for the "
                    "time-aware selector; from v1.5.0 tested on the earliest-completion cover (ECC), the selector "
                    "H1b also tests, with the time-aware selector reported (series time_aware_*)."),
        "primary_series": "earliest_completion",
        # v1.5.1: one window's onsets share most of their passes, so they are not
        # independent trials; within a run the status and its bound describe this
        # window. The inferential test treats windows as the unit (confirmatory
        # protocol, arf_confirmatory_analysis.py).
        "inference": ("descriptive for this window: consecutive onsets share passes, so the Clopper-Pearson bound, "
                      "which assumes independent trials, overstates the evidence; test across independent windows"),
        "statement": (f"At >= {_H1MIN:g}% of onsets the earliest-completion cover reaches {_H1R:g}% of "
                      f"reachable tiles within {H1_FLOOR_MARGIN_MIN:.1f} min of the physical floor, "
                      f"using at most {MAX_SELECTION_STEPS} images (clear sky)."),
        "reachable_coverage_pct": _H1R, "floor_margin_min": H1_FLOOR_MARGIN_MIN,
        "floor_margin_source": H1_MARGIN_SOURCE, "min_onset_pct": _H1MIN,
        "budget_k": MAX_SELECTION_STEPS, **_evaluability,
        "series": {f"{name}_{sky}": _h1_revised(name, sky)
                   for sky in ("clear", "cloud")
                   for name in ("time_aware", "coverage_first", "fcfs", "earliest_completion")},
    }

    def _num(r, col):
        v = r.get(col)
        return None if v is None or pd.isna(v) else float(v)

    # H1b: image economy of ECC over FCFS at the level, ECC within the margin.
    def _h1b(sky):
        n_eval = n_ok = n_one = 0
        n_e_all, n_f_all, ratios = [], [], []
        for _, r in onset_df.iterrows():
            fl = _num(r, f"best_possible_{sky}_t{_H1R:g}r_min")
            n_f = _num(r, f"fcfs_unbudgeted_{sky}_images{_H1R:g}r")
            if not r.get(f"reachable_{sky}_tiles") or fl is None or n_f is None:
                continue
            if n_f < 2:
                n_one += 1                              # one image: no room for economy
                continue
            n_eval += 1
            t_e = _num(r, f"earliest_completion_{sky}_t{_H1R:g}r_min")
            n_e = _num(r, f"earliest_completion_{sky}_images{_H1R:g}r")
            n_f_all.append(n_f)
            if t_e is not None and n_e is not None:
                n_e_all.append(n_e)
                ratios.append(n_e / n_f)
                n_ok += int(t_e - fl <= H1_FLOOR_MARGIN_MIN + 1e-9 and n_e < n_f)
        lcb = cp_lower_bound(n_ok, n_eval)
        return {"n_onsets_evaluated": n_eval, "n_onsets_one_image": n_one, "n_success": n_ok,
                "pct_success": (100.0 * n_ok / n_eval) if n_eval else None,
                "pct_success_lcb95": None if lcb is None else 100.0 * lcb,
                "evaluable": n_eval >= H1_MIN_EVALUABLE_ONSETS,
                "median_images_ecc": float(np.median(n_e_all)) if n_e_all else None,
                "median_images_fcfs": float(np.median(n_f_all)) if n_f_all else None,
                "median_image_ratio_ecc_to_fcfs": float(np.median(ratios)) if ratios else None,
                "max_images_ecc": float(np.max(n_e_all)) if n_e_all else None}

    onset_robustness_summary["h1b_image_economy"] = {
        "framing": "H1b: image economy (added 2026-10-12, after the 15 runs of that day; confirmatory on fresh data only)",
        # v1.5.1: one window's onsets share most of their passes, so they are not
        # independent trials; within a run the status and its bound describe this
        # window. The inferential test treats windows as the unit (confirmatory
        # protocol, arf_confirmatory_analysis.py).
        "inference": ("descriptive for this window: consecutive onsets share passes, so the Clopper-Pearson bound, "
                      "which assumes independent trials, overstates the evidence; test across independent windows"),

        "statement": (f"At >= {_H1MIN:g}% of onsets where first-come-first-served needs two or more images to "
                      f"reach {_H1R:g}% of reachable tiles, the earliest-completion cover reaches that level within "
                      f"{H1_FLOOR_MARGIN_MIN:.1f} min of the physical floor using fewer images than first-come-"
                      f"first-served, with at most {MAX_SELECTION_STEPS} images (clear sky)."),
        "reachable_coverage_pct": _H1R, "floor_margin_min": H1_FLOOR_MARGIN_MIN, "min_onset_pct": _H1MIN,
        "budget_k": MAX_SELECTION_STEPS, **_evaluability,
        "series": {sky: _h1b(sky) for sky in ("clear", "cloud")},
    }

    # Mosaic quality when each policy first reaches the level (reported, not tested).
    def _quality(sky):
        out = {}
        for name in ("time_aware", "coverage_first", "fcfs_unbudgeted", "earliest_completion"):
            q = pd.to_numeric(onset_df[f"{name}_{sky}_quality{_H1R:g}r"], errors="coerce").dropna()
            out[name] = {"median": float(q.median()) if len(q) else None, "n_onsets": int(len(q))}
        for name in ("time_aware", "earliest_completion"):
            a = pd.to_numeric(onset_df[f"{name}_{sky}_quality{_H1R:g}r"], errors="coerce")
            b = pd.to_numeric(onset_df[f"fcfs_unbudgeted_{sky}_quality{_H1R:g}r"], errors="coerce")
            ok = a.notna() & b.notna() & (b > 0)
            out[f"{name}_vs_fcfs"] = {
                "n_onsets": int(ok.sum()),
                "median_difference": float((a[ok] - b[ok]).median()) if ok.any() else None,
                "median_relative_difference_pct": float((100.0 * (a[ok] - b[ok]) / b[ok]).median()) if ok.any() else None,
                "pct_onsets_higher": float(100.0 * (a[ok] > b[ok] + 1e-12).mean()) if ok.any() else None}
        return out

    onset_robustness_summary["quality_at_level"] = {
        "note": (f"Mean utility of the tiles covered when each policy's delivered mosaic first reaches {_H1R:g}% "
                 f"of reachable tiles; first-come-first-served without a budget. Reported, not tested."),
        "series": {sky: _quality(sky) for sky in ("clear", "cloud")},
    }
    with open(os.path.join(_out_dir, "arf_onset_robustness_summary.json"), "w") as f:
        json.dump(onset_robustness_summary, f, indent=2)

    # Figure 5: empirical CDF, across the evaluated onsets, of the time from the
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
            for name, (lab, col, ls, lw) in [("best_possible", (_floor_lab, "0.72", "-", 4.0))] + \
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
            _notes = (["At the floor at every onset:"] + at_floor if at_floor else []) + \
                     (["Never at the level:"] + never if never else [])
            if _notes:
                ax.text(0.98, 0.04, "\n".join(_notes), transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
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
                ax.text(0.97 * H1_FLOOR_MARGIN_MIN, 0.04, "Gap to the floor, per selector\n" + "\n".join(notes),
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
else:
    logger.info("Onset robustness analysis disabled (onset_robustness_enabled=False in Cell 0).")

# =============================================================================
# 3. OPTIMALITY GAP (exact ILP on small sub-problems, clear sky)
# =============================================================================
ilp_rows, ilp_summary = [], {}
_ILP_MIN_HORIZON_MIN = 30
if ILP_ENABLED:
    _span = max(0, WINDOW_HOURS - ILP_HORIZON_HOURS) * 60
    for d in sorted(set(int(x) for x in np.linspace(0, _span, ILP_INSTANCES))):
        # Large ROIs: halve the window until the instance fits ILP_MAX_PAIRS
        # (down to _ILP_MIN_HORIZON_MIN), so the gap stays measurable; the
        # window actually used is recorded per instance.
        h_min = ILP_HORIZON_HOURS * 60
        sub, _ = subset_by_time(event_data, None, d, d + h_min)
        while int((sub.util > 0.0).sum()) > ILP_MAX_PAIRS and h_min > _ILP_MIN_HORIZON_MIN:
            h_min = max(_ILP_MIN_HORIZON_MIN, h_min // 2)
            sub, _ = subset_by_time(event_data, None, d, d + h_min)
        if len(sub) == 0:
            continue
        res = solve_ilp_max_utility(sub, None, ILP_BUDGET_K, ILP_TIME_LIMIT_S, ILP_MAX_PAIRS)
        row = {"shift_min": d, "horizon_min": h_min, "n_events": len(sub), "n_pairs": res["n_pairs"],
               "ilp_status": res["status"],
               "ilp_seconds": round(res["seconds"], 2), "ilp_objective": res["objective"], "ilp_bound": res["bound"],
               "ilp_latency_min": compute_completion_latency(res["selected"], sub) if res["selected"] else None}
        _sel_ta = []
        for name, rate in (("coverage_first", 0.0), ("time_aware", LATENCY_DECAY_RATE)):
            _sel, cov, _, _ = run_submodular_greedy_selection(sub, ILP_BUDGET_K, gap_fill=False,
                                                              time_aware=rate > 0, latency_decay_rate=rate)
            row[f"{name}_objective"] = float(sum(cov.values()))
            if rate > 0:
                _sel_ta = _sel
        row["fcfs_objective"] = float(select_fcfs(sub, None, ILP_BUDGET_K)[1].sum())
        # v1.5.0: the time-aware selector against the exact optimum of ITS OWN
        # objective, the time-weighted coverage (ILP with u replaced by w * u).
        _w = time_aware_weights(sub, LATENCY_DECAY_RATE)
        _rw = solve_ilp_max_utility(sub, None, ILP_BUDGET_K, ILP_TIME_LIMIT_S, ILP_MAX_PAIRS, weight=_w)
        row.update({"tw_ilp_status": _rw["status"], "tw_ilp_objective": _rw["objective"],
                    "time_aware_tw_objective": facility_value(sub, None, _sel_ta, _w)})
        row["time_aware_tw_ratio"] = (row["time_aware_tw_objective"] / _rw["objective"]
                                      if _rw["status"] == "optimal" and _rw["objective"] else None)
        # Ratios only against a PROVEN optimum. A time-limited solve's best-found
        # solution can be arbitrarily poor and its dual bound loose, so such
        # instances are recorded but excluded from the gap statistics.
        _optimal = res["status"] == "optimal" and res["objective"]
        for name in ("coverage_first", "time_aware", "fcfs"):
            row[f"{name}_ratio"] = row[f"{name}_objective"] / res["objective"] if _optimal else None
        # ECC's greedy cover step vs the exact minimum number of images, for the
        # revised-H1 coverage level of this window's reachable tiles.
        _lvl = float(globals().get("H1_REACHABLE_COVERAGE_PCT", 95.0))
        _R = reachable_tile_count(sub, None)
        _need = int(np.ceil(_lvl / 100.0 * _R - 1e-9))
        _tg, _, _ng = earliest_completion(sub, None, len(sub), _lvl, denominator=_R) if _R else (None, [], 0)
        _ex = solve_ilp_min_images(sub, None, _need, ILP_TIME_LIMIT_S, ILP_MAX_PAIRS) if _R else {"status": "empty", "images": None}
        row.update({"ecc_target_pct": _lvl, "ecc_greedy_images": _ng if _tg is not None else None,
                    "ecc_exact_images": _ex["images"], "ecc_exact_status": _ex["status"]})
        ilp_rows.append(row)
    ilp_df = pd.DataFrame(ilp_rows)
    ilp_df.to_csv(os.path.join(_out_dir, "arf_ilp_gap.csv"), index=False)
    _solved = ilp_df[ilp_df["ilp_status"] == "optimal"] if len(ilp_df) else ilp_df
    ilp_summary = {
        "instances": int(len(ilp_df)), "solved": int(len(_solved)),
        "optimal": int((ilp_df["ilp_status"] == "optimal").sum()) if len(ilp_df) else 0,
        "time_limited": int(ilp_df["ilp_status"].astype(str).str.startswith("time limit").sum()) if len(ilp_df) else 0,
        "horizon_hours": ILP_HORIZON_HOURS, "budget_k": ILP_BUDGET_K,
        "horizons_used_min": sorted(int(h) for h in ilp_df["horizon_min"].unique()) if len(ilp_df) else [],
        **{f"{n}_ratio_{s}": (float(getattr(_solved[f"{n}_ratio"], s)()) if len(_solved) else None)
           for n in ("coverage_first", "time_aware", "fcfs") for s in ("median", "min")},
        # v1.5.0: time-aware greedy / exact optimum of the time-weighted coverage
        "time_aware_objective": globals().get("TIME_AWARE_OBJECTIVE", "facility"),
        "tw_solved": int((ilp_df["tw_ilp_status"] == "optimal").sum()) if len(ilp_df) else 0,
        **{f"time_aware_tw_ratio_{s}": (float(getattr(ilp_df["time_aware_tw_ratio"].dropna(), s)())
                                        if len(ilp_df) and ilp_df["time_aware_tw_ratio"].notna().any() else None)
           for s in ("median", "min")},
        "guarantee_note": ("facility objective: the k-step greedy attains >= 1 - 1/e of the optimum of the time-"
                           "weighted coverage (Nemhauser, Wolsey & Fisher 1978); the ratio above is the measured value"
                           if globals().get("TIME_AWARE_OBJECTIVE", "facility") == "facility" else
                           "legacy objective: no guarantee; the ratio is measured against the time-weighted coverage "
                           "optimum for comparison only"),
        **(lambda _e: {"ecc_instances_exact": int(len(_e)),
                       "ecc_matches_exact": int((_e["ecc_greedy_images"] == _e["ecc_exact_images"]).sum()),
                       "ecc_max_extra_images": (int((_e["ecc_greedy_images"] - _e["ecc_exact_images"]).max())
                                                if len(_e) else None)})(
            ilp_df[(ilp_df.get("ecc_exact_status") == "optimal") & ilp_df["ecc_greedy_images"].notna()]
            if len(ilp_df) and "ecc_exact_status" in ilp_df else pd.DataFrame(columns=["ecc_greedy_images", "ecc_exact_images"])),
    }

# =============================================================================
# 4. CLOUD ROBUSTNESS OVER CLOUD REALISATIONS (H2, v1.5.0)
# =============================================================================
# One cloud realisation is one random draw, so a comparison of policies under
# cloud needs many. For CLOUD_SEEDS realisations (seeds RANDOM_SEED,
# RANDOM_SEED + 1, ...; the model of Cell 0), full window, for each greedy
# objective (time-aware -- the primary -- and coverage-first), four plans with
# the SAME number of images n:
#   adaptive          the Fail-Signal selector as configured (select_passes,
#                     CLOUD_PLANNING, gap-fill on): learns each image's cloud
#                     when it arrives and re-plans; n = the images it uses
#   open-loop (exp.)  n images chosen in advance by the SAME planner without
#                     any cloud observation (plan_expected_coverage, expected
#                     coverage from the tiles' cloud probabilities, time-
#                     weighted for the time-aware objective)
#   open-loop (clear) n images chosen in advance by the clear-sky objective
#                     (no cloud information at all)
#   oracle            n images chosen by the greedy knowing every image's cloud
#                     (clairvoyant reference; greedy, so not a strict upper bound)
# and, for the record, the v1.4.0 Fail-Signal (clear-sky planning, re-plan on
# failure) with its own image count and an expected-coverage open-loop plan of
# that size.
# H2 (v1.5.0): the mean paired difference in coverage, adaptive minus
# open-loop (expected), has a one-sided 95% lower confidence bound above zero
# (Student t over the realisations) -- the value of acting on what the images
# show, the planner held fixed. NOT EVALUABLE when the selector never saw a
# failure or compensated in any realisation, when even the oracle never
# beats the open-loop plan (no room for adaptation), or under
# CLOUD_KNOWLEDGE="oracle". Each realisation also records the floor, ECC and
# FCFS under cloud at the revised-H1 level (H1a/H1b under cloud, secondary).
cloud_seed_rows, cloud_robustness_summary = [], {}
cloud_sensitivity_rows = []
_cs_n = int(globals().get("CLOUD_SEEDS", 20))
_ptile = globals().get("CLOUD_P_TILE")
if _ptile is None:
    _ptile = np.full(len(_tile_ids), float(CLOUD_PROBABILITY))
_OBJECTIVES = (("time_aware", float(LATENCY_DECAY_RATE)), ("coverage_first", 0.0))
_PLANS = ("adaptive", "openloop_expected", "openloop_clear", "oracle")
_PLAN_LABELS = {"adaptive": "Fail-Signal (adaptive)", "openloop_expected": "Open-loop (expected coverage)",
                "openloop_clear": "Open-loop (clear-sky plan)", "oracle": "Greedy with the cloud known"}


def _cloud_trial(m, objectives=None, with_h1=True):
    """One cloud realisation: the four plans for each objective (default
    _OBJECTIVES), plus the floor, ECC and FCFS at the revised-H1 level of the
    tiles reachable under that cloud."""
    objectives = _OBJECTIVES if objectives is None else objectives
    row = {"cloudy_rows_pct": 100.0 * m.n_cloudy / max(1, len(m))}
    reach = reachable_tile_count(event_data, m)
    row["reachable_tiles"] = reach
    for obj, rate in objectives:
        ta = rate > 0
        w = time_aware_weights(event_data, rate) if ta else None
        sel_a, _, _, _ = select_passes(event_data, MAX_SELECTION_STEPS, cloud_mask=m, gap_fill=True,
                                       time_aware=ta, latency_decay_rate=rate)
        stats = dict(last_adaptive_run)
        n = len(sel_a)
        plans = {"adaptive": sel_a,
                 "openloop_expected": plan_expected_coverage(event_data, n, _ptile, w),
                 "openloop_clear": run_submodular_greedy_selection(event_data, n, cloud_mask=None, gap_fill=False,
                                                                   time_aware=ta, latency_decay_rate=rate)[0],
                 "oracle": run_submodular_greedy_selection(event_data, n, cloud_mask=m, gap_fill=False,
                                                           time_aware=ta, latency_decay_rate=rate)[0]}
        if globals().get("CLOUD_KNOWLEDGE", "realistic") != "oracle":
            sel_f = simulate_adaptive_selection(event_data, MAX_SELECTION_STEPS, m, gap_fill=True, time_aware=ta,
                                                latency_decay_rate=rate, planning="clear")[0]
            plans["failsignal_v140"] = sel_f
            plans["openloop_expected_v140n"] = plan_expected_coverage(event_data, len(sel_f), _ptile, w)
        for tag, sel in plans.items():
            c = realized_coverage(event_data, m, sel)
            row[f"{obj}_{tag}_images"] = len(sel)
            row[f"{obj}_{tag}_cov_pct"] = 100.0 * c / total_target_tiles
            row[f"{obj}_{tag}_reach_pct"] = (100.0 * c / reach) if reach else None
        row[f"{obj}_replans"] = int(stats.get("replans", 0))
        row[f"{obj}_compensatory"] = int(stats.get("compensatory", 0))
    if with_h1:
        need = int(np.ceil(_H1R / 100.0 * reach - 1e-9)) if reach else 0
        fl = coverage_times(event_data, m, (_H1R,), denominator=reach)[_H1R] if reach else None
        te, _, ne = earliest_completion(event_data, m, MAX_SELECTION_STEPS, _H1R, denominator=reach) if reach else (None, [], 0)
        tf, nf, _ = (level_stats(event_data, m, select_fcfs(event_data, m, len(event_data))[0], need) if reach
                     else (None, None, None))
        row.update({f"floor_t{_H1R:g}r_min": fl, f"ecc_t{_H1R:g}r_min": te,
                    f"ecc_images{_H1R:g}r": ne if te is not None else None,
                    f"fcfs_t{_H1R:g}r_min": tf, f"fcfs_images{_H1R:g}r": nf})
    return row


def _numcol(df, c):
    return (pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float) if c in df.columns
            else np.full(len(df), np.nan))


def _h2_summary(df, obj):
    A, X, C, O = (_numcol(df, f"{obj}_{t}_cov_pct") for t in _PLANS)
    adapted = int(np.nansum(_numcol(df, f"{obj}_replans")) + np.nansum(_numcol(df, f"{obj}_compensatory")))
    room = bool(np.any(O > X + 1e-9))
    d_x, d_c = paired_summary(A - X), paired_summary(A - C)
    captured = [(a - x) / (o - x) for a, x, o in zip(A, X, O) if np.isfinite(o - x) and o - x > 1e-9]
    reasons = []
    if globals().get("CLOUD_KNOWLEDGE", "realistic") == "oracle":
        reasons.append("cloud knowledge is 'oracle', so the selector does not adapt")
    if adapted == 0:
        reasons.append("the selector never re-planned or compensated in any realisation (nothing to adapt)")
    if not room:
        reasons.append("even the greedy that knows the cloud never beats the open-loop plan with the same images "
                       "(no room to adapt)")
    evaluable = not reasons and d_x["n"] >= 2
    status = ("NOT EVALUABLE" if not evaluable else
              ("SUPPORTED" if d_x["lcb"] is not None and d_x["lcb"] > 0 else "NOT SUPPORTED"))
    F, XF = _numcol(df, f"{obj}_failsignal_v140_cov_pct"), _numcol(df, f"{obj}_openloop_expected_v140n_cov_pct")
    return {"status": status, "evaluable": evaluable, "not_evaluable_reasons": reasons,
            "n_realisations": int(len(A)),
            "planning": globals().get("CLOUD_PLANNING", "clear"),
            "coverage_pct_mean": {t: (float(np.nanmean(v)) if np.isfinite(v).any() else None)
                                  for t, v in zip(_PLANS, (A, X, C, O))},
            "images_mean": float(np.nanmean(_numcol(df, f"{obj}_adaptive_images"))) if len(df) else None,
            "difference_vs_openloop_expected": d_x, "difference_vs_openloop_clear": d_c,
            "v140_failsignal": ({"coverage_pct_mean": float(np.nanmean(F)),
                                 "images_mean": float(np.nanmean(_numcol(df, f"{obj}_failsignal_v140_images"))),
                                 "difference_vs_openloop_expected_same_images": paired_summary(F - XF)}
                                if np.isfinite(F).any() else None),
            "share_of_oracle_gain_captured_median": float(np.median(captured)) if captured else None,
            "replans_total": int(np.nansum(_numcol(df, f"{obj}_replans"))),
            "compensatory_total": int(np.nansum(_numcol(df, f"{obj}_compensatory")))}


def _h1_under_cloud(df):
    fl, te = _numcol(df, f"floor_t{_H1R:g}r_min"), _numcol(df, f"ecc_t{_H1R:g}r_min")
    ne, nf = _numcol(df, f"ecc_images{_H1R:g}r"), _numcol(df, f"fcfs_images{_H1R:g}r")
    ok = np.isfinite(fl)
    within = ok & np.isfinite(te) & (te - fl <= H1_FLOOR_MARGIN_MIN + 1e-9)
    elig = ok & np.isfinite(nf) & (nf >= 2)
    econ = elig & within & np.isfinite(ne) & (ne < nf)
    return {"n_realisations": int(ok.sum()),
            "pct_ecc_within_margin": float(100.0 * within.sum() / ok.sum()) if ok.sum() else None,
            "n_realisations_fcfs_two_plus_images": int(elig.sum()),
            "pct_ecc_fewer_images_within_margin": float(100.0 * econ.sum() / elig.sum()) if elig.sum() else None,
            "median_images_ecc": float(np.median(ne[elig & np.isfinite(ne)])) if (elig & np.isfinite(ne)).any() else None,
            "median_images_fcfs": float(np.median(nf[elig])) if elig.any() else None}


if _cs_n >= 2 and len(event_data):
    logger.info(f"Cloud robustness: {_cs_n} cloud realisations ({globals().get('CLOUD_MODEL', 'iid')} model) ...")
    for _i in range(_cs_n):
        _m = build_cloud_mask(event_data, CLOUD_PROBABILITY, RANDOM_SEED + _i)
        cloud_seed_rows.append({"seed": RANDOM_SEED + _i, **_cloud_trial(_m)})
    cloud_seed_df = pd.DataFrame(cloud_seed_rows)
    cloud_seed_df.to_csv(os.path.join(_out_dir, "arf_cloud_robustness.csv"), index=False)
    cloud_robustness_summary = {
        "design": ("Per cloud realisation, four plans with the same number of images: the Fail-Signal adaptive "
                   "selector; the same expected-coverage planner run once in advance with no cloud observation "
                   "(open-loop); an open-loop plan on clear-sky geometry; and the greedy given the cloud in advance. Also "
                   "the v1.4.0 Fail-Signal (clear-sky planning). Coverage = tiles with a usable cloud-free "
                   "observation, % of ROI tiles."),
        "statement": ("H2 (v1.5.0): with the same number of images, the Fail-Signal adaptive selector covers more of "
                      "the region under cloud than the same planner run in advance without seeing any image (open-"
                      "loop, maximising expected coverage from the tiles' cloud probabilities): the mean paired "
                      "difference over cloud realisations has a one-sided 95% lower confidence bound above zero."),
        "cloud_planning": globals().get("CLOUD_PLANNING", "clear"),
        "primary_objective": "time_aware",
        "cloud_model": globals().get("cloud_model_info", {"model": globals().get("CLOUD_MODEL", "iid")}),
        "seeds": [int(RANDOM_SEED), int(RANDOM_SEED + _cs_n - 1)],
        "h2": {obj: _h2_summary(cloud_seed_df, obj) for obj, _ in _OBJECTIVES},
        "h1_under_cloud": _h1_under_cloud(cloud_seed_df),
    }

    # =========================================================================
    # 5. CLOUD SENSITIVITY: the correlation length and time (one at a time)
    # =========================================================================
    _combos = [(float(L), None if float(L) == 0.0 else float(CLOUD_CORR_TIME_H))
               for L in globals().get("CLOUD_SWEEP_LENGTHS_KM", ())]
    _combos += [(float(CLOUD_CORR_LENGTH_KM), float(T)) for T in globals().get("CLOUD_SWEEP_TIMES_H", ())
                if float(T) != float(CLOUD_CORR_TIME_H)]
    _combos = list(dict.fromkeys(_combos))
    _sw_n = int(globals().get("CLOUD_SWEEP_SEEDS", 10))
    if _combos and _sw_n >= 2:
        logger.info(f"Cloud sensitivity: {len(_combos)} settings x {_sw_n} realisations ...")
        for L, T in _combos:
            _rows = []
            for _i in range(_sw_n):
                _m = (build_cloud_mask(event_data, CLOUD_PROBABILITY, RANDOM_SEED + _i, model="iid") if L == 0.0 else
                      build_cloud_mask(event_data, CLOUD_PROBABILITY, RANDOM_SEED + _i, model="correlated",
                                       length_km=L, time_h=T))
                _rows.append(_cloud_trial(_m, objectives=_OBJECTIVES[:1], with_h1=False))
            _d = pd.DataFrame(_rows)
            _A, _X, _C, _O = (_numcol(_d, f"time_aware_{t}_cov_pct") for t in _PLANS)
            _ps = paired_summary(_A - _X)
            cloud_sensitivity_rows.append({
                "model": "iid" if L == 0.0 else "correlated", "length_km": L, "time_h": T,
                "primary": (globals().get("CLOUD_MODEL") == "correlated" and L == float(CLOUD_CORR_LENGTH_KM)
                            and T == float(CLOUD_CORR_TIME_H)),
                "n_realisations": int(len(_d)), "images_mean": float(np.nanmean(_numcol(_d, "time_aware_adaptive_images"))),
                "adaptive_cov_pct": float(np.nanmean(_A)), "openloop_expected_cov_pct": float(np.nanmean(_X)),
                "openloop_clear_cov_pct": float(np.nanmean(_C)), "oracle_cov_pct": float(np.nanmean(_O)),
                "diff_mean": _ps["mean"], "diff_lcb95": _ps["lcb"], "diff_pct_positive": _ps["pct_positive"],
                "replans_total": int(np.nansum(_numcol(_d, "time_aware_replans"))),
                "compensatory_total": int(np.nansum(_numcol(_d, "time_aware_compensatory")))})
        pd.DataFrame(cloud_sensitivity_rows).to_csv(os.path.join(_out_dir, "arf_cloud_sensitivity.csv"), index=False)
        cloud_robustness_summary["sensitivity"] = cloud_sensitivity_rows
    with open(os.path.join(_out_dir, "arf_cloud_robustness_summary.json"), "w") as f:
        json.dump(cloud_robustness_summary, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))

    # Figure 8 (IEEE double-column width): (a) coverage of the four plans with
    # the same images over the realisations; (b), (c) the H2 difference against
    # the correlation length and time (mean, with its one-sided 95% lower bound).
    _cols8 = {"adaptive": "#2a78d6", "openloop_expected": "#eb6834", "openloop_clear": "#1baf7a",
              "oracle": "#8a8986"}
    _n_pan = 3 if cloud_sensitivity_rows else 1
    with plt.rc_context({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
                         "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
                         "axes.spines.right": False}):
        # v1.5.3: 7.0 in (was 7.16): with the tight bounding box the saved figure grew to 7.18 in when
        # panels (b) and (c) carry longer tick labels (open data alone, 15 October 2026)
        fig8, axs8 = plt.subplots(1, _n_pan, figsize=(7.0, 2.6),
                                  gridspec_kw={"width_ratios": [1.4, 1, 1][:_n_pan]})
        axs8 = np.atleast_1d(axs8)
        ax = axs8[0]
        _rng8 = np.random.default_rng(0)
        for j, tag in enumerate(_PLANS):
            v = _numcol(cloud_seed_df, f"time_aware_{tag}_cov_pct")
            v = v[np.isfinite(v)]
            if not v.size:
                continue
            ax.boxplot([v], positions=[j], widths=0.5, showfliers=False, medianprops={"color": "#0b0b0b", "lw": 1.2},
                       boxprops={"color": "#52514e", "lw": 0.8}, whiskerprops={"color": "#52514e", "lw": 0.8},
                       capprops={"color": "#52514e", "lw": 0.8})
            ax.scatter(j + _rng8.uniform(-0.13, 0.13, v.size), v, s=9, color=_cols8[tag], alpha=0.85, lw=0, zorder=3)
        ax.set_xticks(range(len(_PLANS)))
        ax.set_xticklabels(["Fail-Signal\n(adaptive)", "Open-loop\n(expected)", "Open-loop\n(clear sky)",
                            "Greedy,\ncloud known"])
        ax.set_ylabel("Coverage under cloud (% of tiles)")
        _h2p = cloud_robustness_summary["h2"]["time_aware"]
        ax.set_title(f"(a) Same images, {_cs_n} cloud realisations", loc="left")
        ax.grid(axis="y", alpha=0.3, lw=0.5)
        if _n_pan == 3:
            _sens = pd.DataFrame(cloud_sensitivity_rows)
            for ax, (key, sub_df, xlab) in zip(axs8[1:], (
                    ("length_km", _sens[(_sens.model == "iid") | (_sens.time_h == float(CLOUD_CORR_TIME_H))],
                     f"Correlation length (km), time {CLOUD_CORR_TIME_H:g} h"),
                    ("time_h", _sens[(_sens.model == "correlated") & (_sens.length_km == float(CLOUD_CORR_LENGTH_KM))],
                     f"Correlation time (h), length {CLOUD_CORR_LENGTH_KM:g} km"))):
                sub_df = sub_df.sort_values(key, na_position="first")
                labels = ["iid" if (key == "length_km" and r.model == "iid") else f"{r[key]:g}" for _, r in sub_df.iterrows()]
                xs = np.arange(len(sub_df))
                mean = sub_df["diff_mean"].to_numpy(dtype=float)
                lcb = sub_df["diff_lcb95"].to_numpy(dtype=float)
                ax.axhline(0.0, color="#52514e", lw=0.8)
                ax.vlines(xs, lcb, mean, color="#2a78d6", lw=1.6)
                ax.plot(xs, mean, "o", color="#2a78d6", ms=4)
                prim = sub_df["primary"].to_numpy(dtype=bool)
                if prim.any():
                    ax.plot(xs[prim], mean[prim], "o", mfc="none", mec="#0b0b0b", ms=7, mew=1.0)
                ax.set_xticks(xs)
                ax.set_xticklabels(labels)
                ax.set_xlabel(xlab)
                ax.grid(axis="y", alpha=0.3, lw=0.5)
            axs8[1].set_ylabel("Adaptive minus open-loop\n(coverage points)")
            axs8[1].set_title("(b) Sensitivity to cloud size", loc="left")
            axs8[2].set_title("(c) Sensitivity to cloud persistence", loc="left")
        fig8.tight_layout()
        for fmt in ("png", "pdf"):
            fig8.savefig(os.path.join(_out_dir, f"arf_fig8_cloud_robustness.{fmt}"), dpi=300, bbox_inches="tight")
        plt.show()

if 'scenario_results' in globals():
    scenario_results.setdefault(_scenario_key(), {})
    scenario_results[_scenario_key()]["decay_sweep"] = decay_sweep_rows
    scenario_results[_scenario_key()]["onset_robustness"] = onset_robustness_summary
    scenario_results[_scenario_key()]["policy_comparison"] = policy_rows
    scenario_results[_scenario_key()]["ilp_gap"] = ilp_summary
    scenario_results[_scenario_key()]["cost_latency_frontier"] = frontier_rows
    scenario_results[_scenario_key()]["cloud_robustness"] = cloud_robustness_summary
    scenario_results[_scenario_key()]["cloud_model"] = globals().get("cloud_model_info")

print("\n" + "=" * 70)
print("  CELL 7C: Objective Comparison & Onset Robustness (H1)")
print("=" * 70)
print(decay_sweep_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\n  Policy comparison (full window, budget k, no gap-fill):")
print(policy_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
if ilp_summary:
    print(f"\n  Optimality gap ({ilp_summary['solved']}/{ilp_summary['instances']} instances solved, "
          f"{ilp_summary['optimal']} to proven optimality; {ILP_HORIZON_HOURS} h of passes, k={ILP_BUDGET_K}, clear sky):")
    for n in ("coverage_first", "time_aware", "fcfs"):
        mi, md = ilp_summary[f"{n}_ratio_min"], ilp_summary[f"{n}_ratio_median"]
        print(f"    {n:16s} utility / optimum: median {md:.4f}, worst {mi:.4f}" if md is not None else f"    {n}: n/a")
if onset_robustness_summary:
    print(f"\n  {len(onset_df)} onsets | superseded H1 target {_H1_TARGET_MIN:g} min (kept for the record)")
    print(f"  {'series':26s} {'median':>8s} {'p90':>8s} {'max':>8s} {'% <10min':>8s}")
    for k, s in onset_robustness_summary["series"].items():
        f = lambda v: f"{v:8.0f}" if v is not None else "     n/a"
        print(f"  {k:26s} {f(s['median_min'])} {f(s['p90_min'])} {f(s['max_min'])} {s['pct_onsets_under_h1']:7.1f}%")
    _a = onset_robustness_summary["h1_revised"]["series"]["earliest_completion_clear"]
    _at = onset_robustness_summary["h1_revised"]["series"]["time_aware_clear"]
    _b = onset_robustness_summary["h1b_image_economy"]["series"]["clear"]
    _p = lambda v: "n/a" if v is None else f"{v:.1f}%"
    print(f"\n  H1a (clear): ECC within {H1_FLOOR_MARGIN_MIN:.1f} min of the floor at {_a['n_within_margin']}/"
          f"{_a['n_onsets_evaluated']} onsets ({_p(_a['pct_within_margin'])}; 95% lower bound "
          f"{_p(_a['pct_within_margin_lcb95'])}; evaluable from {H1_MIN_EVALUABLE_ONSETS} onsets); "
          f"time-aware {_p(_at['pct_within_margin'])}")
    print(f"  H1b (clear): ECC fewer images than FCFS, within the margin, at {_b['n_success']}/{_b['n_onsets_evaluated']} "
          f"onsets ({_p(_b['pct_success'])}); median images ECC {_b['median_images_ecc']} vs FCFS "
          f"{_b['median_images_fcfs']}; {_b['n_onsets_one_image']} onsets need one image")
if cloud_robustness_summary:
    for _obj in ("time_aware", "coverage_first"):
        _h = cloud_robustness_summary["h2"][_obj]
        _dx = _h["difference_vs_openloop_expected"]
        _f = lambda v, n=2: "n/a" if v is None else f"{v:.{n}f}"
        print(f"  H2 ({_obj}{', primary' if _obj == 'time_aware' else ''}): {_h['status']} -- adaptive minus "
              f"open-loop (expected), same images: mean {_f(_dx['mean'])} points, one-sided 95% lower bound "
              f"{_f(_dx['lcb'])} over {_dx['n']} realisations; mean images {_f(_h['images_mean'], 1)}")
print("=" * 70)
print("\n=== CELL 7C COMPLETE ===")

