"""PlanetScope acquisition-fraction replay (proposed post hoc amendment): how the
floor, the reachable share and the ECC/FCFS results change when only a fraction
f of PlanetScope's predicted daylight passes acquire an image.

Why: the model assumes every PlanetScope pass acquires; Planet describes that as
its default ("image all land areas at every opportunity"), subject to
proprietary priorities, so it is an upper bound. This replay bounds the effect
of the assumption without Planet's catalogue.

How:
  * Inputs are a completed run folder of the release of record: arf_refined_events.csv
    (and, for the identity check, arf_onset_robustness.csv; the H1a margin from
    arf_cell2_summary.json or --margin-min).
  * A PlanetScope pass is a maximal run of consecutive minutes of one FLOCK
    satellite. Each pass gets one uniform draw u per seed; it acquires when u < f,
    so the fractions are nested (every pass kept at f = 0.25 is kept at 0.5).
    Sentinel-2 and Landsat passes are always kept (validated one for one).
  * Per onset (0 to 24 h in 20-min steps, 24-h horizon), the floor, the reachable
    count, ECC (k = 30, 95% of reachable, bisection over delivery times with a
    greedy probe) and FCFS's image count are recomputed with the rules of the
    pipeline's Cell 7C (earliest_completion, select_fcfs + level_stats,
    coverage_times), ported to arrays.
  * Identity check: at f = 1 every onset's floor, reachable count, ECC time, ECC
    images and FCFS images must equal the run's arf_onset_robustness.csv, or the
    script stops before any thinned result is written.
Outputs: OUT_PREFIX_onsets.csv (one row per fraction, seed, onset) and
OUT_PREFIX_summary.json (per fraction: means over seeds of the window's H1a and
H1b shares by the protocol's rules, median floor, reachable share and images).

Usage: arf_planet_thinning.py RUN_DIR OUT_PREFIX [--fractions 1,0.75,0.5,0.25]
         [--seeds 20] [--k 30] [--margin-min M]
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

LEVEL = 95.0
HORIZON, SHIFT_H, STEP = 1440, 24, 20
PLANET = re.compile(r"^FLOCK", re.I)          # MISSION_TABLE name pattern of PlanetScope (release 1.5.3)


def load(run_dir):
    ev = pd.read_csv(os.path.join(run_dir, "arf_refined_events.csv"))
    util = (np.clip(1.0 - ev["A_nadir_deg"].to_numpy(float) / ev["theta_max_deg"].to_numpy(float), 0.0, None)
            * ev["D_downlink"].to_numpy(float))
    ev = ev.assign(util=util)
    tiles = pd.Index(pd.unique(ev["tile"]))
    ev = ev.assign(tid=tiles.get_indexer(ev["tile"]))
    keys, cells, wait = [], [], []
    for (sat, t), g in ev.groupby(["sat", "t_min"], sort=True):          # the EventStore order: (sat, t_min)
        keys.append((sat, int(t)))
        cells.append(g.loc[g["util"] > 0.0, "tid"].to_numpy(np.int64))
        w = g["T_wait_min"].to_numpy(float)[0]
        wait.append(np.inf if np.isnan(w) else float(w))
    return keys, cells, np.array(wait), len(tiles)


def planet_passes(keys):
    """Pass id per event for PlanetScope (consecutive minutes of one satellite), -1 otherwise."""
    pid = np.full(len(keys), -1, dtype=np.int64)
    n, prev = -1, None
    for e, (sat, t) in enumerate(keys):
        if not PLANET.match(sat):
            prev = None
            continue
        if prev is None or prev[0] != sat or t - prev[1] > 1:
            n += 1
        pid[e] = n
        prev = (sat, t)
    return pid, n + 1


def onset_metrics(evs, keys, cells, wait, n_tiles, d, k):
    """Floor, reachable, ECC (time, images) and FCFS images for onset d over events `evs` (store order)."""
    sub = [e for e in evs if d <= keys[e][1] < d + HORIZON]
    if not sub:
        return None
    reach = np.unique(np.concatenate([cells[e] for e in sub])).size if sub else 0
    if reach == 0:
        return None
    need = int(math.ceil(LEVEL / 100.0 * reach - 1e-9))
    dl = np.array([keys[e][1] - d + wait[e] for e in sub])
    order = np.argsort(dl, kind="stable")
    order = order[np.isfinite(dl[order])]
    times = dl[order]
    # floor (coverage_times with the reachable denominator) and FCFS images (select_fcfs + level_stats)
    seen = np.zeros(n_tiles, bool); n_seen = 0; floor = None; fcfs = None; n_kept = 0; pos_floor = None
    for pos, o in enumerate(order):
        c = cells[sub[o]]
        new = c[~seen[c]]
        if new.size:
            seen[new] = True; n_seen += new.size; n_kept += 1
            if floor is None and 100.0 * n_seen / reach >= LEVEL - 1e-9:
                floor = float(times[pos])
            if fcfs is None and n_seen >= need:
                fcfs = n_kept; pos_floor = pos
        if floor is not None and fcfs is not None:
            break
    # ECC: earliest_completion of the pipeline (bisection; greedy probe on the images delivered by times[pos])
    def probe(pos):
        end = int(np.searchsorted(times, times[pos], side="right"))
        cand = [sub[o] for o in order[:end]]
        covered = np.zeros(n_tiles, bool); avail = np.ones(len(cand), bool); chosen = []; n_cov = 0
        for _ in range(k):
            gain = np.array([(~covered[cells[e]]).sum() if avail[i] else -1 for i, e in enumerate(cand)])
            j = int(np.argmax(gain))
            if gain[j] <= 0:
                break
            avail[j] = False; chosen.append(cand[j])
            c = cells[cand[j]]
            n_cov += int((~covered[c]).sum()); covered[c] = True
            if n_cov >= need:
                break
        return n_cov >= need, chosen
    ecc_t = ecc_n = None
    if pos_floor is not None:
        hi = len(order) - 1
        ok, best = probe(hi)
        if ok:
            a, b = pos_floor, hi
            while a < b:
                mid = (a + b) // 2
                okm, ch = probe(mid)
                if okm:
                    b, best = mid, ch
                else:
                    a = mid + 1
            ecc_t = float(max(keys[e][1] - d + wait[e] for e in best)); ecc_n = len(best)
    return {"reachable": reach, "floor_min": floor, "ecc_min": ecc_t, "ecc_images": ecc_n, "fcfs_images": fcfs}


def window_stats(df, margin):
    """H1a and H1b shares of one window by the protocol's rules (Cell 7C), plus medians."""
    ok = df.dropna(subset=["floor_min"])
    h1a = ((ok["ecc_min"] - ok["floor_min"]) <= margin + 1e-9) & ok["ecc_min"].notna()
    elig = ok[ok["fcfs_images"] >= 2]
    h1b = (((elig["ecc_min"] - elig["floor_min"]) <= margin + 1e-9) & (elig["ecc_images"] < elig["fcfs_images"]))
    return {"evaluable_onsets": int(len(ok)), "h1a_pct": 100.0 * float(h1a.mean()) if len(ok) else None,
            "h1b_eligible": int(len(elig)), "h1b_pct": 100.0 * float(h1b.mean()) if len(elig) else None,
            "floor_h_median": float(ok["floor_min"].median()) / 60 if len(ok) else None,
            "ecc_images_median": float(ok["ecc_images"].median()) if len(ok) else None,
            "fcfs_images_median": float(ok["fcfs_images"].median()) if len(ok) else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir"); ap.add_argument("out_prefix")
    ap.add_argument("--fractions", default="1,0.75,0.5,0.25")
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--margin-min", type=float, default=None)
    a = ap.parse_args()
    margin = a.margin_min
    c2 = os.path.join(a.run_dir, "arf_cell2_summary.json")
    if margin is None and os.path.exists(c2):
        margin = float(json.load(open(c2))["t_max_wait_min"])
    if margin is None:
        sys.exit("H1a margin unknown: no arf_cell2_summary.json; pass --margin-min")
    keys, cells, wait, n_tiles = load(a.run_dir)
    pid, n_pass = planet_passes(keys)
    total_tiles = None
    mt = os.path.join(a.run_dir, "arf_mosaic_tiles.csv")
    if os.path.exists(mt):
        total_tiles = len(pd.read_csv(mt, usecols=["tile"]))
    shifts = list(range(0, SHIFT_H * 60 + 1, STEP))

    # ---------------------------------------------------------------- identity check at f = 1
    full = list(range(len(keys)))
    ref = pd.read_csv(os.path.join(a.run_dir, "arf_onset_robustness.csv")).set_index("onset_shift_min")
    cols = {"reachable": "reachable_clear_tiles", "floor_min": "best_possible_clear_t95r_min",
            "ecc_min": "earliest_completion_clear_t95r_min", "ecc_images": "earliest_completion_clear_images95r",
            "fcfs_images": "fcfs_unbudgeted_clear_images95r"}
    bad = []
    for d in shifts:
        m = onset_metrics(full, keys, cells, wait, n_tiles, d, a.k) or {}
        for mine, theirs in cols.items():
            x, y = m.get(mine), ref.loc[d, theirs]
            if (x is None) != pd.isna(y) or (x is not None and not np.isclose(float(x), float(y))):
                bad.append((d, mine, x, y))
    if bad:
        print("IDENTITY CHECK FAILED (first 10):", bad[:10]); sys.exit(1)
    print(f"identity check: {len(shifts)} onsets x {len(cols)} columns equal the run's arf_onset_robustness.csv")

    # ---------------------------------------------------------------- thinned replays
    fr = [float(x) for x in a.fractions.split(",")]
    rows = []
    for seed in range(1, a.seeds + 1):
        u = np.random.default_rng(seed).random(n_pass)
        for f in fr:
            if f >= 1.0 and seed > 1:
                continue                                 # nothing is thinned at f = 1: one replay suffices
            keep = [e for e in full if pid[e] < 0 or u[pid[e]] < f]
            for d in shifts:
                m = onset_metrics(keep, keys, cells, wait, n_tiles, d, a.k)
                rows.append({"fraction": f, "seed": seed, "onset_min": d, **(m or {})})
    out = pd.DataFrame(rows)
    out.to_csv(a.out_prefix + "_onsets.csv", index=False)
    summ = {"run_dir": os.path.basename(os.path.normpath(a.run_dir)), "planetscope_passes": int(n_pass),
            "events": len(keys), "margin_min": margin, "k": a.k, "seeds": a.seeds, "identity_check": "passed",
            "fractions": {}}
    for f in fr:
        per_seed = [window_stats(g, margin) for _, g in out[out.fraction == f].groupby("seed")]
        reach_pct = (out[out.fraction == f]["reachable"] / total_tiles * 100).median() if total_tiles else None
        summ["fractions"][str(f)] = {k_: (float(np.mean([s[k_] for s in per_seed if s[k_] is not None]))
                                          if any(s[k_] is not None for s in per_seed) else None)
                                     for k_ in per_seed[0]}
        summ["fractions"][str(f)]["reachable_pct_median"] = None if reach_pct is None else float(reach_pct)
    print(json.dumps(summ, indent=1))
    json.dump(summ, open(a.out_prefix + "_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
