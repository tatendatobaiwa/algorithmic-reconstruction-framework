# =============================================================================
# ARF EXACT TIME AT TIGHT BUDGETS AND THE FAIL-SIGNAL IN HOURS
# (protocol amendment A9, POST HOC; read-only)
# =============================================================================
# Two questions a reviewer asks, answered on a replay of Cells 6-8 of release
# 1.5.3 on each primary run's own saved inputs (the procedure of
# arf_exact_economy_a4.py: the run's own altitude and period; the replay must
# reproduce the run's arf_onset_robustness.csv and arf_cloud_robustness.csv,
# exactly in the columns used, else the script stops). Windows: the 8 dry
# primary windows and, when present, the 8 wet ones of A5 (p = 0.35 runs).
#
# PART 1, exact earliest time at a tight budget (clear sky, 95% of reachable).
#   For k in (5, 10, 15, 20) and a fixed sample of onsets (every 80 min; all 73
#   with --all-onsets): ECC (earliest_completion, budget k) against the exact
#   earliest time T*(k) = the first delivery time t at which some k images
#   delivered by t cover the level. Feasibility at t is decided exactly by the
#   release's solve_ilp_min_images (HiGHS) on the events delivered by t; it is
#   monotone in t, so T*(k) is found by bisection over the delivery times.
#     ECC at the floor            -> optimal, no solve
#     ECC later than the floor    -> solve at the floor; if feasible, ECC's
#                                    delay was avoidable; else bisection
#     ECC found no cover          -> solve at the horizon end: infeasible means
#                                    no k-cover exists (a property of the passes)
#   A solve that hits the time limit is reported as unknown, never as a result.
#
# PART 2, the Fail-Signal in hours (time-aware objective, the H2 design).
#   For each of the run's cloud realisations (build_cloud_mask, seeds 42..61),
#   the adaptive Fail-Signal plan and the open-loop expected-coverage plan of
#   the same size are rebuilt with the release's functions (their image counts
#   and coverage must equal the run's arf_cloud_robustness.csv). Then, from
#   coverage_times with each plan as the selection and the tiles reachable
#   under that cloud as the denominator:
#     L      = the open-loop plan's final realised coverage (% of reachable)
#     dT(L)  = open-loop's time to L minus the adaptive plan's time to L (h):
#              how much sooner the Fail-Signal reaches what open-loop achieves
#     and the time each reaches 50% and 80% of reachable tiles.
#   Times are delivery times from the window start (48-h window, as in H2).
#
# Writes ARF_batch/exact_time_a9/ (a9_budget.csv, a9_failsignal.csv,
# a9_summary.json, tab_a9_budget.tex, tab_a9_failsignal.tex). Resumable: each
# finished window is cached in exact_time_a9/cache.
#
# Use (Colab):  !pip install -q h3==4.5.0
#               !python /content/drive/MyDrive/ARF_code/arf_exact_time_a9.py
# =============================================================================
import os
import json
import time
import runpy
import argparse

import numpy as np
import pandas as pd

LEVEL = 95.0
BUDGETS = (5, 10, 15, 20)
SAMPLE_MIN = 80
TIME_LIMIT_S = 300.0
MAX_PAIRS = 3_000_000
MAX_BISECT = 14
FS_LEVELS = (50.0, 80.0)
INDEXES = (("confirmatory_v153_primary_index.csv", "dry"), ("wetseason_v153_a5_index.csv", "wet"))
CLOUD_EXACT = ("seed", "reachable_tiles", "time_aware_adaptive_images", "time_aware_openloop_expected_images")


def feasible(ns, ex, sub, t, need, k):
    """True / False / None (unknown): some <= k images delivered by t cover `need` tiles."""
    dl = ex["delivered_by"](ns, sub, t)
    if len(dl) == 0:
        return False, 0, "empty"
    r = ns["solve_ilp_min_images"](dl, None, need, TIME_LIMIT_S, MAX_PAIRS)
    if r["status"] == "infeasible":
        return False, None, r["status"]
    if r["status"] != "optimal":
        return None, None, r["status"]
    return r["images"] <= k, r["images"], "optimal"


def exact_time(ns, ex, sub, need, k, floor, t_e):
    """(T*(k) or None, verdict, solves, unknown)."""
    d = ns["event_delivery_times"](sub)
    times = sorted({float(x) for x in d if np.isfinite(x) and x >= floor - 1e-9})
    solves = 0
    if t_e is not None and t_e <= floor + 1e-9:
        return floor, "ECC at floor", 0, False
    if t_e is None:
        ok, _, st = feasible(ns, ex, sub, times[-1], need, k); solves += 1
        if ok is None:
            return None, f"unknown ({st})", solves, True
        if not ok:
            return None, "no k-cover exists", solves, False
        hi = len(times) - 1
    else:
        hi = max(i for i, t in enumerate(times) if t <= t_e + 1e-9)
    ok, _, st = feasible(ns, ex, sub, times[0], need, k); solves += 1
    if ok is None:
        return None, f"unknown ({st})", solves, True
    if ok:
        return times[0], ("ECC delay avoidable" if t_e is not None else "ECC missed a cover"), solves, False
    lo = 0                                           # infeasible at times[lo], feasible at times[hi]
    while hi - lo > 1:
        if solves >= MAX_BISECT:
            return None, "unknown (bisection cap)", solves, True
        mid = (lo + hi) // 2
        ok, _, st = feasible(ns, ex, sub, times[mid], need, k); solves += 1
        if ok is None:
            return None, f"unknown ({st})", solves, True
        lo, hi = (lo, mid) if ok else (mid, hi)
    T = times[hi]
    if t_e is None:
        return T, "ECC missed a cover", solves, False
    return T, ("ECC at exact time" if abs(T - t_e) < 1e-9 else "ECC later than exact"), solves, False


def part_budget(ns, ex, all_onsets=False, label=""):
    H = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
    rows = []
    for d in range(0, ns["ONSET_SHIFT_HOURS"] * 60 + 1, ns["ONSET_SHIFT_STEP_MIN"]):
        if not all_onsets and d % SAMPLE_MIN:
            continue
        sub, _ = ns["subset_by_time"](ns["event_data"], None, d, d + H)
        R = ns["reachable_tile_count"](sub, None) if len(sub) else 0
        if not R:
            continue
        need = int(np.ceil(LEVEL / 100.0 * R - 1e-9))
        floor = ns["coverage_times"](sub, None, (LEVEL,), denominator=R)[LEVEL]
        for k in BUDGETS:
            t_e, _, n_e = ns["earliest_completion"](sub, None, k, LEVEL, denominator=R)
            t0 = time.time()
            T, verdict, solves, unk = exact_time(ns, ex, sub, need, k, floor, t_e)
            rows.append({"onset_shift_min": d, "k": k, "reachable_tiles": R, "floor_min": floor, "ecc_min": t_e,
                         "ecc_images": n_e if t_e is not None else None, "exact_min": T, "verdict": verdict,
                         "solves": solves, "unknown": unk, "seconds": round(time.time() - t0, 1),
                         "ecc_excess_min": (t_e - T) if (t_e is not None and T is not None) else None})
            if solves:
                print(f"  {label} +{d} min k={k}: {verdict} (ECC {t_e}, exact {T}, {solves} solves)", flush=True)
    return pd.DataFrame(rows)


def realised_time(ns, m, sel, level, reach):
    return ns["coverage_times"](ns["event_data"], m, (level,), selected=sel, denominator=reach)[level]


def part_failsignal(ns, run_cloud_csv):
    run = pd.read_csv(run_cloud_csv)
    ed, rate = ns["event_data"], float(ns["LATENCY_DECAY_RATE"])
    w = ns["time_aware_weights"](ed, rate)
    rows = []
    for i in range(int(ns["CLOUD_SEEDS"])):
        seed = int(ns["RANDOM_SEED"]) + i
        m = ns["build_cloud_mask"](ed, ns["CLOUD_PROBABILITY"], seed)
        reach = ns["reachable_tile_count"](ed, m)
        a = ns["select_passes"](ed, ns["MAX_SELECTION_STEPS"], cloud_mask=m, gap_fill=True, time_aware=True,
                                latency_decay_rate=rate)[0]
        o = ns["plan_expected_coverage"](ed, len(a), ns["_ptile"], w)
        r = run[run.seed == seed].iloc[0]
        ca, co = ns["realized_coverage"](ed, m, a), ns["realized_coverage"](ed, m, o)
        n_t = ns["total_target_tiles"]
        assert (len(a) == r.time_aware_adaptive_images and len(o) == r.time_aware_openloop_expected_images
                and reach == r.reachable_tiles
                and abs(100.0 * ca / n_t - r.time_aware_adaptive_cov_pct) < 1e-9
                and abs(100.0 * co / n_t - r.time_aware_openloop_expected_cov_pct) < 1e-9), \
            f"seed {seed}: rebuilt plans differ from the run's arf_cloud_robustness.csv"
        row = {"seed": seed, "reachable_tiles": reach, "images": len(a),
               "adaptive_reach_pct": 100.0 * ca / reach if reach else None,
               "openloop_reach_pct": 100.0 * co / reach if reach else None}
        if reach and co:
            L = 100.0 * co / reach
            to, ta = realised_time(ns, m, o, L, reach), realised_time(ns, m, a, L, reach)
            row.update(level_L_pct=L, openloop_t_L_h=None if to is None else to / 60.0,
                       adaptive_t_L_h=None if ta is None else ta / 60.0,
                       dT_L_h=None if (to is None or ta is None) else (to - ta) / 60.0)
        for lv in FS_LEVELS:
            to, ta = realised_time(ns, m, o, lv, reach), realised_time(ns, m, a, lv, reach)
            row[f"openloop_t{lv:g}_h"] = None if to is None else to / 60.0
            row[f"adaptive_t{lv:g}_h"] = None if ta is None else ta / 60.0
        rows.append(row)
    return pd.DataFrame(rows)


def windows(batch, code):
    out = []
    for name, season in INDEXES:
        p = os.path.join(batch, name)
        if not os.path.exists(p):
            continue
        idx = pd.read_csv(p, dtype=str).fillna("")
        for _, r in idx[idx.status == "done"].iterrows():
            d = r["out_dir"]
            mp = os.path.join(d, "arf_reproducibility_manifest.json")
            if not os.path.exists(mp):
                continue
            man = json.load(open(mp))
            cfg = man.get("config_overrides") or {}
            if (cfg.get("active_roi"), cfg.get("fleet"), cfg.get("planet_swath_case", "conservative")) != \
                    ("ckgr_okavango", "earth-resources+planet", "conservative") or cfg.get("exclude_satellites") \
                    or cfg.get("cloud_probability_source", "constant") != "constant":
                continue
            out.append({"season": season, "window": str(cfg["analysis_epoch"])[:10], "run_dir": d, "man": man})
    seen, res = set(), []
    for w in sorted(out, key=lambda w: w["window"]):
        if w["window"] not in seen:
            seen.add(w["window"]); res.append(w)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description="Exact time at tight budgets and the Fail-Signal in hours (A9).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    ap.add_argument("--all-onsets", action="store_true")
    ap.add_argument("--only", default="", help="comma-separated windows (YYYY-MM-DD), for a partial run")
    a, _ = ap.parse_known_args(argv)
    ex = runpy.run_path(os.path.join(a.code, "arf_exact_economy_a4.py"), run_name="arf_exact_economy_a4")
    par = ex["load"](os.path.join(a.code, "arf_colab_parity.py"), "arf_colab_parity")
    pipe = os.path.join(a.code, "algorithmic_reconstruction_framework.py")
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    assert par.sha256_file(pipe) == proto["pipeline_sha256"], "ARF_code does not hold the release of record"
    out = os.path.join(a.batch_dir, "exact_time_a9")
    cache = os.path.join(out, "cache")
    os.makedirs(cache, exist_ok=True)
    only = {s.strip() for s in a.only.split(",") if s.strip()}
    B, F, info = [], [], []
    for w in windows(a.batch_dir, a.code):
        if only and w["window"] not in only:
            continue
        cp = os.path.join(cache, f"{w['season']}_{w['window']}.json")
        if os.path.exists(cp):
            c = json.load(open(cp))
            B.append(pd.DataFrame(c["budget"])); F.append(pd.DataFrame(c["fail"])); info.append(c["info"])
            print(f"{w['window']}: from cache", flush=True)
            continue
        d, man = w["run_dir"], w["man"]
        inp = man["inputs"]
        ov = man["output_validation"]["files"]
        refined, tiles_csv = os.path.join(d, "arf_refined_events.csv"), os.path.join(d, "arf_mosaic_tiles.csv")
        for p in (refined, tiles_csv):
            assert par.sha256_file(p) == ov[os.path.basename(p)]["sha256"], f"{p} differs from what Cell 9 recorded"
        key = list(man["scenario_results"])[0]
        orbit = man["scenario_results"][key]["orbital_period_min"]
        c2 = json.load(open(os.path.join(d, "arf_cell2_summary.json")))
        alt, torb = float(c2["median_alt_km"]), float(c2["t_max_wait_min"])
        assert torb == float(orbit), f"{d}: Cell 2's period differs from the manifest's"
        tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
        t0 = time.time()
        ns = ex["replay_exact"](par, alt, torb, pipe, refined, tiles, orbit, inp["active_roi"], inp["fleet"],
                                inp["config"].get("planet_swath_case", "conservative"), inp["n_satellites"])
        ident, ok, rep = ex["compare_onsets"](os.path.join(ns["_out_dir"], "arf_onset_robustness.csv"),
                                              os.path.join(d, "arf_onset_robustness.csv"))
        assert ok, f"{d}: the replay does not reproduce the run's onset results: {rep[:3]}"
        _, okc, repc = ex["compare_onsets"](os.path.join(ns["_out_dir"], "arf_cloud_robustness.csv"),
                                            os.path.join(d, "arf_cloud_robustness.csv"))
        assert okc, f"{d}: the replay does not reproduce the run's cloud realisations: {repc[:3]}"
        b = part_budget(ns, ex, a.all_onsets, label=w["window"])
        f = part_failsignal(ns, os.path.join(d, "arf_cloud_robustness.csv"))
        for df in (b, f):
            df.insert(0, "window", w["window"]); df.insert(0, "season", w["season"])
        inf = {"season": w["season"], "window": w["window"], "replay_identical": ident, "seconds": round(time.time() - t0, 1)}
        json.dump({"budget": b.to_dict("records"), "fail": f.to_dict("records"), "info": inf}, open(cp, "w"), default=str)
        B.append(b); F.append(f); info.append(inf)
        print(json.dumps(inf), flush=True)
    assert B, "no primary-design window found"
    B, F = pd.concat(B, ignore_index=True), pd.concat(F, ignore_index=True)
    for c in ("ecc_min", "exact_min", "ecc_excess_min", "floor_min"):
        B[c] = pd.to_numeric(B[c], errors="coerce")
    B.to_csv(os.path.join(out, "a9_budget.csv"), index=False)
    F.to_csv(os.path.join(out, "a9_failsignal.csv"), index=False)
    summ = {"label": "POST HOC (amendment A9)", "level_pct_of_reachable": LEVEL, "budgets": BUDGETS,
            "onset_sample_min": None if a.all_onsets else SAMPLE_MIN, "windows": info, "budget": {}, "failsignal": {}}
    for (s, k), g in B.groupby(["season", "k"]):
        v = g.verdict.value_counts().to_dict()
        late = g[g.ecc_excess_min.notna() & (g.ecc_excess_min > 1e-9)]
        summ["budget"][f"{s}_k{k}"] = {"onsets": int(len(g)), "verdicts": {kk: int(vv) for kk, vv in v.items()},
                                       "ecc_excess_min_max": float(late.ecc_excess_min.max()) if len(late) else 0.0,
                                       "unknown": int(g.unknown.sum())}
    for c in ("dT_L_h", "openloop_t_L_h", "adaptive_t_L_h"):
        F[c] = pd.to_numeric(F.get(c), errors="coerce")
    # A realisation in which the adaptive plan never reaches the open-loop plan's final coverage (it ends with
    # less) is the worst outcome, not a missing one: dT = -inf, so it ranks below every finite delay. (Revision
    # of 30 Oct: the first version dropped these, which flattered the Fail-Signal.)
    F["dT_rank_h"] = np.where(F.dT_L_h.notna(), F.dT_L_h, np.where(F.openloop_t_L_h.notna(), -np.inf, np.nan))
    for s, g in F.groupby("season"):
        x = g.dT_rank_h.dropna().to_numpy(dtype=float)
        per_w = g.groupby("window").dT_rank_h.median()
        q = np.quantile(x, [0.25, 0.5, 0.75]) if len(x) else [None] * 3
        fin = lambda v: None if v is None or not np.isfinite(v) else float(v)
        summ["failsignal"][s] = {"realisations": int(len(g)), "with_dT": int(len(x)),
                                 "adaptive_never_reaches_openloop_level": int(np.isinf(x).sum()),
                                 "dT_L_h_median": fin(q[1]), "dT_L_h_iqr": [fin(q[0]), fin(q[2])],
                                 "dT_L_h_iqr_lower_is_never": bool(len(x) and not np.isfinite(q[0])),
                                 "pct_adaptive_sooner": float(100 * (x > 1e-9).mean()) if len(x) else None,
                                 "pct_same_time": float(100 * (np.abs(x) <= 1e-9).mean()) if len(x) else None,
                                 "pct_adaptive_later_or_never": float(100 * (x < -1e-9).mean()) if len(x) else None,
                                 "windows_median_dT_positive": int((per_w > 1e-9).sum()), "windows": int(per_w.notna().sum()),
                                 "adaptive_reach_pct_mean": float(g.adaptive_reach_pct.mean()),
                                 "openloop_reach_pct_mean": float(g.openloop_reach_pct.mean())}
    json.dump(summ, open(os.path.join(out, "a9_summary.json"), "w"), indent=1, default=str)
    write_tables(summ, out)
    print(json.dumps({k: summ[k] for k in ("budget", "failsignal")}, indent=1, default=str))
    print(f"Files: {out}/a9_budget.csv, a9_failsignal.csv, a9_summary.json, tab_a9_budget.tex, tab_a9_failsignal.tex")
    return B, F, summ


def write_tables(summ, out):
    V = ("ECC at floor", "ECC at exact time", "ECC delay avoidable", "ECC later than exact", "ECC missed a cover",
         "no k-cover exists")
    L = ["\\begin{table*}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         "\\setlength{\\tabcolsep}{3.5pt}",
         "\\caption{Post hoc (amendment A9): ECC at a tight budget $k$ against the exact earliest time to 95\\% of the "
         "reachable tiles (HiGHS feasibility over the delivered images, bisection over delivery times), sampled "
         "onsets. Columns count onsets by verdict; Excess is ECC's largest delay beyond the exact time (min).}",
         "\\label{tab:a9_budget}", "\\begin{tabular}{|l|c|c|c|c|c|c|c|c|c|}", "\\hline",
         "Season, $k$ & Onsets & At floor & At exact & Delay avoidable & Later & Missed cover & No cover & "
         "Unknown & Excess \\\\", "\\hline"]
    for key, s in summ["budget"].items():
        v = s["verdicts"]
        L += [f"{key.replace('_k', ', k = ')} & {s['onsets']} & " + " & ".join(str(v.get(x, 0)) for x in V)
              + f" & {s['unknown']} & {s['ecc_excess_min_max']:.0f} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table*}"]
    open(os.path.join(out, "tab_a9_budget.tex"), "w").write("\n".join(L) + "\n")
    f = lambda x: "--" if x is None else f"{x:.2f}"
    g0 = lambda x: "--" if x is None else f"{x:.0f}"
    L = ["\\begin{table}[!t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.2}",
         "\\setlength{\\tabcolsep}{3pt}",
         "\\caption{Post hoc (amendment A9): the Fail-Signal in hours. $\\Delta T$: open-loop plan's time to its own "
         "final coverage minus the adaptive plan's time to the same coverage, same number of images, per cloud "
         "realisation ($p = 0.35$, correlated); median and interquartile range in hours; sooner/same/later in \\% "
         "of realisations. A realisation in which the "
         "adaptive plan ends with less coverage counts as never reaching it (later). Win.: windows whose median "
         "is positive.}",
         "\\label{tab:a9_failsignal}", "\\begin{tabular}{|l|c|c|c|c|c|}", "\\hline",
         "Season & $n$ & $\\Delta T$ & IQR & Sooner/same/later & Win.\\ $>0$ \\\\", "\\hline"]
    for s, v in summ["failsignal"].items():
        iqr = v["dT_L_h_iqr"] or [None, None]
        lo = "never" if v.get("dT_L_h_iqr_lower_is_never") else f(iqr[0])
        L += [f"{s} & {v['with_dT']} & {f(v['dT_L_h_median'])} & {lo}, {f(iqr[1])} & "
              f"{g0(v['pct_adaptive_sooner'])}/{g0(v['pct_same_time'])}/{g0(v['pct_adaptive_later_or_never'])} & "
              f"{v['windows_median_dT_positive']} of {v['windows']} \\\\", "\\hline"]
    L += ["\\end{tabular}", "\\end{table}"]
    open(os.path.join(out, "tab_a9_failsignal.tex"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
