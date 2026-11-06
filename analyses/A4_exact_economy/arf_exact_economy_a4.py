# =============================================================================
# ARF EXACT IMAGE ECONOMY (protocol amendment A4, POST HOC; read-only)
# =============================================================================
# H1b asks ECC to reach the H1a level with FEWER images than first-come-first-
# served (FCFS). At an onset where ECC finishes at the physical floor, no
# selector can finish earlier, and the fewest images any selector could use to
# finish then is the exact minimum number of images, among those delivered by
# that time, that cover the level. This script measures that minimum at every
# onset of every primary window (clear sky), so each H1b outcome can be read as
#   ECC images = exact minimum  -> ECC is optimal at that onset
#   exact minimum = FCFS images -> no selector could beat FCFS there: an H1b
#                                  "tie" is then a property of the passes, not
#                                  of the algorithm
# Procedure: for each primary run folder, Cells 6-8 of release 1.5.3 are
# replayed on the run's own saved inputs (arf_refined_events.csv and
# arf_mosaic_tiles.csv, checked against the SHA-256 recorded by Cell 9) by
# arf_colab_parity.replay, the procedure that reproduced the Colab runs of 15
# October 2026 value for value, with ONE change: the fleet's median altitude and
# orbital period are set to the run's own values (arf_cell2_summary.json:
# median_alt_km, t_max_wait_min). arf_colab_parity.replay rebuilds both from the
# rounded period, which moves the period by about 1e-13 min; that is invisible
# at the parity check's tolerance (1e-9) but enough to change the last digits
# of the time-weighted columns, so the first version of this script (byte
# check) stopped on window 1 in Colab, 26 October 2026. The replayed
# arf_onset_robustness.csv is then compared with the run's: it must be
# byte-identical, or else equal in every column (exactly in the columns this
# analysis uses: floor, ECC and FCFS times and images, reachable tiles; to a
# relative 1e-9 in the others); otherwise the script stops and prints, per
# column, how far apart the two are. Then, per onset:
# the onset's passes (subset_by_time), ECC (earliest_completion, k = 30) at
# 95% of reachable tiles, the events delivered by ECC's completion time, and
# the exact minimum images covering the level from them (solve_ilp_min_images,
# HiGHS). All functions are the release's own; nothing is reimplemented
# except the selection of delivered events.
# Writes <BATCH_DIR>/sensitivity_a4/a4_exact_economy.csv and a summary JSON.
# Exact solves (about 13 s each on the 12 October data): at every onset where
# ECC does not use fewer images than FCFS (each H1b tie or loss), and at the
# onsets 0, 80, 160, ... min (a fixed sample of 19 for ECC's optimality);
# --all-onsets solves all 73. Colab time: Cells 6-8 per window (about 5 min)
# plus the solves, about 10 min per window, 1.5 h for eight.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_exact_economy_a4.py
# =============================================================================
import os
import json
import time
import argparse
import importlib.util

import numpy as np
import pandas as pd

LEVEL = 95.0


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


EXACT_COLUMNS = ("onset_shift_min", "reachable_clear_tiles", "best_possible_clear_t95r_min",
                 "earliest_completion_clear_t95r_min", "earliest_completion_clear_images95r",
                 "fcfs_unbudgeted_clear_t95r_min", "fcfs_unbudgeted_clear_images95r")


def replay_exact(par, alt_km, t_orbit_min, *args):
    """arf_colab_parity.replay with the run's own CONSTELLATION_ALTITUDE_KM and T_MAX_WAIT_MINUTES."""
    import inspect
    src = inspect.getsource(par.replay)
    hook = '    ns["h3"] = h3\n'
    assert src.count(hook) == 1 and src.count("def replay(pipeline_path,") == 1, "unexpected arf_colab_parity.replay"
    src = src.replace(hook, hook + '    ns.update(CONSTELLATION_ALTITUDE_KM=_ALT, T_MAX_WAIT_MINUTES=_TORB)\n')
    src = src.replace("def replay(pipeline_path,", "def replay_exact_(_ALT, _TORB, pipeline_path,")
    g = dict(vars(par))
    exec(compile(src, par.__file__, "exec"), g)
    return g["replay_exact_"](alt_km, t_orbit_min, *args)


def compare_onsets(a_path, b_path, rtol=1e-9):
    """(identical bytes, acceptable, report): report lists the columns that differ and by how much."""
    import hashlib
    if hashlib.sha256(open(a_path, "rb").read()).digest() == hashlib.sha256(open(b_path, "rb").read()).digest():
        return True, True, []
    a, b = pd.read_csv(a_path), pd.read_csv(b_path)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        return False, False, [f"shape or columns differ: {a.shape} {b.shape}"]
    rep, ok = [], True
    for c in a.columns:
        x, y = pd.to_numeric(a[c], errors="coerce"), pd.to_numeric(b[c], errors="coerce")
        if x.notna().any() or y.notna().any():
            na = (x.isna() != y.isna())
            both = x.notna() & y.notna()
            dx = (x[both] - y[both]).abs()
            rel = (dx / np.maximum(y[both].abs(), 1e-300)).max() if both.any() else 0.0
            if na.any() or (dx > 0).any():
                exact = c in EXACT_COLUMNS
                bad = na.any() or (exact and (dx > 0).any()) or (not exact and rel > rtol)
                ok &= not bad
                rep.append(f"{c}: {int(na.sum() + (dx > 0).sum())} rows differ, max |diff| {dx.max() if len(dx) else 0:.3g}, "
                           f"max rel {rel:.3g}{' (must be exact)' if exact else ''}{' FAIL' if bad else ''}")
        elif not a[c].fillna("").astype(str).equals(b[c].fillna("").astype(str)):
            ok = False
            rep.append(f"{c}: text values differ FAIL")
    return False, ok, rep


def delivered_by(ns, sub, t):
    """EventStore of the events of `sub` delivered by time t (same rows, utilities and times)."""
    d = ns["event_delivery_times"](sub)
    keep = [e for e in range(len(sub)) if d[e] <= t + 1e-9]
    ES = ns["EventStore"]
    if not keep:
        return ES([], np.zeros(1, np.int64), np.zeros(0, np.int32), np.zeros(0), [])
    rows = np.concatenate([np.arange(sub.offsets[e], sub.offsets[e + 1]) for e in keep])
    counts = sub.counts[keep]
    return ES([sub.keys_[e] for e in keep], np.r_[0, np.cumsum(counts)], sub.tile_idx[rows], sub.util[rows],
              [sub.wait_lookup[sub.keys_[e]] for e in keep])


def onset_exact(ns, time_limit_s=300.0, max_pairs=3_000_000, every_min=80, all_onsets=False, label=""):
    """Per onset (clear sky): floor, ECC time and images, FCFS images; the exact minimum images by ECC's
    time at every onset where ECC does not use fewer images than FCFS (each H1b tie or loss), and at the
    onsets whose shift is a multiple of `every_min` (a fixed sample for ECC's optimality), or at all
    onsets with all_onsets=True. About 13 s per exact solve on the 12 October 2026 CKGR+Okavango data."""
    H = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
    k = int(ns["MAX_SELECTION_STEPS"])
    out = []
    for d in range(0, ns["ONSET_SHIFT_HOURS"] * 60 + 1, ns["ONSET_SHIFT_STEP_MIN"]):
        sub, _ = ns["subset_by_time"](ns["event_data"], None, d, d + H)
        R = ns["reachable_tile_count"](sub, None) if len(sub) else 0
        row = {"onset_shift_min": d, "reachable_tiles": R}
        if not R:
            out.append(row)
            continue
        need = int(np.ceil(LEVEL / 100.0 * R - 1e-9))
        floor = ns["coverage_times"](sub, None, (LEVEL,), denominator=R)[LEVEL]
        t_e, _, n_e = ns["earliest_completion"](sub, None, k, LEVEL, denominator=R)
        _, n_f, _ = ns["level_stats"](sub, None, ns["select_fcfs"](sub, None, len(sub))[0], need)
        row.update(floor_min=floor, ecc_min=t_e, ecc_images=n_e if t_e is not None else None, fcfs_images=n_f)
        tie_or_loss = t_e is not None and n_f is not None and n_e >= n_f
        row["solved_because"] = ("tie or loss" if tie_or_loss else "sample" if d % every_min == 0 or all_onsets
                                 else None)
        if t_e is not None and row["solved_because"]:
            dl = delivered_by(ns, sub, t_e)
            t0 = time.time()
            ex = ns["solve_ilp_min_images"](dl, None, need, time_limit_s, max_pairs)
            row.update(exact_images=ex["images"], exact_status=ex["status"], exact_pairs=ex["n_pairs"],
                       exact_seconds=round(time.time() - t0, 2), delivered_events=len(dl))
            print(f"  {label} onset +{d} min: ECC {n_e}, FCFS {n_f}, exact {ex['images']} ({ex['status']}, "
                  f"{row['exact_seconds']} s; {row['solved_because']})", flush=True)
        out.append(row)
    df = pd.DataFrame(out)
    ok = df.get("exact_status", pd.Series(dtype=str)) == "optimal"
    df["ecc_optimal"] = np.where(ok, df.get("ecc_images") == df.get("exact_images"), None)
    df["no_economy_possible"] = np.where(ok, df.get("exact_images") >= df.get("fcfs_images"), None)
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description="Post-hoc exact image economy at every onset (A4).")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    ap.add_argument("--index", default="confirmatory_v153_primary_index.csv")
    ap.add_argument("--all-onsets", action="store_true", help="exact solve at every onset (about 3 h in Colab)")
    a, _ = ap.parse_known_args(argv)
    par = load(os.path.join(a.code, "arf_colab_parity.py"), "arf_colab_parity")
    pipe = os.path.join(a.code, "algorithmic_reconstruction_framework.py")
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    assert par.sha256_file(pipe) == proto["pipeline_sha256"], "ARF_code does not hold the release of record"
    idx = pd.read_csv(os.path.join(a.batch_dir, a.index), dtype=str).fillna("")
    frames, runs = [], []
    for _, r in idx[(idx.status == "done") & (idx.release == proto["release"])].iterrows():
        d = r["out_dir"]
        man = json.load(open(os.path.join(d, "arf_reproducibility_manifest.json")))
        cfg, inp = man.get("config_overrides") or {}, man["inputs"]
        if cfg.get("analysis_epoch") not in proto["windows"]:
            continue
        ov = man["output_validation"]["files"]
        refined, tiles_csv = os.path.join(d, "arf_refined_events.csv"), os.path.join(d, "arf_mosaic_tiles.csv")
        for p in (refined, tiles_csv):
            assert par.sha256_file(p) == ov[os.path.basename(p)]["sha256"], f"{p} differs from what Cell 9 recorded"
        key = [k for k in man["scenario_results"]][0]
        orbit = man["scenario_results"][key]["orbital_period_min"]
        c2 = json.load(open(os.path.join(d, "arf_cell2_summary.json")))
        alt, torb = float(c2["median_alt_km"]), float(c2["t_max_wait_min"])
        assert torb == float(orbit), f"{d}: Cell 2's period {torb!r} differs from the manifest's {orbit!r}"
        tiles = list(pd.read_csv(tiles_csv, usecols=["tile"])["tile"])
        win = str(cfg["analysis_epoch"])[:10]
        t0 = time.time()
        ns = replay_exact(par, alt, torb, pipe, refined, tiles, orbit, inp["active_roi"], inp["fleet"],
                          inp["config"].get("planet_swath_case", "conservative"), inp["n_satellites"])
        ident, acceptable, rep = compare_onsets(os.path.join(ns["_out_dir"], "arf_onset_robustness.csv"),
                                                os.path.join(d, "arf_onset_robustness.csv"))
        print(f"{win}: replay {'byte-identical to' if ident else 'differs from'} the run's arf_onset_robustness.csv"
              + ("" if ident else ":\n    " + "\n    ".join(rep)), flush=True)
        assert acceptable, f"{d}: the replay does not reproduce the run's arf_onset_robustness.csv (see above)"
        df = onset_exact(ns, all_onsets=a.all_onsets, label=win)
        df.insert(0, "window", win)
        frames.append(df)
        e = df[df.get("exact_status") == "optimal"] if "exact_status" in df else df.iloc[:0]
        elig = e[e.fcfs_images >= 2]
        ties = elig[elig.solved_because == "tie or loss"]
        runs.append({"window": win, "replay_identical": ident, "replay_differences": rep,
                     "seconds": round(time.time() - t0, 1),
                     "onsets_solved": int(len(e)), "ecc_optimal": int(e.ecc_optimal.sum()),
                     "ecc_extra_images_max": float((e.ecc_images - e.exact_images).max()) if len(e) else None,
                     "onsets_tie_or_loss": int(((df.ecc_images >= df.fcfs_images) & (df.fcfs_images >= 2)).sum()),
                     "tie_or_loss_solved": int(len(ties)),
                     "tie_or_loss_no_economy_possible": int(ties.no_economy_possible.sum()),
                     "tie_or_loss_ecc_above_exact": int((ties.ecc_images > ties.exact_images).sum()),
                     "sample_solved": int((e.solved_because == "sample").sum()),
                     "sample_ecc_optimal": int(e[e.solved_because == "sample"].ecc_optimal.sum())})
        print(json.dumps(runs[-1]), flush=True)
    out = os.path.join(a.batch_dir, "sensitivity_a4")
    os.makedirs(out, exist_ok=True)
    all_df = pd.concat(frames) if frames else pd.DataFrame()
    all_df.to_csv(os.path.join(out, "a4_exact_economy.csv"), index=False)
    summ = {"label": "POST HOC (amendment A4)", "level_pct_of_reachable": LEVEL, "windows": runs}
    json.dump(summ, open(os.path.join(out, "a4_exact_economy_summary.json"), "w"), indent=1)
    print(pd.DataFrame(runs).to_string(index=False))
    print(f"Files: {out}/a4_exact_economy.csv, a4_exact_economy_summary.json")
    return all_df, summ


if __name__ == "__main__":
    main()
