# =============================================================================
# ARF EO FOLLOW-UPS (v1.5.4 package): three read-only analyses of the real-fleet
# protocol runs, labelled POST HOC (not part of the pre-registered decisions)
# =============================================================================
# A. H1a gaps: per primary window, how far ECC's time to 95% of reachable tiles
#    lies above the physical floor (minutes), against the H1a margin of one
#    orbital period. Shows whether H1a's 100% is a near thing or a wide one.
# B. H1b failures: per primary window, the eligible onsets (FCFS needs >= 2
#    images) split into: success (ECC within the margin AND fewer images),
#    tie (same images), ECC more images, ECC late (beyond the margin), ECC
#    missing; with the onsets that fail, listed for the window(s) under 90%.
# C. Real fleet against the dense arm: per window, the median over onsets of
#    the physical floor and of ECC (95% of reachable tiles, clear sky), beside
#    the dense arm's medians for the same window (dense_summary/dense_runs.csv).
# Every rule is the pipeline's own (Cell 7C of release 1.5.3: _h1_revised,
# _h1b); each window's recomputed H1a and H1b shares are checked against the
# run's arf_onset_robustness_summary.json and the script stops on a mismatch.
# Writes CSV files to <BATCH_DIR>/eo_followups/; reads only.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_eo_followups.py
# =============================================================================
import os
import json
import argparse

import numpy as np
import pandas as pd


def onset_tables(run_dir):
    on = pd.read_csv(os.path.join(run_dir, "arf_onset_robustness.csv"))
    sm = json.load(open(os.path.join(run_dir, "arf_onset_robustness_summary.json")))
    man = json.load(open(os.path.join(run_dir, "arf_reproducibility_manifest.json")))
    return on, sm, man


def num(v):
    return None if v is None or pd.isna(v) else float(v)


def window_analysis(run_dir, sky="clear"):
    on, sm, man = onset_tables(run_dir)
    h1 = sm["h1_revised"]
    lvl, margin = float(h1["reachable_coverage_pct"]), float(h1["floor_margin_min"])
    L = f"{lvl:g}"
    cfg = man.get("config_overrides") or {}
    gaps, n_eval_a, n_ok_a = [], 0, 0
    b = {"eligible": 0, "success": 0, "tie": 0, "ecc_more_images": 0, "ecc_late": 0, "ecc_missing": 0, "one_image": 0}
    fails = []
    for _, r in on.iterrows():
        if not r.get(f"reachable_{sky}_tiles"):
            continue
        fl = num(r.get(f"best_possible_{sky}_t{L}r_min"))
        if fl is None:
            continue
        te = num(r.get(f"earliest_completion_{sky}_t{L}r_min"))
        # A: the pipeline's _h1_revised rule
        n_eval_a += 1
        if te is not None:
            gaps.append(te - fl)
            n_ok_a += int(te - fl <= margin + 1e-9)
        # B: the pipeline's _h1b rule
        nf = num(r.get(f"fcfs_unbudgeted_{sky}_images{L}r"))
        if nf is None:
            continue
        if nf < 2:
            b["one_image"] += 1
            continue
        b["eligible"] += 1
        ne = num(r.get(f"earliest_completion_{sky}_images{L}r"))
        if te is None or ne is None:
            b["ecc_missing"] += 1
            kind = "ECC missing"
        elif te - fl > margin + 1e-9:
            b["ecc_late"] += 1
            kind = "ECC late"
        elif ne < nf:
            b["success"] += 1
            continue
        elif ne == nf:
            b["tie"] += 1
            kind = "tie"
        else:
            b["ecc_more_images"] += 1
            kind = "ECC more images"
        fails.append({"onset_shift_min": r.get("onset_shift_min"), "kind": kind, "floor_min": fl, "ecc_min": te,
                      "ecc_images": ne, "fcfs_images": nf})
    # self-check against the pipeline's own summary
    ref_a = h1["series"][f"earliest_completion_{sky}"]
    ref_b = sm["h1b_image_economy"]["series"][sky]
    pa = 100.0 * n_ok_a / n_eval_a if n_eval_a else None
    pb = 100.0 * b["success"] / b["eligible"] if b["eligible"] else None
    for mine, ref, what in ((n_eval_a, ref_a["n_onsets_evaluated"], "H1a onsets"), (pa, ref_a["pct_within_margin"], "H1a %"),
                            (b["eligible"], ref_b["n_onsets_evaluated"], "H1b onsets"), (pb, ref_b["pct_success"], "H1b %")):
        assert (mine is None and ref is None) or (mine is not None and ref is not None and abs(mine - ref) < 1e-9), \
            f"{run_dir}: recomputed {what} {mine} differs from the run's summary {ref}"
    g = np.asarray(gaps, float)
    fl_all = pd.to_numeric(on.get(f"best_possible_{sky}_t{L}r_min"), errors="coerce").dropna()
    ecc_all = pd.to_numeric(on.get(f"earliest_completion_{sky}_t{L}r_min"), errors="coerce").dropna()
    row = {"window": str(cfg.get("analysis_epoch") or "")[:10], "roi": cfg.get("active_roi"), "fleet": cfg.get("fleet"),
           "case": cfg.get("planet_swath_case"), "n_satellites": (man.get("inputs") or {}).get("n_satellites"),
           "margin_min": round(margin, 1), "h1a_onsets": n_eval_a, "h1a_pct": pa,
           "gap_median_min": float(np.median(g)) if g.size else None,
           "gap_p90_min": float(np.percentile(g, 90)) if g.size else None,
           "gap_max_min": float(g.max()) if g.size else None,
           "gap_max_share_of_margin_pct": float(100.0 * g.max() / margin) if g.size and margin else None,
           "h1b_eligible": b["eligible"], "h1b_pct": pb, **{f"h1b_{k}": v for k, v in b.items() if k != "eligible"},
           "fcfs_images_median": float(pd.to_numeric(on.get(f"fcfs_unbudgeted_{sky}_images{L}r"), errors="coerce").median()),
           "ecc_images_median": float(pd.to_numeric(on.get(f"earliest_completion_{sky}_images{L}r"), errors="coerce").median()),
           "floor_median_min": float(fl_all.median()) if fl_all.size else None,
           "ecc_median_min": float(ecc_all.median()) if ecc_all.size else None,
           "candidate_events_median": float(pd.to_numeric(on.get("n_events"), errors="coerce").median()),
           "run_dir": run_dir}
    return row, fails


def main(argv=None):
    ap = argparse.ArgumentParser(description="Post-hoc, read-only follow-ups of the real-fleet protocol runs.")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    ap.add_argument("--index", default="confirmatory_v153_primary_index.csv")
    ap.add_argument("--h1b-threshold", type=float, default=90.0)
    a, _ = ap.parse_known_args(argv)
    idx = pd.read_csv(os.path.join(a.batch_dir, a.index), dtype=str).fillna("")
    rows, fails = [], {}
    for _, r in idx[idx.status == "done"].iterrows():
        row, f = window_analysis(r["out_dir"])
        rows.append(row)
        fails[row["run_dir"]] = f
    df = pd.DataFrame(rows).sort_values("window")
    out = os.path.join(a.batch_dir, "eo_followups")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "eo_windows.csv"), index=False)
    pd.set_option("display.width", 220)
    print(f"POST-HOC follow-ups, {len(df)} runs of {a.index} (every window's H1a and H1b shares re-derived and "
          "matched to the run's own summary)\n")
    print("A. H1a: ECC minus physical floor, minutes, at 95% of reachable tiles (clear sky)")
    print(df[["window", "n_satellites", "margin_min", "h1a_onsets", "h1a_pct", "gap_median_min", "gap_p90_min",
              "gap_max_min", "gap_max_share_of_margin_pct"]].round(1).to_string(index=False))
    print("\nB. H1b: eligible onsets by outcome (success = ECC within margin AND fewer images than FCFS)")
    print(df[["window", "h1b_eligible", "h1b_pct", "h1b_success", "h1b_tie", "h1b_ecc_more_images", "h1b_ecc_late",
              "h1b_ecc_missing", "h1b_one_image", "fcfs_images_median", "ecc_images_median"]].round(1).to_string(index=False))
    low = df[df.h1b_pct.notna() & (df.h1b_pct < a.h1b_threshold)]
    rows_f = []
    for _, w in low.iterrows():
        f = pd.DataFrame(fails[w["run_dir"]])
        rows_f.append(f.assign(window=w["window"], run_dir=w["run_dir"]))
        print(f"\n   Window {w['window'] or os.path.basename(w['run_dir'])}: {len(f)} failing onsets of {w['h1b_eligible']}; by kind "
              f"{f['kind'].value_counts().to_dict()}; onset shifts (min) {f['onset_shift_min'].min():g} to "
              f"{f['onset_shift_min'].max():g}; images ECC/FCFS (median) {f['ecc_images'].median():g}/{f['fcfs_images'].median():g}")
    if rows_f:
        pd.concat(rows_f).to_csv(os.path.join(out, "eo_h1b_failing_onsets.csv"), index=False)
    print("\nC. Median time to 95% of reachable tiles (clear sky), minutes: real fleet against the dense arm")
    dp = os.path.join(a.batch_dir, "dense_summary", "dense_runs.csv")
    comp = df[["window", "n_satellites", "floor_median_min", "ecc_median_min"]].rename(
        columns={"n_satellites": "real_fleet_satellites", "floor_median_min": "real_floor", "ecc_median_min": "real_ecc"})
    if os.path.exists(dp):
        dn = pd.read_csv(dp)
        piv = dn.pivot_table(index="window", columns="N", values="floor_t95r_median_min", aggfunc="first")
        piv.columns = [f"dense_floor_N{int(c)}" for c in piv.columns]
        comp = comp.merge(piv.reset_index(), on="window", how="left")
    print(comp.round(1).to_string(index=False))
    comp.to_csv(os.path.join(out, "eo_vs_dense_floor.csv"), index=False)
    print(f"\nFiles: {out}/eo_windows.csv, eo_h1b_failing_onsets.csv, eo_vs_dense_floor.csv")
    return df


if __name__ == "__main__":
    main()
