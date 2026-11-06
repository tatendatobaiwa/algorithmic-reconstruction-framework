"""Cross-run summary of the v1.5.2 replays (same 40 columns as
arf_cross_run_summary_v150_replay.csv, plus refined rows, candidate events and
the input used), and a v1.5.0 -> v1.5.2 comparison of the key results.
Usage: make_cross_run_v152.py REPLAY_JSONL V150_CSV OUT_CSV [DELTA_CSV]
Each run's statuses come from its own arf_hypothesis_validation.json (the
record's out_dir; for records without one, the single output folder of the
same scenario label that no other record claims)."""
import sys, os, json, glob
import pandas as pd

jl, v150_csv, out_csv = sys.argv[1:4]
delta_csv = sys.argv[4] if len(sys.argv) > 4 else None
ROI = {"ckgr": ("CKGR + Okavango", "ckgr_okavango"), "karoo": ("Karoo NR", "karoo_nature_reserve"),
       "lz": ("Lower Zambezi NP", "lower_zambezi")}
FLEET = {"open": ("Landsat 8/9 + Sentinel-2", "earth-resources"),
         "openPC": ("open data + PlanetScope (conservative)", "earth-resources_planet-conservative"),
         "openPO": ("open data + PlanetScope (optimistic)", "earth-resources_planet-optimistic"),
         "PC": ("PlanetScope (conservative)", "planet-conservative"),
         "PO": ("PlanetScope (optimistic)", "planet-optimistic")}
ORDER = [f"{r}_{f}" for r in ("ckgr", "karoo", "lz") for f in ("open", "openPC", "openPO", "PC", "PO")]

recs = {}
for line in open(jl):
    r = json.loads(line)
    recs[r["run"]] = r                                  # the latest record of a run wins
# output folders recorded by any replay log next to this one (e.g. the runs
# without Sentinel-2A share scenario labels with the full-fleet runs)
claimed = set()
for other in glob.glob(os.path.join(os.path.dirname(os.path.abspath(jl)), "replay152_*.jsonl")):
    for line in open(other):
        o = json.loads(line).get("out_dir")
        if o:
            claimed.add(o)
by_label = {}
for d in glob.glob("/tmp/R150_OUT_*"):
    f = os.path.join(d, "arf_hypothesis_validation.json")
    if os.path.exists(f) and d not in claimed:
        lab = list(json.load(open(f))["H1_revised"]["assessments"])[0]
        by_label.setdefault(lab, []).append(d)


def hv_of(run, rec):
    d = rec.get("out_dir")
    if not d:
        roi, fl = run.split("_", 1)
        lab = f"ROI={ROI[roi][1]}, fleet={FLEET[fl][1]}"
        ds = by_label.get(lab, [])
        if len(ds) != 1:
            raise SystemExit(f"{run}: {len(ds)} unclaimed output folders for {lab}; record its out_dir")
        d = ds[0]
    hv = json.load(open(os.path.join(d, "arf_hypothesis_validation.json")))
    lab = list(hv["H1_revised"]["assessments"])[0]
    return {k: (v["assessments"][lab] if isinstance(v, dict) and "assessments" in v else v) for k, v in hv.items()}, d


def g(d, *ks):
    for k in ks:
        d = d.get(k) if isinstance(d, dict) else None
    return d


rows = []
for run in ORDER:
    if run not in recs:
        continue
    r = recs[run]
    if "error" in r:
        rows.append({"ROI": ROI[run.split("_")[0]][0], "Fleet": FLEET[run.split("_", 1)[1]][0],
                     "Run": run, "Replay error": r["error"].strip().splitlines()[-1]})
        continue
    hv, d = hv_of(run, r)
    roi, fl = run.split("_", 1)
    h1a, h1b, h2, cf, il = r["h1a"] or {}, r["h1b"] or {}, r["h2_ta"] or {}, r["h2_cf"] or {}, r["ilp"] or {}
    dif = h2.get("difference_vs_openloop_expected") or {}
    cm = h2.get("coverage_pct_mean") or {}
    orig = hv.get("H2_original_ratio") or {}
    mp = orig.get("coverage_maintenance_pct")
    rows.append({
        "ROI": ROI[roi][0], "Fleet": FLEET[fl][0],
        "Drive folder (12 Oct run)": f"ARF_Outputs_{ROI[roi][1]}_{FLEET[fl][1]}" + (
            f" (replayed without the {', '.join(r['dropped_sats'])} rows)" if r.get("dropped_sats") else ""),
        "Decay time scale tau (h)": round(r["tau_min"] / 60.0, 2) if r.get("tau_min") is not None else None,
        "Reachable usable tiles, clear (%)": r.get("reach_pct"),
        "Coverage-first coverage, clear (%)": r.get("cov_clear"),
        "Time-aware coverage, clear (%)": r.get("ta_cov_clear"),
        "Coverage-first coverage, cloud seed 42 (%)": r.get("cov_cloud"),
        "Time-aware coverage, cloud seed 42 (%)": r.get("ta_cov_cloud"),
        "H1 (H1a and H1b)": hv["H1_reframed"]["status"],
        "H1a status (ECC)": hv["H1_revised"]["status"],
        "H1a onsets": h1a.get("n_onsets_evaluated"),
        "H1a ECC within margin (%)": h1a.get("pct_within_margin"),
        "H1a ECC 95% lower bound (%)": h1a.get("pct_within_margin_lcb95"),
        "H1a status, time-aware selector": hv["H1_revised"].get("status_time_aware"),
        "H1b status": hv["H1b_image_economy"]["status"],
        "H1b eligible onsets": h1b.get("n_onsets_evaluated"),
        "H1b ECC fewer images within margin (%)": h1b.get("pct_success"),
        "H1b 95% lower bound (%)": h1b.get("pct_success_lcb95"),
        "Median images to level, ECC": h1b.get("median_images_ecc"),
        "Median images to level, FCFS": h1b.get("median_images_fcfs"),
        "H2 status (time-aware; adaptive vs open-loop, same images)": hv["H2"]["status"],
        "H2 mean difference (coverage points)": dif.get("mean"),
        "H2 one-sided 95% lower bound": dif.get("lcb"),
        "H2 cloud realisations": dif.get("n"),
        "H2 mean images (adaptive)": h2.get("images_mean"),
        "H2 mean coverage adaptive (%)": cm.get("adaptive"),
        "H2 mean coverage open-loop expected (%)": cm.get("openloop_expected"),
        "H2 mean coverage open-loop clear-sky (%)": cm.get("openloop_clear"),
        "H2 mean coverage greedy with cloud known (%)": cm.get("oracle"),
        "H2 not evaluable because": "; ".join(h2.get("not_evaluable_reasons") or []),
        "H2 status, coverage-first": hv["H2"].get("status_coverage_first"),
        "H2 original ratio test (>50% of clear)": (None if mp is None else
                                                   ("SUPPORTED" if float(mp) > 50.0 else "NOT SUPPORTED")),
        "H3 status": hv["H3"]["status"],
        "ILP instances solved (time-weighted)": il.get("tw_solved"),
        "Greedy/optimum, time-weighted objective (worst)": il.get("time_aware_tw_ratio_min"),
        "Greedy/optimum, coverage objective, coverage-first (worst)": il.get("coverage_first_ratio_min"),
        "Coverage objective ratio, time-aware selector (worst)": il.get("time_aware_ratio_min"),
        "Coverage objective ratio, FCFS (worst)": il.get("fcfs_ratio_min"),
        "ECC images = exact minimum (windows)": (f"{il.get('ecc_matches_exact')}/{il.get('ecc_instances_exact')}"
                                                 if il.get("ecc_instances_exact") is not None else None),
        "Refined rows": r.get("n_rows"), "Candidate events": r.get("n_events"),
        "Input (refined events)": ("v1.5.2 Cell 5 with real Orekit on the saved elements" if "exact152" in
                                   str(r.get("refined_dir")) and not r.get("dedupe") else
                                   "Colab v1.5.0 Cell 5 output with the end-cap rule applied"),
        "Satellites dropped": ",".join(r.get("dropped_sats") or []),
        "Run": run,
    })
df = pd.DataFrame(rows)
for c in df.columns:
    if df[c].dtype == float:
        df[c] = df[c].round(4)
df.to_csv(out_csv, index=False)
print(f"{len(df)} runs -> {out_csv}")

if delta_csv:
    old = pd.read_csv(v150_csv)
    old_txt = pd.read_csv(v150_csv, dtype=str)
    key = ["ROI", "Fleet"]

    def _decimals(col):
        """Decimals the v1.5.0 table used for a column (its values were rounded)."""
        s = old_txt[col].dropna()
        return max((len(v.split(".")[1]) if "." in v else 0) for v in s) if len(s) else 6
    cols = ["Reachable usable tiles, clear (%)", "Time-aware coverage, clear (%)", "H1 (H1a and H1b)",
            "H1a status (ECC)", "H1a ECC within margin (%)", "H1b status", "H1b eligible onsets",
            "H1b ECC fewer images within margin (%)", "Median images to level, ECC", "Median images to level, FCFS",
            "H2 status (time-aware; adaptive vs open-loop, same images)", "H2 mean difference (coverage points)",
            "H2 one-sided 95% lower bound", "H3 status", "Decay time scale tau (h)"]
    m = old[key + cols].merge(df[key + cols + ["Refined rows", "Candidate events"]], on=key,
                              suffixes=(" v1.5.0", " v1.5.2"), how="outer")
    out = m[key].copy()
    # rows of the saved Colab v1.5.0 inputs (12 October 2026), for the end-cap effect
    inv = {(ROI[k.split("_")[0]][0], FLEET[k.split("_", 1)[1]][0]): k for k in ORDER}
    out["Refined rows v1.5.0 (Colab)"] = [
        sum(1 for _ in open(f"/home/claude/realdata/{inv[(a, b)]}_refined_events.csv")) - 1
        for a, b in zip(m.ROI, m.Fleet)]
    out["Refined rows v1.5.2"] = m["Refined rows"]
    out["Candidate events v1.5.2"] = m["Candidate events"]
    changed = []
    for c in cols:
        a, b = m[c + " v1.5.0"], m[c + " v1.5.2"]
        out[c + " v1.5.0"], out[c + " v1.5.2"] = a, b
        num = pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b)
        if num:                                   # compare at the v1.5.0 table's rounding
            nd = _decimals(c)
            diff = ((a - b.round(nd)).abs() > 0.5 * 10 ** (-nd) + 1e-9) | (a.isna() != b.isna())
        else:
            diff = a.astype(str) != b.astype(str)
        changed.append(diff.rename(c))
    ch = pd.concat(changed, axis=1)
    out["Fields changed"] = ch.apply(lambda s: "; ".join(k for k, v in s.items() if v), axis=1)
    out.to_csv(delta_csv, index=False)
    print(f"comparison -> {delta_csv}")
    for _, x in out.iterrows():
        print(f"  {x.ROI} / {x.Fleet}: {x['Fields changed'] or 'no change'}")
