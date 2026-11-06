"""Summaries of post hoc amendments A12 (floor decomposition) and A13 (PlanetScope
acquisition fraction) over the 8 dry-season primary windows and the 8 wet-season
A5 windows, from the result files written in Colab on 4 November 2026.

Means over windows carry percentile 95% intervals from the protocol's
window-block bootstrap (B = 10,000, seed 20261030, whole windows resampled).
Cross-checks (the script stops if any fails):
  * A12: every window's decomposition passed the pipeline check; W + G + L = floor.
  * A13: at f = 1, per-window H1b equals the confirmatory per-window values (dry)
    and the A5 mean (wet); per-window reachable share and share of onsets reaching
    95% of the region equal A7.
Usage: summarise_a12_a13.py A12_DIR A13_DIR OUT_JSON
"""
import sys, glob, json, math
import numpy as np
import pandas as pd

A12, A13, OUT = sys.argv[1:4]
NT = 12231
B, SEED = 10_000, 20261030
FACTS = "/home/claude/ms/facts.json"
A7 = "/home/claude/arf/a7results/a7_windows.csv"
A5 = "/home/claude/ms/data/a5_summary.json"


def boot(values):
    v = np.asarray(values, float)
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(v), size=(B, len(v)))
    m = v[idx].mean(axis=1)
    return {"mean": float(v.mean()), "ci95": [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))],
            "min": float(v.min()), "max": float(v.max()), "windows": int(len(v))}


out = {"label": "POST HOC (amendments A12 and A13): descriptive, no decision", "bootstrap": {"B": B, "seed": SEED}}

# ------------------------------------------------------------------ A12
rows = []
for p in sorted(glob.glob(f"{A12}/*_decomp.csv")):
    season, win = p.split("/")[-1].split("_")[:2]
    s = json.load(open(p.replace(".csv", "_summary.json")))
    assert s["pipeline_check"]["checked"] and s["pipeline_check"]["floors_equal"] and s["pipeline_check"]["reachable_equal"], p
    d = pd.read_csv(p).assign(season=season, window=win)
    assert d["floor_min"].notna().all()
    assert np.allclose(d.W_wait_min + d.G_geometry_min + d.L_latency_min, d.floor_min)
    rows.append(d)
D = pd.concat(rows)
assert len(D) == 16 * 73
D["first_lst_h"] = (D.local_mean_solar_time_h + D.W_wait_min / 60) % 24
D["complete_lst_h"] = (D.local_mean_solar_time_h + (D.W_wait_min + D.G_geometry_min) / 60) % 24
a12 = {}
for season, g in D.groupby("season"):
    w = g.groupby("window").agg(floor=("floor_min", "median"), W=("W_wait_min", "median"), G=("G_geometry_min", "median"))
    a12[season] = {
        "floor_h_median": boot(w.floor / 60), "W_h_median": boot(w.W / 60), "G_h_median": boot(w.G / 60),
        "L_max_min": float(g.L_latency_min.max()),
        "onsets": int(len(g)),
        "pooled_median_W_share_of_floor": float((g.W_wait_min / g.floor_min).median()),
        "share_onsets_G_le_1h": float((g.G_geometry_min <= 60).mean()),
        "first_pass_lst_h_range": [float(g.first_lst_h.min()), float(g.first_lst_h.max())],
        "completion_lst_h_range": [float(g.complete_lst_h.min()), float(g.complete_lst_h.max())],
        "completion_lst_h_median": float(g.complete_lst_h.median()),
        "pooled_floor_h_median_onset_sun_below_horizon": float(g[g.sun_elev_centroid_deg < 0].floor_min.median() / 60),
        "pooled_floor_h_median_onset_sun_up": float(g[g.sun_elev_centroid_deg >= 0].floor_min.median() / 60),
        "completing_pass_planetscope_share": float(g.completing_sat.str.match(r"^FLOCK").mean())}
out["A12"] = a12

# ------------------------------------------------------------------ A13
need = math.ceil(0.95 * NT - 1e-9)
rows, summ = [], []
for p in sorted(glob.glob(f"{A13}/*_thin_onsets.csv")):
    season, win = p.split("/")[-1].split("_")[:2]
    rows.append(pd.read_csv(p).rename(columns={"fraction": "f"}).assign(season=season, window=win))
    s = json.load(open(p.replace("_onsets.csv", "_summary.json")))
    assert s["identity_check"] == "passed", p
    for fr, v in s["fractions"].items():
        summ.append({"season": season, "window": win, "f": float(fr), **v})
O, S = pd.concat(rows), pd.DataFrame(summ)
O["region95"] = O.reachable >= need
O["reach_pct"] = O.reachable / NT * 100
per = (O.groupby(["season", "window", "f", "seed"])
         .agg(reach=("reach_pct", "median"), region95=("region95", "mean"),
              ecc=("ecc_images", "mean"), fcfs=("fcfs_images", "mean"), floor=("floor_min", "median"))
         .groupby(["season", "window", "f"]).mean().reset_index())
per = per.merge(S[["season", "window", "f", "h1a_pct", "h1b_pct"]], on=["season", "window", "f"])
per["region95"] *= 100
per["saving"] = per.fcfs - per.ecc
per["floor"] = per["floor"] / 60                                   # hours
# cross-checks at f = 1
facts = json.load(open(FACTS))
dry1 = per[(per.season == "dry") & (per.f == 1.0)].sort_values("window")
assert np.allclose(np.round(dry1.h1b_pct.values, 2), facts["h1b_by_window"]), "A13 f=1 H1b differs from the confirmatory windows"   # facts.json holds 2 decimals
a5 = json.load(open(A5))["by_variant"]["p035"]["h1b_pct"]["mean"]
assert abs(per[(per.season == "wet") & (per.f == 1.0)].h1b_pct.mean() - a5) < 1e-9, "A13 f=1 wet H1b differs from A5"
a7 = pd.read_csv(A7)
a7 = a7[a7.arm.isin(["Primary: CKGR+Okavango, all 96, conservative", "Primary design, wet season (A5)"])]
chk = per[per.f == 1.0].merge(a7[["window", "region_t95_reached_pct", "reach_pct_median"]], on="window")
assert len(chk) == 16 and np.allclose(chk.region95, chk.region_t95_reached_pct) and np.allclose(chk.reach, chk.reach_pct_median)
out["A13_crosschecks"] = "at f = 1: per-window H1b equals the confirmatory values (dry) and the A5 mean (wet); reachable share and region-95 share equal A7 in all 16 windows"
a13 = {}
for (season, f), g in per.groupby(["season", "f"]):
    a13.setdefault(season, {})[f"{f:g}"] = {k: boot(g[k]) for k in
                                            ("reach", "region95", "h1a_pct", "h1b_pct", "ecc", "fcfs", "saving", "floor")}
out["A13"] = a13
per.to_csv(OUT.replace(".json", "_a13_windows.csv"), index=False)
D.to_csv(OUT.replace(".json", "_a12_onsets.csv"), index=False)
json.dump(out, open(OUT, "w"), indent=1)


def show(x):
    return f"{x['mean']:.2f} [{x['ci95'][0]:.2f}, {x['ci95'][1]:.2f}]"


for s in ("dry", "wet"):
    a = a12[s]
    print(f"A12 {s}: floor {show(a['floor_h_median'])} h = W {show(a['W_h_median'])} + G {show(a['G_h_median'])}; "
          f"completion LST {a['completion_lst_h_range'][0]:.2f}-{a['completion_lst_h_range'][1]:.2f} h")
    for f in ("1", "0.75", "0.5", "0.25"):
        b = a13[s][f]
        print(f"  A13 {s} f={f}: reach {show(b['reach'])}%, region95 {show(b['region95'])}%, H1a {b['h1a_pct']['mean']:.1f}, "
              f"H1b {show(b['h1b_pct'])}, ECC {show(b['ecc'])}, FCFS {show(b['fcfs'])}, saving {show(b['saving'])}")
