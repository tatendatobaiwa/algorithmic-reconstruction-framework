"""facts_a12a13.json: every value the text may quote from amendments A12 and A13,
formatted as printed, with its source key in data/a12_a13_summary.json."""
import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..")) if _o.path.basename(_o.path.dirname(_o.path.abspath(__file__))) in ("figscripts", "checks") else None; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
from paths import D as _D, A as _A, FACTS as _FACTS
import json
s = json.load(open(_D + "a12_a13_summary.json"))
R = {}


def put(v, src, fmts=("{:.0f}", "{:.1f}", "{:.2f}")):
    for f in fmts:
        R.setdefault(f.format(v), []).append(src)


def hhmm(h, src):
    hh, mm = int(h), int(round((h - int(h)) * 60))
    if mm == 60:
        hh, mm = hh + 1, 0
    R.setdefault(f"{hh:02d}", []).append(src + " (hh)")
    R.setdefault(f"{mm:02d}", []).append(src + " (mm)")
    R.setdefault(f"{hh:02d}:{mm:02d}", []).append(src)


for season, a in s["A12"].items():
    for k in ("floor_h_median", "W_h_median", "G_h_median"):
        put(a[k]["mean"], f"A12 {season} {k} mean")
        for i, v in enumerate(a[k]["ci95"]):
            put(v, f"A12 {season} {k} ci{i}")
    put(100 * a["pooled_median_W_share_of_floor"], f"A12 {season} W share of floor, %")
    put(100 * a["completing_pass_planetscope_share"], f"A12 {season} PlanetScope completes, %")
    put(a["pooled_floor_h_median_onset_sun_below_horizon"], f"A12 {season} floor, onset in darkness")
    put(a["pooled_floor_h_median_onset_sun_up"], f"A12 {season} floor, onset after sunrise")
    for i, v in enumerate(a["completion_lst_h_range"]):
        hhmm(v, f"A12 {season} completion LST {'min' if i == 0 else 'max'}")
    for i, v in enumerate(a["first_pass_lst_h_range"]):
        hhmm(v, f"A12 {season} first pass LST {'min' if i == 0 else 'max'}")
for season, fr in s["A13"].items():
    for f, x in fr.items():
        put(float(f), f"A13 fraction {f}", ("{:g}",))
        for k, b in x.items():
            put(b["mean"], f"A13 {season} f={f} {k} mean")
            for i, v in enumerate(b["ci95"]):
                put(v, f"A13 {season} f={f} {k} ci{i}")
put(sum(v["onsets"] for v in s["A12"].values()), "A12 onsets, both seasons", ("{:,}",))
json.dump(R, open("facts_a12a13.json", "w"), indent=0)
print(len(R), "values")
