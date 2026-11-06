"""Checks of release 1.5.5 (proxy_swath_km) on synthetic Starlink-like orbits, real Orekit.
  1. proxy_swath_km unset: Cells 0-5 of 1.5.5 give the same refined events as 1.5.4, row for row.
  2. proxy_swath_km = 32.5: every proxy satellite carries the swath; each pass's edge angle equals
     the off-nadir angle of a point 16.25 km from the ground track at that satellite's altitude
     (within the altitude spread of the strip); every kept row lies inside it and in daylight;
     search radius and folder tag as specified.
  3. Cell 0 refuses the setting for a real fleet and outside [0.5, 400] km.
Usage: test_v155.py S_DIR PIPE_155
"""
import sys, os, re, json, math
import numpy as np
import pandas as pd

S, PIPE = sys.argv[1], sys.argv[2]
res = {}
Re = 6378.137


def inc_angle(ground_km, rs_km):            # the pipeline's _incidence_angle_deg_vec, scalar form
    g = ground_km / Re
    slant = math.sqrt(rs_km ** 2 + Re ** 2 - 2 * rs_km * Re * math.cos(g))
    return math.degrees(math.asin(Re * math.sin(g) / slant))


# 1 ---------------------------------------------------------------------------------------------
a = pd.read_csv(os.path.join(S, "o154", "refined_events_cell5.csv"))
b = pd.read_csv(os.path.join(S, "o155none", "refined_events_cell5.csv"))
same = a.shape == b.shape and list(a.columns) == list(b.columns) and bool((a.astype(str).values == b.astype(str).values).all())
res["1_unset_equals_154"] = {"rows": int(len(a)), "identical": same}
ta = json.load(open(os.path.join(S, "o154", "timing.json")))
tb = json.load(open(os.path.join(S, "o155none", "timing.json")))
res["1_unset_equals_154"]["fleet_tag_same"] = ta["fleet_tag"] == tb["fleet_tag"]

# 2 ---------------------------------------------------------------------------------------------
c = pd.read_csv(os.path.join(S, "o155sw", "refined_events_cell5.csv"))
tc = json.load(open(os.path.join(S, "o155sw", "timing.json")))
con = json.load(open(os.path.join(S, "o155sw", "constellation.json")))
sen = json.load(open(os.path.join(S, "o155sw", "sensor.json")))
alt = {s["id"]: s["alt_km"] for s in con}
th = c.groupby(["sat", "t_min"])["theta_max_deg"].first().reset_index()
exp_lo = th["sat"].map(lambda s: inc_angle(16.25, Re + alt[s] + 15.0))   # higher -> smaller angle
exp_hi = th["sat"].map(lambda s: inc_angle(16.25, Re + alt[s] - 15.0))
inside = ((th["theta_max_deg"] >= exp_lo - 1e-9) & (th["theta_max_deg"] <= exp_hi + 1e-9))
res["2_swath_32.5"] = {
    "fleet_tag": tc["fleet_tag"], "tag_ok": tc["fleet_tag"].endswith("_sw32.5km"),
    "all_satellites_carry_swath": all(s["swath_km"] == 32.5 for s in con),
    "passes": int(len(th)), "edge_angle_range_deg": [round(float(th.theta_max_deg.min()), 4), round(float(th.theta_max_deg.max()), 4)],
    "edge_angle_matches_altitude": bool(inside.all()),
    "rows": int(len(c)), "rows_inside_edge": bool((c["A_nadir_deg"] <= c["theta_max_deg"] + 1e-12).all()),
    "rows_in_daylight_sun_ge_5": bool((c["sun_elev_deg"] >= 5.0).all()),
    "footprint_radius_km": round(sen["footprint_radius_km"], 3),
    "rows_vs_unset": round(len(c) / len(b), 4),
    "cells_per_pass_unset_vs_32.5": [round(float(b.groupby(["sat", "t_min"]).size().mean()), 1),
                                     round(float(c.groupby(["sat", "t_min"]).size().mean()), 1)],
    "max_off_nadir_unset_deg": round(float(b["A_nadir_deg"].max()), 3)}
# expected search radius: 1.05 * swath/2 + half the along-track sample spacing at the fastest ground speed
mu = 398600.4418
h = max(150.0, sen["min_altitude_km"] - 15.0)          # the pipeline's lowest altitude over the window (Cell 4)
v = math.sqrt(mu / (Re + h)) * 60.0 * Re / (Re + h) + 2 * math.pi * Re / 1436.068
_exp_r = 1.05 * 32.5 / 2 + 0.5 * v * 2.0 / 60.0
res["2_swath_32.5"]["footprint_radius_expected_km"] = round(_exp_r, 3)
res["2_swath_32.5"]["footprint_radius_ok"] = abs(_exp_r - sen["footprint_radius_km"]) < 1e-6   # unrounded

# 3 ---------------------------------------------------------------------------------------------
src = open(PIPE, encoding="utf-8").read()
hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
ms = list(hdr.finditer(src))
cell0 = src[ms[0].start(): ms[1].start()]
refusals = {}
for label, ov in (("real_fleet", {"fleet": "earth-resources", "proxy_swath_km": 30.0}),
                  ("too_small", {"fleet": "starlink", "proxy_swath_km": 0.1}),
                  ("too_large", {"fleet": "starlink", "proxy_swath_km": 500.0}),
                  ("accepted", {"fleet": "starlink", "active_n": 100, "proxy_swath_km": 8.0})):
    ns = {"__name__": "__main__", "display": lambda *a, **k: None, "ARF_SKIP_INSTALL": True,
          "ARF_ALLOW_MIXED_VERSIONS": True, "ARF_CONFIG_OVERRIDES": ov, "ARF_BOUNDARY_DIR": "/nonexistent",
          "ARF_SHARED_CACHE_DIR": os.path.join(S, "cell0_cache")}
    old = os.getcwd(); os.makedirs(os.path.join(S, "cell0_cwd"), exist_ok=True); os.chdir(os.path.join(S, "cell0_cwd"))
    try:
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(cell0, "cell0", "exec"), ns)
        refusals[label] = {"raised": False, "tag": ns.get("FLEET_TAG"), "PROXY_SWATH_KM": ns.get("PROXY_SWATH_KM")}
    except Exception as e:
        refusals[label] = {"raised": True, "type": type(e).__name__, "msg": str(e)[:120]}
    finally:
        os.chdir(old)
res["3_cell0"] = refusals

checks = [res["1_unset_equals_154"]["identical"], res["1_unset_equals_154"]["fleet_tag_same"],
          res["2_swath_32.5"]["tag_ok"], res["2_swath_32.5"]["all_satellites_carry_swath"],
          res["2_swath_32.5"]["edge_angle_matches_altitude"], res["2_swath_32.5"]["rows_inside_edge"],
          res["2_swath_32.5"]["rows_in_daylight_sun_ge_5"], res["2_swath_32.5"]["footprint_radius_ok"],
          refusals["real_fleet"]["raised"], refusals["too_small"]["raised"], refusals["too_large"]["raised"],
          (not refusals["accepted"]["raised"]) and refusals["accepted"]["tag"] == "N100_sw8km"]
res["passed"] = f"{sum(map(bool, checks))} of {len(checks)}"
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(S, "test_v155_result.json"), "w"), indent=1)
sys.exit(0 if all(checks) else 1)
