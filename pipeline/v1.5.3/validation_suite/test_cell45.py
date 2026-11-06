"""Cells 4 and 5 (access detection, refinement), v1.5.1 end-cap rule.
Usage: test_cell45.py PIPELINE_v151 PIPELINE_v150
A. _min_track_distance_km on synthetic tracks: abeam points get the cross-track
   distance; points whose closest approach lies beyond either end get +inf
   (v1.5.1) or the end-sample distance (legacy); merged clipped samples.
B. Full Cells 4-5 on a synthetic fleet over the real Lower Zambezi tiles,
   against a brute-force closest approach (0.5-s sampling of the same orbits,
   then a continuous ternary-search refinement on the analytic orbit):
   every (satellite, tile, pass) appears in exactly ONE strip, the strip that
   contains its closest approach, with the same off-nadir angle (v1.5.1);
   v1.5.0 shows the duplicates; no tile is lost."""
import sys, os, json, math, ast, re
import numpy as np, pandas as pd, h3
sys.path.insert(0, "/home/claude/validate")
import e2e15

P151, P150 = sys.argv[1], sys.argv[2]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


def defs(path, cell):
    C = e2e15.cells_of(open(path, encoding="utf-8").read())
    tree = ast.parse(e2e15.strip_org(C[cell]))
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    import typing
    ns = {"np": np, "pd": pd, "Re_km": 6378.137, **{n: getattr(typing, n) for n in ("List", "Dict", "Tuple", "Any", "Optional", "Set")},
          "refine_config": type("RC", (), {"downlink_min_elevation_deg": 15.0, "downlink_coarse_step_min": 5,
                                           "t_max_wait_minutes": 95.0, "swath_half_angle_deg": 9.22,
                                           "wgs84_re_km": 6378.137})()}
    exec(compile(tree, "defs", "exec"), ns)
    return ns


ns = defs(P151, "5")
f = ns["_min_track_distance_km"]
R = 6378.137
# track along the equator, 0.0 .. 3.6 deg longitude in 31 samples (2-s spacing of a ~6.7 km/s ground track)
s_lon = np.linspace(0.0, 3.6, 31); s_lat = np.zeros(31); s_rs = np.full(31, R + 786.0)
pts = [(0.5, 1.8, "abeam"), (-1.0, 0.0, "abeam at start"), (0.8, 3.6, "abeam at end"),
       (0.3, -0.4, "beyond start"), (-0.2, 4.0, "beyond end"), (1.2, 2.2, "abeam")]
c_lat = np.array([p[0] for p in pts]); c_lon = np.array([p[1] for p in pts])
d_new, _ = f(s_lat, s_lon, s_rs, c_lat, c_lon, R, exclude_caps=True)
d_old, _ = f(s_lat, s_lon, s_rs, c_lat, c_lon, R, exclude_caps=False)
ok = True
for (la, lo, kind), dn, do in zip(pts, d_new, d_old):
    xt = abs(la) * math.pi / 180 * R
    if kind.startswith("abeam"):
        ok &= abs(dn - xt) < 1e-6 and abs(do - xt) < 1e-6
    else:
        end_lo = 0.0 if "start" in kind else 3.6
        p1, p2 = 0.0, math.radians(la)
        a = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lo - end_lo) / 2) ** 2
        de = 2 * R * math.asin(math.sqrt(a))
        ok &= math.isinf(dn) and abs(do - de) < 1e-6
check(ok, "track distance: abeam -> cross-track distance (both rules); beyond an end -> +inf (v1.5.1) or end-sample distance (legacy)")
# clipped samples at a window edge (repeated first sample) behave as a shorter strip
s_lon2 = np.r_[np.zeros(15), np.linspace(0.0, 1.8, 16)]
d2, _ = f(np.zeros(31), s_lon2, np.full(31, R + 786.0), np.array([0.5, 0.5]), np.array([-0.3, 1.0]), R)
check(math.isinf(d2[0]) and abs(d2[1] - 0.5 * math.pi / 180 * R) < 1e-6,
      "clipped strip (repeated samples at t = 0): a point before the window start is excluded, an abeam point kept")

# ------------------------------------------------------------------ B: full Cells 4-5 vs brute force
tiles = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
sats = e2e15.eo_fleet(seed=11, n_dove=30)
runs = {}
for tag, path in (("v151", P151), ("v150", P150)):
    nsr = e2e15.run(path, tiles, cells=("4", "5"), sats=sats, overrides={"window_hours": 24, "onset_shift_hours": 12})
    runs[tag] = nsr["refined_events_df"].copy()
    runs[tag]["sat"] = runs[tag]["sat"].astype(str); runs[tag]["tile"] = runs[tag]["tile"].astype(str)
    if tag == "v151":
        cell5 = json.load(open(os.path.join(nsr["_out_dir"], "arf_cell5_summary.json")))
        stray = os.listdir(nsr["_work_dir"])
ll = np.array([h3.cell_to_latlng(t) for t in tiles])
tlat, tlon = ll[:, 0], ll[:, 1]


def dup_pairs(df):
    d = df.sort_values(["sat", "tile", "t_min"])
    same = (d.sat.values[1:] == d.sat.values[:-1]) & (d.tile.values[1:] == d.tile.values[:-1]) & \
           (np.diff(d.t_min.values) == 1)
    return int(same.sum())


check(dup_pairs(runs["v151"]) == 0 and dup_pairs(runs["v150"]) > 0,
      f"consecutive-strip duplicates: v1.5.0 {dup_pairs(runs['v150'])}, v1.5.1 {dup_pairs(runs['v151'])} "
      f"(Cell 5 summary: {cell5.get('rows_end_cap_excluded')} end-cap rows excluded)")

# brute force: per satellite, per pass (contiguous visibility), per tile: min off-nadir angle and its time
def off_nadir(dkm, rs):
    g = dkm / R
    sl = np.sqrt(rs ** 2 + R ** 2 - 2 * rs * R * np.cos(g))
    return np.degrees(np.arcsin(np.clip(R * np.sin(g) / sl, -1, 1)))


bf_rows = []
win = 24 * 60
for s in sats:
    t = np.arange(0, win + 1e-9, 0.5 / 60.0)
    la, lo = e2e15.orbit(s, t)
    near = np.hypot(la - tlat.mean(), (lo - tlon.mean()) * np.cos(np.radians(tlat.mean()))) < 12
    if not near.any():
        continue
    idx = np.flatnonzero(near)
    for g in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
        best = np.full(len(tiles), np.inf); bt = np.zeros(len(tiles))
        for k in range(0, len(g), 300):
            sl = g[k:k + 300]
            p1, p2 = np.radians(la[sl])[:, None], np.radians(tlat)[None, :]
            a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(tlon[None, :] - lo[sl][:, None]) / 2) ** 2
            ang = off_nadir(2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1))), R + s["h"])
            j = ang.argmin(axis=0); m = ang[j, np.arange(len(tiles))]
            upd = m < best; best[upd] = m[upd]; bt[upd] = t[sl][j[upd]]
        theta = off_nadir(np.array(s["swath_km"] / 2.0), R + s["h"])
        ci = np.flatnonzero(best < theta + 0.5)
        if ci.size == 0:
            continue
        # Continuous CPA: 0.5-s point samples are up to 0.25 s (~1.8 km) from the
        # true closest approach, an off-nadir error of ~0.2 deg near nadir for a
        # 475 km Dove. Refine by ternary search on the analytic orbit in
        # [bt - 0.5 s, bt + 0.5 s] (the distance is unimodal there).
        lo_t, hi_t = bt[ci] - 0.5 / 60.0, bt[ci] + 0.5 / 60.0
        def dist_at(tt):
            la_, lo_ = e2e15.orbit(s, tt)
            p1, p2 = np.radians(la_), np.radians(tlat[ci])
            a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(tlon[ci] - lo_) / 2) ** 2
            return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
        for _ in range(60):
            m1, m2 = lo_t + (hi_t - lo_t) / 3.0, hi_t - (hi_t - lo_t) / 3.0
            closer1 = dist_at(m1) < dist_at(m2)
            hi_t = np.where(closer1, m2, hi_t)
            lo_t = np.where(closer1, lo_t, m1)
        t_ref = 0.5 * (lo_t + hi_t)
        a_ref = off_nadir(dist_at(t_ref), R + s["h"])
        for i, a_i, t_i in zip(ci, a_ref, t_ref):
            bf_rows.append((s["id"], tiles[i], float(a_i), float(t_i), float(theta)))
bf = pd.DataFrame(bf_rows, columns=["sat", "tile", "angle", "t_cpa", "theta"])
ev = runs["v151"]
m = ev.merge(bf, on=["sat", "tile"], how="left")
m = m[(m.t_min - m.t_cpa).abs() <= 1.5]
unmatched = len(ev) - len(m.drop_duplicates(["sat", "tile", "t_min"]))
minute_ok = (m.t_min == np.floor(m.t_cpa + 0.5)) | ((m.t_cpa + 0.5) % 1.0 < 0.02) | ((m.t_cpa + 0.5) % 1.0 > 0.98)
angle_err = (m.A_nadir_deg - m.angle).abs()
# completeness: brute-force in-swath tiles (margin) that are sunlit must all appear (sun computed by the pipeline itself)
from_ev = set(zip(ev.sat, ev.tile))
must = bf[bf.angle < bf.theta - 0.02]
t_utc = e2e15.T0 + (must.t_cpa.to_numpy() * 60e6).astype("timedelta64[us]")
sun_ns = defs(P151, "5")["solar_elevation_deg"]
ll_m = np.array([h3.cell_to_latlng(t) for t in must.tile])
sun = np.array([float(np.ravel(sun_ns(np.datetime64(tt, "s"), la_, lo_))[0]) for tt, (la_, lo_) in zip(t_utc, ll_m)])
lim = must.sat.map({s["id"]: s["min_sun_elev_deg"] for s in sats}).to_numpy()
must = must[sun >= lim + 0.5]
missing = [k for k in zip(must.sat, must.tile) if k not in from_ev]
check(unmatched == 0 and minute_ok.all() and angle_err.max() < 0.02,
      f"v1.5.1 rows: each in the strip containing its closest approach ({int((~minute_ok).sum())} off), off-nadir "
      f"angle vs brute force max error {angle_err.max():.4f} deg, unmatched rows {unmatched}")
check(not missing, f"completeness: all {len(must)} brute-force in-swath sunlit (satellite, tile) passes present "
                   f"({len(missing)} missing)")
tiles150 = set(zip(runs["v150"].sat, runs["v150"].tile)); tiles151 = set(zip(ev.sat, ev.tile))
check(tiles151 <= tiles150 and len(tiles150 - tiles151) <= 0.002 * len(tiles150),
      f"no (satellite, tile) lost by the end-cap rule except at the window edges: {len(tiles150 - tiles151)} of {len(tiles150)}")
check(not stray, f"no files written outside OUT_DIR ({stray})")
_ver = re.search(r'^__version__ = "([^"]+)"', open(P151, encoding="utf-8").read(), re.M).group(1)
if tuple(int(x) for x in _ver.split(".")) >= (1, 5, 2):
    n_old, n_new = len(runs["v150"]), len(runs["v151"])
    check(cell5.get("rows_end_cap_excluded") == n_old - n_new,
          f"v1.5.2 summary: end-cap rows not counted ({cell5.get('rows_end_cap_excluded')}) = v1.5.0 rows minus "
          f"this release's rows ({n_old} - {n_new})")
n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell45.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
