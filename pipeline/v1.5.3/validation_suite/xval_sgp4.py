"""Independent cross-check of Cells 2, 4 and 5 on a REAL run (CKGR + Okavango,
open data + PlanetScope optimistic, 12 Oct 2026), open-data satellites only.

Independent of the pipeline: propagation with the `sgp4` package (Vallado's
SGP4, initialised from the same CelesTrak OMM records Cell 2 cached), TEME ->
Earth-fixed by GMST (IAU-82; UT1 ~ UTC, polar motion ignored: ~0.5 km at most),
WGS84 geodetic latitude by Bowring's iteration, tiles sampled every 0.25 s for
the closest point of approach (CPA), sun elevation by astropy.

Checks, per satellite pass over the ROI:
  A. completeness: every tile whose independent CPA off-nadir angle is below
     the swath-edge angle (minus a small tolerance) and is sunlit appears in
     arf_refined_events.csv for that satellite and pass;
  B. no false positives: no tile with CPA angle above the edge angle (plus
     tolerance) appears;
  C. the CPA angle and minute agree with the refined row with the SMALLEST
     angle for that (satellite, tile, pass): angle within 0.1 deg (the
     frame model here ignores polar motion and UT1-UTC, ~0.5 km), strip
     minute = round(CPA minute) except where the independent CPA lies within
     0.5 s of a strip boundary (0.25-s sampling here, plus the frame-model
     difference: such a tile is abeam at the boundary, in either strip);
  D. how many tiles also appear in a NEIGHBOURING minute strip (the end-cap
     duplication found in Cell 5).
Data files are only read; nothing in them is executed.
Usage: xval_sgp4.py [refined_events.csv] [report.json]"""
import json, sys, math
import numpy as np, pandas as pd, h3
from sgp4.api import Satrec, jday
from sgp4 import omm

RD = "/home/claude/realdata/"
RE = 6378.137                        # Cell 5's sphere (WGS84 equatorial radius)
A_WGS, F_WGS = 6378.137, 1 / 298.257223563
SIM0 = np.datetime64("2026-10-12T08:27:09.082080")        # Cell 2 sim_epoch (cell2 summary)
WIN_MIN = 48 * 60
SWATH = {"SENTINEL": 290.0, "LANDSAT": 185.0}              # MISSION_TABLE (sourced)
MIN_SUN = {"SENTINEL": 5.0, "LANDSAT": 5.0}

gp = json.load(open(RD + "ckgr_openPO_arf_gp_resource.json"))["records"]
open_sats = {r["OBJECT_NAME"]: r for r in gp
             if r["OBJECT_NAME"] in ("SENTINEL-2A", "SENTINEL-2B", "SENTINEL-2C", "LANDSAT 8", "LANDSAT 9")}
tiles = pd.read_csv(RD + "ckgr_openPO_mosaic_tiles.csv", usecols=["tile"])["tile"].to_numpy()
ll = np.array([h3.cell_to_latlng(t) for t in tiles])
tlat, tlon = ll[:, 0], ll[:, 1]
EV_PATH = sys.argv[1] if len(sys.argv) > 1 else RD + "ckgr_openPO_refined_events.csv"
REPORT = sys.argv[2] if len(sys.argv) > 2 else "/home/claude/validate/xval_sgp4_report.json"
ev = pd.read_csv(EV_PATH)
ev = ev[ev.sat.isin(list(open_sats))]


def gmst_rad(jd_ut1):
    T = (jd_ut1 - 2451545.0) / 36525.0
    sec = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + 0.093104 * T ** 2 - 6.2e-6 * T ** 3
    return np.radians((sec % 86400.0) / 240.0)


def geodetic(x, y, z):
    lon = np.degrees(np.arctan2(y, x))
    p = np.hypot(x, y)
    e2 = F_WGS * (2 - F_WGS)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(6):
        N = A_WGS / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        lat = np.arctan2(z + e2 * N * np.sin(lat), p)
    return np.degrees(lat), lon


def propagate(rec, t_min):
    sat = Satrec()
    omm.initialize(sat, rec)
    t = SIM0 + (np.asarray(t_min) * 60e6).astype("timedelta64[us]")
    s = (t - np.datetime64("2000-01-01T12:00:00")).astype("timedelta64[us]").astype(np.float64) / 86400e6
    jd = np.full(s.shape, 2451545.0)
    fr = s                                          # jd + fr = JD(UTC)
    e, r, v = sat.sgp4_array(jd, fr)
    assert not np.any(e), "SGP4 error code"
    th = gmst_rad(jd + fr)
    x = np.cos(th) * r[:, 0] + np.sin(th) * r[:, 1]
    y = -np.sin(th) * r[:, 0] + np.cos(th) * r[:, 1]
    z = r[:, 2]
    lat, lon = geodetic(x, y, z)
    return lat, lon, np.sqrt(x * x + y * y + z * z), t


def off_nadir(d_km, rs):
    g = d_km / RE
    slant = np.sqrt(rs ** 2 + RE ** 2 - 2 * rs * RE * np.cos(g))
    return np.degrees(np.arcsin(np.clip(RE * np.sin(g) / slant, -1, 1)))


def gc_km(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
import astropy.units as u
from astropy.utils import iers
iers.conf.auto_download = False

report = {"satellites": {}}
all_ok = True
for name, rec in sorted(open_sats.items()):
    fam = "SENTINEL" if name.startswith("SENTINEL") else "LANDSAT"
    t1 = np.arange(0, WIN_MIN + 1e-9, 1.0 / 60.0)                 # 1-s scan
    la, lo, rs, _ = propagate(rec, t1)
    near = gc_km(la, lo, tlat.mean(), tlon.mean()) < 1500.0
    if not near.any():
        report["satellites"][name] = {"passes": 0}
        continue
    idx = np.flatnonzero(near)
    groups = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    sat_rows = []
    spans = []
    for p_id, g in enumerate(groups):
        spans.append((t1[g[0]], t1[g[-1]]))
        tf = np.arange(t1[g[0]], t1[g[-1]] + 1e-9, 0.25 / 60.0)    # 0.25-s fine samples
        fla, flo, frs, ft = propagate(rec, tf)
        best = np.full(len(tiles), np.inf); bt = np.zeros(len(tiles)); brs = np.zeros(len(tiles))
        for k in range(0, len(tf), 400):
            sl = slice(k, k + 400)
            d = gc_km(fla[sl, None], flo[sl, None], tlat[None, :], tlon[None, :])
            a = off_nadir(d, frs[sl, None])
            j = a.argmin(axis=0)
            m = a[j, np.arange(len(tiles))]
            upd = m < best
            best[upd] = m[upd]; bt[upd] = tf[sl][j[upd]]; brs[upd] = frs[sl][j[upd]]
        cand = np.isfinite(best) & (best < 25.0)
        if not cand.any():
            continue
        # swath-edge angle at the CPA radius (Cell 5 uses the strip's mean radius)
        theta = off_nadir(np.full(len(tiles), SWATH[fam] / 2.0), brs)
        # sun elevation at the CPA instant (astropy, independent of Cell 5's NOAA formula)
        ci = np.flatnonzero(cand)
        tt = Time((SIM0 + (bt[ci] * 60e6).astype("timedelta64[us]")).astype("datetime64[ms]").astype(str), scale="utc")
        loc = EarthLocation.from_geodetic(tlon[ci] * u.deg, tlat[ci] * u.deg, 0 * u.m)
        sun = get_sun(tt).transform_to(AltAz(obstime=tt, location=loc)).alt.deg
        sun_all = np.full(len(tiles), np.nan); sun_all[ci] = sun
        for i in ci:
            sat_rows.append((p_id, i, best[i], bt[i], theta[i], sun_all[i]))
    if not sat_rows:
        report["satellites"][name] = {"passes": len(groups), "tiles_in_view": 0}
        continue
    ind = pd.DataFrame(sat_rows, columns=["pass_id", "i", "angle", "t_cpa", "theta", "sun"])
    ind["tile"] = tiles[ind.i]
    margin = 0.02                                     # deg; 2-s sampling and frame-model differences
    must = ind[(ind.angle < ind.theta - margin) & (ind.sun >= MIN_SUN[fam] + 0.5)]
    must_not = ind[(ind.angle > ind.theta + margin) | (ind.sun < MIN_SUN[fam] - 0.5)]
    e = ev[ev.sat == name].copy()
    # assign each refined row to the independent pass whose time span contains its strip minute
    e["pass_id"] = -1
    for p_id, (a0, a1) in enumerate(spans):
        e.loc[(e.t_min >= a0 - 1) & (e.t_min <= a1 + 1), "pass_id"] = p_id
    unassigned = int((e.pass_id < 0).sum())
    e_min = e.loc[e.groupby(["pass_id", "tile"]).A_nadir_deg.idxmin()] if len(e) else e
    got = set(zip(e.pass_id, e.tile))
    missing = sorted(set(zip(must.pass_id, must.tile)) - got)
    false_pos = sorted(set(zip(must_not.pass_id, must_not.tile)) & got)
    j = e_min.merge(ind, on=["pass_id", "tile"], how="inner")
    da = (j.A_nadir_deg - j.angle).abs()
    dt = (j.t_min - np.round(j.t_cpa)).abs()
    # seconds from the independent CPA to the nearest strip boundary (m + 0.5 min)
    _fr = np.mod(j.t_cpa.to_numpy() - 0.5, 1.0)
    to_boundary_s = 60.0 * np.minimum(_fr, 1.0 - _fr)
    mism = (dt > 0).to_numpy()
    ds = (j.sun_elev_deg - j.sun).abs()
    dup = e.groupby(["pass_id", "tile"]).t_min.nunique()
    res = {"passes": len(groups), "passes_with_refined_rows": int(e.pass_id[e.pass_id >= 0].nunique()),
           "refined_rows_unassigned": unassigned,
           "independent_in_fov_sunlit": int(len(must)), "refined_tile_passes": int(len(got)),
           "missing_from_refined": len(missing), "false_positives": len(false_pos),
           "angle_abs_diff_deg_max": float(da.max()) if len(da) else None,
           "angle_abs_diff_deg_p99": float(da.quantile(0.99)) if len(da) else None,
           "cpa_minute_mismatch_rows": int((dt > 0).sum()), "rows_compared": int(len(j)),
           "mismatch_max_seconds_from_strip_boundary": float(to_boundary_s[mism].max()) if mism.any() else None,
           "sun_abs_diff_deg_max": float(ds.max()) if len(ds) else None,
           "tile_passes_in_two_strips": int((dup > 1).sum()),
           "theta_refined_minus_independent_max_deg": float((j.theta_max_deg - j.theta).abs().max()) if len(j) else None}
    ok = res["missing_from_refined"] == 0 and res["false_positives"] == 0 and (res["angle_abs_diff_deg_max"] or 0) < 0.1 \
        and (not mism.any() or to_boundary_s[mism].max() <= 0.5) and (res["sun_abs_diff_deg_max"] or 0) < 0.6 \
        and unassigned == 0
    res["ok"] = ok
    all_ok &= ok
    report["satellites"][name] = res
    print(name, json.dumps(res), flush=True)
report["all_ok"] = all_ok
report["events_file"] = EV_PATH
json.dump(report, open(REPORT, "w"), indent=1)
print("ALL OK" if all_ok else "SEE REPORT")
