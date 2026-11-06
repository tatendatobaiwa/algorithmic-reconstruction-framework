# =============================================================================
# CELL 5: PHYSICAL CONSTRAINTS & REFINEMENT (Incidence Angle & Downlink)
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Implements physical Earth Observation (EO) constraints for the Spatio-Temporal
#   Selection Function (Equation 4 of the research proposal):
#     1. A_nadir(i, j): Precise off-nadir incidence angle calculation & FOV filter.
#     2. D_downlink(i): Ground station line-of-sight visibility & wait-time reward.
#
# Kinematic Strip Adaptation:
#   Each Cell 4 event is a one-minute push-broom strip. A tile is imaged when
#   the ground track passes abeam of it (closest point of approach, CPA), so
#   its off-nadir angle comes from its great-circle distance to the track:
#   real fleets use the distance to the arcs between the strip's 2-s samples
#   (exact between samples), and each tile observation belongs to the ONE
#   strip that contains its CPA (v1.5.1; checked against an independent SGP4
#   propagation, xval_sgp4.py); the Starlink proxy keeps the nearest sample.
#   The samples are reused from Cell 4 rather than re-propagated. Sunlit tiles
#   only (each mission's limit), and passes inside a sourced acquisition
#   exclusion (Cell 0, v1.5.2) are dropped.
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 5 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["5"] = _ARF_CELL_VERSION

import os
import sys
import time
import logging
from typing import Dict, List, Tuple, Set, Optional, Any
from dataclasses import dataclass

import numpy as np
import pandas as pd
from org.orekit.frames import TopocentricFrame
from org.orekit.bodies import GeodeticPoint
from org.hipparchus.util import FastMath
from org.orekit.utils import Constants
import h3

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logger = logging.getLogger("arf.cell5")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# =============================================================================
# CUSTOM EXCEPTIONS
# =============================================================================

class RefinementError(Exception):
    """Raised when geometric or downlink refinement fails."""
    pass

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()

# =============================================================================
# REFINEMENT CONFIGURATION
# =============================================================================

@dataclass
class RefinementConfig:
    """Physical refinement parameters."""
    downlink_min_elevation_deg: float = 15.0
    downlink_coarse_step_min: int = 5
    t_max_wait_minutes: float = 95.0
    swath_half_angle_deg: float = 9.22
    wgs84_re_km: float = Constants.WGS84_EARTH_EQUATORIAL_RADIUS / 1000.0

    def __post_init__(self):
        if 'DOWNLINK_MIN_ELEVATION_DEG' in globals():
            self.downlink_min_elevation_deg = DOWNLINK_MIN_ELEVATION_DEG
        if 'T_MAX_WAIT_MINUTES' in globals():
            self.t_max_wait_minutes = T_MAX_WAIT_MINUTES
        if 'SWATH_HALF_ANGLE_DEG' in globals():
            self.swath_half_angle_deg = SWATH_HALF_ANGLE_DEG

refine_config = RefinementConfig()
_real_fleet5 = globals().get("FLEET", "starlink") != "starlink"
logger.info(f"Refinement Parameters: Min Elevation = {refine_config.downlink_min_elevation_deg}°, "
            f"T_max = {refine_config.t_max_wait_minutes:.1f} min, "
            + ("swath edge: each satellite's own (MISSION_TABLE)" if _real_fleet5 else
               f"Swath Half-Angle = {refine_config.swath_half_angle_deg}° (Starlink proxy)"))

# =============================================================================
# SATELLITE ACCESS HELPERS
# =============================================================================

def get_satellite_id(sat: Any) -> str:
    if hasattr(sat, 'name'):
        return sat.name
    elif isinstance(sat, dict):
        return sat.get('id', sat.get('name', 'Unknown'))
    return str(sat)

def get_satellite_propagator(sat: Any) -> Any:
    if hasattr(sat, 'propagator'):
        return sat.propagator
    elif isinstance(sat, dict):
        return sat.get('propagator')
    raise RefinementError(f"No propagator found in {type(sat).__name__}")

# =============================================================================
# GROUND STATION TOPOCENTRIC FRAMES
# =============================================================================

gs_frames: Dict[str, TopocentricFrame] = {}
for name, gs in GROUND_STATIONS.items():
    gp = GeodeticPoint(
        FastMath.toRadians(gs["lat"]),
        FastMath.toRadians(gs["lng"]),
        gs["alt_m"]
    )
    gs_frames[name] = TopocentricFrame(earth, gp, name)

# =============================================================================
# DOWNLINK WINDOW DETECTION
# =============================================================================

def _visibility_radius_km(altitude_km: float, min_elev_deg: float) -> float:
    """
    Ground (great-circle) distance from a station within which a satellite at
    altitude_km appears at or above min_elev_deg:
        lambda = arccos(Re * cos(E) / (Re + h)) - E.
    At E = 0 this is the horizon distance (~2,559 km at 550 km altitude).
    """
    E = np.radians(min_elev_deg)
    return float(Re_km * (np.arccos(Re_km * np.cos(E) / (Re_km + altitude_km)) - E))


def detect_downlink_windows(
    constellation: List[Any],
    sim_epoch: Any,
    window_hours: int,
    min_elev_deg: float = refine_config.downlink_min_elevation_deg,
    coarse_step_min: int = refine_config.downlink_coarse_step_min,
    fine_step_sec: int = 10,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Records every satellite-to-ground-station contact minute at or above the
    elevation mask.

    Stage 1 (coarse, every coarse_step_min): flag a station when the
    sub-satellite point lies within the visibility radius of the HIGHEST
    satellite (+15 km altitude pad) plus the distance any satellite can travel
    in half a coarse step. The nearest coarse sample to any contact instant is
    then guaranteed to be flagged, so no pass can fall between samples. (A
    previous elevation-slack test -- flag at >= mask - 5 deg -- missed 2-3% of
    passes: low passes that clear the mask but never exceed mask - 5 deg at a
    5-minute sample.)

    Stage 2 (fine, every fine_step_sec within half a coarse step of a flag): exact
    Orekit topocentric elevation. A contact at time t is recorded against the
    next whole minute, ceil(t), so sub-minute passes are kept and the wait
    derived from this table is never under-estimated. Each satellite is
    propagated once per time step for all stations. Resolution limit: a contact
    above the mask for less than fine_step_sec can be missed -- validated
    residuals were 3-9 s grazing passes, too brief to be usable downlinks.
    """
    window_minutes = window_hours * 60
    win_sec = window_minutes * 60
    step_sec = coarse_step_min * 60
    max_alt = sensor_config.max_altitude_km + 15.0
    reach_km = (_visibility_radius_km(max_alt, min_elev_deg)
                + sensor_config.max_ground_speed_km_min * coarse_step_min / 2.0)
    gs_ll = {n: (float(GROUND_STATIONS[n]["lat"]), float(GROUND_STATIONS[n]["lng"])) for n in gs_frames}

    t0 = time.time()
    n_coarse_hits = n_fine_samples = 0
    best = {}                                    # (sat, gs, minute) -> max elevation
    coarse_secs = list(range(0, win_sec + 1, step_sec))
    if coarse_secs[-1] != win_sec:
        coarse_secs.append(win_sec)

    for sat in constellation:
        sat_id = get_satellite_id(sat)
        propagator = get_satellite_propagator(sat)

        # Stage 1: one propagation per coarse step, tested against all stations.
        windows = {n: [] for n in gs_frames}
        for ts in coarse_secs:
            lat, lon, _ = _subsatellite(propagator, ts / 60.0)
            for n, (gla, glo) in gs_ll.items():
                if _hav_km(lat, lon, gla, glo) <= reach_km:
                    # Every contact instant lies within step/2 of a flagged coarse
                    # sample, so +/-(step/2 + one fine step) is sufficient.
                    half = step_sec // 2 + fine_step_sec
                    windows[n].append((max(0, ts - half), min(win_sec, ts + half)))
                    n_coarse_hits += 1
        windows = {n: w for n, w in windows.items() if w}
        if not windows:
            continue

        # Stage 2: exact elevation, one propagation per fine sample.
        fine_secs = sorted({s for w in windows.values() for lo, hi in w
                            for s in range(lo, hi + 1, fine_step_sec)})
        n_fine_samples += len(fine_secs)
        for ts in fine_secs:
            t = sim_epoch.shiftedBy(float(ts))
            pos = propagator.getPVCoordinates(t, itrf).getPosition()
            for n, w in windows.items():
                if not any(lo <= ts <= hi for lo, hi in w):
                    continue
                elev = FastMath.toDegrees(gs_frames[n].getElevation(pos, itrf, t))
                if elev >= min_elev_deg:
                    key = (sat_id, n, -(-ts // 60))       # ceil to the next whole minute
                    if elev > best.get(key, -90.0):
                        best[key] = elev

    if best:
        df = pd.DataFrame([{"sat": s, "gs": g, "t_min": int(m), "elev": e}
                           for (s, g, m), e in sorted(best.items())])
    else:
        df = pd.DataFrame(columns=["sat", "gs", "t_min", "elev"])

    if verbose:
        logger.info(f"Downlink scan: {n_coarse_hits} coarse hits, {n_fine_samples} fine samples, "
                    f"{time.time() - t0:.1f}s | reach radius {reach_km:.0f} km")
        logger.info(f"Downlink Opportunities: {len(df)} contact minutes identified across "
                    f"{df['sat'].nunique() if len(df) > 0 else 0} satellites.")

    return df

logger.info("Scanning for ground station downlink contact windows...")
_window_hours = WINDOW_HOURS if 'WINDOW_HOURS' in globals() else 48

downlink_events_df = detect_downlink_windows(
    constellation=constellation,
    sim_epoch=sim_epoch,
    window_hours=_window_hours
)

_downlink_lookup = (
    downlink_events_df.groupby("sat")["t_min"]
    .apply(lambda x: sorted(x.unique()))
    .to_dict()
    if len(downlink_events_df) > 0 else {}
)

def time_to_next_downlink(sat_id: str, event_t_min: int) -> Optional[int]:
    """Minutes from event_t_min to the satellite's first contact minute at or
    after it (None if none in the window). Binary search (v1.5.2; was a
    linear scan, same result) over the sorted contact minutes."""
    import bisect
    minutes = _downlink_lookup.get(sat_id, [])
    i = bisect.bisect_left(minutes, event_t_min)
    return int(minutes[i] - event_t_min) if i < len(minutes) else None

def downlink_reward(t_wait_min: Optional[int], t_max_min: float = refine_config.t_max_wait_minutes) -> float:
    if t_wait_min is None:
        return 0.0
    return max(0.0, 1.0 - t_wait_min / t_max_min)

# =============================================================================
# VECTORISED SPHERICAL TRIGONOMETRY (A_nadir)
# =============================================================================

def _incidence_angle_deg_vec(ground_distance_km: np.ndarray, Re_km: float, Rs_km: float) -> np.ndarray:
    """
    Sensor off-nadir angle (deg) to a ground point at the given great-circle
    distance from the sub-satellite point, for a satellite at geocentric
    radius Rs_km over a sphere of radius Re_km. Points beyond the horizon are
    returned as +inf: without that guard, arcsin(Re*sin(gamma)/slant) turns
    small again on the far side of the Earth and would report a spurious
    near-nadir angle for an unobservable point.
    """
    gamma = np.asarray(ground_distance_km, dtype=float) / Re_km
    slant_range = np.sqrt(Rs_km**2 + Re_km**2 - 2.0 * Rs_km * Re_km * np.cos(gamma))
    sin_eta = np.clip(Re_km * np.sin(gamma) / slant_range, -1.0, 1.0)
    eta = np.degrees(np.arcsin(sin_eta))
    return np.where(gamma < np.arccos(Re_km / Rs_km), eta, np.inf)


def solar_elevation_deg(t_utc, lat_deg, lon_deg):
    """
    Solar elevation (degrees above the horizon, geometric, no refraction) at
    ground points (lat_deg, lon_deg, east-positive) and UTC datetime(s) t_utc,
    by NOAA GML "General Solar Position Calculations"
    (https://gml.noaa.gov/grad/solcalc/solareqns.PDF): fractional year,
    equation of time, declination, true solar time, hour angle, zenith.
    Vectorised over arrays of points; t_utc a datetime or an array of them.
    """
    t = np.atleast_1d(np.asarray(t_utc, dtype="datetime64[s]"))
    years = t.astype("datetime64[Y]")
    doy = (t.astype("datetime64[D]") - years).astype(int) + 1
    yr = years.astype(int) + 1970
    ndays = np.where((yr % 4 == 0) & ((yr % 100 != 0) | (yr % 400 == 0)), 366, 365)
    sec = (t - t.astype("datetime64[D]")).astype(int)
    hour = sec / 3600.0
    g = 2.0 * np.pi / ndays * (doy - 1 + (hour - 12.0) / 24.0)
    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                       - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = hour * 60.0 + eqtime + 4.0 * np.asarray(lon_deg)           # timezone 0 (UTC)
    ha = np.radians(tst / 4.0 - 180.0)
    lat = np.radians(lat_deg)
    cos_zen = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(ha)
    return 90.0 - np.degrees(np.arccos(np.clip(cos_zen, -1.0, 1.0)))


# Per-satellite sensor geometry and daylight limit (real fleets); the Starlink
# proxy has neither and keeps the uniform half-angle with no daylight limit.
_sat_meta = {s["id"]: s for s in constellation} if 'constellation' in globals() else {}


def _min_track_distance_km(s_lat, s_lon, s_rs, c_lat, c_lon, re_km, exclude_caps=True):
    """
    Minimum great-circle distance (km) from each ground point (c_lat, c_lon) to
    the sub-satellite GROUND TRACK through the strip samples -- to the arcs
    between consecutive samples, not just the samples. Point sampling misses a
    tile lying between samples (2-s samples are ~15 km apart), which is
    negligible for a ~80 km footprint but loses whole tiles for an 8 km swath.

    exclude_caps (v1.5.1): a push-broom or frame imager sees a ground point when
    the track passes ABEAM of it (its closest point of approach, CPA). A point
    whose CPA lies before the strip's first sample or after its last is imaged
    in the neighbouring strip, so its distance here is +inf: the foot of its
    perpendicular on the great circle of the first (last) arc lies before the
    first (after the last) sample. Up to v1.5.0 such points took their distance
    to the end sample, so every strip also claimed a half-disc ("end cap") of
    radius half the swath beyond each end: 73-80% of the Sentinel-2 tile
    observations of the 12 October CKGR runs appeared in two consecutive strips
    (checked against an independent SGP4 propagation). Samples clipped to the
    window edges are merged first, so a strip at t = 0 or at the window end ends
    there.
    Returns (distance_km, geocentric radius of the nearest arc) per point;
    distance +inf where the CPA lies outside the strip (exclude_caps).
    """
    def unit(lat, lon):
        la, lo = np.radians(lat), np.radians(lon)
        return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    tol = 1e-12                                                              # rad: rounding at arc ends
    P = unit(np.asarray(s_lat, float), np.asarray(s_lon, float))              # (k, 3) samples
    rs = np.asarray(s_rs, float)
    if exclude_caps and len(P) > 1:
        # merge repeated samples (strips clipped at the window edges)
        keep = np.r_[True, np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-15]
        P, rs = P[keep], rs[keep]
    C = unit(np.asarray(c_lat, float), np.asarray(c_lon, float))              # (m, 3) points
    best = np.arccos(np.clip(C @ P.T, -1.0, 1.0))                            # (m, k) to samples
    best_ang, best_idx = best.min(axis=1), best.argmin(axis=1)
    rs_best = rs[best_idx]
    for i in range(len(P) - 1):
        A, B = P[i], P[i + 1]
        n = np.cross(A, B)
        nn = np.linalg.norm(n)
        if nn < 1e-15:
            continue
        n = n / nn
        sin_xt = C @ n                                                        # cross-track
        F = C - np.outer(sin_xt, n)                                           # foot on the great circle
        t_arc = tol if exclude_caps else 0.0                                  # legacy: exact v1.5.0 test
        on_arc = ((np.cross(A, F) @ n) >= -t_arc) & ((np.cross(F, B) @ n) >= -t_arc)
        ang = np.abs(np.arcsin(np.clip(sin_xt, -1.0, 1.0)))
        better = on_arc & (ang < best_ang)
        best_ang = np.where(better, ang, best_ang)
        rs_best = np.where(better, 0.5 * (rs[i] + rs[i + 1]), rs_best)
    if exclude_caps and len(P) > 1:
        def _outside(A, B, before):
            n = np.cross(A, B)
            nn = np.linalg.norm(n)
            if nn < 1e-15:
                return np.zeros(len(C), dtype=bool)
            n = n / nn
            F = C - np.outer(C @ n, n)
            return ((np.cross(A, F) @ n) < -tol) if before else ((np.cross(F, B) @ n) < -tol)
        outside = _outside(P[0], P[1], True) | _outside(P[-2], P[-1], False)
        best_ang = np.where(outside, np.inf, best_ang)
    return best_ang * re_km, rs_best


def _sim_epoch_datetime64():
    """sim_epoch (Orekit AbsoluteDate) as numpy datetime64 UTC."""
    from org.orekit.time import TimeScalesFactory
    return np.datetime64(sim_epoch.toString(TimeScalesFactory.getUTC()).replace("Z", "")[:26])


def _acquisition_exclusion_windows():
    """v1.5.2: ACQUISITION_EXCLUSIONS (Cell 0) as {satellite id: [(start_min,
    end_min, entry index)]}, minutes from the analysis start. A satellite id
    may carry a ' [NORAD id]' suffix (Cell 2, duplicate names); the name before
    it is matched."""
    import re as _re5
    out = {}
    excl = globals().get("ACQUISITION_EXCLUSIONS") or []
    if not excl or not _sat_meta or globals().get("FLEET", "starlink") == "starlink":
        return out
    t0 = _sim_epoch_datetime64()
    for i, ex in enumerate(excl):
        a = float((np.datetime64(ex["start_utc"]) - t0) / np.timedelta64(1, "s")) / 60.0
        b = float((np.datetime64(ex["end_utc"]) - t0) / np.timedelta64(1, "s")) / 60.0
        for sid in _sat_meta:
            if _re5.match(ex["name_pattern"], str(sid).split(" [")[0].strip(), flags=_re5.IGNORECASE):
                out.setdefault(sid, []).append((a, b, i))
    return out


def refine_access_events(
    access_events_df: pd.DataFrame,
    constellation: List[Any],
    sim_epoch: Any,
    swath_half_angle_deg: float = refine_config.swath_half_angle_deg,
) -> pd.DataFrame:
    """
    Exact incidence-angle FOV filter and downlink reward for each strip event.
    Uses the strip samples Cell 4 recorded (strip_samples); constellation and
    sim_epoch are kept in the signature for API compatibility. Distances use
    the same sphere (Re_km = WGS84 equatorial) as the incidence geometry.
    v1.5.2: passes inside a sourced acquisition exclusion (Cell 0) are dropped
    after the geometric and daylight filters and counted; the end-cap count of
    the summary is taken after the daylight filter, so it equals the rows the
    v1.5.0 rule would have added to this table.
    """
    # Accumulated as arrays (one chunk per event), not one dict per row: at
    # N=1000 country scale this table reaches tens of millions of rows.
    ev_sat, ev_t, ev_n, ev_d, ev_w, ev_theta = [], [], [], [], [], []
    tile_chunks, angle_chunks, sun_chunks = [], [], []
    global refine_daylight_stats
    refine_daylight_stats = {"events_dark": 0, "rows_dark": 0, "rows_lit": 0, "rows_end_cap_excluded": 0,
                             "rows_acquisition_excluded": 0, "passes_acquisition_excluded": 0,
                             "acquisition_exclusions_hit": {}}
    exclude_caps = not globals().get("LEGACY_STRIP_END_CAPS", False)        # v1.5.1 (Cell 0)
    _t0_utc = _sim_epoch_datetime64() if any(m.get("min_sun_elev_deg") is not None for m in _sat_meta.values()) else None
    _excl_win = _acquisition_exclusion_windows()
    logger.info(f"Evaluating physical constraints across {len(access_events_df):,} continuous strip events...")
    for row in access_events_df.itertuples(index=False):
        sat_id, t_min = row.sat, int(row.t_min)
        s_lat, s_lon, s_rs = strip_samples[(sat_id, t_min)]
        meta = _sat_meta.get(sat_id, {})
        excluded_by = next((i for a, b, i in _excl_win.get(sat_id, ()) if a <= t_min < b), None)
        cap_ok = None
        if meta.get("swath_km"):
            # Real fleets: exact distance to the ground-track arcs. The search
            # radius (Cell 4) includes half the sample spacing, so a tile
            # between samples is still a candidate.
            cand = _tiles_near(s_lat, s_lon, sensor_config.footprint_radius_km)
            if cand.size == 0:
                continue
            c_lat, c_lon = _tile_lat[cand], _tile_lon[cand]
            d, rs_near = _min_track_distance_km(s_lat, s_lon, s_rs, c_lat, c_lon, refine_config.wgs84_re_km,
                                                exclude_caps=exclude_caps)
            a_nadir = np.full(cand.size, np.inf)            # +inf: closest approach outside this strip
            fin = np.isfinite(d)
            if fin.any():
                a_nadir[fin] = _incidence_angle_deg_vec(d[fin], refine_config.wgs84_re_km, rs_near[fin])
            if exclude_caps and (~fin).any():
                # The end-cap observations the v1.5.0 rule would have kept,
                # found with that rule itself (legacy distance and the radius of
                # its nearest sample); counted below, after the daylight filter.
                d_old, rs_old = _min_track_distance_km(s_lat, s_lon, s_rs, c_lat[~fin], c_lon[~fin],
                                                       refine_config.wgs84_re_km, exclude_caps=False)
                th_end = float(_incidence_angle_deg_vec(meta["swath_km"] / 2.0, refine_config.wgs84_re_km,
                                                        float(np.mean(s_rs))))
                cap_ok = np.zeros(cand.size, dtype=bool)
                cap_ok[~fin] = _incidence_angle_deg_vec(d_old, refine_config.wgs84_re_km, rs_old) <= th_end
        else:
            # Starlink proxy: nearest strip sample, unchanged so earlier results
            # reproduce exactly (swath-edge error ~0.3 km against a ~80 km footprint).
            cand = _tiles_near(s_lat, s_lon, sensor_config.footprint_radius_km)
            if cand.size == 0:
                continue
            c_lat, c_lon = _tile_lat[cand], _tile_lon[cand]
            a_nadir = np.full(cand.size, np.inf)
            for la, lo, rs in zip(s_lat, s_lon, s_rs):
                d = _hav_km(la, lo, c_lat, c_lon)
                a_nadir = np.minimum(a_nadir, _incidence_angle_deg_vec(d, refine_config.wgs84_re_km, rs))
        if meta.get("swath_km"):
            # Off-nadir angle of this satellite's swath edge (half the swath on
            # the ground) at its mean geocentric radius over the strip.
            theta = float(_incidence_angle_deg_vec(meta["swath_km"] / 2.0, refine_config.wgs84_re_km,
                                                   float(np.mean(s_rs))))
        else:
            theta = swath_half_angle_deg
        within_fov = a_nadir <= theta
        if cap_ok is not None and cap_ok.any() and excluded_by is None:
            # v1.5.2: counted only where the v1.5.0 rule would also have kept
            # the row (sunlit), so the count equals the rows it would add.
            if meta.get("min_sun_elev_deg") is not None:
                _t_cap = _t0_utc + np.timedelta64(int(round(t_min * 60)), "s")
                _sun_cap = solar_elevation_deg(_t_cap, c_lat[cap_ok], c_lon[cap_ok])
                refine_daylight_stats["rows_end_cap_excluded"] += int((_sun_cap >= meta["min_sun_elev_deg"]).sum())
            else:
                refine_daylight_stats["rows_end_cap_excluded"] += int(cap_ok.sum())
        if not within_fov.any():
            continue
        sun = np.full(cand.size, np.nan)
        if meta.get("min_sun_elev_deg") is not None:
            # Optical imagers: drop tiles where the Sun is below the mission's
            # limit at the event minute (NOAA solar position, max error 0.46 deg).
            t_utc = _t0_utc + np.timedelta64(int(round(t_min * 60)), "s")
            sun[within_fov] = solar_elevation_deg(t_utc, c_lat[within_fov], c_lon[within_fov])
            lit = within_fov & (sun >= meta["min_sun_elev_deg"])
            refine_daylight_stats["rows_dark"] += int(within_fov.sum() - lit.sum())
            refine_daylight_stats["rows_lit"] += int(lit.sum())
            if not lit.any():
                refine_daylight_stats["events_dark"] += 1
                continue
            within_fov = lit
        if excluded_by is not None:
            # v1.5.2: inside a sourced acquisition exclusion (Cell 0): the
            # satellite is tracked but did not acquire, so the pass is dropped.
            _n_ex = int(within_fov.sum())
            refine_daylight_stats["rows_acquisition_excluded"] += _n_ex
            refine_daylight_stats["passes_acquisition_excluded"] += 1
            _hit = refine_daylight_stats["acquisition_exclusions_hit"].setdefault(str(excluded_by), {"passes": 0, "rows": 0})
            _hit["passes"] += 1
            _hit["rows"] += _n_ex
            continue

        # Downlink Reward (constant per event)
        d_wait = time_to_next_downlink(sat_id, t_min)
        d_reward = downlink_reward(d_wait)
        tile_chunks.append(cand[within_fov].astype(np.int32))
        angle_chunks.append(a_nadir[within_fov])
        sun_chunks.append(sun[within_fov])
        ev_theta.append(theta)
        ev_sat.append(sat_id); ev_t.append(t_min); ev_n.append(int(within_fov.sum()))
        ev_d.append(float(d_reward)); ev_w.append(np.nan if d_wait is None else float(d_wait))

    if len(tile_chunks) == 0:
        logger.warning("No tile events survived the physical FOV filter.")
        return pd.DataFrame(columns=["sat", "tile", "t_min", "A_nadir_deg", "D_downlink", "T_wait_min",
                                     "theta_max_deg", "sun_elev_deg"])

    # Compact columns: categorical sat/tile (tile categories = Cell 4's _tile_ids,
    # so codes are tile indices), int32 minutes, float64 physics values.
    counts = np.array(ev_n, dtype=np.int64)
    sat_cats = sorted(set(ev_sat))
    sat_code = {s: i for i, s in enumerate(sat_cats)}
    return pd.DataFrame({
        "sat": pd.Categorical.from_codes(np.repeat(np.array([sat_code[s] for s in ev_sat], dtype=np.int32), counts),
                                         categories=sat_cats),
        "tile": pd.Categorical.from_codes(np.concatenate(tile_chunks), categories=_tile_ids),
        "t_min": np.repeat(np.array(ev_t, dtype=np.int32), counts),
        "A_nadir_deg": np.concatenate(angle_chunks).astype(np.float64),
        "D_downlink": np.repeat(np.array(ev_d, dtype=np.float64), counts),
        "T_wait_min": np.repeat(np.array(ev_w, dtype=np.float64), counts),
        "theta_max_deg": np.repeat(np.array(ev_theta, dtype=np.float64), counts),   # swath-edge angle, this satellite
        "sun_elev_deg": np.concatenate(sun_chunks).astype(np.float64),             # NaN: no daylight limit
    })

# =============================================================================
# EXECUTION & POST-REFINEMENT ANALYSIS
# =============================================================================

assert 'access_events_df' in globals(), "access_events_df missing. Execute Cell 4."
assert 'constellation' in globals(), "constellation missing. Execute Cell 2."

logger.info("\nRefining continuous access strips with CPA incidence angles and downlink rewards...")
t0 = time.time()
refined_events_df = refine_access_events(
    access_events_df=access_events_df,
    constellation=constellation,
    sim_epoch=sim_epoch
)
t_refine = time.time() - t0
logger.info(f"Refinement complete in {t_refine:.2f}s ({len(refined_events_df):,} verified tile-access records).")

# =============================================================================
# EXPORT & SUMMARY
# =============================================================================

_refined_csv = os.path.join(_out_dir, 'arf_refined_events.csv')
_downlink_csv = os.path.join(_out_dir, 'arf_downlink_events.csv')
_summary_json = os.path.join(_out_dir, 'arf_cell5_summary.json')

try:
    refined_events_df.to_csv(_refined_csv, index=False)
    downlink_events_df.to_csv(_downlink_csv, index=False)
    logger.info(f"Refined datasets exported to {_refined_csv} and {_downlink_csv}")
except Exception as e:
    logger.warning(f"Could not save refined event CSVs: {e}")

_summary = {
    "n_access_events_input": len(access_events_df),
    "n_refined_events": len(refined_events_df),
    "n_unique_tiles": refined_events_df["tile"].nunique() if len(refined_events_df) > 0 else 0,
    "n_unique_sats": refined_events_df["sat"].nunique() if len(refined_events_df) > 0 else 0,
    "a_nadir_mean": float(refined_events_df['A_nadir_deg'].mean()) if len(refined_events_df) > 0 else 0.0,
    "d_downlink_mean": float(refined_events_df['D_downlink'].mean()) if len(refined_events_df) > 0 else 0.0,
    "n_downlink_windows": len(downlink_events_df),
    "refinement_time_s": t_refine,
    # v1.5.1: each tile observation belongs to the strip of its closest approach
    "strip_end_caps": "legacy (v1.5.0 and earlier)" if globals().get("LEGACY_STRIP_END_CAPS", False) else "excluded",
    # v1.5.2: counted after the daylight filter = rows the v1.5.0 rule would add
    "rows_end_cap_excluded": int(refine_daylight_stats.get("rows_end_cap_excluded", 0)),
    "rows_below_sun_limit": int(refine_daylight_stats.get("rows_dark", 0)),
    # v1.5.2: passes dropped by sourced acquisition exclusions (Cell 0)
    "rows_acquisition_excluded": int(refine_daylight_stats.get("rows_acquisition_excluded", 0)),
    "passes_acquisition_excluded": int(refine_daylight_stats.get("passes_acquisition_excluded", 0)),
    "acquisition_exclusions_hit": [
        {**{k: globals()["ACQUISITION_EXCLUSIONS"][int(i)][k] for k in ("name_pattern", "start_utc", "end_utc",
                                                                         "source_url")}, **v}
        for i, v in sorted(refine_daylight_stats.get("acquisition_exclusions_hit", {}).items())],
    "theta_max_deg_range": ([float(refined_events_df["theta_max_deg"].min()), float(refined_events_df["theta_max_deg"].max())]
                            if len(refined_events_df) > 0 else None),
}

try:
    with open(_summary_json, 'w') as f:
        json.dump(_summary, f, indent=2, default=str)
    logger.info(f"Cell 5 summary exported to {_summary_json}")
except Exception as e:
    logger.warning(f"Could not save Cell 5 summary: {e}")

print("\n" + "=" * 60)
print("  CELL 5: Physical Constraints & Refinement Complete")
print("=" * 60)
print(f"  Input Access Strips:       {len(access_events_df):,}")
print(f"  Refined Tile Records:      {len(refined_events_df):,}")
print(f"  Unique Valid Tiles Hit:    {_summary['n_unique_tiles']:,} / {len(target_tiles):,}")
# v1.5.2: real fleets have one swath-edge angle per satellite (Sentinel-2 ~10.3
# deg, Landsat ~7.5 deg, Doves ~1.5-2 deg); 9.22 deg is the Starlink proxy's.
_theta_txt = (f"swath edge {_summary['theta_max_deg_range'][0]:.2f}-{_summary['theta_max_deg_range'][1]:.2f}° by satellite"
              if _real_fleet5 and _summary["theta_max_deg_range"] else f"Max: {refine_config.swath_half_angle_deg}°")
print(f"  Mean Incidence Angle:      {_summary['a_nadir_mean']:.2f}° ({_theta_txt})")
print(f"  Mean Downlink Reward:      {_summary['d_downlink_mean']:.3f} / 1.000")
print(f"  Downlink Contact Windows:  {len(downlink_events_df):,}")
if _summary["rows_end_cap_excluded"]:
    print(f"  End-cap rows not counted:  {_summary['rows_end_cap_excluded']:,} (v1.5.0 rule; each tile kept in the "
          f"strip of its closest approach)")
if _summary["passes_acquisition_excluded"]:
    print(f"  Acquisition exclusions:    {_summary['passes_acquisition_excluded']} passes, "
          f"{_summary['rows_acquisition_excluded']:,} rows dropped (sourced, Cell 0)")
print("=" * 60)
print("\n=== CELL 5 COMPLETE ===")

