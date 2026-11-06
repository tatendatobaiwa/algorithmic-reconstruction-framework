# =============================================================================
# CELL 4: SPATIO-TEMPORAL ACCESS DETECTION (Coarse-to-Fine Scan)
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Implements Section 3.3.1 of the research proposal (Footprint Projection).
#   Discretises continuous SGP4 orbital ground tracks onto the H3 hexagonal grid
#   to construct the candidate Spatio-Temporal Graph event set Φ.
#
# Algorithmic Efficiency & Kinematic Correction:
#   Avoids the O(N x T x M) brute-force scan with a two-stage search:
#     1. Stage 1 (Coarse, 5-min): flags a satellite whose sub-satellite point
#        falls inside the ROI's bounding box inflated by the footprint radius
#        plus the distance a satellite can travel in HALF a coarse step (~1,100
#        km at LEO). The nearest coarse sample is then guaranteed to fall inside
#        the box for every pass, so no pass can slip between samples. Longitude
#        margins are scaled by cos(latitude); boxes are antimeridian-safe.
#     2. Stage 2 (Strips, 1-min): each whole minute near a coarse hit is a
#        candidate event (sat, m) = the continuous push-broom swath swept over
#        [m-0.5, m+0.5) min, sampled every 2 s. Tile distances are computed
#        directly against tile centres (vectorised), not by enumerating H3 rings.
#   Search radii are deliberate supersets (highest satellite's footprint,
#   lowest satellite's ground speed, padded); Cell 5 applies the exact
#   incidence-angle filter. Validated against a 0.25-s ground truth on
#   synthetic 340-570 km constellations (43/53/70/97.6 deg): 99.9% tile-access
#   recall, zero false positives; residual misses lie within ~0.07 deg of
#   the swath edge (2-s sampling tolerance).
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 4 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["4"] = _ARF_CELL_VERSION

import os
import sys
import time
import logging
from typing import Dict, List, Tuple, Set, Optional, Any
from dataclasses import dataclass
import numpy as np
import pandas as pd
import h3
from org.hipparchus.util import FastMath

logger = logging.getLogger("arf.cell4")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

def _footprint_ground_radius_km(altitude_km, half_angle_deg):
    """
    Ground-range radius of a nadir-pointing sensor's footprint on a
    spherical Earth: Re * (asin((Re+h)/Re * sin(eta)) - eta). Same
    spherical-Earth model as compute_incidence (Cell 5), so this search
    radius and Cell 5's exact filter agree on where the swath edge is.
    """
    eta = np.radians(half_angle_deg)
    s = min(1.0, (Re_km + altitude_km) / Re_km * np.sin(eta))
    return Re_km * (np.arcsin(s) - eta)


@dataclass
class SensorConfig:
    swath_half_angle_deg: float = SWATH_HALF_ANGLE_DEG
    max_altitude_km: float = CONSTELLATION_ALTITUDE_KM
    min_altitude_km: float = CONSTELLATION_ALTITUDE_KM
    window_hours: int = WINDOW_HOURS
    coarse_step_minutes: int = 5
    # 2 s along-track sampling (~14 km) caps swath-edge erosion from
    # sampling at ~0.3 km -- about a fifth of one Res-7 tile edge.
    strip_sample_seconds: float = 2.0
    # Real fleets: the WIDEST swath in the fleet (km on the ground). The
    # search radius is then a superset for every satellite; Cell 5 applies
    # each satellite's own swath exactly.
    max_swath_km: Optional[float] = None

    def __post_init__(self):
        # Fastest ground-relative speed: the LOWEST satellite's ground-track
        # speed plus Earth's equatorial surface rotation speed.
        h = max(150.0, self.min_altitude_km - 15.0)
        v_orb = np.sqrt(mu_earth / (Re_km + h)) * 60.0 * Re_km / (Re_km + h)
        v_rot = 2.0 * np.pi * Re_km / 1436.068
        self.max_ground_speed_km_min = v_orb + v_rot
        if self.max_swath_km:
            # Widest half-swath +5%, plus HALF the along-track sample spacing:
            # a tile lying between two 2-s samples (~15 km apart) can be up to
            # that far from both while still on the track. Without this term
            # an 8 km Dove swath loses most tiles between samples (measured).
            self.footprint_radius_km = (1.05 * self.max_swath_km / 2.0
                                        + 0.5 * self.max_ground_speed_km_min * self.strip_sample_seconds / 60.0)
        else:
            # Footprint of the HIGHEST satellite, +15 km altitude for orbital
            # eccentricity over the window, +5%.
            self.footprint_radius_km = 1.05 * _footprint_ground_radius_km(
                self.max_altitude_km + 15.0, self.swath_half_angle_deg)
        self.coarse_margin_km = (self.footprint_radius_km
                                 + self.max_ground_speed_km_min * self.coarse_step_minutes / 2.0)
        self.strip_margin_km = self.footprint_radius_km + self.max_ground_speed_km_min * 0.5


_alts = pos_list if ('pos_list' in globals() and len(pos_list)) else [CONSTELLATION_ALTITUDE_KM]
_swaths = [s["swath_km"] for s in constellation if s.get("swath_km")] if 'constellation' in globals() else []
sensor_config = SensorConfig(max_altitude_km=float(np.max(_alts)), min_altitude_km=float(np.min(_alts)),
                             max_swath_km=(float(max(_swaths)) if _swaths else None))
logger.info(f"Footprint search radius: {sensor_config.footprint_radius_km:.1f} km | "
            f"coarse-box margin: {sensor_config.coarse_margin_km:.0f} km | "
            f"strip margin: {sensor_config.strip_margin_km:.0f} km")

# --- Geometry helpers --------------------------------------------------------
_KM_PER_DEG_LAT_MIN = 110.574   # WGS84 meridian, at the equator (smallest -> conservative)
_KM_PER_DEG_LON_EQ = 111.320    # WGS84 parallel, at the equator


def _wrap_dlon(lon, lon_ref):
    """Longitude difference wrapped into [-180, 180) -- antimeridian-safe."""
    return (np.asarray(lon, dtype=float) - lon_ref + 180.0) % 360.0 - 180.0


def _hav_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(_wrap_dlon(lon2, 0) - _wrap_dlon(lon1, 0)) / 2) ** 2)
    return 2.0 * Re_km * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


class _GeoBox:
    """
    Lat/lon box around a set of points, inflated by margin_km. The
    longitude margin is scaled by cos(latitude) at the box's widest-
    latitude edge (where a degree of longitude is shortest), longitudes are
    handled relative to a centre so a box may straddle +/-180, and a box
    reaching within 1 deg of a pole accepts all longitudes.
    """
    def __init__(self, lats, lons, margin_km):
        lats = np.asarray(lats, dtype=float)
        mlat = margin_km / _KM_PER_DEG_LAT_MIN
        self.lat_lo, self.lat_hi = lats.min() - mlat, lats.max() + mlat
        lat_ref = max(abs(self.lat_lo), abs(self.lat_hi))
        ref = float(np.asarray(lons, dtype=float).flat[0])
        d = _wrap_dlon(lons, ref)
        self.lon_c = ref + (d.min() + d.max()) / 2.0
        self.all_lon = lat_ref >= 89.0
        self.lon_half = (d.max() - d.min()) / 2.0
        if not self.all_lon:
            self.lon_half += margin_km / (_KM_PER_DEG_LON_EQ * np.cos(np.radians(lat_ref)))
            self.all_lon = self.lon_half >= 180.0

    def contains(self, lat, lon):
        lat = np.asarray(lat, dtype=float)
        ok = (lat >= self.lat_lo) & (lat <= self.lat_hi)
        return ok if self.all_lon else ok & (np.abs(_wrap_dlon(lon, self.lon_c)) <= self.lon_half)


# ROI tile centres, sorted by latitude so candidate tiles near a strip can be
# found with a binary search instead of a scan of every tile.
_tile_ids = np.array(sorted(target_tiles))
_ll = np.array([h3.cell_to_latlng(t) for t in _tile_ids])
_order = np.argsort(_ll[:, 0], kind="stable")
_tile_ids, _tile_lat, _tile_lon = _tile_ids[_order], _ll[_order, 0], _ll[_order, 1]
del _ll, _order


def _tiles_near(lats, lons, radius_km):
    """Indices of ROI tiles whose centre lies within radius_km (ground
    range) of at least one of the given sub-satellite points."""
    box = _GeoBox(lats, lons, radius_km)
    lo = np.searchsorted(_tile_lat, box.lat_lo, side="left")
    hi = np.searchsorted(_tile_lat, box.lat_hi, side="right")
    idx = np.arange(lo, hi)
    if idx.size and not box.all_lon:
        idx = idx[np.abs(_wrap_dlon(_tile_lon[idx], box.lon_c)) <= box.lon_half]
    if idx.size == 0:
        return idx
    d = _hav_km(np.asarray(lats)[:, None], np.asarray(lons)[:, None],
                _tile_lat[idx][None, :], _tile_lon[idx][None, :]).min(axis=0)
    return idx[d <= radius_km]


def _subsatellite(prop, t_min):
    """(geodetic lat deg, lon deg, geocentric radius km) at t_min."""
    t = sim_epoch.shiftedBy(float(t_min * 60.0))
    pos = prop.getPVCoordinates(t, itrf).getPosition()
    geo = earth.transform(pos, itrf, t)
    return (FastMath.toDegrees(geo.getLatitude()), FastMath.toDegrees(geo.getLongitude()),
            pos.getNorm() / 1000.0)


_roi_box = _GeoBox(_tile_lat, _tile_lon, 0.0)
coarse_box = _GeoBox(_tile_lat, _tile_lon, sensor_config.coarse_margin_km)
strip_box = _GeoBox(_tile_lat, _tile_lon, sensor_config.strip_margin_km)

# (sat, t_min) -> (lats, lons, radii_km) of the strip's sub-satellite samples;
# Cell 5 reuses these instead of re-propagating.
strip_samples = {}

_ACCESS_COLUMNS = ["sat", "t_min", "lat", "lon", "tiles_hit"]


def detect_access_events(constellation, config):
    win_min = config.window_hours * 60
    step = config.coarse_step_minutes
    coarse_times = list(range(0, win_min + 1, step))
    if coarse_times[-1] != win_min:
        coarse_times.append(win_min)
    offsets = np.arange(-0.5, 0.5, config.strip_sample_seconds / 60.0)
    if config.max_swath_km:
        # Real fleets: include the strip's end point (= next strip's first
        # sample) so the track arcs of consecutive strips join with no gap.
        offsets = np.append(offsets, 0.5)

    events = []
    n_coarse_hits = n_candidate_min = n_strips = 0
    for sat in constellation:
        prop = sat["propagator"]
        # Stage 1: coarse scan against the inflated box.
        hits = []
        for m in coarse_times:
            lat, lon, _ = _subsatellite(prop, m)
            if coarse_box.contains(lat, lon):
                hits.append(m)
        n_coarse_hits += len(hits)
        if not hits:
            continue
        # Stage 2: every minute within one coarse step of a hit.
        minutes = sorted({mm for m_c in hits
                          for mm in range(max(0, m_c - step), min(win_min, m_c + step) + 1)})
        n_candidate_min += len(minutes)
        for m in minutes:
            mid_lat, mid_lon, _ = _subsatellite(prop, m)
            if not strip_box.contains(mid_lat, mid_lon):
                continue            # no part of this minute's strip can reach the ROI
            n_strips += 1
            ts = np.clip(m + offsets, 0.0, float(win_min))
            s = np.array([_subsatellite(prop, t) for t in ts])
            near = _tiles_near(s[:, 0], s[:, 1], config.footprint_radius_km)
            if near.size:
                strip_samples[(sat["id"], m)] = (s[:, 0], s[:, 1], s[:, 2])
                events.append({"sat": sat["id"], "t_min": m, "lat": mid_lat, "lon": mid_lon,
                               "tiles_hit": int(near.size)})

    logger.info(f"Stage 1: {n_coarse_hits} coarse hits | Stage 2: {n_candidate_min} candidate "
                f"minutes, {n_strips} strips sampled, {len(events)} events")
    return pd.DataFrame(events, columns=_ACCESS_COLUMNS)


logger.info(f"Starting Access Scan for {len(constellation)} satellites...")
t0 = time.time()
access_events_df = detect_access_events(constellation, sensor_config)
_t_scan = time.time() - t0
logger.info(f"Scan Complete in {_t_scan:.2f}s. "
            f"Found {len(access_events_df)} continuous access strips.")
display(access_events_df.head())

import os, json
_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()
access_events_df.to_csv(os.path.join(_out_dir, "arf_access_events.csv"), index=False)
with open(os.path.join(_out_dir, "arf_cell4_summary.json"), "w") as f:
    json.dump({
        "n_satellites": len(constellation),
        "n_access_strips": int(len(access_events_df)),
        "footprint_search_radius_km": float(sensor_config.footprint_radius_km),
        "coarse_box_margin_km": float(sensor_config.coarse_margin_km),
        "strip_margin_km": float(sensor_config.strip_margin_km),
        "coarse_step_minutes": sensor_config.coarse_step_minutes,
        "strip_sample_seconds": sensor_config.strip_sample_seconds,
        "altitude_range_km": [float(sensor_config.min_altitude_km), float(sensor_config.max_altitude_km)],
        "scan_seconds": round(_t_scan, 2),
    }, f, indent=2)

