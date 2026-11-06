# =============================================================================
# CELL 3: SPATIAL VALIDATION & GROUND TRACK VISUALISATION
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Study-area map of the active ROI, drawn to journal standard:
#     (a) Regional setting -- the ROI, the ground-station network,
#         international borders and sample ground tracks of the fleet.
#     (b) The ROI itself -- its H3 tiling (the grid the mosaic is built on)
#         and its official boundary, with a scale bar and a north arrow.
#   Panel (b) is fitted to the ROI, so the mosaic is legible at every ROI
#   size, from a 32-tile reserve to an 18,000-tile complex; ground stations
#   far from the ROI appear in panel (a) only. Large tile sets get a zoom
#   inset showing the individual hexagons. Both panels are in geographic
#   coordinates (WGS 84), drawn with the equirectangular aspect of their
#   centre latitude (1/cos(latitude)), so shapes and distances are not
#   stretched.
#
#   The map helpers defined here are reused by Cell 7's mosaic map.
#
# Research Proposal Mapping:
#   - Section 3.2.1: Geodesic Grid Model (Uber H3 Tiles)
#   - Section 3.2.3: Ground Station Telemetry Anchors
#   - Section 4.1: Simulation Environment Validation
#
# Dependencies:
#   - Cell 0 (ACTIVE_N, ACTIVE_ROI, OUT_DIR, FLEET)
#   - Cell 1 (target_tiles, GROUND_STATIONS, ROI_OUTLINE_GEOMS, itrf, earth)
#   - Cell 2 (constellation, sim_epoch)
#
# Outputs (names unchanged):
#   - arf_coverage_map.png (300 dpi)
#   - arf_coverage_map.pdf (vector; tile layers above MAP_RASTER_TILES tiles
#     are embedded as a 300 dpi image to keep the file small)
#   - arf_cell3_metadata.json (reproducibility metadata and a draft caption)
#
# Data drawn on the map:
#   - ROI boundaries: World Database on Protected Areas (UNEP-WCMC and IUCN).
#     The WDPA/WDOECM Manual 1.6 (section 5.3, p. 30) asks that its citation
#     be reproduced in any publication using the WDPA; the draft caption in
#     arf_cell3_metadata.json carries it. Boundaries are drawn on the figure
#     but never written to the output files (the licence forbids
#     redistributing the data).
#   - International borders: Natural Earth 1:50m Admin 0 - Countries, v5.1.2,
#     public domain ("No permission is needed to use Natural Earth",
#     https://www.naturalearthdata.com/about/terms-of-use/). Downloaded once
#     and cached; if it cannot be downloaded the map is drawn without borders.
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

"""
ARF Spatial Validation Visualisation Module
============================================

1. **ROI mosaic grid**: the H3 cells of the ROI, filled as one dissolved
   shape (no seams between cells) with the cell outlines drawn on top while
   they are larger than a pixel; a zoom inset shows them for large ROIs.
2. **Ground stations**: the downlink stations of Cell 1, with labels.
3. **Ground tracks**: SGP4 ground tracks of a deterministic sample of the
   fleet over the first MAP_TRACK_HOURS hours, keeping the first
   MAP_TRACK_MAX_SEGMENTS crossings of the regional frame (arrows give the
   direction of flight).
"""

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 3 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["3"] = _ARF_CELL_VERSION

import os
import sys
import json
import math
import hashlib
import logging
from datetime import datetime, timezone
from typing import Dict, List, Tuple, Set, Optional, Any

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms
import matplotlib.lines as mlines
import matplotlib.patheffects as patheffects
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.path import Path as MplPath
import h3
from org.hipparchus.util import FastMath

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logger = logging.getLogger("arf.cell3")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()

# =============================================================================
# MATPLOTLIB PUBLICATION CONFIGURATION
# =============================================================================
# Times New Roman where installed, otherwise a metric-compatible Times clone
# (e.g. Liberation Serif); DejaVu Serif is the last resort. pdf.fonttype 42
# embeds the fonts in PDFs as TrueType (Type 42) instead of matplotlib's
# default Type 3, so figure text stays embedded, selectable text.

ARF_SERIF_FONTS = ['Times New Roman', 'Liberation Serif', 'Nimbus Roman', 'FreeSerif', 'DejaVu Serif']

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ARF_SERIF_FONTS,
    'font.size': 10,
    'mathtext.fontset': 'stix',
    'axes.linewidth': 1.0,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
    'figure.dpi': 100,
    'figure.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# =============================================================================
# MAP HELPERS (shared with Cell 7's mosaic map)
# =============================================================================

MAP_FIG_WIDTH_IN = 7.0          # double-column figure width (about 178 mm)
MAP_ASPECT_RANGE = (0.8, 1.5)   # panel width/height is kept inside this range
MAP_RASTER_TILES = 3000         # tile layers with more tiles are rasterised in the PDF
MAP_EDGE_TILES = 3000           # cell outlines are drawn up to this many tiles (beyond: sub-pixel)
MAP_ZOOM_TILES = 1500           # above this many tiles, a zoom inset shows the hexagons
MAP_TRACK_HOURS = 24            # ground tracks: hours from the analysis start ...
MAP_TRACK_MAX_SATS = 40         # ... for an evenly spaced sample of this many satellites ...
MAP_TRACK_MAX_SEGMENTS = 8      # ... keeping the first crossings of the regional frame

MAP_COLORS = {
    "tile": "#41ab5d",          # ROI cells
    "tile_edge": "white",
    "boundary": "#1a1a1a",
    "station": "#b2182b",
    "track": "#525252",
    "extent": "#1a1a1a",
    "land": "#f0efe9",
    "water": "#e4eef6",
    "border": "#a0a0a0",
    "country_label": "#808080",
}

# Map text sizes at the printed figure width (journals ask for >= 6-7 pt).
_MAP_RC = {
    "font.family": "serif",
    "font.serif": ARF_SERIF_FONTS,
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

_HALO = [patheffects.withStroke(linewidth=2.2, foreground="white")]

NATURAL_EARTH_VERSION = "v5.1.2"
NATURAL_EARTH_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
                     + NATURAL_EARTH_VERSION + "/50m_cultural/ne_50m_admin_0_countries.{ext}")
_NATURAL_EARTH_BASE = "ne_50m_admin_0_countries"
_country_shape_cache: Dict[str, Any] = {}


def _km_per_deg_lon(lat_deg: float) -> float:
    """Length of one degree of longitude along the parallel at lat_deg on the
    WGS 84 ellipsoid (a = 6378.137 km, e^2 = 6.69437999014e-3): (pi/180) N cos(lat)."""
    a_km, e2 = 6378.137, 6.69437999014e-3
    s = math.sin(math.radians(lat_deg))
    return math.radians(1.0) * a_km * math.cos(math.radians(lat_deg)) / math.sqrt(1.0 - e2 * s * s)


_KM_PER_DEG_LAT = 110.574       # WGS 84 meridian degree at the equator (shortest; used only for padding)


def _cos_lat(ext) -> float:
    return max(math.cos(math.radians(0.5 * (ext[2] + ext[3]))), 0.05)


def panel_aspect(ext) -> float:
    """Displayed width/height of a map panel with extent [lon0, lon1, lat0, lat1]."""
    return (ext[1] - ext[0]) * _cos_lat(ext) / (ext[3] - ext[2])


def pad_extent(ext, frac: float = 0.06, min_pad_km: float = 1.0):
    """Grow an extent on all sides by frac of its larger side (at least min_pad_km)."""
    lon0, lon1, lat0, lat1 = ext
    mid = 0.5 * (lat0 + lat1)
    w_km = (lon1 - lon0) * _km_per_deg_lon(mid)
    h_km = (lat1 - lat0) * _KM_PER_DEG_LAT
    pad_km = max(frac * max(w_km, h_km), min_pad_km)
    dlon = pad_km / max(_km_per_deg_lon(mid), 1e-6)
    dlat = pad_km / _KM_PER_DEG_LAT
    return [lon0 - dlon, lon1 + dlon, max(lat0 - dlat, -89.0), min(lat1 + dlat, 89.0)]


def fit_aspect(ext, aspect: float):
    """Grow the shorter side of an extent (about its centre) to the given displayed aspect."""
    lon0, lon1, lat0, lat1 = ext
    if panel_aspect(ext) < aspect:
        half = 0.5 * (lat1 - lat0) * aspect / _cos_lat(ext)
        cx = 0.5 * (lon0 + lon1)
        return [cx - half, cx + half, lat0, lat1]
    half = 0.5 * (lon1 - lon0) * _cos_lat(ext) / aspect
    cy = 0.5 * (lat0 + lat1)
    return [lon0, lon1, cy - half, cy + half]


def _ring_area(ring) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def geo_polygons(geoms) -> List[Tuple[np.ndarray, List[np.ndarray]]]:
    """[(exterior, [holes])] as (lon, lat) arrays from GeoJSON-like Polygon /
    MultiPolygon dicts (pyshp __geo_interface__, h3.cells_to_geo)."""
    out = []
    for g in geoms or []:
        if not g:
            continue
        if g["type"] == "Polygon":
            polys = [g["coordinates"]]
        elif g["type"] == "MultiPolygon":
            polys = g["coordinates"]
        else:
            continue
        for poly in polys:
            rings = []
            for r in poly:
                if len(r) < 3:
                    continue
                r = np.asarray(r, dtype=float)[:, :2]
                if not np.allclose(r[0], r[-1]):          # close open rings (e.g. ROI corner lists)
                    r = np.vstack([r, r[:1]])
                rings.append(r)
            if rings:
                out.append((rings[0], rings[1:]))
    return out


def geo_lines(geoms) -> List[np.ndarray]:
    """Every ring (exteriors and holes) of the geometries, for outlines."""
    return [r for ext_ring, holes in geo_polygons(geoms) for r in [ext_ring] + list(holes)]


def geo_path(geoms) -> Optional[MplPath]:
    """One compound Path for filling (exteriors anticlockwise, holes clockwise)."""
    verts, codes = [], []
    for ext_ring, holes in geo_polygons(geoms):
        for i, ring in enumerate([ext_ring] + list(holes)):
            if (_ring_area(ring) > 0) != (i == 0):
                ring = ring[::-1]
            verts.extend(ring.tolist())
            verts.append(ring[0].tolist())
            codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(ring) - 1) + [MplPath.CLOSEPOLY])
    return MplPath(verts, codes) if verts else None


def tiles_outline(tiles) -> Optional[dict]:
    """Dissolved outline of a set of H3 cells (GeoJSON-like, lon/lat)."""
    tiles = list(tiles)
    return h3.cells_to_geo(tiles) if tiles else None


def tile_polygons(tiles) -> List[np.ndarray]:
    """Hexagon vertices (lon, lat) of each H3 cell, in the order given."""
    polys = []
    for t in tiles:
        p = np.array([(lng, lat) for lat, lng in h3.cell_to_boundary(t)], dtype=float)
        if np.ptp(p[:, 0]) > 180.0:                      # cell straddles the antimeridian
            p[:, 0] = np.where(p[:, 0] < 0.0, p[:, 0] + 360.0, p[:, 0])
        polys.append(p)
    return polys


def geo_bounds(geoms=None, points=()):
    """[lon0, lon1, lat0, lat1] of geometries and (lat, lon) points."""
    xs, ys = [], []
    for r in geo_lines(geoms):
        xs.append(r[:, 0])
        ys.append(r[:, 1])
    for lat, lon in points:
        xs.append([lon])
        ys.append([lat])
    x, y = np.concatenate(xs), np.concatenate(ys)
    return [float(x.min()), float(x.max()), float(y.min()), float(y.max())]


def roi_map_extent(tiles, outline_geoms=None, pad_frac: float = 0.06):
    """Map extent fitted to the ROI (its cells and boundary), padded, with the
    displayed aspect clamped to MAP_ASPECT_RANGE."""
    ext = pad_extent(geo_bounds([tiles_outline(tiles)] + list(outline_geoms or [])), pad_frac, min_pad_km=1.0)
    lo, hi = MAP_ASPECT_RANGE
    return fit_aspect(ext, min(max(panel_aspect(ext), lo), hi))


def regional_extent(detail_ext, stations, aspect: float, context: float = 3.0,
                    min_span_km: float = 900.0, pad_frac: float = 0.12):
    """Extent for the regional panel: the ROI panel and every ground station,
    padded, at least `context` times the ROI panel and min_span_km across,
    fitted to the same displayed aspect as the ROI panel."""
    xs = [detail_ext[0], detail_ext[1]] + [gs["lng"] for gs in stations.values()]
    ys = [detail_ext[2], detail_ext[3]] + [gs["lat"] for gs in stations.values()]
    ext = pad_extent([min(xs), max(xs), min(ys), max(ys)], pad_frac, min_pad_km=50.0)
    mid = 0.5 * (ext[2] + ext[3])
    need_lon = max(context * (detail_ext[1] - detail_ext[0]), min_span_km / max(_km_per_deg_lon(mid), 1e-6))
    need_lat = max(context * (detail_ext[3] - detail_ext[2]), min_span_km / _KM_PER_DEG_LAT)
    cx, cy = 0.5 * (ext[0] + ext[1]), 0.5 * (ext[2] + ext[3])
    half_lon = 0.5 * max(ext[1] - ext[0], need_lon)
    half_lat = 0.5 * max(ext[3] - ext[2], need_lat)
    ext = fit_aspect([cx - half_lon, cx + half_lon, cy - half_lat, cy + half_lat], aspect)
    if ext[2] < -85.0 or ext[3] > 85.0:                  # keep clear of the poles
        shift = (-85.0 - ext[2]) if ext[2] < -85.0 else (85.0 - ext[3])
        ext = [ext[0], ext[1], max(ext[2] + shift, -85.0), min(ext[3] + shift, 85.0)]
    return ext


_TICK_STEPS_DEG = (1 / 60, 2 / 60, 5 / 60, 10 / 60, 15 / 60, 20 / 60, 30 / 60, 1, 2, 5, 10, 15, 20, 30)


def _fmt_deg(v: float, step: float, pos: str, neg: str) -> str:
    """24°E, 24°30′E, 32°10′S ... (degrees, or degrees and minutes when step < 1°)."""
    hemi = pos if v > 1e-9 else (neg if v < -1e-9 else "")
    a = abs(v)
    if step >= 1.0 - 1e-9:
        return f"{a:.0f}°{hemi}"
    d = int(math.floor(a + 1e-9))
    m = int(round((a - d) * 60.0))
    if m == 60:
        d, m = d + 1, 0
    return f"{d}°{m:02d}′{hemi}" if m else f"{d}°{hemi}"


def _ticks(lo: float, hi: float, max_intervals: int) -> Tuple[np.ndarray, float]:
    """Graticule ticks at a 'nice' step, with at least two ticks inside [lo, hi]."""
    steps = [s for s in _TICK_STEPS_DEG if (hi - lo) / s <= max_intervals] or [30.0]
    step = steps[0]
    for s in reversed(_TICK_STEPS_DEG[:_TICK_STEPS_DEG.index(step) + 1] if step in _TICK_STEPS_DEG else ()):
        t = np.arange(math.ceil(lo / s - 1e-9) * s, hi + 1e-9, s)
        if t.size >= 2:
            return t, s
    return np.arange(math.ceil(lo / step - 1e-9) * step, hi + 1e-9, step), step


def setup_map_axes(ax, ext, max_intervals=(4, 4), water: Optional[str] = None, labels=(True, True)):
    """Limits, equirectangular aspect, graticule ticks with degree labels."""
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect(1.0 / _cos_lat(ext), adjustable="box")
    xt, sx = _ticks(ext[0], ext[1], max_intervals[0])
    yt, sy = _ticks(ext[2], ext[3], max_intervals[1])
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_xticklabels([_fmt_deg(v, sx, "E", "W") for v in xt] if labels[0] else [])
    ax.set_yticklabels([_fmt_deg(v, sy, "N", "S") for v in yt] if labels[1] else [])
    ax.grid(True, color="0.5", linewidth=0.3, linestyle=(0, (2, 2)), alpha=0.6)
    ax.set_axisbelow(True)
    if water:
        ax.set_facecolor(water)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)


def _nice_length(target: float) -> float:
    for v in (1000, 500, 250, 200, 100, 50, 25, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1):
        if v <= target:
            return float(v)
    return 0.1


_SCALE_BAR_BOX = (0.36, 0.11)    # axes-fraction footprint of a scale bar with its labels
_NORTH_BOX = (0.09, 0.17)        # ... of the north arrow

# Places for map furniture: the four corners and the four edge midpoints.
MAP_ANCHORS = ("lower left", "lower right", "upper left", "upper right",
               "lower center", "upper center", "center left", "center right")
MAP_FURNITURE_MAX_GROW = 0.5     # the map frame may grow by up to this share of its span to clear the ROI


def _anchor_origin(anchor: str, size) -> Tuple[float, float]:
    """Lower-left (axes fraction) of a box of `size` at the named anchor, 2% from
    the frame: a corner ('lower left' ...) or an edge midpoint ('lower center',
    'center left' ...)."""
    sx, sy = size
    x = 0.02 if "left" in anchor else (0.98 - sx if "right" in anchor else 0.5 - 0.5 * sx)
    y = 0.02 if "lower" in anchor else (0.98 - sy if "upper" in anchor else 0.5 - 0.5 * sy)
    return x, y


_corner_origin = _anchor_origin      # name used up to v1.3.1


def map_obstacle(tiles, outline_geoms=None):
    """What map furniture must not cover, as one shapely geometry in lon/lat:
    the ROI's cells and its boundary polygons (None if there is nothing)."""
    from shapely.geometry import shape
    from shapely.ops import unary_union
    parts = []
    g = tiles_outline(tiles)
    if g:
        parts.append(shape(g).buffer(0))
    for og in outline_geoms or []:
        try:
            parts.append(shape(og).buffer(0))
        except Exception:
            pass
    return unary_union(parts) if parts else None


def _grow_extent(ext, side: str, g: float):
    """Extent grown on one side by g times its span on that axis, then fitted
    back into MAP_ASPECT_RANGE (growing the other axis about its centre if needed)."""
    lon0, lon1, lat0, lat1 = ext
    if side == "bottom":
        lat0 = max(lat0 - g * (lat1 - lat0), -89.0)
    elif side == "top":
        lat1 = min(lat1 + g * (lat1 - lat0), 89.0)
    elif side == "left":
        lon0 -= g * (lon1 - lon0)
    elif side == "right":
        lon1 += g * (lon1 - lon0)
    new = [lon0, lon1, lat0, lat1]
    lo, hi = MAP_ASPECT_RANGE
    return fit_aspect(new, min(max(panel_aspect(new), lo), hi))


def plan_furniture(ext, obstacle, items, points=(), max_grow: float = MAP_FURNITURE_MAX_GROW,
                   clearance: float = 0.015):
    """Positions for map furniture that hide none of the ROI. items: [(name,
    size in axes fraction, preferred anchors)]. Every combination of
    MAP_ANCHORS with no two boxes overlapping is scored by the area of
    `obstacle` (lon/lat geometry) plus small discs around `points` (lon, lat),
    widened by `clearance` (axes fraction) so furniture keeps a visible gap,
    that the boxes would cover; ties go to the preferred anchors, in item order.
    If every combination hides part of the ROI, the frame is grown on one side
    (bottom, right, top, left; by 10%, 20%, ... up to max_grow of its span) and
    the first clear arrangement is used; failing that, the least-hiding one.
    Returns (extent, {name: anchor}, report) with report = {"grown_side",
    "grown_by", "roi_hidden_pct" (share of the ROI area, with the station
    discs, inside the frame that the furniture boxes hide)}."""
    import itertools
    from shapely.geometry import box as _box, Point
    from shapely.ops import unary_union
    from shapely import affinity

    def evaluate(e):
        sx, sy = 1.0 / (e[1] - e[0]), 1.0 / (e[3] - e[2])
        obs = []
        if obstacle is not None and not obstacle.is_empty:
            obs.append(affinity.affine_transform(obstacle, [sx, 0.0, 0.0, sy, -e[0] * sx, -e[2] * sy]))
        for lon, lat in points:
            fx, fy = (lon - e[0]) * sx, (lat - e[2]) * sy
            if 0.0 <= fx <= 1.0 and 0.0 <= fy <= 1.0:
                obs.append(Point(fx, fy).buffer(0.02))
        geom = unary_union(obs).intersection(_box(0.0, 0.0, 1.0, 1.0)) if obs else None
        total = geom.area if geom is not None else 0.0
        near = geom.buffer(clearance) if geom is not None else None     # keep a small gap, too
        boxes, hidden, hidden_roi = {}, {}, {}
        for name, size, _ in items:
            for a in MAP_ANCHORS:
                x0, y0 = _anchor_origin(a, size)
                boxes[name, a] = (x0, y0, x0 + size[0], y0 + size[1])
                b = _box(*boxes[name, a])
                hidden[name, a] = b.intersection(near).area if near is not None else 0.0
                hidden_roi[name, a] = b.intersection(geom).area if geom is not None else 0.0

        def overlap(p, q):
            return min(p[2], q[2]) - max(p[0], q[0]) > 1e-9 and min(p[3], q[3]) - max(p[1], q[1]) > 1e-9
        best = None
        for combo in itertools.product(MAP_ANCHORS, repeat=len(items)):
            bs = [boxes[it[0], a] for it, a in zip(items, combo)]
            if any(overlap(bs[i], bs[j]) for i in range(len(bs)) for j in range(i + 1, len(bs))):
                continue
            cost = sum(hidden[it[0], a] for it, a in zip(items, combo))
            ranks = tuple(it[2].index(a) if a in it[2] else len(it[2]) + MAP_ANCHORS.index(a)
                          for it, a in zip(items, combo))
            key = (round(cost, 9), ranks)
            if best is None or key < best[0]:
                best = (key, combo)
        roi = sum(hidden_roi[it[0], a] for it, a in zip(items, best[1]))
        return best[0][0], (100.0 * roi / total if total > 0 else 0.0), dict(zip([it[0] for it in items], best[1]))

    tried = []
    for g in [0.0] + [round(0.1 * i, 1) for i in range(1, int(round(max_grow * 10)) + 1)]:
        for side in (("none",) if g == 0.0 else ("bottom", "right", "top", "left")):
            e = ext if g == 0.0 else _grow_extent(ext, side, g)
            cost, pct, anchors = evaluate(e)
            tried.append((round(cost, 9), pct, g, e, anchors, side))
            if cost <= 1e-7:
                return e, anchors, {"grown_side": None if g == 0.0 else side, "grown_by": g, "roi_hidden_pct": 0.0,
                                    "clearance": clearance}
    cost, pct, g, e, anchors, side = min(tried, key=lambda t: (round(t[1], 6), t[0], t[2]))
    logger.warning(f"Map furniture hides {pct:.2f}% of the ROI at best (frame grown by {g:g} on the {side}).")
    return e, anchors, {"grown_side": None if g == 0.0 else side, "grown_by": g, "roi_hidden_pct": pct,
                        "clearance": clearance}


def choose_corners(ext, lonlat: np.ndarray, items) -> Dict[str, str]:
    """Assign each map furniture item to a corner, in the order given, choosing
    the free corner whose box hides the fewest of the points `lonlat` (N x 2)
    -- ties go to the item's preferred order. items: [(name, size, preferred corners)]."""
    lonlat = np.asarray(lonlat, dtype=float).reshape(-1, 2)
    fx = (lonlat[:, 0] - ext[0]) / (ext[1] - ext[0])
    fy = (lonlat[:, 1] - ext[2]) / (ext[3] - ext[2])
    free = ["lower left", "lower right", "upper left", "upper right"]
    out = {}
    for name, size, preferred in items:
        def cost(c):
            x0, y0 = _corner_origin(c, size)
            hidden = int(((fx >= x0) & (fx <= x0 + size[0]) & (fy >= y0) & (fy <= y0 + size[1])).sum())
            return (hidden, preferred.index(c) if c in preferred else len(preferred))
        best = min(free, key=cost)
        out[name] = best
        free.remove(best)
    return out


def add_scale_bar(ax, ext, frac: float = 0.22, corner: str = "lower left", fontsize: float = 6.5,
                  compact: bool = False):
    """Two-segment scale bar, true along the parallel it is drawn on, on a
    white backing box, inside the _SCALE_BAR_BOX footprint at `corner` (any of
    MAP_ANCHORS). compact: a single bar labelled with its length only."""
    bx, by = _anchor_origin(corner, _SCALE_BAR_BOX)
    lat_b = ext[2] + (by + 0.035) * (ext[3] - ext[2])
    kpd = _km_per_deg_lon(lat_b)
    length_km = _nice_length(frac * (ext[1] - ext[0]) * kpd)
    dlon = length_km / kpd
    bar_frac = dlon / (ext[1] - ext[0])
    span = bar_frac + (0.0 if compact else 0.07)              # bar plus its right-hand label
    x0 = (bx + 0.03 if "left" in corner else
          (bx + _SCALE_BAR_BOX[0] - 0.03 - span if "right" in corner else 0.5 - 0.5 * span))
    lon_b = ext[0] + x0 * (ext[1] - ext[0])
    h = 0.012 * (ext[3] - ext[2])
    box_w = dlon + (0.03 if compact else 0.10) * (ext[1] - ext[0])
    backing = ax.add_patch(mpatches.FancyBboxPatch(
        (lon_b - 0.015 * (ext[1] - ext[0]), lat_b - 0.022 * (ext[3] - ext[2])), box_w, 0.085 * (ext[3] - ext[2]),
        boxstyle="square,pad=0", facecolor="white", edgecolor="none", alpha=0.85, zorder=8))
    segments = (("black",),) if compact else (("black",), ("white",))
    for i, (color,) in enumerate(segments):
        w = dlon / len(segments)
        ax.add_patch(mpatches.Rectangle((lon_b + i * w, lat_b), w, h, facecolor=color,
                                        edgecolor="black", linewidth=0.5, zorder=9))
    labels = (((lon_b + dlon / 2, f"{length_km:g} km", "center"),) if compact else
              ((lon_b, "0", "center"), (lon_b + dlon / 2, f"{length_km / 2:g}", "center"),
               (lon_b + dlon, f"{length_km:g} km", "left")))
    for x, txt, ha in labels:
        ax.text(x, lat_b + 1.8 * h, txt, ha=ha, va="bottom", fontsize=fontsize, zorder=10, path_effects=_HALO)
    return length_km, backing


def add_north_arrow(ax, corner: str = "upper right", fontsize: float = 7.5):
    """North arrow on a white backing box (north is up in these lat/lon panels),
    at `corner` (any of MAP_ANCHORS)."""
    x0, y0 = _anchor_origin(corner, _NORTH_BOX)
    x, y, length = x0 + 0.5 * _NORTH_BOX[0], y0 + 0.02, 0.09
    backing = ax.add_patch(mpatches.FancyBboxPatch((x0 + 0.01, y0), _NORTH_BOX[0] - 0.02, _NORTH_BOX[1],
                                                   boxstyle="square,pad=0", transform=ax.transAxes,
                                                   facecolor="white", edgecolor="none", alpha=0.85, zorder=8))
    ax.annotate("", xy=(x, y + length), xytext=(x, y), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", color="black", lw=0.9),
                zorder=10)
    ax.text(x, y + length + 0.008, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold", zorder=10)
    return backing


class LabelPlacer:
    """Greedy collision-avoiding placement of map labels: each label takes the
    first candidate position whose text box stays inside the panel and clear
    of everything placed or registered before it."""

    def __init__(self, ax):
        self.ax = ax
        ax.apply_aspect()
        self.renderer = ax.figure.canvas.get_renderer()
        self.frame = ax.get_window_extent(self.renderer)
        self.boxes = []

    def occupy_point(self, lon: float, lat: float, radius_pt: float = 4.0):
        from matplotlib.transforms import Bbox
        x, y = self.ax.transData.transform((lon, lat))
        r = radius_pt * self.ax.figure.dpi / 72.0
        self.boxes.append(Bbox([[x - r, y - r], [x + r, y + r]]))

    def occupy_artist(self, artist):
        self.boxes.append(artist.get_window_extent(self.renderer))

    def _fits(self, bb) -> bool:
        f = self.frame
        return (f.x0 <= bb.x0 and bb.x1 <= f.x1 and f.y0 <= bb.y0 and bb.y1 <= f.y1
                and not any(bb.overlaps(b) for b in self.boxes))

    def place(self, text: str, candidates, required: bool = True, **kw):
        """candidates: [(lon, lat, dx_pt, dy_pt, ha, va)]. Returns the text, or
        None if required is False and no candidate fits."""
        for lon, lat, dx, dy, ha, va in candidates:
            t = self.ax.annotate(text, (lon, lat), xytext=(dx, dy), textcoords="offset points",
                                 ha=ha, va=va, **kw)
            bb = t.get_window_extent(self.renderer).expanded(1.04, 1.10)
            if self._fits(bb):
                self.boxes.append(bb)
                return t
            t.remove()
        if not required:
            return None
        lon, lat, dx, dy, ha, va = candidates[0]
        t = self.ax.annotate(text, (lon, lat), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, **kw)
        self.boxes.append(t.get_window_extent(self.renderer))
        return t


def _point_candidates(lon, lat, r: float = 5.0):
    """Label positions around a point marker of radius r points: right, left, above, below, diagonals."""
    return [(lon, lat, r, 0, "left", "center"), (lon, lat, -r, 0, "right", "center"),
            (lon, lat, 0, r, "center", "bottom"), (lon, lat, 0, -r, "center", "top"),
            (lon, lat, r, r, "left", "bottom"), (lon, lat, -r, r, "right", "bottom"),
            (lon, lat, r, -r, "left", "top"), (lon, lat, -r, -r, "right", "top")]


def _box_candidates(ext, gap: float = 4.0):
    """Label positions around a rectangle [lon0, lon1, lat0, lat1]."""
    xc, yc = 0.5 * (ext[0] + ext[1]), 0.5 * (ext[2] + ext[3])
    return [(ext[1], yc, gap, 0, "left", "center"), (ext[0], yc, -gap, 0, "right", "center"),
            (xc, ext[3], 0, gap, "center", "bottom"), (xc, ext[2], 0, -gap, "center", "top")]


def load_country_shapes(search_dirs=(), download_dir: Optional[str] = None) -> List[dict]:
    """Natural Earth 1:50m Admin 0 country polygons, [{'name', 'label', 'bbox',
    'rings'}]; [] if no local copy exists and the download fails. A local copy
    in any of search_dirs is used first; otherwise the files are downloaded
    from the Natural Earth repository (tag NATURAL_EARTH_VERSION) into
    download_dir."""
    base = _NATURAL_EARTH_BASE
    path = None
    for d in [x for x in search_dirs if x]:
        if all(os.path.exists(os.path.join(d, f"{base}.{e}")) for e in ("shp", "shx", "dbf")):
            path = os.path.join(d, base + ".shp")
            break
    if path is None and download_dir:
        try:
            import urllib.request
            os.makedirs(download_dir, exist_ok=True)
            for e in ("shp", "shx", "dbf"):
                part = os.path.join(download_dir, f"{base}.{e}.part")
                with urllib.request.urlopen(NATURAL_EARTH_URL.format(ext=e), timeout=60) as r, open(part, "wb") as fh:
                    fh.write(r.read())
                os.replace(part, os.path.join(download_dir, f"{base}.{e}"))
            path = os.path.join(download_dir, base + ".shp")
            logger.info(f"Natural Earth {NATURAL_EARTH_VERSION} borders downloaded to {download_dir}")
        except Exception as exc:
            logger.warning(f"Natural Earth borders unavailable ({exc}); maps are drawn without borders.")
            return []
    if path is None:
        return []
    if path in _country_shape_cache:
        return _country_shape_cache[path]
    try:
        import shapefile   # pyshp
        reader = shapefile.Reader(path)
        fields = [f[0] for f in reader.fields[1:]]
        shapes = []
        for sr in reader.iterShapeRecords():
            pts = np.asarray(sr.shape.points, dtype=float)
            if pts.size == 0:
                continue
            rec = dict(zip(fields, sr.record))
            parts = list(sr.shape.parts) + [len(pts)]
            label = (rec.get("LABEL_X"), rec.get("LABEL_Y"))
            shapes.append({"name": str(rec.get("NAME") or rec.get("ADMIN") or ""),
                           "label": label if None not in label else None,
                           "bbox": list(sr.shape.bbox),
                           "rings": [pts[parts[i]:parts[i + 1]] for i in range(len(parts) - 1)]})
        with open(path, "rb") as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
        _country_shape_cache[path] = shapes
        _country_shape_cache["_meta"] = {"path": path, "sha256": sha, "version": NATURAL_EARTH_VERSION}
        return shapes
    except Exception as exc:
        logger.warning(f"Could not read Natural Earth borders at {path} ({exc}); maps are drawn without borders.")
        return []


def _in_ext(ext, lon, lat, margin: float = 0.0) -> bool:
    mx, my = margin * (ext[1] - ext[0]), margin * (ext[3] - ext[2])
    return ext[0] + mx <= lon <= ext[1] - mx and ext[2] + my <= lat <= ext[3] - my


def draw_countries(ax, shapes, ext) -> List[dict]:
    """Land fill and borders/coastlines inside ext; returns the countries drawn."""
    vis = [s for s in shapes if not (s["bbox"][2] < ext[0] or s["bbox"][0] > ext[1]
                                     or s["bbox"][3] < ext[2] or s["bbox"][1] > ext[3])]
    rings = [r for s in vis for r in s["rings"]]
    if rings:
        ax.add_collection(PolyCollection(rings, facecolor=MAP_COLORS["land"], edgecolor="none", zorder=0.2))
        ax.add_collection(LineCollection(rings, colors=MAP_COLORS["border"], linewidths=0.45, zorder=0.7))
    return vis


def label_countries(placer: "LabelPlacer", countries, ext, fontsize: float = 6.0) -> int:
    """Country names at their Natural Earth label points, for countries with a
    sizeable part of the frame, skipping any name that would collide."""
    frame = (ext[1] - ext[0]) * (ext[3] - ext[2])
    n = 0
    for s in countries:
        if not s["label"] or not _in_ext(ext, *s["label"], margin=0.04):
            continue
        bx0, by0, bx1, by1 = s["bbox"]
        overlap = max(0.0, min(bx1, ext[1]) - max(bx0, ext[0])) * max(0.0, min(by1, ext[3]) - max(by0, ext[2]))
        if overlap < 0.03 * frame:
            continue
        t = placer.place(s["name"].upper(), [(*s["label"], 0, 0, "center", "center")], required=False,
                         fontsize=fontsize, color=MAP_COLORS["country_label"], style="italic", zorder=1.7,
                         path_effects=[patheffects.withStroke(linewidth=2.0, foreground=MAP_COLORS["land"])])
        n += t is not None
    return n


def draw_station_markers(ax, stations, ext, margin: float = 0.0) -> List[Tuple[str, float, float]]:
    """Ground-station triangles inside ext (at least `margin` of the frame from its edges)."""
    drawn = []
    for name, gs in stations.items():
        lon, lat = gs["lng"], gs["lat"]
        if _in_ext(ext, lon, lat, margin):
            ax.plot(lon, lat, marker="^", ms=5.5, mfc=MAP_COLORS["station"], mec="black", mew=0.6,
                    ls="none", zorder=6)
            drawn.append((name, lon, lat))
    return drawn


def label_stations(placer: "LabelPlacer", drawn, fontsize: float = 6.5):
    for name, lon, lat in drawn:
        placer.place(name.replace("_", " "), _point_candidates(lon, lat, 4.5), fontsize=fontsize,
                     zorder=7, path_effects=_HALO)


def map_figure(aspects, fig_w: float = MAP_FIG_WIDTH_IN, left: float = 0.52, gap: float = 0.62,
               right: float = 0.12, top: float = 0.30, bottom: float = 0.80):
    """Figure with equal-height map panels side by side whose boxes match the
    panels' displayed aspects exactly (so set_aspect never shrinks them).
    Returns (fig, axes, layout) with layout sizes in inches."""
    h = (fig_w - left - right - gap * (len(aspects) - 1)) / float(sum(aspects))
    fig_h = top + h + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))
    axes, x = [], left
    for a in aspects:
        axes.append(fig.add_axes([x / fig_w, bottom / fig_h, a * h / fig_w, h / fig_h]))
        x += a * h + gap
    return fig, axes, {"fig_w": fig_w, "fig_h": fig_h, "panel_h": h, "bottom": bottom, "top": top}


def draw_tile_grid(ax, tiles, color=None, edges: Optional[bool] = None, zorder: float = 2.0,
                   rasterize: Optional[bool] = None):
    """ROI cells as one dissolved fill plus (while visible) white cell outlines."""
    tiles = list(tiles)
    path = geo_path([tiles_outline(tiles)])
    if path is not None:
        patch = mpatches.PathPatch(path, facecolor=color or MAP_COLORS["tile"], edgecolor="none", zorder=zorder)
        ax.add_patch(patch)
    if edges if edges is not None else len(tiles) <= MAP_EDGE_TILES:
        lw = 0.9 if len(tiles) <= 100 else (0.5 if len(tiles) <= 800 else 0.25)
        lines = LineCollection([np.vstack([p, p[:1]]) for p in tile_polygons(tiles)],
                               colors=MAP_COLORS["tile_edge"], linewidths=lw, zorder=zorder + 0.1)
        if rasterize if rasterize is not None else len(tiles) > MAP_RASTER_TILES:
            lines.set_rasterized(True)
        ax.add_collection(lines)


def draw_outline(ax, geoms, color=None, lw: float = 0.8, zorder: float = 4.0, ls="-", rasterized=None):
    """ROI outline. rasterized (v1.5.0): default True when the outline is a
    WDPA boundary (ROI_OUTLINE_SOURCE == "WDPA"). The WDPA licence forbids
    redistributing the data, so saved figures carry the boundary as pixels --
    as a PNG always does -- and a PDF holds no vector copy of the polygon."""
    lines = geo_lines(geoms)
    if lines:
        lc = LineCollection(lines, colors=color or MAP_COLORS["boundary"], linewidths=lw, linestyles=ls,
                            zorder=zorder)
        lc.set_gid("roi_outline")
        if globals().get("ROI_OUTLINE_SOURCE") == "WDPA" if rasterized is None else rasterized:
            lc.set_rasterized(True)
        ax.add_collection(lc)


_ZOOM_BOX = (0.34, 0.34)


def zoom_centre(tiles, n_across: int = 9) -> str:
    """The cell deepest inside the ROI (hexagon rings to the nearest cell
    outside, capped at 3x the zoom radius), ties going to the cell nearest the
    centroid of the ROI's cells -- so the zoom shows a full patch of the grid,
    not the ROI's rim, even for ROIs made of separate parts."""
    tset = tiles if isinstance(tiles, (set, frozenset, dict)) else set(tiles)
    cap = 3 * (n_across // 2 + 1)
    depth, frontier = {}, []
    for t in tset:
        if any(n not in tset for n in h3.grid_disk(t, 1)):
            depth[t] = 0
            frontier.append(t)
    d = 0
    while frontier and d < cap:
        d += 1
        nxt = []
        for t in frontier:
            for n in h3.grid_disk(t, 1):
                if n in tset and n not in depth:
                    depth[n] = d
                    nxt.append(n)
        frontier = nxt
    tl = list(tset)
    cent = np.array([h3.cell_to_latlng(t) for t in tl])
    dist2 = ((cent - cent.mean(axis=0)) ** 2).sum(axis=1)
    deep = np.array([min(depth.get(t, cap), cap) for t in tl])
    return tl[int(np.lexsort((dist2, -deep))[0])]


def add_zoom_inset(ax, tiles, ext, corner: str, outline_geoms=None, n_across: int = 9, color=None,
                   colors_by_tile: Optional[Dict[str, Any]] = None):
    """Magnified view (n_across cells wide) of the cells around zoom_centre(),
    in the given corner of `ax`, with a zoom indicator and its own scale bar.
    colors_by_tile: optional {tile: facecolor} (Cell 7's mosaic map)."""
    centre = zoom_centre(tiles, n_across)
    lat_c, lon_c = h3.cell_to_latlng(centre)
    edge_km = h3.average_hexagon_edge_length(h3.get_resolution(centre), unit="km")
    half_km = 0.5 * n_across * edge_km * math.sqrt(3.0)
    dlon, dlat = half_km / _km_per_deg_lon(lat_c), half_km / _KM_PER_DEG_LAT
    zext = fit_aspect([lon_c - dlon, lon_c + dlon, lat_c - dlat, lat_c + dlat], 1.0)
    axins = ax.inset_axes([*_anchor_origin(corner, _ZOOM_BOX), *_ZOOM_BOX])
    tset = tiles if isinstance(tiles, (set, frozenset, dict)) else set(tiles)
    near = [t for t in h3.grid_disk(centre, n_across) if t in tset]
    if colors_by_tile is None:
        draw_tile_grid(axins, near, color=color, edges=True, rasterize=False)
    else:
        polys = tile_polygons(near)
        axins.add_collection(PolyCollection(polys, facecolors=[colors_by_tile[t] for t in near],
                                            edgecolors="white", linewidths=0.6, zorder=2))
    draw_outline(axins, outline_geoms, lw=0.8)
    axins.set_xlim(zext[0], zext[1])
    axins.set_ylim(zext[2], zext[3])
    axins.set_aspect(1.0 / _cos_lat(zext), adjustable="box")
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_facecolor("white")
    for spine in axins.spines.values():
        spine.set_linewidth(0.6)
    add_scale_bar(axins, zext, frac=0.35, corner="lower left", fontsize=5.5, compact=True)
    ax.indicate_inset_zoom(axins, edgecolor="black", linewidth=0.6, alpha=1.0)
    return axins


def _wdpa_citation(reports) -> Tuple[Optional[str], Optional[str]]:
    """(full citation, short credit) for the WDPA release the boundaries came from,
    in the form the WDPA/WDOECM Manual 1.6 (section 5.3) asks for."""
    if not reports:
        return None, None
    rel = max(_wdpa_release(r.get("file", "")) for r in reports) if '_wdpa_release' in globals() else (0, 0)
    if rel == (0, 0):
        when, short = "[month/year of the version downloaded]", "WDPA"
    else:
        months = ("January", "February", "March", "April", "May", "June", "July", "August",
                  "September", "October", "November", "December")
        when = f"{months[rel[1] - 1]} {rel[0]}"
        short = f"WDPA, {months[rel[1] - 1][:3]} {rel[0]}"
    year = rel[0] if rel[0] else "[year]"
    return (f"UNEP-WCMC and IUCN ({year}), Protected Planet: The World Database on Protected Areas (WDPA) "
            f"[On-line], {when}, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.",
            short)


def roi_outline_for_maps():
    """(geometries, legend label, source) of the ROI outline: the official
    boundary kept in memory by Cell 1, the ROI polygon, or -- if Cell 1 predates
    ROI_OUTLINE_GEOMS -- the dissolved outline of the ROI's cells."""
    geoms = globals().get("ROI_OUTLINE_GEOMS")
    if geoms:
        src = globals().get("ROI_OUTLINE_SOURCE", "ROI definition")
        return list(geoms), ("WDPA boundary" if src == "WDPA" else "ROI boundary"), src
    if globals().get("roi_coords"):
        return ([{"type": "Polygon", "coordinates": [[(lon, lat) for lat, lon in roi_coords]]}],
                "ROI boundary", "ROI definition")
    return [tiles_outline(target_tiles)], "ROI outline (cells)", "H3 cells"


def roi_short_name() -> str:
    spec = globals().get("ROI_DEFINITIONS", {}).get(ACTIVE_ROI, {})
    return spec.get("short_name") or ACTIVE_ROI.replace("_", " ").title()


# =============================================================================
# GROUND TRACK COMPUTATION
# =============================================================================

def compute_track(sat: Dict, epoch: Any, window_min: int, step_min: int) -> Tuple[List[float], List[float]]:
    """Propagates satellite state vectors and transforms them into geodetic coordinates."""
    lats, lngs = [], []
    prop = sat["propagator"]
    for m in range(0, window_min + 1, step_min):
        try:
            t = epoch.shiftedBy(float(m * 60))
            pv = prop.getPVCoordinates(t, itrf)
            geo = earth.transform(pv.getPosition(), itrf, t)
            lats.append(FastMath.toDegrees(geo.getLatitude()))
            lngs.append(FastMath.toDegrees(geo.getLongitude()))
        except Exception as e:
            logger.warning(f"Propagation anomaly at t={m}m for {sat['id']}: {e}")
            lats.append(float('nan'))
            lngs.append(float('nan'))
    return lats, lngs


def _break_at_antimeridian(lats: List[float], lngs: List[float]) -> Tuple[List[float], List[float]]:
    """Inserts NaN to lift the plotter pen when crossing the International Date Line."""
    if not lngs:
        return lats, lngs
    out_lat, out_lng = [lats[0]], [lngs[0]]
    for i in range(1, len(lngs)):
        if abs(lngs[i] - lngs[i - 1]) > 180.0:
            out_lat.append(np.nan)
            out_lng.append(np.nan)
        out_lat.append(lats[i])
        out_lng.append(lngs[i])
    return out_lat, out_lng


def frame_crossings(sats, epoch, ext, hours: float = MAP_TRACK_HOURS, max_sats: int = MAP_TRACK_MAX_SATS,
                    max_segments: int = MAP_TRACK_MAX_SEGMENTS):
    """Ground-track segments crossing ext during the first `hours` hours, for an
    evenly spaced, deterministic sample of up to max_sats satellites (1-minute
    steps). Returns the earliest max_segments as dicts (sat, t_min, lon, lat,
    inside = index range within the frame), and the sample size."""
    n = len(sats)
    if n == 0:
        return [], 0
    idx = sorted(set(np.linspace(0, n - 1, min(max_sats, n)).round().astype(int).tolist()))
    segs = []
    for i in idx:
        lats, lons = compute_track(sats[i], epoch, int(round(hours * 60)), 1)
        lats, lons = np.asarray(lats, dtype=float), np.asarray(lons, dtype=float)
        inside = ((lons >= ext[0]) & (lons <= ext[1]) & (lats >= ext[2]) & (lats <= ext[3]))
        j = 0
        while j < inside.size:
            if not inside[j]:
                j += 1
                continue
            k = j
            while k + 1 < inside.size and inside[k + 1]:
                k += 1
            a, b = max(j - 1, 0), min(k + 1, inside.size - 1)      # one step past each edge
            seg_lon, seg_lat = lons[a:b + 1], lats[a:b + 1]
            if b > a and np.all(np.isfinite(seg_lon)) and np.all(np.abs(np.diff(seg_lon)) < 180.0):
                segs.append({"sat": sats[i]["id"], "t_min": int(j), "lon": seg_lon, "lat": seg_lat,
                             "inside": (j - a, k - a)})
            j = k + 1
    segs.sort(key=lambda s: (s["t_min"], s["sat"]))
    return segs[:max_segments], len(idx)


def draw_track(ax, seg, color=None, lw: float = 0.7):
    color = color or MAP_COLORS["track"]
    ax.plot(seg["lon"], seg["lat"], color=color, lw=lw, alpha=0.85, zorder=1.5, solid_capstyle="butt")
    i0, i1 = seg["inside"]
    m = max(i0, min((i0 + i1) // 2, len(seg["lon"]) - 2))
    ax.annotate("", xy=(seg["lon"][m + 1], seg["lat"][m + 1]), xytext=(seg["lon"][m], seg["lat"][m]),
                arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.22", color=color, lw=lw),
                zorder=1.6)


# =============================================================================
# FIGURE GENERATION
# =============================================================================

logger.info("Constructing study-area map...")
_roi_name = roi_short_name()
_outline_geoms, _outline_label, _outline_source = roi_outline_for_maps()
_n_tiles_map = len(target_tiles)
_roi_area_km2 = float(sum(h3.cell_area(t, unit="km^2") for t in target_tiles))
_wdpa_cite, _wdpa_short = _wdpa_citation(globals().get("roi_boundary_report"))
_tile_lonlat = np.array([h3.cell_to_latlng(t)[::-1] for t in target_tiles])        # (lon, lat)
_roi_centroid = tuple(_tile_lonlat.mean(axis=0))                                   # (lon, lat)

_detail_ext = roi_map_extent(target_tiles, _outline_geoms)
# Furniture of panel (b) goes where it hides none of the ROI (cells, boundary,
# stations); the frame grows on one side if no such place exists.
_furniture = [("scale", _SCALE_BAR_BOX, ["lower left", "lower right", "upper left", "upper right"]),
              ("north", _NORTH_BOX, ["upper right", "upper left", "lower right", "lower left"])]
if _n_tiles_map > MAP_ZOOM_TILES:
    _furniture.insert(0, ("zoom", _ZOOM_BOX, ["lower right", "upper left", "upper right", "lower left"]))
_detail_ext, _corners, _furniture_report = plan_furniture(
    _detail_ext, map_obstacle(target_tiles, _outline_geoms), _furniture,
    points=[(gs["lng"], gs["lat"]) for gs in GROUND_STATIONS.values()])
_map_aspect = panel_aspect(_detail_ext)
_regional_ext = regional_extent(_detail_ext, GROUND_STATIONS, _map_aspect)

_ne_dirs = [os.path.join(ARF_BOUNDARY_DIR, "natural_earth") if 'ARF_BOUNDARY_DIR' in globals() else None,
            os.path.join(ARF_CACHE_DIR, "natural_earth") if 'ARF_CACHE_DIR' in globals() else None]
_ne_download = next((d for d in _ne_dirs if d and os.path.isdir(os.path.dirname(d))), None)
_country_shapes = load_country_shapes(_ne_dirs, _ne_download)

with plt.rc_context(_MAP_RC):
    fig, (ax_a, ax_b), _layout = map_figure([_map_aspect, _map_aspect])

    # --- (a) Regional setting ------------------------------------------------
    setup_map_axes(ax_a, _regional_ext, max_intervals=(5, 5),
                   water=MAP_COLORS["water"] if _country_shapes else None)
    _countries = draw_countries(ax_a, _country_shapes, _regional_ext)
    try:
        _segments, _n_track_sample = frame_crossings(constellation, sim_epoch, _regional_ext)
    except Exception as exc:          # a map without tracks beats no map
        logger.warning(f"Ground tracks not drawn ({exc}).")
        _segments, _n_track_sample = [], 0
    for seg in _segments:
        draw_track(ax_a, seg)
    _roi_patch_path = geo_path([tiles_outline(target_tiles)])
    if _roi_patch_path is not None:
        ax_a.add_patch(mpatches.PathPatch(_roi_patch_path, facecolor=MAP_COLORS["tile"],
                                          edgecolor=MAP_COLORS["tile"], linewidth=0.4, zorder=3))
    draw_outline(ax_a, _outline_geoms, lw=0.5, zorder=3.5)
    _small_roi = (_detail_ext[1] - _detail_ext[0]) < 0.05 * (_regional_ext[1] - _regional_ext[0])
    if _small_roi:
        _roi_marker = ax_a.plot(_roi_centroid[0], _roi_centroid[1], marker="o", ms=9, mfc="none", mec="black",
                                mew=0.9, zorder=5)[0]
    else:
        _roi_marker = ax_a.add_patch(mpatches.Rectangle(
            (_detail_ext[0], _detail_ext[2]), _detail_ext[1] - _detail_ext[0], _detail_ext[3] - _detail_ext[2],
            fill=False, edgecolor=MAP_COLORS["extent"], linewidth=0.7, zorder=5))
    _stations_a = draw_station_markers(ax_a, GROUND_STATIONS, _regional_ext)
    _, _bar_a = add_scale_bar(ax_a, _regional_ext, corner="lower left")
    _placer_a = LabelPlacer(ax_a)
    _placer_a.occupy_artist(_roi_marker)
    _placer_a.occupy_artist(_bar_a)
    for _, lon, lat in _stations_a:
        _placer_a.occupy_point(lon, lat)
    label_stations(_placer_a, _stations_a)
    _placer_a.place(_roi_name, _point_candidates(*_roi_centroid, 7.0) if _small_roi else _box_candidates(_detail_ext),
                    fontsize=6.5, fontweight="bold", zorder=7, path_effects=_HALO)
    _n_country_labels = label_countries(_placer_a, _countries, _regional_ext)
    ax_a.set_title("(a) Regional setting", loc="left", pad=4)

    # --- (b) ROI and its H3 tiling --------------------------------------------
    setup_map_axes(ax_b, _detail_ext)
    ax_b.set_facecolor("white")
    draw_tile_grid(ax_b, target_tiles)
    draw_outline(ax_b, _outline_geoms, lw=0.8)
    _, _bar_b = add_scale_bar(ax_b, _detail_ext, corner=_corners["scale"])
    _north_b = add_north_arrow(ax_b, corner=_corners["north"])
    _zoom_corner = None
    if "zoom" in _corners:
        add_zoom_inset(ax_b, target_tiles, _detail_ext, _corners["zoom"], _outline_geoms)
        _zoom_corner = _corners["zoom"]
    _stations_b = draw_station_markers(ax_b, GROUND_STATIONS, _detail_ext, margin=0.03)
    _placer_b = LabelPlacer(ax_b)
    for artist in (_bar_b, _north_b):
        _placer_b.occupy_artist(artist)
    if _zoom_corner:
        x0, y0 = _corner_origin(_zoom_corner, _ZOOM_BOX)
        _placer_b.boxes.append(ax_b.transAxes.transform_bbox(
            matplotlib.transforms.Bbox.from_bounds(x0, y0, *_ZOOM_BOX)))
    for _, lon, lat in _stations_b:
        _placer_b.occupy_point(lon, lat)
    label_stations(_placer_b, _stations_b)
    ax_b.set_title(f"(b) {_roi_name}: H3 tiling", loc="left", pad=4)

    # --- Legend and credits ----------------------------------------------------
    _res = H3_RESOLUTION if 'H3_RESOLUTION' in globals() else h3.get_resolution(next(iter(target_tiles)))
    _handles = [mpatches.Patch(facecolor=MAP_COLORS["tile"], edgecolor="white",
                               label=f"H3 res-{_res} cells (n = {_n_tiles_map:,}; {_roi_area_km2:,.0f} km²)"),
                mlines.Line2D([], [], color=MAP_COLORS["boundary"], lw=0.8, label=_outline_label),
                mlines.Line2D([], [], color=MAP_COLORS["station"], marker="^", mec="black", mew=0.6, ms=5.5,
                              ls="none", label="Ground station (notional)")]
    if _segments:
        _handles.append(mlines.Line2D([], [], color=MAP_COLORS["track"], lw=0.7,
                                      label=f"Ground track, first {MAP_TRACK_HOURS} h (arrow: direction)"))
    _handles.append(mlines.Line2D([], [], color="black", marker="o", ms=7, mfc="none", ls="none", mew=0.9,
                                  label="Location of (b)") if _small_roi else
                    mpatches.Patch(facecolor="none", edgecolor=MAP_COLORS["extent"], lw=0.7, label="Extent of (b)"))
    if _countries:
        _handles.append(mlines.Line2D([], [], color=MAP_COLORS["border"], lw=0.6, label="Border / coastline"))
    _fw, _fh, _bot = _layout["fig_w"], _layout["fig_h"], _layout["bottom"]
    fig.legend(handles=_handles, loc="upper center", bbox_to_anchor=(0.5, (_bot - 0.24) / _fh),
               ncol=3, frameon=False, handlelength=1.8, columnspacing=1.4, handletextpad=0.5)
    _credits = [f"Boundary: {_wdpa_short} (UNEP-WCMC & IUCN)" if _wdpa_short else None,
                "Borders: Natural Earth" if _countries else None,
                "Coordinates: WGS 84"]
    fig.text(0.52 / _fw, 0.04 / _fh, " · ".join(c for c in _credits if c), fontsize=6, color="0.4",
             ha="left", va="bottom")

# =============================================================================
# EXPORT & SUMMARY
# =============================================================================

_png_path = os.path.join(_out_dir, "arf_coverage_map.png")
_pdf_path = os.path.join(_out_dir, "arf_coverage_map.pdf")
_meta_path = os.path.join(_out_dir, "arf_cell3_metadata.json")

try:
    fig.savefig(_png_path, dpi=300, bbox_inches='tight')
    fig.savefig(_pdf_path, dpi=300, bbox_inches='tight')
    logger.info(f"Cartographic outputs exported to {_out_dir}")
except Exception as e:
    logger.warning(f"Failed to save figures: {e}")

plt.show()
plt.close(fig)

_fleet_label = (f"proxy constellation (N={ACTIVE_N})" if globals().get("FLEET", "starlink") == "starlink"
                else f"{FLEET} fleet ({len(constellation)} satellites)")
# v1.5.2: the stations are notional receiving sites at approximate coordinates
# (Cell 1), not the operators' ground networks; the caption says so.
_caption = (
    f"Study area. (a) Regional setting of {_roi_name}, with the notional ground stations used for downlink in the "
    f"model ({', '.join(n.replace('_', ' ') for n in GROUND_STATIONS)}; approximate coordinates, not the "
    f"operators' ground networks)"
    + (f" and the first {len(_segments)} crossings of the frame by ground tracks of a sample of "
       f"{_n_track_sample} satellites of the {_fleet_label} during the first {MAP_TRACK_HOURS} h "
       f"from {sim_epoch} (arrows: direction of flight)" if _segments else "")
    + (". Borders and coastlines: Natural Earth 1:50m." if _countries else ".")
    + f" (b) The ROI discretised into {_n_tiles_map:,} H3 resolution-{_res} cells "
      f"({_roi_area_km2:,.0f} km²; mean cell area {_roi_area_km2 / _n_tiles_map:.2f} km²)"
    + (f" and its boundary from the World Database on Protected Areas" if _outline_source == "WDPA" else "")
    + (f"; the inset magnifies the cells near the ROI centre" if _zoom_corner else "")
    + ". Geographic coordinates (WGS 84), drawn with the equirectangular aspect of each panel's centre latitude."
)

_metadata = {
    "figure_type": "study_area_map",
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "roi": ACTIVE_ROI,
    "roi_label": _roi_name,
    "h3_resolution": _res,
    "n_tiles": _n_tiles_map,
    "roi_area_km2": _roi_area_km2,
    "n_ground_stations": len(GROUND_STATIONS),
    "constellation_size": len(constellation),
    "sim_epoch": str(sim_epoch),
    "panel_b_extent_lon_lat": [round(v, 5) for v in _detail_ext],
    "panel_a_extent_lon_lat": [round(v, 4) for v in _regional_ext],
    "outline_source": _outline_source,
    "zoom_inset": _zoom_corner,
    "furniture": dict(_furniture_report, anchors=_corners),
    "ground_tracks": {"hours": MAP_TRACK_HOURS, "satellites_sampled": _n_track_sample,
                      "segments_shown": [{"sat": s["sat"], "t_min": s["t_min"]} for s in _segments]},
    "borders": dict(_country_shape_cache.get("_meta", {}), countries_drawn=len(_countries)) if _countries else None,
    "wdpa_citation": _wdpa_cite,
    "draft_caption": _caption,
}

try:
    with open(_meta_path, 'w') as f:
        json.dump(_metadata, f, indent=2, default=str)
    logger.info(f"Metadata exported to {_meta_path}")
except Exception as e:
    logger.warning(f"Could not save metadata: {e}")

print("\n" + "=" * 60)
print("  CELL 3: Spatial Validation & Visualisation Complete")
print("=" * 60)
print(f"  ROI:                       {_roi_name} ({_n_tiles_map:,} tiles, {_roi_area_km2:,.0f} km²)")
print(f"  Panel (b) extent:          {_detail_ext[1] - _detail_ext[0]:.3f}° x {_detail_ext[3] - _detail_ext[2]:.3f}°")
print(f"  Ground-track crossings:    {len(_segments)} (sample of {_n_track_sample} satellites)")
print(f"  Borders:                   {'Natural Earth ' + NATURAL_EARTH_VERSION if _countries else 'not available'}")
print(f"  Exported PNG:              arf_coverage_map.png (300 dpi)")
print(f"  Exported PDF:              arf_coverage_map.pdf")
print("=" * 60)
print("\n=== CELL 3 COMPLETE ===")

