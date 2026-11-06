# =============================================================================
# CELL 1: ENVIRONMENT SETUP, OREKIT INITIALISATION, & SPATIAL DISCRETISATION
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Establishes the astrodynamics engine (Orekit) and discretises the continuous
#   Region of Interest (ROI) into a Discrete Global Grid System (DGGS) using
#   the H3 hexagonal index at resolution 7 (average cell 5.161 km2, average edge
#   1.406 km; H3 resolution table, https://h3geo.org/docs/core-library/restable/
#   -- the 1.22 km figure of the H3 v3 table is a different edge definition).
#
# Research Proposal Mapping:
#   - Section 3.2.1: Geodesic Grid Model (The Target Manifold R)
#   - Section 3.2.3: Ground Station Network (Telemetry Anchors)
#   - Gap 2 & Justification 2.4.1: Hexagonal uniform adjacency over Lat/Lon
#
# Dependencies:
#   - Cell 0 must be executed (provides ACTIVE_ROI, H3_RESOLUTION, OUT_DIR)
#
# Outputs (Globals for Downstream Pipeline):
#   - itrf: Orekit ITRF reference frame (IERS 2010 conventions)
#   - earth: Orekit WGS84 ellipsoid model
#   - target_tiles: Set[str] of H3 cell indices covering the ROI
#   - ground_station_points: Dict[str, GeodeticPoint]
#   - GROUND_STATIONS: Dict[str, Dict] of ground station metadata
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

"""
ARF Environment Setup & Spatial Discretisation
===============================================

This module executes three foundational initialisation tasks:

1. **JVM & Astrodynamics Engine Initialisation**: Boots the Java Virtual Machine
   (JVM) and loads high-precision physical models (IERS Earth Orientation Parameters,
   leap seconds, and gravitational harmonics).

2. **Geodetic Reference Frame Setup**: Instantiates the International Terrestrial
   Reference Frame (ITRF 2010) and WGS84 ellipsoid for SGP4 orbital propagation.

3. **Discrete Global Grid System (DGGS) Discretisation**: Projects the continuous
   target manifold into H3 hexagonal cells (near-uniform area; H3 cells are not
   exactly equal-area), enabling O(1) set-theoretic spatial operations. ROIs are
   official WDPA protected-area boundaries (the Botswana, Zambia and South
   Africa study areas), or, for the legacy runs, a bounding box (Botswana) or a
   buffered polyline (Nepal disaster corridor).
"""

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 1 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["1"] = _ARF_CELL_VERSION

import os
import sys
import platform
import logging
import json
from typing import Dict, List, Tuple, Set, Optional, Any
from dataclasses import dataclass

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logger = logging.getLogger("arf.cell1")
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

class EnvironmentSetupError(Exception):
    """Raised when JVM, Orekit physical data, or dependencies fail to initialise."""
    pass

class SpatialDiscretisationError(Exception):
    """Raised when the H3 polygon discretisation fails or yields empty tile sets."""
    pass

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()

# =============================================================================
# 1. DEPENDENCY INSTALLATION & JVM INITIALISATION
# =============================================================================

logger.info("Verifying and installing astrodynamics & geospatial packages...")

# v1.5.0: a batch driver that runs many configurations in one session sets
# ARF_SKIP_INSTALL = True after the first run (pip, and the orekit-data clone
# from git, would otherwise repeat for every run).
import subprocess
# v1.5.2: the packages whose behaviour the results depend on are pinned to the
# versions the release was validated with (Colab, 12 October 2026, and the
# test environment): orekit-jpype 13.1.9.0 (with jdk4py 25.0.2.1), h3 4.5.0,
# pyshp 3.1.6. numpy, pandas, shapely, scipy and matplotlib are Colab's own
# (the refined events reproduced to 1e-13 deg across numpy 2.1 and 2.5); every
# version actually used is recorded in the manifest (Cell 9). orekit-data
# (Earth orientation and leap seconds) is installed from its repository head on
# purpose -- a pinned copy would lack the Earth orientation data of new windows
# -- and its version is recorded.
ARF_VALIDATED_VERSIONS = {"orekit-jpype": "13.1.9.0", "jdk4py": "25.0.2.1", "h3": "4.5.0", "pyshp": "3.1.6"}
if globals().get("ARF_SKIP_INSTALL"):
    logger.info("ARF_SKIP_INSTALL: package installation skipped (installed earlier in this session).")
else:
    try:
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q",
             f"orekit-jpype[jdk4py]=={ARF_VALIDATED_VERSIONS['orekit-jpype']}",
             f"jdk4py=={ARF_VALIDATED_VERSIONS['jdk4py']}", f"h3=={ARF_VALIDATED_VERSIONS['h3']}",
             f"pyshp=={ARF_VALIDATED_VERSIONS['pyshp']}", "numpy", "pandas", "shapely"],
            check=True
        )
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q",
             "git+https://gitlab.orekit.org/orekit/orekit-data.git"],
            check=True
        )
        logger.info("Library installation complete.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Package installation failed: {e}")
        raise EnvironmentSetupError("Failed to install required dependencies.") from e

import jdk4py
import orekit_jpype

# Auto-detect JVM Path
try:
    _jvm_path = os.path.join(jdk4py.JAVA_HOME, "lib", "server", "libjvm.so") if platform.system() == "Linux" else None
    if _jvm_path and os.path.exists(_jvm_path):
        orekit_jpype.initVM(jvmpath=_jvm_path)
    else:
        orekit_jpype.initVM()
    logger.info("Orekit JVM initialised successfully.")
except Exception as e:
    logger.info(f"JVM already running or init handled: {e}")

# =============================================================================
# 2. OREKIT PHYSICAL DATA CONTEXT
# =============================================================================

from orekit_jpype.pyhelpers import setup_orekit_curdir

_orekit_data_loaded = False
try:
    setup_orekit_curdir(from_pip_library=True)
    _orekit_data_loaded = True
    logger.info("Orekit physical data context loaded successfully (pip library).")
except Exception:
    try:
        import orekitdata
        from orekit_jpype.pyhelpers import setup_orekit_data
        setup_orekit_data(filenames=orekitdata.__path__[0], from_pip_library=False)
        _orekit_data_loaded = True
        logger.info("Orekit physical data loaded via fallback package path.")
    except Exception as e:
        logger.error(f"Failed to load physical data context: {e}")

if not _orekit_data_loaded:
    raise EnvironmentSetupError("Could not load Orekit physical data (EOP/Leap Seconds).")

# =============================================================================
# 3. SPATIAL & ASTRODYNAMICS IMPORTS
# =============================================================================

import h3
import numpy as np
import pandas as pd

if not h3.__version__.startswith("4."):
    raise EnvironmentSetupError(f"Detected H3 version {h3.__version__}. ARF requires H3 v4.x.")
logger.info(f"H3 DGGS library version confirmed: {h3.__version__}")
# v1.5.2: report any package that differs from the validated version (a warning,
# not an error: the run is still recorded with the versions it used).
from importlib import metadata as _md1
ARF_PACKAGE_VERSIONS = {}
for _pkg, _want in ARF_VALIDATED_VERSIONS.items():
    try:
        ARF_PACKAGE_VERSIONS[_pkg] = _md1.version(_pkg)
    except _md1.PackageNotFoundError:
        ARF_PACKAGE_VERSIONS[_pkg] = None
    if ARF_PACKAGE_VERSIONS[_pkg] != _want:
        logger.warning(f"{_pkg} {ARF_PACKAGE_VERSIONS[_pkg]} differs from the validated {_want}; results may "
                       f"differ from the release's reference values (arf_colab_parity.py checks them).")

from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions, Constants
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.hipparchus.util import FastMath

# Establish ITRF Reference Frame and WGS84 Oblate Spheroid
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    itrf
)
logger.info("Astrodynamics reference frame (ITRF 2010) & WGS84 ellipsoid initialised.")

# =============================================================================
# 4. REGION OF INTEREST (ROI) GEOMETRY DEFINITION
# =============================================================================

def _buffer_waypoints_to_polygon(waypoints_latlon: List[Tuple[float, float]], buffer_km: float) -> List[Tuple[float, float]]:
    """
    Buffer a lat/lon waypoint polyline by buffer_km and return the resulting polygon.
    Uses a flat-earth degree/km approximation local to the waypoints' average latitude.
    """
    import math
    from shapely.geometry import LineString

    lat_avg = sum(w[0] for w in waypoints_latlon) / len(waypoints_latlon)
    km_per_deg_lat = 111.0
    km_per_deg_lon = 111.0 * math.cos(math.radians(lat_avg))
    buffer_deg = buffer_km / min(km_per_deg_lat, km_per_deg_lon)

    line_lonlat = LineString([(lon, lat) for lat, lon in waypoints_latlon])
    corridor_lonlat = line_lonlat.buffer(buffer_deg)
    return [(lat, lon) for lon, lat in corridor_lonlat.exterior.coords]

ROI_DEFINITIONS = {
    "botswana_full": {
        "short_name": "Botswana (bounding box)",  # label on the maps (Cells 3 and 7)
        "description": "Botswana full national boundary (bounding rectangle)",
        "roi_coords": [
            (-17.6618156877, 19.8954577979), (-17.6618156877, 29.4321883481),
            (-26.8285429827, 29.4321883481), (-26.8285429827, 19.8954577979)
        ],
        "ground_stations": {
            "BIUST": {"lat": -22.59, "lng": 27.14, "alt_m": 994.0},
            "Kasane": {"lat": -17.80, "lng": 25.15, "alt_m": 930.0},
            "Kebonye": {"lat": -24.84, "lng": 25.84, "alt_m": 994.0},
        },
        "tle_source": "celestrak_live",
    },
    "nepal_trishuli_corridor": {
        "short_name": "Trishuli corridor",  # label on the maps (Cells 3 and 7)
        "description": (
            "2026 Nepal Trishuli/Lende Khola debris flow and flash flood corridor. "
            "Buffered 7.5km either side of the waypoint polyline."
        ),
        "waypoints_latlon": [
            (28.27778, 85.37778),      # Rasuwagadhi
            (28.1635, 85.3363),        # Syaphru / Syabrubesi
            (27.9833, 85.1833),        # Betrawati
            (27.8881907, 85.1340051),  # Devighat
        ],
        "buffer_km": 7.5,
        "historical_epoch": "2026-08-26T02:52:00",
        "tle_source": "spacetrack_history",
        "ground_stations": {
            "NAST_Khumaltar": {"lat": 27.65479, "lng": 85.32882, "alt_m": 1335.0},
        },
    },
    "ckgr_okavango": {
        "short_name": "CKGR + Okavango Delta",  # label on the maps (Cells 3 and 7)
        "description": (
            "Central Kalahari Game Reserve and the Okavango Delta World Heritage Site, "
            "mosaicked together from their official WDPA boundaries."
        ),
        # Expected GIS areas are the WDPA's own polygon areas (GIS_AREA), used
        # to flag a revised boundary in a later WDPA release.
        "wdpa_sites": {
            7510: {"name": "Central Kalahari Game Reserve",
                   "expected_gis_area_km2": 52225.35, "designation_contains": "Game Reserve"},
            555577555: {"name": "Okavango Delta (World Heritage Site)",
                        "expected_gis_area_km2": 20551.25, "designation_contains": "World Heritage"},
        },
        "ground_stations": {
            "BIUST": {"lat": -22.59, "lng": 27.14, "alt_m": 994.0},
            "Kasane": {"lat": -17.80, "lng": 25.15, "alt_m": 930.0},
            "Kebonye": {"lat": -24.84, "lng": 25.84, "alt_m": 994.0},
        },
        "tle_source": "celestrak_live",
    },
    # Site IDs, names, countries and reported areas below are as shown on each
    # site's Protected Planet page (source_url), accessed 2026-10-08. The pages do
    # not show the designation, so designation_contains is None: the designation
    # is read from the file and recorded, and the site's identity is instead
    # checked by its reported area (REP_AREA) against the page value.
    # Ground stations: the Botswana network (BIUST, Kasane, Kebonye), so downlink
    # waits reflect a Botswana-operated ground segment. ASSUMPTION: these are
    # notional receiving sites (approximate coordinates), not the operators'
    # ground networks; Sentinel-2, Landsat and PlanetScope downlink to their own
    # stations. All of them see every pass over these regions, so every wait is
    # 0 and "delivery" times are acquisition times (Cells 6, 8 say so).
    # Replaces Middle Zambezi (WDPA 555547578), whose WDPA record is a point
    # with no boundary polygon (WDPA/WDOECM Manual 1.6, p.15), so it cannot be tiled.
    "lower_zambezi": {
        "short_name": "Lower Zambezi NP",  # label on the maps (Cells 3 and 7)
        "description": "Lower Zambezi protected area (Zambia), from its official WDPA boundary.",
        "wdpa_sites": {
            7962: {"name": "Lower Zambezi", "designation_contains": None,
                   "expected_rep_area_km2": 4143.0,
                   "source_url": "https://www.protectedplanet.net/7962", "accessed": "2026-10-08"},
        },
        "ground_stations": {
            "BIUST": {"lat": -22.59, "lng": 27.14, "alt_m": 994.0},
            "Kasane": {"lat": -17.80, "lng": 25.15, "alt_m": 930.0},
            "Kebonye": {"lat": -24.84, "lng": 25.84, "alt_m": 994.0},
        },
        "tle_source": "celestrak_live",
    },
    "karoo_nature_reserve": {
        "short_name": "Karoo NR",  # label on the maps (Cells 3 and 7)
        "description": ("Karoo Nature Reserve (South Africa), from its official WDPA boundary. "
                        "32 H3 res-7 tiles with the October 2026 WDPA release (167 km2 tiled for a "
                        "143 km2 site: small sites are quantised to whole tiles)."),
        "wdpa_sites": {
            555563860: {"name": "Karoo Nature Reserve", "designation_contains": None,
                        "expected_rep_area_km2": 143.15,
                        "source_url": "https://www.protectedplanet.net/555563860", "accessed": "2026-10-08"},
        },
        "ground_stations": {
            "BIUST": {"lat": -22.59, "lng": 27.14, "alt_m": 994.0},
            "Kasane": {"lat": -17.80, "lng": 25.15, "alt_m": 930.0},
            "Kebonye": {"lat": -24.84, "lng": 25.84, "alt_m": 994.0},
        },
        "tle_source": "celestrak_live",
    },
}

# -----------------------------------------------------------------------------
# WDPA boundary loading (ROIs defined by "wdpa_sites")
# -----------------------------------------------------------------------------
# Official protected-area boundaries come from the World Database on Protected
# Areas (UNEP-WCMC and IUCN, Protected Planet). The WDPA licence forbids
# redistributing the data, so the boundary files are NOT shipped with this
# framework: download each site from https://www.protectedplanet.net/<SITE_ID>
# as a shapefile ("SHP") and put the download in ARF_BOUNDARY_DIR -- as the zip,
# or extracted, at any folder depth. Files are matched by WDPA site ID, split
# parts (_shp_0, _shp_1, ...) are combined, and the newest release is used. Cite as: UNEP-WCMC and IUCN (<year>), Protected
# Planet: The World Database on Protected Areas (WDPA) [On-line], <month/year
# of the version downloaded>, Cambridge, UK: UNEP-WCMC and IUCN.

WDPA_AREA_TOLERANCE = 0.05   # tiled area must be within 5% of the WDPA GIS_AREA ...
# ... or within WDPA_TILE_TOLERANCE tiles' worth, whichever is larger. Small sites
# are quantised to whole tiles (one res-7 tile = 5.16 km2 = 3.6% of a 143 km2
# site), so a fixed 5% fails correct boundaries. Empirical setting: for 400 random
# irregular 143 km2 polygons the largest tiling error was 4.6 tiles; 5 tiles
# covered all of them (ARF validation, 2026-10-08). Above ~3,600 km2 the 5% term
# dominates (largest error seen there: 2.1%).
WDPA_TILE_TOLERANCE = 5
WDPA_REP_AREA_TOLERANCE = 0.005   # file REP_AREA vs the value on the site page


def _default_boundary_dir():
    drive_dir = "/content/drive/MyDrive/ARF_boundaries"
    if os.path.isdir(drive_dir):
        return drive_dir
    return os.path.join(os.getcwd(), "ARF_boundaries")


ARF_BOUNDARY_DIR = globals().get("ARF_BOUNDARY_DIR") or _default_boundary_dir()


_MONTHS = {m: i for i, m in enumerate(("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 1)}


def _wdpa_release(path):
    """(year, month) of the WDPA release in a filename like ..._Oct2026_..., or (0, 0)."""
    import re
    m = re.search(r"_([A-Z][a-z]{2})(\d{4})_", os.path.basename(path))
    return (int(m.group(2)), _MONTHS.get(m.group(1), 0)) if m else (0, 0)


def _find_wdpa_files(site_id, boundary_dir):
    """
    Every polygon shapefile part for a WDPA site under boundary_dir, at ANY depth.
    Protected Planet splits downloads into numbered parts (..._shp_0, _shp_1, ...)
    nested inside the download folder, and may deliver them as zips inside a zip;
    zips are unpacked (only this site's members) into a temporary folder. If
    several monthly releases are present, only the newest release's parts are
    returned. Returns [] if nothing is found.
    """
    import glob, zipfile, tempfile
    target = f"_{site_id}_shp-polygons.shp"
    found = [p for p in glob.glob(os.path.join(boundary_dir, "**", "*.shp"), recursive=True)
             if p.endswith(target)]
    zips = [p for p in glob.glob(os.path.join(boundary_dir, "**", "*.zip"), recursive=True)
            if f"_{site_id}_" in os.path.basename(p)]
    for zp in zips:
        out = tempfile.mkdtemp(prefix=f"wdpa_{site_id}_")
        def _unpack(zf, dest):
            for member in zf.namelist():
                base = os.path.basename(member)
                if member.endswith(".zip"):
                    import io
                    with zipfile.ZipFile(io.BytesIO(zf.read(member))) as inner:
                        _unpack(inner, os.path.join(dest, os.path.splitext(base)[0]))
                elif f"_{site_id}_shp-polygons." in base:
                    os.makedirs(dest, exist_ok=True)
                    with open(os.path.join(dest, base), "wb") as fh:
                        fh.write(zf.read(member))
        with zipfile.ZipFile(zp) as zf:
            _unpack(zf, out)
        found += [p for p in glob.glob(os.path.join(out, "**", "*.shp"), recursive=True) if p.endswith(target)]
    if not found:
        return []
    newest = max(_wdpa_release(p) for p in found)
    # De-duplicate by content: the same part may exist both extracted and inside
    # its zip; reading it twice would double-count the site's area.
    import hashlib
    unique = {}
    for p in sorted(p for p in found if _wdpa_release(p) == newest):
        with open(p, "rb") as fh:
            unique.setdefault(hashlib.sha256(fh.read()).hexdigest(), p)
    return sorted(unique.values())


def _find_wdpa_file(site_id, boundary_dir):
    """First matching file (kept for the manifest in Cell 9); see _find_wdpa_files."""
    files = _find_wdpa_files(site_id, boundary_dir)
    return files[0] if files else None


def _wdpa_sites_present(boundary_dir):
    """WDPA site IDs for which polygon shapefiles exist anywhere under boundary_dir."""
    import glob, re
    ids = set()
    for p in glob.glob(os.path.join(boundary_dir, "**", "*_shp-polygons.shp"), recursive=True):
        m = re.search(r"_(\d+)_shp-polygons\.shp$", p)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def load_wdpa_site(site_id, spec, boundary_dir, h3_res, geoms_out=None):
    """
    Read one WDPA site's boundary, check it is the intended site, and tile it.
    Returns (tile_set, report_dict). If geoms_out is a list, the boundary
    polygons (GeoJSON-like) are appended to it for the maps; they are kept in
    memory only and never written to the outputs (WDPA licence). Raises ConfigurationError with download
    instructions if the file is missing, or if its designation does not match
    (e.g. a Ramsar file supplied where the World Heritage Site is required).
    """
    import shapefile   # pyshp
    paths = _find_wdpa_files(site_id, boundary_dir)
    if not paths:
        import glob
        points = [p for p in glob.glob(os.path.join(boundary_dir, "**", f"*_{site_id}_shp-points.shp"), recursive=True)]
        if points:
            # WDPA Manual 1.6, p.15: "Where boundary data are unavailable, the
            # latitude and longitude of the centremost point of the site is
            # requested ... instead." A point has no shape, so it cannot be tiled.
            raise ConfigurationError(
                f"WDPA site {site_id} ({spec['name']}) has NO boundary polygon in the WDPA: the download "
                f"contains only a point record ({os.path.basename(points[0])}). WDPA stores a point where "
                f"the data provider supplied no boundary (WDPA/WDOECM Manual 1.6, p.15), so this site's "
                f"shape is unknown and it cannot be tiled. Choose a site whose download contains "
                f"'_shp-polygons' files, or remove '{spec['name']}' from ROI_DEFINITIONS.")
        present = _wdpa_sites_present(boundary_dir) if os.path.isdir(boundary_dir) else []
        raise ConfigurationError(
            f"WDPA boundary for site {site_id} ({spec['name']}) not found anywhere under {boundary_dir}. "
            f"Download it as a shapefile from https://www.protectedplanet.net/{site_id} and place the "
            f"download (zip, or extracted folders at any depth) in that folder. "
            f"WDPA sites found there: {present if present else 'none'}.")
    path = paths[0]
    tiles, gis_area, desig, name, rep_areas = set(), 0.0, set(), set(), set()
    shp_records = []
    for p in paths:                                   # every split part of this site
        reader = shapefile.Reader(p)
        fields = [f[0] for f in reader.fields[1:]]
        shp_records += [(fields, sr) for sr in reader.iterShapeRecords()]
    for fields, shp_rec in shp_records:
        attrs = dict(zip(fields, shp_rec.record))
        rec_id = int(attrs.get("SITE_ID", attrs.get("WDPAID", -1)))
        if rec_id != int(site_id):
            raise ConfigurationError(f"{os.path.basename(path)} holds WDPA site {rec_id}, expected {site_id}.")
        desig.add(str(attrs.get("DESIG_ENG", "")))
        name.add(str(attrs.get("NAME_ENG", attrs.get("NAME", ""))))
        gis_area += float(attrs.get("GIS_AREA", 0.0))
        rep_areas.add(round(float(attrs.get("REP_AREA", 0.0)), 4))
        tiles |= set(h3.geo_to_cells(shp_rec.shape.__geo_interface__, h3_res))
        if geoms_out is not None:
            geoms_out.append(shp_rec.shape.__geo_interface__)
    # Identity check 1: reported area (site-level, repeated on every part) must
    # match the value on the site's Protected Planet page, when one is given.
    exp_rep = spec.get("expected_rep_area_km2")
    if exp_rep:
        rep = max(rep_areas) if rep_areas else 0.0
        if abs(rep - exp_rep) / exp_rep > WDPA_REP_AREA_TOLERANCE:
            raise ConfigurationError(
                f"WDPA site {site_id}: file reports REP_AREA {rep:,.2f} km2, but {spec.get('source_url', 'the site page')} "
                f"shows {exp_rep:,.2f} km2. This is probably not {spec['name']}.")
    # Identity check 2: designation, when the ROI specifies one. Otherwise the
    # file's designation is only recorded.
    if spec.get("designation_contains") is None:
        logger.info(f"  WDPA {site_id}: designation {sorted(desig)} (recorded, not pre-specified).")
    elif not any(spec["designation_contains"].lower() in d.lower() for d in desig):
        raise ConfigurationError(
            f"WDPA site {site_id} is designated {sorted(desig)}, but this ROI requires "
            f"'{spec['designation_contains']}' ({spec['name']}). Check the site ID: Protected Planet "
            f"lists several overlapping designations under similar names.")
    tiled_area = float(sum(h3.cell_area(c, unit="km^2") for c in tiles))
    deviation = (tiled_area - gis_area) / gis_area if gis_area else float("nan")
    tile_km2 = tiled_area / len(tiles) if tiles else 0.0
    tolerance = max(WDPA_AREA_TOLERANCE, WDPA_TILE_TOLERANCE * tile_km2 / gis_area if gis_area else 0.0)
    report = {"site_id": int(site_id), "name": sorted(name), "designation": sorted(desig),
              "file": os.path.basename(path), "n_parts": len(paths), "wdpa_gis_area_km2": gis_area,
              "expected_gis_area_km2": spec.get("expected_gis_area_km2"),
              "tiled_area_km2": tiled_area, "area_deviation_pct": 100 * deviation, "n_tiles": len(tiles),
              "area_tolerance_pct": 100 * tolerance, "wdpa_rep_area_km2": max(rep_areas) if rep_areas else None,
              "source_url": spec.get("source_url"), "accessed": spec.get("accessed")}
    if not abs(deviation) <= tolerance:
        raise SpatialDiscretisationError(
            f"WDPA site {site_id}: tiled area {tiled_area:,.0f} km2 differs from WDPA GIS_AREA "
            f"{gis_area:,.0f} km2 by {100 * deviation:+.1f}% (tolerance {100 * tolerance:.1f}%).")
    exp = spec.get("expected_gis_area_km2")
    if exp and abs(gis_area - exp) / exp > WDPA_AREA_TOLERANCE:
        logger.warning(f"WDPA site {site_id}: GIS_AREA {gis_area:,.0f} km2 differs from the expected "
                       f"{exp:,.0f} km2 -- the boundary may have been revised in this WDPA release.")
    return tiles, report


def discretise_wdpa_roi(roi_def, h3_res, boundary_dir=None, geoms_out=None):
    """Union of the H3 tiles of every WDPA site in roi_def['wdpa_sites'].
    geoms_out: optional list collecting the boundary polygons (see load_wdpa_site)."""
    boundary_dir = boundary_dir or ARF_BOUNDARY_DIR
    all_tiles, reports = set(), []
    for site_id, spec in roi_def["wdpa_sites"].items():
        tiles, rep = load_wdpa_site(site_id, spec, boundary_dir, h3_res, geoms_out=geoms_out)
        all_tiles |= tiles
        reports.append(rep)
        logger.info(f"  WDPA {site_id} ({', '.join(rep['name'])}; {', '.join(rep['designation'])}): "
                    f"{rep['n_tiles']:,} tiles, {rep['tiled_area_km2']:,.0f} km2 vs WDPA "
                    f"{rep['wdpa_gis_area_km2']:,.0f} km2 ({rep['area_deviation_pct']:+.2f}%)")
    return all_tiles, reports


def resolve_roi_coords(roi_name: str) -> List[Tuple[float, float]]:
    """Resolve a named ROI to its (lat, lon) polygon."""
    if roi_name not in ROI_DEFINITIONS:
        raise ConfigurationError(f"ACTIVE_ROI={roi_name!r} is not defined in ROI_DEFINITIONS.")

    roi_def = ROI_DEFINITIONS[roi_name]
    if "roi_coords" in roi_def:
        return roi_def["roi_coords"]
    elif "waypoints_latlon" in roi_def:
        return _buffer_waypoints_to_polygon(roi_def["waypoints_latlon"], roi_def["buffer_km"])
    elif "wdpa_sites" in roi_def:
        return None          # boundary polygons are loaded from WDPA files (see section 5)
    else:
        raise ConfigurationError(f"ROI {roi_name!r} lacks coordinate definitions.")

# Resolve active ROI from Cell 0
_active_roi = ACTIVE_ROI if 'ACTIVE_ROI' in globals() else "botswana_full"
roi_coords = resolve_roi_coords(_active_roi)

# =============================================================================
# 5. H3 DGGS SPATIAL DISCRETISATION
# =============================================================================

_h3_res = H3_RESOLUTION if 'H3_RESOLUTION' in globals() else 7

roi_boundary_report = None
# ROI outline for the maps (Cells 3 and 7): the official boundary polygons, or
# the ROI polygon. Kept in memory only -- WDPA boundaries must never be written
# to the outputs (the WDPA licence forbids redistributing the data).
ROI_OUTLINE_GEOMS: List[dict] = []
if "wdpa_sites" in ROI_DEFINITIONS[_active_roi]:
    # Configuration/discretisation errors here carry their own instructions.
    target_tiles, roi_boundary_report = discretise_wdpa_roi(ROI_DEFINITIONS[_active_roi], _h3_res,
                                                            geoms_out=ROI_OUTLINE_GEOMS)
    ROI_OUTLINE_SOURCE = "WDPA"
else:
    ROI_OUTLINE_GEOMS = [{"type": "Polygon", "coordinates": [[(lon, lat) for lat, lon in roi_coords]]}]
    ROI_OUTLINE_SOURCE = "ROI definition"
    try:
        geo_polygon = h3.LatLngPoly(roi_coords)
        target_tiles: Set[str] = set(h3.polygon_to_cells(geo_polygon, _h3_res))
    except Exception as e:
        raise SpatialDiscretisationError(f"Failed to discretise ROI '{_active_roi}' at H3 Res {_h3_res}: {e}") from e

_n_tiles = len(target_tiles)
if _n_tiles == 0:
    raise SpatialDiscretisationError(f"Discretisation of '{_active_roi}' returned 0 tiles.")

logger.info(f"━ ROI: {_active_roi} discretised into {_n_tiles:,} tiles at Resolution {_h3_res}.")

# =============================================================================
# 6. GROUND STATION NETWORK (Downlink Physical Constraints)
# =============================================================================

GROUND_STATIONS = ROI_DEFINITIONS[_active_roi].get("ground_stations", {})
if not GROUND_STATIONS:
    logger.warning(f"ROI {_active_roi!r} defines no ground stations; D_downlink will be zero.")

ground_station_points: Dict[str, GeodeticPoint] = {
    name: GeodeticPoint(FastMath.toRadians(d['lat']), FastMath.toRadians(d['lng']), d['alt_m'])
    for name, d in GROUND_STATIONS.items()
}

# =============================================================================
# SUMMARY & PERSISTENCE
# =============================================================================

# v1.5.2: the smallest cell ID, not next(iter(set)): string hashing is salted
# per Python process, so a set's first element (and this JSON) changed between
# identical runs.
_sample_tile = min(target_tiles)
_sample_lat, _sample_lng = h3.cell_to_latlng(_sample_tile)

_validation_results = {
    "jvm_initialised": True,
    "orekit_data_loaded": _orekit_data_loaded,
    "h3_version": h3.__version__,
    "roi_active": _active_roi,
    "roi_tiles": _n_tiles,
    "sample_tile": _sample_tile,
    "sample_tile_centroid": (_sample_lat, _sample_lng),
    "package_versions": ARF_PACKAGE_VERSIONS,
    "validated_versions": ARF_VALIDATED_VERSIONS,
    "ground_stations": len(GROUND_STATIONS),
    "wdpa_boundaries": roi_boundary_report,
    "h3_resolution": _h3_res,
}

print("\n" + "=" * 60)
print("  CELL 1: Environment Setup & Spatial Discretisation")
print("=" * 60)
print(f"  JVM Status:                Initialised ({platform.system()})")
print(f"  Orekit Data:               Loaded Successfully")
print(f"  H3 Version:                {h3.__version__}")
print(f"  Reference Frames:          ITRF 2010, WGS84 Ellipsoid")
print(f"  ─────────────────────────────────────────────────────────")
print(f"  Target ROI:                {_active_roi}")
print(f"  Discretised Tiles:         {_n_tiles:,} (Res {_h3_res})")
print(f"  Ground Stations Configured:{len(GROUND_STATIONS)} ({', '.join(GROUND_STATIONS.keys())})")
print("=" * 60)

validation_filepath = os.path.join(_out_dir, 'arf_cell1_validation.json')
try:
    with open(validation_filepath, 'w') as f:
        json.dump(_validation_results, f, indent=2, default=str)
    logger.info(f"Validation results saved to {validation_filepath}")
except Exception as e:
    logger.warning(f"Could not save validation results: {e}")

print("\n=== CELL 1 COMPLETE ===")

