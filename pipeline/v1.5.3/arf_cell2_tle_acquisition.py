# =============================================================================
# CELL 2: TLE ACQUISITION & DETERMINISTIC SATELLITE SELECTION
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Ingests orbital elements for the imaging fleet and builds SGP4 propagators
#   (Orekit TLEPropagator).
#
#   Real fleets (default): CelesTrak GP data as OMM JSON (catalogue numbers
#   above 99999 included), live from the newest element sets or from
#   REFERENCE_EPOCH_UTC; past windows (ANALYSIS_EPOCH) use Space-Track
#   gp_history elements as of that instant. EXCLUDE_SATELLITES (v1.5.2) leaves
#   named satellites out for a fleet sensitivity run.
#   Starlink proxy (legacy): TLE text, live (celestrak_live) or at a disaster
#   onset (spacetrack_history, e.g. the Nepal flash flood).
#
# Reproducibility:
#   Real fleets: every operational satellite with a sourced MISSION_TABLE entry,
#   from CelesTrak GP data cached for 6 h (ARF_cache), or from Space-Track
#   gp_history as of ANALYSIS_EPOCH (retrospective). The fleet list, the element
#   files and the window start are recorded (Cell 9). Runs to be compared must
#   share the window start (REFERENCE_EPOCH_UTC or ANALYSIS_EPOCH, Cell 0).
#   Starlink proxy (legacy): NORAD ID pinning of a seeded sample of N.
#
# Dependencies:
#   - Cell 0 (provides ACTIVE_N, ACTIVE_ROI, OUT_DIR)
#   - Cell 1 (provides ROI_DEFINITIONS, itrf, earth)
#
# Outputs (Globals for Downstream Pipeline):
#   - constellation: List[Dict] of active satellite objects and TLE propagators
#   - sim_epoch: Orekit AbsoluteDate marking the simulation origin
#   - CONSTELLATION_ALTITUDE_KM: Derived physical median altitude
#   - T_MAX_WAIT_MINUTES: Dynamically derived orbital period
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 2 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["2"] = _ARF_CELL_VERSION

import os
import sys
import json
import time
import logging
import urllib.request
import urllib.parse
import numpy as np
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.hipparchus.util import FastMath

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logger = logging.getLogger("arf.cell2")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# =============================================================================
# CACHE & SOURCE CONFIGURATION
# =============================================================================

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()
# v1.5.0: ARF_SHARED_CACHE_DIR (set by a batch driver) shares the CelesTrak group
# lists and Space-Track element sets between runs, so a batch does not download
# the same data once per run; every cache file name carries its ROI, fleet or
# epoch, so runs cannot read each other's selections.
ARF_CACHE_DIR = globals().get("ARF_SHARED_CACHE_DIR") or os.path.join(_out_dir, "ARF_cache")
os.makedirs(ARF_CACHE_DIR, exist_ok=True)

CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle"
TLE_RAW_CACHE_PATH = os.path.join(ARF_CACHE_DIR, f"arf_raw_cache_{ACTIVE_ROI}.json")
SELECTION_CACHE_PATH = os.path.join(ARF_CACHE_DIR, f"arf_selection_{ACTIVE_ROI}_N{ACTIVE_N}.json")
TLE_CACHE_MAX_AGE_HOURS = 6

_roi_def = ROI_DEFINITIONS[ACTIVE_ROI]
TLE_SOURCE = _roi_def.get("tle_source", "celestrak_live")

# =============================================================================
# 1A. LIVE FETCH (CELESTRAK)
# =============================================================================

def fetch_tle_data(url: str = CELESTRAK_URL) -> List[Tuple[str, str, str]]:
    logger.info(f"Fetching live TLE data from: {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "ARF-Research-Pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw_text = resp.read().decode("utf-8")

    lines = [l.strip() for l in raw_text.splitlines() if l.strip()]
    records = []
    for i in range(0, len(lines), 3):
        records.append((lines[i], lines[i+1], lines[i+2]))
    return records

def load_tle_records(url=CELESTRAK_URL, cache_path=TLE_RAW_CACHE_PATH, max_age_hours=TLE_CACHE_MAX_AGE_HOURS, force_refresh=False):
    if os.path.exists(cache_path) and not force_refresh:
        with open(cache_path) as f:
            cache = json.load(f)
        age_hours = (time.time() - cache["fetched_at"]) / 3600.0
        if age_hours < max_age_hours:
            logger.info(f"Using cached TLE fetch ({age_hours:.1f}h old): {len(cache['records'])} records.")
            return [tuple(r) for r in cache["records"]]

    records = fetch_tle_data(url)
    logger.info(f"Total records fetched: {len(records)}")
    try:
        with open(cache_path, "w") as f:
            json.dump({"fetched_at": time.time(), "url": url, "records": records}, f)
    except Exception as e:
        logger.warning(f"Could not write TLE raw cache: {e}")
    return records

# =============================================================================
# 1B. HISTORICAL FETCH (SPACE-TRACK)
# =============================================================================

def _get_spacetrack_credentials():
    try:
        from google.colab import userdata
        username = userdata.get('SPACETRACK_USERNAME')
        password = userdata.get('SPACETRACK_PASSWORD')
    except Exception as e:
        raise RuntimeError("Could not read Space-Track credentials from Colab userdata. "
                           "Set SPACETRACK_USERNAME and SPACETRACK_PASSWORD in the Secrets tab.")
    if not username or not password:
        raise RuntimeError("SPACETRACK_USERNAME or SPACETRACK_PASSWORD is empty.")
    return username, password

def _spacetrack_login(username, password, base_url="https://www.space-track.org"):
    import http.cookiejar
    cookie_jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cookie_jar))
    login_url = f"{base_url}/ajaxauth/login"
    payload = urllib.parse.urlencode({"identity": username, "password": password}).encode()
    req = urllib.request.Request(login_url, data=payload, headers={"User-Agent": "ARF-Research-Pipeline/1.0"})

    resp = opener.open(req, timeout=30)
    if resp.getcode() != 200:
        raise RuntimeError(f"Space-Track auth failed (HTTP {resp.getcode()})")
    # Space-Track signals bad credentials with HTTP 200 and a body of exactly
    # {"Login": "Failed"} (success returns an empty body plus the session
    # cookie), so the status code alone cannot detect a failed login.
    body = resp.read().decode("utf-8").strip()
    if body:
        try:
            parsed = json.loads(body)
            if isinstance(parsed, dict) and parsed.get("Login") == "Failed":
                raise RuntimeError("Space-Track authentication rejected: "
                                   "check SPACETRACK_USERNAME / SPACETRACK_PASSWORD.")
        except json.JSONDecodeError:
            pass
    return opener

def fetch_starlink_gp_history(opener, epoch_date, base_url="https://www.space-track.org", lookback_days=30):
    from datetime import datetime, timedelta
    end_dt = datetime.strptime(epoch_date, "%Y-%m-%dT%H:%M:%S.%f") if "." in epoch_date else datetime.strptime(epoch_date, "%Y-%m-%dT%H:%M:%S")
    start_dt = end_dt - timedelta(days=lookback_days)
    start_str = start_dt.strftime("%Y-%m-%dT%H:%M:%S")

    query_path = f"/basicspacedata/query/class/gp_history/OBJECT_NAME/~~STARLINK/EPOCH/{start_str}--{epoch_date}/orderby/EPOCH%20desc/format/json"
    req = urllib.request.Request(base_url + query_path, headers={"User-Agent": "ARF-Research-Pipeline/1.0"})

    with opener.open(req, timeout=120) as resp:
        parsed = json.loads(resp.read().decode("utf-8"))

    if not parsed or (isinstance(parsed, dict) and "error" in parsed):
        raise RuntimeError(f"Space-Track gp_history query error or empty result.")
    return parsed

def _parse_gp_history_epoch(epoch_str):
    from datetime import datetime
    try:
        return datetime.strptime(epoch_str, "%Y-%m-%dT%H:%M:%S.%f")
    except ValueError:
        return datetime.strptime(epoch_str, "%Y-%m-%dT%H:%M:%S")

def _dedupe_gp_history_to_latest_per_satellite(gp_history_records):
    by_epoch_desc = sorted(gp_history_records, key=lambda r: _parse_gp_history_epoch(r["EPOCH"]), reverse=True)
    by_id_stable = sorted(by_epoch_desc, key=lambda r: r["NORAD_CAT_ID"])
    seen, out = set(), []
    for rec in by_id_stable:
        nid = rec["NORAD_CAT_ID"]
        if nid not in seen:
            seen.add(nid)
            out.append((rec["OBJECT_NAME"], rec["TLE_LINE1"], rec["TLE_LINE2"]))
    return out

def load_historical_tle_records(epoch_date, cache_path=TLE_RAW_CACHE_PATH, force_refresh=False):
    if os.path.exists(cache_path) and not force_refresh:
        with open(cache_path) as f:
            cache = json.load(f)
        logger.info(f"Using cached Space-Track gp_history fetch (epoch <= {epoch_date}).")
        return [tuple(r) for r in cache["records"]]

    username, password = _get_spacetrack_credentials()
    opener = _spacetrack_login(username, password)
    raw_records = fetch_starlink_gp_history(opener, epoch_date)
    records = _dedupe_gp_history_to_latest_per_satellite(raw_records)
    logger.info(f"Deduplicated to {len(records)} satellites at historical epoch.")

    try:
        with open(cache_path, "w") as f:
            json.dump({"fetched_at": time.time(), "epoch_date": epoch_date, "records": records}, f)
    except Exception as e:
        logger.warning(f"Could not write Space-Track cache: {e}")
    return records

# =============================================================================
# 2. EXECUTE FETCH & DETERMINISTIC PINNING
# =============================================================================

if FLEET != "starlink" and TLE_SOURCE == "spacetrack_history" and not globals().get("ANALYSIS_EPOCH"):
    # A historical ROI with a real fleet is a retrospective run from the ROI's onset.
    from datetime import datetime as _dt_hist
    ANALYSIS_EPOCH = _dt_hist.fromisoformat(_roi_def["historical_epoch"]).strftime("%Y-%m-%dT%H:%M:%S.%f")
    logger.info(f"Historical ROI {ACTIVE_ROI!r} with fleet {FLEET!r}: retrospective run from {ANALYSIS_EPOCH} UTC.")
if FLEET != "starlink":
    raw_records = []        # real fleets use CelesTrak GP JSON (section 2B)
elif TLE_SOURCE == "spacetrack_history":
    _historical_epoch = _roi_def["historical_epoch"]
    logger.info(f"━ ACTIVE_ROI={ACTIVE_ROI!r} requires historical TLEs (epoch <= {_historical_epoch}).")
    raw_records = load_historical_tle_records(_historical_epoch)
else:
    raw_records = load_tle_records()

import re

# =============================================================================
# 2B. REAL EARTH-OBSERVATION FLEETS (CelesTrak GP data, JSON)
# =============================================================================
# JSON (OMM) instead of TLE text: CelesTrak cannot deliver catalogue numbers
# above 99999 as TLEs, so newly catalogued satellites would silently vanish
# from a TLE fetch. Orekit TLE objects are built from the numeric elements.

def omm_to_tle(omm, utc):
    """
    Orekit TLE from a CelesTrak OMM JSON record (GP data, FORMAT=json). Built
    from the numeric elements, not TLE text, so catalogue numbers above 99999
    are supported. OMM units: MEAN_MOTION rev/day; MEAN_MOTION_DOT and
    MEAN_MOTION_DDOT as in the TLE fields (rev/day^2, rev/day^3); angles in
    degrees; EPOCH UTC ISO-8601. Orekit expects rad/s, rad/s^2, rad/s^3, radians.
    """
    from org.orekit.propagation.analytical.tle import TLE
    from org.orekit.time import AbsoluteDate
    from jpype import JChar
    import math
    two_pi, day = 2.0 * math.pi, 86400.0
    obj_id = str(omm.get("OBJECT_ID", "") or "")            # e.g. "1998-067A"
    try:
        launch_year, rest = int(obj_id[:4]), obj_id[5:]
        launch_number = int("".join(ch for ch in rest if ch.isdigit()) or 0)
        launch_piece = "".join(ch for ch in rest if ch.isalpha()) or "A"
    except ValueError:
        launch_year, launch_number, launch_piece = 2000, 0, "A"
    return TLE(int(omm["NORAD_CAT_ID"]), JChar(str(omm.get("CLASSIFICATION_TYPE", "U"))[0]),
               launch_year, launch_number, launch_piece,
               int(omm.get("EPHEMERIS_TYPE", 0)), int(omm.get("ELEMENT_SET_NO", 999)),
               AbsoluteDate(str(omm["EPOCH"]).replace("Z", ""), utc),
               float(omm["MEAN_MOTION"]) * two_pi / day,
               # TLE/OMM derivative fields hold n-dot/2 and n-ddot/6; Orekit stores
               # the derivatives themselves (verified against Orekit's text parser).
               2.0 * float(omm["MEAN_MOTION_DOT"]) * two_pi / day ** 2,
               6.0 * float(omm["MEAN_MOTION_DDOT"]) * two_pi / day ** 3,
               float(omm["ECCENTRICITY"]), math.radians(float(omm["INCLINATION"])),
               math.radians(float(omm["ARG_OF_PERICENTER"])), math.radians(float(omm["RA_OF_ASC_NODE"])),
               math.radians(float(omm["MEAN_ANOMALY"])),
               int(omm.get("REV_AT_EPOCH", 0)), float(omm["BSTAR"]), utc)


def load_gp_json(groups, cache_dir, max_age_hours=TLE_CACHE_MAX_AGE_HOURS):
    """OMM records for CelesTrak GROUPs, cached per group; de-duplicated by
    NORAD_CAT_ID across groups (a satellite may be listed in more than one)."""
    by_id = {}
    for g in groups:
        path = os.path.join(cache_dir, f"arf_gp_{g}.json")
        recs = None
        if os.path.exists(path):
            with open(path) as f:
                cache = json.load(f)
            if (time.time() - cache["fetched_at"]) / 3600.0 < max_age_hours:
                recs = cache["records"]
                logger.info(f"Using cached GP JSON for GROUP={g}: {len(recs)} records.")
        if recs is None:
            url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP={g}&FORMAT=json"
            logger.info(f"Fetching {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "ARF-Research-Pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                recs = json.loads(resp.read().decode("utf-8"))
            with open(path, "w") as f:
                json.dump({"fetched_at": time.time(), "url": url, "records": recs}, f)
        for r in recs:
            r = dict(r, _GROUP=g)
            by_id.setdefault(int(r["NORAD_CAT_ID"]), r)
    return list(by_id.values())


def load_eo_gp_history(live_records, epoch_iso, cache_dir, lookback_days=7.0, base_url="https://www.space-track.org"):
    """
    Orbital elements of a real fleet AS OF epoch_iso, for retrospective runs.
    Satellites: those in the current CelesTrak list (live_records) whose name
    matches a MISSION_TABLE entry. Elements: Space-Track gp_history, newest
    element set per satellite with lookback_days <= age <= 0 before epoch_iso
    (no element set from after the epoch is used). The CelesTrak name is kept,
    so mission matching is identical to live runs. Cached per fleet and epoch.
    Satellites launched after epoch_iso have no elements and drop out;
    satellites retired since then are not in the current list (stated limit).
    """
    from datetime import datetime, timedelta
    names = {int(r["NORAD_CAT_ID"]): str(r.get("OBJECT_NAME", "")).strip() for r in live_records}
    groups = {int(r["NORAD_CAT_ID"]): r.get("_GROUP") for r in live_records}
    ids = sorted(i for i, n in names.items() if _match_mission(n) is not None)
    end = datetime.fromisoformat(epoch_iso.replace("Z", ""))
    start = end - timedelta(days=lookback_days)
    path = os.path.join(cache_dir, f"arf_gp_history_{FLEET_TAG}_{end:%Y%m%dT%H%M%S}.json")
    globals()["GP_HISTORY_PATH"] = path                   # recorded in the manifest (Cell 9)
    if os.path.exists(path):
        with open(path) as f:
            recs = json.load(f)["records"]
        logger.info(f"Using cached Space-Track gp_history for {len(ids)} satellites (epoch <= {epoch_iso}).")
    else:
        if not ids:
            raise ConfigurationError("No satellite in the current CelesTrak list matches a MISSION_TABLE entry.")
        username, password = _get_spacetrack_credentials()
        opener = _spacetrack_login(username, password)
        recs = []
        for i in range(0, len(ids), 100):                  # keep each URL short
            chunk = ",".join(str(x) for x in ids[i:i + 100])
            q = (f"/basicspacedata/query/class/gp_history/NORAD_CAT_ID/{chunk}/EPOCH/"
                 f"{start:%Y-%m-%dT%H:%M:%S}--{end:%Y-%m-%dT%H:%M:%S}/orderby/EPOCH%20desc/format/json")
            req = urllib.request.Request(base_url + q, headers={"User-Agent": "ARF-Research-Pipeline/1.0"})
            with opener.open(req, timeout=120) as resp:
                part = json.loads(resp.read().decode("utf-8"))
            if isinstance(part, dict) and "error" in part:
                raise RuntimeError(f"Space-Track gp_history error: {part}")
            recs.extend(part)
        with open(path, "w") as f:
            json.dump({"fetched_at": time.time(), "epoch": epoch_iso, "norad_ids": ids, "records": recs}, f)
    latest = {}
    for r in recs:
        nid = int(r["NORAD_CAT_ID"])
        ep = _omm_epoch(r)
        if ep.replace(tzinfo=None) > end:                  # defensive: never use a later element set
            continue
        if nid not in latest or ep > _omm_epoch(latest[nid]):
            latest[nid] = r
    out = [dict(r, OBJECT_NAME=names.get(nid, r.get("OBJECT_NAME", "")), _GROUP=groups.get(nid))
           for nid, r in sorted(latest.items())]
    logger.info(f"Retrospective elements at {epoch_iso}: {len(out)} of {len(ids)} sourced satellites have an "
                f"element set within {lookback_days:g} days before the epoch.")
    return out


def _omm_orbit_km(omm, mu=398600.4418, re=6378.137):
    """(perigee_alt_km, apogee_alt_km) from OMM mean motion (rev/day) and eccentricity."""
    import math
    n = float(omm["MEAN_MOTION"]) * 2.0 * math.pi / 86400.0
    a = (mu / n ** 2) ** (1.0 / 3.0)
    e = float(omm["ECCENTRICITY"])
    return a * (1.0 - e) - re, a * (1.0 + e) - re


def swath_from_fov_km(altitude_km, fov_deg, re=6378.137):
    """Ground swath (km) of a nadir-pointing sensor with full field of view
    fov_deg at altitude_km over a sphere of radius re: twice the Earth-central
    angle to the edge of the field, lambda = arcsin((re+h)/re * sin(fov/2)) - fov/2."""
    import math
    eta = math.radians(fov_deg / 2.0)
    s = (re + altitude_km) / re * math.sin(eta)
    if s >= 1.0:
        raise ValueError(f"field of view {fov_deg} deg reaches beyond the horizon at {altitude_km} km")
    return 2.0 * re * (math.asin(s) - eta)


def _omm_epoch(omm):
    from datetime import datetime, timezone
    return datetime.fromisoformat(str(omm["EPOCH"]).replace("Z", "")).replace(tzinfo=timezone.utc)


def _match_mission(name):
    for m in MISSION_TABLE:
        if re.match(m["name_pattern"], str(name).strip(), flags=re.IGNORECASE):
            return m
    return None


def build_eo_constellation(omm_records, reference_iso, min_perigee_km=300.0, max_age_days=7.0, min_age_days=0.0,
                           exclude_names=()):
    """
    Every satellite that (1) matches a MISSION_TABLE entry, (2) has a perigee of
    at least min_perigee_km and (3) an element set at most max_age_days older
    than the reference epoch (and at most -min_age_days newer: 0 by default, -2
    when Cell 0 sets a live window start, v1.5.1). No sampling: N is the whole
    operational fleet. exclude_names (v1.5.2): OBJECT_NAMEs left out on
    purpose (EXCLUDE_SATELLITES, Cell 0), counted in stats. Returns
    (constellation, stats); stats lists unmatched names, never guessed.
    """
    import collections, jpype
    from datetime import datetime, timezone
    from org.orekit.time import TimeScalesFactory
    utc = TimeScalesFactory.getUTC()
    ref = datetime.fromisoformat(reference_iso).replace(tzinfo=timezone.utc)
    excl = {str(n).strip().upper() for n in exclude_names}
    stats = {"fetched": len(omm_records), "malformed": 0, "stale_or_future": 0, "low_perigee": 0,
             "high_apogee": 0, "unmatched_mission": 0, "excluded_by_config": 0, "operational": 0,
             "criteria": {"min_perigee_km": min_perigee_km, "max_elset_age_days": max_age_days,
                          "min_elset_age_days": min_age_days, "exclude_satellites": sorted(excl)}, "seed": None}
    unmatched, per_mission, constellation, seen_names = collections.Counter(), collections.Counter(), [], set()
    unmatched_names, excluded_found = [], []
    for omm in sorted(omm_records, key=lambda r: int(r["NORAD_CAT_ID"])):
        name = str(omm.get("OBJECT_NAME", "")).strip()
        m = _match_mission(name)
        if m is None:
            stats["unmatched_mission"] += 1
            unmatched_names.append(name)
            unmatched[re.sub(r"[\s\-]*\d+[A-Z]?$", "", name) or name] += 1
            continue
        if name.upper() in excl:
            stats["excluded_by_config"] += 1
            excluded_found.append(name)
            continue
        try:
            age = (ref - _omm_epoch(omm)).total_seconds() / 86400.0
            perigee, apogee = _omm_orbit_km(omm)
            tle = omm_to_tle(omm, utc)
        except (KeyError, ValueError, TypeError, ZeroDivisionError, jpype.JException):
            # Bad or missing fields in the record itself. Programming errors
            # (NameError etc.) are deliberately NOT caught, so they surface.
            stats["malformed"] += 1
            continue
        if not min_age_days <= age <= max_age_days:
            stats["stale_or_future"] += 1
            continue
        if perigee < min_perigee_km:
            stats["low_perigee"] += 1
            continue
        sid = name if name not in seen_names else f"{name} [{int(omm['NORAD_CAT_ID'])}]"
        seen_names.add(sid)
        per_mission[m["mission"]] += 1
        if "swath_km" in m:
            swath, fov = float(m["swath_km"]["value"]), None
        else:      # field of view given: swath at this satellite's mean altitude
            fov = float(m["fov_deg"]["value"])
            swath = swath_from_fov_km(0.5 * (perigee + apogee), fov)
        constellation.append({"id": sid, "tle": tle, "propagator": TLEPropagator.selectExtrapolator(tle),
                              "norad_id": int(omm["NORAD_CAT_ID"]), "group": omm.get("_GROUP"),
                              "mission": m["mission"], "swath_km": swath, "fov_deg": fov,
                              "min_sun_elev_deg": float(m["min_sun_elev_deg"]["value"])})
    stats["operational"] = len(constellation)
    stats["per_mission"] = dict(per_mission)
    stats["unmatched_families"] = dict(unmatched.most_common(40))
    stats["unmatched_names"] = sorted(set(unmatched_names))[:200]
    stats["excluded_names"] = sorted(excluded_found)
    stats["exclusions_not_found"] = sorted(excl - {n.upper() for n in excluded_found})
    return constellation, stats


def _parse_norad_id(line1):
    try:
        return int(line1[2:7])
    except (ValueError, IndexError):
        return None

# -----------------------------------------------------------------------------
# Historical runs: satellites OPERATIONAL AT THE ONSET, seeded random sample
# -----------------------------------------------------------------------------
# gp_history returns every Starlink object with an element set in the lookback
# window, including satellites that are deorbiting or no longer tracked.
# Taking the first N by catalogue number would pick the OLDEST satellites
# (2019-2022 launches), many of them decaying. Instead: keep satellites whose
# latest element set describes an operational orbit at the onset, then draw a
# reproducible random sample of N from them.
#   perigee >= 300 km : excludes deorbiting / re-entering objects
#   apogee  <= 700 km : excludes anomalous or non-Starlink-shell orbits
#   element set <= 7 days older than the onset : excludes stale / untracked
# Orbit size comes from the element set's own mean motion and eccentricity
# (semi-major axis a = (mu / n^2)^(1/3)), not from a single propagated point.
# Checked against Orekit's SGP4: propagated altitude dips ~7 km below the
# mean-element perigee (short-period J2 terms), so the 300 km bound admits
# orbits whose true minimum altitude is down to ~293 km.

OPERATIONAL_CRITERIA = {"min_perigee_km": 300.0, "max_apogee_km": 700.0, "max_elset_age_days": 7.0}


def _elset_epoch_utc(line1):
    """TLE line 1 epoch (cols 19-32, YYDDD.DDDDDDDD) as a UTC datetime."""
    from datetime import datetime, timedelta, timezone
    yy = int(line1[18:20])
    year = 2000 + yy if yy < 57 else 1900 + yy
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=float(line1[20:32]) - 1.0)


def _elset_orbit_km(line2, mu=398600.4418, re=6378.137):
    """(perigee_alt_km, apogee_alt_km) from TLE line 2: eccentricity (cols 27-33,
    implied leading decimal point) and mean motion (cols 53-63, rev/day)."""
    import math
    e = float("0." + line2[26:33].strip())
    n = float(line2[52:63]) * 2.0 * math.pi / 86400.0          # rad/s
    a = (mu / n ** 2) ** (1.0 / 3.0)
    return a * (1.0 - e) - re, a * (1.0 + e) - re


def select_operational_sample(records, onset_iso, seed, criteria=OPERATIONAL_CRITERIA):
    """
    Filter records to satellites operational at onset_iso (see criteria), put
    them in catalogue-number order, and shuffle with numpy's seeded generator.
    Returns (shuffled_records, stats); the first N are the sample.
    """
    from datetime import datetime, timezone
    onset = datetime.fromisoformat(onset_iso.replace("Z", "")).replace(tzinfo=timezone.utc)
    stats = {"fetched": len(records), "malformed": 0, "stale_or_future": 0, "low_perigee": 0,
             "high_apogee": 0, "operational": 0, "criteria": dict(criteria), "seed": seed}
    keep = []
    for rec in records:
        name, l1, l2 = rec
        try:
            age_days = (onset - _elset_epoch_utc(l1)).total_seconds() / 86400.0
            perigee, apogee = _elset_orbit_km(l2)
            nid = int(l1[2:7])
        except (ValueError, IndexError, ZeroDivisionError):
            stats["malformed"] += 1
            continue
        if not (0.0 <= age_days <= criteria["max_elset_age_days"]):
            stats["stale_or_future"] += 1
        elif perigee < criteria["min_perigee_km"]:
            stats["low_perigee"] += 1
        elif apogee > criteria["max_apogee_km"]:
            stats["high_apogee"] += 1
        else:
            keep.append((nid, rec))
    keep.sort(key=lambda x: x[0])
    stats["operational"] = len(keep)
    order = np.random.default_rng(seed).permutation(len(keep))
    return [keep[i][1] for i in order], stats


def select_pinned_satellites(records, active_n, cache_path):
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            pinned_ids = json.load(f)["norad_ids"]
        by_id = {nid: (name, l1, l2) for name, l1, l2 in records if (nid := _parse_norad_id(l1)) is not None}

        selected = [by_id[nid] for nid in pinned_ids if nid in by_id]
        if len(selected) < len(pinned_ids):
            logger.warning(f"Only {len(selected)}/{len(pinned_ids)} pinned satellites found in current fetch.")
        logger.info(f"Loaded {len(selected)} pinned satellites from {cache_path}.")
        return selected
    else:
        chosen = records[:active_n]
        norad_ids = [nid for name, l1, l2 in chosen if (nid := _parse_norad_id(l1)) is not None]
        try:
            with open(cache_path, "w") as f:
                json.dump({"active_n": active_n, "norad_ids": norad_ids}, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not write selection cache: {e}")
        logger.info(f"Pinned {len(norad_ids)} satellites to {cache_path} for deterministic re-runs.")
        return chosen

def _newest_elset_epoch_iso(records):
    """Newest element-set epoch among the fetched records (UTC, ISO 8601)."""
    epochs = []
    for _, l1, _ in records:
        try:
            epochs.append(_elset_epoch_utc(l1))
        except (ValueError, IndexError):
            continue
    if not epochs:
        raise RuntimeError("No parseable element sets in the TLE fetch.")
    return max(epochs).strftime("%Y-%m-%dT%H:%M:%S.%f")


if FLEET == "starlink":
    # Reference epoch for both the operational filter and the simulation start:
    #   historical runs -> the disaster onset (ROI_DEFINITIONS historical_epoch)
    #   live runs       -> the newest element-set epoch in the fetch (effectively
    #                      the fetch time; reproducible because the fetch is cached)
    _reference_epoch = (_historical_epoch if TLE_SOURCE == "spacetrack_history"
                        else _newest_elset_epoch_iso(raw_records))
    REFERENCE_EPOCH_SOURCE = ("ROI historical_epoch (Space-Track gp_history)" if TLE_SOURCE == "spacetrack_history"
                              else "newest element set in the Starlink fetch (live)")
    GP_INPUT_PATHS = [TLE_RAW_CACHE_PATH]

    # Same selection method for every run: satellites operational at the reference
    # epoch, seeded random sample of N. Pin file named for this method, so a pin
    # written by the earlier first-N-by-catalogue-number method is never reused.
    SELECTION_CACHE_PATH = os.path.join(
        ARF_CACHE_DIR, f"arf_selection_{ACTIVE_ROI}_N{ACTIVE_N}_opsample_seed{RANDOM_SEED}.json")
    _candidates, _op_stats = select_operational_sample(raw_records, _reference_epoch, RANDOM_SEED)
    logger.info(f"Operational at {_reference_epoch} UTC: {_op_stats['operational']} of {_op_stats['fetched']} fetched "
                f"(excluded: {_op_stats['low_perigee']} low perigee, {_op_stats['high_apogee']} high apogee, "
                f"{_op_stats['stale_or_future']} stale, {_op_stats['malformed']} malformed). "
                f"Seeded random sample (seed={RANDOM_SEED}) of up to N={ACTIVE_N}.")
    if _op_stats["operational"] < ACTIVE_N:
        logger.warning(f"Only {_op_stats['operational']} operational satellites available; using all of them.")

    selected_records = select_pinned_satellites(_candidates, ACTIVE_N, SELECTION_CACHE_PATH)

    # =============================================================================
    # 3. OREKIT SGP4 PROPAGATOR CONSTRUCTION
    # =============================================================================

    constellation = []
    for name, l1, l2 in selected_records:
        try:
            tle = TLE(l1, l2)
            constellation.append({
                "id": name,
                "tle": tle,
                "propagator": TLEPropagator.selectExtrapolator(tle)
            })
        except Exception as e:
            logger.warning(f"Skipping {name}: {e}")

else:
    # Real Earth-observation fleet: all operational satellites with a sourced
    # mission entry (section 2B). Reference epoch = newest element set (live),
    # or ANALYSIS_EPOCH with elements as of then (retrospective).
    _gp_records = load_gp_json(CELESTRAK_GROUPS[FLEET], ARF_CACHE_DIR)
    GP_INPUT_PATHS = [os.path.join(ARF_CACHE_DIR, f"arf_gp_{g}.json") for g in CELESTRAK_GROUPS[FLEET]]
    _min_age_days = 0.0
    if globals().get("ANALYSIS_EPOCH"):
        _gp_records = load_eo_gp_history(_gp_records, ANALYSIS_EPOCH, ARF_CACHE_DIR)
        _reference_epoch = ANALYSIS_EPOCH
        REFERENCE_EPOCH_SOURCE = "analysis_epoch (retrospective: elements as of then, Space-Track gp_history)"
        GP_INPUT_PATHS.append(globals().get("GP_HISTORY_PATH"))
    elif globals().get("REFERENCE_EPOCH_UTC"):
        _reference_epoch = REFERENCE_EPOCH_UTC
        REFERENCE_EPOCH_SOURCE = "reference_epoch_utc (set in Cell 0; live elements)"
        _min_age_days = -2.0          # live elements may be up to 2 days newer than a set start
    else:
        # v1.5.1: the newest element set over BOTH real-fleet groups, whichever
        # fleet is run, so the fleets of one session share their window start.
        _ref_groups = sorted({g for gs in CELESTRAK_GROUPS.values() for g in gs})
        _ref_records = load_gp_json(_ref_groups, ARF_CACHE_DIR)
        GP_INPUT_PATHS = sorted(set(GP_INPUT_PATHS) | {os.path.join(ARF_CACHE_DIR, f"arf_gp_{g}.json")
                                                       for g in _ref_groups})
        _reference_epoch = max(_omm_epoch(r) for r in _ref_records).strftime("%Y-%m-%dT%H:%M:%S.%f")
        REFERENCE_EPOCH_SOURCE = f"newest element set over CelesTrak groups {_ref_groups} (live)"
    constellation, _op_stats = build_eo_constellation(_gp_records, _reference_epoch, min_age_days=_min_age_days,
                                                      exclude_names=globals().get("EXCLUDE_SATELLITES", ()))
    if _op_stats["exclusions_not_found"]:
        logger.warning(f"EXCLUDE_SATELLITES names not found in the fleet: {_op_stats['exclusions_not_found']} "
                       f"(names are matched exactly against CelesTrak OBJECT_NAME).")
    TLE_RAW_CACHE_PATH = os.path.join(ARF_CACHE_DIR, f"arf_gp_{CELESTRAK_GROUPS[FLEET][0]}.json")
    SELECTION_CACHE_PATH = os.path.join(ARF_CACHE_DIR, f"arf_fleet_{FLEET_TAG}.json")
    with open(SELECTION_CACHE_PATH, "w") as f:
        json.dump({"fleet": FLEET, "reference_epoch": _reference_epoch,
                   "satellites": [{"id": s["id"], "norad_id": s["norad_id"], "mission": s["mission"],
                                   "swath_km": round(s["swath_km"], 3), "fov_deg": s.get("fov_deg"),
                                   "min_sun_elev_deg": s["min_sun_elev_deg"]}
                                  for s in constellation]}, f, indent=2)
    logger.info(f"Fleet {FLEET!r}: {_op_stats['operational']} of {_op_stats['fetched']} satellites used "
                f"{_op_stats['per_mission']}; excluded {_op_stats['unmatched_mission']} with no sourced mission "
                f"entry, {_op_stats['low_perigee']} low perigee, {_op_stats['stale_or_future']} stale, "
                f"{_op_stats['malformed']} malformed"
                + (f", {_op_stats['excluded_by_config']} left out by EXCLUDE_SATELLITES "
                   f"({', '.join(_op_stats['excluded_names'])})" if _op_stats['excluded_by_config'] else "") + ".")
    if _op_stats["unmatched_families"]:
        logger.info(f"Excluded (no sourced swath -- add to MISSION_TABLE with a source to include): "
                    f"{_op_stats['unmatched_families']}; names: {_op_stats['unmatched_names'][:60]}")
    if not constellation:
        raise ConfigurationError(f"No satellite in fleet {FLEET!r} matched a MISSION_TABLE entry.")

real_epochs = [sat["tle"].getDate() for sat in constellation]
# The simulation window starts at the reference epoch: the onset for a
# disaster scenario (so only acquisitions after the event count), the newest
# element set for a live run. Each satellite is propagated from its own
# element-set epoch (at most 7 days earlier) to this instant.
from org.orekit.time import AbsoluteDate, TimeScalesFactory
sim_epoch = AbsoluteDate(_reference_epoch, TimeScalesFactory.getUTC())

# =============================================================================
# 4. COMPUTE DYNAMIC PARAMETERS (Altitude & Orbital Period)
# =============================================================================

mu_earth = 398600.4418
Re_km = 6378.137

pos_list = []
for sat in constellation:
    pv = sat["propagator"].getPVCoordinates(sim_epoch, itrf)
    r = np.sqrt(pv.getPosition().getX()**2 + pv.getPosition().getY()**2 + pv.getPosition().getZ()**2) / 1000.0
    pos_list.append(r - Re_km)

CONSTELLATION_ALTITUDE_KM = float(np.median(pos_list))
T_MAX_WAIT_MINUTES = 2 * np.pi * np.sqrt((Re_km + CONSTELLATION_ALTITUDE_KM)**3 / mu_earth) / 60.0

# =============================================================================
# SUMMARY & PERSISTENCE
# =============================================================================

_summary = {
    "n_satellites": len(constellation),
    "median_alt_km": CONSTELLATION_ALTITUDE_KM,
    "t_max_wait_min": T_MAX_WAIT_MINUTES,
    "sim_epoch": str(sim_epoch),
    "tle_source": TLE_SOURCE,
    "analysis_epoch": globals().get("ANALYSIS_EPOCH"),
    "reference_epoch": _reference_epoch,
    "reference_epoch_source": REFERENCE_EPOCH_SOURCE,
    "fleet": FLEET,
    "satellite_selection": dict(method=("operational at reference epoch, seeded random sample" if FLEET == "starlink"
                                        else "all operational satellites with a sourced MISSION_TABLE entry"),
                                **_op_stats),
    "mission_table": MISSION_TABLE if FLEET != "starlink" else None,
    "acquisition_exclusions": globals().get("ACQUISITION_EXCLUSIONS") if FLEET != "starlink" else None,   # v1.5.2
}

print("\n" + "=" * 60)
print("  CELL 2: TLE Acquisition & SGP4 Dynamics")
print("=" * 60)
# v1.5.2: real fleets read CelesTrak GP data (OMM JSON) or Space-Track history,
# and use every operational satellite -- "(Pinned)" described the Starlink
# proxy's seeded sample only.
print(f"  Orbital elements:          " + (TLE_SOURCE if FLEET == "starlink" else
                                          ("Space-Track gp_history (as of the analysis epoch)"
                                           if globals().get("ANALYSIS_EPOCH") else "CelesTrak GP data (OMM JSON)")))
print(f"  Window start:              {_reference_epoch} UTC [{REFERENCE_EPOCH_SOURCE}]")
print(f"  Simulation Epoch:          {sim_epoch}")
print(f"  Constellation Size (N):    {len(constellation)} "
      + ("(seeded sample, pinned)" if FLEET == "starlink" else "(every operational satellite with a sourced swath)"))
if _op_stats:
    print(f"  Operational at reference:  {_op_stats['operational']} of {_op_stats['fetched']} fetched "
          + (f"(seeded sample, seed={RANDOM_SEED})" if FLEET == "starlink" else f"(fleet {FLEET}: {_op_stats.get('per_mission')})"))
print(f"  Fleet Median Altitude:     {CONSTELLATION_ALTITUDE_KM:.2f} km")
print(f"  Derived Orbital Period:    {T_MAX_WAIT_MINUTES:.2f} min")
print("=" * 60)

try:
    with open(os.path.join(_out_dir, 'arf_cell2_summary.json'), 'w') as f:
        json.dump(_summary, f, indent=2)
    logger.info("Summary exported to arf_cell2_summary.json")
except Exception as e:
    logger.warning(f"Could not save Cell 2 summary: {e}")

print("\n=== CELL 2 COMPLETE ===")

