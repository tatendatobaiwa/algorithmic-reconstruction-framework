"""Build release 1.5.5 from release 1.5.4 by exact, counted text replacements.

1.5.5 serves the Paper 2 swath correction (proposed amendment A14); the real
Earth-observation fleets and the 1.5.4 proxy default are untouched:
  Cell 0  new field proxy_swath_km (None = the off-nadir half-angle
          swath_half_angle_deg, 9.22 deg, as in every earlier proxy run; a value =
          a ground swath in km given to every proxy satellite, which is then
          treated exactly as a real-fleet satellite with a sourced swath: its
          swath-edge angle at its own altitude, distance to the ground track
          rather than to the nearest sample); folder tag _sw<km>km when set;
          console line.
  Cell 2  the swath attached to each proxy satellite when set; recorded in the
          Cell 2 summary.
  Cell 5  log and console lines name the per-satellite swath when set.
  All     version strings.
No line computes anything differently when proxy_swath_km is None.

Usage: make_v155.py IN_1.5.4 OUT_1.5.5
"""
import sys
import hashlib

SRC, DST = sys.argv[1], sys.argv[2]
src = open(SRC, encoding="utf-8").read()
EXPECTED_SHA_IN = "2b618ec442e71cf896e1104c6da89c94f2d6b6785a04c4a9236b92c21115a1cb"
assert hashlib.sha256(src.encode("utf-8")).hexdigest() == EXPECTED_SHA_IN, "input is not the shipped 1.5.4"


def rep(old, new, count=1):
    global src
    n = src.count(old)
    assert n == count, f"expected {count} occurrence(s), found {n}: {old[:90]!r}"
    src = src.replace(old, new)


# ----------------------------------------------------------------------------- version strings
rep('__version__ = "1.5.4"   # 1.5.4: dense-constellation arm',
    '__version__ = "1.5.5"   # 1.5.5: Starlink proxy ground swath in km (proxy_swath_km; None = the 1.5.4 '
    'half-angle, unchanged), real fleets unchanged (Cells 0, 2, 5); 1.5.4: dense-constellation arm')
rep('_ARF_CELL_VERSION = "1.5.4"', '_ARF_CELL_VERSION = "1.5.5"', count=13)

# ----------------------------------------------------------------------------- Cell 0: field
rep("    proxy_min_sun_elev_deg: Optional[float] = None\n\n",
    "    proxy_min_sun_elev_deg: Optional[float] = None\n"
    "    # Starlink proxy only (v1.5.5): None = the off-nadir half-angle\n"
    "    # swath_half_angle_deg below (9.22 deg, about 155-180 km on the ground at\n"
    "    # 475-560 km; every proxy run before 1.5.5, value undocumented); a value =\n"
    "    # the GROUND SWATH in km of every proxy satellite, applied as for a real\n"
    "    # satellite with a sourced swath (edge angle at its own altitude, exact\n"
    "    # distance to the ground track). An ASSUMPTION (Paper 2 swath sensitivity).\n"
    "    proxy_swath_km: Optional[float] = None\n\n")
rep("            if not -10.0 <= float(self.proxy_min_sun_elev_deg) < 90.0:\n"
    "                raise ConfigurationError(f\"PROXY_MIN_SUN_ELEV_DEG={self.proxy_min_sun_elev_deg} must be in [-10, 90).\")\n",
    "            if not -10.0 <= float(self.proxy_min_sun_elev_deg) < 90.0:\n"
    "                raise ConfigurationError(f\"PROXY_MIN_SUN_ELEV_DEG={self.proxy_min_sun_elev_deg} must be in [-10, 90).\")\n"
    "        if self.proxy_swath_km is not None:                       # v1.5.5\n"
    "            if self.fleet != \"starlink\":\n"
    "                raise ConfigurationError(\"PROXY_SWATH_KM applies to the Starlink proxy only (real fleets take \"\n"
    "                                         \"theirs from MISSION_TABLE).\")\n"
    "            if not 0.5 <= float(self.proxy_swath_km) <= 400.0:\n"
    "                raise ConfigurationError(f\"PROXY_SWATH_KM={self.proxy_swath_km} must be in [0.5, 400] km.\")\n")
rep("    proxy_min_sun_elev_deg=None,   # Starlink proxy daylight limit (v1.5.4); None = none, as before\n",
    "    proxy_min_sun_elev_deg=None,   # Starlink proxy daylight limit (v1.5.4); None = none, as before\n"
    "    proxy_swath_km=None,           # Starlink proxy ground swath, km (v1.5.5); None = 9.22 deg half-angle, as before\n")
rep("                          else float(config.proxy_min_sun_elev_deg))     # v1.5.4\n",
    "                          else float(config.proxy_min_sun_elev_deg))     # v1.5.4\n"
    "PROXY_SWATH_KM = None if config.proxy_swath_km is None else float(config.proxy_swath_km)   # v1.5.5\n")
rep("    FLEET_TAG += f\"_sun{PROXY_MIN_SUN_ELEV_DEG:g}\"\n",
    "    FLEET_TAG += f\"_sun{PROXY_MIN_SUN_ELEV_DEG:g}\"\n"
    "if FLEET == \"starlink\" and PROXY_SWATH_KM is not None:           # proxy ground swath (v1.5.5)\n"
    "    FLEET_TAG += f\"_sw{PROXY_SWATH_KM:g}km\"\n")
rep("if FLEET == \"starlink\":\n    print(f\"  Sensor Half-Angle:         {SWATH_HALF_ANGLE_DEG}°\")\n",
    "if FLEET == \"starlink\":\n"
    "    print(f\"  Sensor Half-Angle:         {SWATH_HALF_ANGLE_DEG}°\" if PROXY_SWATH_KM is None else\n"
    "          f\"  Proxy ground swath:        {PROXY_SWATH_KM:g} km per satellite (assumption; v1.5.5)\")\n")

# ----------------------------------------------------------------------------- Cell 2: attach swath
rep("            if PROXY_MIN_SUN_ELEV_DEG is not None:  # v1.5.4: daylight limit for the proxy (assumption)\n"
    "                constellation[-1][\"min_sun_elev_deg\"] = PROXY_MIN_SUN_ELEV_DEG\n",
    "            if PROXY_MIN_SUN_ELEV_DEG is not None:  # v1.5.4: daylight limit for the proxy (assumption)\n"
    "                constellation[-1][\"min_sun_elev_deg\"] = PROXY_MIN_SUN_ELEV_DEG\n"
    "            if globals().get(\"PROXY_SWATH_KM\") is not None:   # v1.5.5: ground swath (assumption)\n"
    "                constellation[-1][\"swath_km\"] = PROXY_SWATH_KM\n")
rep("    \"acquisition_exclusions\": globals().get(\"ACQUISITION_EXCLUSIONS\") if FLEET != \"starlink\" else None,   # v1.5.2\n}\n",
    "    \"acquisition_exclusions\": globals().get(\"ACQUISITION_EXCLUSIONS\") if FLEET != \"starlink\" else None,   # v1.5.2\n"
    "    \"proxy_swath_km\": globals().get(\"PROXY_SWATH_KM\") if FLEET == \"starlink\" else None,   # v1.5.5\n}\n")

# ----------------------------------------------------------------------------- Cell 5: log and console
rep("_real_fleet5 = globals().get(\"FLEET\", \"starlink\") != \"starlink\"\n",
    "_real_fleet5 = globals().get(\"FLEET\", \"starlink\") != \"starlink\"\n"
    "_proxy_swath5 = (not _real_fleet5) and globals().get(\"PROXY_SWATH_KM\") is not None   # v1.5.5\n")
rep("            + (\"swath edge: each satellite's own (MISSION_TABLE)\" if _real_fleet5 else\n",
    "            + (\"swath edge: each satellite's own (MISSION_TABLE)\" if _real_fleet5 else\n"
    "               f\"Starlink proxy ground swath {globals().get('PROXY_SWATH_KM'):g} km, edge angle per satellite\"\n"
    "               if _proxy_swath5 else\n")
rep("              if _real_fleet5 and _summary[\"theta_max_deg_range\"] else",
    "              if (_real_fleet5 or _proxy_swath5) and _summary[\"theta_max_deg_range\"] else")

out = src
open(DST, "w", encoding="utf-8").write(out)
print(DST, hashlib.sha256(out.encode("utf-8")).hexdigest())
