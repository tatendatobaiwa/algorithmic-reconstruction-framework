"""Cell 2 (orbital elements and imaging fleet), with real Orekit.
Usage: test_cell2.py PIPELINE
A. omm_to_tle: every CelesTrak OMM record of the 12 October 2026 cache (both
   groups) is propagated with Orekit's TLEPropagator from the TLE that
   omm_to_tle builds, and with the independent `sgp4` package (Vallado's
   SGP4, initialised from the same OMM record): TEME positions agree to
   centimetres over 0, 1, 12 and 24 h (unit conversions of mean motion, its
   derivatives, angles and BSTAR).
B. build_eo_constellation on synthetic records: stale, future, low-perigee,
   unmatched, malformed and duplicate-name records are counted, never used;
   exclude_names (v1.5.2) leaves the named satellites out and reports names
   not found; the real cache gives the fleet of the Colab run (97 = 2 Landsat
   + 3 Sentinel-2 + 92 Doves at the run's reference epoch).
C. swath_from_fov_km: PlanetScope 3.0 / 4.0 deg at 475 km against the
   PSScene figures quoted in MISSION_TABLE (24.9 km; ~33 km) and an
   independent ray-sphere intersection.
D. load_eo_gp_history's selection: never an element set newer than the epoch,
   the newest one at or before it, per satellite (synthetic cache file)."""
import sys, os, ast, json, math, copy, tempfile, logging, re
import numpy as np
sys.path.insert(0, "/home/claude/validate")
sys.path.insert(0, "/home/claude")
import e2e15
import orekit_env                                   # JVM + Orekit data (/home/claude/orekit-data)
from org.orekit.frames import FramesFactory
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from sgp4.api import Satrec, jday
from sgp4 import omm as sgp4_omm

P = sys.argv[1]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


C = e2e15.cells_of(open(P, encoding="utf-8").read())
src0 = C["0"]
ns0 = {"__name__": "__main__", "ARF_CONFIG_OVERRIDES": {"fleet": "earth-resources+planet", "planet_swath_case": "optimistic"}}
import io, contextlib
d0 = tempfile.mkdtemp()
_old = os.getcwd(); os.chdir(d0)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src0, "cell0", "exec"), ns0)
os.chdir(_old)
tree = ast.parse(C["2"])
tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef)]


class ConfigurationError(Exception):
    pass


ns = {"os": os, "json": json, "np": np, "re": re, "time": __import__("time"), "logger": logging.getLogger("t2"),
      "MISSION_TABLE": ns0["MISSION_TABLE"], "ConfigurationError": ConfigurationError, "TLEPropagator": TLEPropagator,
      "FLEET_TAG": ns0["FLEET_TAG"], "List": list, "Tuple": tuple, "Optional": object, "Any": object, "Dict": dict,
      "CELESTRAK_URL": "unused", "TLE_RAW_CACHE_PATH": "unused", "TLE_CACHE_MAX_AGE_HOURS": 6,
      "OPERATIONAL_CRITERIA": {"min_perigee_km": 300.0, "max_apogee_km": 700.0, "max_elset_age_days": 7.0}}
exec(compile(tree, "cell2_defs", "exec"), ns)
utc = TimeScalesFactory.getUTC()
teme = FramesFactory.getTEME()
RD = "/home/claude/realdata/"
recs = []
for g in ("resource", "planet"):
    recs += [dict(r, _GROUP=g) for r in json.load(open(RD + f"ckgr_openPO_arf_gp_{g}.json"))["records"]]

# A
worst, n = 0.0, 0
for r in recs:
    tle = ns["omm_to_tle"](r, utc)
    prop = TLEPropagator.selectExtrapolator(tle)
    sat = Satrec()
    sgp4_omm.initialize(sat, r)
    ep = tle.getDate()
    for dt_min in (0.0, 60.0, 720.0, 1440.0):
        t = ep.shiftedBy(dt_min * 60.0)
        p = prop.getPVCoordinates(t, teme).getPosition()
        # sgp4 time: minutes since its own epoch (no UTC/TT ambiguity)
        e, rr, _ = sat.sgp4_tsince(dt_min)
        if e != 0:
            continue
        d = math.dist((p.getX() / 1000.0, p.getY() / 1000.0, p.getZ() / 1000.0), rr)
        worst = max(worst, d)
        n += 1
check(n > 1000 and worst < 0.01, f"omm_to_tle + Orekit SGP4 vs the sgp4 package: {len(recs)} records x 4 times "
                                 f"({n} states), worst TEME position difference {worst * 1e6:.3f} mm")
# A2: SGP4 itself ignores n-dot and n-ddot, so A cannot check their units. Build
# TLE text from each OMM record (sgp4's exporter, an independent writer of the
# TLE format), parse it with Orekit's own TLE parser, and compare every element
# with omm_to_tle's TLE (the derivative fields hold n-dot/2 and n-ddot/6).
# TLE text keeps 5 significant digits of BSTAR (OMM JSON keeps more), so BSTAR
# is compared to a relative 1e-4; the other elements to their text precision.
from sgp4.exporter import export_tle
bad = []
for r in recs:
    sat = Satrec()
    sgp4_omm.initialize(sat, r)
    l1, l2 = export_tle(sat)
    if int(r["NORAD_CAT_ID"]) > 99999:
        continue                                    # no TLE text exists for 6-digit catalogue numbers
    a, b = TLE(l1, l2), ns["omm_to_tle"](r, utc)
    for name, fa, fb, tol in (("n", a.getMeanMotion(), b.getMeanMotion(), 1e-12),
                              ("ndot", a.getMeanMotionFirstDerivative(), b.getMeanMotionFirstDerivative(), 1e-16),
                              ("nddot", a.getMeanMotionSecondDerivative(), b.getMeanMotionSecondDerivative(), 1e-20),
                              ("e", a.getE(), b.getE(), 1e-7), ("i", a.getI(), b.getI(), 1e-6),
                              ("raan", a.getRaan(), b.getRaan(), 1e-6), ("argp", a.getPerigeeArgument(), b.getPerigeeArgument(), 1e-6),
                              ("M", a.getMeanAnomaly(), b.getMeanAnomaly(), 1e-6), ("bstar", a.getBStar() / max(abs(b.getBStar()), 1e-12), b.getBStar() / max(abs(b.getBStar()), 1e-12), 1e-4),
                              ("epoch", 0.0, a.getDate().durationFrom(b.getDate()), 1e-3)):
        if abs(fa - fb) > tol * max(1.0, abs(fa)):
            bad.append((r["OBJECT_NAME"], name, fa, fb))
check(not bad, f"omm_to_tle vs Orekit's parser of TLE text written from the same OMM record: every element equal, "
               f"including n-dot and n-ddot ({len(recs)} records)" + (f"; differences {bad[:4]}" if bad else ""))

# B
ref = "2026-10-12T08:27:09.082080"
cons, st = ns["build_eo_constellation"](recs, ref)
pm = st["per_mission"]
check(len(cons) == 97 and pm == {"Landsat 8/9": 2, "Sentinel-2": 3, "PlanetScope": 92},
      f"real cache at the Colab reference epoch: {len(cons)} satellites {pm} (Colab run: 97 = 2 + 3 + 92)")
cons2, st2 = ns["build_eo_constellation"](recs, ref, exclude_names=("sentinel-2a", "NOT-A-SAT"))
check(len(cons2) == 96 and not any(s["id"] == "SENTINEL-2A" for s in cons2) and st2["excluded_by_config"] == 1
      and st2["exclusions_not_found"] == ["NOT-A-SAT"],
      "exclude_names: SENTINEL-2A left out (case-insensitive), a name not in the fleet reported")
base = next(r for r in recs if r["OBJECT_NAME"] == "LANDSAT 9")
synth = []
from datetime import datetime, timedelta


def mk(name, nid, epoch, **kw):
    r = copy.deepcopy(base)
    r.update(OBJECT_NAME=name, NORAD_CAT_ID=nid, EPOCH=epoch, **kw)
    return r


refdt = datetime.fromisoformat(ref)
iso = lambda d: d.strftime("%Y-%m-%dT%H:%M:%S.%f")
synth += [mk("LANDSAT 9", 90001, iso(refdt - timedelta(days=1))),                      # ok
          mk("LANDSAT 9", 90002, iso(refdt - timedelta(days=2))),                      # ok, duplicate name
          mk("LANDSAT 8", 90003, iso(refdt - timedelta(days=8))),                      # stale
          mk("LANDSAT 8", 90004, iso(refdt + timedelta(hours=3))),                     # future
          mk("SENTINEL-2B", 90005, iso(refdt - timedelta(days=1)), MEAN_MOTION=16.4),  # low perigee (~260 km)
          mk("SKYSAT-C1", 90006, iso(refdt - timedelta(days=1))),                      # unmatched
          mk("FLOCK 4X-1", 90007, iso(refdt - timedelta(days=1)), ECCENTRICITY="bad")] # malformed
cons3, st3 = ns["build_eo_constellation"](synth, ref)
ids = sorted(s["id"] for s in cons3)
check(ids == ["LANDSAT 9", "LANDSAT 9 [90002]"] and st3["stale_or_future"] == 2 and st3["low_perigee"] == 1
      and st3["unmatched_mission"] == 1 and st3["malformed"] == 1,
      f"synthetic records: kept {ids}; stale/future {st3['stale_or_future']}, low perigee {st3['low_perigee']}, "
      f"unmatched {st3['unmatched_mission']}, malformed {st3['malformed']}")

# C
sw3 = ns["swath_from_fov_km"](475.0, 3.0)
sw4 = ns["swath_from_fov_km"](475.0, 4.0)


def indep(h, fov, re_=6378.137):
    """Ray from the satellite at fov/2 off nadir, intersected with the sphere
    (vector geometry, no closed-form formula); swath = 2 x Re x central angle."""
    eta = math.radians(fov / 2)
    sx, sy = 0.0, re_ + h
    dx, dy = math.sin(eta), -math.cos(eta)
    b = 2 * (sx * dx + sy * dy)
    c = sx * sx + sy * sy - re_ * re_
    t = (-b - math.sqrt(b * b - 4 * c)) / 2                      # nearer intersection
    ex, ey = sx + t * dx, sy + t * dy
    return 2 * re_ * math.atan2(ex, ey)


check(abs(sw3 - 24.9) < 0.1 and abs(sw4 - indep(475.0, 4.0)) < 1e-9 and abs(sw3 - indep(475.0, 3.0)) < 1e-9
      and 32.0 < sw4 < 34.5,
      f"FOV swath at 475 km: 3.0 deg -> {sw3:.2f} km (spec 24.9 km), 4.0 deg -> {sw4:.2f} km; matches a "
      f"ray-sphere intersection")

# D
tmp = tempfile.mkdtemp()
epoch_iso = "2026-08-17T00:00:00"
end = datetime.fromisoformat(epoch_iso)
live = [dict(base, OBJECT_NAME="LANDSAT 9", NORAD_CAT_ID=49260, _GROUP="resource")]
hist = [mk("LANDSAT 9", 49260, iso(end - timedelta(days=3))), mk("LANDSAT 9", 49260, iso(end - timedelta(hours=5))),
        mk("LANDSAT 9", 49260, iso(end + timedelta(hours=1)))]                  # the last is AFTER the epoch
path = os.path.join(tmp, f"arf_gp_history_{ns0['FLEET_TAG']}_{end:%Y%m%dT%H%M%S}.json")
json.dump({"fetched_at": 0, "epoch": epoch_iso, "norad_ids": [49260], "records": hist}, open(path, "w"))
out = ns["load_eo_gp_history"](live, epoch_iso, tmp)
check(len(out) == 1 and out[0]["EPOCH"] == hist[1]["EPOCH"],
      "retrospective elements: the newest set at or before the epoch is used, never a later one")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell2.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
