"""End-to-end harness for v1.5.x: Cell 0 runs for real (ARF_CONFIG_OVERRIDES);
Cells 1-3 are stood in (target tiles, a synthetic Earth-observation fleet on
circular orbits over a rotating spherical Earth, Orekit stand-ins); Cells 4-9
(not 7D, which queries the online archive) run from the pipeline file with
their `from org.` imports replaced by the stand-ins. OUT_DIR differs from the
working directory, so a stray write is visible."""
import sys, os, io, math, json, tempfile, contextlib, logging, typing, re
import numpy as np, pandas as pd, h3, matplotlib
matplotlib.use("Agg")
Re_km, mu_earth = 6378.137, 398600.4418
WE = 360.0 / 1436.068
T0 = np.datetime64("2026-10-07T06:00:00")


def orbit(sat, t_min):
    h, inc = sat["h"], math.radians(sat["inc"])
    per = 2 * math.pi * math.sqrt((Re_km + h) ** 3 / mu_earth) / 60.0
    u = sat["phase"] + 2 * np.pi * np.asarray(t_min) / per
    lat = np.degrees(np.arcsin(np.sin(inc) * np.sin(u)))
    lon = np.degrees(np.arctan2(np.cos(inc) * np.sin(u), np.cos(u))) - WE * np.asarray(t_min) + sat["raan"]
    return lat, (lon + 180) % 360 - 180


class T:
    def __init__(self, s): self.s = s
    def shiftedBy(self, dt): return T(self.s + dt)
    def toString(self, scale=None): return str(T0 + np.timedelta64(int(round(self.s * 1e6)), "us"))
    def __str__(self): return self.toString()
class Pos:
    def __init__(self, lat, lon, r): self.lat, self.lon, self.r = lat, lon, r
    def getNorm(self): return self.r * 1000.0
class PV:
    def __init__(self, p): self.p = p
    def getPosition(self): return self.p
class Prop:
    def __init__(self, sat): self.sat = sat
    def getPVCoordinates(self, t, frame):
        lat, lon = orbit(self.sat, t.s / 60.0)
        return PV(Pos(float(lat), float(lon), Re_km + self.sat["h"]))
class Geo:
    def __init__(self, p): self.p = p
    def getLatitude(self): return math.radians(self.p.lat)
    def getLongitude(self): return math.radians(self.p.lon)
class Earth:
    def transform(self, pos, frame, t): return Geo(pos)
class FastMath:
    toDegrees = staticmethod(math.degrees); toRadians = staticmethod(math.radians)
class Constants:
    WGS84_EARTH_EQUATORIAL_RADIUS = Re_km * 1000.0
class GeodeticPoint:
    def __init__(self, lat, lon, alt): self.lat, self.lon, self.alt = lat, lon, alt
class TopocentricFrame:
    def __init__(self, earth, gp, name): self.gp = gp
    def getElevation(self, pos, frame, t):
        p1, p2 = self.gp.lat, math.radians(pos.lat)
        dl = math.radians(pos.lon) - self.gp.lon
        g = math.acos(max(-1, min(1, math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl))))
        return math.atan2(math.cos(g) - Re_km / pos.r, math.sin(g))
class TimeScalesFactory:
    @staticmethod
    def getUTC(): return None


def eo_fleet(seed=3, n_dove=40):
    """Sentinel-2-like (786 km, 290 km swath), Landsat-like (705 km, 185 km),
    Dove-like (475 km, 25 km), sun-synchronous-like 98 deg orbits."""
    rng = np.random.default_rng(seed)
    sats = []
    for fam, mission, h, sw, sun, k in (("SENTINEL-2", "Sentinel-2", 786, 290.0, 5.0, 3),
                                         ("LANDSAT", "Landsat 8/9", 705, 185.0, 5.0, 2),
                                         ("FLOCK", "PlanetScope", 475, 25.0, 10.0, n_dove)):
        for i in range(k):
            s = {"id": f"{fam}-{i}", "h": float(h), "inc": 98.0, "phase": float(rng.uniform(0, 2 * np.pi)),
                 "raan": float(rng.uniform(15, 40)), "swath_km": sw, "min_sun_elev_deg": sun, "mission": mission}
            s["propagator"] = Prop(s)
            sats.append(s)
    return sats


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


def strip_org(code):
    return "\n".join(l for l in code.splitlines() if not l.strip().startswith("from org."))


def run(path, tiles, overrides=None, cells=("4", "5", "6", "6B", "7", "7B", "7C", "8", "9"), sats=None,
        quiet=True, extra=None, window_hours=48):
    C = cells_of(open(path, encoding="utf-8").read())
    sats = sats if sats is not None else eo_fleet()
    out_dir, work_dir = tempfile.mkdtemp(prefix="E2E15_OUT_"), tempfile.mkdtemp(prefix="E2E15_CWD_")
    ov = {"active_roi": "lower_zambezi", "fleet": "earth-resources+planet", "planet_swath_case": "optimistic",
          "window_hours": window_hours, "cloud_seeds": 3, "cloud_sweep_lengths_km": (0.0, 50.0),
          "cloud_sweep_times_h": (3.0,), "cloud_sweep_seeds": 2, "onset_shift_step_min": 240,
          "random_baseline_draws": 2, "ilp_instances": 2}
    ov.update(overrides or {})
    ns = {"__name__": "__main__", "ARF_CONFIG_OVERRIDES": ov, "display": lambda *a, **k: None}
    sink = io.StringIO()
    old = os.getcwd()
    alts = [s["h"] for s in sats]
    gs = {"BIUST": {"lat": -22.59, "lng": 27.14, "alt_m": 994.0}, "Kasane": {"lat": -17.80, "lng": 25.15, "alt_m": 930.0},
          "Kebonye": {"lat": -24.84, "lng": 25.84, "alt_m": 994.0}}
    try:
        with contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext():
            os.chdir(out_dir)
            exec(compile(C["0"], "cell0", "exec"), ns)
            os.chdir(work_dir)
            ver = ns["__version__"]
            for name in ("arf_cell1_validation.json", "arf_cell2_summary.json", "arf_cell3_metadata.json",
                         "arf_coverage_map.png", "arf_coverage_map.pdf", "arf_acquisition_validation.csv",
                         "arf_acquisition_validation_summary.json"):
                open(os.path.join(ns["OUT_DIR"], name), "w").write("{}")
            ns.update(np=np, pd=pd, h3=h3, logging=logging, os=os, json=json, math=math,
                      **{n: getattr(typing, n) for n in ("List", "Dict", "Tuple", "Any", "Optional", "Set")},
                      Re_km=Re_km, mu_earth=mu_earth, TopocentricFrame=TopocentricFrame, GeodeticPoint=GeodeticPoint,
                      TimeScalesFactory=TimeScalesFactory, pos_list=alts, CONSTELLATION_ALTITUDE_KM=float(np.median(alts)),
                      T_MAX_WAIT_MINUTES=2 * np.pi * np.sqrt((Re_km + float(np.median(alts))) ** 3 / mu_earth) / 60.0,
                      target_tiles=set(tiles), constellation=sats, sim_epoch=T(0.0), itrf=None, earth=Earth(),
                      FastMath=FastMath, Constants=Constants, GROUND_STATIONS=gs, ROI_DEFINITIONS={
                          ov["active_roi"]: {"short_name": "test ROI", "tle_source": "celestrak_live"}},
                      roi_boundary_report=None, _reference_epoch=str(T0), REFERENCE_EPOCH_SOURCE="synthetic",
                      _ARF_CELLS_RUN={**ns["_ARF_CELLS_RUN"], **{c: ver for c in ("1", "2", "3", "7D")}})
            ns.update(extra or {})
            for c in cells:
                exec(compile(strip_org(C[c]), f"cell{c}", "exec"), ns)
    finally:
        os.chdir(old)
    ns["_stdout"], ns["_out_dir"], ns["_work_dir"] = sink.getvalue(), out_dir, work_dir
    return ns
