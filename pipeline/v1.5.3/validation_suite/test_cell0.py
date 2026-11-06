"""Cell 0 (configuration): every validation rule rejects what it should and
accepts the defaults; overrides; folder tags; UTC handling of epochs; the
mission table. Usage: test_cell0.py PIPELINE"""
import sys, os, io, re, json, contextlib, tempfile
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
C = e2e15.cells_of(open(P, encoding="utf-8").read())
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


def cell0(overrides=None, extra=None):
    d = tempfile.mkdtemp(prefix="c0_")
    ns = {"__name__": "__main__"}
    if overrides is not None:
        ns["ARF_CONFIG_OVERRIDES"] = overrides
    ns.update(extra or {})
    old = os.getcwd()
    os.chdir(d)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(C["0"], "cell0", "exec"), ns)
    finally:
        os.chdir(old)
    return ns, d


ns, d = cell0()
check(ns["__version__"] and ns["_ARF_CELLS_RUN"] == {"0": ns["__version__"]} and os.path.exists(os.path.join(d, "arf_config.json")),
      f"defaults accepted; version {ns['__version__']}; arf_config.json written; cell registry reset")
cfg = json.load(open(os.path.join(d, "arf_config.json")))
check(cfg == ns["config"].to_dict() | {k: list(v) if isinstance(v, tuple) else v for k, v in ns["config"].to_dict().items()}
      or json.loads(ns["config"].to_json()) == cfg, "arf_config.json equals the config used")

bad = {"active_roi": "mars", "active_n": 7, "cloud_probability": 1.5, "max_selection_steps": 0, "random_seed": 4.2,
       "window_hours": 0, "fleet": "landsat", "planet_swath_case": "medium", "assumed_min_sun_elev_deg": 95.0,
       "swath_half_angle_deg": 0.0, "h3_resolution": 16, "downlink_min_elevation_deg": 90.0, "t_max_wait_minutes": 0,
       "latency_decay_rate": -1, "decay_sweep": (), "onset_shift_hours": 48, "onset_shift_step_min": 0,
       "random_baseline_draws": 0, "ilp_instances": 0, "ilp_horizon_hours": 0, "ilp_time_limit_s": 0,
       "coverage_targets_pct": (101.0,), "h1_reachable_coverage_pct": 0, "h1_floor_margin_min": -1,
       "h1_min_onset_pct": 0, "h1_min_evaluable_onsets": 0, "latency_decay_scale": "slow", "cloud_knowledge": "psychic",
       "time_aware_objective": "fast", "cloud_model": "fractal", "cloud_corr_length_km": 0.0, "cloud_corr_time_h": -1,
       "cloud_planning": "guess", "cloud_probability_source": "era5", "cloud_climatology_month": 13, "cloud_seeds": 1,
       "cloud_sweep_lengths_km": (-5.0,), "analysis_epoch": "next tuesday", "reference_epoch_utc": "soon"}
rejected = []
for k, v in bad.items():
    try:
        extra = {"fleet": "earth-resources"} if k in ("analysis_epoch", "reference_epoch_utc") else {}
        cell0({**extra, k: v})
        rejected.append((k, False))
    except Exception as e:
        rejected.append((k, type(e).__name__ == "ConfigurationError"))
fails = [k for k, ok in rejected if not ok]
check(not fails, f"{len(bad)} invalid values each raise ConfigurationError" + (f"; not rejected: {fails}" if fails else ""))
for k, v in (("analysis_epoch", "2026-08-17T00:00:00"), ("reference_epoch_utc", "2026-10-12T08:27:09Z")):
    try:
        cell0({"fleet": "starlink", k: v}); ok = False
    except Exception as e:
        ok = type(e).__name__ == "ConfigurationError"
    check(ok, f"{k} with the Starlink proxy is rejected")
try:
    cell0({"fleet": "earth-resources", "analysis_epoch": "2099-01-01T00:00:00"}); ok = False
except Exception as e:
    ok = type(e).__name__ == "ConfigurationError"
check(ok, "a future analysis_epoch is rejected")
try:
    cell0({"fleet": "earth-resources", "analysis_epoch": "2026-08-17T00:00:00", "reference_epoch_utc": "2026-08-17T00:00:00"}); ok = False
except Exception as e:
    ok = type(e).__name__ == "ConfigurationError"
check(ok, "analysis_epoch and reference_epoch_utc together are rejected")
try:
    cell0({"no_such_field": 1}); ok = False
except Exception as e:
    ok = type(e).__name__ == "ConfigurationError" and "no_such_field" in str(e)
check(ok, "unknown override field is rejected by name")

# UTC handling: an offset is converted to UTC (2026-08-17T02:00+02:00 == 2026-08-17T00:00Z)
ns, _ = cell0({"fleet": "earth-resources", "analysis_epoch": "2026-08-17T02:00:00+02:00"})
check(ns["ANALYSIS_EPOCH"] == "2026-08-17T00:00:00.000000" and ns["FLEET_TAG"].endswith("_from20260817"),
      f"analysis_epoch with +02:00 -> {ns['ANALYSIS_EPOCH']} UTC, folder tag {ns['FLEET_TAG']}")
ns, _ = cell0({"fleet": "earth-resources+planet", "planet_swath_case": "optimistic", "reference_epoch_utc": "2026-10-12T08:27:09Z",
               "window_hours": 384})
check(ns["REFERENCE_EPOCH_UTC"] == "2026-10-12T08:27:09.000000" and ns["FLEET_TAG"] ==
      "earth-resources_planet-optimistic_384h_at20261012T0827",
      f"reference_epoch_utc -> {ns['REFERENCE_EPOCH_UTC']}, folder tag {ns['FLEET_TAG']}")
# lists from JSON overrides become tuples
ns, _ = cell0({"cloud_sweep_times_h": [1.0, 6.0], "coverage_targets_pct": [99.0, 95.0]})
check(ns["CLOUD_SWEEP_TIMES_H"] == (1.0, 6.0) and ns["COVERAGE_TARGETS_PCT"] == (95.0, 99.0) and
      ns["ARF_CONFIG_OVERRIDES_APPLIED"] == {"cloud_sweep_times_h": [1.0, 6.0], "coverage_targets_pct": [99.0, 95.0]},
      "list overrides become tuples (targets sorted); the applied overrides are recorded")
# mission table and planet case
for case, fov in (("conservative", 3.0), ("optimistic", 4.0)):
    ns, _ = cell0({"fleet": "planet", "planet_swath_case": case})
    ps = [m for m in ns["MISSION_TABLE"] if m["mission"] == "PlanetScope"][0]
    s2 = [m for m in ns["MISSION_TABLE"] if m["mission"] == "Sentinel-2"][0]
    ls = [m for m in ns["MISSION_TABLE"] if m["mission"] == "Landsat 8/9"][0]
    check(ps["fov_deg"]["value"] == fov and s2["swath_km"]["value"] == 290.0 and ls["swath_km"]["value"] == 185.0 and
          ps["min_sun_elev_deg"]["value"] == 10.0 and ls["min_sun_elev_deg"]["value"] == 5.0 and
          s2["min_sun_elev_deg"]["status"] == "assumption" and ns["FLEET_TAG"] == f"planet-{case}",
          f"mission table ({case}): PlanetScope FOV {fov} deg, S2 290 km (sun 5 deg ASSUMPTION), Landsat 185 km "
          f"(sun 5 deg sourced), Planet sun 10 deg; tag {ns['FLEET_TAG']}")
import re as _re
names_ok = all(_re.match(ns["MISSION_TABLE"][i]["name_pattern"], n, _re.I) for i, n in
               ((0, "SENTINEL-2A"), (0, "SENTINEL-2C"), (1, "LANDSAT 8"), (1, "LANDSAT 9"), (2, "FLOCK 4G-9")))
names_no = not any(_re.match(m["name_pattern"], n, _re.I) for m in ns["MISSION_TABLE"]
                   for n in ("SENTINEL-1A", "SENTINEL-3A", "LANDSAT 7", "SKYSAT-C1", "PELICAN-2"))
check(names_ok and names_no, "name patterns: Sentinel-2A-C, Landsat 8/9 and FLOCK match; Sentinel-1/3, Landsat 7, SkySat, Pelican do not")
# overrides do not persist into a later run in the same namespace unless set again
ns = {"__name__": "__main__", "ARF_CONFIG_OVERRIDES": {"fleet": "planet"}}
d = tempfile.mkdtemp(); old = os.getcwd(); os.chdir(d)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(C["0"], "cell0", "exec"), ns)
    ns.pop("ARF_CONFIG_OVERRIDES")
    exec(compile(C["0"], "cell0", "exec"), ns)
os.chdir(old)
_v = tuple(int(x) for x in ns["__version__"].split("."))
_default_fleet = "earth-resources+planet" if _v >= (1, 5, 2) else "starlink"
check(ns["FLEET"] == _default_fleet and ns["ARF_CONFIG_OVERRIDES_APPLIED"] == {},
      f"Cell 0 re-run without overrides returns to the cell's own settings (fleet {ns['FLEET']}; overrides applied: none)")
if _v >= (1, 5, 2):
    # v1.5.2: default scenario, exclude_satellites, scenario key, acquisition exclusions
    ns, _ = cell0()
    check(ns["ACTIVE_ROI"] == "ckgr_okavango" and ns["FLEET"] == "earth-resources+planet" and
          ns["PLANET_SWATH_CASE"] == "conservative" and ns["FLEET_TAG"] == "earth-resources_planet-conservative"
          and ns["SCENARIO_N"] == ns["FLEET_TAG"] and ns["EXCLUDE_SATELLITES"] == (),
          f"v1.5.2 default: {ns['ACTIVE_ROI']}, {ns['FLEET']}, {ns['PLANET_SWATH_CASE']}; scenario key = fleet tag "
          f"{ns['SCENARIO_N']}")
    ns, _ = cell0({"fleet": "starlink", "active_n": 20, "active_roi": "botswana_full"})
    check(ns["SCENARIO_N"] == 20 and ns["FLEET_TAG"] == "N20", "Starlink proxy keeps (ROI, N) as its scenario key")
    ns, _ = cell0({"exclude_satellites": ["sentinel-2a", "LANDSAT 8"]})
    check(ns["EXCLUDE_SATELLITES"] == ("SENTINEL-2A", "LANDSAT 8") and
          ns["FLEET_TAG"] == "earth-resources_planet-conservative_without-LANDSAT8-SENTINEL2A",
          f"exclude_satellites from a JSON list: names upper-cased, folder tag {ns['FLEET_TAG']}")
    bad_ex = []
    for v in ("SENTINEL-2A", ("",), (3,)):
        try:
            cell0({"exclude_satellites": v}); bad_ex.append(False)
        except Exception as e:
            bad_ex.append(type(e).__name__ == "ConfigurationError")
    try:
        cell0({"fleet": "starlink", "exclude_satellites": ("STARLINK-1",)}); bad_ex.append(False)
    except Exception as e:
        bad_ex.append(type(e).__name__ == "ConfigurationError")
    check(all(bad_ex), "exclude_satellites: a bare string, an empty name, a non-string, or the Starlink proxy is rejected")
    from datetime import datetime as _dt
    ex = ns["ACQUISITION_EXCLUSIONS"]
    ok = len(ex) >= 1
    for e in ex:
        ok &= (_dt.fromisoformat(e["start_utc"]) < _dt.fromisoformat(e["end_utc"]) and e["status"] == "sourced"
               and e["source_url"].startswith("https://") and len(e["quote"]) > 20)
    s2a = ex[0]
    ok &= (bool(_re.match(s2a["name_pattern"], "SENTINEL-2A", _re.I)) and not _re.match(s2a["name_pattern"], "SENTINEL-2B", _re.I)
           and not _re.match(s2a["name_pattern"], "SENTINEL-2C", _re.I) and s2a["start_utc"] == "2026-01-13T00:00:00"
           and s2a["end_utc"] == "2026-03-09T00:00:00")
    check(ok, "acquisition exclusions: sourced entries with ISO times; Sentinel-2A tandem pause 13 Jan - 8 Mar 2026 "
              "(end inclusive), pattern matches 2A only")
n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell0.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
