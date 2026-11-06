# =============================================================================
# ARF A5 PREFETCH (Colab): Space-Track element sets for the eight wet-season
# windows, fetched once while the notebook tab is active (secrets available)
# =============================================================================
# The A5 runs execute the pipeline with google.colab hidden (so each run writes
# to its own folder), which also hides the Colab secrets; they must therefore
# find their element sets already cached. This cell fetches them with the
# release's OWN functions (load_gp_json, load_eo_gp_history and their helpers,
# taken from algorithmic_reconstruction_framework.py of release 1.5.3 by AST,
# unchanged), into ARF_batch/shared_cache, under the exact file names Cell 2
# reads (arf_gp_history_<fleet tag>_<window>.json). Windows already cached are
# not fetched again. One gp_history query per 100 satellites per window.
# Secrets are read by the release's own code; nothing is printed of them.
#
# Use: run in a cell (not !python), with the Space-Track secrets set.
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
CACHE = f"{ROOT}/ARF_batch/shared_cache"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
WINDOWS = ("2026-01-05T00:00:00", "2026-01-12T00:00:00", "2026-01-19T00:00:00", "2026-01-26T00:00:00",
           "2026-02-02T00:00:00", "2026-02-09T00:00:00", "2026-02-16T00:00:00", "2026-02-23T00:00:00")
DESIGN = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative"}

import os, re, ast, json, time, hashlib, logging, urllib.request, urllib.parse
from datetime import datetime

FUNCS = ("_get_spacetrack_credentials", "_spacetrack_login", "load_gp_json", "load_eo_gp_history", "_omm_epoch",
         "_match_mission")
ASSIGNS = ("CELESTRAK_GROUPS", "_ACCESSED", "MISSION_TABLE", "TLE_CACHE_MAX_AGE_HOURS")


def fleet_tag(cfg):
    """Cell 0's FLEET_TAG of a past-window real-fleet run (as in arf_colab_rebuild_20261026.py)."""
    fl = cfg["fleet"]
    t = fl.replace("+", "_") + (f"-{cfg.get('planet_swath_case', 'conservative')}" if "planet" in fl else "")
    return t + "_from" + cfg["analysis_epoch"][:10].replace("-", "")


def release_namespace(pipe, sha_expected=SHA_153):
    """The release's fetch functions and the constants they read, by AST, unchanged."""
    src = open(pipe, encoding="utf-8").read()
    assert hashlib.sha256(src.encode("utf-8")).hexdigest() == sha_expected, f"{pipe} is not release 1.5.3"
    tree = ast.parse(src)
    keep = []
    for n in tree.body:
        if isinstance(n, ast.FunctionDef) and n.name in FUNCS:
            keep.append(n)
        elif isinstance(n, ast.ClassDef) and n.name == "ConfigurationError":
            keep.append(n)
        elif isinstance(n, ast.Assign) and any(getattr(t, "id", None) in ASSIGNS for t in n.targets):
            keep.append(n)
    got = {n.name for n in keep if hasattr(n, "name")} | {t.id for n in keep if isinstance(n, ast.Assign)
                                                           for t in n.targets}
    assert set(FUNCS) | set(ASSIGNS) | {"ConfigurationError"} <= got, f"missing in the release: {set(FUNCS) - got}"
    # MISSION_TABLE reads two Cell 0 settings that do not affect the name patterns used for matching;
    # they are given the primary design's values.
    ns = {"os": os, "re": re, "json": json, "time": time, "urllib": urllib, "logger": logging.getLogger("arf.a5"),
          "ASSUMED_MIN_SUN_ELEV_DEG": 5.0, "PLANET_SWATH_CASE": DESIGN["planet_swath_case"]}
    exec(compile(ast.Module(body=keep, type_ignores=[]), pipe, "exec"), ns)
    return ns


def history_path(cache, cfg):
    epoch = datetime.fromisoformat(cfg["analysis_epoch"])
    return os.path.join(cache, f"arf_gp_history_{fleet_tag(cfg)}_{epoch:%Y%m%dT%H%M%S}.json")


def prefetch(code=CODE, cache=CACHE, windows=WINDOWS):
    ns = release_namespace(os.path.join(code, "algorithmic_reconstruction_framework.py"))
    os.makedirs(cache, exist_ok=True)
    out = []
    for w in windows:
        cfg = {**DESIGN, "analysis_epoch": w}
        path = history_path(cache, cfg)
        was = os.path.exists(path)
        ns["FLEET_TAG"] = fleet_tag(cfg)
        live = ns["load_gp_json"](ns["CELESTRAK_GROUPS"][cfg["fleet"]], cache)
        epoch = datetime.fromisoformat(w).strftime("%Y-%m-%dT%H:%M:%S.%f")       # as Cell 0 formats ANALYSIS_EPOCH
        recs = ns["load_eo_gp_history"](live, epoch, cache)
        assert os.path.exists(path), f"expected cache file {path} was not written"
        out.append((w[:10], "cached already" if was else "fetched now", len(recs)))
        print(f"{w[:10]}: {out[-1][1]}, {len(recs)} satellites with an element set in the 7 days before the window")
    return out


if __name__ == "__main__" or "get_ipython" in globals():
    if not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    logging.basicConfig(level=logging.INFO)
    prefetch()
    print("All eight wet-season windows cached; the A5 runs need no secrets.")
