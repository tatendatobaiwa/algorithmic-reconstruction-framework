"""Test of arf_colab_prefetch_wet_a5.py with the REAL release 1.5.3 functions (taken by AST) and the REAL
CelesTrak lists saved on 12 October 2026 (made fresh); Space-Track is replaced by a stand-in that returns one
element set per requested satellite, dated one day before the window.
  - file names: identical to the names Cell 2 wrote for the confirmatory windows on Drive (checked against the
    real shared_cache listing of 26 October for 2026-08-03 and 2026-09-21)
  - the satellites requested are exactly those of the current lists whose names match a MISSION_TABLE entry
  - a second call uses the cache (no Space-Track call); a release file other than 1.5.3 is refused
Usage: test_prefetch_a5.py PREFETCH PIPELINE_V153 GP_RESOURCE GP_PLANET DRIVE_LISTING"""
import sys, os, json, time, re, tempfile, shutil, importlib.util, io, contextlib
from datetime import datetime, timedelta

PF, P, GR, GPL, LIST = sys.argv[1:6]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


spec = importlib.util.spec_from_file_location("pf", PF); pf = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(pf)
listing = set(open(LIST).read().split())
for w in ("2026-08-03T00:00:00", "2026-09-21T00:00:00"):
    n = os.path.basename(pf.history_path("/x", {**pf.DESIGN, "analysis_epoch": w}))
    check(n in listing, f"name for {w[:10]} equals the real Drive cache file: {n}")

root = tempfile.mkdtemp(prefix="pfa5_")
code, cache = os.path.join(root, "ARF_code"), os.path.join(root, "cache")
os.makedirs(code); os.makedirs(cache)
shutil.copy(P, os.path.join(code, "algorithmic_reconstruction_framework.py"))
for src, g in ((GR, "resource"), (GPL, "planet")):
    d = json.load(open(src)); d["fetched_at"] = time.time()
    json.dump(d, open(os.path.join(cache, f"arf_gp_{g}.json"), "w"))
calls = []


class Resp:
    def __init__(self, b): self.b = b
    def read(self): return self.b
    def getcode(self): return 200
    def __enter__(self): return self
    def __exit__(self, *a): return False


class Opener:
    def open(self, req, timeout=None):
        url = req.full_url
        calls.append(url)
        ids = [int(x) for x in re.search(r"NORAD_CAT_ID/([0-9,]+)/", url).group(1).split(",")]
        end = re.search(r"--([0-9T:\-]+)/orderby", url).group(1)
        ep = (datetime.fromisoformat(end) - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S")
        return Resp(json.dumps([{"NORAD_CAT_ID": str(i), "OBJECT_NAME": "X", "EPOCH": ep, "MEAN_MOTION": "15",
                                 "ECCENTRICITY": "0.001"} for i in ids]).encode())


real_ns = pf.release_namespace


def patched(pipe, sha_expected=pf.SHA_153):
    ns = real_ns(pipe, sha_expected)
    ns["_get_spacetrack_credentials"] = lambda: ("user", "pw")      # stand-in; the real one reads Colab secrets
    ns["_spacetrack_login"] = lambda u, p, base_url=None: Opener()
    return ns


pf.release_namespace = patched
with contextlib.redirect_stdout(io.StringIO()):
    out = pf.prefetch(code=code, cache=cache, windows=pf.WINDOWS[:2])
check([o[1] for o in out] == ["fetched now", "fetched now"] and len(calls) >= 2, f"two windows fetched ({len(calls)} queries)")
ns = real_ns(os.path.join(code, "algorithmic_reconstruction_framework.py"))
live = ns["load_gp_json"](ns["CELESTRAK_GROUPS"]["earth-resources+planet"], cache)
want = sorted(int(r["NORAD_CAT_ID"]) for r in live if ns["_match_mission"](str(r.get("OBJECT_NAME", "")).strip()))
got = json.load(open(pf.history_path(cache, {**pf.DESIGN, "analysis_epoch": pf.WINDOWS[0]})))["norad_ids"]
check(got == want and len(want) > 50, f"requested satellites = those matching MISSION_TABLE ({len(want)})")
n0 = len(calls)
with contextlib.redirect_stdout(io.StringIO()):
    out2 = pf.prefetch(code=code, cache=cache, windows=pf.WINDOWS[:2])
check(len(calls) == n0 and [o[1] for o in out2] == ["cached already"] * 2, "second call: from the cache, no query")
open(os.path.join(code, "algorithmic_reconstruction_framework.py"), "a").write("\n")
try:
    real_ns(os.path.join(code, "algorithmic_reconstruction_framework.py")); ok = False
except AssertionError as e:
    ok = "not release 1.5.3" in str(e)
check(ok, "a changed release file is refused")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
