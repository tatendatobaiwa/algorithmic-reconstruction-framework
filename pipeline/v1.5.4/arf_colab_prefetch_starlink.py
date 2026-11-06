# =============================================================================
# ARF 1.5.4 PREFETCH CELL (Colab): Space-Track Starlink histories for every
# window of the dense-constellation arm, fetched once while you are at the tab
# =============================================================================
# Colab hands out the Space-Track secrets only while the notebook tab is active,
# so a window whose first run starts unattended fails ("Could not read
# Space-Track credentials"). This cell fetches the history of each protocol
# window that is not yet cached, with the pipeline's own function
# (load_starlink_history_at of release 1.5.4, taken from the release file, so
# the cache files are exactly what Cell 2 would write), into the batch's shared
# cache. After it, every run of the arm reads the cache and needs no secrets.
# One query per window (at most 8), within Space-Track's published limits.
# It also counts the element sets with Alpha-5 catalogue numbers, which the
# release's operational filter leaves out (reported, not changed).
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
CACHE = f"{ROOT}/ARF_batch/shared_cache"
SHA_154 = "2b618ec442e71cf896e1104c6da89c94f2d6b6785a04c4a9236b92c21115a1cb"

import os, ast, json, time, hashlib, logging, runpy, urllib.request, urllib.parse

if ROOT.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
    from google.colab import drive
    drive.mount("/content/drive")

pipe = os.path.join(CODE, "algorithmic_reconstruction_framework_v154.py")
src = open(pipe, encoding="utf-8").read()
assert hashlib.sha256(src.encode("utf-8")).hexdigest() == SHA_154, "ARF_code does not hold release 1.5.4"
# the four functions the fetch uses, exactly as in the release file
wanted = ("_get_spacetrack_credentials", "_spacetrack_login", "_dedupe_gp_history_to_latest_per_satellite",
          "_parse_gp_history_epoch", "load_starlink_history_at")
tree = ast.parse(src)
defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
assert sorted(n.name for n in defs) == sorted(wanted), f"functions not found: {set(wanted) - {n.name for n in defs}}"
ns = {"os": os, "json": json, "time": time, "urllib": urllib, "logger": logging.getLogger("arf.prefetch")}
exec(compile(ast.Module(body=defs, type_ignores=[]), pipe, "exec"), ns)
logging.basicConfig(level=logging.INFO)

windows = runpy.run_path(os.path.join(CODE, "arf_batch_runner_v154.py"), run_name="arf_batch_runner_v154")["PROTOCOL_WINDOWS"]
os.makedirs(CACHE, exist_ok=True)
for w in windows:
    path = os.path.join(CACHE, f"arf_gp_history_starlink_{w.replace('-', '').replace(':', '')}_7d.json")
    was_cached = os.path.exists(path)
    recs, p = ns["load_starlink_history_at"](w, CACHE)
    assert p == path, (p, path)
    alpha5 = sum(1 for _, l1, _ in recs if not str(l1)[2:7].strip().isdigit())
    print(f"{w}: {'cached already' if was_cached else 'fetched now'}, {len(recs)} satellites, "
          f"{alpha5} with Alpha-5 catalogue numbers ({100.0 * alpha5 / max(len(recs), 1):.1f}%)")
print("All 8 windows cached: the arm's runs no longer need the Space-Track secrets.")
