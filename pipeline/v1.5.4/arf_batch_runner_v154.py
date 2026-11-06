# =============================================================================
# ARF BATCH RUNNER (v1.5.4): the dense-constellation arm (protocol amendment A3)
# =============================================================================
# The batch machinery of arf_batch_runner.py v1.5.3 (fresh namespace per run,
# resumable index, a run marked "done" is skipped only when its folder holds a
# manifest of the same release with the same settings), run on release 1.5.4,
# which is installed beside 1.5.3 as algorithmic_reconstruction_framework_v154.py
# so the real-fleet batch keeps reading 1.5.3.
#
# Runs (27, batch name exploratory_v154_dense): the Starlink proxy standing in
# for a hypothetical dense Earth-observation constellation over CKGR + Okavango,
# daylight limit 5 deg (assumption), cloud-sensitivity sweep off; every other
# setting as in the protocol. Window 1 (3 August 2026): N = 100, 500, 1000,
# 2000, 4000, 8000; windows 2-8: N = 100, 500, 1000. N = 8000 runs last (about
# 12 h on its own, estimated). About 43 h for the 27 runs. Exploratory: no
# hypothesis decision is drawn from it.
#
# Use (Colab): through arf_colab_launch_v154_dense.py (Space-Track secrets
# SPACETRACK_USERNAME and SPACETRACK_PASSWORD, tab kept open).
# =============================================================================
import os
import io
import re
import gc
import json
import time
import hashlib
import traceback
import contextlib
from datetime import datetime, timedelta, timezone

PIPELINE_PATH = "/content/drive/MyDrive/ARF_code/algorithmic_reconstruction_framework_v154.py"
BATCH_DIR = "/content/drive/MyDrive/ARF_batch"
SHARED_CACHE_DIR = os.path.join(BATCH_DIR, "shared_cache")


def windows(first_utc, n, step_days):
    """ISO start times of n analysis windows, step_days apart (UTC)."""
    t0 = datetime.fromisoformat(first_utc.replace("Z", ""))
    return [(t0 + timedelta(days=step_days * i)).strftime("%Y-%m-%dT%H:%M:%S") for i in range(n)]


def run_id(cfg):
    """Readable, stable identifier of a configuration."""
    parts = [str(cfg.get("active_roi", "roi")), str(cfg.get("fleet", "fleet"))]
    if "planet" in str(cfg.get("fleet", "")):
        parts.append(str(cfg.get("planet_swath_case", "conservative")))
    if cfg.get("analysis_epoch"):
        parts.append(str(cfg["analysis_epoch"]).replace(":", "").replace("-", "")[:13])
    rest = {k: v for k, v in sorted(cfg.items()) if k not in ("active_roi", "fleet", "planet_swath_case", "analysis_epoch")}
    if rest:
        parts.append(hashlib.sha256(json.dumps(rest, sort_keys=True, default=str).encode()).hexdigest()[:8])
    return "__".join(parts)


def _normalise(d):
    return json.loads(json.dumps(d, sort_keys=True, default=str))


def _manifest_matches(out_dir, release, cfg):
    p = os.path.join(out_dir or "", "arf_reproducibility_manifest.json")
    if not out_dir or not os.path.exists(p):
        return False
    try:
        m = json.load(open(p))
    except Exception:
        return False
    return (m.get("framework_version") == release and
            _normalise(m.get("config_overrides") or {}) == _normalise({k: (list(v) if isinstance(v, tuple) else v)
                                                                       for k, v in cfg.items()}))


def _load_index(path):
    if not os.path.exists(path):
        return {}
    import pandas as pd
    df = pd.read_csv(path, dtype=str).fillna("")
    return {r["run_id"]: dict(r) for _, r in df.iterrows()}


def _save_index(path, index):
    import pandas as pd
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(list(index.values())).to_csv(path + ".part", index=False)
    os.replace(path + ".part", path)


def run_batch(runs, batch_name, pipeline_path=PIPELINE_PATH, batch_dir=BATCH_DIR, shared_cache_dir=SHARED_CACHE_DIR,
              extra_globals=None, quiet=True):
    """Run each configuration (dict of Cell 0 field -> value) through the whole
    pipeline. Returns the batch index (dict run_id -> record)."""
    import matplotlib.pyplot as plt
    src = open(pipeline_path).read()
    m = re.search(r'^__version__ = "([^"]+)"', src, flags=re.M)
    release = m.group(1) if m else "unknown"
    code = compile(src, pipeline_path, "exec")
    src_sha = hashlib.sha256(src.encode()).hexdigest()
    index_path = os.path.join(batch_dir, f"{batch_name}_index.csv")
    index = _load_index(index_path)
    os.makedirs(shared_cache_dir, exist_ok=True)
    print(f"ARF batch '{batch_name}': {len(runs)} runs, pipeline {os.path.basename(pipeline_path)} "
          f"(release {release}, sha256 {src_sha[:16]}), index {index_path}")
    installed = False
    _show = plt.show
    try:
        if quiet:
            plt.show = lambda *a, **k: None
        for i, cfg in enumerate(runs, 1):
            rid = run_id(cfg)
            rec = index.get(rid, {})
            if rec.get("status") == "done" and _manifest_matches(rec.get("out_dir"), release, cfg):
                print(f"[{i}/{len(runs)}] {rid}: done earlier, skipped")
                continue
            ns = {"__name__": "__main__", "ARF_CONFIG_OVERRIDES": dict(cfg), "ARF_SHARED_CACHE_DIR": shared_cache_dir,
                  "ARF_SKIP_INSTALL": installed}
            if quiet:
                ns["display"] = lambda *a, **k: None
            ns.update(extra_globals or {})
            buf = io.StringIO()
            t0 = time.time()
            status, err = "done", ""
            try:
                with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                    exec(code, ns)
            except KeyboardInterrupt:
                status, err = "interrupted", "KeyboardInterrupt"
                raise
            except BaseException:                       # record and continue with the next run
                status, err = "error", traceback.format_exc()[-1500:]
            finally:
                out_dir = ns.get("OUT_DIR", "")
                if "_ARF_CELLS_RUN" in ns and "1" in ns["_ARF_CELLS_RUN"]:
                    installed = True                    # Cell 1 ran once: skip pip from now on
                if out_dir and os.path.isdir(out_dir):
                    with open(os.path.join(out_dir, "arf_batch_console.txt"), "w") as f:
                        f.write(buf.getvalue())
                index[rid] = {"run_id": rid, "status": status, "release": release, "pipeline_sha256": src_sha,
                              "out_dir": out_dir, "minutes": f"{(time.time() - t0) / 60.0:.1f}",
                              "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S"),
                              "config": json.dumps(cfg, sort_keys=True, default=str), "error": err}
                _save_index(index_path, index)
                plt.close("all")
                del ns
                gc.collect()
            print(f"[{i}/{len(runs)}] {rid}: {status} in {index[rid]['minutes']} min -> {out_dir}"
                  + (f"\n    {err.strip().splitlines()[-1]}" if err else ""))
    finally:
        plt.show = _show
    n_done = sum(r.get("status") == "done" for r in index.values())
    print(f"Batch '{batch_name}': {n_done}/{len(index)} runs done. Index: {index_path}")
    return index


# =============================================================================
# THE DENSE-CONSTELLATION ARM (amendment A3, release 1.5.4)
# =============================================================================
PROTOCOL_WINDOWS = windows("2026-08-03T00:00:00", 8, 7)
DENSE_SETTINGS = {"active_roi": "ckgr_okavango", "fleet": "starlink", "proxy_min_sun_elev_deg": 5.0,
                  "cloud_sweep_lengths_km": [], "cloud_sweep_times_h": []}
DENSE_N_ALL = (100, 500, 1000, 2000, 4000, 8000)       # window 1
DENSE_N_SMALL = (100, 500, 1000)                        # windows 2-8
DENSE = ([dict(DENSE_SETTINGS, active_n=n, analysis_epoch=PROTOCOL_WINDOWS[0]) for n in DENSE_N_ALL if n != 8000]
         + [dict(DENSE_SETTINGS, active_n=n, analysis_epoch=w) for w in PROTOCOL_WINDOWS[1:] for n in DENSE_N_SMALL]
         + [dict(DENSE_SETTINGS, active_n=8000, analysis_epoch=PROTOCOL_WINDOWS[0])])

if __name__ == "__main__" and os.path.isdir("/content/drive"):
    run_batch(DENSE, "exploratory_v154_dense")
