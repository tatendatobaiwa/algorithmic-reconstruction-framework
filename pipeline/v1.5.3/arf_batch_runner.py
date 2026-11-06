# =============================================================================
# ARF BATCH RUNNER (v1.5.3): many configurations in one Colab session
# =============================================================================
# Runs the whole pipeline (Cells 0-9 of algorithmic_reconstruction_framework.py)
# once per configuration, each run in a FRESH namespace, so nothing carries over
# from one run to the next. Each run writes to its own Drive folder, which Cell 0
# names from the ROI, the fleet and, for a past window, its start date
# (ARF_Outputs_<roi>_<fleet tag>[_from<YYYYMMDD>][_without-...]); from v1.5.2 a
# folder holding a complete run of another release is kept and the run writes
# to <folder>_v<release digits> (the index records the folder actually used).
#
# Resumable: the batch index (one row per run) is saved after every run; a run
# already marked "done" whose folder holds a manifest from the same release with
# the same settings is skipped. A run that fails is recorded with the end of its
# traceback, and the batch moves on.
#
# Use (Colab):
#   1. Upload the release file to MyDrive/ARF_code/ (PIPELINE_PATH).
#   2. Past windows read orbital elements from Space-Track: add the Colab secrets
#      SPACETRACK_USERNAME and SPACETRACK_PASSWORD (as for a single past run).
#   3. Choose the runs (below: the confirmatory protocol's list) and run this cell.
# Console output of each run: <run folder>/arf_batch_console.txt. Figures are not
# shown on screen during a batch (they are saved in each run folder as usual).
# Etiquette: CelesTrak lists and Space-Track element sets are cached in
# SHARED_CACHE_DIR and reused by every run (each cache file is keyed by its
# group, fleet or epoch); runs are sequential.
# v1.5.2: batch names carry the release; protocol amendment A1 (14 October 2026,
# before any window was run) adds the exploratory runs without Sentinel-2A.
# v1.5.3: release 1.5.3 (protocol amendment A2, 18 October 2026, before any
# window was run: figures and report labels only, every computed value as in
# 1.5.2); batch names confirmatory_v153_*. A protocol window has no folder yet,
# so each run writes to its plain folder; a folder holding a complete run of
# another release would be kept and the run would write to <folder>_v153.
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

PIPELINE_PATH = "/content/drive/MyDrive/ARF_code/algorithmic_reconstruction_framework.py"
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
# THE CONFIRMATORY PROTOCOL'S RUNS (frozen 2026-10-12, amended 2026-10-14 before
# any window was run; see arf_confirmatory_protocol.json): eight 48-hour windows
# starting Mondays 00:00 UTC, 3 August - 21 September 2026, none analysed before.
# =============================================================================
PROTOCOL_WINDOWS = windows("2026-08-03T00:00:00", 8, 7)
PRIMARY = [{"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative",
            "analysis_epoch": w} for w in PROTOCOL_WINDOWS]
SECONDARY = [{"active_roi": roi, "fleet": fl, "planet_swath_case": case, "analysis_epoch": w}
             for w in PROTOCOL_WINDOWS
             for roi, fl, case in (("ckgr_okavango", "earth-resources+planet", "optimistic"),
                                   ("ckgr_okavango", "earth-resources", "conservative"),
                                   ("lower_zambezi", "earth-resources+planet", "conservative"),
                                   ("karoo_nature_reserve", "earth-resources+planet", "conservative"))]
# Amendment A1 (exploratory): the primary design and the open-data fleet without
# Sentinel-2A, whose extension-campaign plan images a region on a reduced schedule.
EXPLORATORY = [{"active_roi": "ckgr_okavango", "fleet": fl, "planet_swath_case": "conservative",
                "analysis_epoch": w, "exclude_satellites": ["SENTINEL-2A"]}
               for w in PROTOCOL_WINDOWS for fl in ("earth-resources+planet", "earth-resources")]

if __name__ == "__main__" and os.path.isdir("/content/drive"):
    # The 8 primary runs first; then the 32 secondary and 16 exploratory runs (resumable).
    run_batch(PRIMARY, "confirmatory_v153_primary")
    run_batch(SECONDARY, "confirmatory_v153_secondary")
    run_batch(EXPLORATORY, "confirmatory_v153_exploratory")
