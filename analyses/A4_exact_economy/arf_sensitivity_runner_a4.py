# =============================================================================
# ARF SENSITIVITY BATCH (protocol amendment A4, POST HOC): release 1.5.3
# unchanged, primary design, the eight protocol windows, one setting changed
# per run. Labelled post hoc: decided after the confirmatory results were seen,
# so none of it can change a confirmatory decision.
# =============================================================================
# Variants (each against the confirmatory run of the same window):
#   k05, k10, k15, k20   max_selection_steps = 5, 10, 15, 20 (protocol: 30):
#                        where the image budget starts to bind, i.e. where ECC
#                        leaves the physical floor (H1a is met trivially at 30)
#   earthenv             cloud_probability_source = "earthenv": per-tile
#                        monthly cloud frequency (Wilson & Jetz 2016) instead
#                        of the constant p = 0.35 (H2 under a realistic,
#                        dry-season cloud probability)
#   ilp48                the exact integer optimum on ONE full 48-h instance
#                        (ilp_horizon_hours = 48, ilp_instances = 1,
#                        ilp_max_pairs = 3,000,000, ilp_time_limit_s = 900):
#                        ECC's greedy cover against the exact minimum images at
#                        full scale. Every other output must equal the
#                        confirmatory run's (a determinism check).
# 48 runs, about 3.5 to 4 min each in Colab (the confirmatory runs took 2.9 to
# 4.0 min); ilp48 runs take longer by the solve time.
#
# SAFETY: Cell 0 of release 1.5.3 names its Drive folder from the ROI, fleet
# and window only, so a changed budget would write into the CONFIRMATORY
# folder of the same window. This runner therefore runs each configuration
# with google.colab hidden from the pipeline (Cell 0 then writes to the
# working directory) and the working directory set to its own folder under
# MyDrive/ARF_sensitivity_A4/<variant>/<window>. The pipeline file is not
# modified. Before the batch, every file in every confirmatory run folder is
# recorded (size, modification time, and the SHA-256 of each JSON and CSV);
# after each run all are checked, and the batch stops at the first change. Boundaries (ARF_boundaries), the climatology
# (ARF_climatology) and the shared caches are found by their absolute Drive
# paths, as in a normal run.
# The past windows' element sets come from ARF_batch/shared_cache (written by
# the confirmatory runs), so no Space-Track query is made; the runner checks
# that the cache is present before it starts.
#
# Use (Colab): mount Drive, then, as its own process (not runpy in the
# notebook kernel, whose own Colab modules must stay importable):
#   !python /content/drive/MyDrive/ARF_code/arf_sensitivity_runner_a4.py
# Resumable: completed runs are skipped.
# =============================================================================
import os
import sys
import glob
import json
import runpy
import hashlib
import contextlib

# Run as its own process (!python ...), so hiding google.colab affects the pipeline runs only, never the
# notebook kernel; figures are written to files with the non-interactive backend.
os.environ["MPLBACKEND"] = "Agg"

ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
BATCH = f"{ROOT}/ARF_batch"
SENS_ROOT = f"{ROOT}/ARF_sensitivity_A4"
BATCH_NAME = "sensitivity_v153_a4"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
CONFIRMATORY_INDEXES = ("confirmatory_v153_primary_index.csv", "confirmatory_v153_secondary_index.csv",
                        "confirmatory_v153_exploratory_index.csv")
VARIANTS = {
    "k05": {"max_selection_steps": 5},
    "k10": {"max_selection_steps": 10},
    "k15": {"max_selection_steps": 15},
    "k20": {"max_selection_steps": 20},
    "earthenv": {"cloud_probability_source": "earthenv"},
    "ilp48": {"ilp_horizon_hours": 48, "ilp_instances": 1, "ilp_max_pairs": 3_000_000, "ilp_time_limit_s": 900.0},
}
PRIMARY_DESIGN = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def confirmatory_folders(batch_dir):
    import pandas as pd
    out = []
    for name in CONFIRMATORY_INDEXES:
        p = os.path.join(batch_dir, name)
        if os.path.exists(p):
            idx = pd.read_csv(p, dtype=str).fillna("")
            out += [d for d in idx.loc[idx.status == "done", "out_dir"] if d and os.path.isdir(d)]
    return sorted(set(out))


def snapshot(folders):
    """Every file (recursively) in the given folders: size and modification time, plus the SHA-256 of
    the JSON and CSV files (the results); cheap enough to repeat after every run on Drive."""
    snap = {}
    for d in folders:
        for p in glob.glob(os.path.join(d, "**", "*"), recursive=True):
            if os.path.isfile(p):
                st = os.stat(p)
                snap[p] = (st.st_size, st.st_mtime_ns, sha(p) if p.endswith((".json", ".csv")) else None)
    return snap


def changed(before, folders):
    now = snapshot(folders)
    diff = sorted(set(before) ^ set(now)) + sorted(p for p in before if p in now and before[p] != now[p])
    return diff


@contextlib.contextmanager
def outside_colab_folder(work_dir):
    """Hide google.colab from the pipeline (Cell 0 then writes to the working directory) and work in
    work_dir. Restores sys.modules and the working directory afterwards."""
    saved = {k: sys.modules[k] for k in list(sys.modules) if k == "google.colab" or k.startswith("google.colab.")}
    cwd = os.getcwd()
    for k in saved:
        del sys.modules[k]
    sys.modules["google.colab"] = None            # `import google.colab` raises ImportError
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)
    try:
        yield
    finally:
        os.chdir(cwd)
        sys.modules.pop("google.colab", None)
        sys.modules.update(saved)


def build_runs(windows, variants=VARIANTS, order=("ilp48", "k10", "k20", "k05", "k15", "earthenv")):
    runs = []
    for v in order:
        for w in windows:
            runs.append((v, w, {**PRIMARY_DESIGN, "analysis_epoch": w, **variants[v]}))
    return runs


def main(argv=None, root=ROOT, code=CODE, batch=BATCH, sens_root=SENS_ROOT, only=None, pipeline=None):
    pipeline = pipeline or os.path.join(code, "algorithmic_reconstruction_framework.py")
    if root.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    assert sha(pipeline) == SHA_153, f"{pipeline} is not release 1.5.3"
    R = runpy.run_path(os.path.join(code, "arf_batch_runner.py"), run_name="arf_batch_runner")
    W = R["PROTOCOL_WINDOWS"]
    cache = os.path.join(batch, "shared_cache")
    hist = glob.glob(os.path.join(cache, "arf_gp_history_*"))
    assert len(hist) >= len(W), (f"{cache} holds {len(hist)} element-set histories; the confirmatory runs should "
                                 f"have left at least {len(W)} (one per window)")
    runs = build_runs(W)
    if only:
        runs = [r for r in runs if r[0] in only]
    folders = confirmatory_folders(batch)
    assert folders, f"no completed confirmatory run folder found through {CONFIRMATORY_INDEXES} in {batch}"
    print(f"A4 sensitivity batch: {len(runs)} runs of release 1.5.3; {len(folders)} confirmatory folders protected")
    before = snapshot(folders)
    for i, (v, w, cfg) in enumerate(runs, 1):
        work = os.path.join(sens_root, v, w[:10])
        print(f"--- [{i}/{len(runs)}] {v} {w[:10]} -> {work}")
        with outside_colab_folder(work):
            R["run_batch"]([cfg], BATCH_NAME, pipeline_path=pipeline, batch_dir=batch,
                           shared_cache_dir=cache)
        bad = changed(before, folders)
        assert not bad, f"STOPPED: a confirmatory folder changed during run {v} {w}: {bad[:5]}"
    print(f"Done. Index: {os.path.join(batch, BATCH_NAME + '_index.csv')}; confirmatory folders unchanged.")


if __name__ == "__main__":
    main()
