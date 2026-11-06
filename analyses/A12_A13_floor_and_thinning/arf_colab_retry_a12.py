# =============================================================================
# ARF AMENDMENT A12 RETRY (Colab): the floor decomposition only
# =============================================================================
# Run after arf_colab_run_a12_a13.py, in the same or a new runtime. It
#   1. prints the end of the earlier A12 log of 3 August, so the cause is on record;
#   2. checks the installed protocol and that ARF_code/arf_floor_decomposition.py
#      is the registered A12 script (SHA-256 in the protocol);
#   3. installs h3 4.5.0 (the version the pipeline pins) if this runtime lacks it;
#   4. runs A12 READ-ONLY on the 16 windows into a NEW folder
#      MyDrive/ARF_posthoc_A12_results_<time>, zipped to MyDrive/ARF_A12_results_<time>.zip.
# A13's results from the first cell are complete and are not rerun.
# Nothing is deleted or overwritten.
# =============================================================================
ROOT = "/content/drive/MyDrive"
H3_VERSION = "4.5.0"
DRY = ["2026-08-03", "2026-08-10", "2026-08-17", "2026-08-24", "2026-08-31", "2026-09-07", "2026-09-14", "2026-09-21"]
WET = ["2026-01-05", "2026-01-12", "2026-01-19", "2026-01-26", "2026-02-02", "2026-02-09", "2026-02-16", "2026-02-23"]
NEED = ("arf_refined_events.csv", "arf_onset_robustness.csv", "arf_mosaic_tiles.csv", "arf_cell2_summary.json")

import os, sys, glob, json, hashlib, zipfile, subprocess, datetime


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run_folder(path):
    hits = [r for r, _, fs in os.walk(path) if all(n in fs for n in NEED)]
    assert len(hits) == 1, f"{path}: expected one run folder with {NEED}, found {len(hits)}: {hits[:3]}"
    return hits[0]


def h3_version(py):
    r = subprocess.run([py, "-I", "-c", "import h3; print(h3.__version__)"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def main(root=ROOT, py=sys.executable):
    if root.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    code = os.path.join(root, "ARF_code")
    # 1. the earlier failure, on record
    prev = sorted(glob.glob(os.path.join(root, "ARF_posthoc_A12A13_results_*", "dry_2026-08-03_A12.log")))
    if prev:
        print(f"--- end of the earlier A12 log ({prev[-1]}):")
        print("\n".join(open(prev[-1]).read().strip().splitlines()[-4:]))
        print("---")
    # 2. protocol and script
    prot = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
    a12 = [a for a in prot["amendments"] if a["id"] == "A12"]
    assert a12, "the installed protocol has no A12: run arf_colab_run_a12_a13.py first"
    script = os.path.join(code, "arf_floor_decomposition.py")
    assert sha(script) == a12[0]["script"]["arf_floor_decomposition.py"], "ARF_code/arf_floor_decomposition.py is not the registered A12 script"
    # 3. h3
    v = h3_version(py)
    if v != H3_VERSION:
        print(f"h3 in this runtime: {v or 'not installed'}; installing h3=={H3_VERSION} (the pipeline's pinned version)")
        r = subprocess.run([py, "-m", "pip", "install", "-q", f"h3=={H3_VERSION}"], capture_output=True, text=True)
        assert r.returncode == 0, "pip install h3 failed:\n" + r.stderr[-1500:]
        v = h3_version(py)
        assert v == H3_VERSION, f"h3 still {v} after install"
    print(f"h3 {v}")
    # 4. run
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = os.path.join(root, f"ARF_posthoc_A12_results_{stamp}")
    os.makedirs(out)
    runs = [("dry", d, os.path.join(root, "ARF_Outputs_ckgr_okavango_earth-resources_planet-conservative_from"
                                     + d.replace("-", ""))) for d in DRY]
    runs += [("wet", d, os.path.join(root, "ARF_wetseason_A5", "p035", d)) for d in WET]
    log = {"stamp": stamp, "h3": v, "script_sha256": sha(script), "windows": []}
    for season, d, path in runs:
        rd = run_folder(path)
        tag = f"{season}_{d}"
        r = subprocess.run([py, "-I", script, rd, d + "T00:00:00", os.path.join(out, tag + "_decomp.csv")],
                           capture_output=True, text=True)
        open(os.path.join(out, f"{tag}_A12.log"), "w").write(r.stdout + r.stderr)
        status = "ok" if r.returncode == 0 else f"FAILED (exit {r.returncode})"
        print(f"{tag} A12: {status}")
        if r.returncode != 0:
            print("   " + "\n   ".join((r.stdout + r.stderr).strip().splitlines()[-3:]))
        log["windows"].append({"season": season, "window": d, "run_dir": rd, "A12": status})
    json.dump(log, open(os.path.join(out, "run_log.json"), "w"), indent=1)
    zpath = os.path.join(root, f"ARF_A12_results_{stamp}.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(out)):
            z.write(os.path.join(out, f), f)
    n_ok = sum(w["A12"] == "ok" for w in log["windows"])
    print(f"\nDone: A12 complete for {n_ok} of {len(runs)} windows. Upload {zpath} and the earlier "
          f"ARF_A12A13_results_<time>.zip (it holds A13) to the chat.")
    return zpath


if __name__ == "__main__" or "google.colab" in sys.modules:
    main()
