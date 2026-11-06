# =============================================================================
# ARF AMENDMENTS A12 + A13 (Colab): install, then run on the 16 primary windows
# =============================================================================
# Upload ARF_posthoc_A12A13.zip to MyDrive (do NOT unzip it there) and run this
# cell once. It
#   1. checks the zip against SHA256SUMS.txt;
#   2. checks that ARF_code holds release 1.5.3 and the A11 protocol, and that the
#      zip's protocol is that protocol with A12 and A13 appended and nothing else
#      changed; keeps the A11 protocol as arf_confirmatory_protocol.a11.json and
#      installs the new one with the two scripts;
#   3. runs A12 (floor decomposition) and A13 (PlanetScope acquisition fraction,
#      20 seeds) READ-ONLY on the 8 dry-season primary runs and the 8 wet-season
#      A5 runs (p = 0.35), each script refusing a run whose floors it cannot
#      reproduce exactly;
#   4. writes everything to a NEW folder MyDrive/ARF_posthoc_A12A13_results_<time>
#      and zips it to MyDrive/ARF_A12A13_results_<time>.zip for upload.
# Nothing is deleted or overwritten. Expected run time: a few minutes per window.
# Do NOT run an earlier install cell or any launch cell again afterwards.
# =============================================================================
ROOT = "/content/drive/MyDrive"
ZIP_NAME = "ARF_posthoc_A12A13.zip"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
DRY = ["2026-08-03", "2026-08-10", "2026-08-17", "2026-08-24", "2026-08-31", "2026-09-07", "2026-09-14", "2026-09-21"]
WET = ["2026-01-05", "2026-01-12", "2026-01-19", "2026-01-26", "2026-02-02", "2026-02-09", "2026-02-16", "2026-02-23"]
NEED = ("arf_refined_events.csv", "arf_onset_robustness.csv", "arf_mosaic_tiles.csv", "arf_cell2_summary.json")

import os, sys, json, shutil, hashlib, zipfile, subprocess, datetime


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run_folder(path):
    """The one folder at or under `path` holding a complete run's outputs."""
    hits = [r for r, _, fs in os.walk(path) if all(n in fs for n in NEED)]
    assert len(hits) == 1, f"{path}: expected one run folder with {NEED}, found {len(hits)}: {hits[:3]}"
    return hits[0]


def main(root=ROOT, work="/content/arf_a12a13"):
    if root.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    code = os.path.join(root, "ARF_code")
    # 1. package
    shutil.rmtree(work, ignore_errors=True)
    with zipfile.ZipFile(os.path.join(root, ZIP_NAME)) as z:
        for n in z.namelist():
            assert not n.startswith("/") and ".." not in n.split("/"), f"unsafe path in zip: {n}"
        z.extractall(work)
    sums = {}
    for line in open(os.path.join(work, "SHA256SUMS.txt")).read().splitlines():
        h, n = line.split(maxsplit=1)
        sums[n.strip()] = h
    for n, h in sums.items():
        assert sha(os.path.join(work, n)) == h, f"{n}: checksum mismatch (damaged upload?)"
    print(f"package: {len(sums)} files, checksums OK")
    # 2. install
    assert sha(os.path.join(code, "algorithmic_reconstruction_framework.py")) == SHA_153, "ARF_code: release 1.5.3 changed"
    cur_p = os.path.join(code, "arf_confirmatory_protocol.json")
    cur, new = json.load(open(cur_p)), json.load(open(os.path.join(work, "arf_confirmatory_protocol.json")))
    ids = [a["id"] for a in cur["amendments"]]
    if cur == new:
        print("the A12/A13 protocol is already installed")
    else:
        assert ids == [f"A{i}" for i in range(1, 12)], f"ARF_code protocol has {ids}, expected A1-A11"
        strip = lambda d: {k: v for k, v in d.items() if k not in ("amendments", "amended_utc")}
        assert strip(cur) == strip(new) and new["amendments"][:11] == cur["amendments"] and \
            [a["id"] for a in new["amendments"]] == [f"A{i}" for i in range(1, 14)], \
            "the zip's protocol is not the installed A11 protocol with A12 and A13 appended"
        keep = os.path.join(code, "arf_confirmatory_protocol.a11.json")
        assert not os.path.exists(keep), f"{keep} exists already; not overwritten"
        shutil.copy2(cur_p, keep)
        shutil.copy2(os.path.join(work, "arf_confirmatory_protocol.json"), cur_p)
        print("installed: protocol with A12 and A13 (A11 protocol kept as arf_confirmatory_protocol.a11.json)")
    for s in ("arf_floor_decomposition.py", "arf_planet_thinning.py"):
        dst = os.path.join(code, s)
        if os.path.exists(dst):
            assert sha(dst) == sha(os.path.join(work, s)), f"ARF_code/{s} differs from the package; not overwritten"
        else:
            shutil.copy2(os.path.join(work, s), dst)
    # 3. run
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = os.path.join(root, f"ARF_posthoc_A12A13_results_{stamp}")
    os.makedirs(out)
    runs = [("dry", d, os.path.join(root, "ARF_Outputs_ckgr_okavango_earth-resources_planet-conservative_from"
                                     + d.replace("-", ""))) for d in DRY]
    runs += [("wet", d, os.path.join(root, "ARF_wetseason_A5", "p035", d)) for d in WET]
    log = {"stamp": stamp, "protocol_sha256": sha(cur_p), "windows": []}
    for season, d, path in runs:
        rd = run_folder(path)
        tag = f"{season}_{d}"
        rec = {"season": season, "window": d, "run_dir": rd}
        for name, args in (("A12", [os.path.join(work, "arf_floor_decomposition.py"), rd, d + "T00:00:00",
                                    os.path.join(out, tag + "_decomp.csv")]),
                           ("A13", [os.path.join(work, "arf_planet_thinning.py"), rd, os.path.join(out, tag + "_thin"),
                                    "--seeds", "20"])):
            r = subprocess.run([sys.executable, "-I"] + args, capture_output=True, text=True)
            open(os.path.join(out, f"{tag}_{name}.log"), "w").write(r.stdout + r.stderr)
            rec[name] = "ok" if r.returncode == 0 else f"FAILED (exit {r.returncode}; see {tag}_{name}.log)"
            print(f"{tag} {name}: {rec[name]}")
        log["windows"].append(rec)
    json.dump(log, open(os.path.join(out, "run_log.json"), "w"), indent=1)
    zpath = os.path.join(root, f"ARF_A12A13_results_{stamp}.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(out)):
            z.write(os.path.join(out, f), f)
    failed = [w for w in log["windows"] if w["A12"] != "ok" or w["A13"] != "ok"]
    print(f"\nDone: {len(runs) - len(failed)} of {len(runs)} windows complete. Upload {zpath} to the chat."
          + (f" FAILED: {[(w['season'], w['window']) for w in failed]}" if failed else ""))
    return zpath


if __name__ == "__main__" or "google.colab" in sys.modules:
    main()
