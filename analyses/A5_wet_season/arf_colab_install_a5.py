# =============================================================================
# ARF AMENDMENT A5 INSTALL CELL (Colab): the wet-season arm
# =============================================================================
# Upload ARF_wetseason_A5.zip to MyDrive (do not unzip it there) and run this
# cell once. It checks the zip against SHA256SUMS.txt; checks that ARF_code holds
# release 1.5.3 and the A4 protocol and A4 runner; checks that the zip's protocol
# is the installed A4 protocol with amendment A5 appended and nothing else
# changed; keeps the A4 protocol as arf_confirmatory_protocol.a4.json; installs
# the A5 protocol and the three A5 scripts. Nothing is run.
# Do NOT run the A4 install cell or any launch cell again afterwards.
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
ZIP = f"{ROOT}/ARF_wetseason_A5.zip"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SCRIPTS = ("arf_colab_prefetch_wet_a5.py", "arf_wetseason_runner_a5.py", "arf_wetseason_summary_a5.py")

import os, json, shutil, hashlib, zipfile


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def install(root=ROOT, code=CODE, zip_path=ZIP, work="/content/arf_wetseason_a5"):
    if root.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    shutil.rmtree(work, ignore_errors=True)
    with zipfile.ZipFile(zip_path) as z:
        for n in z.namelist():
            assert not n.startswith("/") and ".." not in n.split("/"), f"unsafe path in zip: {n}"
        z.extractall(work)
    sums = dict(reversed(l.split(maxsplit=1)) for l in open(os.path.join(work, "SHA256SUMS.txt")).read().splitlines())
    sums = {k.strip(): v for k, v in sums.items()}
    for n, h in sums.items():
        assert sha(os.path.join(work, n)) == h, f"{n}: checksum mismatch (damaged upload?)"
    print(f"package: {len(sums)} files, checksums OK")
    assert sha(os.path.join(code, "algorithmic_reconstruction_framework.py")) == SHA_153, "ARF_code: release 1.5.3 changed"
    for need in ("arf_batch_runner.py", "arf_confirmatory_analysis.py", "arf_sensitivity_runner_a4.py"):
        assert os.path.exists(os.path.join(code, need)), f"ARF_code/{need} missing"
    cur_p = os.path.join(code, "arf_confirmatory_protocol.json")
    cur, new = json.load(open(cur_p)), json.load(open(os.path.join(work, "arf_confirmatory_protocol.json")))
    ids = [a["id"] for a in cur["amendments"]]
    if cur == new:
        print("the A5 protocol is already installed")
    else:
        assert ids in (["A1", "A2", "A3", "A4"], ["A1", "A2", "A3", "A4", "A5"]), f"ARF_code protocol has {ids}"
        strip = lambda d: {k: v for k, v in d.items() if k not in ("amendments", "amended_utc")}
        assert strip(cur) == strip(new) and new["amendments"][:4] == cur["amendments"][:4] and \
            [a["id"] for a in new["amendments"]] == ["A1", "A2", "A3", "A4", "A5"], \
            "the zip's protocol is not the installed A4 protocol with A5 appended"
        keep = os.path.join(code, "arf_confirmatory_protocol.a4.json" if len(ids) == 4 else
                            "arf_confirmatory_protocol.a5_prev.json")
        shutil.copy2(cur_p, keep)
        shutil.copy2(os.path.join(work, "arf_confirmatory_protocol.json"), cur_p)
        print(f"protocol: A5 installed (previous kept as {os.path.basename(keep)})")
    for s in SCRIPTS:
        shutil.copy2(os.path.join(work, s), os.path.join(code, s))
        assert sha(os.path.join(code, s)) == sums[s]
    print("scripts installed:", ", ".join(SCRIPTS))
    print("\nNext:\n  1. download MODCF_monthlymean_01.tif and _02.tif from https://www.earthenv.org/cloud into "
          "MyDrive/ARF_climatology\n  2. run arf_colab_prefetch_wet_a5.py in a cell (Space-Track secrets set, tab open)\n"
          f"  3. !pip install -q h3==4.5.0\n     !python {code}/arf_wetseason_runner_a5.py\n"
          f"  4. !python {code}/arf_wetseason_summary_a5.py")


if __name__ == "__main__" or "get_ipython" in globals():
    install()
