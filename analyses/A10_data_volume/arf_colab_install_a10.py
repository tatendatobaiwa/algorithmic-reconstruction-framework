# =============================================================================
# ARF AMENDMENT A10 INSTALL CELL (Colab): the corrected data-volume analysis
# =============================================================================
# Upload ARF_datavolume_A10.zip to MyDrive (do not unzip it there) and run this
# cell once. It checks the zip against SHA256SUMS.txt; checks that ARF_code holds
# release 1.5.3, the A9 protocol, arf_colab_parity.py and arf_exact_economy_a4.py; checks that the zip's protocol
# is the installed A9 protocol with amendment A10 appended and nothing else
# changed; keeps the A9 protocol as arf_confirmatory_protocol.a9.json; installs
# the A10 protocol and the A10 script. Nothing is run.
# Do NOT run an earlier install cell or any launch cell again afterwards.
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
ZIP = f"{ROOT}/ARF_datavolume_A10.zip"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SCRIPTS = ("arf_data_volume_a10.py",)

import os, json, shutil, hashlib, zipfile


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def install(root=ROOT, code=CODE, zip_path=ZIP, work="/content/arf_datavolume_a10"):
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
    for need in ("arf_colab_parity.py", "arf_exact_economy_a4.py"):
        assert os.path.exists(os.path.join(code, need)), f"ARF_code/{need} missing"
    cur_p = os.path.join(code, "arf_confirmatory_protocol.json")
    cur, new = json.load(open(cur_p)), json.load(open(os.path.join(work, "arf_confirmatory_protocol.json")))
    ids = [a["id"] for a in cur["amendments"]]
    if cur == new:
        print("the A10 protocol is already installed")
    else:
        assert ids in (["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"], ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"]), f"ARF_code protocol has {ids}"
        strip = lambda d: {k: v for k, v in d.items() if k not in ("amendments", "amended_utc")}
        assert strip(cur) == strip(new) and new["amendments"][:9] == cur["amendments"][:9] and \
            [a["id"] for a in new["amendments"]] == ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"], \
            "the zip's protocol is not the installed A9 protocol with A10 appended"
        keep = os.path.join(code, "arf_confirmatory_protocol.a9.json" if len(ids) == 9 else
                            "arf_confirmatory_protocol.a10_prev.json")
        shutil.copy2(cur_p, keep)
        shutil.copy2(os.path.join(work, "arf_confirmatory_protocol.json"), cur_p)
        print(f"protocol: A10 installed (previous kept as {os.path.basename(keep)})")
    for s in SCRIPTS:
        shutil.copy2(os.path.join(work, s), os.path.join(code, s))
        assert sha(os.path.join(code, s)) == sums[s]
    print("scripts installed:", ", ".join(SCRIPTS))
    print(f"\nNext:\n  !pip install -q h3==4.5.0\n  !python {code}/arf_data_volume_a10.py   (about 30 min)")


if __name__ == "__main__" or "get_ipython" in globals():
    install()
