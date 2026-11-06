# =============================================================================
# ARF AMENDMENT A8 INSTALL CELL (Colab): the real Sentinel-2 cloud check
# =============================================================================
# Upload ARF_realcloud_A8.zip to MyDrive (do not unzip it there) and run this
# cell once. It checks the zip against SHA256SUMS.txt; checks that ARF_code holds
# release 1.5.3, the A7 protocol; checks that the zip's protocol
# is the installed A7 protocol with amendment A8 appended and nothing else
# changed; keeps the A7 protocol as arf_confirmatory_protocol.a7.json; installs
# the A8 protocol and the A8 script. Nothing is run.
# Do NOT run an earlier install cell or any launch cell again afterwards.
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
ZIP = f"{ROOT}/ARF_realcloud_A8.zip"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SCRIPTS = ("arf_real_cloud_a8.py",)

import os, json, shutil, hashlib, zipfile


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def install(root=ROOT, code=CODE, zip_path=ZIP, work="/content/arf_realcloud_a8"):
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
    for need in ():
        assert os.path.exists(os.path.join(code, need)), f"ARF_code/{need} missing"
    cur_p = os.path.join(code, "arf_confirmatory_protocol.json")
    cur, new = json.load(open(cur_p)), json.load(open(os.path.join(work, "arf_confirmatory_protocol.json")))
    ids = [a["id"] for a in cur["amendments"]]
    if cur == new:
        print("the A8 protocol is already installed")
    else:
        assert ids in (["A1", "A2", "A3", "A4", "A5", "A6", "A7"], ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8"]), f"ARF_code protocol has {ids}"
        strip = lambda d: {k: v for k, v in d.items() if k not in ("amendments", "amended_utc")}
        assert strip(cur) == strip(new) and new["amendments"][:7] == cur["amendments"][:7] and \
            [a["id"] for a in new["amendments"]] == ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8"], \
            "the zip's protocol is not the installed A7 protocol with A8 appended"
        keep = os.path.join(code, "arf_confirmatory_protocol.a7.json" if len(ids) == 7 else
                            "arf_confirmatory_protocol.a8_prev.json")
        shutil.copy2(cur_p, keep)
        shutil.copy2(os.path.join(work, "arf_confirmatory_protocol.json"), cur_p)
        print(f"protocol: A8 installed (previous kept as {os.path.basename(keep)})")
    for s in SCRIPTS:
        shutil.copy2(os.path.join(work, s), os.path.join(code, s))
        assert sha(os.path.join(code, s)) == sums[s]
    print("scripts installed:", ", ".join(SCRIPTS))
    print(f"\nNext:\n  !pip install -q rasterio h3==4.5.0\n  !python {code}/arf_real_cloud_a8.py   (reads Sentinel-2 scenes over the internet; resumable)")


if __name__ == "__main__" or "get_ipython" in globals():
    install()
