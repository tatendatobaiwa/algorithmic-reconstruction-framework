# =============================================================================
# ARF PAPER 1 / AMENDMENT A4 INSTALL CELL (Colab)
# =============================================================================
# Upload ARF_paper1_A4.zip to MyDrive (do not unzip it there), then run this
# cell once. It:
#   1. extracts the zip to /content/arf_paper1_a4 (outside Drive) and checks
#      every file against SHA256SUMS.txt
#   2. checks ARF_code: release 1.5.3 (and 1.5.4 beside it) untouched, and the
#      protocol in place is the A3 protocol
#   3. checks that the A4 protocol in the zip equals the installed one with
#      amendment A4 appended and nothing else changed
#   4. keeps the A3 protocol as ARF_code/arf_confirmatory_protocol.a3.json,
#      installs the A4 protocol and the five A4 scripts
# Nothing is run; the cells to run next are printed. Secrets are not needed.
# Do NOT run the 1.5.3 or 1.5.4 launch cells again afterwards: each would put
# an older protocol back over A4.
# =============================================================================
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
ZIP = f"{ROOT}/ARF_paper1_A4.zip"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SHA_154 = "2b618ec442e71cf896e1104c6da89c94f2d6b6785a04c4a9236b92c21115a1cb"
SCRIPTS = ("arf_paper1_evidence.py", "arf_exact_economy_a4.py", "arf_sensitivity_runner_a4.py",
           "arf_sensitivity_summary_a4.py")

import os, json, shutil, hashlib, zipfile

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def install(root=ROOT, code=CODE, zip_path=ZIP, work="/content/arf_paper1_a4"):
    if root.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
        from google.colab import drive
        drive.mount("/content/drive")
    shutil.rmtree(work, ignore_errors=True)
    with zipfile.ZipFile(zip_path) as z:
        for n in z.namelist():
            assert not n.startswith("/") and ".." not in n.split("/"), f"unsafe path in zip: {n}"
        z.extractall(work)
    sums = {}
    for line in open(os.path.join(work, "SHA256SUMS.txt")):
        h, n = line.split(maxsplit=1)
        sums[n.strip()] = h
    for n, h in sums.items():
        assert sha(os.path.join(work, n)) == h, f"{n}: checksum mismatch (damaged upload?)"
    print(f"package: {len(sums)} files, checksums OK")

    assert sha(os.path.join(code, "algorithmic_reconstruction_framework.py")) == SHA_153, "ARF_code: release 1.5.3 changed"
    p154 = os.path.join(code, "algorithmic_reconstruction_framework_v154.py")
    assert not os.path.exists(p154) or sha(p154) == SHA_154, "ARF_code: release 1.5.4 changed"
    for need in ("arf_batch_runner.py", "arf_confirmatory_analysis.py", "arf_colab_parity.py"):
        assert os.path.exists(os.path.join(code, need)), f"ARF_code/{need} missing (1.5.3 release files)"
    cur_p = os.path.join(code, "arf_confirmatory_protocol.json")
    cur = json.load(open(cur_p))
    new = json.load(open(os.path.join(work, "arf_confirmatory_protocol.json")))
    ids = [a["id"] for a in cur["amendments"]]
    if ids == ["A1", "A2", "A3", "A4"] and cur == new:
        print("the A4 protocol is already installed")
    else:
        # from the A3 protocol, or from an earlier version of A4 (kept as .a4_prev.json)
        assert ids in (["A1", "A2", "A3"], ["A1", "A2", "A3", "A4"]), \
            f"ARF_code protocol has amendments {ids}, expected A1-A3 (or an earlier A4)"
        a = {k: v for k, v in new.items() if k not in ("amendments", "amended_utc")}
        b = {k: v for k, v in cur.items() if k not in ("amendments", "amended_utc")}
        assert a == b and new["amendments"][:3] == cur["amendments"][:3] and new["amendments"][3]["id"] == "A4" \
            and len(new["amendments"]) == 4, "the zip's protocol is not the installed A3 protocol with A4 appended"
        keep = os.path.join(code, "arf_confirmatory_protocol.a3.json" if len(ids) == 3 else
                            "arf_confirmatory_protocol.a4_prev.json")
        if len(ids) == 4 or not os.path.exists(keep):
            shutil.copy2(cur_p, keep)
        shutil.copy2(os.path.join(work, "arf_confirmatory_protocol.json"), cur_p)
        print(f"protocol: A4 installed (previous protocol kept as {keep})")
    for s in SCRIPTS:
        shutil.copy2(os.path.join(work, s), os.path.join(code, s))
        assert sha(os.path.join(code, s)) == sums[s]
    print("scripts installed:", ", ".join(SCRIPTS))
    print("\nNext, one per cell, in this order:\n"
          f"  !python {code}/arf_paper1_evidence.py          # minutes; read-only\n"
          f"  !python {code}/arf_exact_economy_a4.py          # about 1.5 h; read-only\n"
          f"  !python {code}/arf_sensitivity_runner_a4.py     # about 3.5-4 h, 48 runs; resumable\n"
          f"  !python {code}/arf_sensitivity_summary_a4.py    # minutes; read-only")


if __name__ == "__main__" or "get_ipython" in globals():
    install()
