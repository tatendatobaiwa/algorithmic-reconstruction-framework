# =============================================================================
# ARF 1.5.4 LAUNCH CELL (Colab): the dense-constellation arm (protocol amendment A3)
# =============================================================================
# Before running: finish the real-fleet batch first (release 1.5.3, launch cell
# arf_colab_launch_v153.py, "56/56 protocol runs PASS" in the audit), then upload
# ARF_v1.5.4_dense_release.zip to MyDrive (top level), or its unpacked folder as
# MyDrive/ARF_v1.5.4_dense_release/ARF_v1.5.4. Colab secrets SPACETRACK_USERNAME
# and SPACETRACK_PASSWORD as before. Fresh runtime, ONE notebook only, tab kept open.
#   STAGE = "check"  installs 1.5.4 beside 1.5.3 and checks it reproduces the
#                    real-fleet reference (parity, about 10 min)
#   STAGE = "batch"  the same, then the 27 runs of the arm (about 43 h in total,
#                    the N = 8000 run alone about 12 h; resumable between runs: if
#                    Colab disconnects, run the cell again; a run cut off part-way
#                    starts again from its beginning)
# 1.5.3 stays in ARF_code untouched (the real-fleet release of record); 1.5.4 is
# installed as algorithmic_reconstruction_framework_v154.py. The protocol file in
# ARF_code is replaced by the one with amendment A3 (the A2 file is kept as
# arf_confirmatory_protocol.a2.json). Secrets are tested for presence only.
# =============================================================================
STAGE = "check"                                   # "check" or "batch"
ROOT = "/content/drive/MyDrive"
ZIP = f"{ROOT}/ARF_v1.5.4_dense_release.zip"
FOLDER = f"{ROOT}/ARF_v1.5.4_dense_release/ARF_v1.5.4"
CODE = f"{ROOT}/ARF_code"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SHA_154 = "2b618ec442e71cf896e1104c6da89c94f2d6b6785a04c4a9236b92c21115a1cb"

import os, sys, json, shutil, zipfile, hashlib, subprocess, runpy

if ROOT.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
    from google.colab import drive
    drive.mount("/content/drive")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# 1. unpack (or take the unpacked folder) and verify every file against SHA256SUMS.txt
work = "/tmp/ARF_v1.5.4_unpacked"
shutil.rmtree(work, ignore_errors=True)
if os.path.exists(ZIP):
    zipfile.ZipFile(ZIP).extractall(work)
    rel = os.path.join(work, "ARF_v1.5.4")
elif os.path.isdir(FOLDER):
    rel = shutil.copytree(FOLDER, os.path.join(work, "ARF_v1.5.4"))
else:
    raise FileNotFoundError(f"neither {ZIP} nor {FOLDER} found: upload the 1.5.4 package to MyDrive first")
bad = []
for line in open(os.path.join(rel, "SHA256SUMS.txt")):
    h, name = line.split(None, 1)
    if sha(os.path.join(rel, name.strip())) != h:
        bad.append(name.strip())
assert not bad, f"checksum mismatch: {bad}"
assert sha(os.path.join(rel, "algorithmic_reconstruction_framework_v154.py")) == SHA_154, "not the 1.5.4 release file"
proto = json.load(open(os.path.join(rel, "arf_confirmatory_protocol.json")))
a3 = [a for a in proto["amendments"] if a["id"] == "A3"]
assert a3 and a3[0]["pipeline_sha256_of_arm"] == SHA_154 and proto["pipeline_sha256"] == SHA_153, \
    "protocol file without amendment A3 for this release"
print("1/5 package unpacked; every checksum OK")

# 2. the real-fleet release 1.5.3 must already be in ARF_code (never replaced here)
p153 = os.path.join(CODE, "algorithmic_reconstruction_framework.py")
assert os.path.exists(p153) and sha(p153) == SHA_153, \
    f"{p153} is not release 1.5.3: install the real-fleet release first (arf_colab_launch_v153.py)"
print("2/5 release 1.5.3 in ARF_code untouched (the real-fleet release of record)")

# 3. install 1.5.4 beside it, the arm's runner, the audit and the A3 protocol
pp = os.path.join(CODE, "arf_confirmatory_protocol.json")
if os.path.exists(pp) and sha(pp) != sha(os.path.join(rel, "arf_confirmatory_protocol.json")):
    keep = os.path.join(CODE, "arf_confirmatory_protocol.a2.json")
    if not os.path.exists(keep):
        shutil.copy2(pp, keep)
for f in ("algorithmic_reconstruction_framework_v154.py", "arf_batch_runner_v154.py", "arf_batch_audit.py",
          "arf_confirmatory_protocol.json", "arf_dense_summary.py"):
    dst = os.path.join(CODE, f)
    if f == "arf_batch_audit.py" and os.path.exists(dst) and sha(dst) != sha(os.path.join(rel, f)) \
            and not os.path.exists(dst + ".v153"):
        shutil.copy2(dst, dst + ".v153")
    shutil.copy2(os.path.join(rel, f), dst)
print(f"3/5 release 1.5.4 installed beside 1.5.3 in {CODE} (protocol with amendment A3)")

# 4. parity: 1.5.4 must reproduce the real-fleet reference of the five 15 October runs
try:
    import h3
    h3_ok = h3.__version__ == "4.5.0"
except ImportError:
    h3_ok = False
if not h3_ok:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "h3==4.5.0"])
exp = os.path.join(CODE, "arf_parity_expected_v153.json")
assert os.path.exists(exp), f"{exp} missing: it is installed by the 1.5.3 launch cell"
report_path = os.path.join(CODE, "arf_parity_report.json")
if os.path.exists(report_path):
    os.replace(report_path, report_path + ".previous")
r = subprocess.run([sys.executable, os.path.join(CODE, "arf_colab_parity.py"),
                    "--pipeline", os.path.join(CODE, "algorithmic_reconstruction_framework_v154.py"),
                    "--expected", exp, "--drive", ROOT], capture_output=True, text=True)
print("\n".join(l for l in r.stdout.splitlines() if not l.startswith("INFO")))
if r.returncode != 0 or not os.path.exists(report_path):
    print("parity script error output (last lines):\n" + "\n".join(r.stderr.splitlines()[-25:]))
    raise RuntimeError(f"parity script failed (exit code {r.returncode}); see the error output above")
rep = json.load(open(report_path))
st = {k: v["status"] for k, v in rep["runs"].items()}
assert st and all(s == "MATCH" for s in st.values()), f"parity failed: {st}"
print(f"4/5 parity of 1.5.4 on the real fleets: {len(st)}/{len(st)} MATCH")

# 5. the dense-constellation arm
if STAGE != "batch":
    print("5/5 skipped (STAGE = 'check'). Set STAGE = 'batch' and run again to start the 27 runs of the arm.")
else:
    try:
        from google.colab import userdata
        missing = [k for k in ("SPACETRACK_USERNAME", "SPACETRACK_PASSWORD") if not userdata.get(k)]
    except Exception as e:
        missing = [f"Colab secrets unavailable ({type(e).__name__})"]
    assert not missing, f"missing Space-Track secrets: {missing}"
    print("5/5 starting the dense-constellation arm: 27 runs (N = 8000 last)")
    runpy.run_path(os.path.join(CODE, "arf_batch_runner_v154.py"), run_name="__main__")
