# =============================================================================
# ARF 1.5.3 LAUNCH CELL (Colab): install the release, check parity, run the batch
# =============================================================================
# Before running: upload ARF_v1.5.3_release.zip to MyDrive (top level), or the unpacked
# folder as MyDrive/ARF_v1.5.3_release/ARF_v1.5.3 (both are checked the same way), and add the
# Colab secrets SPACETRACK_USERNAME and SPACETRACK_PASSWORD (key icon, left bar;
# "Notebook access" on). Then run this cell in a fresh runtime.
#   STAGE = "check"  installs and runs the parity check only (about 10 min)
#   STAGE = "batch"  the same, then the 56 protocol runs (about 4-5 h); the batch
#                    is resumable: if Colab disconnects, run the cell again
# The cell stops if a checksum fails, a parity run does not MATCH, or a secret
# is missing. Secrets are only tested for presence, never printed.
# =============================================================================
STAGE = "check"                                   # "check" or "batch"
ROOT = "/content/drive/MyDrive"
ZIP = f"{ROOT}/ARF_v1.5.3_release.zip"
CODE = f"{ROOT}/ARF_code"
EXPECTED_SHA = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"

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


# 1. unpack (or take the folder Drive unpacked on upload) and verify every file
#    against SHA256SUMS.txt
FOLDER = f"{ROOT}/ARF_v1.5.3_release/ARF_v1.5.3"
work = "/tmp/ARF_v1.5.3_unpacked"
shutil.rmtree(work, ignore_errors=True)
if os.path.exists(ZIP):
    zipfile.ZipFile(ZIP).extractall(work)
    rel = os.path.join(work, "ARF_v1.5.3")
elif os.path.isdir(FOLDER):
    rel = shutil.copytree(FOLDER, os.path.join(work, "ARF_v1.5.3"))
else:
    raise FileNotFoundError(f"neither {ZIP} nor {FOLDER} found: upload the release to MyDrive first")
bad = []
for line in open(os.path.join(rel, "SHA256SUMS.txt")):
    h, name = line.split(None, 1)
    if sha(os.path.join(rel, name.strip())) != h:
        bad.append(name.strip())
assert not bad, f"checksum mismatch: {bad}"
assert sha(os.path.join(rel, "algorithmic_reconstruction_framework.py")) == EXPECTED_SHA, "not the 1.5.3 release file"
print("1/4 release unpacked; every checksum OK")

# 2. install into ARF_code (earlier files kept with a .v152 suffix, never deleted)
os.makedirs(CODE, exist_ok=True)
for f in ("algorithmic_reconstruction_framework.py", "arf_batch_runner.py", "arf_colab_parity.py",
          "arf_parity_expected_v153.json", "arf_parity_expected_v153_12oct.json",
          "arf_confirmatory_analysis.py", "arf_confirmatory_protocol.json"):
    dst = os.path.join(CODE, f)
    if os.path.exists(dst) and sha(dst) != sha(os.path.join(rel, f)) and not os.path.exists(dst + ".v152"):
        shutil.copy2(dst, dst + ".v152")
    shutil.copy2(os.path.join(rel, f), dst)
print(f"2/4 release 1.5.3 installed in {CODE}")

# 3. parity: Cells 6-8 on the five 15 October runs must MATCH the reference.
#    A fresh runtime has no h3; install the version Cell 1 validates (4.5.0).
try:
    import h3
    h3_ok = h3.__version__ == "4.5.0"
except ImportError:
    h3_ok = False
if not h3_ok:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "h3==4.5.0"])
report_path = os.path.join(CODE, "arf_parity_report.json")
if os.path.exists(report_path):                   # never read a report left by an earlier run
    os.replace(report_path, report_path + ".previous")
r = subprocess.run([sys.executable, os.path.join(CODE, "arf_colab_parity.py"),
                    "--pipeline", os.path.join(CODE, "algorithmic_reconstruction_framework.py"),
                    "--expected", os.path.join(CODE, "arf_parity_expected_v153.json"), "--drive", ROOT],
                   capture_output=True, text=True)
print("\n".join(l for l in r.stdout.splitlines() if not l.startswith("INFO")))
if r.returncode != 0 or not os.path.exists(report_path):
    print("parity script error output (last lines):\n" + "\n".join(r.stderr.splitlines()[-25:]))
    raise RuntimeError(f"parity script failed (exit code {r.returncode}); see the error output above")
rep = json.load(open(report_path))
st = {k: v["status"] for k, v in rep["runs"].items()}
assert r.returncode == 0 and st and all(s == "MATCH" for s in st.values()), f"parity failed: {st}"
print(f"3/4 parity: {len(st)}/{len(st)} MATCH")

# 4. the protocol batch
if STAGE != "batch":
    print("4/4 skipped (STAGE = 'check'). Set STAGE = 'batch' and run again to start the 56 protocol runs.")
else:
    try:
        from google.colab import userdata
        missing = [k for k in ("SPACETRACK_USERNAME", "SPACETRACK_PASSWORD") if not userdata.get(k)]
    except Exception as e:
        missing = [f"Colab secrets unavailable ({type(e).__name__})"]
    assert not missing, f"missing Space-Track secrets: {missing}"
    print("4/4 starting the batch: PRIMARY (8), SECONDARY (32), EXPLORATORY (16)")
    runpy.run_path(os.path.join(CODE, "arf_batch_runner.py"), run_name="__main__")
