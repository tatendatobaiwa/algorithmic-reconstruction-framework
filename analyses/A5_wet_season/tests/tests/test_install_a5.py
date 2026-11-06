"""Test of arf_colab_install_a5.py against a MyDrive copy with the REAL release files, the REAL A4 protocol and
A4 runner: installs A5 (A4 kept), the three scripts; second run no-op; refuses a damaged zip, a changed 1.5.3,
and a protocol without A4. Usage: test_install_a5.py INSTALL ZIP SHIP153 SHIP_P1"""
import sys, os, json, shutil, tempfile, hashlib, zipfile, importlib.util, io, contextlib
INST, ZIP, S153, SP1 = sys.argv[1:5]
RES = []
def check(c, m): RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
spec = importlib.util.spec_from_file_location("i", INST); inst = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(inst)
def fresh():
    root = tempfile.mkdtemp(prefix="i5_"); code = os.path.join(root, "ARF_code"); os.makedirs(code)
    for f in ("algorithmic_reconstruction_framework.py", "arf_batch_runner.py", "arf_confirmatory_analysis.py"):
        shutil.copy(os.path.join(S153, f), code)
    for f in ("arf_sensitivity_runner_a4.py", "arf_confirmatory_protocol.json"):
        shutil.copy(os.path.join(SP1, f), code)
    shutil.copy(ZIP, os.path.join(root, "ARF_wetseason_A5.zip")); return root, code
def run(root, code, zp=None):
    b = io.StringIO()
    with contextlib.redirect_stdout(b):
        inst.install(root=root, code=code, zip_path=zp or os.path.join(root, "ARF_wetseason_A5.zip"), work=os.path.join(root, "w"))
    return b.getvalue()
root, code = fresh(); a4 = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
before = {f: sha(os.path.join(code, f)) for f in os.listdir(code)}
run(root, code)
p = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
check([a["id"] for a in p["amendments"]] == ["A1", "A2", "A3", "A4", "A5"] and p["amendments"][:4] == a4["amendments"], "A5 appended, A1-A4 unchanged")
check(json.load(open(os.path.join(code, "arf_confirmatory_protocol.a4.json"))) == a4, "A4 protocol kept")
check(all(os.path.exists(os.path.join(code, s)) for s in inst.SCRIPTS), "three scripts installed")
check(all(sha(os.path.join(code, f)) == h for f, h in before.items() if f != "arf_confirmatory_protocol.json"), "release and A4 files untouched")
check("already installed" in run(root, code), "second run: nothing to do")
def refused(fn, t):
    try: fn(); return False
    except AssertionError as e: return t in str(e)
root, code = fresh(); bad = os.path.join(root, "bad.zip")
with zipfile.ZipFile(ZIP) as zi, zipfile.ZipFile(bad, "w") as zo:
    for n in zi.namelist(): zo.writestr(n, zi.read(n) + (b"\n#" if n == "arf_wetseason_runner_a5.py" else b""))
check(refused(lambda: run(root, code, bad), "checksum mismatch"), "a damaged package is refused")
root, code = fresh(); open(os.path.join(code, "algorithmic_reconstruction_framework.py"), "a").write("\n")
check(refused(lambda: run(root, code), "1.5.3 changed"), "a changed release file is refused")
root, code = fresh(); pp = os.path.join(code, "arf_confirmatory_protocol.json"); q = json.load(open(pp)); q["amendments"] = q["amendments"][:3]; json.dump(q, open(pp, "w"))
check(refused(lambda: run(root, code), "ARF_code protocol has"), "a protocol without A4 is refused")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
