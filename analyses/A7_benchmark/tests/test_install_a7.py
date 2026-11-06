"""Test of arf_colab_install_a7.py against a MyDrive copy with the REAL release 1.5.3 files, the REAL analysis module
 and the REAL A6 protocol: installs A7 (A6 kept) and the script; second run no-op; refuses a damaged zip, a
changed 1.5.3, a missing analysis module, and a protocol without A6.
Usage: test_install_a7.py INSTALL ZIP SHIP153 SHIP_P1 SHIP_A6"""
import sys, os, json, shutil, tempfile, hashlib, zipfile, importlib.util, io, contextlib
INST, ZIP, S153, SP1, SA6 = sys.argv[1:6]
RES = []
def check(c, m): RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
spec = importlib.util.spec_from_file_location("i", INST); inst = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(inst)
def fresh():
    root = tempfile.mkdtemp(prefix="i7_"); code = os.path.join(root, "ARF_code"); os.makedirs(code)
    for f in ("algorithmic_reconstruction_framework.py", "arf_batch_runner.py", "arf_colab_parity.py"):
        shutil.copy(os.path.join(S153, f), code)
    shutil.copy(os.path.join(S153, "arf_confirmatory_analysis.py"), code)
    shutil.copy(os.path.join(SA6, "arf_confirmatory_protocol.json"), code)
    shutil.copy(ZIP, os.path.join(root, "ARF_benchmark_A7.zip")); return root, code
def run(root, code, zp=None):
    b = io.StringIO()
    with contextlib.redirect_stdout(b):
        inst.install(root=root, code=code, zip_path=zp or os.path.join(root, "ARF_benchmark_A7.zip"), work=os.path.join(root, "w"))
    return b.getvalue()
root, code = fresh(); a5 = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
before = {f: sha(os.path.join(code, f)) for f in os.listdir(code)}
out = run(root, code)
p = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
check([a["id"] for a in p["amendments"]] == ["A1", "A2", "A3", "A4", "A5", "A6", "A7"] and p["amendments"][:6] == a5["amendments"], "A7 appended, A1-A6 unchanged")
check(json.load(open(os.path.join(code, "arf_confirmatory_protocol.a6.json"))) == a5, "A6 protocol kept")
check(all(os.path.exists(os.path.join(code, s)) for s in inst.SCRIPTS), "script installed")
check(all(sha(os.path.join(code, f)) == h for f, h in before.items() if f != "arf_confirmatory_protocol.json"), "release and analysis files untouched")
check("already installed" in run(root, code), "second run: nothing to do")
def refused(fn, t):
    try: fn(); return False
    except AssertionError as e: return t in str(e)
root, code = fresh(); bad = os.path.join(root, "bad.zip")
with zipfile.ZipFile(ZIP) as zi, zipfile.ZipFile(bad, "w") as zo:
    for n in zi.namelist(): zo.writestr(n, zi.read(n) + (b"\n#" if n == "arf_benchmark_a7.py" else b""))
check(refused(lambda: run(root, code, bad), "checksum mismatch"), "a damaged package is refused")
root, code = fresh(); open(os.path.join(code, "algorithmic_reconstruction_framework.py"), "a").write("\n")
check(refused(lambda: run(root, code), "1.5.3 changed"), "a changed release file is refused")
root, code = fresh(); os.remove(os.path.join(code, "arf_confirmatory_analysis.py"))
check(refused(lambda: run(root, code), "arf_confirmatory_analysis.py missing"), "a missing analysis module is refused")
root, code = fresh(); pp = os.path.join(code, "arf_confirmatory_protocol.json"); q = json.load(open(pp)); q["amendments"] = q["amendments"][:5]; json.dump(q, open(pp, "w"))
check(refused(lambda: run(root, code), "ARF_code protocol has"), "a protocol without A6 is refused")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
