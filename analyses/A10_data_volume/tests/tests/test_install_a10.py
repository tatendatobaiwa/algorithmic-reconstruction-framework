"""Test of arf_colab_install_a10.py against a MyDrive copy with the REAL release 1.5.3 files, the REAL analysis module
 and the REAL A8 protocol: installs A9 (A8 kept) and the script; second run no-op; refuses a damaged zip, a
changed 1.5.3, and a protocol without A8.
Usage: test_install_a9.py INSTALL ZIP SHIP153 SHIP_P1 SHIP_A8"""
import sys, os, json, shutil, tempfile, hashlib, zipfile, importlib.util, io, contextlib
INST, ZIP, S153, SP1, SA8 = sys.argv[1:6]
RES = []
def check(c, m): RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
spec = importlib.util.spec_from_file_location("i", INST); inst = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(inst)
def fresh():
    root = tempfile.mkdtemp(prefix="i10_"); code = os.path.join(root, "ARF_code"); os.makedirs(code)
    for f in ("algorithmic_reconstruction_framework.py", "arf_batch_runner.py", "arf_colab_parity.py"):
        shutil.copy(os.path.join(S153, f), code)
    shutil.copy(os.path.join(S153, "arf_confirmatory_analysis.py"), code)
    shutil.copy(os.path.join(SP1, "arf_exact_economy_a4.py"), code)
    shutil.copy(os.path.join(SA8, "arf_confirmatory_protocol.json"), code)
    shutil.copy(ZIP, os.path.join(root, "ARF_datavolume_A10.zip")); return root, code
def run(root, code, zp=None):
    b = io.StringIO()
    with contextlib.redirect_stdout(b):
        inst.install(root=root, code=code, zip_path=zp or os.path.join(root, "ARF_datavolume_A10.zip"), work=os.path.join(root, "w"))
    return b.getvalue()
root, code = fresh(); a5 = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
before = {f: sha(os.path.join(code, f)) for f in os.listdir(code)}
out = run(root, code)
p = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
check([a["id"] for a in p["amendments"]] == ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"] and p["amendments"][:9] == a5["amendments"], "A10 appended, A1-A9 unchanged")
check(json.load(open(os.path.join(code, "arf_confirmatory_protocol.a9.json"))) == a5, "A9 protocol kept")
check(all(os.path.exists(os.path.join(code, s)) for s in inst.SCRIPTS), "script installed")
check(all(sha(os.path.join(code, f)) == h for f, h in before.items() if f != "arf_confirmatory_protocol.json"), "release and analysis files untouched")
check("already installed" in run(root, code), "second run: nothing to do")
def refused(fn, t):
    try: fn(); return False
    except AssertionError as e: return t in str(e)
root, code = fresh(); bad = os.path.join(root, "bad.zip")
with zipfile.ZipFile(ZIP) as zi, zipfile.ZipFile(bad, "w") as zo:
    for n in zi.namelist(): zo.writestr(n, zi.read(n) + (b"\n#" if n == "arf_data_volume_a10.py" else b""))
check(refused(lambda: run(root, code, bad), "checksum mismatch"), "a damaged package is refused")
root, code = fresh(); open(os.path.join(code, "algorithmic_reconstruction_framework.py"), "a").write("\n")
check(refused(lambda: run(root, code), "1.5.3 changed"), "a changed release file is refused")
root, code = fresh(); os.remove(os.path.join(code, "arf_exact_economy_a4.py"))
check(refused(lambda: run(root, code), "arf_exact_economy_a4.py missing"), "a missing A4 exact-economy script is refused")
root, code = fresh(); pp = os.path.join(code, "arf_confirmatory_protocol.json"); q = json.load(open(pp)); q["amendments"] = q["amendments"][:8]; json.dump(q, open(pp, "w"))
check(refused(lambda: run(root, code), "ARF_code protocol has"), "a protocol without A9 is refused")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
