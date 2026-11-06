"""Test of arf_colab_install_a4.py against a MyDrive copy holding the REAL release files (1.5.3, 1.5.4, the
A3 protocol, arf_batch_runner.py, arf_confirmatory_analysis.py, arf_colab_parity.py):
  installs the A4 protocol (A3 kept), the four scripts; second run is a no-op; refuses a damaged zip,
  a changed 1.5.3 file, and a protocol that is not A3; leaves 1.5.3 and 1.5.4 untouched.
Usage: test_install_a4.py INSTALL ZIP SHIP153 SHIP154 [FIRST_A4_ZIP]"""
import sys, os, json, shutil, tempfile, hashlib, zipfile, importlib.util, io, contextlib

INST, ZIP, S153, S154 = sys.argv[1:5]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def fresh():
    root = tempfile.mkdtemp(prefix="inst_")
    code = os.path.join(root, "ARF_code"); os.makedirs(code)
    for f in ("algorithmic_reconstruction_framework.py", "arf_batch_runner.py", "arf_confirmatory_analysis.py",
              "arf_colab_parity.py"):
        shutil.copy(os.path.join(S153, f), code)
    shutil.copy(os.path.join(S154, "algorithmic_reconstruction_framework_v154.py"), code)
    shutil.copy(os.path.join(S154, "arf_confirmatory_protocol.json"), code)
    shutil.copy(ZIP, os.path.join(root, "ARF_paper1_A4.zip"))
    return root, code


spec = importlib.util.spec_from_file_location("inst", INST)
inst = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(inst)


def run(root, code, zp=None):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        inst.install(root=root, code=code, zip_path=zp or os.path.join(root, "ARF_paper1_A4.zip"),
                     work=os.path.join(root, "work"))
    return buf.getvalue()


root, code = fresh()
before = {f: sha(os.path.join(code, f)) for f in os.listdir(code)}
a3 = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
out = run(root, code)
p = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
check([a["id"] for a in p["amendments"]] == ["A1", "A2", "A3", "A4"], "A4 protocol installed")
check(json.load(open(os.path.join(code, "arf_confirmatory_protocol.a3.json"))) == a3, "A3 protocol kept as .a3.json")
check(all(os.path.exists(os.path.join(code, s)) for s in inst.SCRIPTS), "four scripts installed")
check(all(sha(os.path.join(code, f)) == h for f, h in before.items() if f != "arf_confirmatory_protocol.json"),
      "1.5.3, 1.5.4 and the other release files untouched")
check(p["amendments"][3]["runs"]["pipeline_sha256"] == inst.SHA_153, "A4 registers runs of release 1.5.3")
out2 = run(root, code)
check("already installed" in out2, "second run: nothing to do")

# damaged zip
root, code = fresh()
bad = os.path.join(root, "bad.zip")
with zipfile.ZipFile(ZIP) as zi, zipfile.ZipFile(bad, "w") as zo:
    for n in zi.namelist():
        data = zi.read(n)
        zo.writestr(n, data + b"\n# x" if n == "arf_sensitivity_runner_a4.py" else data)
try:
    run(root, code, bad); ok = False
except AssertionError as e:
    ok = "checksum mismatch" in str(e)
check(ok, "a damaged package is refused")
check([a["id"] for a in json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))["amendments"]][-1] == "A3",
      "... and nothing installed")
# changed 1.5.3
root, code = fresh()
open(os.path.join(code, "algorithmic_reconstruction_framework.py"), "a").write("\n")
try:
    run(root, code); ok = False
except AssertionError as e:
    ok = "1.5.3 changed" in str(e)
check(ok, "a changed release 1.5.3 file is refused")
# protocol not A3
root, code = fresh()
pp = os.path.join(code, "arf_confirmatory_protocol.json")
q = json.load(open(pp)); q["amendments"] = q["amendments"][:2]; json.dump(q, open(pp, "w"))
try:
    run(root, code); ok = False
except AssertionError as e:
    ok = "expected A1-A3" in str(e)
check(ok, "a protocol other than A3 is refused")
# upgrade from the first version of the package (its A4) to this one
if len(sys.argv) > 5:
    root, code = fresh()
    run(root, code, sys.argv[5])
    first = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
    out = run(root, code)
    now = json.load(open(os.path.join(code, "arf_confirmatory_protocol.json")))
    check([a["id"] for a in now["amendments"]] == ["A1", "A2", "A3", "A4"] and now != first
          and json.load(open(os.path.join(code, "arf_confirmatory_protocol.a4_prev.json"))) == first
          and json.load(open(os.path.join(code, "arf_confirmatory_protocol.a3.json"))) == a3,
          "upgrade from the first A4 package: new A4 installed, first A4 and A3 kept")
    check(sha(os.path.join(code, "arf_exact_economy_a4.py")) == sha(os.path.join(root, "work", "arf_exact_economy_a4.py")),
          "upgrade: revised scripts installed")
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
