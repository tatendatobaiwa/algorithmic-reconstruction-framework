"""Cell 0 output folder on Colab (v1.5.2 guard): a complete run of another
release in the scenario's folder is kept and this run writes to
<folder>_v<release digits>; a folder holding a run of the same release, or no
manifest, is reused. Colab is simulated: a stub google.colab.drive whose
mount() does nothing, and a scratch /content/drive/MyDrive (removed after).
Usage: test_cell0_outdir.py PIPELINE"""
import sys, os, io, json, shutil, hashlib, tempfile, contextlib
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
C = e2e15.cells_of(open(P, encoding="utf-8").read())
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


assert not os.path.exists("/content"), "/content exists: refusing to touch it"
stub = tempfile.mkdtemp(prefix="colabstub_")
os.makedirs(os.path.join(stub, "google", "colab"))
open(os.path.join(stub, "google", "__init__.py"), "w").close()
open(os.path.join(stub, "google", "colab", "__init__.py"), "w").close()
with open(os.path.join(stub, "google", "colab", "drive.py"), "w") as f:
    f.write("def mount(path, *a, **k):\n    pass\n")
sys.path.insert(0, stub)
ROOT = "/content/drive/MyDrive"
os.makedirs(ROOT)


def cell0():
    for m in [m for m in sys.modules if m == "google" or m.startswith("google.")]:
        del sys.modules[m]
    ns = {"__name__": "__main__"}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(C["0"], "cell0", "exec"), ns)
    return ns, buf.getvalue()


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


try:
    ns, out = cell0()
    base = f"{ROOT}/ARF_Outputs_{ns['ACTIVE_ROI']}_{ns['FLEET_TAG']}"
    rel = ns["__version__"]
    suffixed = base + "_v" + rel.replace(".", "")
    check(ns["IN_COLAB"] and ns["OUT_DIR"] == base and os.path.exists(os.path.join(base, "arf_config.json")),
          f"no earlier run: writes to {os.path.basename(base)}")

    man = os.path.join(base, "arf_reproducibility_manifest.json")
    json.dump({"framework_version": "1.5.0"}, open(man, "w"))
    open(os.path.join(base, "arf_refined_events.csv"), "w").write("sat,tile\nX,1\n")
    before = {f: sha(os.path.join(base, f)) for f in os.listdir(base)}
    ns, out = cell0()
    after = {f: sha(os.path.join(base, f)) for f in os.listdir(base)}
    check(ns["OUT_DIR"] == suffixed and os.path.exists(os.path.join(suffixed, "arf_config.json"))
          and before == after and "release 1.5.0" in out,
          f"complete run of release 1.5.0 present: kept byte for byte; this run writes to {os.path.basename(suffixed)}")

    ns2, _ = cell0()
    check(ns2["OUT_DIR"] == suffixed, "repeat run of this release with the 1.5.0 run still present: same suffixed folder")

    json.dump({"framework_version": rel}, open(man, "w"))
    ns, _ = cell0()
    check(ns["OUT_DIR"] == base, f"complete run of this release ({rel}): folder reused, as before")

    open(man, "w").write("{not json")
    ns, out = cell0()
    check(ns["OUT_DIR"] == suffixed and "unreadable manifest" in out, "unreadable manifest: treated as another release")

    os.remove(man)
    ns, _ = cell0()
    check(ns["OUT_DIR"] == base, "folder without a manifest (incomplete run): reused, as before")
finally:
    shutil.rmtree("/content")
    shutil.rmtree(stub)
check(not os.path.exists("/content"), "scratch /content removed")
n_fail = sum(not ok for ok, _ in RESULTS)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump({"pipeline": P, "results": [{"ok": ok, "check": m} for ok, m in RESULTS]},
          open("/home/claude/validate/logs/test_cell0_outdir.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
