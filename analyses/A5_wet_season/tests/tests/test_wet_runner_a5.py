"""Test of arf_wetseason_runner_a5.py with the REAL arf_batch_runner.py (1.5.3), the REAL A4 runner (its
isolation helpers) and the REAL A5 prefetch module (cache file names), and a stub pipeline whose output-folder
code is Cell 0's Colab block of release 1.5.3 verbatim; a fake google.colab stands for Colab.
  - 16 runs (2 cloud settings x 8 windows), each in ARF_wetseason_A5/<variant>/<window>
  - nothing under /content; confirmatory and A4 folders unchanged; resumable
  - refuses: missing element-set cache, missing EarthEnv file, a changed release file, a changed A4 runner;
    stops if a run touches an A4 sensitivity folder
Usage: test_wet_runner_a5.py RUNNER_A5 PREFETCH_A5 RUNNER_A4 BATCH_RUNNER_V153 CELL0_BLOCK"""
import sys, os, json, types, shutil, tempfile, importlib.util, io, contextlib, runpy

RUN5, PF5, RUN4, BR, BLOCK = sys.argv[1:6]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


fake = types.ModuleType("google.colab"); fake.drive = types.SimpleNamespace(mount=lambda *a, **k: None)
g = types.ModuleType("google"); g.colab = fake
sys.modules["google"], sys.modules["google.colab"] = g, fake
root = tempfile.mkdtemp(prefix="a5run_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
cache, clim = os.path.join(batch, "shared_cache"), os.path.join(root, "ARF_climatology")
for x in (code, cache, clim):
    os.makedirs(x)
shutil.copy(BR, os.path.join(code, "arf_batch_runner.py"))
shutil.copy(RUN4, os.path.join(code, "arf_sensitivity_runner_a4.py"))
shutil.copy(PF5, os.path.join(code, "arf_colab_prefetch_wet_a5.py"))
STUB = '''__version__ = "1.5.3"
import os, json
cfg = dict(ARF_CONFIG_OVERRIDES)
ACTIVE_ROI = cfg["active_roi"]
FLEET_TAG = cfg["fleet"].replace("+", "_") + "-" + cfg["planet_swath_case"] + "_from" + cfg["analysis_epoch"][:10].replace("-", "")
''' + open(BLOCK).read() + '''
with open(os.path.join(OUT_DIR, "arf_reproducibility_manifest.json"), "w") as f:
    json.dump({"framework_version": __version__, "config_overrides": cfg}, f)
with open(os.path.join(OUT_DIR, "arf_result.json"), "w") as f:
    json.dump({"cloud": cfg.get("cloud_probability_source", "constant")}, f)
_ARF_CELLS_RUN = {"1": __version__}
'''
stub = os.path.join(code, "algorithmic_reconstruction_framework.py"); open(stub, "w").write(STUB)
spec = importlib.util.spec_from_file_location("r5", RUN5); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
r5.SHA_153 = r5.sha(stub)
with contextlib.redirect_stdout(io.StringIO()):
    pf = runpy.run_path(PF5, run_name="x")
for w in r5.WINDOWS:
    open(pf["history_path"](cache, {**pf["DESIGN"], "analysis_epoch": w}), "w").write("{}")
for m in (1, 2):
    open(os.path.join(clim, f"MODCF_monthlymean_{m:02d}.tif"), "w").write("x")
import pandas as pd
conf = []
for i in range(3):
    d = os.path.join(root, f"ARF_Outputs_conf{i}"); os.makedirs(d)
    json.dump({"framework_version": "1.5.3"}, open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    conf.append({"run_id": str(i), "status": "done", "out_dir": d})
pd.DataFrame(conf).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
a4dir = os.path.join(root, "ARF_sensitivity_A4", "k10", "2026-08-03"); os.makedirs(a4dir)
open(os.path.join(a4dir, "arf_onset_robustness.csv"), "w").write("a,b\n1,2\n")
wet = os.path.join(root, "ARF_wetseason_A5")


def run(**k):
    with contextlib.redirect_stdout(io.StringIO()):
        r5.main(root=root, code=code, batch=batch, wet_root=wet, **k)


shutil.rmtree("/content", ignore_errors=True)
run()
idx = pd.read_csv(os.path.join(batch, "wetseason_v153_a5_index.csv"), dtype=str)
check(len(idx) == 16 and (idx.status == "done").all(), "16 runs done")
ok = all(json.load(open(os.path.join(wet, v, w[:10], "arf_result.json")))["cloud"] ==
         ("earthenv" if v == "earthenv" else "constant") for v in r5.VARIANTS for w in r5.WINDOWS)
check(ok, "each run in ARF_wetseason_A5/<variant>/<window> with its own cloud setting")
check(not os.path.exists("/content"), "nothing written under /content")
check(open(os.path.join(a4dir, "arf_onset_robustness.csv")).read() == "a,b\n1,2\n", "A4 folder unchanged")
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    r5.main(root=root, code=code, batch=batch, wet_root=wet)
check(buf.getvalue().count("done earlier, skipped") == 16, "second call: all 16 skipped")


def refused(fn, text):
    try:
        fn(); return False
    except AssertionError as e:
        return text in str(e)


hp = pf["history_path"](cache, {**pf["DESIGN"], "analysis_epoch": r5.WINDOWS[3]})
os.rename(hp, hp + ".x")
check(refused(run, "element sets not cached for 1 windows"), "a missing element-set cache is refused")
os.rename(hp + ".x", hp)
os.rename(os.path.join(clim, "MODCF_monthlymean_02.tif"), os.path.join(clim, "x.tif"))
check(refused(run, "month(s) [2]"), "a missing EarthEnv file is refused")
os.rename(os.path.join(clim, "x.tif"), os.path.join(clim, "MODCF_monthlymean_02.tif"))
check(refused(lambda: run(sha_a4="0" * 64), "differs from the A4 release"), "a changed A4 runner is refused")
bad = os.path.join(code, "bad.py")
open(bad, "w").write(STUB + f'\nopen({os.path.join(a4dir, "arf_onset_robustness.csv")!r}, "a").write("3,4\\n")\n')
os.remove(os.path.join(batch, "wetseason_v153_a5_index.csv"))
r5.SHA_153 = r5.sha(bad)
check(refused(lambda: run(pipeline=bad), "STOPPED"), "a run that touches an A4 folder stops the batch")
r5.SHA_153 = "1" * 64
check(refused(run, "not release 1.5.3"), "a release file other than 1.5.3 is refused")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
