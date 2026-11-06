"""Tests of arf_sensitivity_runner_a4.py with the REAL arf_batch_runner.py (1.5.3) and a stub pipeline
whose output-folder code is Cell 0's Colab block of release 1.5.3, verbatim (cell0_colab_block_v153.txt).
A fake google.colab module (drive.mount = no-op) stands for Colab.
  1 negative control: without the runner's guard the stub writes to /content/drive/MyDrive/ARF_Outputs_...
    (the confirmatory folder of the window)
  2 with the runner: every run writes to ARF_sensitivity_A4/<variant>/<window>, the confirmatory folders
    are unchanged, google.colab and the working directory are restored, the index records 48 runs
  3 resumable: a second call skips all 48
  4 a pipeline that touches a confirmatory folder stops the batch
  5 run list: 6 variants x 8 windows, primary design, one setting changed per variant
Usage: test_sensitivity_runner.py RUNNER BATCH_RUNNER_V153 CELL0_BLOCK"""
import sys, os, json, types, shutil, tempfile, hashlib, importlib.util, glob

RUNNER, BR, BLOCK = sys.argv[1:4]
RES = []


def check(c, m):
    RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m, flush=True)


fake = types.ModuleType("google.colab")
fake.drive = types.SimpleNamespace(mount=lambda *a, **k: None)
g = types.ModuleType("google"); g.colab = fake
sys.modules["google"], sys.modules["google.colab"] = g, fake

root = tempfile.mkdtemp(prefix="a4test_")
code, batch = os.path.join(root, "ARF_code"), os.path.join(root, "ARF_batch")
os.makedirs(code); os.makedirs(os.path.join(batch, "shared_cache"))
shutil.copy(BR, os.path.join(code, "arf_batch_runner.py"))
block = open(BLOCK).read()
STUB = '''__version__ = "1.5.3"
import os, json
cfg = dict(ARF_CONFIG_OVERRIDES)
ACTIVE_ROI = cfg["active_roi"]
FLEET_TAG = cfg["fleet"].replace("+", "_") + "-" + cfg["planet_swath_case"] + "_from" + cfg["analysis_epoch"][:10].replace("-", "")
''' + block + '''
with open(os.path.join(OUT_DIR, "arf_reproducibility_manifest.json"), "w") as f:
    json.dump({"framework_version": __version__, "config_overrides": cfg}, f)
with open(os.path.join(OUT_DIR, "arf_result.json"), "w") as f:
    json.dump({"k": cfg.get("max_selection_steps", 30)}, f)
_ARF_CELLS_RUN = {"1": __version__}
'''
stub = os.path.join(code, "algorithmic_reconstruction_framework.py")
open(stub, "w").write(STUB)

spec = importlib.util.spec_from_file_location("runner", RUNNER)
rn = importlib.util.module_from_spec(spec); spec.loader.exec_module(rn)
rn.SHA_153 = rn.sha(stub)                      # the stub stands for the release in this test
R = __import__("runpy").run_path(os.path.join(code, "arf_batch_runner.py"), run_name="x")
W = R["PROTOCOL_WINDOWS"]
for w in W:
    open(os.path.join(batch, "shared_cache", f"arf_gp_history_x_{w[:10]}.json"), "w").write("{}")

# 5 run list
runs = rn.build_runs(W)
check(len(runs) == 48 and len({R["run_id"](c) for _, _, c in runs}) == 48, "48 distinct runs (6 variants x 8 windows)")
check(all({k: c[k] for k in rn.PRIMARY_DESIGN} == rn.PRIMARY_DESIGN for _, _, c in runs), "every run is the primary design")
check(all(set(c) - set(rn.PRIMARY_DESIGN) - {"analysis_epoch"} == set(rn.VARIANTS[v]) for v, _, c in runs),
      "each run changes only its variant's settings")

# 1 negative control: the real batch runner without the guard writes into the confirmatory folder
neg_cfg = {**rn.PRIMARY_DESIGN, "analysis_epoch": W[0], "max_selection_steps": 10}
conf0 = "/content/drive/MyDrive/ARF_Outputs_ckgr_okavango_earth-resources_planet-conservative_from20260803"
shutil.rmtree("/content", ignore_errors=True)
R["run_batch"]([neg_cfg], "negative_control", pipeline_path=stub, batch_dir=batch,
               shared_cache_dir=os.path.join(batch, "shared_cache"))
check(os.path.exists(os.path.join(conf0, "arf_result.json")),
      "negative control: without the guard a k = 10 run lands in the confirmatory folder of its window")
shutil.rmtree("/content", ignore_errors=True)

# confirmatory folders (as the real batch made them: in the index, status done)
import pandas as pd
rows = []
for w in W:
    d = os.path.join(root, "ARF_Outputs_conf_" + w[:10])
    os.makedirs(d)
    json.dump({"framework_version": "1.5.3", "config_overrides": {}}, open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    open(os.path.join(d, "fig.png"), "wb").write(b"\x89PNG" + os.urandom(64))
    rows.append({"run_id": w, "status": "done", "out_dir": d})
pd.DataFrame(rows).to_csv(os.path.join(batch, "confirmatory_v153_primary_index.csv"), index=False)
before = rn.snapshot(rn.confirmatory_folders(batch))
cwd0 = os.getcwd()
sens = os.path.join(root, "ARF_sensitivity_A4")

# 2 the guarded batch
rn.main(root=root, code=code, batch=batch, sens_root=sens)
idx = pd.read_csv(os.path.join(batch, "sensitivity_v153_a4_index.csv"), dtype=str)
check(len(idx) == 48 and (idx.status == "done").all(), "index: 48 runs done")
ok = all(os.path.normpath(r.out_dir) == os.path.normpath(os.path.join(sens, v, w[:10]))
         for (v, w, c), r in zip(sorted(runs, key=lambda x: R["run_id"](x[2])),
                                 idx.sort_values("run_id").itertuples()))
check(ok, "every run wrote to ARF_sensitivity_A4/<variant>/<window> (index out_dir)")
check(all(json.load(open(os.path.join(sens, v, w[:10], "arf_result.json")))["k"] == c.get("max_selection_steps", 30)
          for v, w, c in runs), "each folder holds its own variant's result")
check(not os.path.exists("/content"), "nothing written under /content/drive")
check(rn.changed(before, rn.confirmatory_folders(batch)) == [], "confirmatory folders unchanged")
check(sys.modules.get("google.colab") is fake and os.getcwd() == cwd0, "google.colab and working directory restored")

# 3 resumable
import io, contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rn.main(root=root, code=code, batch=batch, sens_root=sens)
check(buf.getvalue().count("done earlier, skipped") == 48, "second call: all 48 skipped")

# 4 a pipeline that touches a confirmatory folder stops the batch
bad = STUB + f'\nopen({os.path.join(rows[3]["out_dir"], "arf_result_extra.json")!r}, "w").write("x")\n'
bad_path = os.path.join(code, "bad_pipeline.py"); open(bad_path, "w").write(bad)
os.remove(os.path.join(batch, "sensitivity_v153_a4_index.csv"))   # fresh index, so the runs execute
rn.SHA_153 = rn.sha(bad_path)
try:
    rn.main(root=root, code=code, batch=batch, sens_root=os.path.join(root, "sens_bad"), only={"k05"}, pipeline=bad_path)
    stopped = False
except AssertionError as e:
    stopped = "STOPPED" in str(e)
check(stopped, "a run that writes into a confirmatory folder stops the batch")

# release guard
rn.SHA_153 = "0" * 64
try:
    rn.main(root=root, code=code, batch=batch, sens_root=sens); g_ok = False
except AssertionError as e:
    g_ok = "not release 1.5.3" in str(e)
check(g_ok, "refuses a pipeline file other than release 1.5.3")
print(f"{sum(RES)}/{len(RES)} passed")
sys.exit(0 if all(RES) else 1)
