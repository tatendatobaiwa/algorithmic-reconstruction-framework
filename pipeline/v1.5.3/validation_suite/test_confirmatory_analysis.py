"""arf_confirmatory_analysis.py (v1.5.2, v1.5.3) on synthetic run folders: designs are
matched on ROI, fleet, swath case AND exclude_satellites (a run without
Sentinel-2A is never pooled with the full-fleet runs); primary design at
alpha/3, others at alpha; runs of another release are excluded; the lower
bound equals an independent Student-t computation; NOT EVALUABLE below
min_windows. Usage: test_confirmatory_analysis.py ANALYSIS_SCRIPT PROTOCOL_JSON"""
import sys, os, json, math, tempfile, importlib.util
import numpy as np
from scipy import stats

A, PROT = sys.argv[1:3]
spec = importlib.util.spec_from_file_location("arf_confirmatory_analysis", A)
an = importlib.util.module_from_spec(spec); spec.loader.exec_module(an)
protocol = json.load(open(PROT))
REL = protocol["release"]          # the release of record (1.5.2 after A1, 1.5.3 after A2)
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


root = tempfile.mkdtemp(prefix="confirm_")
rng = np.random.default_rng(7)


def make_run(name, cfg, release, h1a, h1b, h2):
    d = os.path.join(root, name)
    os.makedirs(d)
    json.dump({"framework_version": release, "config_overrides": cfg,
               "cell_versions": {"consistent": True}}, open(os.path.join(d, "arf_reproducibility_manifest.json"), "w"))
    json.dump({"h1_revised": {"series": {"earliest_completion_clear": {"n_onsets_evaluated": 73, "pct_within_margin": h1a}}},
               "h1b_image_economy": {"series": {"clear": {"n_onsets_evaluated": 73, "pct_success": h1b,
                                                          "median_images_ecc": 2, "median_images_fcfs": 5}}}},
              open(os.path.join(d, "arf_onset_robustness_summary.json"), "w"))
    json.dump({"h2": {"time_aware": {"evaluable": True, "n_realisations": 20,
                                     "difference_vs_openloop_expected": {"mean": h2}}},
               "cloud_model": {"model": "correlated"}}, open(os.path.join(d, "arf_cloud_robustness_summary.json"), "w"))
    json.dump({"status": "ok", "precision": 0.9, "recall": 0.8,
               "by_platform": {"sentinel-2a": {"precision": 0.2, "recall": 1.0}}},
              open(os.path.join(d, "arf_acquisition_validation_summary.json"), "w"))
    return d


dirs, truth = [], {"full": [], "noS2A": []}
for w in protocol["windows"]:
    base = {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet", "planet_swath_case": "conservative",
            "analysis_epoch": w}
    v = (float(rng.uniform(92, 100)), float(rng.uniform(90, 100)), float(rng.uniform(0.2, 1.5)))
    dirs.append(make_run("full_" + w[:10], base, REL, *v)); truth["full"].append(v)
    u = (float(rng.uniform(60, 90)), float(rng.uniform(40, 80)), float(rng.uniform(-0.5, 0.5)))
    dirs.append(make_run("noS2A_" + w[:10], {**base, "exclude_satellites": ["sentinel-2a"]}, REL, *u))
    truth["noS2A"].append(u)
dirs.append(make_run("old_release", {"active_roi": "ckgr_okavango", "fleet": "earth-resources+planet",
                                     "planet_swath_case": "conservative", "analysis_epoch": protocol["windows"][0]},
                     "1.5.0", 0.0, 0.0, -9.0))
# a design with only 3 windows: NOT EVALUABLE
for w in protocol["windows"][:3]:
    dirs.append(make_run("lz_" + w[:10], {"active_roi": "lower_zambezi", "fleet": "earth-resources+planet",
                                          "planet_swath_case": "conservative", "analysis_epoch": w},
                         REL, 99.0, 99.0, 1.0))

res, df = an.analyse(dirs, protocol)
designs = {(d["design"]["active_roi"], d["design"]["fleet"], d["design"].get("planet_swath_case"),
            tuple(d["design"].get("exclude_satellites") or ())): d for d in res["designs"]}
prim = designs[("ckgr_okavango", "earth-resources+planet", "conservative", ())]
expl = designs[("ckgr_okavango", "earth-resources+planet", "conservative", ("SENTINEL-2A",))]
check(len(prim["windows_found"]) == 8 and len(expl["windows_found"]) == 8,
      "primary and exploratory (without Sentinel-2A) designs each find their own 8 windows")
check(res["runs_excluded_release_or_mixed"] == [os.path.join(root, "old_release")],
      "the release-1.5.0 run is excluded (not pooled)")
check(abs(prim["alpha"] - protocol["alpha"] / 3) < 1e-12 and abs(expl["alpha"] - protocol["alpha"]) < 1e-12,
      f"alpha: primary {prim['alpha']:.4f} (Bonferroni over 3), exploratory {expl['alpha']:.4f}")


def lcb(x, alpha):
    x = np.asarray(x)
    return x.mean() - stats.t.ppf(1 - alpha, len(x) - 1) * x.std(ddof=1) / math.sqrt(len(x))


ok = True
for dsg, key in ((prim, "full"), (expl, "noS2A")):
    for i, h in enumerate(("H1a", "H1b", "H2")):
        x = [t[i] for t in truth[key]]
        got = [v for k, v in dsg[h].items() if k.startswith("lcb_")][0]
        ok &= abs(got - lcb(x, dsg["alpha"])) < 1e-9 and abs(dsg[h]["mean"] - np.mean(x)) < 1e-9
check(ok, "means and one-sided lower bounds equal an independent Student-t computation (both designs, 3 hypotheses)")
exp_status = []
for dsg, key in ((prim, "full"), (expl, "noS2A")):
    for i, h in enumerate(("H1a", "H1b", "H2")):
        b = lcb([t[i] for t in truth[key]], dsg["alpha"])
        want = ("SUPPORTED" if (b >= protocol["threshold_pct"] if h != "H2" else b > 0) else "NOT SUPPORTED")
        exp_status.append(dsg[h]["status"] == want)
check(all(exp_status), "statuses follow the frozen decision rules: " +
      ", ".join(f"{h} {d[h]['status']}" for d in (prim, expl) for h in ("H1a", "H1b", "H2")))
lz = designs[("lower_zambezi", "earth-resources+planet", "conservative", ())]
check(all(lz[h]["status"] == "NOT EVALUABLE" for h in ("H1a", "H1b", "H2")) and len(lz["windows_found"]) == 3,
      "a design with 3 of 8 windows (min_windows 5) is NOT EVALUABLE")
check("archive_sentinel-2a_precision" in df.columns and set(df["exclude_satellites"]) == {"", "SENTINEL-2A"},
      "window table carries per-platform archive precision and the exclusion key")
_ids = [a["id"] for a in protocol.get("amendments", [])]
check(res.get("amendments") == _ids and _ids[:1] == ["A1"],
      f"results record the protocol amendments ({', '.join(_ids)})")
n_fail = sum(not ok for ok, _ in RESULTS)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS],
          open("/home/claude/validate/logs/test_confirmatory_analysis.json", "w"), indent=1)
sys.exit(1 if n_fail else 0)
