"""Append amendment A5 (EXPLORATORY wet-season arm) to the A4 protocol, before any A5 run.
Usage: make_protocol_amendment_a5.py PROTOCOL_A4 RUNNER_A5 PREFETCH_A5 SUMMARY_A5 OUT_JSON"""
import sys, json, copy, hashlib, runpy, datetime
src_p, run_p, pf_p, sum_p, out = sys.argv[1:6]
p = json.load(open(src_p))
assert p["release"] == "1.5.3" and [a["id"] for a in p["amendments"]] == ["A1", "A2", "A3", "A4"], "input is not the A4 protocol"
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
R = runpy.run_path(run_p, run_name="x")
assert R["SHA_153"] == p["pipeline_sha256"]
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
new = copy.deepcopy(p)
new["amended_utc"] = now
new["amendments"].append({
    "id": "A5", "date_utc": now,
    "status_when_amended": ("EXPLORATORY. Written after the confirmatory results and the A4 post hoc results were read, "
                            "and before any run of the windows below existed."),
    "scope": ("Adds an exploratory arm. The confirmatory plan and its decisions are unchanged. No A5 result is a "
              "hypothesis decision; every A5 number is reported as exploratory, beside the confirmatory results."),
    "reason": ("The eight confirmatory windows (August and September 2026) lie in months whose EarthEnv climatology over "
               "the region is 4.9% (August) and 8.3% (September) cloud frequency (A4 runs), and the archived scenes of "
               "the first window were almost cloud-free (Cell 7D); the value of acting on cloud observations (H2) is therefore "
               "tested where cloud is rare. A5 repeats the primary design in eight windows of January and February "
               "2026, run retrospectively, with the protocol's constant cloud probability and with the climatology."),
    "runs": {"release": "1.5.3", "pipeline_sha256": p["pipeline_sha256"], "design": R["PRIMARY_DESIGN"] if "PRIMARY_DESIGN" in R else None,
             "windows": list(R["WINDOWS"]), "window_hours": 48, "variants": R["VARIANTS"],
             "n_runs": len(R["WINDOWS"]) * len(R["VARIANTS"]), "batch_name": R["BATCH_NAME"],
             "output_folders": "MyDrive/ARF_wetseason_A5/<variant>/<window> (never an earlier run folder)",
             "element_sets": "Space-Track gp_history as of each window start (7-day lookback), fetched by the release's own functions (arf_colab_prefetch_wet_a5.py)",
             "climatology": "EarthEnv MODCF_monthlymean_01.tif and _02.tif (Wilson & Jetz 2016; CC BY-NC 4.0)",
             "stated_limit": "the fleet is built from the current CelesTrak list, so satellites retired since the window drop out (Cell 2)"},
    "scripts": {"arf_wetseason_runner_a5.py": sha(run_p), "arf_colab_prefetch_wet_a5.py": sha(pf_p),
                "arf_wetseason_summary_a5.py": sha(sum_p)},
    "analysis": ("descriptive: per window H1a, H1b, H2 by the protocol's own window rules; over windows the mean, SD, "
                 "one-sided Student-t lower bounds at alpha and alpha/3, exact one-sided sign and Wilcoxon tests; "
                 "observed scene cloud of the archived acquisitions; beside the dry-season values"),
})
json.dump(new, open(out, "w"), indent=1)
print("A5 written", now)
