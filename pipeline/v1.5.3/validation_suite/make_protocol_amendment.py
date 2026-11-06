"""Write the v1.5.2 confirmatory protocol: the frozen protocol of 12 October 2026
with amendment A1 (14 October 2026, before any protocol window was run).
Usage: make_protocol_amendment.py PROTOCOL_V150 PIPELINE_V152 OUT_JSON"""
import sys, json, hashlib, re, copy, datetime

src_p, pipe, out = sys.argv[1:4]
p = json.load(open(src_p))
raw = open(pipe, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
rel = re.search(rb'^__version__ = "([^"]+)"', raw, re.M).group(1).decode()
assert rel == "1.5.2", rel
new = copy.deepcopy(p)
new["release"] = rel
new["pipeline_sha256"] = sha
new["amended_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
explor = [{"role": "exploratory", "active_roi": "ckgr_okavango", "fleet": fl, "planet_swath_case": "conservative",
           "exclude_satellites": ["SENTINEL-2A"]} for fl in ("earth-resources+planet", "earth-resources")]
new["designs"] = p["designs"] + explor
new["amendments"] = [{
    "id": "A1",
    "date_utc": new["amended_utc"],
    "status_when_amended": ("No window listed in this protocol had been run or analysed: the batch "
                            "(arf_batch_runner.py) had not been started. Decided before any confirmatory data existed."),
    "release_before": p["release"], "pipeline_sha256_before": p["pipeline_sha256"],
    "release_after": rel, "pipeline_sha256_after": sha,
    "reason": ("Per-cell audit of the pipeline before the confirmatory runs (14 October 2026). Release 1.5.2 "
               "corrects defects found in 1.5.0 and 1.5.1 and is the release the confirmatory runs will use; "
               "arf_confirmatory_analysis.py excludes runs of any other release."),
    "unchanged": ("hypotheses H1a, H1b, H2; thresholds; windows; evaluability rules; decision rules; alpha and the "
                  "Bonferroni adjustment; the primary and secondary designs; all settings"),
    "changes_that_can_affect_window_outcomes": [
        ("Cell 5 (1.5.1): each tile observation belongs to the one-minute strip that contains its closest approach. "
         "Up to 1.5.0 every strip also claimed an end cap beyond each end, so most Sentinel-2 observations appeared "
         "in two consecutive strips (12 October 2026, primary design: 12,110 of 50,075 rows, 24.2%, removed; in "
         "the optimistic-swath run, 12,911 of 57,467: every removed row has a twin in the adjacent strip of the "
         "same pass and no pass loses a tile; checked with real Orekit and an independent SGP4 propagation). "
         "Coverage and floor times change by at most one strip minute; image counts (H1b), quality and the "
         "candidate-event set can change."),
        ("Cell 6 (1.5.1): time-aware weights on delivery time instead of pass time (identical when every downlink "
         "wait is 0, as in all modelled runs) and floored at 1e-300 (only matters beyond about 745 time scales)."),
        ("Cell 5 (1.5.2): sourced acquisition exclusions (Sentinel-2A tandem pause, 13 January - 8 March 2026) are "
         "applied per pass time. No protocol window overlaps it, so no effect on the protocol runs."),
    ],
    "changes_that_cannot_affect_window_outcomes": [
        "Cell 7C (1.5.1): exact ILP tolerances (ILP gaps are reported, not tested)",
        "Cell 7D (1.5.2): per-platform archive lags and precision/recall (reported, not tested)",
        "Cell 8 (1.5.2): H3 overall status rule (H3 is reported, not tested, in this protocol)",
        "Cells 7, 7C (1.5.2): figures (IEEE widths, hours, H1a-gap figure)",
        "Cell 0 (1.5.2): default scenario and scenario-ledger key (the batch passes every design field explicitly; "
        "the analysis reads the per-run summary files)",
        "Cell 0 (1.5.2): output folder on Colab: a folder holding a complete run of another release is kept and the "
        "run writes to <folder>_v152 (changes where results are written, not what they are; the batch index records "
        "the folder used)",
        "Cells 1, 2, 9 (1.5.2): package pins and version records; deterministic sample tile; printouts",
    ],
    "added_designs": {
        "role": "exploratory (unadjusted alpha, reported as exploratory, not part of the confirmatory claims)",
        "designs": explor,
        "reason": ("Sentinel-2A is tracked (it has orbital elements) but flies an extension campaign with its own "
                   "acquisition scenario, and over the CKGR it images less often than its orbit allows: in the "
                   "Earth Search catalogue it imaged the CKGR every 20 days in August-October 2026, and a "
                   "predicted Sentinel-2A pass over the CKGR on 13 October 2026 has no archive scene. The full-fleet designs therefore carry an optimistic Sentinel-2A contribution; "
                   "the designs without it bracket the open-data results. Source for the campaign: Copernicus Data "
                   "Space Ecosystem news, 'Sentinel-2A Extension Campaign to Continue Until the End of 2027', "
                   "9 October 2026."),
    },
    "runner": "arf_batch_runner.py v1.5.2: PRIMARY, SECONDARY, then EXPLORATORY (batch names confirmatory_v152_*)",
    "deviation_report": (
        "Under deviations_policy this amendment is reported as a deviation from the plan frozen on "
        f"{p['frozen_utc']}, with the reason above. No result under the frozen release exists, because no window "
        "had been run. Release 1.5.0 counts most Sentinel-2 observations in two strips (the Cell 5 defect above), so "
        "it is not the release of record. If results under the frozen release are also wanted, run the PRIMARY list "
        "with release 1.5.0 (pipeline_sha256_before) BEFORE the 1.5.2 batch, and analyse those runs with the frozen "
        "protocol file: release 1.5.0 has no output-folder guard and would overwrite a 1.5.2 run of the same window, "
        "whereas 1.5.2 keeps a 1.5.0 run and writes to <folder>_v152."),
}]
new["status_when_frozen"] = p["status_when_frozen"]
json.dump(new, open(out, "w"), indent=2)
print(f"protocol {new['id']}: release {rel}, sha256 {sha[:16]}..., {len(new['designs'])} designs "
      f"({len(explor)} exploratory added), amendment A1 dated {new['amended_utc']}")
