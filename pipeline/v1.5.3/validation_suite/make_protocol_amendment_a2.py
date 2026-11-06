"""Write the v1.5.3 confirmatory protocol: the protocol as amended on 14 October
2026 (A1, release 1.5.2) with amendment A2 (18 October 2026, before any protocol
window was run): release 1.5.3, presentation only.
Usage: make_protocol_amendment_a2.py PROTOCOL_V152 PIPELINE_V153 IDENTITY_JSONL OUT_JSON"""
import sys, json, hashlib, re, copy, datetime

src_p, pipe, ident, out = sys.argv[1:5]
p = json.load(open(src_p))
assert p["release"] == "1.5.2" and [a["id"] for a in p["amendments"]] == ["A1"], "input is not the A1 protocol"
raw = open(pipe, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
rel = re.search(rb'^__version__ = "([^"]+)"', raw, re.M).group(1).decode()
assert rel == "1.5.3", rel
recs = [json.loads(l) for l in open(ident)]
assert recs and all(r.get("ok") for r in recs), "the 1.5.2 / 1.5.3 identity check did not pass for every run"
n12 = sum(r["mode"] == "oct12" for r in recs)
n15 = sum(r["mode"] == "oct15" for r in recs)
n_csv = sum(len(r["csv"]) for r in recs)
n_json = sum(v["values"] for r in recs for v in r["json"].values())
n_sr = sum(r["scenario_results_fields"] for r in recs)

new = copy.deepcopy(p)
new["release"] = rel
new["pipeline_sha256"] = sha
new["amended_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
new["amendments"].append({
    "id": "A2",
    "date_utc": new["amended_utc"],
    "status_when_amended": ("No window listed in this protocol had been run or analysed: the batch "
                            "(arf_batch_runner.py) had not been started; on 18 October 2026 the Drive held no "
                            "ARF_batch folder and no run folder of a protocol window. Decided before any "
                            "confirmatory data existed."),
    "release_before": p["release"], "pipeline_sha256_before": p["pipeline_sha256"],
    "release_after": rel, "pipeline_sha256_after": sha,
    "reason": ("Evaluation of the five complete Colab runs of 15 October 2026 (release 1.5.2) found presentation "
               "defects only: Figure 5 empty when a window cannot reach 99% of the region; Figures 3 and 9 "
               "unreadable (curves pressed against an axis); Figure 7 7.75 in wide, above two IEEE columns; the "
               "budget panels (Figures 3 and 4) empty without a word when a window offers fewer than five images; "
               "a run's report headline pooling every scenario run in the same Colab kernel; superseded "
               "hypotheses printed without a label. Release 1.5.3 corrects only these."),
    "unchanged": ("hypotheses H1a, H1b, H2; thresholds; windows; evaluability rules; decision rules; alpha and the "
                  "Bonferroni adjustment; the primary, secondary and exploratory designs; all settings; every "
                  "computed value"),
    "changes_that_can_affect_window_outcomes": [],
    "changes_that_cannot_affect_window_outcomes": [
        "Cell 7 (1.5.3): Figure 3 and composite panel (c) span their data with the legend under the axes; Figures "
        "3 and 4 state why they are empty when the window offers fewer images than the smallest budget swept "
        "(k = 5); whole-number k ticks; Figure 7 panel (b) title shortened, its cloud settings moved to the footer",
        "Cell 7C (1.5.3): Figure 5 at the H1a level of reachable tiles (was 99% of all ROI tiles), hours on a linear "
        "axis; Figure 9 in minutes on a linear axis with each selector's largest gap; Figure 8 drawn 7.0 in wide",
        "Cells 7, 7C, 8, 9 (1.5.3): superseded hypotheses (original H1, its time-aware variant, original H2 ratio "
        "test) labelled in the console, the report, the JSON (superseded flags) and the pipeline log",
        "Cell 8 (1.5.3): each run's own statuses printed first and saved as hypothesis_validation.this_run; "
        "ledger-wide statuses labelled as pooled when a kernel holds more than one scenario (the batch runs each "
        "configuration in a fresh namespace, so its runs hold one scenario each)",
        "All cells (1.5.3): version strings; Cell 0 workflow note",
    ],
    "evidence_no_value_changed": (
        f"Cells 6-8 of 1.5.2 and of 1.5.3 replayed on the same saved inputs, {n12} runs of 12 October 2026 and "
        f"{n15} complete Colab runs of 15 October 2026: {n_csv} CSV files identical (byte for byte, except the "
        f"solver's wall-clock column ilp_seconds of arf_ilp_gap.csv, which differs from run to run whatever the "
        f"release); {n_json} JSON values "
        f"identical apart from the keys 1.5.3 adds (hypothesis_validation.this_run; superseded flags); {n_sr} "
        "scenario_results fields identical; the report differs in label lines only; the figures 1.5.3 does not "
        "touch are byte-identical PNG files (compare_v152_v153.py; results/arf_v152_v153_identity.csv). Cells 0-5 "
        "differ from 1.5.2 in version strings and comments only, and with real Orekit reproduce Colab's refined "
        "events of the five 15 October runs row for row."),
    "runner": "arf_batch_runner.py v1.5.3: PRIMARY, SECONDARY, then EXPLORATORY (batch names confirmatory_v153_*)",
    "deviation_report": ("Under deviations_policy this amendment is reported as a deviation from the plan as amended "
                         "on 2026-10-14 (A1), with the reason above. No protocol window has been run under release "
                         "1.5.2. Release 1.5.3 computes every value as 1.5.2 does, so results under either release "
                         "are the same; arf_confirmatory_analysis.py takes the release of record from this file and "
                         "excludes runs of any other release."),
})
json.dump(new, open(out, "w"), indent=1, ensure_ascii=False)
print(f"{out}: release {rel}, sha256 {sha[:16]}..., amendments {[a['id'] for a in new['amendments']]}")
