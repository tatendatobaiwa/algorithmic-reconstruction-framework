"""Write amendment A4 to the confirmatory protocol: POST HOC sensitivity runs and
manuscript evidence, recorded AFTER the confirmatory results were read. Nothing in
the confirmatory plan changes (release of record 1.5.3, designs, hypotheses,
windows, rules); A4 only registers what was added afterwards and why, so that a
reader can tell the post hoc analyses from the pre-registered ones.

Usage: make_protocol_amendment_a4.py PROTOCOL_A3 RUNNER_A4 SUMMARY_A4 EVIDENCE EXACT_ECONOMY OUT_JSON
"""
import sys, json, hashlib, copy, runpy, datetime

src_p, runner_p, summ_p, ev_p, ex_p, out = sys.argv[1:7]
p = json.load(open(src_p))
assert p["release"] == "1.5.3" and [a["id"] for a in p["amendments"]] == ["A1", "A2", "A3"], "input is not the A3 protocol"
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
R = runpy.run_path(runner_p, run_name="arf_sensitivity_runner_a4")
assert R["SHA_153"] == p["pipeline_sha256"], "the runner does not run the release of record"
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
new = copy.deepcopy(p)
new["amended_utc"] = now
new["amendments"].append({
    "id": "A4",
    "date_utc": now,
    "status_when_amended": (
        "POST HOC. Written after the confirmatory analysis had been run and read (primary design: H1a SUPPORTED, "
        "H1b NOT SUPPORTED, H2 SUPPORTED at alpha/3) and after the post hoc follow-ups of the real-fleet runs "
        "(arf_eo_followups.py). Everything registered here was chosen with those results known."),
    "scope": ("Registers post hoc sensitivity runs and read-only manuscript evidence. The confirmatory plan and its "
              "decisions are unchanged: release of record 1.5.3, designs, hypotheses, windows, thresholds, "
              "evaluability and decision rules, alpha and its adjustment. No A4 result can change a confirmatory "
              "decision, and every A4 result is to be reported as post hoc."),
    "reason": (
        "Points a referee is expected to raise. (1) H1a is met at every onset with the protocol's budget of 30 "
        "images because earliest-completion cover (ECC) searches from the physical floor, so H1a shows that the "
        "budget suffices; smaller budgets show where it starts to bind. (2) The archive validation (Cell 7D) "
        "recorded near-zero scene cloud cover in these dry-season windows, against the protocol's constant cloud "
        "probability of 0.35; the per-tile climatology option of the release shows H2 under a realistic "
        "probability. (3) H1b failed in one window (14 September) through ties with first-come-first-served; the exact minimum number of images at each onset shows whether any selector could have used fewer. (4) The exact optimum was measured on 6-h sub-problems only; one full-window instance per "
        "window measures ECC's greedy cover against the exact minimum at full scale and repeats the confirmatory run "
        "otherwise (determinism check)."),
    "runs": {
        "release": "1.5.3", "pipeline_sha256": p["pipeline_sha256"],
        "design": R["PRIMARY_DESIGN"], "windows": p["windows"],
        "variants": R["VARIANTS"],
        "n_runs": len(R["build_runs"](p["windows"])),
        "batch_name": R["BATCH_NAME"],
        "output_folders": "MyDrive/ARF_sensitivity_A4/<variant>/<window> (never the confirmatory folders)",
        "runner": "arf_sensitivity_runner_a4.py", "runner_sha256": sha(runner_p)},
    "analyses": {
        "arf_sensitivity_summary_a4.py": {"sha256": sha(summ_p),
                                          "content": "budget, cloud-probability and full-window optimum tables"},
        "arf_exact_economy_a4.py": {"sha256": sha(ex_p),
                                    "content": ("at every onset of every primary window (clear sky), the exact minimum "
                                                "number of images, among those delivered by ECC's completion time, "
                                                "that reach the H1a level (integer programme), against ECC's and "
                                                "FCFS's image counts; from a replay of Cells 6-8 that must reproduce "
                                                "the run's arf_onset_robustness.csv byte for byte")},
        "arf_paper1_evidence.py": {"sha256": sha(ev_p),
                                   "content": ("read-only evidence from the confirmatory primary runs: selectors at "
                                               "the H1a level, the existing 6-h exact optima, the four-plan H2 bracket "
                                               "and cloud-scale sweep, archive validation per platform with observed "
                                               "scene cloud, and exact one-sided sign and Wilcoxon tests over windows "
                                               "(complements to the protocol's Student-t bounds)")}},
    "revisions": [
        ("2026-10-26 (no A4 result had been analysed): arf_exact_economy_a4.py stopped on window 1 in Colab because "
         "its replay (arf_colab_parity.replay) rebuilds the fleet's altitude and orbital period from the rounded "
         "period, about 1e-13 min off the run's; the replay now takes both from the run's arf_cell2_summary.json, "
         "and a replay that is not byte-identical is accepted only if the columns the analysis uses are exactly "
         "equal and every other column agrees to a relative 1e-9, with the differences printed")],
    "reporting_rule": ("Every A4 number is labelled post hoc in the thesis and in any manuscript, and is reported "
                       "beside, never instead of, the confirmatory decisions."),
})
json.dump(new, open(out, "w"), indent=1)
print(f"A4 written to {out} ({now})")
