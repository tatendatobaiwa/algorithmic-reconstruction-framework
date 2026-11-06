# =============================================================================
# CELL 8: HYPOTHESIS VALIDATION (H1/H2/H3)
# =============================================================================
# Statuses (SUPPORTED / NOT SUPPORTED / NOT EVALUABLE) for every scenario in
# the ledger, from the analyses of Cells 6 and 7C:
#   H1 = H1a and H1b (Cell 7C onset analysis, clear sky; level = 95% of the
#        tiles reachable in each onset's horizon):
#     H1a latency non-inferiority: the earliest-completion cover (ECC, k images)
#        reaches the level within one orbital period of the physical floor at
#        >= 90% of onsets. ECC aims at the floor by construction, so H1a checks
#        that the budget allows it; it does not compare policies (the time-aware
#        selector's result is reported alongside).
#     H1b image economy: at >= 90% of the onsets where first-come-first-served
#        needs two or more images, ECC meets the H1a criterion with fewer images.
#     The onsets of one window share most of their passes, so these statuses
#     describe the window; the confirmatory test treats windows as the unit.
#   H2 Fail-Signal adaptivity (Cell 7C, cloud realisations): with the same
#        number of images, the adaptive selector covers more under cloud than
#        the same planner run in advance with no cloud observation (one-sided
#        95% lower bound of the mean paired difference > 0).
#   H3 redundancy below 1.5 when 95% coverage is first reached (clear sky);
#        NOT EVALUABLE when the passes cannot reach 95%.
#   Overall status of H1a, H1b, H2 and (from v1.5.2) H3: SUPPORTED only if
#   every evaluable scenario in the ledger is; NOT EVALUABLE if none is.
# Kept for the record: the original H1 (completion < 10 min) and the original
# H2 (cloud coverage > 50% of clear-sky coverage, which a plan with no
# Fail-Signal passes whenever p < 0.5, so it could not fail).
# Provenance: rebuilt in October 2026 when the working notebook was found to
# have lost this cell; two original criteria were corrected then (H2 tested
# coverage > 0; H3 used OR instead of AND).
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 8 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["8"] = _ARF_CELL_VERSION

logger = logging.getLogger("arf.cell8")

H1_LATENCY_THRESHOLD_MIN = 10.0
H2_MAINTENANCE_THRESHOLD_PCT = 50.0
H3_COVERAGE_THRESHOLD_PCT = 95.0
H3_REDUNDANCY_THRESHOLD = 1.5


def _scenario_label(key):
    """Format a scenario_results key -- (active_roi, active_n) -- as a
    readable label for validator assessments. Falls back to str() for
    any legacy int-only key, so this doesn't hard-crash against an
    older-format ledger loaded from a previous run."""
    if isinstance(key, tuple) and len(key) == 2:
        roi, n = key
        return f"ROI={roi}, N={n}" if isinstance(n, int) else f"ROI={roi}, fleet={n}"
    return f"N={key}"


def validate_h1_latency(scenario_results, threshold_min=H1_LATENCY_THRESHOLD_MIN,
                        latency_key='latency_clear_min', label="H1 (Latency Reduction)"):
    """H1 (Latency Reduction): SUPPORTED if latency_clear_min -- the mosaic
    completion time, i.e. when delivered coverage reaches its final value --
    is < threshold_min for at least one tested scenario."""
    assessments = {}
    any_under = False
    for scenario_k, data in scenario_results.items():
        latency = data.get(latency_key)
        under = latency is not None and latency < threshold_min
        assessments[_scenario_label(scenario_k)] = {"latency_clear_min": latency, "under_threshold": under}
        if under:
            any_under = True
    return {
        "hypothesis": label,
        "overall_status": "SUPPORTED ✅" if any_under else "NOT SUPPORTED ❌",
        "threshold_min": threshold_min,
        "assessments": assessments,
    }


def validate_h2_reliability_fail_signal(scenario_results, maintenance_threshold_pct=H2_MAINTENANCE_THRESHOLD_PCT):
    """H2 (Reliability via Fail-Signal): SUPPORTED if cloud-scenario coverage
    maintains > maintenance_threshold_pct of clear-sky coverage for at
    least one tested scenario."""
    assessments = {}
    any_maintained = False
    for scenario_k, data in scenario_results.items():
        cc = data.get('coverage_clear_pct', 0.0)
        cl = data.get('coverage_cloud_pct', 0.0)
        maintenance_pct = (cl / cc * 100) if cc > 0 else 0.0
        maintained = maintenance_pct > maintenance_threshold_pct
        assessments[_scenario_label(scenario_k)] = {
            "coverage_clear_pct": cc, "coverage_cloud_pct": cl,
            "coverage_maintenance_pct": maintenance_pct,
        }
        if maintained:
            any_maintained = True
    return {
        "hypothesis": "H2 (Reliability via Fail-Signal)",
        "overall_status": "SUPPORTED ✅" if any_maintained else "NOT SUPPORTED ❌",
        "maintenance_threshold_pct": maintenance_threshold_pct,
        "assessments": assessments,
    }


def validate_h3_redundancy_minimisation(scenario_results,
                                         coverage_threshold_pct=H3_COVERAGE_THRESHOLD_PCT,
                                         redundancy_threshold=H3_REDUNDANCY_THRESHOLD):
    """H3 (Redundancy Minimisation): SUPPORTED if coverage_clear_pct >
    coverage_threshold_pct AND h3_redundancy_clear -- the redundancy at the
    point coverage_threshold_pct is first reached (Cell 6) -- is below
    redundancy_threshold,
    TOGETHER, for at least one tested scenario.
    NOT EVALUABLE (v1.4.0) for a scenario whose passes cannot reach
    coverage_threshold_pct at all (usable clear-sky reach below it): the
    redundancy at that coverage is then undefined, whatever the selector does.
    Overall (v1.5.2): SUPPORTED only if every evaluable scenario is; NOT
    EVALUABLE if none is -- the same rule as H1a, H1b and H2 (until v1.5.1 a
    single supported scenario sufficed, the rule of the original H1/H2)."""
    assessments = {}
    any_supported = any_evaluable = False
    for scenario_k, data in scenario_results.items():
        cov = data.get('coverage_clear_pct', 0.0)
        if data.get('h3_coverage_pct') not in (None, coverage_threshold_pct):
            logger.warning(f"H3: Cell 6 measured redundancy at {data.get('h3_coverage_pct')}% coverage but "
                           f"H3_COVERAGE_THRESHOLD_PCT is {coverage_threshold_pct}%; update _H3_COVERAGE_PCT in Cell 6.")
        reach = data.get('reachable_usable_clear_pct', data.get('reachable_tiles_pct'))
        red = data.get('h3_redundancy_clear')
        evaluable = reach is None or reach >= coverage_threshold_pct - 1e-9
        supported = evaluable and red is not None and (cov > coverage_threshold_pct) and (red < redundancy_threshold)
        assessments[_scenario_label(scenario_k)] = {
            "status": ("NOT EVALUABLE" if not evaluable else ("SUPPORTED" if supported else "NOT SUPPORTED")),
            "coverage_clear_pct": cov, "reachable_clear_pct": reach,
            "redundancy_clear": float('inf') if red is None else red,
            "passes_to_coverage": data.get('h3_passes_clear'),
            "redundancy_all_passes": data.get('redundancy_clear'), "and_condition": supported,
        }
        any_evaluable |= evaluable
        any_supported |= supported
    _ev = [a["status"] == "SUPPORTED" for a in assessments.values() if a["status"] != "NOT EVALUABLE"]
    return {
        "hypothesis": "H3 (Redundancy Minimisation)",
        "overall_status": ("NOT EVALUABLE" if not _ev else ("SUPPORTED ✅" if all(_ev) else "NOT SUPPORTED ❌")),
        "overall_rule": "SUPPORTED only if every evaluable scenario is (v1.5.2)",
        "coverage_threshold_pct": coverage_threshold_pct,
        "redundancy_threshold": redundancy_threshold,
        "assessments": assessments,
    }


h1_result = validate_h1_latency(scenario_results)
h1_ta_result = validate_h1_latency(scenario_results, latency_key='ta_latency_clear_min',
                                   label="H1 (Latency Reduction) -- time-aware objective")
h2_result = validate_h2_reliability_fail_signal(scenario_results)
h3_result = validate_h3_redundancy_minimisation(scenario_results)
# v1.5.3: the superseded hypotheses carry the label wherever they are printed or
# saved. Until v1.5.2 the console and the pipeline log printed, e.g., "H2
# (Reliability via Fail-Signal): SUPPORTED" with nothing to say that this test
# was replaced (and passes whenever p < 0.5).
for _r, _note in ((h1_result, "superseded by H1a and H1b: the original < 10 min completion target, falsified "
                              "by the real-fleet runs of 2026-10-09"),
                  (h1_ta_result, "superseded by H1a and H1b: the < 10 min completion target for the "
                                 "time-aware selector"),
                  (h2_result, "superseded by H2 (v1.5.0, comparative): a plan with no Fail-Signal passes "
                              "this ratio test whenever p < 0.5")):
    _r["superseded"], _r["superseded_note"] = True, _note


def _share_status(n_eval, pct, min_pct, min_n):
    """Status of an 'at >= min_pct of onsets' test: NOT EVALUABLE below min_n
    evaluated onsets (see Cell 7C, evaluability), else SUPPORTED / NOT SUPPORTED."""
    if not n_eval or (min_n is not None and n_eval < min_n):
        return "NOT EVALUABLE"
    return "SUPPORTED" if pct is not None and pct >= min_pct - 1e-9 else "NOT SUPPORTED"


def _overall(statuses):
    """SUPPORTED only if every evaluable scenario is; NOT EVALUABLE if none is."""
    verdicts = [s == "SUPPORTED" for s in statuses if s in ("SUPPORTED", "NOT SUPPORTED")]
    return "NOT EVALUABLE" if not verdicts else ("SUPPORTED ✅" if all(verdicts) else "NOT SUPPORTED ❌")


def validate_h1_revised(scenario_results):
    """H1a, latency non-inferiority: SUPPORTED for a scenario if the primary
    selector -- the earliest-completion cover from v1.5.0 (the time-aware
    selector in v1.4.0 runs, whose summaries carry no primary_series) -- is
    within the floor margin at >= min_onset_pct of the evaluated onsets, clear
    sky. NOT EVALUABLE if fewer onsets were evaluated than min_evaluable_onsets.
    The time-aware selector's status is reported alongside. Overall: SUPPORTED
    only if every evaluable scenario is supported."""
    assessments = {}
    for k, d in scenario_results.items():
        h = (d.get("onset_robustness") or {}).get("h1_revised")
        if not h:
            assessments[_scenario_label(k)] = {"status": "NOT RUN"}
            continue
        prim = h.get("primary_series", "time_aware")
        s = h["series"][f"{prim}_clear"]
        st = _share_status(s["n_onsets_evaluated"], s["pct_within_margin"], h["min_onset_pct"],
                           h.get("min_evaluable_onsets"))
        s_ta = h["series"]["time_aware_clear"]
        st_ta = _share_status(s_ta["n_onsets_evaluated"], s_ta["pct_within_margin"], h["min_onset_pct"],
                              h.get("min_evaluable_onsets"))
        assessments[_scenario_label(k)] = {"status": st, "primary_selector": prim, "status_time_aware": st_ta,
                                           "statement": h["statement"],
                                           "floor_margin_source": h["floor_margin_source"],
                                           "min_evaluable_onsets": h.get("min_evaluable_onsets"),
                                           **{f"{n}": h["series"][n] for n in h["series"]}}
    return {"hypothesis": "H1a: latency non-inferiority (earliest-completion cover vs physical floor, reachable tiles)",
            "overall_status": _overall(a["status"] for a in assessments.values()), "assessments": assessments,
            "note": "Revised after the original < 10 min H1 was falsified by real-fleet runs (2026-10-09, then for "
                    "the time-aware selector); tested on the earliest-completion cover from v1.5.0 (2026-10-12), so "
                    "that H1a and H1b concern one selector. Confirmatory only on data not used to write it. Every "
                    "policy meets the floor margin when the budget allows (first-come-first-served by construction), "
                    "so H1a shows the selector loses no time, not that it beats the baselines; H1b carries the "
                    "comparison. Baseline rates are listed."}


def validate_h1b_image_economy(scenario_results):
    """H1b, image economy (added 2026-10-12): SUPPORTED for a scenario if, at
    >= min_onset_pct of the onsets where first-come-first-served needs two or
    more images to reach the level, the earliest-completion cover reaches it
    within the floor margin with fewer images (clear sky). NOT EVALUABLE below
    min_evaluable_onsets such onsets."""
    assessments = {}
    for k, d in scenario_results.items():
        h = (d.get("onset_robustness") or {}).get("h1b_image_economy")
        if not h:
            assessments[_scenario_label(k)] = {"status": "NOT RUN"}
            continue
        s = h["series"]["clear"]
        st = _share_status(s["n_onsets_evaluated"], s["pct_success"], h["min_onset_pct"], h.get("min_evaluable_onsets"))
        assessments[_scenario_label(k)] = {"status": st, "statement": h["statement"],
                                           "min_evaluable_onsets": h.get("min_evaluable_onsets"),
                                           "clear": h["series"]["clear"], "cloud": h["series"]["cloud"]}
    return {"hypothesis": "H1b: image economy (earliest-completion cover vs first-come-first-served)",
            "overall_status": _overall(a["status"] for a in assessments.values()), "assessments": assessments,
            "note": "Formulated after the 15 real-fleet runs of 2026-10-12: confirmatory only on fresh data."}


def combine_h1(h1a, h1b):
    """H1 as reframed in v1.4.0 = H1a AND H1b, per scenario: NOT SUPPORTED if
    either part is; SUPPORTED if both are; otherwise NOT EVALUABLE."""
    assessments = {}
    labels = list(h1a["assessments"]) + [lab for lab in h1b["assessments"] if lab not in h1a["assessments"]]
    for lab in labels:                                            # scenario order kept, as in the ledger
        a = h1a["assessments"].get(lab, {}).get("status", "NOT RUN")
        b = h1b["assessments"].get(lab, {}).get("status", "NOT RUN")
        st = ("NOT SUPPORTED" if "NOT SUPPORTED" in (a, b) else
              ("SUPPORTED" if a == b == "SUPPORTED" else "NOT EVALUABLE"))
        assessments[lab] = {"status": st, "H1a": a, "H1b": b}
    return {"hypothesis": "H1 (v1.4.0): latency non-inferiority (H1a) and image economy (H1b)",
            "overall_status": _overall(a["status"] for a in assessments.values()), "assessments": assessments}


def validate_h2_comparative(scenario_results):
    """H2 (v1.5.0), Fail-Signal adaptivity: per scenario, the status computed in
    Cell 7C (section 4) for the primary (time-aware) objective -- SUPPORTED if,
    with the same number of images, the adaptive selector covers more under
    cloud than an open-loop plan maximising expected coverage (one-sided 95%
    lower bound of the mean paired difference > 0 over the cloud
    realisations); NOT EVALUABLE when there was nothing to adapt or no room to
    gain. Coverage-first is reported alongside. Overall: SUPPORTED only if every
    evaluable scenario is supported."""
    assessments = {}
    for k, d in scenario_results.items():
        c = (d.get("cloud_robustness") or {})
        if not c.get("h2"):
            assessments[_scenario_label(k)] = {"status": "NOT RUN"}
            continue
        h, hc = c["h2"]["time_aware"], c["h2"].get("coverage_first", {})
        assessments[_scenario_label(k)] = {"status": h["status"], "status_coverage_first": hc.get("status"),
                                           "time_aware": h, "coverage_first": hc,
                                           "cloud_model": c.get("cloud_model"), "seeds": c.get("seeds")}
    return {"hypothesis": "H2 (v1.5.0): Fail-Signal adaptivity (adaptive vs open-loop plan with the same images, "
                          "over cloud realisations)",
            "overall_status": _overall(a["status"] for a in assessments.values()), "assessments": assessments,
            "note": "Replaces the original H2 (cloud coverage > 50% of clear-sky coverage), which a plan with no "
                    "Fail-Signal passes whenever p < 0.5. Formulated 2026-10-12 after the 15 runs of that day: "
                    "confirmatory only on fresh data."}


h1_revised_result = validate_h1_revised(scenario_results)
h1b_result = validate_h1b_image_economy(scenario_results)
h1_reframed_result = combine_h1(h1_revised_result, h1b_result)
h2c_result = validate_h2_comparative(scenario_results)

hypothesis_validation = {
    "H1_reframed": h1_reframed_result,
    "H1_revised": h1_revised_result,          # = H1a
    "H1b_image_economy": h1b_result,
    "H1": h1_result,
    "H2": h2c_result,                         # v1.5.0: comparative
    "H2_original_ratio": h2_result,           # cloud coverage > 50% of clear (kept for the record)
    "H3": h3_result,
    "H1_time_aware": h1_ta_result,
    "H1_onset_robustness": {_scenario_label(k): d.get("onset_robustness", {}) for k, d in scenario_results.items()},
    "approximation_guarantee": (
        "The k-step greedy attains at least 1 - 1/e of the optimum of a monotone submodular objective under a "
        "cardinality budget (Nemhauser, Wolsey & Fisher, 1978). This covers the coverage-first policy (the best "
        "utility per tile, summed) and, from v1.5.0, the time-aware policy, whose objective is the time-weighted "
        "coverage sum_t max w_e u_e,t -- also a weighted facility-location function. Cell 7C measures both against "
        "proven ILP optima of their own objectives. No guarantee is claimed for the Fail-Signal phases (re-planning "
        "on arrival and compensation beyond the budget), the earliest-completion cover's time (its greedy cover step "
        "has the classical logarithmic bound on the number of images), or the legacy v1.4.0 time-aware weighting."
        if globals().get("TIME_AWARE_OBJECTIVE", "facility") == "facility" else
        "(1-1/e) approximation guarantee (Nemhauser, Wolsey & Fisher, 1978) applies to the coverage-first policy. "
        "The legacy time-aware policy re-weights marginal gains per event; no guarantee is claimed for it."
    ),
}

# v1.5.3: this run's own statuses first. An "overall" status pools every
# scenario in the kernel's ledger (SUPPORTED only if every evaluable one is),
# so a run made after others in one kernel could read "H2 NOT SUPPORTED"
# overall while its own scenario was supported (15 October 2026).
_this_label = _scenario_label(_scenario_key()) if "_scenario_key" in globals() else None
_this_run = {"scenario": _this_label, "ledger_scenarios": [_scenario_label(k) for k in scenario_results],
             "statuses": {_hk: ((_hr.get("assessments") or {}).get(_this_label) or {}).get("status", "NOT RUN")
                          for _hk, _hr in (("H1a", h1_revised_result), ("H1b", h1b_result),
                                           ("H1 (H1a and H1b)", h1_reframed_result), ("H2", h2c_result),
                                           ("H3", h3_result))},
             "note": "statuses of this run's own scenario; the overall_status of each hypothesis pools "
                     "every scenario in ledger_scenarios"}
hypothesis_validation["this_run"] = _this_run
_ov = (f" (overall: pools the {len(scenario_results)} scenarios in this kernel)"
       if len(scenario_results) > 1 else "")
print("=" * 70)
print("HYPOTHESIS VALIDATION")
print("=" * 70)
print(f"\nThis run ({_this_label}): " + ", ".join(f"{k} {v}" for k, v in _this_run["statuses"].items()))
if len(scenario_results) > 1:
    print(f"NOTE: this kernel's ledger holds {len(scenario_results)} scenarios; each 'overall' status below "
          f"pools them all. Quote this run's line above or the per-scenario lines.")
print(f"\n{h1_reframed_result['hypothesis']}: {h1_reframed_result['overall_status']}{_ov}")
for n_label, a in h1_reframed_result['assessments'].items():
    print(f"  {n_label}: {a['status']} (H1a {a['H1a']}, H1b {a['H1b']})")
print(f"\n{h1_revised_result['hypothesis']}: {h1_revised_result['overall_status']}{_ov}")
for n_label, a in h1_revised_result['assessments'].items():
    _s = a.get(f"{a.get('primary_selector', 'time_aware')}_clear") or {}
    print(f"  {n_label}: {a['status']} -- {a.get('primary_selector')} within margin at "
          f"{_s.get('n_within_margin')}/{_s.get('n_onsets_evaluated')} onsets "
          f"(evaluable from {a.get('min_evaluable_onsets')}); time-aware: {a.get('status_time_aware')}")
print(f"\n{h2c_result['hypothesis']}: {h2c_result['overall_status']}{_ov}")
for n_label, a in h2c_result['assessments'].items():
    _h = a.get("time_aware") or {}
    _dx = _h.get("difference_vs_openloop_expected") or {}
    _f = lambda v: "n/a" if v is None else f"{v:.2f}"
    print(f"  {n_label}: {a['status']} -- adaptive minus open-loop (expected), same images: mean "
          f"{_f(_dx.get('mean'))} points, one-sided 95% lower bound {_f(_dx.get('lcb'))} over {_dx.get('n')} "
          f"realisations; coverage-first: {a.get('status_coverage_first')}"
          + (f"; not evaluable: {'; '.join(_h.get('not_evaluable_reasons') or [])}" if a['status'] == 'NOT EVALUABLE' else ""))
print(f"\n{h1b_result['hypothesis']}: {h1b_result['overall_status']}{_ov}")
for n_label, a in h1b_result['assessments'].items():
    _s = a.get("clear") or {}
    print(f"  {n_label}: {a['status']} -- fewer images than FCFS within the margin at "
          f"{_s.get('n_success')}/{_s.get('n_onsets_evaluated')} onsets; median images ECC "
          f"{_s.get('median_images_ecc')} vs FCFS {_s.get('median_images_fcfs')}")
for h_key, h_data in [("H1 (original, < 10 min)", h1_result), ("H1 time-aware", h1_ta_result),
                      ("H2 (original ratio)", h2_result), ("H3", h3_result)]:
    _sup = "SUPERSEDED (kept for the record): " if h_data.get("superseded") else ""
    print(f"\n{_sup}{h_data['hypothesis']}: {h_data['overall_status']}{_ov}")
    if h_data.get("superseded"):
        print(f"  ({h_data['superseded_note']})")
    for n_label, a in h_data['assessments'].items():
        print(f"  {n_label}: {a}")

print(f"\n{hypothesis_validation['approximation_guarantee']}")

try:
    with open(os.path.join(OUT_DIR if 'OUT_DIR' in globals() else os.getcwd(), 'arf_hypothesis_validation.json'), 'w') as f:
        json.dump(hypothesis_validation, f, indent=2, default=str)
    logger.info("Hypothesis validation saved to arf_hypothesis_validation.json")
except Exception as e:
    logger.warning(f"Could not save hypothesis validation: {e}")

# =============================================================================
# RESULTS SUMMARY (CSV) & FINAL REPORT (TXT)
# =============================================================================
# One row per scenario in the ledger, so a kernel that has run several
# (ROI, N) scenarios produces a cross-scenario comparison.
import textwrap
from datetime import datetime, timezone
_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()


def _r(v, nd=2):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else round(float(v), nd)


def _frontier_summary(rows, k_budget):
    """Per sky, for the reachable-tile target: the smallest image budget whose
    earliest-completion time equals the floor, and the time at the run's budget."""
    out = {}
    for sky in ("clear", "cloud"):
        f = [r for r in rows if r["sky"] == sky and r["target"] == "reachable"]
        if not f:
            continue
        fl = f[0]["floor_min"]
        at_floor = [r["k"] for r in f if r["t_complete_min"] is not None and fl is not None
                    and abs(r["t_complete_min"] - fl) < 1e-9]
        at_k = [r for r in f if r["k"] == k_budget]
        out[f"ECC Images to Match Floor - {sky}"] = min(at_floor) if at_floor else None
        out[f"ECC Time at Budget k (min) - {sky}"] = _r(at_k[0]["t_complete_min"], 1) if at_k else None
        out[f"Floor to Reachable Target (min) - {sky}"] = _r(fl, 1)
    return out


def _h1_v140_columns(key, ob):
    """Summary columns of the v1.4.0 H1 framing (H1a, H1b, quality at the level)."""
    lab = _scenario_label(key)
    out = {"H1 Status (v1.4.0: H1a and H1b)": h1_reframed_result["assessments"].get(lab, {}).get("status"),
           "H1a Status (latency non-inferiority)": h1_revised_result["assessments"].get(lab, {}).get("status"),
           "H1b Status (image economy)": h1b_result["assessments"].get(lab, {}).get("status")}
    a, b, q = ob.get("h1_revised"), ob.get("h1b_image_economy"), ob.get("quality_at_level")
    lvl = f"{(a or b or {}).get('reachable_coverage_pct', 95.0):g}%r"
    if a:
        out.update({"H1 Min Evaluable Onsets": a.get("min_evaluable_onsets"),
                    "H1a Primary Selector": a.get("primary_series", "time_aware"),
                    "H1a Within Margin 95% Lower Bound (%) - ECC clear":
                        _r(a["series"]["earliest_completion_clear"].get("pct_within_margin_lcb95"), 1),
                    "H1a Within Margin 95% Lower Bound (%) - time-aware clear":
                        _r(a["series"]["time_aware_clear"].get("pct_within_margin_lcb95"), 1)})
    if b:
        s = b["series"]["clear"]
        out.update({"H1b Onsets Evaluated - clear": s["n_onsets_evaluated"],
                    "H1b Onsets Needing One Image - clear": s["n_onsets_one_image"],
                    "H1b ECC Fewer Images Within Margin (%) - clear": _r(s["pct_success"], 1),
                    "H1b 95% Lower Bound (%) - clear": _r(s["pct_success_lcb95"], 1),
                    f"Median Images to {lvl} - ECC clear": _r(s["median_images_ecc"], 1),
                    f"Median Images to {lvl} - FCFS clear": _r(s["median_images_fcfs"], 1)})
    if q:
        s = q["series"]["clear"]
        out.update({f"Quality at {lvl} (median) - {lab_} clear": _r(s[k_]["median"], 3)
                    for k_, lab_ in (("time_aware", "time-aware"), ("fcfs_unbudgeted", "FCFS"),
                                     ("earliest_completion", "ECC"))})
        out[f"Quality at {lvl} time-aware vs FCFS (median relative difference %) - clear"] = _r(
            s["time_aware_vs_fcfs"]["median_relative_difference_pct"], 1)
    return out


_summary_rows = []
for _key, _d in scenario_results.items():
    _roi, _n = _key if isinstance(_key, tuple) and len(_key) == 2 else ("?", _key)
    _summary_rows.append({
        "ROI": _roi,
        "N": _n,
        "Candidate Events": _d.get("candidate_events"),
        "Selection Budget (k)": _d.get("selection_budget_k"),
        "Coverage - Clear Sky (%)": _r(_d.get("coverage_clear_pct")),
        "Coverage - Cloud/Fail-Signal (%)": _r(_d.get("coverage_cloud_pct")),
        "Cloud Knowledge": _d.get("cloud_knowledge", "oracle"),
        "Oracle Coverage - Cloud (%) (perfect cloud knowledge)": _r(_d.get("oracle_coverage_cloud_pct")),
        "Oracle H1 Latency - Cloud (min) (perfect cloud knowledge)": _r(_d.get("oracle_latency_cloud_min"), 1),
        "Cloud Re-plans": _d.get("cloud_replans"),
        "Cloud Compensatory Passes": _d.get("cloud_compensatory_passes"),
        "H3 Redundancy at 95% - Clear": _r(_d.get("h3_redundancy_clear"), 3),
        "H3 Redundancy at 95% - Cloud": _r(_d.get("h3_redundancy_cloud"), 3),
        "H3 Passes to 95% - Clear": _d.get("h3_passes_clear"),
        "Cloud Redundancy Floor (model)": _r(_d.get("cloud_redundancy_floor"), 3),
        "Redundancy all passes - Clear (budget-driven)": _r(_d.get("redundancy_clear"), 3),
        "Redundancy all passes - Cloud (budget-driven)": _r(_d.get("redundancy_cloud"), 3),
        "H1 Latency (completion) - Clear (min)": _r(_d.get("latency_clear_min"), 1),
        "H1 Latency (completion) - Cloud (min)": _r(_d.get("latency_cloud_min"), 1),
        "Last Selected Delivery - Clear (min)": _r(_d.get("last_delivery_clear_min"), 1),
        "Last Selected Delivery - Cloud (min)": _r(_d.get("last_delivery_cloud_min"), 1),
        "Reachable Tiles (%)": _r(_d.get("reachable_tiles_pct")),
        **(lambda _h: {} if not _h else {
            "H1r Margin (min)": _r(_h["floor_margin_min"], 1),
            **{f"H1r Onsets Within Margin (%) - {lab} {sky}": _r(_h["series"][f"{key}_{sky}"]["pct_within_margin"], 1)
               for key, lab in (("time_aware", "time-aware"), ("coverage_first", "coverage-first"), ("fcfs", "FCFS"),
                                ("earliest_completion", "ECC"))
               for sky in ("clear", "cloud")},
            "H1r Onsets Evaluated - clear": _h["series"]["time_aware_clear"]["n_onsets_evaluated"],
            "H1r Median Gap to Floor (min) - time-aware clear": _r(_h["series"]["time_aware_clear"]["median_gap_min"], 1),
        })((_d.get("onset_robustness") or {}).get("h1_revised")),
        **_h1_v140_columns(_key, _d.get("onset_robustness") or {}),
        **(lambda _c: {} if not _c.get("h2") else {
            "H2 Status (v1.5.0: adaptive vs open-loop, same images)": _c["h2"]["time_aware"]["status"],
            "H2 Status - coverage-first": _c["h2"]["coverage_first"]["status"],
            "H2 Cloud Realisations": _c["h2"]["time_aware"]["n_realisations"],
            "H2 Cloud Model": (_c.get("cloud_model") or {}).get("model"),
            "H2 Mean Images (adaptive, time-aware)": _r(_c["h2"]["time_aware"]["images_mean"], 1),
            **{f"H2 Mean Coverage (%) - {lab}": _r(_c["h2"]["time_aware"]["coverage_pct_mean"].get(t), 2)
               for t, lab in (("adaptive", "adaptive"), ("openloop_expected", "open-loop expected"),
                              ("openloop_clear", "open-loop clear-sky"), ("oracle", "oracle"))},
            "H2 Adaptive minus Open-loop (points, mean)": _r(_c["h2"]["time_aware"]["difference_vs_openloop_expected"]["mean"], 3),
            "H2 Adaptive minus Open-loop (points, one-sided 95% LB)": _r(_c["h2"]["time_aware"]["difference_vs_openloop_expected"]["lcb"], 3),
            "H2 Share of Oracle Gain Captured (median)": _r(_c["h2"]["time_aware"]["share_of_oracle_gain_captured_median"], 3),
        })(_d.get("cloud_robustness") or {}),
        "H2 Status (original ratio > 50%)": (lambda _st: _st)(
            "SUPPORTED" if (h2_result["assessments"].get(_scenario_label(_key), {}).get("coverage_maintenance_pct") or 0)
            > h2_result["maintenance_threshold_pct"] else "NOT SUPPORTED"),
        "Time-aware Objective": _d.get("time_aware_objective"),
        "Greedy/optimal time-weighted coverage - time-aware (median)": _r((_d.get("ilp_gap") or {}).get("time_aware_tw_ratio_median"), 4),
        "Greedy/optimal time-weighted coverage - time-aware (worst)": _r((_d.get("ilp_gap") or {}).get("time_aware_tw_ratio_min"), 4),
        "Decay Time Scale": _d.get("latency_decay_scale"),
        "Decay tau (min)": _r(_d.get("latency_decay_tau_min"), 1),
        "Reachable Usable Tiles - Clear (%)": _r(_d.get("reachable_usable_clear_pct")),
        "H3 Status": h3_result["assessments"].get(_scenario_label(_key), {}).get("status"),
        "Selected Events - Clear": _d.get("selected_events_clear"),
        "Selected Events - Cloud": _d.get("selected_events_cloud"),
        "Elbow k": _d.get("real_elbow_k"),
        "Saturation k": _d.get("saturation_k"),
        "Gap-Fill Enabled": _d.get("gap_fill_enabled"),
        "Decay Rate (time-aware)": _d.get("latency_decay_rate"),
        "H1 Latency time-aware - Clear (min)": _r(_d.get("ta_latency_clear_min"), 1),
        "H1 Latency time-aware - Cloud (min)": _r(_d.get("ta_latency_cloud_min"), 1),
        "Mean Quality - Clear (coverage-first)": _r(_d.get("mean_quality_clear"), 3),
        "Mean Quality - Clear (time-aware)": _r(_d.get("ta_mean_quality_clear"), 3),
        **{f"Onsets <H1 (%) - {lab}": _r(((_d.get("onset_robustness") or {}).get("series", {})
                                          .get(key, {}) or {}).get("pct_onsets_under_h1"), 1)
           for lab, key in (("best possible clear", "best_possible_clear"),
                            ("coverage-first clear", "coverage_first_clear"),
                            ("time-aware clear", "time_aware_clear"),
                            ("FCFS clear", "fcfs_clear"), ("random clear", "random_clear"))},
        **{f"Onset median (min) - {lab}": _r(((_d.get("onset_robustness") or {}).get("series", {})
                                            .get(key, {}) or {}).get("median_min"), 1)
           for lab, key in (("best possible clear", "best_possible_clear"),
                            ("coverage-first clear", "coverage_first_clear"),
                            ("time-aware clear", "time_aware_clear"),
                            ("FCFS clear", "fcfs_clear"), ("random clear", "random_clear"))},
        **{f"Onset median to {tg:g}% (min) - {lab} {sky}": _r(((_d.get("onset_robustness") or {})
                .get("equal_coverage_series", {}).get(f"{key}_{sky}_t{tg:g}", {}) or {}).get("median_min"), 1)
           for tg in COVERAGE_TARGETS_PCT for sky in ("clear", "cloud")
           for lab, key in (("best possible", "best_possible"), ("time-aware", "time_aware"), ("FCFS", "fcfs"),
                            ("ECC", "earliest_completion"))},
        **{f"Onsets reaching {tg:g}% (%) - {lab} {sky}": _r(((_d.get("onset_robustness") or {})
                .get("equal_coverage_series", {}).get(f"{key}_{sky}_t{tg:g}", {}) or {}).get("pct_onsets_reached"), 1)
           for tg in COVERAGE_TARGETS_PCT for sky in ("clear", "cloud")
           for lab, key in (("time-aware", "time_aware"), ("FCFS", "fcfs"), ("ECC", "earliest_completion"))},
        **{f"Greedy/optimal utility - {lab} (median)": _r((_d.get("ilp_gap") or {}).get(f"{key}_ratio_median"), 4)
           for lab, key in (("coverage-first", "coverage_first"), ("time-aware", "time_aware"), ("FCFS", "fcfs"))},
        **_frontier_summary(_d.get("cost_latency_frontier") or [], _d.get("selection_budget_k")),
        **(lambda _v: {"Validation Status": (_v.get("status") or "")[:60],
                       "Validation Precision": _r(_v.get("precision"), 3),
                       "Validation Recall": _r(_v.get("recall"), 3),
                       "Validation Predicted Tiles Confirmed (%)": _r(_v.get("predicted_tiles_confirmed_pct"), 1),
                       "Validation Real Cloud Cover Mean (%)": _r(_v.get("actual_cloud_cover_mean_pct"), 1)}
           )(_d.get("acquisition_validation") or {}),
    })
results_summary_df = pd.DataFrame(_summary_rows)
results_summary_df.to_csv(os.path.join(_out_dir, "arf_results_summary.csv"), index=False)

_labels = {_scenario_label(k): d for k, d in scenario_results.items()}
_bar, _sub = "=" * 70, "-" * 70
_lines = [_bar, "ALGORITHMIC RECONSTRUCTION FRAMEWORK (ARF)", "FINAL RESULTS REPORT", _bar, "",
          f"Date (UTC): {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}",
          f"Scenarios Evaluated: {list(_labels)}",
          f"Framework Version: {globals().get('__version__', 'unknown')}", "",
          _sub, "THIS RUN", _sub, "", f"Scenario: {_this_label}",
          *[f"  {k}: {v}" for k, v in _this_run["statuses"].items()],
          *([f"NOTE: the ledger holds {len(scenario_results)} scenarios (this kernel ran more than one); each "
             f"'overall' status under HYPOTHESIS VALIDATION pools them all. This folder's run is the one above."]
            if len(scenario_results) > 1 else []), "",
          _sub, "CROSS-SCENARIO COMPARISON", _sub, "",
          results_summary_df.to_string(index=False) if len(results_summary_df) else "(no scenarios in ledger)", "",
          _sub, "HYPOTHESIS VALIDATION", _sub, ""]
# v1.5.0: what the times mean in this run
_acq = {lab: d.get("acquisition_time_only") for lab, d in _labels.items()}
if _acq and all(v is True for v in _acq.values()):
    _lines += ["Times: delivery time = acquisition time + downlink wait. Every wait in this run is 0 (the modelled "
               "ground stations see each imaging pass), so all latencies below are ACQUISITION times; operator "
               "downlink, processing and distribution delays are not modelled.", ""]
elif _acq:
    _lines += [f"Times: delivery time = acquisition time + downlink wait to the modelled ground stations "
               f"(acquisition-only scenarios: {[k for k, v in _acq.items() if v]}). Operator processing and "
               f"distribution delays are not modelled.", ""]
# v1.5.1: times are measured from each scenario's own window start
_ref = {lab: d.get("reference_epoch") for lab, d in _labels.items()}
if len({v for v in _ref.values() if v}) > 1:
    _lines += ["WARNING: the scenarios in this ledger start at different times "
               + "; ".join(f"{k}: {v}" for k, v in _ref.items())
               + ". Their latencies are measured from different instants and are not comparable across "
                 "scenarios; set REFERENCE_EPOCH_UTC or ANALYSIS_EPOCH (Cell 0) to one value for runs to compare.", ""]
_lines += ["Onset statistics (H1a, H1b) describe one window: consecutive onsets share most of their passes, so "
           "their binomial bounds are optimistic; confirmatory tests use independent windows.", ""]

_pc = lambda v: "n/a" if v is None else f"{v:.1f}%"
_lines.append(f"{h1_reframed_result['hypothesis']}: {h1_reframed_result['overall_status']}{_ov}")
for _lab, _a in h1_reframed_result["assessments"].items():
    _lines.append(f"  {_lab}: {_a['status']} (H1a {_a['H1a']}, H1b {_a['H1b']})")
_lines.append("")
_lines.append(f"{h1_revised_result['hypothesis']}: {h1_revised_result['overall_status']}{_ov}")
for _lab, _a in h1_revised_result["assessments"].items():
    if "statement" not in _a:
        _lines.append(f"  {_lab}: {_a['status']}")
        continue
    _lines.append(f"  {_lab}: {_a['status']}")
    _lines.append(f"    {_a['statement']}  [margin {_a['floor_margin_source']}; evaluable from "
                  f"{_a.get('min_evaluable_onsets')} onsets]")
    _lines.append(f"    primary selector: {_a.get('primary_selector')}; time-aware selector: {_a.get('status_time_aware')}")
    for _n in ("earliest_completion_clear", "time_aware_clear", "coverage_first_clear", "fcfs_clear",
               "earliest_completion_cloud", "time_aware_cloud"):
        _s = _a[_n]
        _g = "n/a" if _s["median_gap_min"] is None else f"{_s['median_gap_min']:.0f} min"
        _lines.append(f"    {_n:22s} within margin at {_s['n_within_margin']}/{_s['n_onsets_evaluated']} onsets "
                      f"({_pc(_s['pct_within_margin'])}; 95% lower bound {_pc(_s.get('pct_within_margin_lcb95'))}); "
                      f"median gap to floor {_g}")
_lines.append(f"  Note: {h1_revised_result['note']}")
_lines.append("")
_lines.append(f"{h1b_result['hypothesis']}: {h1b_result['overall_status']}{_ov}")
for _lab, _a in h1b_result["assessments"].items():
    if "statement" not in _a:
        _lines.append(f"  {_lab}: {_a['status']}")
        continue
    _lines.append(f"  {_lab}: {_a['status']}")
    _lines.append(f"    {_a['statement']}  [evaluable from {_a.get('min_evaluable_onsets')} onsets]")
    for _sky in ("clear", "cloud"):
        _s = _a[_sky]
        _lines.append(f"    {_sky:5s} fewer images within the margin at {_s['n_success']}/{_s['n_onsets_evaluated']} "
                      f"onsets ({_pc(_s['pct_success'])}; 95% lower bound {_pc(_s['pct_success_lcb95'])}); median "
                      f"images ECC {_s['median_images_ecc']} vs FCFS {_s['median_images_fcfs']}; "
                      f"{_s['n_onsets_one_image']} onsets need one image (not evaluated)")
    _q = ((_labels.get(_lab, {}).get("onset_robustness") or {}).get("quality_at_level") or {}).get("series", {})
    for _sky in ("clear", "cloud"):
        _s = _q.get(_sky)
        if _s:
            _f3 = lambda v: "n/a" if v is None else f"{v:.3f}"
            _rd = _s["time_aware_vs_fcfs"]["median_relative_difference_pct"]
            _lines.append(f"    {_sky:5s} mosaic quality at the level (median; reported, not tested): time-aware "
                          f"{_f3(_s['time_aware']['median'])}, FCFS {_f3(_s['fcfs_unbudgeted']['median'])}, ECC "
                          f"{_f3(_s['earliest_completion']['median'])}; time-aware vs FCFS median "
                          f"{'n/a' if _rd is None else format(_rd, '+.1f') + '%'} "
                          f"(higher at {_pc(_s['time_aware_vs_fcfs']['pct_onsets_higher'])} of onsets)")
_lines.append(f"  Note: {h1b_result['note']}")
_lines.append("")
_lines.append("SUPERSEDED HYPOTHESES (kept for the record; not tested by the protocol, not results to quote):")
_lines.append(f"SUPERSEDED {h1_result['hypothesis']}: {h1_result['overall_status']}{_ov} "
              f"({h1_result['superseded_note']})")
_lines.append(f"  Latency = mosaic completion time (when delivered coverage reaches its final value);"
              f" target < {h1_result['threshold_min']:g} min.")
for _lab, _a in h1_result["assessments"].items():
    _lat = _a.get("latency_clear_min")
    _last = _labels.get(_lab, {}).get("last_delivery_clear_min")
    _lines.append(f"  {_lab}: " + ("N/A" if _lat is None else
                  f"{_lat:.1f} min ({'UNDER' if _a.get('under_threshold') else 'OVER'}-{h1_result['threshold_min']:g})")
                  + ("" if _last is None else f"  [last selected delivery: {_last:.1f} min]"))
_lines.append("")
_lines.append(f"SUPERSEDED {h1_ta_result['hypothesis']}: {h1_ta_result['overall_status']}{_ov} "
              f"({h1_ta_result['superseded_note']})")
for _lab, _a in h1_ta_result["assessments"].items():
    _lat = _a.get("latency_clear_min")
    _rate = _labels.get(_lab, {}).get("latency_decay_rate")
    _lines.append(f"  {_lab}: " + ("N/A" if _lat is None else
                  f"{_lat:.1f} min ({'UNDER' if _a.get('under_threshold') else 'OVER'}-{h1_ta_result['threshold_min']:g})")
                  + ("" if _rate is None else f"  [decay rate {_rate:g}; full sweep in arf_decay_sweep.csv]"))
_lines.append("")
_lines.append("H1 onset robustness (completion time across shifted onsets; see arf_onset_robustness.csv):")
for _lab, _d in _labels.items():
    _ob = _d.get("onset_robustness") or {}
    if not _ob:
        _lines.append(f"  {_lab}: not run")
        continue
    _lines.append(f"  {_lab}: {_ob['series']['best_possible_clear']['n_onsets']} onsets over "
                  f"{_ob['onset_shift_hours']} h, {_ob['horizon_hours']:g} h of passes each")
    _fm = lambda v: "n/a" if v is None else f"{v:.0f}"
    _lines.append("    Completion time (to each policy's OWN final coverage -- compare across policies only "
                  "where coverage is equal):")
    for _key, _s in _ob["series"].items():
        _lines.append(f"    {_key:22s} median {_fm(_s['median_min']):>5} | p90 {_fm(_s['p90_min']):>5} | "
                      f"max {_fm(_s['max_min']):>5} min | under the superseded 10-min H1 in {_s['pct_onsets_under_h1']:.1f}% of onsets")
    for _tg in _ob.get("coverage_targets_pct", []):
        _lines.append(f"    Equal coverage -- time to {_tg:g}% of ROI tiles (like-for-like across policies):")
        for _key, _s in _ob["equal_coverage_series"].items():
            if _key.endswith(f"_t{_tg:g}"):
                _lines.append(f"    {_key[:-len(f'_t{_tg:g}')]:22s} median {_fm(_s['median_min']):>5} | p90 {_fm(_s['p90_min']):>5} | "
                              f"reached in {_s['pct_onsets_reached']:5.1f}% | under the superseded 10-min H1 in {_s['pct_onsets_under_h1']:.1f}% of onsets")
for _lab, _d in _labels.items():
    _pcmp = _d.get("policy_comparison") or []    # v1.5.2: was _pc, which shadowed the % formatter above
    if _pcmp:
        _lines.append("")
        _lines.append(f"Policy comparison, {_lab} (full window, budget k, no gap-fill; see arf_policy_comparison.csv):")
        _fmt_t = lambda v: "never" if v is None else f"{v:.0f} min"
        for _p in _pcmp:
            _tt = " | ".join(f"to {tg:g}%: {_fmt_t(_p.get(f't{tg:g}_min')):>9}" for tg in COVERAGE_TARGETS_PCT)
            if _p["policy"] == "best_possible":
                _lines.append(f"    {_p['sky']:5s} {'best possible':16s} {'':>13} | {'':>16} | {'':>18} | {_tt}")
                continue
            _lat = "n/a" if _p["h1_latency_min"] is None else f"{_p['h1_latency_min']:.0f} min"
            _lines.append(f"    {_p['sky']:5s} {_p['policy']:16s} H1 {_lat:>9} | coverage {_p['coverage_pct']:6.2f}% "
                          f"| mean quality {_p['mean_quality']:.3f} | {_tt}")
    _fr = _d.get("cost_latency_frontier") or []
    if _fr:
        _lines.append(f"Cost-latency frontier, {_lab} (earliest-completion cover; see arf_cost_latency_frontier.csv):")
        for _sky in ("clear", "cloud"):
            _f = [r for r in _fr if r["sky"] == _sky and r["target"] == "reachable"]
            if not _f:
                continue
            _fl = _f[0]["floor_min"]
            _pts = ", ".join(f"k={r['k']}: {'never' if r['t_complete_min'] is None else format(r['t_complete_min'], '.0f') + ' min'}"
                             for r in _f)
            _lines.append(f"    {_sky:5s} to {_f[0]['target_pct']:g}% of {_f[0]['denominator_tiles']} reachable tiles "
                          f"(floor {'never' if _fl is None else format(_fl, '.0f') + ' min'}): {_pts}")
    _ig = _d.get("ilp_gap") or {}
    if _ig.get("solved"):
        _hz = _ig.get("horizons_used_min") or [int(_ig["horizon_hours"] * 60)]
        _lines.append(f"Optimality gap, {_lab}: {_ig['solved']}/{_ig['instances']} ILP instances solved to proven "
                      f"optimality ({_ig.get('time_limited', 0)} hit the time limit and are excluded); windows of "
                      f"{', '.join(f'{h} min' for h in _hz)} of passes, k={_ig['budget_k']}, clear sky")
        for _n in ("coverage_first", "time_aware", "fcfs"):
            _lines.append(f"    {_n:16s} utility / optimum: median {_ig[_n + '_ratio_median']:.4f}, "
                          f"worst {_ig[_n + '_ratio_min']:.4f}")
        if _ig.get("ecc_instances_exact"):
            _lines.append(f"    earliest-completion cover: greedy image count equalled the exact minimum in "
                          f"{_ig['ecc_matches_exact']}/{_ig['ecc_instances_exact']} windows "
                          f"(worst case {_ig['ecc_max_extra_images']} extra image(s))")
_lines.append("")
_lines.append(f"{h2c_result['hypothesis']}: {h2c_result['overall_status']}{_ov}")
for _lab, _a in h2c_result["assessments"].items():
    if _a.get("status") == "NOT RUN":
        _lines.append(f"  {_lab}: NOT RUN")
        continue
    _cm = _a.get("cloud_model") or {}
    _lines.append(f"  {_lab}: {_a['status']} (coverage-first: {_a.get('status_coverage_first')}); cloud model "
                  f"{_cm.get('model')}" + (f", L = {_cm.get('length_km'):g} km, T = {_cm.get('time_h'):g} h"
                                           if _cm.get('model') == 'correlated' else "")
                  + f"; seeds {_a.get('seeds')}")
    for _obj in ("time_aware", "coverage_first"):
        _h = _a.get(_obj) or {}
        if not _h:
            continue
        _cv = _h.get("coverage_pct_mean") or {}
        _dx, _dc = _h["difference_vs_openloop_expected"], _h["difference_vs_openloop_clear"]
        _f = lambda v, n=2: "n/a" if v is None else f"{v:.{n}f}"
        _lines.append(f"    {_obj}: mean images {_f(_h.get('images_mean'), 1)}; mean coverage adaptive "
                      f"{_f(_cv.get('adaptive'))}%, open-loop expected {_f(_cv.get('openloop_expected'))}%, open-loop "
                      f"clear-sky {_f(_cv.get('openloop_clear'))}%, oracle {_f(_cv.get('oracle'))}%")
        _lines.append(f"      adaptive minus open-loop (expected): mean {_f(_dx['mean'], 3)} points, one-sided 95% lower "
                      f"bound {_f(_dx['lcb'], 3)}, positive in {_f(_dx['pct_positive'], 0)}% of realisations; minus "
                      f"open-loop (clear-sky): mean {_f(_dc['mean'], 3)}; share of the oracle's gain captured "
                      f"(median) {_f(_h.get('share_of_oracle_gain_captured_median'), 3)}; {_h.get('replans_total')} "
                      f"re-plans, {_h.get('compensatory_total')} compensatory images in all")
        if _h.get("not_evaluable_reasons"):
            _lines.append(f"      not evaluable: {'; '.join(_h['not_evaluable_reasons'])}")
_lines.append(f"  Note: {h2c_result['note']}")
_lines.append("")
_lines.append(f"SUPERSEDED original {h2_result['hypothesis']} (kept for the record): "
              f"{h2_result['overall_status']}{_ov} (requires cloud coverage to maintain > "
              f"{h2_result['maintenance_threshold_pct']:g}% of clear-sky coverage; {h2_result['superseded_note']})")
for _lab, _a in h2_result["assessments"].items():
    _d = _labels.get(_lab, {})
    _lines.append(f"  {_lab}: Clear={_a['coverage_clear_pct']:.2f}%, Cloud={_a['coverage_cloud_pct']:.2f}% "
                  f"(maintained {_a['coverage_maintenance_pct']:.1f}%) [cloud knowledge: {_d.get('cloud_knowledge', 'oracle')}]")
    if _d.get("cloud_knowledge") == "realistic":
        _oc = _d.get("oracle_coverage_cloud_pct")
        _plan_txt = ("passes planned on expected coverage from each tile's cloud probability, each image's "
                     "cloud revealed on arrival, re-planned after every image"
                     if _d.get("cloud_planning") == "expected" else
                     "passes planned on clear-sky geometry, cloud revealed on arrival")
        _lines.append(f"    Realistic: {_plan_txt}; {_d.get('cloud_replans')} Fail-Signals, "
                      f"{_d.get('cloud_compensatory_passes')} compensatory passes. "
                      + ("" if _oc is None else f"Greedy given the cloud in advance: {_oc:.2f}%."))
_lines.append("")
_lines.append(f"{h3_result['hypothesis']}: {h3_result['overall_status']}{_ov} (requires coverage "
              f">{h3_result['coverage_threshold_pct']:g}% AND redundancy <{h3_result['redundancy_threshold']:g} "
              f"together, per the hypothesis as stated)")
for _lab, _a in h3_result["assessments"].items():
    _d = _labels.get(_lab, {})
    _ra = _a.get("redundancy_all_passes")
    if _a.get("status") == "NOT EVALUABLE":
        _lines.append(f"  {_lab}: NOT EVALUABLE -- the passes in the window can reach only "
                      f"{_a['reachable_clear_pct']:.2f}% of ROI tiles, below the {h3_result['coverage_threshold_pct']:g}% "
                      f"coverage at which redundancy is measured (selection covered {_a['coverage_clear_pct']:.2f}%)"
                      + ("" if _ra is None else f"  [redundancy over all selected passes: {_ra:.2f}, budget-driven]"))
        continue
    _red = _a['redundancy_clear']
    _lines.append(f"  {_lab}: {_a.get('status', '')} -- Coverage={_a['coverage_clear_pct']:.2f}%, Redundancy at "
                  f"{h3_result['coverage_threshold_pct']:g}% coverage="
                  + (f"{_red:.2f} " if np.isfinite(_red) else "never reached within the budget ")
                  + f"({_a.get('passes_to_coverage')} pass{'' if _a.get('passes_to_coverage') == 1 else 'es'})"
                  + ("" if _ra is None else f"  [over all selected passes: {_ra:.2f}, budget-driven]"))
    _rc, _fl = _d.get("h3_redundancy_cloud"), _d.get("cloud_redundancy_floor")
    if _rc is not None:
        _lines.append(f"    Cloud (reported, not part of the test): redundancy at {h3_result['coverage_threshold_pct']:g}% "
                      f"coverage={_rc:.2f} ({_d.get('h3_passes_cloud')} pass{'' if _d.get('h3_passes_cloud') == 1 else 'es'})"
                      + ("" if _fl is None else f"; model floor ~{_fl:.2f} at this cloud probability"))
_lines += ["", _sub, "VALIDATION AGAINST REAL ACQUISITIONS (Cell 7D)", _sub, ""]
for _lab, _d in _labels.items():
    _v = _d.get("acquisition_validation") or {}
    if not _v:
        continue
    _lines.append(f"  {_lab}: {_v.get('status')}")
    if _v.get("status") == "ok":
        _f = lambda x, nd=1: "n/a" if x is None else f"{x:.{nd}f}"
        _lines.append(f"    {_v['predicted_overpasses']} predicted and {_v['actual_overpasses']} real overpasses "
                      f"({', '.join(_v['platforms'])}) from {_v['window_start_utc']} to {_v['evaluated_until_utc']} UTC; "
                      f"{_v['matched']} matched within {_v['match_tolerance_min']:g} min.")
        _lines.append(f"    Precision {_f(_v['precision'], 3)}, recall {_f(_v['recall'], 3)}, median time offset "
                      f"{_f(_v['median_time_offset_min'], 2)} min; {_f(_v['predicted_tiles_confirmed_pct'])}% of predicted "
                      f"tiles lie inside the real scenes, which cover {_f(_v['actual_tiles_predicted_pct'])}% predicted.")
        _lines.append(f"    Real scene cloud cover: mean {_f(_v['actual_cloud_cover_mean_pct'])}%, median "
                      f"{_f(_v['actual_cloud_cover_median_pct'])}% (cloud model: {_v['model_cloud_probability_pct']:g}%). "
                      f"{_v['note']}")
        for _plat, _bp in (_v.get("by_platform") or {}).items():          # v1.5.2
            if not _bp.get("evaluated", True):
                _lines.append(f"    {_plat}: not evaluated yet (archive lag {_bp.get('archive_lag_h'):g} h)")
                continue
            _lines.append(f"    {_plat}: {_bp['matched']} of {_bp['predicted_overpasses']} predicted overpasses "
                          f"matched (precision {_f(_bp['precision'], 3)}), {_bp['actual_overpasses']} real "
                          f"(recall {_f(_bp['recall'], 3)}), evaluated to {_bp['evaluated_until_utc']} UTC")
_lines += ["", _sub, "THEORETICAL GUARANTEE", _sub, "",
           textwrap.fill(hypothesis_validation["approximation_guarantee"], width=70), ""]

with open(os.path.join(_out_dir, "arf_final_report.txt"), "w") as f:
    f.write("\n".join(_lines))
logger.info("Results summary and final report saved to arf_results_summary.csv / arf_final_report.txt")

print("\n=== CELL 8 COMPLETE ===")

