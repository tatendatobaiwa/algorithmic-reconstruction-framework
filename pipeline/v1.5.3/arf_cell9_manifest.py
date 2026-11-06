# =============================================================================
# CELL 9: REPRODUCIBILITY MANIFEST & RUN INSTRUCTIONS
# =============================================================================
#
# Records everything needed to reproduce this scenario's results:
#   - the execution environment (Python, platform, installed package versions)
#   - the scenario inputs (ROI, N, TLE source and epoch, and SHA-256 checksums
#     of the TLE cache and pinned satellite selection the run actually used)
#   - every output file the pipeline writes, with size and SHA-256 checksum
#   - the scenario_results ledger and the hypothesis validation (Cell 8)
#
# All outputs of Cells 0-9 are written to OUT_DIR (Cell 0): a per-scenario
# Google Drive folder on Colab (ARF_Outputs_<roi>_<fleet tag>[_from<date>]
# [_at<time>][_without-...], plus _v<release> when the folder already holds a
# complete run of another release), or the working directory when run locally.
# This cell validates files there.
#
# This cell does not execute other cells: notebook cells cannot re-run one
# another, so multi-scenario runs follow SCENARIO_RUN_INSTRUCTIONS below.
#
# Outputs (in OUT_DIR):
#   - arf_reproducibility_manifest.json (machine-readable)
#   - arf_pipeline_log.txt (human-readable)
# =============================================================================
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 9 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["9"] = _ARF_CELL_VERSION

import os
import sys
import json
import logging
import platform
import hashlib
from datetime import datetime, timezone
from importlib import metadata as _importlib_metadata
from typing import Dict, List, Any

logger = logging.getLogger("arf.cell9")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()
_framework_version = globals().get("__version__", "unknown")

# Every file Cells 0-8 write into OUT_DIR, grouped by the cell that writes it.
EXPECTED_OUTPUTS: Dict[str, List[str]] = {
    "Cell 0": ["arf_config.json"],
    "Cell 1": ["arf_cell1_validation.json"],
    "Cell 2": ["arf_cell2_summary.json"],
    "Cell 3": ["arf_cell3_metadata.json", "arf_coverage_map.png", "arf_coverage_map.pdf"],
    "Cell 4": ["arf_cell4_summary.json", "arf_access_events.csv"],
    "Cell 5": ["arf_cell5_summary.json", "arf_refined_events.csv", "arf_downlink_events.csv"],
    "Cell 6": ["arf_cell6_summary.json", "arf_selection_history_clear.csv",
               "arf_selection_history_cloud.csv", "arf_coverage_sweep.csv"],
    "Cell 7": ["arf_cell7_metadata.json",
               "arf_fig1_coverage_progression.png", "arf_fig1_coverage_progression.pdf",
               "arf_fig2_quality_progression.png", "arf_fig2_quality_progression.pdf",
               "arf_fig3_latency_analysis.png", "arf_fig3_latency_analysis.pdf",
               "arf_fig4_coverage_budget.png", "arf_fig4_coverage_budget.pdf",
               "arf_manuscript_figures.png", "arf_manuscript_figures.pdf",
               "arf_fig7_mosaic_map.png", "arf_fig7_mosaic_map.pdf", "arf_mosaic_tiles.csv"],
    "Cell 7B": ["arf_coverage_vs_time.csv"],                       # v1.5.2
    "Cell 7C": ["arf_decay_sweep.csv", "arf_policy_comparison.csv", "arf_cost_latency_frontier.csv",
                "arf_fig6_cost_latency.png", "arf_fig6_cost_latency.pdf"]
        + (["arf_ilp_gap.csv"] if globals().get("ILP_ENABLED", True) else []) + (
        ["arf_onset_robustness.csv", "arf_onset_robustness_summary.json",
         "arf_fig5_onset_robustness.png", "arf_fig5_onset_robustness.pdf",
         "arf_fig9_h1a_gap.png", "arf_fig9_h1a_gap.pdf"]                # fig9: v1.5.2
        if globals().get("ONSET_ROBUSTNESS_ENABLED", True) else [])
        # v1.5.0: cloud realisations (H2) and the cloud-scale sensitivity sweep
        + ["arf_cloud_robustness.csv", "arf_cloud_robustness_summary.json",
           "arf_fig8_cloud_robustness.png", "arf_fig8_cloud_robustness.pdf"]
        + (["arf_cloud_sensitivity.csv"] if (globals().get("CLOUD_SWEEP_LENGTHS_KM") or
                                             globals().get("CLOUD_SWEEP_TIMES_H")) else []),
    "Cell 7D": ["arf_acquisition_validation.csv", "arf_acquisition_validation_summary.json"],
    "Cell 8": ["arf_hypothesis_validation.json", "arf_results_summary.csv", "arf_final_report.txt"],
}

# Installed distributions to record. The orekit-data package is installed from
# git; its distribution name is tried under each common spelling.
_PACKAGE_DISTRIBUTIONS = {
    "numpy": ["numpy"], "pandas": ["pandas"], "matplotlib": ["matplotlib"],
    "h3": ["h3"], "shapely": ["shapely"], "orekit-jpype": ["orekit-jpype"],
    "jdk4py": ["jdk4py"], "orekit-data": ["orekitdata", "orekit-data", "orekit_data"],
    "scipy": ["scipy"], "rasterio": ["rasterio"],          # v1.5.0: ILP, cloud field; climatology reader (optional)
    "pyshp": ["pyshp"],                                    # v1.5.1: WDPA and Natural Earth shapefile reader (Cells 1, 3)
    "JPype1": ["JPype1", "jpype1"],                        # v1.5.2: the Python-Java bridge under orekit-jpype
}


def _java_version():
    """Version of the JVM running Orekit (v1.5.2), or None if no JVM is up."""
    try:
        import jpype
        if not jpype.isJVMStarted():
            return None
        from java.lang import System as _JSystem
        return {"java.version": str(_JSystem.getProperty("java.version")),
                "java.vendor": str(_JSystem.getProperty("java.vendor")),
                "java.vm.name": str(_JSystem.getProperty("java.vm.name"))}
    except Exception:
        return None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _file_record(path: str) -> Dict[str, Any]:
    if os.path.exists(path):
        return {"status": "FOUND", "size_bytes": os.path.getsize(path), "sha256": _sha256(path)}
    return {"status": "MISSING"}


def capture_environment() -> Dict[str, Any]:
    """Python, platform and installed package versions (via importlib.metadata,
    which reports the version actually installed, unlike module __version__)."""
    env_info = {
        "timestamp_utc": _utc_now(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
        "packages": {},
        "cwd": os.getcwd(),
        "out_dir": _out_dir,
    }
    for label, candidates in _PACKAGE_DISTRIBUTIONS.items():
        env_info["packages"][label] = "not installed"
        for dist in candidates:
            try:
                env_info["packages"][label] = _importlib_metadata.version(dist)
                break
            except _importlib_metadata.PackageNotFoundError:
                continue
    env_info["java"] = _java_version()                                       # v1.5.2
    env_info["validated_versions"] = globals().get("ARF_VALIDATED_VERSIONS")   # v1.5.2 (Cell 1)
    env_info["differs_from_validated"] = {
        k: {"validated": v, "installed": env_info["packages"].get(k)}
        for k, v in (globals().get("ARF_VALIDATED_VERSIONS") or {}).items()
        if env_info["packages"].get(k) not in (v, None) and k in env_info["packages"]}
    return env_info


def record_inputs() -> Dict[str, Any]:
    """The scenario definition and the exact input files this run used."""
    roi_def = ROI_DEFINITIONS.get(ACTIVE_ROI, {}) if 'ROI_DEFINITIONS' in globals() else {}
    inputs = {
        "active_roi": globals().get("ACTIVE_ROI"),
        "active_n": globals().get("ACTIVE_N"),
        "fleet": globals().get("FLEET", "starlink"),
        "fleet_tag": globals().get("FLEET_TAG"),
        "mission_table": globals().get("MISSION_TABLE") if globals().get("FLEET", "starlink") != "starlink" else None,
        "tle_source": roi_def.get("tle_source", "celestrak_live"),
        "historical_epoch": roi_def.get("historical_epoch"),
        "sim_epoch": str(sim_epoch) if 'sim_epoch' in globals() else None,
        # v1.5.1: where the window starts and why (runs are comparable only with
        # the same reference epoch)
        "reference_epoch": globals().get("_reference_epoch"),
        "reference_epoch_source": globals().get("REFERENCE_EPOCH_SOURCE"),
        "analysis_epoch": globals().get("ANALYSIS_EPOCH"),
        "n_satellites": len(constellation) if 'constellation' in globals() else None,
        "config": config.to_dict() if 'config' in globals() and hasattr(config, "to_dict") else None,
        "files": {},
    }
    for name in ("TLE_RAW_CACHE_PATH", "SELECTION_CACHE_PATH"):
        path = globals().get(name)
        if path:
            inputs["files"][name] = {"path": path, **_file_record(path)}
    # v1.5.1: every orbital-element file the run read (both CelesTrak groups of a
    # combined fleet, the reference-epoch groups, a Space-Track history file)
    for i, path in enumerate(globals().get("GP_INPUT_PATHS") or []):
        inputs["files"][f"GP_INPUT_{i}"] = {"path": path, **_file_record(path)}
    # WDPA boundary files (licence: not redistributed; recorded by checksum only)
    if roi_def.get("wdpa_sites") and "_find_wdpa_file" in globals():
        for sid in roi_def["wdpa_sites"]:
            p = _find_wdpa_file(sid, ARF_BOUNDARY_DIR)
            inputs["files"][f"WDPA_{sid}"] = ({"path": p, **_file_record(p)} if p
                                              else {"path": None, "status": "MISSING"})
    return inputs


def validate_all_outputs(out_dir: str = None) -> Dict[str, Any]:
    """Check every expected output (EXPECTED_OUTPUTS) exists in out_dir, with
    size and SHA-256 checksum."""
    out_dir = out_dir or _out_dir
    results = {"timestamp": _utc_now(), "out_dir": out_dir, "files": {},
               "n_expected": 0, "n_found": 0, "n_missing": 0}
    for cell, names in EXPECTED_OUTPUTS.items():
        for name in names:
            rec = {"cell": cell, **_file_record(os.path.join(out_dir, name))}
            results["files"][name] = rec
            results["n_expected"] += 1
            results["n_found" if rec["status"] == "FOUND" else "n_missing"] += 1
    results["all_present"] = results["n_missing"] == 0
    return results


# Every cell carries a version guard and reports in (v1.5.1; v1.5.0 guarded only
# the cells it changed). A cell that did not report was not run in this session,
# or is a copy from a release without the guard.
ARF_GUARDED_CELLS = ("0", "1", "2", "3", "4", "5", "6", "6B", "7", "7B", "7C", "7D", "8", "9")


def check_cell_versions() -> Dict[str, Any]:
    """Every guarded cell must have run in this session from this release.
    Raises unless ARF_ALLOW_MIXED_VERSIONS = True (the manifest then records
    the mismatch)."""
    run = dict(globals().get("_ARF_CELLS_RUN") or {})
    missing = [c for c in ARF_GUARDED_CELLS if c not in run]
    other = {c: v for c, v in run.items() if v != _framework_version}
    res = {"release": _framework_version, "cells_run": run, "missing": missing, "other_release": other,
           "consistent": not missing and not other}
    if not res["consistent"]:
        msg = (f"Mixed or incomplete notebook: cells not run from ARF {_framework_version}: {missing or 'none'}; "
               f"cells from another release: {other or 'none'}. Re-run Cells 0-9 from one release.")
        if not globals().get("ARF_ALLOW_MIXED_VERSIONS", False):
            raise RuntimeError(msg + " (Set ARF_ALLOW_MIXED_VERSIONS = True to write the manifest anyway; it "
                                     "will record the mismatch.)")
        logger.warning(msg)
    return res


def generate_reproducibility_manifest() -> Dict[str, Any]:
    """Build and save the manifest (JSON) and the pipeline log (TXT) in OUT_DIR."""
    logger.info("Generating reproducibility manifest...")
    cell_versions = check_cell_versions()

    def _json_scenario_key(k):
        return f"{k[0]}__N{k[1]}" if isinstance(k, tuple) and len(k) == 2 else str(k)

    manifest = {
        "manifest_version": "1.1.0",
        "framework": "Algorithmic Reconstruction Framework (ARF)",
        "framework_version": _framework_version,
        "cell_versions": cell_versions,                  # v1.5.0 version guard
        "config_overrides": dict(globals().get("ARF_CONFIG_OVERRIDES_APPLIED") or {}),   # v1.5.0 batch runs
        "manifest_timestamp": _utc_now(),
        "environment": capture_environment(),
        "inputs": record_inputs(),
        "output_validation": validate_all_outputs(),
        "scenario_results": ({_json_scenario_key(k): v for k, v in scenario_results.items()}
                             if 'scenario_results' in globals() else {}),
        "hypothesis_validation": hypothesis_validation if 'hypothesis_validation' in globals() else {},
    }

    manifest_path = os.path.join(_out_dir, "arf_reproducibility_manifest.json")
    try:
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2, default=str)
        logger.info(f"Manifest saved to {manifest_path}")
    except Exception as e:
        logger.warning(f"Could not save manifest: {e}")

    log_path = os.path.join(_out_dir, "arf_pipeline_log.txt")
    rule = "-" * 70 + "\n"
    try:
        with open(log_path, "w") as f:
            f.write("=" * 70 + "\nALGORITHMIC RECONSTRUCTION FRAMEWORK (ARF)\n"
                    "REPRODUCIBILITY MANIFEST & PIPELINE LOG\n" + "=" * 70 + "\n\n")
            f.write(f"Generated: {manifest['manifest_timestamp']}\n"
                    f"Framework Version: {_framework_version}\n\n")

            f.write(rule + "EXECUTION ENVIRONMENT\n" + rule + "\n")
            env = manifest["environment"]
            f.write(f"Python: {env['python_version']}\n"
                    f"Platform: {env['platform']['system']} {env['platform']['release']}\n"
                    f"Machine: {env['platform']['machine']}\n\nPackage Versions:\n")
            for pkg, version in env["packages"].items():
                f.write(f"  {pkg}: {version}\n")

            f.write("\n" + rule + "SCENARIO INPUTS\n" + rule + "\n")
            inp = manifest["inputs"]
            f.write(f"ROI: {inp['active_roi']}   fleet: {inp['fleet']}   "
                    f"satellites propagated: {inp['n_satellites']}\n"
                    f"Window start (reference epoch): {inp.get('reference_epoch')} UTC "
                    f"[{inp.get('reference_epoch_source')}]\n"
                    f"TLE source: {inp['tle_source']}   historical epoch: {inp['historical_epoch']}\n"
                    f"Simulation epoch: {inp['sim_epoch']}\n")
            for name, rec in inp["files"].items():
                f.write(f"  [{rec['status']}] {name}: {rec['path']}"
                        + (f" ({rec['size_bytes']} bytes, sha256 {rec['sha256'][:16]}...)" if rec["status"] == "FOUND" else "")
                        + "\n")

            f.write("\n" + rule + "OUTPUT FILE VALIDATION\n" + rule + "\n")
            val = manifest["output_validation"]
            f.write(f"Directory: {val['out_dir']}\nExpected: {val['n_expected']}   Found: {val['n_found']}   "
                    f"Missing: {val['n_missing']}   All present: {val['all_present']}\n\n")
            for name, rec in val["files"].items():
                size = rec.get("size_bytes", "N/A")
                f.write(f"  [{rec['status']}] {rec['cell']}: {name} ({size} bytes)\n")

            if manifest["scenario_results"]:
                f.write("\n" + rule + "SCENARIO RESULTS\n" + rule + "\n")
                for label, data in manifest["scenario_results"].items():
                    f.write(f"{label}:\n")
                    for key, value in data.items():
                        if key not in ("sweep_results", "real_k_sweep"):   # detailed sweeps: see CSV/JSON
                            f.write(f"  {key}: {value}\n")
                    f.write("\n")

            if manifest["hypothesis_validation"]:
                f.write(rule + "HYPOTHESIS VALIDATION\n" + rule + "\n")
                for h_key, h_data in manifest["hypothesis_validation"].items():
                    if isinstance(h_data, dict) and "overall_status" in h_data:
                        _sup = "SUPERSEDED (kept for the record): " if h_data.get("superseded") else ""
                        f.write(f"{_sup}{h_data.get('hypothesis', h_key)}: {h_data['overall_status']}\n")
                        for label, a in h_data.get("assessments", {}).items():
                            f.write(f"  {label}: {a}\n")
                        f.write("\n")
        logger.info(f"Pipeline log saved to {log_path}")
    except Exception as e:
        logger.warning(f"Could not save pipeline log: {e}")

    return manifest


# =============================================================================
# RUN INSTRUCTIONS
# =============================================================================

SCENARIO_RUN_INSTRUCTIONS = """
=============================================================================
ARF -- HOW TO RUN SCENARIOS (v1.5.2)
=============================================================================
A scenario is a region (active_roi) and a fleet (fleet, planet_swath_case),
live or for a past analysis window (analysis_epoch). Default: CKGR + Okavango
Delta with the open-data + PlanetScope fleet (conservative Dove swath).

STEP 1: Configure (Cell 0)
   - Set active_roi (e.g. "ckgr_okavango"), fleet ("earth-resources",
     "planet", "earth-resources+planet") and, with PlanetScope,
     planet_swath_case ("conservative" / "optimistic").
   - For a past window set analysis_epoch (ISO UTC, e.g. "2026-08-17T00:00:00");
     its orbital elements come from Space-Track gp_history and need the Colab
     secrets SPACETRACK_USERNAME / SPACETRACK_PASSWORD.
   - Runs you will COMPARE (one region, several fleets) must share their
     window start: the same analysis_epoch, or the same reference_epoch_utc for
     live runs (by default live runs start at the newest element set over both
     CelesTrak groups, which runs fetched in one session share).
   - Fleet sensitivity: exclude_satellites=("SENTINEL-2A",) leaves a tracked
     satellite out (Sentinel-2A flies a reduced extension-campaign plan, so
     its predicted passes are an upper bound); run with and without it.
   - Run Cell 0. Outputs go to the Drive folder
     ARF_Outputs_<roi>_<fleet tag>[_from<YYYYMMDD>][_at<time>][_without-...],
     one per scenario. A folder holding a complete run of another release is
     kept; this release then writes to <folder>_v<release digits> (v1.5.2).

STEP 2: Run Cells 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 6B -> 7 -> 7B -> 7C -> 7D -> 8 -> 9
   in order, all from ONE release: every cell checks Cell 0's version and
   stops on a mismatch, and Cell 9 refuses a manifest from a mixed or
   incomplete notebook.

STEP 3: Check the results
   - Cell 6:  selections (coverage-first, time-aware on the time-weighted
              coverage objective), cloud model, H1 latency, H3 redundancy
   - Cell 6B: data-driven selection budget (elbow and saturation k)
   - Cell 7:  figures and the mosaic map
   - Cell 7B: time to 25/50/75/90% coverage
   - Cell 7C: decay sweep, baselines, onset robustness (H1a/H1b on the
              earliest-completion cover), exact ILP gaps (both objectives),
              cloud realisations (H2: adaptive vs open-loop with the same
              images) and the cloud-scale sensitivity sweep (Fig. 8)
   - Cell 7D: predicted Sentinel-2 / Landsat passes vs the real archive
   - Cell 8:  hypothesis validation, arf_results_summary.csv, final report
   - Cell 9:  this manifest -- every expected file FOUND, cell versions consistent

STEP 4: Further scenarios and batches
   - Return to Cell 0, change the settings, and repeat. Without restarting
     the kernel, the scenario_results ledger keeps every scenario so far.
   - Many scenarios or windows: arf_batch_runner.py runs each configuration
     through all cells in a fresh namespace, resumably; the confirmatory
     protocol's runs are listed there, and arf_confirmatory_analysis.py
     applies the protocol's decision rules to the finished runs.
   - arf_colab_parity.py re-runs Cells 6-8 of saved runs and compares them
     with the reference results computed outside Colab (environment check).
=============================================================================
"""


# =============================================================================
# EXECUTE
# =============================================================================

manifest = generate_reproducibility_manifest()
validation = manifest["output_validation"]

print("\n" + "=" * 70)
print("  CELL 9: Reproducibility Manifest")
print("=" * 70)
# v1.5.2: real fleets are identified by their fleet tag (N is the Starlink proxy's)
print(f"  Scenario:                  ROI={manifest['inputs']['active_roi']}, "
      + (f"N={manifest['inputs']['active_n']}" if manifest['inputs']['fleet'] == "starlink"
         else f"fleet {manifest['inputs']['fleet_tag']} ({manifest['inputs']['n_satellites']} satellites)"))
print(f"  Output directory:          {validation['out_dir']}")
print(f"  Output files:              {validation['n_found']} / {validation['n_expected']} found"
      f"{'' if validation['all_present'] else '  (see missing below)'}")
for name, rec in validation["files"].items():
    if rec["status"] == "MISSING":
        print(f"    - MISSING ({rec['cell']}): {name}")
for name, rec in manifest["inputs"]["files"].items():
    print(f"  Input {name:18s} {rec['status']}")
print(f"  Manifest:                  arf_reproducibility_manifest.json")
print(f"  Pipeline log:              arf_pipeline_log.txt")
print("=" * 70)
print(SCENARIO_RUN_INSTRUCTIONS)
print(f"\n=== CELL 9 COMPLETE ===")
