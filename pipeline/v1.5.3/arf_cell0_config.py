# =============================================================================
# CELL 0: FRAMEWORK METADATA, CONFIGURATION, & SCENARIO SELECTION
# Algorithmic Reconstruction Framework (ARF)
# =============================================================================
#
# Description:
#   Single source of truth for every experimental parameter. One run is one
#   scenario: a region of interest (active_roi), a real Earth-observation
#   fleet (fleet, planet_swath_case) and an analysis window (window_hours,
#   live from the newest orbital elements, from a common reference_epoch_utc,
#   or retrospective from analysis_epoch). The Starlink proxy (fleet="starlink",
#   N in VALID_N_VALUES) and the Botswana bounding box / Nepal corridor ROIs are
#   kept only to reproduce the earliest runs.
#
# Theoretical basis:
#   - Cardinality-constrained weighted maximum coverage (NP-hard). The k-step
#     greedy attains at least 1 - 1/e (~0.632) of the optimum of a monotone
#     submodular objective (Nemhauser, Wolsey & Fisher 1978, Math. Program. 14,
#     265-294, doi:10.1007/BF01588971). Both selection objectives of Cell 6 are
#     of this kind; Cell 7C measures their exact gaps with an ILP.
#
# Execution workflow:
#   1. Set active_roi, fleet, planet_swath_case (and the window) in the
#      ARFConfig(...) call below, or pass ARF_CONFIG_OVERRIDES from a batch
#      driver (arf_batch_runner.py), and run this cell.
#   2. Run Cells 1 -> 9 in order, all from this release (each cell checks).
#   3. Results accumulate in the scenario_results ledger until the kernel
#      restarts; each scenario also writes its own output folder. From v1.5.3
#      each run's report leads with that run's own statuses; an "overall"
#      status pools every scenario in the ledger and says so.
#
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 ARF Contributors
# Author: Tatenda Allen Tobaiwa (BIUST)
# =============================================================================

"""
ARF Scenario Configuration Module
==================================

Defines and validates every parameter; later cells only read them.

Regions of interest (official WDPA boundaries, Cell 1):
    - ckgr_okavango: Central Kalahari Game Reserve + Okavango Delta (Botswana)
    - lower_zambezi: Lower Zambezi National Park (Zambia)
    - karoo_nature_reserve: Karoo Nature Reserve (South Africa)
Fleets (CelesTrak GP data, Cell 2; swaths and daylight limits in MISSION_TABLE):
    - earth-resources: Sentinel-2A/B/C and Landsat 8/9 (free and open data)
    - planet: PlanetScope Doves (commercial)
    - earth-resources+planet: both
Legacy: fleet="starlink" (proxy, N in VALID_N_VALUES), botswana_full,
nepal_trishuli_corridor.
"""

import os
import sys
import json
import logging
from dataclasses import dataclass, asdict
from typing import Dict, Any, Set, ClassVar, Tuple, Optional

# =============================================================================
# FRAMEWORK METADATA
# =============================================================================

__version__ = "1.5.3"   # 1.5.3: presentation only, every computed value as in 1.5.2 -- Figure 3 and composite panel (c) zoomed to their data with the legend under the axes, Figure 7 panel (b) title shortened so the map fits two IEEE columns, its cloud settings moved to the footer (Cell 7); Figures 5 and 9 at the H1a level of reachable tiles, Figure 5 in hours and Figure 9 in minutes on linear axes, legends under the panels (Cell 7C); superseded hypotheses labelled as such in the console, the report, the JSON and the pipeline log, and each run's own statuses printed first with ledger-wide statuses labelled as pooled (Cells 7, 7C, 8, 9); 1.5.2: per-cell audit -- sourced acquisition exclusions (Sentinel-2A tandem pause) applied per pass time, and an exclude_satellites option for fleet sensitivity runs (Cells 0, 2, 5); default scenario CKGR + Okavango with the open-data + PlanetScope fleet, and a complete run of another release in the output folder kept rather than overwritten (Cell 0); validated package versions pinned and checked (Cell 1), recorded with JPype and Java (Cell 9); deterministic sample tile (Cell 1); scenario ledger keyed by fleet tag and reset per run (Cells 0, 6); end-cap count after the daylight filter, per-satellite swath-edge print (Cell 5); figures at IEEE column widths, latency figures in hours against the physical floor and the H1a margin instead of the falsified 10-min line, new H1a-gap figure (Cells 7, 7C); coverage-vs-time CSV and undelivered-image-safe ordering (Cell 7B); per-platform archive lags and precision/recall (Cell 7D); H3 overall status by the same rule as H1a, H1b and H2 (Cell 8); documentation (Cells 0-9); 1.5.1: each tile observation assigned to the strip of its closest approach, no end-cap duplicates (Cell 5); fleet-independent reference epoch and optional common window start for live runs (Cells 0, 2); version guard in every cell; exact ILP tolerances (Cell 7C); time-aware weights on delivery time (Cell 6); single-window H1 statistics labelled descriptive; documentation corrected (Cells 0-9); 1.5.0: time-aware selection on a time-weighted coverage objective with the 1-1/e guarantee and its exact ILP gap (Cells 6, 7C); spatially and temporally correlated cloud with optional per-tile climatology (Cell 6); comparative H2 over cloud seeds (Fail-Signal vs open-loop plans with the same images) and H1a/H1b on the earliest-completion cover (Cells 7C, 8); per-cell version check (Cell 9); protected-area outline rasterised in saved maps (Cells 3, 7); config overrides for batch runs (Cell 0); 1.4.0: fleet-aware time-decay scale (Cell 6); H1 as latency non-inferiority (H1a) + image economy (H1b) with evaluability rules, H3 not-evaluable label (Cells 7C, 8); map furniture kept off the ROI (Cells 3, 7); Cell 7 metadata fleet size; NaN-safe delivery order; 1.3.1: journal-standard study-area map (Cell 3) and mosaic map (Cell 7), figures only; 1.3.0: Planet FOV swath, realistic cloud, ECC + frontier, retrospective runs, archive validation; 1.2.0: revised H1; earlier runs report 1.0.0
__framework_name__ = "Algorithmic Reconstruction Framework (ARF)"
__license__ = "MIT"
__authors__ = ["Tatenda Allen Tobaiwa"]
# Cells run in this session, with their release (version guard since v1.5.0, in
# every cell since v1.5.1; Cell 9 refuses a manifest from a mixed notebook).
# Reset whenever Cell 0 runs.
_ARF_CELLS_RUN = {"0": __version__}
__institution__ = "Botswana International University of Science and Technology (BIUST)"

# =============================================================================
# CUSTOM EXCEPTIONS
# =============================================================================

class ConfigurationError(Exception):
    """Raised when an invalid scenario or parameter configuration is detected."""
    pass


def _utc_naive(iso):
    """ISO 8601 time -> naive datetime in UTC. 'Z' and offsets (+02:00) are
    converted to UTC (v1.5.1; an offset previously crashed the future-time check
    or was dropped when formatting)."""
    from datetime import datetime as _dt, timezone as _tz
    d = _dt.fromisoformat(str(iso).strip().replace("Z", "+00:00"))
    return d.astimezone(_tz.utc).replace(tzinfo=None) if d.tzinfo is not None else d


def _utc_now_naive():
    from datetime import datetime as _dt, timezone as _tz
    return _dt.now(_tz.utc).replace(tzinfo=None)

# =============================================================================
# CONFIGURATION DATACLASS
# =============================================================================

@dataclass
class ARFConfig:
    """
    Configuration container for the Algorithmic Reconstruction Framework.
    Enforces runtime type validation and physical parameter sanity checks.
    """

    # --- Allowed Scenarios ---
    VALID_N_VALUES: ClassVar[Set[int]] = {5, 20, 100, 250, 1000}
    VALID_FLEETS: ClassVar[Set[str]] = {"starlink", "earth-resources", "planet", "earth-resources+planet"}
    VALID_ROI_NAMES: ClassVar[Set[str]] = {"botswana_full", "nepal_trishuli_corridor", "ckgr_okavango",
                                           "lower_zambezi", "karoo_nature_reserve"}

    # --- Scenario Parameters ---
    active_n: int = 1000                 # Starlink proxy only (legacy); real fleets use every operational satellite
    active_roi: str = "ckgr_okavango"    # v1.5.2 default (was the legacy Nepal corridor)
    cloud_probability: float = 0.35
    max_selection_steps: int = 30
    random_seed: int = 42

    # Imaging fleet. "starlink" = the Starlink proxy used for all results before
    # Path 1 (uniform swath_half_angle_deg, no daylight limit, N sampled). The
    # others are real Earth-observation fleets from CelesTrak GP data: every
    # operational satellite with a MISSION_TABLE entry, each with its own swath
    # and minimum sun elevation. Satellites without an entry are excluded and
    # listed, never given a guessed swath. v1.5.2 default: the open-data +
    # PlanetScope fleet (was the Starlink proxy).
    fleet: str = "earth-resources+planet"
    # PlanetScope frame orientation is not stated in the sources found, and the
    # generation of each Dove cannot be told from its name, so its swath is run
    # as a bound from Planet's field of view: "conservative" = PS2.SD Dove-R
    # (3.0 deg, ~26-27 km at 500 km) or "optimistic" = PSB.SD SuperDove
    # (4.0 deg, ~35-36 km). See MISSION_TABLE.
    planet_swath_case: str = "conservative"
    # Used only for missions whose sources give no sun-elevation limit
    # (Sentinel-2). An ASSUMPTION, borrowed from the Landsat criterion (5 deg).
    assumed_min_sun_elev_deg: float = 5.0

    # --- Physical & Astrodynamics Parameters ---
    # Simulation window. 48 h is the rapid-response case. For the open-data
    # fleet (fleet="earth-resources") use 384 h (16 days): the longest single-
    # satellite repeat cycle in that fleet, so every satellite flies its whole
    # ground track once. Sources (accessed 2026-10-09): Landsat 9 "Has a 16-day
    # repeat cycle" (https://www.usgs.gov/landsat-missions/landsat-9);
    # Sentinel-2 "The revisit frequency of each single Sentinel-2 satellite is
    # 10 days" (https://sentiwiki.copernicus.eu/web/s2-mission). A window other
    # than 48 h is added to the output folder name (e.g. _384h).
    window_hours: int = 48
    swath_half_angle_deg: float = 9.22
    h3_resolution: int = 7
    downlink_min_elevation_deg: float = 15.0

    # Fallback wait-time bound only. Overwritten dynamically in Cell 2.
    t_max_wait_minutes: float = 95.0

    # Time-aware objective (Cells 6 and 7C): marginal gain weighted by
    # exp(-rate * t_min / tau). Rate 0 == coverage-first.
    # The default 3.0 is where the Nepal decay sweep reached the physical-best
    # latency; it was chosen after seeing that sweep, so the full sweep
    # (decay_sweep, reported in Cell 7C) is the result to cite, not 3.0 alone.
    latency_decay_rate: float = 3.0
    decay_sweep: Tuple[float, ...] = (0.5, 1.0, 2.0, 3.0, 5.0)
    # tau of that weight (v1.4.0). "fleet": the longer of one orbital period
    # and the physical floor of the candidate passes (time for ALL of them to
    # cover h1_reachable_coverage_pct of the reachable tiles, clear-sky
    # geometry), so a sparse fleet that needs a day to cover the ROI is not
    # restricted to its first orbit. "fixed": one orbital period, as in v1.3.1
    # (T_MAX_WAIT_MINUTES). Identical whenever the floor is within one orbit.
    latency_decay_scale: str = "fleet"

    # Onset-robustness analysis (Cell 7C): request time shifted across
    # [0, onset_shift_hours] in onset_shift_step_min steps, each using the
    # following (window_hours - onset_shift_hours) hours of passes.
    onset_robustness_enabled: bool = True
    onset_shift_hours: int = 24
    onset_shift_step_min: int = 20

    # Baselines (Cell 7C): random draws per comparison, and the exact ILP
    # optimum on small sub-problems (clear sky) to measure the greedy's gap.
    random_baseline_draws: int = 10
    ilp_enabled: bool = True
    ilp_instances: int = 5
    ilp_horizon_hours: int = 6
    ilp_budget_k: int = 10
    # Instances are shrunk (shorter window) until they fit ilp_max_pairs; at
    # ~100k pass-tile pairs HiGHS proves optimality well inside the time limit.
    ilp_time_limit_s: float = 120.0
    ilp_max_pairs: int = 100_000

    # Equal-coverage latency (Cell 7C): time at which each policy's delivered
    # mosaic first reaches these coverage levels (% of ROI tiles). Completion
    # time alone is only comparable between policies reaching equal coverage.
    coverage_targets_pct: Tuple[float, ...] = (95.0, 99.0)

    # Revised H1 (pre-specified 2026-10-09, before the fresh-data test; the
    # original "< 10 min" H1 is still evaluated and reported as falsified):
    #   At >= h1_min_onset_pct of onsets, the time-aware selector's delivered
    #   mosaic reaches h1_reachable_coverage_pct of the tiles REACHABLE in that
    #   onset's horizon no later than h1_floor_margin_min after the physical
    #   floor (all passes) reaches the same level, using at most
    #   max_selection_steps passes. Clear sky is the test; cloud is reported.
    # h1_floor_margin_min=None derives the margin as ONE ORBITAL PERIOD at the
    # fleet's median altitude (Kepler's third law): missing one revisit costs
    # about one orbit. h1_min_onset_pct=90 is an ASSUMPTION set by the authors,
    # not derived; change it only before seeing fresh-data results.
    h1_reachable_coverage_pct: float = 95.0
    h1_floor_margin_min: Optional[float] = None
    h1_min_onset_pct: float = 90.0
    # H1 framing from v1.4.0 (Cells 7C and 8). The criterion above is H1a,
    # LATENCY NON-INFERIORITY: first-come-first-served also meets the floor by
    # construction when the budget allows, so H1a shows the selector loses no
    # time, not that it beats the baselines (their pass rates are reported).
    # H1b, IMAGE ECONOMY (added 2026-10-12, after the 15 runs of that day, so
    # confirmatory only on fresh data): at >= h1_min_onset_pct of onsets where
    # first-come-first-served needs two or more images to reach the level, the
    # earliest-completion cover reaches it within the same margin of the floor
    # using FEWER images than first-come-first-served. Mosaic quality at that
    # level (time-aware vs first-come-first-served) is reported, not tested.
    # Evaluability: a share of onsets can only show ">= h1_min_onset_pct" with
    # 95% confidence if enough onsets are evaluated -- with every onset a pass,
    # the exact one-sided 95% lower bound (Clopper & Pearson 1934) is
    # 0.05**(1/n), which reaches 0.90 only from n = 29. Below that the test is
    # NOT EVALUABLE. None derives n from h1_min_onset_pct; an integer overrides.
    h1_min_evaluable_onsets: Optional[int] = None

    # What the selector may know about cloud when it chooses passes.
    #   "realistic" (default): no image's cloud is known before that image
    #     ARRIVES (delivery time). How the remaining budget is planned is set by
    #     cloud_planning below (v1.5.0 default: expected coverage from each
    #     tile's cloud probability, re-planned after every arrival; v1.4.0: clear-
    #     sky geometry, re-planned when a planned observation comes back cloudy,
    #     the Fail-Signal). After the budget, compensation takes each arriving
    #     image that clearly covers a still-missing tile.
    #   "oracle": the selector knows in advance which tiles each pass will see
    #     clear (the behaviour of every run before 2026-10-09). Not achievable
    #     in practice; reported as a clairvoyant reference (greedy, so not a
    #     strict upper bound).
    cloud_knowledge: str = "realistic"

    # Time-aware objective (v1.5.0, Cells 6 and 7C).
    #   "facility" (default): greedy on the time-weighted coverage
    #     f(S) = sum over tiles of max over selected passes e of w_e * u_e,t,
    #     w_e = exp(-rate * t_e / tau). f is a weighted facility-location
    #     function, hence monotone submodular, so the k-step greedy keeps the
    #     1 - 1/e guarantee of Nemhauser, Wolsey & Fisher (1978), Math. Program.
    #     14, 265-294, doi:10.1007/BF01588971; Cell 7C measures its exact gap
    #     with the same ILP as coverage-first, utilities replaced by w * u.
    #   "legacy": v1.4.0 (unweighted marginal gain times w_e). Those scores are
    #     not the marginal gains of any set function, so no guarantee applies.
    time_aware_objective: str = "facility"

    # Cloud model (v1.5.0, Cell 6).
    #   "correlated" (default): a latent Gaussian field over the region,
    #     spatially correlated with e-folding length cloud_corr_length_km and
    #     evolving between pass times as an AR(1) process with e-folding time
    #     cloud_corr_time_h (separable covariance). A tile is cloudy for a pass
    #     when the field there exceeds the threshold that makes its marginal
    #     probability exactly the tile's cloud probability. Passes minutes apart
    #     therefore see nearly the same clouds, and clouds cover neighbouring
    #     tiles together.
    #   "iid": every (pass, tile) pair independent (v1.4.0 and earlier).
    # Why: about half of global cloud cover comes from clouds 200 km or larger
    # and about 15% from clouds under 10 km (Wood & Field 2011, J. Climate 24,
    # 4800-4816, doi:10.1175/2011JCLI4056.1), so one independent draw per
    # ~2.4 km tile is unrealistic. The two scales are ASSUMPTIONS within that
    # range, not measured values: Cell 7C sweeps them (cloud_sweep).
    cloud_model: str = "correlated"
    cloud_corr_length_km: float = 50.0
    cloud_corr_time_h: float = 3.0
    # Each tile's cloud probability.
    #   "constant": cloud_probability everywhere (all runs before v1.5.0).
    #   "earthenv": the 2000-2014 mean cloud frequency of the analysis month
    #     at the tile (EarthEnv Global 1-km Cloud Frequency v1; Wilson & Jetz
    #     2016, PLoS Biol 14(3): e1002415, doi:10.1371/journal.pbio.1002415;
    #     licence CC BY-NC 4.0), read from
    #     CLIMATOLOGY_DIR/MODCF_monthlymean_MM.tif (Cell 6 downloads it there
    #     once if missing). cloud_climatology_month overrides the month.
    cloud_probability_source: str = "constant"
    cloud_climatology_month: Optional[int] = None
    # How the adaptive (Fail-Signal) selector plans the images it has not yet
    # seen (v1.5.0, Cell 6).
    #   "expected" (default): plan the remaining budget to maximise EXPECTED
    #     coverage of the tiles still missing, from each tile's cloud
    #     probability (redundancy where cloud is likely), and re-plan after
    #     every arrival from the coverage actually achieved.
    #   "clear": plan on clear-sky geometry and re-plan only when a planned
    #     observation comes back cloudy (v1.4.0).
    cloud_planning: str = "expected"
    # Cloud realisations for the cloud analyses of Cell 7C (H2, seeds
    # random_seed, random_seed + 1, ...). The maps and the onset analysis use
    # the single realisation random_seed.
    cloud_seeds: int = 20
    # Sensitivity of the cloud results to the two correlation scales (Cell 7C):
    # lengths x times, cloud_sweep_seeds realisations each; empty to skip.
    cloud_sweep_lengths_km: Tuple[float, ...] = (0.0, 10.0, 25.0, 50.0, 100.0, 200.0)
    cloud_sweep_times_h: Tuple[float, ...] = (1.0, 3.0, 6.0, 24.0)
    cloud_sweep_seeds: int = 10

    # Retrospective runs of a REAL fleet (not Starlink): simulate from this past
    # UTC instant (ISO 8601, e.g. "2026-09-15T00:00:00") instead of "now". The
    # fleet's satellites are taken from the current CelesTrak list, and their
    # orbital elements AS OF this instant (newest element set with epoch at or
    # before it, at most 7 days old) from Space-Track gp_history (needs the
    # SPACETRACK_USERNAME / SPACETRACK_PASSWORD Colab secrets). A past window
    # can be checked at once against the real image archive (Cell 7D).
    # None = live run from the newest elements.
    analysis_epoch: Optional[str] = None
    # Live runs (v1.5.1): start the window at this UTC instant (ISO 8601) instead
    # of at the newest element set. Runs that will be COMPARED (one region with
    # different fleets, say) must share one window: give them the same value
    # here, or use analysis_epoch. None = the newest element-set epoch over BOTH
    # CelesTrak groups (resource and planet) whichever fleet is run, so runs
    # fetched together share their start. The 15 runs of 12 October 2026 did
    # not: the CKGR open-data run started 4 h 51 min before the open + PlanetScope
    # run, so their times to coverage are not comparable as they stand.
    reference_epoch_utc: Optional[str] = None
    # Strip assignment (v1.5.1, Cell 5). False (default): each tile observation
    # belongs to the one-minute strip that contains its closest approach.
    # True: the v1.5.0-and-earlier rule, under which every strip also claimed
    # an "end cap" of radius half the swath beyond each end (most Sentinel-2
    # observations appeared in two strips). Only to reproduce earlier runs.
    legacy_strip_end_caps: bool = False
    # Fleet sensitivity (v1.5.2, Cell 2): satellites to leave out, by CelesTrak
    # OBJECT_NAME (case-insensitive, exact), e.g. ("SENTINEL-2A",). Being tracked
    # is not the same as acquiring: Sentinel-2A flies an extended campaign with a
    # reduced acquisition plan (it imaged the CKGR every 20 days in August-October
    # 2026, Earth Search catalogue), so its predicted passes are an upper bound;
    # run with and without it to bracket the open-data results. Added to the
    # output folder name (_without-...).
    exclude_satellites: Tuple[str, ...] = ()

    def __post_init__(self):
        """Validate configuration parameters upon instantiation."""
        self._validate()

    def _validate(self):
        if self.active_roi not in self.VALID_ROI_NAMES:
            raise ConfigurationError(
                f"ACTIVE_ROI={self.active_roi!r} is invalid. Valid ROIs: {sorted(self.VALID_ROI_NAMES)}"
            )
        if self.active_n not in self.VALID_N_VALUES:
            raise ConfigurationError(
                f"ACTIVE_N={self.active_n} is invalid. Valid scenarios: {sorted(self.VALID_N_VALUES)}"
            )
        if not 0.0 <= self.cloud_probability <= 1.0:
            raise ConfigurationError(f"CLOUD_PROBABILITY={self.cloud_probability} must be in [0.0, 1.0].")
        if self.max_selection_steps < 1:
            raise ConfigurationError(f"MAX_SELECTION_STEPS={self.max_selection_steps} must be >= 1.")
        if not isinstance(self.random_seed, int):
            raise ConfigurationError(f"RANDOM_SEED must be an integer, received {type(self.random_seed)}.")
        if self.window_hours < 1:
            raise ConfigurationError(f"WINDOW_HOURS={self.window_hours} must be >= 1.")
        if self.fleet not in self.VALID_FLEETS:
            raise ConfigurationError(f"FLEET={self.fleet!r} is invalid. Valid: {sorted(self.VALID_FLEETS)}")
        if self.planet_swath_case not in ("conservative", "optimistic"):
            raise ConfigurationError(f"PLANET_SWATH_CASE={self.planet_swath_case!r} must be 'conservative' or 'optimistic'.")
        if not -10.0 <= self.assumed_min_sun_elev_deg < 90.0:
            raise ConfigurationError(f"ASSUMED_MIN_SUN_ELEV_DEG={self.assumed_min_sun_elev_deg} must be in [-10, 90).")
        if not 0.0 < self.swath_half_angle_deg < 90.0:
            raise ConfigurationError(f"SWATH_HALF_ANGLE_DEG={self.swath_half_angle_deg} must be strictly in (0, 90).")
        if not 0 <= self.h3_resolution <= 15:
            raise ConfigurationError(f"H3_RESOLUTION={self.h3_resolution} must be an integer in [0, 15].")
        if not 0.0 <= self.downlink_min_elevation_deg < 90.0:
            raise ConfigurationError(f"DOWNLINK_MIN_ELEVATION_DEG={self.downlink_min_elevation_deg} must be in [0, 90).")
        if self.t_max_wait_minutes <= 0:
            raise ConfigurationError(f"T_MAX_WAIT_MINUTES={self.t_max_wait_minutes} must be > 0.")
        if self.latency_decay_rate < 0:
            raise ConfigurationError(f"LATENCY_DECAY_RATE={self.latency_decay_rate} must be >= 0.")
        if not self.decay_sweep or any(r < 0 for r in self.decay_sweep):
            raise ConfigurationError(f"DECAY_SWEEP={self.decay_sweep} must be non-empty with all rates >= 0.")
        if not 0 < self.onset_shift_hours < self.window_hours:
            raise ConfigurationError(f"ONSET_SHIFT_HOURS={self.onset_shift_hours} must be in (0, WINDOW_HOURS).")
        if self.onset_shift_step_min < 1:
            raise ConfigurationError(f"ONSET_SHIFT_STEP_MIN={self.onset_shift_step_min} must be >= 1.")
        if self.random_baseline_draws < 1:
            raise ConfigurationError(f"RANDOM_BASELINE_DRAWS={self.random_baseline_draws} must be >= 1.")
        if self.ilp_instances < 1 or self.ilp_budget_k < 1 or self.ilp_max_pairs < 1:
            raise ConfigurationError("ILP_INSTANCES, ILP_BUDGET_K and ILP_MAX_PAIRS must all be >= 1.")
        if not 0 < self.ilp_horizon_hours <= self.window_hours:
            raise ConfigurationError(f"ILP_HORIZON_HOURS={self.ilp_horizon_hours} must be in (0, WINDOW_HOURS].")
        if self.ilp_time_limit_s <= 0:
            raise ConfigurationError(f"ILP_TIME_LIMIT_S={self.ilp_time_limit_s} must be > 0.")
        if not self.coverage_targets_pct or any(not 0 < t <= 100 for t in self.coverage_targets_pct):
            raise ConfigurationError(f"COVERAGE_TARGETS_PCT={self.coverage_targets_pct} must be non-empty, each in (0, 100].")
        if not 0 < self.h1_reachable_coverage_pct <= 100:
            raise ConfigurationError(f"H1_REACHABLE_COVERAGE_PCT={self.h1_reachable_coverage_pct} must be in (0, 100].")
        if self.h1_floor_margin_min is not None and self.h1_floor_margin_min < 0:
            raise ConfigurationError(f"H1_FLOOR_MARGIN_MIN={self.h1_floor_margin_min} must be >= 0 or None.")
        if not 0 < self.h1_min_onset_pct <= 100:
            raise ConfigurationError(f"H1_MIN_ONSET_PCT={self.h1_min_onset_pct} must be in (0, 100].")
        if self.h1_min_evaluable_onsets is not None and int(self.h1_min_evaluable_onsets) < 1:
            raise ConfigurationError(f"H1_MIN_EVALUABLE_ONSETS={self.h1_min_evaluable_onsets} must be >= 1 or None.")
        if self.latency_decay_scale not in ("fleet", "fixed"):
            raise ConfigurationError(f"LATENCY_DECAY_SCALE={self.latency_decay_scale!r} must be 'fleet' or 'fixed'.")
        if self.cloud_knowledge not in ("realistic", "oracle"):
            raise ConfigurationError(f"CLOUD_KNOWLEDGE={self.cloud_knowledge!r} must be 'realistic' or 'oracle'.")
        if self.time_aware_objective not in ("facility", "legacy"):
            raise ConfigurationError(f"TIME_AWARE_OBJECTIVE={self.time_aware_objective!r} must be 'facility' or 'legacy'.")
        if self.cloud_model not in ("correlated", "iid"):
            raise ConfigurationError(f"CLOUD_MODEL={self.cloud_model!r} must be 'correlated' or 'iid'.")
        if not self.cloud_corr_length_km > 0:
            raise ConfigurationError(f"CLOUD_CORR_LENGTH_KM={self.cloud_corr_length_km} must be > 0 "
                                     "(use cloud_model='iid' for independent tiles).")
        if self.cloud_corr_time_h < 0:
            raise ConfigurationError(f"CLOUD_CORR_TIME_H={self.cloud_corr_time_h} must be >= 0 "
                                     "(0 = independent between pass times).")
        if self.cloud_planning not in ("expected", "clear"):
            raise ConfigurationError(f"CLOUD_PLANNING={self.cloud_planning!r} must be 'expected' or 'clear'.")
        if self.cloud_probability_source not in ("constant", "earthenv"):
            raise ConfigurationError(f"CLOUD_PROBABILITY_SOURCE={self.cloud_probability_source!r} must be "
                                     "'constant' or 'earthenv'.")
        if self.cloud_climatology_month is not None and not 1 <= int(self.cloud_climatology_month) <= 12:
            raise ConfigurationError(f"CLOUD_CLIMATOLOGY_MONTH={self.cloud_climatology_month} must be 1-12 or None.")
        if self.cloud_seeds < 2 or self.cloud_sweep_seeds < 2:
            raise ConfigurationError("CLOUD_SEEDS and CLOUD_SWEEP_SEEDS must be >= 2 (a confidence interval needs "
                                     "at least two realisations).")
        if any(x < 0 for x in self.cloud_sweep_lengths_km) or any(x < 0 for x in self.cloud_sweep_times_h):
            raise ConfigurationError("CLOUD_SWEEP_LENGTHS_KM and CLOUD_SWEEP_TIMES_H must be >= 0 "
                                     "(length 0 = independent tiles).")
        for _name in ("analysis_epoch", "reference_epoch_utc"):
            _v = getattr(self, _name)
            if _v is None:
                continue
            try:
                _ae = _utc_naive(_v)
            except (ValueError, TypeError):
                raise ConfigurationError(f"{_name.upper()}={_v!r} is not an ISO 8601 time.")
            if self.fleet == "starlink":
                raise ConfigurationError(f"{_name.upper()} applies to real fleets only; Starlink historical runs use "
                                         "the ROI's historical_epoch.")
            if _name == "analysis_epoch" and _ae > _utc_now_naive():
                raise ConfigurationError(f"ANALYSIS_EPOCH={_v} is in the future.")
        if self.analysis_epoch is not None and self.reference_epoch_utc is not None:
            raise ConfigurationError("Set ANALYSIS_EPOCH (retrospective, elements as of then) or REFERENCE_EPOCH_UTC "
                                     "(live elements, common window start), not both.")
        if not isinstance(self.legacy_strip_end_caps, bool):
            raise ConfigurationError("LEGACY_STRIP_END_CAPS must be True or False.")
        if (isinstance(self.exclude_satellites, str) or not isinstance(self.exclude_satellites, (tuple, list))
                or any(not isinstance(n, str) or not n.strip() for n in self.exclude_satellites)):
            raise ConfigurationError(f"EXCLUDE_SATELLITES={self.exclude_satellites!r} must be a tuple of satellite "
                                     "names (CelesTrak OBJECT_NAME), e.g. ('SENTINEL-2A',).")
        if self.exclude_satellites and self.fleet == "starlink":
            raise ConfigurationError("EXCLUDE_SATELLITES applies to real fleets only.")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)

# =============================================================================
# CONFIGURATION INITIALISATION
# =============================================================================

# Primary Configuration Instantiation. v1.5.2 default scenario: the Botswana
# study area (CKGR + Okavango Delta) with the open-data + PlanetScope fleet,
# conservative Dove swath -- the primary design of the confirmatory protocol.
# The Starlink proxy (fleet="starlink", active_n) and the Nepal corridor /
# Botswana bounding-box ROIs reproduce the earliest runs only.
config = ARFConfig(
    active_n=1000,
    active_roi="ckgr_okavango",
    cloud_probability=0.35,
    max_selection_steps=30,
    random_seed=42,
    fleet="earth-resources+planet",
    planet_swath_case="conservative",
    assumed_min_sun_elev_deg=5.0,
    window_hours=48,
    swath_half_angle_deg=9.22,
    h3_resolution=7,
    downlink_min_elevation_deg=15.0,
    t_max_wait_minutes=95.0,
    latency_decay_rate=3.0,
    decay_sweep=(0.5, 1.0, 2.0, 3.0, 5.0),
    latency_decay_scale="fleet",   # "fixed" reproduces v1.3.1 (one orbital period)
    onset_robustness_enabled=True,
    onset_shift_hours=24,
    onset_shift_step_min=20,
    random_baseline_draws=10,
    ilp_enabled=True,
    ilp_instances=5,
    ilp_horizon_hours=6,
    ilp_budget_k=10,
    ilp_time_limit_s=120.0,
    ilp_max_pairs=100_000,
    coverage_targets_pct=(95.0, 99.0),
    h1_reachable_coverage_pct=95.0,
    h1_floor_margin_min=None,      # None = one orbital period at the fleet's median altitude
    h1_min_onset_pct=90.0,
    h1_min_evaluable_onsets=None,  # None = derived (29 for 90%): fewer evaluable onsets -> NOT EVALUABLE
    cloud_knowledge="realistic",   # "oracle" reproduces the pre-2026-10-09 (clairvoyant) cloud runs
    time_aware_objective="facility",  # "legacy" reproduces v1.4.0's time-aware selection
    cloud_model="correlated",      # "iid" reproduces v1.4.0's cloud mask
    cloud_corr_length_km=50.0,     # ASSUMPTION, swept in Cell 7C
    cloud_corr_time_h=3.0,         # ASSUMPTION, swept in Cell 7C
    cloud_probability_source="constant",  # "earthenv": per-tile monthly climatology (Wilson & Jetz 2016)
    cloud_climatology_month=None,  # None = month of the analysis start
    cloud_planning="expected",     # "clear" reproduces v1.4.0's Fail-Signal planning
    cloud_seeds=20,
    cloud_sweep_lengths_km=(0.0, 10.0, 25.0, 50.0, 100.0, 200.0),
    cloud_sweep_times_h=(1.0, 3.0, 6.0, 24.0),
    cloud_sweep_seeds=10,
    analysis_epoch=None,           # e.g. "2026-09-15T00:00:00" for a retrospective real-fleet run
    reference_epoch_utc=None,      # live runs to be compared: one common window start (v1.5.1)
    legacy_strip_end_caps=False,   # True reproduces the end-cap strips of v1.5.0 and earlier (Cell 5)
    exclude_satellites=(),         # e.g. ("SENTINEL-2A",): fleet sensitivity run (v1.5.2, Cell 2)
)

# Batch runs (v1.5.0): a driver may set ARF_CONFIG_OVERRIDES = {field: value}
# in the namespace before running the cells; they replace the values above and
# are validated by the same checks. Recorded in the manifest (Cell 9).
ARF_CONFIG_OVERRIDES_APPLIED = dict(globals().get("ARF_CONFIG_OVERRIDES") or {})
if ARF_CONFIG_OVERRIDES_APPLIED:
    from dataclasses import replace as _dc_replace, fields as _dc_fields
    _unknown = sorted(set(ARF_CONFIG_OVERRIDES_APPLIED) - {f.name for f in _dc_fields(ARFConfig)})
    if _unknown:
        raise ConfigurationError(f"ARF_CONFIG_OVERRIDES has unknown fields: {_unknown}")
    config = _dc_replace(config, **{k: (tuple(v) if isinstance(v, list) else v)
                                    for k, v in ARF_CONFIG_OVERRIDES_APPLIED.items()})
    print(f"⚙️  Config overrides applied: {ARF_CONFIG_OVERRIDES_APPLIED}")

# Expose module-level globals for downstream pipeline cells
ACTIVE_N = config.active_n
ACTIVE_ROI = config.active_roi
CLOUD_PROBABILITY = config.cloud_probability
MAX_SELECTION_STEPS = config.max_selection_steps
RANDOM_SEED = config.random_seed
WINDOW_HOURS = config.window_hours
SWATH_HALF_ANGLE_DEG = config.swath_half_angle_deg
FLEET = config.fleet
PLANET_SWATH_CASE = config.planet_swath_case
ASSUMED_MIN_SUN_ELEV_DEG = config.assumed_min_sun_elev_deg
EXCLUDE_SATELLITES = tuple(str(n).strip().upper() for n in config.exclude_satellites)
FLEET_TAG = f"N{ACTIVE_N}" if FLEET == "starlink" else FLEET.replace("+", "_") + (
    f"-{PLANET_SWATH_CASE}" if "planet" in FLEET else "")
# Non-default windows get their own output folder, so a 16-day run never
# overwrites the 48 h run of the same region and fleet.
if config.window_hours != 48:
    FLEET_TAG += f"_{config.window_hours}h"
if config.analysis_epoch is not None:          # retrospective runs get their own folder too
    FLEET_TAG += "_from" + _utc_naive(config.analysis_epoch).strftime("%Y%m%d")
if config.reference_epoch_utc is not None:     # live runs with a common start (v1.5.1)
    FLEET_TAG += "_at" + _utc_naive(config.reference_epoch_utc).strftime("%Y%m%dT%H%M")
if EXCLUDE_SATELLITES:                          # fleet sensitivity runs (v1.5.2)
    import re as _re0
    FLEET_TAG += "_without-" + "-".join(_re0.sub(r"[^A-Z0-9]+", "", n) for n in sorted(EXCLUDE_SATELLITES))
# Scenario identity in the scenario_results ledger: (ROI, N) for the Starlink
# proxy, as before; (ROI, fleet tag) for real fleets (v1.5.2; was (ROI, fleet),
# so the two Planet swath cases, a 384 h window or a past window of the same
# fleet overwrote each other's entry when run in one kernel).
SCENARIO_N = ACTIVE_N if FLEET == "starlink" else FLEET_TAG

# -----------------------------------------------------------------------------
# MISSION TABLE (real fleets): every value carries its status and source.
#   sourced    -- stated on the cited page (quote recorded)
#   assumption -- no source found; labelled, and to be covered by a sweep
# name_pattern is matched against CelesTrak OBJECT_NAME. The patterns
# themselves are unconfirmed until a live fetch: Cell 2 lists every
# satellite left unmatched, so a missed name shows up rather than vanishing.
# -----------------------------------------------------------------------------
CELESTRAK_GROUPS = {
    # GROUP values of https://celestrak.org/NORAD/elements/gp.php ("resource" =
    # Earth Resources, "planet" = Planet), as listed in
    # https://resources.wolframcloud.com/FunctionRepository/resources/CelesTrakTLE
    "earth-resources": ["resource"],
    "planet": ["planet"],
    "earth-resources+planet": ["resource", "planet"],
}
_ACCESSED = "2026-10-08"
MISSION_TABLE = [
    {"mission": "Sentinel-2", "name_pattern": r"^SENTINEL[- ]?2[A-Z]?$",
     "swath_km": {"value": 290.0, "status": "sourced", "accessed": _ACCESSED,
                  "source_url": "https://sentiwiki.copernicus.eu/web/s2-mission",
                  "quote": "The orbital swath width is 290 km."},
     "min_sun_elev_deg": {"value": ASSUMED_MIN_SUN_ELEV_DEG, "status": "assumption",
                          "note": "Source states MSI collects reflected sunlight but gives no sun-elevation "
                                  "limit; the Landsat 5 deg criterion is borrowed (ASSUMED_MIN_SUN_ELEV_DEG)."}},
    {"mission": "Landsat 8/9", "name_pattern": r"^LANDSAT[- ]?[89]$",
     "swath_km": {"value": 185.0, "status": "sourced", "accessed": _ACCESSED,
                  "source_url": "https://science.nasa.gov/mission/landsat-9",
                  "quote": "Swath width: 185 km (115 mi)"},
     "min_sun_elev_deg": {"value": 5.0, "status": "sourced", "accessed": _ACCESSED,
                          "source_url": "https://www.usgs.gov/landsat-missions/landsat-9",
                          "quote": "including islands at solar elevation angles greater than 5 degrees"}},
    # PlanetScope swath is given by Planet as a FIELD OF VIEW ACROSS THE SWATH
    # (PSScene spec, June 2026, Table 2-A): PS2 and PS2.SD "3.0° (swath)",
    # PSB.SD "4.0° (swath)"; the frame sizes in the same table (24 x 8 km,
    # 24 x 16 km, 32.5 x 19.6 km at 475 km) are swath x scene length, which the
    # FOV confirms (3.0 deg at 475 km = 24.9 km; 1.0 deg = 8.3 km). An earlier
    # version of this table used 8 km -- the along-track scene length -- as the
    # swath; that was wrong. PS2 (Dove Classic) imagery ended April 2022, so the
    # current fleet is PS2.SD (Dove-R, 3.0 deg) or PSB.SD (SuperDove, 4.0 deg).
    # Which Dove carries which instrument cannot be told from its name, so the
    # two cases bracket it. Each satellite's swath is computed from the FOV at
    # its own altitude (spherical Earth), not taken at the 475 km reference.
    {"mission": "PlanetScope", "name_pattern": r"^FLOCK",
     "fov_deg": ({"value": 3.0, "status": "sourced (lower bound: PS2.SD Dove-R)", "accessed": "2026-10-09",
                  "source_url": "https://assets.planet.com/docs/Planet_PSScene_Imagery_Product_Spec_letter_screen.pdf",
                  "quote": "Field of View 3.0° (swath) 1.0° (scene length) 3.0° (swath) 2.0° (scene length) "
                           "4.0° (swath) 2.3° (scene length) [PS2, PS2.SD, PSB.SD]"}
                 if PLANET_SWATH_CASE == "conservative" else
                 {"value": 4.0, "status": "sourced (upper bound: PSB.SD SuperDove)", "accessed": "2026-10-09",
                  "source_url": "https://assets.planet.com/docs/Planet_PSScene_Imagery_Product_Spec_letter_screen.pdf",
                  "quote": "Field of View 3.0° (swath) 1.0° (scene length) 3.0° (swath) 2.0° (scene length) "
                           "4.0° (swath) 2.3° (scene length) [PS2, PS2.SD, PSB.SD]"}),
     "min_sun_elev_deg": {"value": 10.0, "status": "sourced", "accessed": _ACCESSED,
                          "source_url": "https://docs.planet.com/data/imagery/planetscope/",
                          "quote": "All of the Doves have a system constraint to image regions that have a "
                                   "sun elevation angle greater than 10°."}},
]
# -----------------------------------------------------------------------------
# ACQUISITION EXCLUSIONS (v1.5.2): periods in which a satellite of the fleet is
# tracked (it has orbital elements, so Cell 2 keeps it) but its operator states
# that it does not acquire nominal imagery. Cell 5 drops the passes whose time
# falls in [start_utc, end_utc) and counts them. Every entry carries its source.
# Not covered here: routine acquisition plans, e.g. the Sentinel-2A extension
# campaign (from 13 March 2025; "through 31 December 2027", Copernicus Data
# Space Ecosystem news, 9 October 2026), under which Sentinel-2A images a given
# region on a reduced schedule (over the CKGR: every 20 days in August-October
# 2026 in the Earth Search catalogue). Use exclude_satellites (above) to
# bracket such a satellite. Checked and not listed: the Sentinel-2A night-time
# campaign from 23 July 2026 (20 days), whose page states that "The S2A
# expected data volume during that period should not be impacted" (daytime
# imaging continued: an S2A scene of the CKGR on 4 August 2026 is in the archive).
# -----------------------------------------------------------------------------
ACQUISITION_EXCLUSIONS = [
    {"name_pattern": r"^SENTINEL[- ]?2A$", "start_utc": "2026-01-13T00:00:00", "end_utc": "2026-03-09T00:00:00",
     "status": "sourced", "accessed": "2026-10-14",
     "source_url": "https://sentinels.copernicus.eu/-/sentinel-2a-and-sentinel-2b-tandem-acquisitions",
     "quote": "From 13 January to 08 March 2026, Sentinel-2A will temporarily pause its nominal acquisition "
              "activities.",
     "note": "Sentinel-2A/2B tandem calibration (page published 18 December 2025). The end date is inclusive, so "
             "the window ends at 9 March 00:00 UTC. Images taken in vicarious (VIC) mode during the drift were "
             "distributed to expert users only and are not modelled."},
]
H3_RESOLUTION = config.h3_resolution
DOWNLINK_MIN_ELEVATION_DEG = config.downlink_min_elevation_deg
T_MAX_WAIT_MINUTES = config.t_max_wait_minutes
LATENCY_DECAY_RATE = config.latency_decay_rate
LATENCY_DECAY_SCALE = config.latency_decay_scale
DECAY_SWEEP = tuple(config.decay_sweep)
ONSET_ROBUSTNESS_ENABLED = config.onset_robustness_enabled
ONSET_SHIFT_HOURS = config.onset_shift_hours
ONSET_SHIFT_STEP_MIN = config.onset_shift_step_min
RANDOM_BASELINE_DRAWS = config.random_baseline_draws
ILP_ENABLED = config.ilp_enabled
ILP_INSTANCES = config.ilp_instances
ILP_HORIZON_HOURS = config.ilp_horizon_hours
ILP_BUDGET_K = config.ilp_budget_k
ILP_TIME_LIMIT_S = config.ilp_time_limit_s
ILP_MAX_PAIRS = config.ilp_max_pairs
COVERAGE_TARGETS_PCT = tuple(sorted(float(t) for t in config.coverage_targets_pct))
H1_REACHABLE_COVERAGE_PCT = config.h1_reachable_coverage_pct
H1_FLOOR_MARGIN_MIN_CONFIG = config.h1_floor_margin_min
H1_MIN_ONSET_PCT = config.h1_min_onset_pct
H1_MIN_EVALUABLE_ONSETS_CONFIG = config.h1_min_evaluable_onsets
CLOUD_KNOWLEDGE = config.cloud_knowledge
TIME_AWARE_OBJECTIVE = config.time_aware_objective
CLOUD_MODEL = config.cloud_model
CLOUD_CORR_LENGTH_KM = config.cloud_corr_length_km
CLOUD_CORR_TIME_H = config.cloud_corr_time_h
CLOUD_PROBABILITY_SOURCE = config.cloud_probability_source
CLOUD_CLIMATOLOGY_MONTH = config.cloud_climatology_month
CLOUD_PLANNING = config.cloud_planning
CLOUD_SEEDS = config.cloud_seeds
CLOUD_SWEEP_LENGTHS_KM = tuple(config.cloud_sweep_lengths_km)
CLOUD_SWEEP_TIMES_H = tuple(config.cloud_sweep_times_h)
CLOUD_SWEEP_SEEDS = config.cloud_sweep_seeds
ANALYSIS_EPOCH = (None if config.analysis_epoch is None else
                  _utc_naive(config.analysis_epoch).strftime("%Y-%m-%dT%H:%M:%S.%f"))
REFERENCE_EPOCH_UTC = (None if config.reference_epoch_utc is None else
                       _utc_naive(config.reference_epoch_utc).strftime("%Y-%m-%dT%H:%M:%S.%f"))
LEGACY_STRIP_END_CAPS = bool(config.legacy_strip_end_caps)

# =============================================================================
# GOOGLE COLAB DRIVE MOUNT & OUTPUT DIRECTORY
# =============================================================================

try:
    import google.colab
    IN_COLAB = True
except ImportError:
    IN_COLAB = False

if IN_COLAB:
    from google.colab import drive
    drive.mount('/content/drive')
    OUT_DIR = f"/content/drive/MyDrive/ARF_Outputs_{ACTIVE_ROI}_{FLEET_TAG}"
    # v1.5.2: a complete run of ANOTHER release in this folder (its Cell 9
    # manifest names the release) is kept, never overwritten: this run writes
    # to <folder>_v<release digits> instead. Protects, e.g., the 12 October
    # 2026 runs of release 1.5.0, whose arf_refined_events.csv files are the
    # inputs of arf_colab_parity.py. A folder without a manifest (an incomplete
    # run) is reused as before; an unreadable manifest counts as another release.
    _prev_manifest = os.path.join(OUT_DIR, "arf_reproducibility_manifest.json")
    if os.path.exists(_prev_manifest):
        try:
            with open(_prev_manifest) as _f:
                _prev_release = json.load(_f).get("framework_version")
        except Exception:
            _prev_release = None
        if _prev_release != __version__:
            print(f"📁 {OUT_DIR} holds a complete run of release {_prev_release or '(unreadable manifest)'}: "
                  f"kept unchanged; this run writes to its own folder.")
            OUT_DIR += "_v" + __version__.replace(".", "")
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"📁 Google Drive mounted. Outputs will be saved to: {OUT_DIR}")
else:
    OUT_DIR = os.getcwd()
    print(f"📁 Running locally. Outputs will be saved to: {OUT_DIR}")

# =============================================================================
# RESULTS LEDGER INITIALISATION
# =============================================================================

# Preserve multi-scenario state in memory across cell runs
if 'scenario_results' not in globals():
    scenario_results: Dict[Tuple[str, Any], Dict[str, Any]] = {}
    print("📊 Initialised empty scenario_results ledger (fresh kernel).")
else:
    print(f"📊 Existing scenario_results preserved for (ROI, N or fleet tag) = "
          f"{sorted(scenario_results.keys(), key=str)}.")

# =============================================================================
# SUMMARY & PERSISTENCE
# =============================================================================

print("\n" + "=" * 60)
print(f"  {__framework_name__} v{__version__}")
print("=" * 60)
print(f"  Active Configuration Summary:")
print(f"  ─────────────────────────────────────────────────────────")
print(f"  Region (ROI):              {ACTIVE_ROI}")
if FLEET == "starlink":
    print(f"  Starlink proxy, N:         {ACTIVE_N}")
else:
    print(f"  Fleet:                     {FLEET}" + (f" (PlanetScope {PLANET_SWATH_CASE})" if "planet" in FLEET else ""))
    print(f"  Window start:              " + (f"{ANALYSIS_EPOCH} UTC (retrospective)" if ANALYSIS_EPOCH else
                                              f"{REFERENCE_EPOCH_UTC} UTC (set)" if REFERENCE_EPOCH_UTC else
                                              "newest orbital elements (live)"))
    if EXCLUDE_SATELLITES:
        print(f"  Left out (sensitivity):    {', '.join(EXCLUDE_SATELLITES)}")
print(f"  Cloud Obscuration Prob:    {CLOUD_PROBABILITY:.2f}")
print(f"  Selection Budget (k):      {MAX_SELECTION_STEPS}")
print(f"  Random Seed:               {RANDOM_SEED}")
print(f"  Simulation Window:         {WINDOW_HOURS} hours")
if FLEET == "starlink":
    print(f"  Sensor Half-Angle:         {SWATH_HALF_ANGLE_DEG}°")
print(f"  H3 DGGS Resolution:        {H3_RESOLUTION}")
print(f"  Downlink Min Elevation:    {DOWNLINK_MIN_ELEVATION_DEG}°")
print(f"  Working Directory:         {OUT_DIR}")
print(f"  ─────────────────────────────────────────────────────────")
print(f"  License:                   {__license__}")
print(f"  Recorded Scenarios:        {len(scenario_results)}")
print("=" * 60)

config_filepath = os.path.join(OUT_DIR, 'arf_config.json')
try:
    with open(config_filepath, 'w') as f:
        f.write(config.to_json())
    print(f"✅ Configuration exported to {config_filepath}")
except Exception as e:
    print(f"⚠️ Warning: Could not write configuration file: {e}")

print("\n=== CELL 0 COMPLETE ===")

