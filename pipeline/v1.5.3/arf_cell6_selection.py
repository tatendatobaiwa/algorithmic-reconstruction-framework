# =============================================================================
# CELL 6: SUBMODULAR GREEDY SELECTION WITH FAIL-SIGNAL GAP-FILLING
# =============================================================================
# Builds the candidate events (one-minute acquisition strips, Cell 5) and
# selects images:
#   - coverage-first: greedy on sum over tiles of the best utility, k images
#   - time-aware: greedy on the time-weighted coverage sum_t max_e w_e u_e,t
#     (TIME_AWARE_OBJECTIVE, Cell 0); both objectives are monotone submodular
#   - under cloud, the causal Fail-Signal selector: each image's cloud is known
#     only once it arrives; the plan for the remaining budget is revised from
#     what has been seen (CLOUD_PLANNING), then compensation (gap-fill)
#   - the clairvoyant greedy (cloud known in advance), reported for reference
# Utility of a tile observation: u = max(0, 1 - A_nadir / theta_max) * D_downlink
# (1 at nadir, 0 at the swath edge; D_downlink from Cell 5). This form was
# chosen for this pipeline (a design choice, not taken from a source).
# Scope of the (1 - 1/e) guarantee (v1.5.2 wording): it holds for the k-step
# greedy on a FIXED objective -- coverage-first and the facility time-aware
# objective, in clear sky or with the cloud known in advance (the oracle
# reference). The causal Fail-Signal selector re-plans after each image and
# adds compensatory images beyond k, so no approximation bound is claimed for
# it; its expected-coverage planner is itself a greedy on a monotone
# submodular function, but with the time-aware weighting it scores gain times
# weight (as the legacy objective did), a heuristic with no bound.
# Provenance: the cell was rebuilt in October 2026 from the contracts of the
# cells that use it (event keys (sat, t_min), cloud mask keyed by (event, tile)),
# when the working notebook was found to have lost it.
# =============================================================================

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 6 is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["6"] = _ARF_CELL_VERSION

logger = logging.getLogger("arf.cell6")

import random
import heapq   # Fail-Signal gap-fill (lazy-greedy max-heap)

class EventStore:
    """
    Compact store of candidate events (compressed sparse rows). Event e owns
    rows offsets[e]:offsets[e+1] of the parallel arrays tile_idx (int32 index
    into _tile_ids, Cell 4) and util (float64 utility). Keys are (sat, t_min)
    in sorted order -- the same order the previous dict-of-dicts had -- and
    each event's downlink wait is kept once rather than once per tile.

    Memory: ~12 bytes per tile-access row instead of ~100+ for a nested
    Python dict, which is what lets N=1000 country-scale runs fit in Colab.

    Read access used by downstream cells still works as before:
    len(event_data), event_data.keys(), key in event_data, and
    event_data[key].items() (materialised on demand for that one event).
    """
    def __init__(self, keys, offsets, tile_idx, util, waits):
        self.keys_ = list(keys)
        self.index = {k: i for i, k in enumerate(self.keys_)}
        self.offsets = np.asarray(offsets, dtype=np.int64)
        self.tile_idx = np.asarray(tile_idx, dtype=np.int32)
        self.util = np.asarray(util, dtype=np.float64)
        self.wait_lookup = dict(zip(self.keys_, waits))
        self.counts = np.diff(self.offsets)
        self.n_rows = int(self.offsets[-1]) if len(self.offsets) else 0
        # v1.5.2: every event owns at least one row. The selectors sum per event
        # with np.add.reduceat(x, offsets[:-1]), which for an empty segment
        # returns the NEXT row's value instead of 0 (NumPy reduceat semantics),
        # so an empty event would silently get another event's gain.
        if len(self.keys_) and not (self.counts.size == len(self.keys_) and (self.counts >= 1).all()
                                    and self.n_rows == self.tile_idx.size == self.util.size):
            raise ValueError("EventStore: every event needs >= 1 row and offsets must match the row arrays.")

    def __len__(self):
        return len(self.keys_)

    def __contains__(self, key):
        return key in self.index

    def keys(self):
        return list(self.keys_)

    def rows(self, e):
        return slice(int(self.offsets[e]), int(self.offsets[e + 1]))

    def __getitem__(self, key):
        s = self.rows(self.index[key])
        return dict(zip(_tile_ids[self.tile_idx[s]], self.util[s].tolist()))


# Rank of each ROI tile in sorted-ID order: reproduces the previous cloud
# mask's "sorted(tiles)" draw order exactly, without sorting strings per event.
_tile_string_rank = np.empty(len(_tile_ids), dtype=np.int64)
_tile_string_rank[np.argsort(_tile_ids, kind="stable")] = np.arange(len(_tile_ids))

# total_target_tiles is the module-level denominator run_submodular_greedy_selection's
# history entries use for coverage_pct. Cell 7 also computes this (as a
# convenience alias, using refined_events_df's tile-reach as a cross-check),
# but Cell 6 needs it available now, from Cell 1's target_tiles.
total_target_tiles = len(target_tiles)


def build_event_data(refined_df, theta_max=SWATH_HALF_ANGLE_DEG):
    """
    Build the EventStore from the refined events table (Cell 5's output).
    Utility combines incidence-angle geometry (1 at nadir, 0 at
    theta_max=SWATH_HALF_ANGLE_DEG -- the same swath-edge cutoff Cell 5
    already filters on) with D_downlink (already a 0..1 reward from Cell 5):
        utility = max(0, 1 - A_nadir/theta_max) * D_downlink
    Events are ordered by (sat, t_min); rows keep their Cell 5 order within
    each event, so greedy tie-breaks and sums match the previous dict form.
    """
    n = len(refined_df)
    if n == 0:
        return EventStore([], np.zeros(1, np.int64), np.zeros(0, np.int32), np.zeros(0), [])

    sat = refined_df["sat"]
    if not isinstance(sat.dtype, pd.CategoricalDtype):
        sat = sat.astype("category")
    sat = sat.cat.reorder_categories(sorted(sat.cat.categories))   # lexicographic key order
    sat_code = sat.cat.codes.to_numpy(np.int64)
    t_min = refined_df["t_min"].to_numpy(np.int64)

    tile = refined_df["tile"]
    if (isinstance(tile.dtype, pd.CategoricalDtype) and len(tile.cat.categories) == len(_tile_ids)
            and bool((np.asarray(tile.cat.categories) == _tile_ids).all())):
        tidx = tile.cat.codes.to_numpy(np.int32)
    else:
        pos = {t: i for i, t in enumerate(_tile_ids)}
        tidx = np.fromiter((pos[t] for t in tile), dtype=np.int32, count=n)

    # Real fleets: normalise to EACH satellite's own swath-edge angle.
    theta = (refined_df["theta_max_deg"].to_numpy(np.float64) if "theta_max_deg" in refined_df.columns
             else theta_max)
    util = ((1.0 - refined_df["A_nadir_deg"].to_numpy(np.float64) / theta).clip(min=0.0)
            * refined_df["D_downlink"].to_numpy(np.float64))
    wait = refined_df["T_wait_min"].to_numpy(np.float64)

    order = np.lexsort((np.arange(n), t_min, sat_code))            # stable within an event
    sc, tm = sat_code[order], t_min[order]
    starts = np.flatnonzero(np.r_[True, (sc[1:] != sc[:-1]) | (tm[1:] != tm[:-1])])
    cats = list(sat.cat.categories)
    keys = [(cats[sc[i]], int(tm[i])) for i in starts]
    w = wait[order][starts]
    waits = [None if np.isnan(x) else float(x) for x in w]
    return EventStore(keys, np.r_[starts, n], tidx[order], util[order], waits)


event_data = build_event_data(refined_events_df)
logger.info(f"Built event_data: {len(event_data)} candidate (sat, t_min) events, "
            f"{event_data.n_rows:,} tile-access rows.")


class CloudMask:
    """
    Per-row cloud flags aligned with an EventStore. .get((event_key, tile))
    is kept for downstream cells (Cell 7B), backed by a one-event cache.
    """
    def __init__(self, store, flags):
        self.store = store
        self.flags = np.asarray(flags, dtype=bool)
        self.n_cloudy = int(self.flags.sum())
        self._cache_e, self._cache = None, {}

    def __len__(self):
        return len(self.flags)

    def get(self, key, default=False):
        event_key, tile = key
        e = self.store.index.get(event_key)
        if e is None:
            return default
        if e != self._cache_e:
            s = self.store.rows(e)
            self._cache = dict(zip(_tile_ids[self.store.tile_idx[s]], self.flags[s].tolist()))
            self._cache_e = e
        return self._cache.get(tile, default)


# =============================================================================
# CLOUD MODEL (v1.5.0): see CLOUD_MODEL and CLOUD_PROBABILITY_SOURCE in Cell 0
# =============================================================================

def _default_climatology_dir():
    if os.path.isdir("/content/drive/MyDrive"):
        return "/content/drive/MyDrive/ARF_climatology"
    return os.path.join(os.getcwd(), "ARF_climatology")


CLIMATOLOGY_DIR = globals().get("CLIMATOLOGY_DIR") or _default_climatology_dir()
EARTHENV_URL = "https://data.earthenv.org/cloud/MODCF_monthlymean_{month:02d}.tif"
EARTHENV_CITATION = ("Wilson, A. M. & Jetz, W. (2016). Remotely sensed high-resolution global cloud dynamics for "
                     "predicting ecosystem and biodiversity distributions. PLoS Biol 14(3): e1002415, "
                     "doi:10.1371/journal.pbio.1002415. Data: EarthEnv Global 1-km Cloud Frequency v1 "
                     "(https://www.earthenv.org/cloud), licence CC BY-NC 4.0.")


def _tile_latlng():
    """(lat, lng) arrays of the ROI tile centres, aligned with _tile_ids (cached)."""
    key = (len(_tile_ids), str(_tile_ids[0]) if len(_tile_ids) else "", str(_tile_ids[-1]) if len(_tile_ids) else "")
    c = globals().get("_tile_latlng_cache")
    if c is None or c[0] != key:
        ll = np.array([h3.cell_to_latlng(t) for t in _tile_ids], dtype=np.float64).reshape(-1, 2)
        c = (key, ll[:, 0].copy(), ll[:, 1].copy())
        globals()["_tile_latlng_cache"] = c
    return c[1], c[2]


def _climatology_month():
    """(month, how it was chosen) for the per-tile cloud climatology."""
    if globals().get("CLOUD_CLIMATOLOGY_MONTH") is not None:
        return int(CLOUD_CLIMATOLOGY_MONTH), "set in Cell 0"
    if globals().get("ANALYSIS_EPOCH"):
        return int(str(ANALYSIS_EPOCH)[5:7]), "month of ANALYSIS_EPOCH"
    try:
        return int(str(_sim_epoch_datetime64())[5:7]), "month of the simulation start"
    except Exception:
        from datetime import datetime as _dt, timezone as _tz
        return _dt.now(_tz.utc).month, "current month (simulation start unavailable)"


def earthenv_cloud_probability(month, lat, lng, path=None):
    """
    Per-tile cloud probability from the EarthEnv 2000-2014 monthly mean cloud
    frequency (Wilson & Jetz 2016): the mean of the valid ~1-km pixels in the
    3 x 3 block around each tile centre (about one H3 res-7 tile, 5.2 km2).
    Stored values 0-10,000 are percent x 100; values above 10,000 are fill
    (EarthEnv data description, Cloud_DataDescription.pdf). Tiles with no valid
    pixel fall back to CLOUD_PROBABILITY and are counted. The file is read from
    CLIMATOLOGY_DIR, downloaded there once if missing. Returns (p, info).
    """
    try:
        import rasterio
    except ImportError:                     # optional dependency, only for cloud_probability_source="earthenv"
        import subprocess, sys as _sys
        subprocess.check_call([_sys.executable, "-m", "pip", "install", "-q", "rasterio"])
        import rasterio
    from rasterio.windows import Window
    path = path or os.path.join(CLIMATOLOGY_DIR, f"MODCF_monthlymean_{month:02d}.tif")
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        import urllib.request
        url = EARTHENV_URL.format(month=month)
        logger.info(f"Downloading the EarthEnv cloud climatology (month {month:02d}) to {path} (one time) ...")
        urllib.request.urlretrieve(url, path + ".part")
        os.replace(path + ".part", path)
    lat, lng = np.asarray(lat, dtype=np.float64), np.asarray(lng, dtype=np.float64)
    with rasterio.open(path) as src:
        epsg = src.crs.to_epsg() if src.crs is not None else 4326
        if epsg != 4326:
            from rasterio.warp import transform as _warp
            xs, ys = _warp("EPSG:4326", src.crs, lng.tolist(), lat.tolist())
            xs, ys = np.asarray(xs), np.asarray(ys)
        else:
            xs, ys = lng, lat
        inv = ~src.transform
        cols, rows = inv * (xs, ys)
        rows, cols = np.floor(rows).astype(np.int64), np.floor(cols).astype(np.int64)
        r0, c0 = max(int(rows.min()) - 1, 0), max(int(cols.min()) - 1, 0)
        r1, c1 = min(int(rows.max()) + 2, src.height), min(int(cols.max()) + 2, src.width)
        if r0 >= r1 or c0 >= c1:
            raise RuntimeError(f"The ROI lies outside {path}.")
        a = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype(np.float64)
        nodata = src.nodata
    bad = (a > 10000) | (a < 0) | ~np.isfinite(a)
    if nodata is not None:
        bad |= a == nodata
    a[bad] = np.nan
    vals = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            rr = np.clip(rows - r0 + di, 0, a.shape[0] - 1)
            cc = np.clip(cols - c0 + dj, 0, a.shape[1] - 1)
            vals.append(a[rr, cc])
    with np.errstate(all="ignore"):
        pct = np.nanmean(np.vstack(vals), axis=0) * 0.01           # % cloudy days
    miss = ~np.isfinite(pct)
    p = np.where(miss, float(CLOUD_PROBABILITY), pct / 100.0)
    p = np.clip(p, 0.001, 0.999)
    info = {"source": "earthenv", "month": int(month), "path": path, "citation": EARTHENV_CITATION,
            "n_tiles": int(p.size), "n_tiles_without_data": int(miss.sum()),
            "mean": float(p.mean()), "min": float(p.min()), "max": float(p.max())}
    return p, info


def cloud_probability_per_tile():
    """Marginal cloud probability of every ROI tile (aligned with _tile_ids) and how it was set."""
    if globals().get("CLOUD_PROBABILITY_SOURCE", "constant") == "earthenv":
        month, how = _climatology_month()
        lat, lng = _tile_latlng()
        p, info = earthenv_cloud_probability(month, lat, lng)
        info["month_from"] = how
        return p, info
    return (np.full(len(_tile_ids), float(CLOUD_PROBABILITY)),
            {"source": "constant", "value": float(CLOUD_PROBABILITY)})


def _gauss_kernel1d(s_cells, truncate=4.0):
    """Normalised 1-D Gaussian kernel, as scipy.ndimage.gaussian_filter1d builds it."""
    radius = int(truncate * s_cells + 0.5)
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    k = np.exp(-0.5 * (x / s_cells) ** 2)
    return k / k.sum()


class LatentCloudField:
    """
    Unit-variance Gaussian field on a local equirectangular plane around the ROI
    with correlation exp(-d^2 / L^2) (e-folding distance L = length_km): white
    noise filtered by a separable Gaussian kernel of standard deviation L/2,
    on a grid of spacing min(L/6, 10 km), padded by 4.5 kernel widths and
    wrapped (opposite ROI edges are then >= 4.5 L apart through the wrap).
    Values at tile centres are bilinear interpolations; their exact variance
    (computed from the discrete kernel) rescales each tile's threshold, so the
    marginal cloud probability of every tile is exactly its p.
    """
    def __init__(self, length_km):
        from scipy.ndimage import correlate1d
        from scipy.special import ndtri
        self._correlate1d, self._ndtri = correlate1d, ndtri
        lat, lng = _tile_latlng()
        R = 6371.0088
        lat0, lng0 = np.radians(lat.mean()), np.radians(lng.mean())
        x = R * (np.radians(lng) - lng0) * np.cos(lat0)
        y = R * (np.radians(lat) - lat0)
        self.length_km = float(length_km)
        sigma = self.length_km / 2.0
        self.dx = min(self.length_km / 6.0, 10.0)
        pad = 4.5 * sigma
        x0, y0 = x.min() - pad, y.min() - pad
        self.nx = int(np.ceil((x.max() + pad - x0) / self.dx)) + 2
        self.ny = int(np.ceil((y.max() + pad - y0) / self.dx)) + 2
        fx, fy = (x - x0) / self.dx, (y - y0) / self.dx
        self.ix, self.iy = np.floor(fx).astype(np.int64), np.floor(fy).astype(np.int64)
        wx, wy = fx - self.ix, fy - self.iy
        self.k1 = _gauss_kernel1d(sigma / self.dx)
        ss = float((self.k1 ** 2).sum())
        self.norm = ss                                       # sd of the 2-D filtered noise = sum(k1^2)
        c1 = float((self.k1[:-1] * self.k1[1:]).sum() / ss)  # correlation at a lag of one cell
        var_x = (1 - wx) ** 2 + wx ** 2 + 2 * (1 - wx) * wx * c1
        var_y = (1 - wy) ** 2 + wy ** 2 + 2 * (1 - wy) * wy * c1
        self.sd = np.sqrt(var_x * var_y)
        self.w00, self.w10 = (1 - wx) * (1 - wy), wx * (1 - wy)
        self.w01, self.w11 = (1 - wx) * wy, wx * wy

    def draw(self, rng):
        """One independent unit-variance field."""
        z = rng.standard_normal((self.ny, self.nx))
        z = self._correlate1d(z, self.k1, axis=0, mode="wrap")
        z = self._correlate1d(z, self.k1, axis=1, mode="wrap")
        return z / self.norm

    def at_tiles(self, z, tiles):
        ix, iy = self.ix[tiles], self.iy[tiles]
        return (self.w00[tiles] * z[iy, ix] + self.w10[tiles] * z[iy, ix + 1]
                + self.w01[tiles] * z[iy + 1, ix] + self.w11[tiles] * z[iy + 1, ix + 1])

    def thresholds(self, p_tile):
        """Per-tile threshold on the interpolated field giving marginal probability p_tile."""
        return self._ndtri(1.0 - np.asarray(p_tile, dtype=np.float64)) * self.sd


def _latent_field(length_km):
    cache = globals().setdefault("_latent_field_cache", {})
    key = (float(length_km), len(_tile_ids), str(_tile_ids[0]) if len(_tile_ids) else "",
           str(_tile_ids[-1]) if len(_tile_ids) else "")
    if key not in cache:
        cache[key] = LatentCloudField(length_km)
    return cache[key]


def build_cloud_mask(event_data, cloud_probability, random_seed, model=None, length_km=None, time_h=None,
                     p_tile=None):
    """
    Deterministic cloud mask over every (event, tile) row (seeded).
    model "iid" (v1.4.0 and earlier): events visited in sorted key order and
      tiles in sorted-ID order, drawing from random.Random(random_seed) -- the
      exact draw sequence of every earlier run when the probability is constant.
    model "correlated" (v1.5.0): the latent field (LatentCloudField) at each
      distinct pass time, evolved between times as AR(1),
      Z(t2) = rho Z(t1) + sqrt(1 - rho^2) E, rho = exp(-(t2 - t1) / time_h),
      E an independent field; time_h = 0 makes pass times independent. Passes at
      the same minute see the same clouds. Drawn from numpy's default_rng(seed).
    p_tile: per-tile probability (aligned with _tile_ids); default from
      CLOUD_PROBABILITY_SOURCE (constant cloud_probability, or the climatology).
    """
    model = model or globals().get("CLOUD_MODEL", "iid")
    store = event_data
    flags = np.zeros(store.n_rows, dtype=bool)
    if p_tile is None:
        p_tile = (CLOUD_P_TILE if globals().get("CLOUD_PROBABILITY_SOURCE", "constant") != "constant"
                  and "CLOUD_P_TILE" in globals() else None)
    if model == "iid":
        rng = random.Random(random_seed)
        for e in sorted(range(len(store)), key=lambda i: store.keys_[i]):
            s = store.rows(e)
            rows = np.arange(s.start, s.stop)
            rows = rows[np.argsort(_tile_string_rank[store.tile_idx[s]], kind="stable")]
            if p_tile is None:
                flags[rows] = [rng.random() < cloud_probability for _ in range(rows.size)]
            else:
                pr = p_tile[store.tile_idx[rows]]
                flags[rows] = [rng.random() < pr[i] for i in range(rows.size)]
        return CloudMask(store, flags)
    if store.n_rows == 0:
        return CloudMask(store, flags)
    length_km = float(CLOUD_CORR_LENGTH_KM if length_km is None else length_km)
    time_h = float(CLOUD_CORR_TIME_H if time_h is None else time_h)
    field = _latent_field(length_km)
    p = np.full(len(_tile_ids), float(cloud_probability)) if p_tile is None else np.asarray(p_tile, dtype=np.float64)
    thr = field.thresholds(p)
    t_event = np.array([k[1] for k in store.keys_], dtype=np.float64)
    t_row = np.repeat(t_event, store.counts)
    times, inv = np.unique(t_row, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order], np.arange(len(times) + 1))
    rng = np.random.default_rng(random_seed)
    z, t_prev = None, None
    for i, t in enumerate(times):
        e = field.draw(rng)
        if z is None:
            z = e
        else:
            rho = 0.0 if time_h <= 0 else float(np.exp(-(t - t_prev) / (60.0 * time_h)))
            z = rho * z + np.sqrt(max(0.0, 1.0 - rho * rho)) * e
        t_prev = t
        rows = order[bounds[i]:bounds[i + 1]]
        tiles = store.tile_idx[rows]
        flags[rows] = field.at_tiles(z, tiles) > thr[tiles]
    return CloudMask(store, flags)


CLOUD_P_TILE, cloud_climatology_info = cloud_probability_per_tile()
cloud_model_info = {"model": CLOUD_MODEL, "probability": cloud_climatology_info}
if CLOUD_MODEL == "correlated":
    _f = _latent_field(CLOUD_CORR_LENGTH_KM)
    cloud_model_info.update(length_km=float(CLOUD_CORR_LENGTH_KM), time_h=float(CLOUD_CORR_TIME_H),
                            grid_spacing_km=float(_f.dx), grid_cells=int(_f.nx * _f.ny),
                            reference="Wood, R. & Field, P. R. (2011). The distribution of cloud horizontal sizes. "
                                      "J. Climate 24, 4800-4816, doi:10.1175/2011JCLI4056.1 (scale range); "
                                      "the length and time scales are assumptions, swept in Cell 7C.")
cloud_mask = build_cloud_mask(event_data, CLOUD_PROBABILITY, RANDOM_SEED)
logger.info(f"Built cloud_mask ({CLOUD_MODEL}"
            + (f", L={CLOUD_CORR_LENGTH_KM:g} km, T={CLOUD_CORR_TIME_H:g} h" if CLOUD_MODEL == "correlated" else "")
            + f"): {cloud_mask.n_cloudy}/{len(cloud_mask)} tile-events cloud-obscured "
            + (f"({CLOUD_PROBABILITY:.0%} target)." if cloud_climatology_info["source"] == "constant" else
               f"(climatology month {cloud_climatology_info['month']:02d}, mean p = {cloud_climatology_info['mean']:.2f})."))


# =============================================================================
# DELIVERY ORDER AND THE TIME SCALE OF THE TIME-AWARE WEIGHT (v1.4.0)
# =============================================================================
# Delivery time = t_min + downlink wait. An image with no downlink contact in
# the window has no delivery time (and zero utility). Sorting on a NaN key
# leaves Python's sort unordered (sorted([3.0, nan, 1.0]) == [3.0, nan, 1.0]),
# so every delivery-order sort uses delivery_sort_key, which puts such images
# last. Real and test runs so far had no such image, so their results are
# unchanged; the fix matters for ground stations far from the ROI.
#
# Time-aware objective: each pass's marginal gain is weighted by
#     exp(-rate * t_min / tau).
# Up to v1.3.1 tau was one orbital period (T_MAX_WAIT_MINUTES). That suits a
# fleet that covers the ROI within about an orbit, but for a sparse fleet whose
# passes need a day to cover it, every pass after the first orbit is weighted
# almost to nothing (exp(-3 * 1440 / 94) ~ 1e-20), so the budget goes on the
# earliest small passes. Real run, CKGR + Okavango, Landsat + Sentinel-2 +
# PlanetScope, 2026-10-12: 91.0% clear-sky coverage at k = 30, against 100%
# for coverage-first. LATENCY_DECAY_SCALE (Cell 0):
#   "fleet" (default): tau = max(T_MAX_WAIT_MINUTES, F), where F is the
#       physical floor of the passes being selected from: the delivery time at
#       which ALL of them together cover H1_REACHABLE_COVERAGE_PCT of the tiles
#       they can reach, from clear-sky geometry (pass predictions, known before
#       any image is taken: no cloud information is used). A pass delivered at
#       the floor is weighted exp(-rate). If the fleet reaches that level within
#       one orbital period, tau is the orbital period, exactly as before.
#   "fixed": tau = T_MAX_WAIT_MINUTES (v1.3.1).
# tau is computed per EventStore, so each onset of Cell 7C's analysis uses the
# floor of its own passes.

def delivery_sort_key(store):
    """Sort key for (sat, t_min) keys of `store` by delivery time; images with
    no downlink in the window (no delivery time) sort last."""
    def key(k):
        w = store.wait_lookup.get(k)
        return np.inf if w is None or np.isnan(w) else k[1] + w
    return key


def floor_time(store, level_pct):
    """Physical floor: with ALL events of `store` in delivery order, the delivery
    time (min) at which they first cover level_pct of the tiles reachable from
    the store (utility > 0; clear-sky geometry). None if nothing is reachable.
    Same rule as Cell 7C's coverage_times(store, None, (level_pct,),
    denominator=reachable_tile_count(store, None))."""
    if len(store) == 0:
        return None
    usable = store.util > 0.0
    reach = int(np.unique(store.tile_idx[usable]).size)
    if reach == 0:
        return None
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n_seen = 0
    for k in sorted(store.keys_, key=delivery_sort_key(store)):
        w = store.wait_lookup.get(k)
        if w is None or np.isnan(w):
            break                                  # undelivered images sort last
        s = store.rows(store.index[k])
        t = store.tile_idx[s][usable[s]]
        new = t[~seen[t]]
        if new.size:
            seen[new] = True
            n_seen += new.size
            if 100.0 * n_seen / reach >= level_pct - 1e-9:
                return float(k[1] + w)
    return None


def decay_time_scale(store):
    """tau (minutes) of the time-aware weight for selections from `store`
    (see the section header). Cached on the store."""
    t_orbit = float(T_MAX_WAIT_MINUTES)
    if globals().get("LATENCY_DECAY_SCALE", "fleet") == "fixed":
        return t_orbit
    level = float(globals().get("H1_REACHABLE_COVERAGE_PCT", 95.0))
    cache = store.__dict__.setdefault("_tau_cache", {})
    if (level, t_orbit) not in cache:
        f = floor_time(store, level)
        cache[(level, t_orbit)] = t_orbit if f is None else max(t_orbit, f)
    return cache[(level, t_orbit)]


def time_aware_weights(store, latency_decay_rate):
    """Per-event weight exp(-rate * t / tau), t = DELIVERY time (t_min + downlink
    wait; v1.5.1 -- v1.5.0 used t_min, identical whenever every wait is 0, as in
    all runs so far), tau = decay_time_scale(store). An image never delivered has
    zero utility, so its weight does not matter (set to 0). Weights are floored
    at 1e-300: exp() underflows to exactly 0 beyond ~745 time scales (e.g. a
    16-day window on the fixed one-orbit scale), which would make a NEW tile
    unselectable under the facility objective."""
    tau = decay_time_scale(store)
    d = event_delivery_times(store)
    w = np.where(np.isfinite(d), np.exp(-latency_decay_rate * np.where(np.isfinite(d), d, 0.0) / tau), 0.0)
    return np.where(np.isfinite(d), np.maximum(w, 1e-300), 0.0)


def run_submodular_greedy_selection(
    events_dict, max_steps, cloud_mask=None, gap_fill=True,
    time_aware=False, latency_decay_rate=1.0,
):
    """
    Weighted maximum coverage, cardinality-constrained greedy, with
    Fail-Signal gap-filling and an optional time-aware objective.

    Implements the Nemhauser-Wolsey-Fisher (1978) greedy algorithm for
    submodular function maximisation. Selects the event with the
    highest marginal utility gain at each step, for exactly max_steps
    iterations, unless every remaining candidate can add zero new
    coverage (true information exhaustion). The (1-1/e) approximation
    guarantee applies directly to this unmodified cardinality-
    constrained weighted maximum coverage formulation.

    TIME-AWARE MODE (time_aware=True), each pass weighted
    w_e = exp(-latency_decay_rate * t_min / tau), tau =
    decay_time_scale(events_dict) (one orbital period, or the fleet's
    physical floor when that is longer -- see the section above):
      TIME_AWARE_OBJECTIVE "facility" (v1.5.0 default): greedy on the
        time-weighted coverage f(S) = sum_t max_{e in S} w_e * u_e,t, a
        weighted facility-location function -- monotone submodular, because
        the marginal gain of e, sum_t max(0, w_e u_e,t - max_{S} w u), can only
        shrink as S grows. The k-step greedy therefore attains at least
        (1 - 1/e) of the best k-pass value of f (Nemhauser, Wolsey & Fisher
        1978), and Cell 7C checks it against the exact ILP optimum of f.
        The greedy stops when no pass adds to f (a new tile always adds,
        since w_e > 0).
      "legacy" (v1.4.0): unweighted marginal gain times w_e. Those scores are
        not marginal gains of any one set function (two passes covering the
        same tile, weights 1 and 0.5, give the pair a value of 1 in one order
        and 0.5 in the other), so no guarantee applies.
    Reported coverage and quality (history, covered utility) are always the
    unweighted best utility per tile.

    FAIL-SIGNAL GAP-FILL: when cloud_mask is not None and gap_fill is
    True, a lazy-greedy set-cover search (Minoux 1978 acceleration --
    a max-heap of marginal new-tile counts, lazily re-validated on
    pop) fills any tiles left uncovered after the main k-step loop,
    always adding the event covering the most still-uncovered tiles,
    stopping when nothing remaining can cover anything outstanding.
    Never fires when cloud_mask is None (a clear-sky run has nothing
    to compensate for).

    IMPLEMENTATION: each step evaluates every event's marginal gain in one
    vectorised pass over the EventStore rows (gather covered utility,
    subtract, clip at zero, segment-sum). Ties go to the first event in
    key order, as in a sequential scan. Only per-tile state (covered
    utility, one float per ROI tile) is updated between steps.

    Args:
        events_dict: EventStore (Cell 6) -- event (sat, t_min) -> tiles, utilities
        max_steps: cardinality constraint k
        cloud_mask: optional CloudMask aligned with events_dict
        gap_fill: whether to run Fail-Signal gap-filling
        time_aware: bias selection toward earlier t_min
        latency_decay_rate: decay strength for time_aware mode

    Returns:
        (selected_events, covered_tiles_utility,
         cumulative_downlinked_tiles, history)
    """
    store = events_dict
    n_events = len(store)
    selected_events = []
    history = []
    cumulative_downlinked_tiles = 0
    covered = np.zeros(len(_tile_ids), dtype=np.float64)       # best utility per ROI tile
    if n_events == 0:
        return selected_events, {}, cumulative_downlinked_tiles, history

    tile_idx, starts = store.tile_idx, store.offsets[:-1]
    util = np.where(cloud_mask.flags, 0.0, store.util) if cloud_mask is not None else store.util
    selected = np.zeros(n_events, dtype=bool)
    weight = time_aware_weights(store, latency_decay_rate) if time_aware else None
    facility = time_aware and globals().get("TIME_AWARE_OBJECTIVE", "facility") == "facility"
    if facility:
        util_w = util * np.repeat(weight, store.counts)               # w_e * u_e,t per row
        covered_w = np.zeros(len(_tile_ids), dtype=np.float64)       # best weighted utility per tile
    buf = np.empty(store.n_rows, dtype=np.float64)

    for step in range(max_steps):
        if facility:
            np.take(covered_w, tile_idx, out=buf)
            np.subtract(util_w, buf, out=buf)
            np.maximum(buf, 0.0, out=buf)
            gains = np.add.reduceat(buf, starts)                       # marginal gain of f
            # A gain counts only if it exceeds 1e-12 in the pass's own utility
            # units (gain / w_e): rounding noise is ignored, while a new tile --
            # gain >= w_e * u -- always counts, however small w_e is.
            score = np.where(gains > 1e-12 * weight, gains, -np.inf)
            tol = 0.0
        else:
            np.take(covered, tile_idx, out=buf)
            np.subtract(util, buf, out=buf)
            np.maximum(buf, 0.0, out=buf)
            gains = np.add.reduceat(buf, starts)
            score = gains * weight if time_aware else gains
            tol = 1e-12
        score = np.where(selected, -np.inf, score)
        e = int(np.argmax(score))

        # Stop ONLY on exhaustion of the objective: coverage-first and legacy
        # time-aware stop when nothing adds real (unweighted) coverage; the
        # facility objective when nothing adds to f.
        if not np.isfinite(score[e]) or gains[e] <= tol:
            logger.info(f"  Stopped at step {step}: no remaining event adds "
                        f"to the objective (true ceiling for this dataset).")
            break

        selected[e] = True
        best_event = store.keys_[e]
        selected_events.append(best_event)
        cumulative_downlinked_tiles += int(store.counts[e])
        s = store.rows(e)
        np.maximum.at(covered, tile_idx[s], util[s])
        if facility:
            np.maximum.at(covered_w, tile_idx[s], util_w[s])

        cov_mask = covered > 0.0
        n_covered = int(cov_mask.sum())
        mean_quality = float(covered[cov_mask].sum() / n_covered) if n_covered else 0.0
        history.append({
            'step': step + 1,
            'selected': best_event,
            'coverage_pct': 100 * n_covered / total_target_tiles,
            'mean_quality': mean_quality,
        })

    # --- Fail-Signal gap-fill (lazy-greedy, heap-accelerated) ---
    if gap_fill and cloud_mask is not None and len(selected_events) > 0:
        in_access = np.zeros(len(_tile_ids), dtype=bool)
        in_access[tile_idx] = True
        remaining = in_access & ~(covered > 0.0)

        if remaining.any():
            n_uncovered_initial = int(remaining.sum())
            logger.info(f"  Gap-Filling: {n_uncovered_initial} tiles uncovered after greedy selection.")

            eligible = (store.util > 0.0) & ~cloud_mask.flags
            counts = np.add.reduceat((eligible & remaining[tile_idx]).astype(np.int64), starts)
            heap = [(-int(c), store.keys_[i]) for i, c in enumerate(counts) if c > 0 and not selected[i]]
            heapq.heapify(heap)

            gap_fill_count = 0
            while remaining.any() and heap:
                neg_cnt, event_key = heapq.heappop(heap)
                e = store.index[event_key]
                if selected[e]:
                    continue
                s = store.rows(e)
                true_cnt = int((eligible[s] & remaining[tile_idx[s]]).sum())
                if true_cnt == 0:
                    continue
                if heap and true_cnt < -heap[0][0]:
                    heapq.heappush(heap, (-true_cnt, event_key))
                    continue

                selected[e] = True
                selected_events.append(event_key)
                gap_fill_count += 1
                cumulative_downlinked_tiles += int(store.counts[e])
                t, u = tile_idx[s], util[s]
                improve = u > covered[t]
                covered[t[improve]] = u[improve]
                remaining[t[improve]] = False

            if remaining.any():
                logger.warning(f"  Gap-fill exhausted candidate events with "
                                f"{int(remaining.sum())} tiles still uncovered.")
            if gap_fill_count:
                logger.info(f"  Added {gap_fill_count} compensatory passes (Fail-Signal gap-fill).")

    covered_tiles_utility = {_tile_ids[i]: float(covered[i]) for i in np.flatnonzero(covered > 0.0)}
    return selected_events, covered_tiles_utility, cumulative_downlinked_tiles, history


# =============================================================================
# REALISTIC CLOUD KNOWLEDGE (non-clairvoyant selection under cloud)
# =============================================================================
# run_submodular_greedy_selection with a cloud mask chooses passes KNOWING
# which tiles each pass will see clear -- information no real system has
# before the image arrives. simulate_adaptive_selection is the causal
# version: plan on clear-sky geometry, reveal cloud when each image arrives,
# re-plan on failure (Fail-Signal). See CLOUD_KNOWLEDGE in Cell 0.
last_adaptive_run = {}


def event_delivery_times(store):
    """Delivery time (t_min + downlink wait, minutes) per event; +inf where
    no ground-station contact follows (such events have zero utility)."""
    return np.array([k[1] + (np.inf if store.wait_lookup.get(k) is None else store.wait_lookup[k])
                     for k in store.keys_], dtype=np.float64)


def _greedy_plan(store, state, candidates, budget, weight, util_plan=None):
    """Greedy utility-weighted max-coverage plan of up to `budget` events from
    `candidates` (bool per event), starting from per-tile best utility
    `state` (not modified). Uses CLEAR-SKY utilities only. Same selection
    rule and tie-breaking as run_submodular_greedy_selection. util_plan:
    per-row utilities to plan on (the time-weighted w * u of the facility
    objective, with `state` the best weighted utility per tile); default the
    clear-sky utilities, scored legacy-style (gain times weight)."""
    tile_idx, starts = store.tile_idx, store.offsets[:-1]
    util = store.util if util_plan is None else util_plan
    tol = 1e-12 if util_plan is None else 0.0
    state = state.copy()
    cand = candidates.copy()
    plan = []
    buf = np.empty(store.n_rows, dtype=np.float64)
    for _ in range(budget):
        np.take(state, tile_idx, out=buf)
        np.subtract(util, buf, out=buf)
        np.maximum(buf, 0.0, out=buf)
        gains = np.add.reduceat(buf, starts)
        if util_plan is None:
            score = gains * weight if weight is not None else gains
        else:                                      # facility objective: see run_submodular_greedy_selection
            score = np.where(gains > 1e-12 * weight, gains, -np.inf)
        score = np.where(cand, score, -np.inf)
        e = int(np.argmax(score))
        if not np.isfinite(score[e]) or gains[e] <= tol:
            break
        cand[e] = False
        plan.append(e)
        s = store.rows(e)
        np.maximum.at(state, tile_idx[s], util[s])
    return plan


def _expected_plan(store, uncovered, candidates, budget, p_tile, weight=None):
    """
    Greedy plan of up to `budget` events from `candidates` (bool per event)
    maximising the EXPECTED number of `uncovered` tiles (bool per ROI tile)
    that end up with a cloud-free usable view, each tile cloudy with
    probability p_tile in every pass independently (the planner's assumption;
    the simulated cloud may be correlated). The gain of a pass is the sum over
    its usable uncovered tiles of P(not yet covered by the plan) * (1 - p) --
    the greedy on sum_t [1 - prod_{e covering t} p_t], a monotone submodular
    function. weight: optional per-event multiplier (time-aware objective):
    the score is then gain times weight, which is not the marginal gain of any
    one set function, so with a weight this is a heuristic without the
    (1 - 1/e) bound (v1.5.2 note).
    """
    tile_idx, starts = store.tile_idx, store.offsets[:-1]
    usable = store.util > 0.0
    p = np.asarray(p_tile, dtype=np.float64)
    q = np.where(uncovered, 1.0, 0.0)
    one_minus_p = 1.0 - p
    cand = candidates.copy()
    plan = []
    for _ in range(int(budget)):
        g = np.where(usable, q[tile_idx] * one_minus_p[tile_idx], 0.0)
        raw = np.add.reduceat(g, starts)
        gains = raw * weight if weight is not None else raw
        gains = np.where(cand & (raw > 1e-12), gains, -np.inf)
        e = int(np.argmax(gains))
        if not np.isfinite(gains[e]):
            break
        cand[e] = False
        plan.append(e)
        s = store.rows(e)
        t = np.unique(tile_idx[s][usable[s]])
        q[t] *= p[t]
    return plan


def _planning_probabilities():
    pt = globals().get("CLOUD_P_TILE")
    return pt if pt is not None else np.full(len(_tile_ids), float(globals().get("CLOUD_PROBABILITY", 0.35)))


def simulate_adaptive_selection(events_dict, max_steps, cloud_mask, gap_fill=True,
                                time_aware=False, latency_decay_rate=1.0, planning=None, p_tile=None):
    """
    Causal (non-clairvoyant) selection under cloud. Same return values as
    run_submodular_greedy_selection. No decision uses the cloud flags of an
    image before it arrives.

    planning "expected" (v1.5.0 default, CLOUD_PLANNING in Cell 0), budget
    max_steps: plan the budget to maximise the EXPECTED coverage of the tiles
    still missing (_expected_plan, from each tile's cloud probability p_tile --
    climatology, not the realised cloud), time-weighted if time_aware. Execute
    in delivery order; each image's cloud is revealed when it arrives; after
    every arrival re-plan the remaining budget over images not yet arrived,
    from the coverage actually achieved. A planned new tile that comes back
    cloudy is a Fail-Signal (counted in `replans`).
    planning "clear" (v1.4.0): plan greedily on clear-sky geometry (the
    time-aware objective if time_aware) and re-plan only on a Fail-Signal (a
    tile a planned observation should have improved stays below its clear-sky
    utility). In clear sky there are no failures, so the selected set equals
    the clear-sky greedy plan.
    Phase 2 (gap_fill), both: after the budget, take each arriving image that
    clearly covers a still-uncovered accessible tile (decided on arrival).
    """
    store = events_dict
    n = len(store)
    if n == 0:
        last_adaptive_run.update(replans=0, compensatory=0)
        return [], {}, 0, []
    planning = planning or globals().get("CLOUD_PLANNING", "clear")
    tile_idx, util, flags = store.tile_idx, store.util, cloud_mask.flags
    realized = np.where(flags, 0.0, util)
    d = event_delivery_times(store)
    arrives = np.isfinite(d)
    weight = time_aware_weights(store, latency_decay_rate) if time_aware else None    # clear-sky geometry only
    facility = time_aware and globals().get("TIME_AWARE_OBJECTIVE", "facility") == "facility"
    covered = np.zeros(len(_tile_ids), dtype=np.float64)
    if facility:                       # plan and re-plan on f = sum_t max w u (Cell 0)
        _wr = np.repeat(weight, store.counts)
        util_plan, realized_w = util * _wr, realized * _wr
        covered_w = np.zeros(len(_tile_ids), dtype=np.float64)
    else:
        util_plan, realized_w, covered_w = None, None, covered
    taken = np.zeros(n, dtype=bool)
    selected, history = [], []
    downlinked, replans, n_comp = 0, 0, 0
    tau, budget = -np.inf, max_steps

    if planning == "expected":
        p_plan = np.asarray(p_tile if p_tile is not None else _planning_probabilities(), dtype=np.float64)
        usable = util > 0.0
        got = np.zeros(len(_tile_ids), dtype=bool)               # tiles with a cloud-free usable view so far
        plan = _expected_plan(store, ~got, arrives, budget, p_plan, weight)
        while budget > 0 and plan:
            i = min(range(len(plan)), key=lambda j: (d[plan[j]], j))   # next planned image to arrive
            e = plan.pop(i)
            tau = d[e]
            s = store.rows(e)
            t = tile_idx[s]
            replans += int(np.any(flags[s] & usable[s] & ~got[t]))    # a planned new tile came back cloudy
            np.maximum.at(covered, t, realized[s])
            got[t[realized[s] > 0.0]] = True
            taken[e] = True
            selected.append(store.keys_[e])
            downlinked += int(store.counts[e])
            budget -= 1
            cov_mask = covered > 0.0
            n_cov = int(cov_mask.sum())
            history.append({"step": len(history) + 1, "selected": store.keys_[e],
                            "coverage_pct": 100 * n_cov / total_target_tiles,
                            "mean_quality": float(covered[cov_mask].sum() / n_cov) if n_cov else 0.0})
            if budget > 0:
                plan = _expected_plan(store, ~got, (~taken) & arrives & (d >= tau), budget, p_plan, weight)

    plan = [] if planning == "expected" else _greedy_plan(store, covered_w, arrives, budget, weight, util_plan)
    while budget > 0 and plan:
        i = min(range(len(plan)), key=lambda j: (d[plan[j]], j))    # next planned image to arrive
        e = plan.pop(i)
        tau = d[e]
        s = store.rows(e)
        t = tile_idx[s]
        if facility:
            failed = bool(np.any(flags[s] & (util_plan[s] > covered_w[t] + 1e-12)))
            np.maximum.at(covered_w, t, realized_w[s])
        else:
            failed = bool(np.any(flags[s] & (util[s] > covered[t] + 1e-12)))
        np.maximum.at(covered, t, realized[s])
        taken[e] = True
        selected.append(store.keys_[e])
        downlinked += int(store.counts[e])
        budget -= 1
        cov_mask = covered > 0.0
        n_cov = int(cov_mask.sum())
        history.append({"step": len(history) + 1, "selected": store.keys_[e],
                        "coverage_pct": 100 * n_cov / total_target_tiles,
                        "mean_quality": float(covered[cov_mask].sum() / n_cov) if n_cov else 0.0})
        if failed and budget > 0:
            replans += 1
            plan = _greedy_plan(store, covered_w, (~taken) & arrives & (d >= tau), budget, weight, util_plan)

    if gap_fill and selected:
        accessible = np.zeros(len(_tile_ids), dtype=bool)
        accessible[tile_idx[util > 0.0]] = True
        remaining = accessible & ~(covered > 0.0)
        if remaining.any():
            for e in np.argsort(d, kind="stable"):
                if taken[e] or not arrives[e] or d[e] < tau:
                    continue
                s = store.rows(e)
                t, r = tile_idx[s], realized[s]
                hit = (r > 0.0) & remaining[t]
                if hit.any():
                    np.maximum.at(covered, t, r)
                    remaining[t[hit]] = False
                    taken[e] = True
                    selected.append(store.keys_[e])
                    downlinked += int(store.counts[e])
                    n_comp += 1
                    if not remaining.any():
                        break
            if remaining.any():
                logger.info(f"  Compensation ran out of arriving images with {int(remaining.sum())} tiles uncovered.")
    last_adaptive_run.update(replans=replans, compensatory=n_comp)
    covered_tiles_utility = {_tile_ids[i]: float(covered[i]) for i in np.flatnonzero(covered > 0.0)}
    return selected, covered_tiles_utility, downlinked, history


def select_passes(events_dict, max_steps, cloud_mask=None, gap_fill=True,
                  time_aware=False, latency_decay_rate=1.0):
    """The selector every cell uses. Clear sky (no mask) or CLOUD_KNOWLEDGE ==
    "oracle": run_submodular_greedy_selection unchanged. Otherwise the causal
    simulate_adaptive_selection."""
    if cloud_mask is None or globals().get("CLOUD_KNOWLEDGE", "realistic") == "oracle":
        return run_submodular_greedy_selection(events_dict, max_steps, cloud_mask=cloud_mask, gap_fill=gap_fill,
                                               time_aware=time_aware, latency_decay_rate=latency_decay_rate)
    return simulate_adaptive_selection(events_dict, max_steps, cloud_mask, gap_fill=gap_fill,
                                       time_aware=time_aware, latency_decay_rate=latency_decay_rate)


logger.info(f"=== Running clear-sky selection (k={MAX_SELECTION_STEPS}) ===")
sel_events_clear, cov_clear, down_clear, hist_clear = run_submodular_greedy_selection(
    event_data, MAX_SELECTION_STEPS, cloud_mask=None, gap_fill=True
)
logger.info(f"Clear-sky: {len(sel_events_clear)} events, "
            f"{hist_clear[-1]['coverage_pct']:.2f}% coverage" if hist_clear else "Clear-sky: no events selected")

logger.info(f"=== Running cloud/Fail-Signal selection (k={MAX_SELECTION_STEPS}, "
            f"cloud knowledge: {globals().get('CLOUD_KNOWLEDGE', 'realistic')}) ===")
sel_events_cloud, cov_cloud, down_cloud, hist_cloud = select_passes(
    event_data, MAX_SELECTION_STEPS, cloud_mask=cloud_mask, gap_fill=True
)
_adaptive_cloud_stats = dict(last_adaptive_run)
logger.info(f"Cloud: {len(sel_events_cloud)} events selected, "
            f"{100 * len(cov_cloud) / total_target_tiles:.2f}% final coverage"
            + (f" ({_adaptive_cloud_stats.get('replans', 0)} re-plans, "
               f"{_adaptive_cloud_stats.get('compensatory', 0)} compensatory passes)"
               if globals().get("CLOUD_KNOWLEDGE", "realistic") != "oracle" else ""))
# Clairvoyant reference, reported for transparency: the same greedy selector
# KNOWING the cloud (not a strict upper bound, since greedy is not optimal).
_oracle_sel_cloud, _oracle_cov_cloud, _, _ = run_submodular_greedy_selection(
    event_data, MAX_SELECTION_STEPS, cloud_mask=cloud_mask, gap_fill=True)


def compute_selection_latency(selected_events, refined_df):
    """
    SECONDARY METRIC -- last selected delivery: the delivery time
    (t_min + T_wait_min) of the LAST-arriving selected event. Reported
    alongside H1 but NOT the H1 metric: once coverage saturates, the greedy
    keeps adding events for quality gains, so this grows with the budget k
    even though the mosaic stopped changing. See compute_completion_latency.
    """
    if not selected_events:
        return None
    # EventStore keeps one wait per event; a DataFrame is still accepted.
    wait_lookup = (refined_df.wait_lookup if hasattr(refined_df, "wait_lookup") else
                   refined_df[["sat", "t_min", "T_wait_min"]]
                   .drop_duplicates(subset=["sat", "t_min"])
                   .set_index(["sat", "t_min"])["T_wait_min"]
                   .to_dict())
    waits = [wait_lookup.get((sat, t_min), 0) for sat, t_min in selected_events]
    delivery_times = [t_min + (np.nan if w is None else w) for (sat, t_min), w in zip(selected_events, waits)]
    return max(delivery_times)


def compute_completion_latency(selected_events, store, cloud_mask=None):
    """
    H1 METRIC -- mosaic completion time: with the selected events ordered by
    delivery time (t_min + T_wait_min), the delivery time at which the
    mosaic reaches its final coverage, i.e. of the last event that adds a
    tile not already delivered. A tile counts as covered under the same
    rule the selection uses: utility > 0 and not cloud-obscured. Events
    that only improve quality on already-covered tiles do not delay it.
    Returns None if no selected event covers anything.
    """
    if not selected_events:
        return None
    util = np.where(cloud_mask.flags, 0.0, store.util) if cloud_mask is not None else store.util
    def _delivery(key):
        w = store.wait_lookup.get(key)
        return key[1] + (np.nan if w is None else w)
    seen = np.zeros(len(_tile_ids), dtype=bool)
    completion = None
    for key in sorted(selected_events, key=delivery_sort_key(store)):   # stable: ties keep selection order
        s = store.rows(store.index[key])
        t = store.tile_idx[s][util[s] > 0.0]
        new = t[~seen[t]]
        if new.size:
            seen[new] = True
            completion = _delivery(key)
    return completion


latency_clear_min = compute_completion_latency(sel_events_clear, event_data)
latency_cloud_min = compute_completion_latency(sel_events_cloud, event_data, cloud_mask)
last_delivery_clear_min = compute_selection_latency(sel_events_clear, event_data)
last_delivery_cloud_min = compute_selection_latency(sel_events_cloud, event_data)
_fmt = lambda v: f"{v:.1f} min" if v is not None else "N/A (no events selected)"
logger.info(f"H1 latency (mosaic completion): clear={_fmt(latency_clear_min)}, cloud={_fmt(latency_cloud_min)} | "
            f"last selected delivery: clear={_fmt(last_delivery_clear_min)}, cloud={_fmt(last_delivery_cloud_min)}")

# Time-aware objective (decay LATENCY_DECAY_RATE), reported alongside the
# coverage-first objective above. Same budget, same cloud mask, same gap-fill.
LATENCY_DECAY_TAU_MIN = decay_time_scale(event_data)
logger.info(f"=== Time-aware selection (decay={LATENCY_DECAY_RATE}, tau={LATENCY_DECAY_TAU_MIN:.1f} min "
            f"[{globals().get('LATENCY_DECAY_SCALE', 'fleet')} scale; orbital period {float(T_MAX_WAIT_MINUTES):.1f} min], "
            f"k={MAX_SELECTION_STEPS}) ===")
sel_ta_clear, cov_ta_clear, down_ta_clear, hist_ta_clear = run_submodular_greedy_selection(
    event_data, MAX_SELECTION_STEPS, cloud_mask=None, gap_fill=True,
    time_aware=True, latency_decay_rate=LATENCY_DECAY_RATE)
sel_ta_cloud, cov_ta_cloud, down_ta_cloud, hist_ta_cloud = select_passes(
    event_data, MAX_SELECTION_STEPS, cloud_mask=cloud_mask, gap_fill=True,
    time_aware=True, latency_decay_rate=LATENCY_DECAY_RATE)
_oracle_sel_ta_cloud, _oracle_cov_ta_cloud, _, _ = run_submodular_greedy_selection(
    event_data, MAX_SELECTION_STEPS, cloud_mask=cloud_mask, gap_fill=True,
    time_aware=True, latency_decay_rate=LATENCY_DECAY_RATE)
ta_latency_clear_min = compute_completion_latency(sel_ta_clear, event_data)
ta_latency_cloud_min = compute_completion_latency(sel_ta_cloud, event_data, cloud_mask)
_mean_q = lambda cov: float(np.mean(list(cov.values()))) if cov else 0.0
logger.info(f"H1 latency, time-aware: clear={_fmt(ta_latency_clear_min)}, cloud={_fmt(ta_latency_cloud_min)} | "
            f"mean quality clear: coverage-first {_mean_q(cov_clear):.3f}, time-aware {_mean_q(cov_ta_clear):.3f}")

_cell6_summary = {
    "n_candidate_events": len(event_data),
    "selected_events_clear": len(sel_events_clear),
    "selected_events_cloud": len(sel_events_cloud),
    "coverage_clear_pct": hist_clear[-1]['coverage_pct'] if hist_clear else 0.0,
    "coverage_cloud_pct": 100 * len(cov_cloud) / total_target_tiles,
}
_out_dir = OUT_DIR if 'OUT_DIR' in globals() else os.getcwd()
with open(os.path.join(_out_dir, 'arf_cell6_summary.json'), 'w') as f:
    json.dump(_cell6_summary, f)


def _selection_table(selected_events, history):
    """One row per selected event, in selection order: greedy steps carry the
    coverage/quality after that step; gap-fill events follow with stage 'gap_fill'."""
    rows = [{"step": h["step"], "stage": "greedy", "sat": h["selected"][0], "t_min": h["selected"][1],
             "coverage_pct": h["coverage_pct"], "mean_quality": h["mean_quality"]} for h in history]
    for sat, t in selected_events[len(history):]:
        rows.append({"step": None, "stage": "gap_fill", "sat": sat, "t_min": t,
                     "coverage_pct": None, "mean_quality": None})
    return pd.DataFrame(rows, columns=["step", "stage", "sat", "t_min", "coverage_pct", "mean_quality"])


_selection_table(sel_events_clear, hist_clear).to_csv(
    os.path.join(_out_dir, "arf_selection_history_clear.csv"), index=False)
_selection_table(sel_events_cloud, hist_cloud).to_csv(
    os.path.join(_out_dir, "arf_selection_history_cloud.csv"), index=False)

# Record in scenario_results (the cross-scenario ledger from Cell 0)
# for the hypothesis validators (Cell 8) and any later comparison
# across runs. Keyed by (ACTIVE_ROI, SCENARIO_N): N for the Starlink proxy,
# the fleet tag for real fleets (v1.5.2: fleet, swath case, window, start and
# left-out satellites), so scenarios run in one kernel never share an entry.
def _scenario_key(active_roi=None, active_n=None):
    """Shared key for scenario_results, used consistently across Cell
    6, 6B, 7, and 8 so every cell keys and reads the same entry."""
    return (active_roi if active_roi is not None else ACTIVE_ROI,
            active_n if active_n is not None else globals().get("SCENARIO_N", ACTIVE_N))

if 'scenario_results' not in globals():
    scenario_results = {}
# v1.5.2: this run's entry starts empty. Cells 6B, 7C and 7D add to it; with
# setdefault, a re-run of the same scenario kept fields of the previous run that
# the new run did not overwrite (e.g. its onset analysis, if 7C was skipped).
scenario_results[_scenario_key()] = {}
_n_covered_clear = len(cov_clear)
_redundancy_clear = (down_clear / _n_covered_clear) if _n_covered_clear else float('nan')
_n_covered_cloud = len(cov_cloud)
_redundancy_cloud = (down_cloud / _n_covered_cloud) if _n_covered_cloud else float('nan')

# H3 redundancy is measured AT THE POINT THE COVERAGE REQUIREMENT IS MET, not
# over all k passes: once coverage is complete the greedy keeps selecting passes
# to refine image quality, so all-k redundancy grows with the budget (e.g. 30.0
# when each of 30 passes covers a whole small ROI) whatever the efficiency of
# reaching coverage. _H3_COVERAGE_PCT matches Cell 8's H3_COVERAGE_THRESHOLD_PCT
# (Cell 8 runs later, so it is restated here and checked there).
_H3_COVERAGE_PCT = 95.0


def redundancy_at_coverage(selected, store, mask, level_pct):
    """
    Tile observations downlinked by the selection, in selection order, up to
    and including the pass at which coverage first reaches level_pct, divided
    by the unique tiles covered at that point. Every observation of a pass
    counts, cloud-obscured ones included (they are downlinked too); a tile
    counts as covered under the usual rule (utility > 0, not obscured).
    Returns (redundancy or None if the level is never reached, passes used).
    """
    util = np.where(mask.flags, 0.0, store.util) if mask is not None else store.util
    seen = np.zeros(len(_tile_ids), dtype=bool)
    n_seen = n_obs = 0
    for i, key in enumerate(selected, 1):
        s = store.rows(store.index[key])
        t = store.tile_idx[s]
        n_obs += t.size
        new = t[(util[s] > 0.0) & ~seen[t]]
        seen[new] = True
        n_seen += new.size
        if n_seen and 100.0 * n_seen / total_target_tiles >= level_pct - 1e-9:
            return n_obs / n_seen, i
    return None, len(selected)


_h3 = {
    "h3_redundancy_clear": redundancy_at_coverage(sel_events_clear, event_data, None, _H3_COVERAGE_PCT),
    "h3_redundancy_cloud": redundancy_at_coverage(sel_events_cloud, event_data, cloud_mask, _H3_COVERAGE_PCT),
    "ta_h3_redundancy_clear": redundancy_at_coverage(sel_ta_clear, event_data, None, _H3_COVERAGE_PCT),
    "ta_h3_redundancy_cloud": redundancy_at_coverage(sel_ta_cloud, event_data, cloud_mask, _H3_COVERAGE_PCT),
}
# Under the independent per-observation cloud model, a tile needs n looks with
# 1 - p^n >= level to be seen clearly with that probability; averaged over
# tiles this bounds the expected redundancy at that coverage from below.
# (v1.5.1) Valid for the independent ("iid") cloud model with a constant
# probability only; correlated cloud needs more looks, so no floor is reported.
_cloud_redundancy_floor = (float(np.log(1.0 - _H3_COVERAGE_PCT / 100.0) / np.log(CLOUD_PROBABILITY))
                           if 0.0 < CLOUD_PROBABILITY < 1.0 and _H3_COVERAGE_PCT < 100.0
                           and globals().get("CLOUD_MODEL", "iid") == "iid"
                           and globals().get("CLOUD_PROBABILITY_SOURCE", "constant") == "constant" else None)
logger.info(f"H3 redundancy at {_H3_COVERAGE_PCT:g}% coverage: clear "
            f"{_h3['h3_redundancy_clear'][0]} ({_h3['h3_redundancy_clear'][1]} passes), cloud "
            f"{_h3['h3_redundancy_cloud'][0]} ({_h3['h3_redundancy_cloud'][1]} passes; "
            f"model floor ~{_cloud_redundancy_floor})")
scenario_results[_scenario_key()].update({
    'candidate_events': len(event_data),
    'selection_budget_k': MAX_SELECTION_STEPS,
    'coverage_clear_pct': hist_clear[-1]['coverage_pct'] if hist_clear else 0.0,
    'coverage_cloud_pct': 100 * _n_covered_cloud / total_target_tiles,
    'redundancy_clear': _redundancy_clear,           # over all selected passes (budget-driven)
    'redundancy_cloud': _redundancy_cloud,
    'h3_coverage_pct': _H3_COVERAGE_PCT,              # H3: redundancy when this coverage is first reached
    'h3_redundancy_clear': _h3['h3_redundancy_clear'][0],
    'h3_passes_clear': _h3['h3_redundancy_clear'][1],
    'h3_redundancy_cloud': _h3['h3_redundancy_cloud'][0],
    'h3_passes_cloud': _h3['h3_redundancy_cloud'][1],
    'ta_h3_redundancy_clear': _h3['ta_h3_redundancy_clear'][0],
    'ta_h3_redundancy_cloud': _h3['ta_h3_redundancy_cloud'][0],
    'cloud_redundancy_floor': _cloud_redundancy_floor,
    'selected_events_clear': len(sel_events_clear),
    'selected_events_cloud': len(sel_events_cloud),
    'reachable_tiles_pct': 100 * int(np.unique(event_data.tile_idx).size) / total_target_tiles,
    # Tiles with at least one USABLE clear-sky observation (utility > 0): the
    # most any selection could cover. H3 is not evaluable below its threshold.
    'reachable_usable_clear_pct': 100 * int(np.unique(event_data.tile_idx[event_data.util > 0.0]).size)
                                  / total_target_tiles,
    'gap_fill_enabled': True,
    'latency_clear_min': latency_clear_min,          # H1: mosaic completion time
    'latency_cloud_min': latency_cloud_min,
    'last_delivery_clear_min': last_delivery_clear_min,   # secondary
    'last_delivery_cloud_min': last_delivery_cloud_min,
    'mean_quality_clear': _mean_q(cov_clear),
    'mean_quality_cloud': _mean_q(cov_cloud),
    'latency_decay_rate': LATENCY_DECAY_RATE,
    'latency_decay_scale': globals().get("LATENCY_DECAY_SCALE", "fleet"),
    'latency_decay_tau_min': LATENCY_DECAY_TAU_MIN,   # time scale of the time-aware weight (full window)
    'orbital_period_min': float(T_MAX_WAIT_MINUTES),
    'ta_latency_clear_min': ta_latency_clear_min,     # H1, time-aware objective
    'ta_latency_cloud_min': ta_latency_cloud_min,
    'ta_coverage_clear_pct': 100 * len(cov_ta_clear) / total_target_tiles,
    'ta_coverage_cloud_pct': 100 * len(cov_ta_cloud) / total_target_tiles,
    'ta_mean_quality_clear': _mean_q(cov_ta_clear),
    'ta_mean_quality_cloud': _mean_q(cov_ta_cloud),
    'cloud_knowledge': globals().get("CLOUD_KNOWLEDGE", "realistic"),
    'time_aware_objective': globals().get("TIME_AWARE_OBJECTIVE", "facility"),   # v1.5.0
    'cloud_planning': globals().get("CLOUD_PLANNING", "clear"),                 # v1.5.0
    'cloud_model': globals().get("cloud_model_info"),                           # v1.5.0
    # v1.5.0: True when no selected-or-candidate image waits for a downlink
    # (every wait 0), so "delivery time" is the acquisition time.
    'acquisition_time_only': bool(np.all(np.nan_to_num(np.array(list(event_data.wait_lookup.values()),
                                                                 dtype=np.float64), nan=1.0) == 0.0))
                             if len(event_data) else None,
    'cloud_replans': _adaptive_cloud_stats.get("replans"),
    'cloud_compensatory_passes': _adaptive_cloud_stats.get("compensatory"),
    # v1.5.1: the window start; scenarios are comparable only if it is the same
    'reference_epoch': globals().get("_reference_epoch"),
    'reference_epoch_source': globals().get("REFERENCE_EPOCH_SOURCE"),
    # Reference: the same greedy selectors given every image's cloud in advance
    # (clairvoyant greedy). Not a strict upper bound: greedy is not optimal (v1.5.0 note).
    'oracle_coverage_cloud_pct': 100 * len(_oracle_cov_cloud) / total_target_tiles,
    'oracle_latency_cloud_min': compute_completion_latency(_oracle_sel_cloud, event_data, cloud_mask),
    'oracle_ta_coverage_cloud_pct': 100 * len(_oracle_cov_ta_cloud) / total_target_tiles,
    'oracle_ta_latency_cloud_min': compute_completion_latency(_oracle_sel_ta_cloud, event_data, cloud_mask),
})


def build_sweep_results(event_data, refined_df, cloud_mask, k_values=None):
    """
    Sweep coverage and latency across a range of selection budgets k,
    independent of MAX_SELECTION_STEPS -- feeds Cell 7's Figure 3
    (Reconstruction Latency Analysis) and Figure 4 (Coverage vs.
    Budget Trade-off), which show how these metrics move as the
    budget itself changes, not just their value at the one k Cell 6
    ran with.

    Uses gap_fill=False so the sweep reflects the raw greedy
    coverage/latency tradeoff at each k, not values inflated by the
    separate Fail-Signal gap-fill step (which exists to guarantee
    full coverage under cloud, not to characterise the tradeoff
    curve).
    """
    if k_values is None:
        _ceiling = len(event_data)
        k_values = sorted(k for k in {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70} if k <= _ceiling)

    results = []
    if not k_values:
        return results
    # Greedy selection is incremental: the events chosen for budget k are
    # exactly the first k of any longer run (and a run that exhausts early
    # stops at the same step either way). One run at the largest k therefore
    # yields every smaller budget's selection and coverage exactly.
    k_max = max(k_values)
    _, _, _, hist_c_full = run_submodular_greedy_selection(
        event_data, max_steps=k_max, cloud_mask=None, gap_fill=False)
    _realistic = globals().get("CLOUD_KNOWLEDGE", "realistic") != "oracle"
    if not _realistic:
        _, _, _, hist_cl_full = run_submodular_greedy_selection(
            event_data, max_steps=k_max, cloud_mask=cloud_mask, gap_fill=False)

    for k in k_values:
        hist_c = hist_c_full[:k]
        # Realistic cloud: re-planning depends on the budget, so a smaller k is
        # NOT a prefix of a larger one -- each k is simulated on its own.
        hist_cl = (select_passes(event_data, k, cloud_mask=cloud_mask, gap_fill=False)[3]
                   if _realistic else hist_cl_full[:k])
        sel_c = [h['selected'] for h in hist_c]
        cov_c_pct = hist_c[-1]['coverage_pct'] if hist_c else 0.0
        cov_cl_pct = hist_cl[-1]['coverage_pct'] if hist_cl else 0.0
        latency_c = compute_completion_latency(sel_c, event_data)
        last_c = compute_selection_latency(sel_c, event_data)

        results.append({
            'k': k,
            'coverage_clear_pct': cov_c_pct,
            'coverage_cloud_pct': cov_cl_pct,
            'latency_clear_min': latency_c if latency_c is not None else 0,       # H1 metric
            'last_delivery_clear_min': last_c if last_c is not None else 0,       # secondary
        })
    return results


sweep_results = build_sweep_results(event_data, refined_events_df, cloud_mask)
pd.DataFrame(sweep_results).to_csv(os.path.join(_out_dir, "arf_coverage_sweep.csv"), index=False)
logger.info(f"Built sweep_results: {len(sweep_results)} k-values "
            f"({sweep_results[0]['k']}-{sweep_results[-1]['k']}) for Cell 7's Figures 3/4."
            if sweep_results else "Built sweep_results: empty (no valid k values for this event_data size).")

print("\n=== CELL 6 COMPLETE ===")

