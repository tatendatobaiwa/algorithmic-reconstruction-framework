"""Replay of Cells 6-8 of a pipeline file (v1.5.x) on the REAL refined events
of a 12 Oct 2026 Colab run. Cell 0 is executed for real (with
ARF_CONFIG_OVERRIDES), so every configuration global is the pipeline's own;
the outputs of Cells 1-5 that Cells 6-8 read are injected from the run's files:
  <run>_refined_events.csv  -> refined_events_df (Cell 5's output, unchanged)
  ROI tile list             -> target_tiles, _tile_ids (Cell 4's latitude order)
  Cell 2 summary values     -> CONSTELLATION_ALTITUDE_KM, T_MAX_WAIT_MINUTES,
                               sim_epoch, constellation size
Options: dedupe_caps (keep, per satellite and tile, only the strip with the
smallest off-nadir angle among consecutive minutes: the end-cap correction of
Cell 5, applied to existing data), family (subset the events to one fleet
family on the same window), cells (which cells to run).
Data files are only read; nothing in them is executed."""
import os, io, re, sys, json, math, tempfile, contextlib, logging, typing, time
import numpy as np, pandas as pd, h3, matplotlib
matplotlib.use("Agg")

RD = "/home/claude/realdata/"
RE_KM, MU = 6378.137, 398600.4418
# run: (refined csv, ROI, fleet, planet case, n satellites, median altitude km, sim epoch UTC)
RUNS = {
    "ckgr_open":   ("ckgr_open",   "ckgr_okavango", "earth-resources", "conservative", 5, 780.9448039106837, "2026-10-12T03:36:03.604608"),
    "ckgr_openPO": ("ckgr_openPO", "ckgr_okavango", "earth-resources+planet", "optimistic", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "ckgr_openPC": ("ckgr_openPC", "ckgr_okavango", "earth-resources+planet", "conservative", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "ckgr_PO":     ("ckgr_PO",     "ckgr_okavango", "planet", "optimistic", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
    "ckgr_PC":     ("ckgr_PC",     "ckgr_okavango", "planet", "conservative", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
    "lz_open":     ("lz_open",     "lower_zambezi", "earth-resources", "conservative", 5, 782.461934128808, "2026-10-12T08:27:09.082080"),
    "lz_openPO":   ("lz_openPO",   "lower_zambezi", "earth-resources+planet", "optimistic", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "lz_openPC":   ("lz_openPC",   "lower_zambezi", "earth-resources+planet", "conservative", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "lz_PO":       ("lz_PO",       "lower_zambezi", "planet", "optimistic", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
    "lz_PC":       ("lz_PC",       "lower_zambezi", "planet", "conservative", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
    "karoo_open":  ("karoo_open",  "karoo_nature_reserve", "earth-resources", "conservative", 5, 782.461934128808, "2026-10-12T08:27:09.082080"),
    "karoo_openPO": ("karoo_openPO", "karoo_nature_reserve", "earth-resources+planet", "optimistic", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "karoo_openPC": ("karoo_openPC", "karoo_nature_reserve", "earth-resources+planet", "conservative", 97, 473.3345109664915, "2026-10-12T08:27:09.082080"),
    "karoo_PO":    ("karoo_PO",    "karoo_nature_reserve", "planet", "optimistic", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
    "karoo_PC":    ("karoo_PC",    "karoo_nature_reserve", "planet", "conservative", 92, 475.44764351422737, "2026-10-12T07:59:24.084096"),
}
TILES = {"ckgr_okavango": RD + "ckgr_openPO_mosaic_tiles.csv", "lower_zambezi": RD + "lz_tiles.csv",
         "karoo_nature_reserve": RD + "karoo_tiles.csv"}
FAMILY_ALT = {"open": 780.9448039106837, "planet": 475.44764351422737}   # median altitude of each family alone


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


def dedupe_caps(df):
    """Keep, per (satellite, tile), one row per pass: among rows in consecutive
    minutes, the one with the smallest off-nadir angle (the strip containing
    the closest point of approach). Approximates Cell 5's end-cap fix on data
    already produced; rows of different passes (minutes apart) are all kept."""
    d = df.sort_values(["sat", "tile", "t_min"], kind="stable").reset_index()
    new_pass = ((d.sat != d.sat.shift()) | (d.tile != d.tile.shift()) | (d.t_min - d.t_min.shift() > 1)).cumsum()
    keep = d.groupby(new_pass)["A_nadir_deg"].idxmin()
    return df.loc[np.sort(d.loc[keep, "index"].to_numpy())]


class _Epoch:
    def __init__(self, iso): self.iso = iso
    def toString(self, scale=None): return self.iso + "Z"
    def __str__(self): return self.iso + "Z"


def run(pipe_path, run_name, overrides=None, cells=("6", "6B", "7", "7B", "7C", "8"), dedupe=False, family=None,
        quiet=True, extra=None):
    csv_name, roi, fleet, case, n_sats, alt, epoch = RUNS[run_name]
    src = open(pipe_path, encoding="utf-8").read()
    C = cells_of(src)
    refined = pd.read_csv(RD + f"{csv_name}_refined_events.csv")
    if family:
        is_open = refined.sat.str.startswith(("SENTINEL", "LANDSAT"))
        refined = refined[is_open] if family == "open" else refined[~is_open]
        alt = FAMILY_ALT[family]
        fleet = "earth-resources" if family == "open" else "planet"
        n_sats = 5 if family == "open" else 92
    if dedupe:
        refined = dedupe_caps(refined)
    refined = refined.reset_index(drop=True)
    tiles = set(pd.read_csv(TILES[roi], usecols=["tile"])["tile"])
    t_orbit = 2 * np.pi * np.sqrt((RE_KM + alt) ** 3 / MU) / 60.0
    out_dir, work_dir = tempfile.mkdtemp(prefix="R150_OUT_"), tempfile.mkdtemp(prefix="R150_CWD_")
    ov = {"active_roi": roi, "fleet": fleet, "planet_swath_case": case}
    ov.update(overrides or {})
    ns = {"__name__": "__main__", "ARF_CONFIG_OVERRIDES": ov, "display": lambda *a, **k: None}
    sink = io.StringIO()
    old = os.getcwd()
    try:
        with contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext():
            os.chdir(out_dir)
            exec(compile(C["0"], "cell0", "exec"), ns)          # OUT_DIR = out_dir (local run)
            os.chdir(work_dir)                                   # later stray writes land here, not in OUT_DIR
            ll = np.array([h3.cell_to_latlng(t) for t in sorted(tiles)])
            ids = np.array(sorted(tiles))
            order = np.argsort(ll[:, 0], kind="stable")          # Cell 4's latitude order
            names = sorted(set(refined.sat))
            fam = lambda s: ("Sentinel-2" if s.startswith("SENTINEL") else "Landsat 8/9" if s.startswith("LANDSAT")
                             else "PlanetScope")
            constellation = [{"id": s, "mission": fam(s)} for s in names]
            constellation += [{"id": f"(no pass) {i}", "mission": "?"} for i in range(max(0, n_sats - len(names)))]
            ns.update(
                np=np, pd=pd, h3=h3, logging=logging, os=os, json=json, math=math,
                **{n: getattr(typing, n) for n in ("List", "Dict", "Tuple", "Any", "Optional", "Set")},
                Re_km=RE_KM, mu_earth=MU, CONSTELLATION_ALTITUDE_KM=alt, T_MAX_WAIT_MINUTES=t_orbit,
                target_tiles=tiles, _tile_ids=ids[order], _tile_lat=ll[order, 0], _tile_lon=ll[order, 1],
                refined_events_df=refined, constellation=constellation, sim_epoch=_Epoch(epoch),
                _sim_epoch_datetime64=lambda: np.datetime64(epoch), GROUND_STATIONS={},
                ROI_DEFINITIONS={roi: {"short_name": roi}}, roi_boundary_report=None,
                _ARF_CELLS_RUN={**ns.get("_ARF_CELLS_RUN", {}), **{c: ns["__version__"] for c in ("1", "2", "3", "4", "5")}},
            )
            ns.update(extra or {})
            t0 = time.time()
            for c in cells:
                exec(compile(C[c], f"cell{c}", "exec"), ns)
            ns["_seconds"] = time.time() - t0
    finally:
        os.chdir(old)
    ns["_stdout"] = sink.getvalue()
    ns["_out_dir"], ns["_work_dir"] = out_dir, work_dir
    return ns


def result(ns):
    return ns["scenario_results"][ns["_scenario_key"]()]


if __name__ == "__main__":
    p = sys.argv[1]
    name = sys.argv[2]
    ns = run(p, name, quiet=True)
    d = result(ns)
    print(json.dumps({k: d.get(k) for k in ("coverage_clear_pct", "coverage_cloud_pct", "ta_coverage_clear_pct",
                                            "ta_coverage_cloud_pct", "latency_decay_tau_min", "reachable_usable_clear_pct")},
                     default=str))
    print("seconds", round(ns["_seconds"], 1), "out", ns["_out_dir"], "stray files:", os.listdir(ns["_work_dir"]))
