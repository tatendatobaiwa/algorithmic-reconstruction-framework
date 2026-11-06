"""Cell 5 (v1.5.2): sourced acquisition exclusions applied per pass time, and
the end-cap count. Usage: test_cell5_exclusions.py PIPELINE_v152
A. Synthetic fleet over the real Lower Zambezi tiles (e2e15 harness): an
   exclusion window for one satellite removes exactly its passes inside the
   window (strip minute in [start, end)), nothing else; the summary counts
   the removed passes and rows.
B. _acquisition_exclusion_windows: name matching ignores the ' [NORAD id]'
   suffix Cell 2 adds to duplicate names; windows in minutes from the start;
   the Starlink proxy gets none.
C. The real Cell 0 table has no effect on a window outside it (all rows kept)."""
import sys, os, json, ast
import numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/validate")
import e2e15

P = sys.argv[1]
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


tiles = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
sats = e2e15.eo_fleet(seed=11, n_dove=30)
ov = {"window_hours": 24, "onset_shift_hours": 12}


def run(excl):
    ns = e2e15.run(P, tiles, cells=("4", "5"), sats=sats, overrides=ov, extra={"ACQUISITION_EXCLUSIONS": excl})
    df = ns["refined_events_df"].copy()
    df["sat"], df["tile"] = df["sat"].astype(str), df["tile"].astype(str)
    return df, json.load(open(os.path.join(ns["_out_dir"], "arf_cell5_summary.json"))), ns


base, s0, _ = run([])
# choose the satellite with the most passes and a window covering about half of them
per = base.groupby("sat").t_min.unique()
sat = max(per.index, key=lambda s: len(per[s]))
mins = np.sort(per[sat])
lo_min, hi_min = int(mins[0]), int(mins[len(mins) // 2]) + 1          # [first pass, middle pass] inclusive
t0 = e2e15.T0
win = [{"name_pattern": rf"^{sat}$", "start_utc": str(t0 + np.timedelta64(lo_min * 60, "s")),
        "end_utc": str(t0 + np.timedelta64(hi_min * 60, "s")), "status": "test", "source_url": "https://example.org",
        "quote": "test window"}]
cut, s1, _ = run(win)
inside = (base.sat == sat) & (base.t_min >= lo_min) & (base.t_min < hi_min)
expect = base[~inside].reset_index(drop=True)
check(cut.reset_index(drop=True).equals(expect),
      f"exclusion window for {sat} over minutes [{lo_min}, {hi_min}) removes exactly its {int(inside.sum())} rows "
      f"in {base[inside].t_min.nunique()} passes; all other rows identical")
check(s1["rows_acquisition_excluded"] == int(inside.sum()) and
      s1["passes_acquisition_excluded"] == int(base[inside].t_min.nunique()) and
      s1["acquisition_exclusions_hit"] and s1["acquisition_exclusions_hit"][0]["rows"] == int(inside.sum()),
      f"summary counts: {s1['passes_acquisition_excluded']} passes, {s1['rows_acquisition_excluded']} rows")
check(s0["rows_acquisition_excluded"] == 0 and s0["passes_acquisition_excluded"] == 0,
      f"no exclusions: nothing removed (end-cap rows not counted: {s0['rows_end_cap_excluded']})")

# B: the window helper on its own
C = e2e15.cells_of(open(P, encoding="utf-8").read())
tree = ast.parse(e2e15.strip_org(C["5"]))
tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_acquisition_exclusion_windows"]


class _T:
    def toString(self, scale=None):
        return "2026-01-10T00:00:00.000Z"


ns = {"np": np, "_sim_epoch_datetime64": lambda: np.datetime64("2026-01-10T00:00:00"), "FLEET": "earth-resources",
      "_sat_meta": {"SENTINEL-2A": {}, "SENTINEL-2A [99999]": {}, "SENTINEL-2B": {}, "LANDSAT 9": {}}}
exec(compile(tree, "c5win", "exec"), ns)
src0 = e2e15.cells_of(open(P, encoding="utf-8").read())["0"]
_cfg = {}
exec(compile(ast.Module(body=[n for n in ast.parse(src0).body if isinstance(n, ast.Assign) and
                              any(getattr(t, "id", "") == "ACQUISITION_EXCLUSIONS" for t in n.targets)],
                        type_ignores=[]), "c0ex", "exec"), _cfg)
ns["ACQUISITION_EXCLUSIONS"] = _cfg["ACQUISITION_EXCLUSIONS"]
w = ns["_acquisition_exclusion_windows"]()
_t0 = np.datetime64("2026-01-10T00:00:00")
exp_a = float((np.datetime64("2026-01-13T00:00:00") - _t0) / np.timedelta64(1, "m"))
exp_b = float((np.datetime64("2026-03-09T00:00:00") - _t0) / np.timedelta64(1, "m"))
check(set(w) == {"SENTINEL-2A", "SENTINEL-2A [99999]"} and all(abs(a - exp_a) < 1e-9 and abs(b - exp_b) < 1e-9
                                                             for v in w.values() for a, b, _ in v),
      f"Sentinel-2A tandem pause: matched with and without the NORAD suffix, 2B/Landsat not; window "
      f"[{exp_a:.0f}, {exp_b:.0f}) min from 10 Jan 2026")
ns["FLEET"] = "starlink"
check(ns["_acquisition_exclusion_windows"]() == {}, "Starlink proxy: no exclusion windows")

# C: the release table on an October window (e2e15 start 2026-10-07) changes nothing
rel, s2, _ = run(ns["ACQUISITION_EXCLUSIONS"])
check(rel.reset_index(drop=True).equals(base.reset_index(drop=True)) and s2["rows_acquisition_excluded"] == 0,
      "release exclusion table on a window outside it: every row kept")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
os.makedirs("/home/claude/validate/logs", exist_ok=True)
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open("/home/claude/validate/logs/test_cell5_exclusions.json", "w"),
          indent=1)
sys.exit(1 if n_fail else 0)
