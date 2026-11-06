# =============================================================================
# CELL 7B: COVERAGE-VS-TIME CURVE (honest check on H1 before writing it off)
# Algorithmic Reconstruction Framework (ARF)
# Run after Cell 6 (needs sel_events_clear, sel_events_cloud, event_data,
# refined_events_df, total_target_tiles already in the kernel namespace).
# =============================================================================
#
# H1 (Cell 6/8) is the mosaic COMPLETION time: when delivered coverage
# reaches its final value. This cell adds the partial milestones (time to
# 25/50/75/90% coverage), useful when a region is usable before it is
# complete. Its final-coverage time equals the H1 figure, because tiles
# are counted under the same rule (utility > 0, not cloud-obscured).
#
# This reorders the SELECTED events by actual delivery time (t_min +
# T_wait_min) — NOT greedy selection order, which has no temporal meaning
# — and tracks cumulative coverage against real elapsed time. Verified
# against a deliberately scrambled test case before being run for real:
# correctly re-sorts regardless of the order events were greedily chosen in.
# Output (v1.5.2): arf_coverage_vs_time.csv (both curves); images that are
# never delivered (no downlink in the window) are left out.

# Version guard (v1.5.0; every cell from v1.5.1): stop if this cell and Cell 0 come
# from different ARF releases -- a mixed notebook would silently run old code under a new label.
_ARF_CELL_VERSION = "1.5.3"
if globals().get("__version__") != _ARF_CELL_VERSION:
    raise RuntimeError(f"Cell 7B is from ARF {_ARF_CELL_VERSION}, but Cell 0 reports ARF "
                       f"{globals().get('__version__', 'unknown')}. Run every cell from one release.")
globals().setdefault("_ARF_CELLS_RUN", {})["7B"] = _ARF_CELL_VERSION

import numpy as np


def compute_coverage_vs_time(selected_events, event_data, refined_events_df,
                               total_target_tiles, cloud_mask=None):
    # Build (sat, t_min) -> T_wait_min lookup from the refined dataset
    # (T_wait_min is constant per event, computed once in Cell 5)
    wait_lookup = (event_data.wait_lookup if hasattr(event_data, "wait_lookup") else
                   refined_events_df[["sat", "t_min", "T_wait_min"]]
                   .drop_duplicates(subset=["sat", "t_min"])
                   .set_index(["sat", "t_min"])["T_wait_min"]
                   .to_dict())

    # v1.5.2: an image with no downlink in the window (wait None or NaN) is
    # never delivered, so it is left out of the curve; before, None crashed
    # the addition and NaN left the sort unordered (Cell 6, delivery_sort_key).
    delivery_times = []
    for event_key in selected_events:
        sat, t_min = event_key
        t_wait = wait_lookup.get((sat, t_min), 0)
        if t_wait is None or (isinstance(t_wait, float) and np.isnan(t_wait)):
            continue
        delivery_times.append((t_min + t_wait, event_key))
    delivery_times.sort(key=lambda x: x[0])       # stable: ties keep selection order

    covered_tiles = set()
    curve = []
    for delivery_time, event_key in delivery_times:
        for tile, util in event_data[event_key].items():
            if util <= 0.0:                 # same coverage rule as the selection (Cell 6)
                continue
            if cloud_mask is not None and cloud_mask.get((event_key, tile), False):
                continue
            covered_tiles.add(tile)
        coverage_pct = len(covered_tiles) / total_target_tiles * 100
        curve.append({"delivery_time_min": delivery_time, "coverage_pct": coverage_pct,
                      "event": event_key})

    return curve


def time_to_reach(curve, target_pct):
    for point in curve:
        if point["coverage_pct"] >= target_pct:
            return point["delivery_time_min"]
    return None  # never reached


print("=== COVERAGE-VS-TIME (Clear-Sky) ===")
curve_clear = compute_coverage_vs_time(sel_events_clear, event_data, refined_events_df, total_target_tiles)
final_cov_clear = curve_clear[-1]["coverage_pct"] if curve_clear else 0.0

print(f"{'delivery_time_min':<20}{'coverage_pct':<15}")
for pt in curve_clear:
    print(f"{pt['delivery_time_min']:<20}{pt['coverage_pct']:.2f}")

print(f"\nFinal coverage: {final_cov_clear:.2f}%")
for pct in [25, 50, 75, 90]:
    if pct <= final_cov_clear:
        t = time_to_reach(curve_clear, pct)
        print(f"Time to reach {pct}% coverage: {t} min ({t/60:.1f} hours)")

print("\n=== COVERAGE-VS-TIME (Cloud/Fail-Signal) ===")
curve_cloud = compute_coverage_vs_time(sel_events_cloud, event_data, refined_events_df,
                                         total_target_tiles, cloud_mask=cloud_mask)
final_cov_cloud = curve_cloud[-1]["coverage_pct"] if curve_cloud else 0.0

print(f"{'delivery_time_min':<20}{'coverage_pct':<15}")
for pt in curve_cloud:
    print(f"{pt['delivery_time_min']:<20}{pt['coverage_pct']:.2f}")

print(f"\nFinal coverage: {final_cov_cloud:.2f}%")
for pct in [25, 50, 75, 90]:
    if pct <= final_cov_cloud:
        t = time_to_reach(curve_cloud, pct)
        print(f"Time to reach {pct}% coverage: {t} min ({t/60:.1f} hours)")

# v1.5.2: the two curves as a file (they were only printed), one row per
# delivered image of the Cell 6 selections (coverage-first, gap-fill included).
import os as _os7b
import pandas as _pd7b
_pd7b.DataFrame(
    [{"sky": sky, "delivery_time_min": p["delivery_time_min"], "coverage_pct": p["coverage_pct"],
      "sat": p["event"][0], "t_min": p["event"][1]}
     for sky, curve in (("clear", curve_clear), ("cloud", curve_cloud)) for p in curve],
    columns=["sky", "delivery_time_min", "coverage_pct", "sat", "t_min"]).to_csv(
    _os7b.path.join(OUT_DIR if 'OUT_DIR' in globals() else _os7b.getcwd(), "arf_coverage_vs_time.csv"), index=False)
print("Saved arf_coverage_vs_time.csv")

print("\n=== CELL 7B COMPLETE ===")

