"""Real-data check of the v1.5.1 end-cap rule (Cell 5) on the CKGR openPO run:
every legacy row removed by v1.5.1 must have a twin (same sat and tile, an adjacent
strip t_min +/- 1) that v1.5.1 keeps, unless its strip touches a window edge."""
import json, sys
import pandas as pd
legacy = pd.read_csv("/home/claude/realdata/ckgr_openPO_refined_events.csv")
new = pd.read_csv("/home/claude/runs/orekit_ckgr_openPO_v151caps/arf_refined_events.csv")
for d in (legacy, new):
    d["sat"] = d["sat"].astype(str); d["tile"] = d["tile"].astype(str)
win_min = 72 * 60
k = ["sat", "tile", "t_min"]
m = legacy.merge(new[k].assign(kept=1), on=k, how="left")
removed = m[m.kept.isna()].copy()
keep_keys = set(map(tuple, new[k].itertuples(index=False, name=None)))
removed["twin"] = [((s, t, tm - 1) in keep_keys) or ((s, t, tm + 1) in keep_keys)
                   for s, t, tm in removed[k].itertuples(index=False, name=None)]
removed["edge"] = (removed.t_min <= 0) | (removed.t_min >= win_min)
orph = removed[~removed.twin & ~removed.edge]
# per-pass tile sets: a pass = run of consecutive strip minutes of one satellite
def passes(df):
    out = {}
    for s, g in df.groupby("sat"):
        mins = sorted(g.t_min.unique()); pid = {}; p = 0
        for i, tm in enumerate(mins):
            if i and tm - mins[i - 1] > 1: p += 1
            pid[tm] = p
        for (tm, t) in g[["t_min", "tile"]].itertuples(index=False, name=None):
            out.setdefault((s, pid[tm]), set()).add(t)
    return out
pl, pn = passes(legacy), passes(new)
lost = {kk: len(pl[kk] - pn.get(kk, set())) for kk in pl if pl[kk] - pn.get(kk, set())}
dup_legacy = int(legacy.duplicated(["sat", "tile"]).sum())
res = {"legacy_rows": len(legacy), "v151_rows": len(new), "removed": len(removed),
       "removed_with_twin": int(removed.twin.sum()), "removed_at_window_edge": int(removed.edge.sum()),
       "removed_orphans": len(orph), "passes_legacy": len(pl), "passes_v151": len(pn),
       "passes_losing_tiles": len(lost), "tiles_lost_total": int(sum(lost.values())),
       "legacy_rows_sharing_sat_tile_with_another_row": dup_legacy,
       "v151_rows_same_sat_tile_adjacent_minute": int(sum(((s, t, tm + 1) in keep_keys)
                                                           for s, t, tm in new[k].itertuples(index=False, name=None))),
       "removed_by_platform": removed.groupby(removed.sat.str.extract(r"^(\D+)")[0]).size().to_dict()}
print(json.dumps(res, indent=1, default=str))
if len(orph):
    print(orph.head(10).to_string())
json.dump(res, open("/home/claude/validate/logs/endcap_real_check.json", "w"), indent=1, default=str)
