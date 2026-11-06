"""ECC (Cell 7C, earliest_completion) finds the earliest delivery time at which a
greedy k-image cover reaches the target by BISECTION over delivery times. That
is exact only if 'the greedy cover succeeds' is monotone in time (adding later
images never makes the greedy fail). Greedy maximum coverage is not monotone in
its ground set in general, so this checks it on the real 12 Oct runs: for every
onset (Cell 7C's subsets) and both skies, the bisection result against a linear
scan over all delivery times (the true earliest success of the same greedy).
Usage: test_ecc_monotone.py PIPELINE run1 [run2 ...]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, "/home/claude/validate")
import replay150 as R

PIPE = sys.argv[1]
runs = sys.argv[2:] or ["ckgr_openPO", "lz_openPO", "karoo_openPO", "ckgr_open"]
report = {}
all_ok = True
for name in runs:
    ns = R.run(PIPE, name, cells=("6",), overrides={"cloud_seeds": 2})
    import ast
    _c7c = R.cells_of(open(PIPE, encoding="utf-8").read())["7C"]
    _tree = ast.parse(_c7c)
    _tree.body = [n for n in _tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    ns.setdefault("H1_CONFIDENCE_ALPHA", 0.05)
    exec(compile(_tree, "cell7c_defs", "exec"), ns)
    st, mask = ns["event_data"], ns["cloud_mask"]
    k = ns["MAX_SELECTION_STEPS"]
    lvl = float(ns["H1_REACHABLE_COVERAGE_PCT"])
    horizon = (ns["WINDOW_HOURS"] - ns["ONSET_SHIFT_HOURS"]) * 60
    n_cmp = n_diff = 0
    worst = 0.0
    diffs = []
    t0 = time.time()
    for d in range(0, ns["ONSET_SHIFT_HOURS"] * 60 + 1, 60):
        for sky in ("clear", "cloud"):
            sub, sm = ns["subset_by_time"](st, mask if sky == "cloud" else None, d, d + horizon)
            if len(sub) == 0:
                continue
            reach = ns["reachable_tile_count"](sub, sm)
            if not reach:
                continue
            tb, selb, nb = ns["earliest_completion"](sub, sm, k, lvl, denominator=reach)
            # linear scan: earliest delivery time at which the same greedy succeeds
            util_eff = ns["_effective_util"](sub, sm)
            dd = ns["event_delivery_times"](sub)
            order = np.argsort(dd, kind="stable")
            order = order[np.isfinite(dd[order])]
            times = dd[order]
            need = int(np.ceil(lvl / 100.0 * reach - 1e-9))
            tl = None
            last_t = None
            for pos in range(len(order)):
                if last_t is not None and times[pos] == last_t:
                    continue
                last_t = times[pos]
                end = int(np.searchsorted(times, times[pos], side="right"))
                n_cov, chosen = ns["_count_cover"](sub, util_eff, order[:end], k, need)
                if n_cov >= need:
                    tl = float(max(dd[e] for e in chosen))
                    break
            n_cmp += 1
            if (tb is None) != (tl is None) or (tb is not None and abs(tb - tl) > 1e-9):
                n_diff += 1
                diffs.append({"onset_min": d, "sky": sky, "bisection": tb, "linear": tl})
                if tb is not None and tl is not None:
                    worst = max(worst, tb - tl)
    ok = n_diff == 0
    all_ok &= ok
    report[name] = {"comparisons": n_cmp, "differences": n_diff, "worst_delay_min": worst,
                    "examples": diffs[:5], "seconds": round(time.time() - t0, 1)}
    print(name, json.dumps(report[name]), flush=True)
os.makedirs("/home/claude/validate/logs", exist_ok=True)
json.dump(report, open("/home/claude/validate/logs/test_ecc_monotone.json", "w"), indent=1)
print("ALL MONOTONE" if all_ok else "BISECTION DIFFERS FROM LINEAR SCAN")
