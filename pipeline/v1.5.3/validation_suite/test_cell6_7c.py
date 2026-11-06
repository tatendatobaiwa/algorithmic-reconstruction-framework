"""Independent checks of Cells 6 and 7C (selection, cloud model, Fail-Signal
causality, exact ILP, statistics). Usage: test_cell6_7c.py PIPELINE [V140_PIPELINE]
Every expected value is computed here by a separate, plain implementation
(brute force, SciPy, or direct probability), never by calling the code under test
twice."""
import ast, re, sys, math, json, heapq, logging, itertools, random, os
import numpy as np, pandas as pd, h3
from scipy import stats
from scipy.stats import multivariate_normal

PIPE = sys.argv[1]
V140 = sys.argv[2] if len(sys.argv) > 2 else None
logging.disable(logging.WARNING)
RESULTS = []


def check(cond, msg):
    RESULTS.append((bool(cond), msg))
    print(("PASS " if cond else "FAIL ") + msg, flush=True)


def cells_of(src):
    hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
    ms = list(hdr.finditer(src))
    return {m.group(1): src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)] for i, m in enumerate(ms)}


def defs_only(code, ns):
    tree = ast.parse(code)
    tree.body = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    exec(compile(tree, "defs", "exec"), ns)
    return ns


def load(path, tile_ids, extra=None):
    C = cells_of(open(path, encoding="utf-8").read())
    tile_ids = np.array(tile_ids)
    rank = np.empty(len(tile_ids), dtype=np.int64)
    rank[np.argsort(tile_ids, kind="stable")] = np.arange(len(tile_ids))
    ns = dict(np=np, pd=pd, math=math, heapq=heapq, random=random, h3=h3, os=os, logger=logging.getLogger("t"),
              _tile_ids=tile_ids, _tile_string_rank=rank, total_target_tiles=len(tile_ids), T_MAX_WAIT_MINUTES=94.0,
              SWATH_HALF_ANGLE_DEG=9.22, H1_CONFIDENCE_ALPHA=0.05, COVERAGE_TARGETS_PCT=(95.0, 99.0),
              MAX_SELECTION_STEPS=30, TIME_AWARE_OBJECTIVE="facility", LATENCY_DECAY_SCALE="fixed",
              CLOUD_PROBABILITY=0.35, CLOUD_MODEL="iid", CLOUD_CORR_LENGTH_KM=50.0, CLOUD_CORR_TIME_H=3.0,
              CLOUD_PLANNING="expected", H1_REACHABLE_COVERAGE_PCT=95.0, CLOUD_PROBABILITY_SOURCE="constant",
              CLOUD_KNOWLEDGE="realistic", last_adaptive_run={})
    ns.update(extra or {})
    defs_only(C["6"], ns)
    defs_only(C["7C"], ns)
    return ns


def random_instance(rng, n_tiles, n_events, p_cover=0.35, waits=False):
    tiles = [f"t{i:03d}" for i in range(n_tiles)]
    rows = []
    for e in range(n_events):
        sat = f"S{rng.integers(0, 4)}"
        t = int(rng.integers(0, 2000))
        cov = [tl for tl in tiles if rng.random() < p_cover] or [tiles[int(rng.integers(0, n_tiles))]]
        theta = 10.0
        for tl in cov:
            ang = float(rng.uniform(0, 11.0))                    # some beyond theta: utility 0
            rows.append({"sat": sat, "tile": tl, "t_min": t, "A_nadir_deg": ang, "D_downlink": 1.0,
                         "T_wait_min": float(rng.integers(0, 30)) if waits else 0.0, "theta_max_deg": theta})
    df = pd.DataFrame(rows).drop_duplicates(["sat", "t_min", "tile"])
    return tiles, df


def rows_of(store, e):
    s = store.rows(e)
    return store.tile_idx[s], store.util[s]


def f_value(store, sel_idx, w, util=None):
    util = store.util if util is None else util
    best = {}
    for e in sel_idx:
        s = store.rows(e)
        for t, u in zip(store.tile_idx[s], util[s]):
            v = w[e] * u
            if v > best.get(t, 0.0):
                best[t] = v
    return sum(best.values())


def plain_greedy(store, k, w, util=None):
    """Greedy on f(S) = sum_t max_{e in S} w_e u_e,t; ties -> lowest event index; stop when no gain.
    A gain counts only above 1e-12 * w_e (rounding noise in the pass's own utility units), the
    rule stated in Cell 6; no absolute floor, since weights can be as small as 1e-300."""
    chosen, cur = [], {}
    util = store.util if util is None else util
    for _ in range(k):
        best_e, best_g = None, -math.inf
        for e in range(len(store)):
            if e in chosen:
                continue
            s = store.rows(e)
            g = sum(max(0.0, w[e] * u - cur.get(t, 0.0)) for t, u in zip(store.tile_idx[s], util[s]))
            if g > 1e-12 * w[e] and g > best_g:
                best_e, best_g = e, g
        if best_e is None:
            break
        chosen.append(best_e)
        s = store.rows(best_e)
        for t, u in zip(store.tile_idx[s], util[s]):
            cur[t] = max(cur.get(t, 0.0), w[best_e] * u)
    return chosen


rng = np.random.default_rng(20261014)

# ---------------------------------------------------------------- A, B, C: facility objective, guarantee, exact ILP
n_inst, min_ratio, ilp_bad, greedy_mismatch, cov_mismatch = 0, 1.0, [], [], []
for it in range(160):
    tiles, df = random_instance(rng, int(rng.integers(8, 25)), int(rng.integers(5, 11)))
    ns = load(PIPE, tiles)
    st = ns["build_event_data"](df)
    if len(st) < 3:
        continue
    n_inst += 1
    rate = float(rng.choice([0.5, 3.0, 5.0]))
    w = ns["time_aware_weights"](st, rate)
    k = int(rng.integers(1, min(5, len(st)) + 1))
    sel, _, _, _ = ns["run_submodular_greedy_selection"](st, k, gap_fill=False, time_aware=True, latency_decay_rate=rate)
    sel_idx = [st.index[s] for s in sel]
    ref = plain_greedy(st, k, w)
    if sel_idx != ref:
        greedy_mismatch.append((it, sel_idx, ref))
    # coverage-first = the same greedy with unit weights
    selc, _, _, _ = ns["run_submodular_greedy_selection"](st, k, gap_fill=False)
    if [st.index[s] for s in selc] != plain_greedy(st, k, np.ones(len(st))):
        cov_mismatch.append(it)
    # brute-force optimum of f over all k-subsets
    opt = max(f_value(st, c, w) for c in itertools.combinations(range(len(st)), min(k, len(st))))
    g = f_value(st, sel_idx, w)
    if opt > 0:
        min_ratio = min(min_ratio, g / opt)
    res = ns["solve_ilp_max_utility"](st, None, k, 60.0, 10 ** 6, weight=w)
    if res["status"] != "optimal" or abs(res["objective"] - opt) > 1e-7 * max(1.0, opt):
        ilp_bad.append((it, res["status"], res["objective"], opt))
check(not greedy_mismatch, f"facility greedy == independent greedy on f(S) = sum_t max w u ({n_inst} random instances)"
      + (f"; mismatches {greedy_mismatch[:3]}" if greedy_mismatch else ""))
check(not cov_mismatch, f"coverage-first greedy == independent unit-weight greedy ({n_inst} instances)")
check(min_ratio >= 1 - 1 / math.e - 1e-12, f"greedy / brute-force optimum >= 1 - 1/e on every instance (worst {min_ratio:.4f})")
check(not ilp_bad, f"weighted ILP (HiGHS) objective == brute-force optimum of f ({n_inst} instances)"
      + (f"; bad {ilp_bad[:3]}" if ilp_bad else ""))

# submodularity of f (random nested sets)
viol = 0
for it in range(300):
    tiles, df = random_instance(rng, 15, 9)
    ns = load(PIPE, tiles)
    st = ns["build_event_data"](df)
    w = ns["time_aware_weights"](st, 3.0)
    n = len(st)
    B = [e for e in range(n) if rng.random() < 0.5]
    A = [e for e in B if rng.random() < 0.5]
    rest = [e for e in range(n) if e not in B]
    if not rest:
        continue
    x = rest[0]
    gA = f_value(st, A + [x], w) - f_value(st, A, w)
    gB = f_value(st, B + [x], w) - f_value(st, B, w)
    viol += gA < gB - 1e-12
check(viol == 0, "f is submodular on 300 random nested pairs A subset B (gain of x at A >= gain at B)")

# ---------------------------------------------------------------- D: legacy reproduces v1.4.0
if V140:
    diffs = 0
    for it in range(80):
        tiles, df = random_instance(rng, 20, 9)
        n15 = load(PIPE, tiles, {"TIME_AWARE_OBJECTIVE": "legacy"})
        n14 = load(V140, tiles)
        s15, s14 = n15["build_event_data"](df), n14["build_event_data"](df)
        a = n15["run_submodular_greedy_selection"](s15, 5, gap_fill=False, time_aware=True, latency_decay_rate=3.0)[0]
        b = n14["run_submodular_greedy_selection"](s14, 5, gap_fill=False, time_aware=True, latency_decay_rate=3.0)[0]
        diffs += a != b
    check(diffs == 0, "TIME_AWARE_OBJECTIVE='legacy' selects exactly as v1.4.0 (80 instances)")
    # iid cloud mask with constant p identical to v1.4.0
    same = True
    for it in range(20):
        tiles, df = random_instance(rng, 30, 12)
        n15 = load(PIPE, tiles, {"CLOUD_MODEL": "iid"})
        n14 = load(V140, tiles)
        s15, s14 = n15["build_event_data"](df), n14["build_event_data"](df)
        m15 = n15["build_cloud_mask"](s15, 0.35, 42 + it, model="iid")
        m14 = n14["build_cloud_mask"](s14, 0.35, 42 + it)
        same &= np.array_equal(m15.flags, m14.flags)
    check(same, "cloud_model='iid' mask is bit-identical to v1.4.0's mask (20 seeds)")

# ---------------------------------------------------------------- E: expected-coverage planner
bad = 0
for it in range(60):
    tiles, df = random_instance(rng, 18, 9)
    ns = load(PIPE, tiles)
    st = ns["build_event_data"](df)
    p = rng.uniform(0.1, 0.8, len(tiles))
    unc = rng.random(len(tiles)) < 0.8
    cand = np.ones(len(st), dtype=bool)
    k = 4
    plan = ns["_expected_plan"](st, unc, cand, k, p)
    # independent greedy on E[#uncovered tiles with >= 1 clear usable view], independence across passes
    def E(sel):
        tot = 0.0
        for t in range(len(tiles)):
            if not unc[t]:
                continue
            q = 1.0
            for e in sel:
                s = st.rows(e)
                if np.any((st.tile_idx[s] == t) & (st.util[s] > 0)):
                    q *= p[t]
            tot += 1.0 - q
        return tot
    ref = []
    for _ in range(k):
        gains = [(E(ref + [e]) - E(ref), -e) for e in range(len(st)) if e not in ref]
        g, ne = max(gains)
        if g <= 1e-12:
            break
        ref.append(-ne)
    bad += plan != ref
check(bad == 0, "expected-coverage planner == independent greedy on E[sum_t 1 - prod p] (60 instances)")

# ---------------------------------------------------------------- F: causality of the adaptive selector
viol, tests = 0, 0
for planning in ("expected", "clear"):
    for ta in (False, True):
        for it in range(25):
            tiles, df = random_instance(rng, 20, 12, waits=True)
            ns = load(PIPE, tiles, {"CLOUD_PLANNING": planning})
            st = ns["build_event_data"](df)
            if len(st) < 4:
                continue
            m = ns["CloudMask"](st, rng.random(st.n_rows) < 0.4)
            sel = ns["simulate_adaptive_selection"](st, 4, m, gap_fill=True, time_aware=ta, latency_decay_rate=3.0,
                                                    planning=planning)[0]
            d = ns["event_delivery_times"](st)
            for j in range(1, len(sel) + 1):
                tj = d[st.index[sel[j - 1]]]
                later = np.repeat(d > tj, st.counts)               # rows of images arriving AFTER decision j
                flags = m.flags.copy()
                flags[later] = rng.random(int(later.sum())) < 0.5
                sel2 = ns["simulate_adaptive_selection"](st, 4, ns["CloudMask"](st, flags), gap_fill=True, time_aware=ta,
                                                         latency_decay_rate=3.0, planning=planning)[0]
                tests += 1
                viol += sel2[:j] != sel[:j]
check(viol == 0, f"adaptive selector is causal: decisions up to each arrival unchanged when the cloud of later "
                 f"images is redrawn ({tests} checks, both planning modes, both objectives)")

# ---------------------------------------------------------------- G: correlated cloud model statistics
tiles_real = sorted(pd.read_csv("/home/claude/realdata/lz_tiles.csv")["tile"])
ll = np.array([h3.cell_to_latlng(t) for t in tiles_real])
order = np.argsort(ll[:, 0], kind="stable")
tid = np.array(tiles_real)[order]
L, T = 50.0, 3.0
ns = load(PIPE, tid, {"CLOUD_MODEL": "correlated", "CLOUD_CORR_LENGTH_KM": L, "CLOUD_CORR_TIME_H": T})
field = ns["LatentCloudField"](L)
# (i) marginal unit variance and spatial correlation exp(-d^2/L^2) of the interpolated field
Z = np.array([field.at_tiles(field.draw(np.random.default_rng(s)), np.arange(len(tid))) / field.sd for s in range(2000)])
var_ok = abs(Z.var(axis=0).mean() - 1.0) < 0.03
lat, lon = ll[order, 0], ll[order, 1]
R = 6371.0088
def dist(i, j):
    p1, p2 = np.radians(lat[i]), np.radians(lat[j])
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon[j] - lon[i]) / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))
pairs = [(int(i), int(j)) for i, j in rng.integers(0, len(tid), (600, 2)) if i != j]
emp, thy = [], []
for i, j in pairs:
    d_ij = dist(i, j)
    if d_ij > 120:
        continue
    emp.append(np.corrcoef(Z[:, i], Z[:, j])[0, 1])
    thy.append(np.exp(-(d_ij / L) ** 2))
emp, thy = np.array(emp), np.array(thy)
err = np.abs(emp - thy)
check(var_ok and err.mean() < 0.03 and err.max() < 0.12,
      f"latent field: unit variance ({Z.var(axis=0).mean():.3f}); spatial correlation vs exp(-d^2/L^2), L = {L:g} km: "
      f"mean |error| {err.mean():.3f}, max {err.max():.3f} over {len(emp)} tile pairs (2,000 draws)")
# (ii) marginal cloud probability per tile and temporal joint probability for two passes dt apart
rows = []
dt = 90
for t in (100, 100 + dt):
    for tl in tid:
        rows.append({"sat": f"S{t}", "tile": tl, "t_min": t, "A_nadir_deg": 1.0, "D_downlink": 1.0, "T_wait_min": 0.0,
                     "theta_max_deg": 10.0})
st = ns["build_event_data"](pd.DataFrame(rows))
F = np.array([ns["build_cloud_mask"](st, 0.35, s, model="correlated", length_km=L, time_h=T).flags for s in range(1500)])
n_t = len(tid)
e0 = st.index[("S100", 100)]
e1 = st.index[(f"S{100 + dt}", 100 + dt)]
r0, r1 = st.rows(e0), st.rows(e1)
a0 = F[:, r0][:, np.argsort(st.tile_idx[r0])]
a1 = F[:, r1][:, np.argsort(st.tile_idx[r1])]
marg = a0.mean()
both = (a0 & a1).mean()
rho_t = math.exp(-dt / (60.0 * T))
thr = stats.norm.ppf(1 - 0.35)
p_both = multivariate_normal(mean=[0, 0], cov=[[1, rho_t], [rho_t, 1]]).cdf([-thr, -thr])   # P(Z1 > thr, Z2 > thr)
check(abs(marg - 0.35) < 0.01 and abs(both - p_both) < 0.01,
      f"cloud flags: marginal {marg:.4f} (target 0.35); P(cloudy at t and t+{dt} min) {both:.4f} vs bivariate-normal "
      f"{p_both:.4f} (rho = exp(-dt/T) = {rho_t:.3f}); 1,500 seeds x {n_t} tiles")
# (iii) iid limit: independent passes
Fi = np.array([ns["build_cloud_mask"](st, 0.35, s, model="iid").flags for s in range(800)])
b0, b1 = Fi[:, r0], Fi[:, r1]
check(abs(b0.mean() - 0.35) < 0.01 and abs((b0 & b1).mean() - 0.35 ** 2) < 0.01,
      f"iid model: marginal {b0.mean():.4f}, joint {(b0 & b1).mean():.4f} vs 0.35^2 = {0.35 ** 2:.4f}")

# ---------------------------------------------------------------- H: statistics helpers
ok_ps = True
for it in range(50):
    d = rng.normal(rng.uniform(-1, 1), rng.uniform(0.1, 2), int(rng.integers(3, 30)))
    ps = ns["paired_summary"](d)
    ref = stats.ttest_1samp(d, 0.0, alternative="greater").confidence_interval(0.95).low
    ok_ps &= abs(ps["lcb"] - ref) < 1e-9 and abs(ps["mean"] - d.mean()) < 1e-12
check(ok_ps, "paired_summary lower bound == SciPy one-sided 95% t interval (50 samples)")
ok_cp = True
for n in range(1, 80):
    for x in range(0, n + 1, max(1, n // 7)):
        lb = ns["cp_lower_bound"](x, n)
        ref = stats.binomtest(x, n).proportion_ci(confidence_level=0.90, method="exact").low
        ok_cp &= abs(lb - ref) < 1e-9
check(ok_cp and ns["min_onsets_for"](90.0) == 29, "Clopper-Pearson one-sided 95% bound == SciPy exact; 29 onsets for 90%")

# ---------------------------------------------------------------- I: open-loop plan, realised coverage
ok_ol = True
for it in range(30):
    tiles, df = random_instance(rng, 20, 10)
    ns2 = load(PIPE, tiles)
    st = ns2["build_event_data"](df)
    p = np.full(len(tiles), 0.35)
    k = 4
    a = ns2["plan_expected_coverage"](st, k, p)
    b = [st.keys_[e] for e in ns2["_expected_plan"](st, np.ones(len(tiles), bool), np.isfinite(ns2["event_delivery_times"](st)), k, p)]
    m = ns2["CloudMask"](st, rng.random(st.n_rows) < 0.35)
    rc = ns2["realized_coverage"](st, m, a)
    seen = set()
    for key in a:
        s = st.rows(st.index[key])
        for t, u, f in zip(st.tile_idx[s], st.util[s], m.flags[s]):
            if u > 0 and not f:
                seen.add(t)
    ok_ol &= a == b and rc == len(seen)
check(ok_ol, "plan_expected_coverage == the adaptive selector's planner run once; realized_coverage == brute force")

n_fail = sum(1 for ok, _ in RESULTS if not ok)
print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
json.dump([{"ok": ok, "check": m} for ok, m in RESULTS], open(os.environ.get("ARF_TEST_OUT", "/home/claude/validate/logs/test_cell6_7c.json"), "w"), indent=1)
sys.exit(1 if n_fail else 0)
