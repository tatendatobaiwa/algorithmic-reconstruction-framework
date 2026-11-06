"""Unit test of exact_time's branches in arf_exact_time_a9.py with a synthetic monotone feasibility (feasible from
time F onwards), delivery times 0..200 step 10, floor 30. Usage: test_exact_time_logic_a9.py SCRIPT"""
import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("a9", sys.argv[1]); a9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a9)
RES = []
def check(c, m): RES.append(bool(c)); print(("PASS " if c else "FAIL ") + m)
times = np.arange(0, 201, 10.0)
ns = {"event_delivery_times": lambda sub: times}
calls = []
def mk(F, unknown_at=None):
    def feas(ns_, ex, sub, t, need, k):
        calls.append(t)
        if unknown_at is not None and abs(t - unknown_at) < 1e-9:
            return None, None, "time limit"
        return (t >= F), None, "optimal"
    return feas
cases = [  # (F, t_e, expected T, verdict)
    (30, 30, 30, "ECC at floor"), (30, 90, 30, "ECC delay avoidable"), (70, 90, 70, "ECC later than exact"),
    (90, 90, 90, "ECC at exact time"), (1e9, None, None, "no k-cover exists"), (120, None, 120, "ECC missed a cover"),
    (40, None, 40, "ECC missed a cover")]
ok = True
for F, te, T, v in cases:
    a9.feasible = mk(F); calls.clear()
    got = a9.exact_time(ns, {}, None, 10, 5, 30.0, te)
    ok &= got[0] == T and got[1] == v and got[3] is False
    ok &= len(calls) <= 2 + int(np.ceil(np.log2(len(times))))
check(ok, "all six verdicts correct, exact time found by bisection within log2(n) + 2 solves")
a9.feasible = mk(70, unknown_at=60); r = a9.exact_time(ns, {}, None, 10, 5, 30.0, 90)
check(r[0] is None and r[3] is True and r[1].startswith("unknown"), "a solve hitting the time limit is reported as unknown")
a9.MAX_BISECT = 2; a9.feasible = mk(170); r = a9.exact_time(ns, {}, None, 10, 5, 30.0, 190)
check(r[0] is None and r[1] == "unknown (bisection cap)", "the bisection cap is reported as unknown, never as a result")
print(f"{sum(RES)}/{len(RES)} passed"); sys.exit(0 if all(RES) else 1)
