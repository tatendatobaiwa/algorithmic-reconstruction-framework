"""Results tables generated from the result files (no number typed by hand)."""
import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..")) if _o.path.basename(_o.path.dirname(_o.path.abspath(__file__))) in ("figscripts", "checks") else None; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
from paths import D as _D, A as _A, FACTS as _FACTS
import json, shutil
D, A = _D, _A
conf = json.load(open(D + "arf_confirmatory_results.json"))
names = {("ckgr_okavango", "earth-resources+planet", "conservative", None): "Primary: CKGR and Okavango, full fleet",
         ("ckgr_okavango", "earth-resources+planet", "optimistic", None): "Optimistic PlanetScope swath",
         ("ckgr_okavango", "earth-resources", "conservative", None): "Open-data fleet only",
         ("lower_zambezi", "earth-resources+planet", "conservative", None): "Lower Zambezi NP",
         ("karoo_nature_reserve", "earth-resources+planet", "conservative", None): "Karoo Nature Reserve",
         ("ckgr_okavango", "earth-resources+planet", "conservative", "SENTINEL-2A"): "Full fleet without Sentinel-2A$^{\\dagger}$",
         ("ckgr_okavango", "earth-resources", "conservative", "SENTINEL-2A"): "Open data without Sentinel-2A$^{\\dagger}$"}
short = {"SUPPORTED": "S", "NOT SUPPORTED": "NS", "NOT EVALUABLE": "NE"}
def cell(x, d=1):
    lcb = [v for k, v in x.items() if k.startswith("lcb")][0]
    b = "--" if lcb is None else f"{lcb:.{d}f}".replace("-", "$-$")
    return f"{x['mean']:.{d}f} [{b}] {short[x['status']]} ({x['evaluable_windows']})"
L = [r"\begin{table}[H]", r"\centering", r"\footnotesize", r"\renewcommand{\arraystretch}{1.25}", r"\setlength{\tabcolsep}{3pt}",
     r"\caption{Confirmatory decisions (protocol ARF-confirmatory-2026-10, release 1.5.3; eight dry-season windows). "
     r"Each cell: mean over evaluable windows [one-sided lower confidence bound], decision (S supported, NS not supported, "
     r"NE not evaluable) and (number of evaluable windows). H1a and H1b in \% of onsets (threshold 90); H2 in percentage "
     r"points of the region's cells (threshold 0). Bound at level $1-\alpha/3$ for the primary design and $1-\alpha$ "
     r"otherwise, $\alpha = 0.05$. $^{\dagger}$Exploratory designs. $^{\S}$Post hoc replication (A5): the primary design in "
     r"eight wet-season windows, 5 January to 23 February 2026, under the same rules and bound ($1-\alpha/3$); "
     r"R: the rule is met, reported without a decision.}",
     r"\label{tab:confirm}", r"\begin{tabular}{|>{\raggedright\arraybackslash}p{3.6cm}|c|c|c|}", r"\hline",
     r"Design & H1a latency & H1b economy & H2 Fail-Signal \\", r"\hline"]
# Post hoc replication (A5): the primary design in eight wet-season windows, same rules, no decision.
a5 = json.load(open(D + "a5_summary.json"))["by_variant"]["p035"]
def rcell(x, d=1):
    return f"{x['mean']:.{d}f} [{x['lcb_alpha_over_3']:.{d}f}] {'R' if x['lcb_alpha_over_3'] >= (90 if d < 2 else 1e-12) else 'n'} ({x['windows']})"
for d in conf["designs"]:
    g = d["design"]; ex = (g.get("exclude_satellites") or [None])[0]
    L += [f"{names[(g['active_roi'], g['fleet'], g['planet_swath_case'], ex)]} & {cell(d['H1a'], 0)} & {cell(d['H1b'], 1)} & {cell(d['H2'], 2)} \\\\", r"\hline"]
    if d["role"] == "primary":
        L += [f"Primary, wet season$^{{\\S}}$ & {rcell(a5['h1a_pct'], 0)} & {rcell(a5['h1b_pct'], 1)} & {rcell(a5['h2_diff'], 2)} \\\\", r"\hline"]
L += [r"\end{tabular}", r"\end{table}"]
open("tables/tab_confirm.tex", "w").write("\n".join(L) + "\n")

a7 = json.load(open(A + "a7results/a7_summary.json"))["arms"]
rows = [("Primary: CKGR+Okavango, all 96, conservative", "Primary (dry season)"),
        ("Primary design, wet season (A5)", "Primary (wet season)$^{\\ddagger}$"),
        ("PlanetScope optimistic swath", "Optimistic swath"), ("Without Sentinel-2A", "Without Sentinel-2A"),
        ("Earth-resources only (no PlanetScope)", "Open-data fleet only"), ("Lower Zambezi", "Lower Zambezi NP"),
        ("Karoo Nature Reserve", "Karoo Nature Reserve")]
def m(v, k, d=1): return f"{v[k]['mean_over_windows']:.{d}f}"
def ci(v, k, d=1):
    lo, hi = v[k]["ci95_window_bootstrap"]; return f"{v[k]['mean_over_windows']:.{d}f} \\newline [{lo:.{d}f}, {hi:.{d}f}]"
L = [r"\begin{table}[H]", r"\centering", r"\footnotesize", r"\renewcommand{\arraystretch}{1.25}", r"\setlength{\tabcolsep}{3pt}",
     r"\caption{What the fleet can deliver within 24~h of an onset, and the images each selector needs (post hoc, amendment A7; clear sky). "
     r"Per window the median over onsets is taken; the table gives the mean over the eight windows with a window-block bootstrap "
     r"95\% interval. Floor: earliest time at which all passes together cover 95\% of the reachable cells. Reach: share of the region "
     r"any pass can image in 24~h. Region: share of onsets at which 95\% of the whole region is covered within 24~h. Images at the H1a level "
     r"for ECC, coverage-first greedy (CF), time-aware greedy (TA) and FCFS. Saving: mean per-onset images saved by ECC. "
     r"$^{\ddagger}$Exploratory (amendment A5).}",
     r"\label{tab:benchmark}", r"\begin{tabular}{|>{\raggedright\arraybackslash}p{2.4cm}|>{\centering\arraybackslash}p{1.45cm}|>{\centering\arraybackslash}p{1.25cm}|>{\centering\arraybackslash}p{1.25cm}|c|c|c|c|>{\centering\arraybackslash}p{1.45cm}|>{\centering\arraybackslash}p{1.45cm}|}", r"\hline",
     r"Design & Floor (h) & Reach (\%) & Region (\%) & ECC & CF & TA & FCFS & Saving vs CF & Saving vs FCFS \\", r"\hline"]
for key, lab in rows:
    v = a7[key]
    L += [f"{lab} & {ci(v, 'floor_h_median')} & {ci(v, 'reach_pct_median', 0)} & {ci(v, 'region_t95_reached_pct', 0)} & "
          f"{m(v, 'ecc_images_median')} & {m(v, 'coverage_first_images_median')} & {m(v, 'time_aware_images_median')} & "
          f"{m(v, 'fcfs_unbudgeted_images_median')} & {ci(v, 'ecc_saving_vs_coverage_first_mean')} & {ci(v, 'ecc_saving_vs_fcfs_mean')} \\\\", r"\hline"]
L += [r"\end{tabular}", r"\end{table}"]
open("tables/tab_benchmark.tex", "w").write("\n".join(L) + "\n")
for src, dst in ((A + "a4results/tab_a4_exact_economy.tex", "tab_a4_exact_economy.tex"), (A + "a9results/tab_a9_budget.tex", "tab_a9_budget.tex"),
                 (A + "a9results/tab_a9_failsignal.tex", "tab_a9_failsignal.tex"), (A + "a8results/tab_a8_pooled.tex", "tab_a8_pooled.tex"),
                 (A + "a11results/tab_a11_volume.tex", "tab_a11_volume.tex")):
    s = open(src).read().replace("table*", "table").replace("[!t]", "[H]").replace("& -0.39, 0.05 &", "& $-$0.39, 0.05 &").replace("cloudy.}", "cloudy. Contains modified Copernicus Sentinel data 2026.}")
    open("tables/" + dst, "w").write(s)
print(open("tables/tab_confirm.tex").read())
