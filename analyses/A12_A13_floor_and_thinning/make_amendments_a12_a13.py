"""Append post hoc amendments A12 (floor decomposition) and A13 (PlanetScope
acquisition fraction) to the A11 protocol. Refuses any other input protocol.
Usage: make_amendments_a12_a13.py IN_A11_PROTOCOL OUT_PROTOCOL"""
import sys, json, hashlib, datetime

SRC, DST = sys.argv[1], sys.argv[2]
raw = open(SRC, "rb").read()
assert hashlib.sha256(raw).hexdigest() == "bb88fa3217c29403929174f838bab8c5cec4bc3f9282357009753ff335ecb3a6", \
    "input is not the A11 protocol"
p = json.loads(raw)
assert [a["id"] for a in p["amendments"]] == [f"A{i}" for i in range(1, 12)]
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
windows = ("the 8 primary dry-season windows (3 Aug to 21 Sep 2026) and the 8 wet-season windows of A5 "
           "(5 Jan to 23 Feb 2026, p = 0.35 runs)")
p["amendments"] += [
    {"id": "A12", "date_utc": now,
     "status_when_amended": ("POST HOC. Written after the confirmatory results and A4 to A11 were read. While the script "
                             "was tested, it was run on the window of 3 August 2026 (sandbox copy of that run): median "
                             "floor 12.1 h = W 10.6 h + G 0.2 h + L 0 h. No other window has been decomposed."),
     "question": "What sets the physical floor: the wait for the first usable pass, the coverage geometry, or delivery?",
     "method": ("Each onset's floor (95% of reachable cells, clear sky, 24-h horizon) split exactly into W (onset to the "
                "first acquisition with a usable cell), G (that acquisition to the acquisition of the pass that completes "
                "the level) and L (acquisition to delivery of that pass); sun elevation at the region's centroid at the "
                "onset from the pipeline's NOAA routine. Every floor and reachable count is recomputed from "
                "arf_refined_events.csv and must equal arf_onset_robustness.csv, or the script stops."),
     "scope": f"Read-only replay of {windows}. Descriptive; no hypothesis, no decision changes.",
     "script": {"arf_floor_decomposition.py": sha("arf_floor_decomposition.py")}},
    {"id": "A13", "date_utc": now,
     "status_when_amended": ("POST HOC. Written after the confirmatory results and A4 to A11 were read, to bound the first "
                             "limitation (every daylight PlanetScope pass assumed to acquire). While the script was "
                             "tested, it was run on the window of 3 August 2026 (sandbox copy of that run, 20 seeds): "
                             "H1a share 100% at every fraction; H1b share 100, 96.5, 84.0 and 52.9% at f = 1, 0.75, "
                             "0.5 and 0.25; median floor unchanged at 12.1 h. No other window has been replayed."),
     "question": ("How do the floor, the reachable share and the H1a and H1b shares change if only a fraction f of "
                  "PlanetScope's predicted daylight passes acquire an image?"),
     "method": ("A PlanetScope pass is a maximal run of consecutive minutes of one FLOCK satellite. Per seed (1 to 20) "
                "each pass gets one uniform draw u and acquires when u < f (nested fractions f = 1, 0.75, 0.5, 0.25); "
                "Sentinel-2 and Landsat passes are always kept. Per onset, the floor, reachable count, ECC (k = 30) "
                "and FCFS image count are recomputed with the rules of Cell 7C; at f = 1 every onset must equal the "
                "run's arf_onset_robustness.csv (floor, reachable, ECC time, ECC images, FCFS images), or the script "
                "stops. H1a and H1b shares by the protocol's rules, margin = the run's t_max_wait_min."),
     "scope": (f"Read-only replay of {windows}. Descriptive sensitivity analysis; no hypothesis, no decision changes. "
               "It bounds the acquisition assumption from the model side; it does not measure Planet's acquisition "
               "fraction, which needs catalogue metadata."),
     "script": {"arf_planet_thinning.py": sha("arf_planet_thinning.py")}},
]
p["amended_utc"] = now
open(DST, "w").write(json.dumps(p, indent=1))
print(DST, now, hashlib.sha256(open(DST, "rb").read()).hexdigest())
