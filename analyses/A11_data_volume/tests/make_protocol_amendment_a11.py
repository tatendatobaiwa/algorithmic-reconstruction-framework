"""Append amendment A11 (POST HOC correction of the A6/A10 per-sensor rates) to the A10 protocol.
Usage: make_protocol_amendment_a11.py PROTOCOL_A10 SCRIPT_A11 OUT_JSON"""
import sys, json, copy, hashlib, datetime
src_p, scr, out = sys.argv[1:4]
p = json.load(open(src_p))
assert [a["id"] for a in p["amendments"]] == [f"A{i}" for i in range(1, 11)], "input is not the A10 protocol"
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
new = copy.deepcopy(p); new["amended_utc"] = now
new["amendments"].append({
    "id": "A11", "date_utc": now,
    "status_when_amended": "POST HOC CORRECTION. Written after the A10 results were read, on checking the A6 rate sources against their documents while drafting the Methods.",
    "corrects": ["A10", "A6"],
    "error": {
        "PlanetScope": ("The upper rate (2.43 MB/km2) divided 700 MB by the 287.5 km2 frame of PS2 (Dove Classic, 25 x 11.5 km), "
                        "whose imagery ended in April 2022 (Planet, PlanetScope Product Specifications, June 2026, Sec. 2.1, Table 2-A). "
                        "The pipeline models the current fleet as PS2.SD or PSB.SD (Cell 0, MISSION_TABLE)."),
        "Landsat": ("The scene area (185 x 180 km) is not stated in the source cited (LSDS-1619 v6.0)."),
    },
    "correction": {
        "rates_mb_per_km2": {"Sentinel-2": [0.066, 0.080], "Landsat": [0.022, 0.038], "PlanetScope": [0.78, 1.22]},
        "PlanetScope": "500-700 MB per 4-band scene (Table 3-B) over the PS2.SD 25 x 23.0 km and PSB.SD 32.5 x 19.6 km frames (Table 3-C): 500/637 and 700/575.",
        "Landsat": ("LSDS-1619 v6.0 Table A-1 sizes (674,478 KB if the seven surface-reflectance bands are one 83,000 KB entry; 1,172,478 KB if per band) "
                    "over the USGS EROS Landsat 8 scene size, 'approximately 170 km north-south by 183 km east-west' (catalog.data.gov/dataset/landsat-8)."),
        "Sentinel-2": "unchanged",
        "also": "per-sensor usable footprint areas (km2) are written for ECC and FCFS, so a later rate change needs no replay",
    },
    "scope": "Read-only replay of the 8 primary windows, as A6 and A10. No decision changes; reported as post hoc. A11 supersedes A10 and A6.",
    "script": {"arf_data_volume_a11.py": hashlib.sha256(open(scr, "rb").read()).hexdigest()},
})
json.dump(new, open(out, "w"), indent=1); print("A11 written", now)
