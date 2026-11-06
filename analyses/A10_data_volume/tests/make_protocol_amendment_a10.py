"""Append amendment A10 (POST HOC correction of A6: per-tile H3 cell areas) to the A9 protocol.
Usage: make_protocol_amendment_a10.py PROTOCOL_A9 SCRIPT_A10 OUT_JSON"""
import sys, json, copy, hashlib, datetime
src_p, scr, out = sys.argv[1:4]
p = json.load(open(src_p))
assert [a["id"] for a in p["amendments"]] == [f"A{i}" for i in range(1, 10)], "input is not the A9 protocol"
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
new = copy.deepcopy(p); new["amended_utc"] = now
new["amendments"].append({
    "id": "A10", "date_utc": now,
    "status_when_amended": "POST HOC CORRECTION. Written after the A6 results were read, on finding an error in A6.",
    "corrects": "A6",
    "error": ("A6 gave every usable tile the global mean area of an H3 resolution-7 cell (5.16 km2). Over the "
              "CKGR+Okavango region the cells average 5.96 km2 (12,231 cells, 72,853 km2 by h3.cell_area), so every "
              "A6 volume was about 15% too low. Shares, percentage reductions and the PlanetScope share are almost "
              "unaffected (areas vary 5.85-6.07 km2 across the region)."),
    "correction": "Each image's footprint area is the sum of its usable tiles' own areas (h3.cell_area, km2). "
                  "All other A6 definitions, rates and sources are unchanged. A6 results are superseded by A10.",
    "scope": "Read-only replay of the 8 primary windows, as A6. No decision changes; reported as post hoc.",
    "script": {"arf_data_volume_a10.py": hashlib.sha256(open(scr, "rb").read()).hexdigest()},
})
json.dump(new, open(out, "w"), indent=1); print("A10 written", now)
