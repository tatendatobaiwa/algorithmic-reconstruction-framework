"""With vs without Sentinel-2A (v1.5.2 replays, same inputs minus the S2A rows).
Usage: make_s2a_compare.py FULL_CSV WITHOUT_CSV OUT_CSV"""
import sys
import pandas as pd
full, wo, out = sys.argv[1:4]
a, b = pd.read_csv(full), pd.read_csv(wo)
key = ["ROI", "Fleet"]
cols = ["Refined rows", "Candidate events", "Reachable usable tiles, clear (%)", "Time-aware coverage, clear (%)",
        "Decay time scale tau (h)", "H1 (H1a and H1b)", "H1a status (ECC)", "H1a ECC within margin (%)",
        "H1a ECC 95% lower bound (%)", "H1b status", "H1b eligible onsets", "H1b ECC fewer images within margin (%)",
        "Median images to level, ECC", "Median images to level, FCFS",
        "H2 status (time-aware; adaptive vs open-loop, same images)", "H2 mean difference (coverage points)",
        "H2 one-sided 95% lower bound", "H3 status"]
m = a[key + cols].merge(b[key + cols], on=key, suffixes=(" with S2A", " without S2A"), how="inner")
order = key + [c + s for c in cols for s in (" with S2A", " without S2A")]
m[order].to_csv(out, index=False)
for _, r in m.iterrows():
    ch = [c for c in cols if str(r[c + " with S2A"]) != str(r[c + " without S2A"])]
    print(f"{r.ROI} / {r.Fleet}: " + "; ".join(f"{c}: {r[c + ' with S2A']} -> {r[c + ' without S2A']}" for c in ch))
