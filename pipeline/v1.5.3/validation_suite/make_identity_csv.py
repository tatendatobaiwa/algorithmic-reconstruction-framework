"""Summarise compare_v152_v153.py (JSON lines) as a CSV, one row per run.
Usage: make_identity_csv.py COMPARE_JSONL OUT_CSV"""
import sys, json
import pandas as pd

rows = []
for line in open(sys.argv[1]):
    r = json.loads(line)
    rows.append({
        "run": r["input"].rstrip("/").split("/")[-1], "inputs": "12 Oct 2026 (exact Cell 5 output of 1.5.2; ckgr_open "
        "from Colab)" if r["mode"] == "oct12" else "15 Oct 2026 Colab run (saved inputs from Drive)",
        "releases compared": " vs ".join(r.get("releases") or []),
        "identical": r.get("ok"),
        "CSV files byte-identical": f"{sum(r.get('csv', {}).values())}/{len(r.get('csv', {}))}",
        "JSON values compared": sum(v["values"] for v in r.get("json", {}).values()),
        "JSON values changed (other than keys 1.5.3 adds)": sum(len(v["unexpected"]) for v in r.get("json", {}).values()),
        "scenario_results fields compared": r.get("scenario_results_fields"),
        "scenario_results fields changed": len(r.get("scenario_results_differences") or []),
        "report lines changed (all label lines)": len(r.get("report_lines_removed") or []),
        "report lines changed otherwise": len(r.get("report_lines_unexplained") or []),
        "PNG figures byte-identical": len(r.get("png_same") or []),
        "PNG figures redrawn": "; ".join(f.replace(".png", "") for f in r.get("png_changed") or []),
        "seconds": r.get("seconds"),
    })
df = pd.DataFrame(rows)
df.to_csv(sys.argv[2], index=False)
print(df[["run", "identical", "CSV files byte-identical", "JSON values compared", "scenario_results fields compared",
          "PNG figures redrawn"]].to_string(index=False))
print(f"{int(df['identical'].sum())}/{len(df)} runs identical")
