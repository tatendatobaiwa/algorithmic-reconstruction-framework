# =============================================================================
# ARF BATCH AUDIT (v1.5.4): read-only check of every protocol run folder
# =============================================================================
# For each of the 56 runs of the confirmatory protocol (PRIMARY, SECONDARY,
# EXPLORATORY, as listed in arf_batch_runner.py) this script checks the run
# folder itself, independently of the batch index:
#   1. the manifest exists and records release 1.5.3 with every cell run from it
#      (cell_versions.consistent);
#   2. its config_overrides equal the protocol configuration of that run;
#   3. Cell 9 found every expected output (output_validation.all_present);
#   4. every output file still has the SHA-256 that Cell 9 recorded, so no file
#      has been rewritten since by another run or session;
#   5. the index row's pipeline_sha256 equals the protocol's release file, and the
#      manifest was written within 15 minutes before the index row's finish time
#      (the indexed run is the one that wrote the folder).
# Input caches (shared element sets) are compared too, but reported only: a shared
# cache may be refreshed legitimately by a later run.
# v1.5.4: when ARF_code holds arf_batch_runner_v154.py and the protocol carries
# amendment A3, the 27 runs of the dense-constellation arm (batch
# exploratory_v154_dense) are checked the same way against the arm's release and
# pipeline file (A3: release_of_arm, pipeline_sha256_of_arm).
# Nothing is written except the audit table <BATCH_DIR>/confirmatory_v153_audit.csv.
#
# Use (Colab):  !python /content/drive/MyDrive/ARF_code/arf_batch_audit.py
#           or paste this file into a cell and run it
# =============================================================================
import os
import sys
import json
import runpy
import hashlib
import argparse
from datetime import datetime, timezone

import pandas as pd

BATCHES = (("PRIMARY", "confirmatory_v153_primary"), ("SECONDARY", "confirmatory_v153_secondary"),
           ("EXPLORATORY", "confirmatory_v153_exploratory"))
MAX_LAG_MIN = 15.0


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def norm(d):
    return json.loads(json.dumps(d, sort_keys=True, default=str))


def ts(s):
    t = datetime.fromisoformat(str(s).replace("Z", "+00:00"))
    return t if t.tzinfo else t.replace(tzinfo=timezone.utc)


def audit_run(cfg, rec, release, pipe_sha):
    """Return (verdict, reasons, notes) for one protocol run."""
    if not rec:
        return "NOT RUN", ["no index row"], []
    if rec.get("status") != "done":
        last = (rec.get("error") or "").strip().splitlines()
        return "NOT DONE", [f"index status {rec.get('status')}: {last[-1][:160] if last else ''}"], []
    out_dir = rec.get("out_dir", "")
    mpath = os.path.join(out_dir, "arf_reproducibility_manifest.json")
    if not os.path.exists(mpath):
        return "FAIL", [f"no manifest in {out_dir}"], []
    m = json.load(open(mpath))
    bad, notes = [], []
    if m.get("framework_version") != release:
        bad.append(f"manifest release {m.get('framework_version')} (expected {release})")
    cv = m.get("cell_versions") or {}
    if not cv.get("consistent") or cv.get("release") != release:
        bad.append(f"cells not all from {release}: missing {cv.get('missing')}, other {cv.get('other_release')}")
    if norm(m.get("config_overrides") or {}) != norm(cfg):
        bad.append(f"config {m.get('config_overrides')} differs from protocol {cfg}")
    if rec.get("pipeline_sha256") != pipe_sha:
        bad.append(f"index pipeline sha {str(rec.get('pipeline_sha256'))[:16]} is not the protocol release file")
    ov = m.get("output_validation") or {}
    if not ov.get("all_present"):
        bad.append(f"Cell 9 found {ov.get('n_missing')} expected outputs missing")
    changed, gone = [], []
    for name, f in (ov.get("files") or {}).items():
        if f.get("status") != "FOUND":
            continue
        p = os.path.join(out_dir, name)
        if not os.path.exists(p):
            gone.append(name)
        elif sha(p) != f.get("sha256"):
            changed.append(name)
    if gone:
        bad.append(f"{len(gone)} outputs deleted since Cell 9: {gone[:5]}")
    if changed:
        bad.append(f"{len(changed)} outputs rewritten since Cell 9: {changed[:5]}")
    try:
        lag = (ts(rec.get("finished_utc")) - ts(m.get("manifest_timestamp"))).total_seconds() / 60.0
        if not (-1.0 <= lag <= MAX_LAG_MIN):
            bad.append(f"manifest written {lag:.1f} min before the index finish time (not this run's manifest?)")
    except Exception as e:
        bad.append(f"timestamps unreadable ({type(e).__name__})")
    old_root = ov.get("out_dir") or out_dir
    for name, f in ((m.get("inputs") or {}).get("files") or {}).items():
        if f.get("status") != "FOUND":
            continue
        p = f.get("path", "")
        if not os.path.exists(p) and p.startswith(old_root):
            p = out_dir + p[len(old_root):]
        if not os.path.exists(p):
            notes.append(f"input {name} absent now")
        elif sha(p) != f.get("sha256"):
            notes.append(f"input {name} changed since the run")
    return ("PASS" if not bad else "FAIL"), bad, notes


def main(argv=None):
    ap = argparse.ArgumentParser(description="Read-only audit of the ARF 1.5.3 protocol run folders.")
    ap.add_argument("--code", default="/content/drive/MyDrive/ARF_code")
    ap.add_argument("--batch-dir", default="/content/drive/MyDrive/ARF_batch")
    # parse_known_args: also works when the code is pasted into a notebook cell,
    # where Jupyter passes its own '-f <kernel.json>' argument
    a, _ = ap.parse_known_args(argv)
    runner = runpy.run_path(os.path.join(a.code, "arf_batch_runner.py"), run_name="arf_batch_runner")
    proto = json.load(open(os.path.join(a.code, "arf_confirmatory_protocol.json")))
    release, pipe_sha = proto["release"], proto["pipeline_sha256"]
    specs = [(runner, list_name, batch, release, pipe_sha) for list_name, batch in BATCHES]
    a3 = [x for x in proto.get("amendments", []) if x.get("id") == "A3"]
    r154 = os.path.join(a.code, "arf_batch_runner_v154.py")
    if a3 and os.path.exists(r154):                     # v1.5.4: the dense-constellation arm
        specs.append((runpy.run_path(r154, run_name="arf_batch_runner_v154"), "DENSE", "exploratory_v154_dense",
                      a3[0]["release_of_arm"], a3[0]["pipeline_sha256_of_arm"]))
    rows = []
    for runner, list_name, batch, release, pipe_sha in specs:
        ip = os.path.join(a.batch_dir, f"{batch}_index.csv")
        idx = {}
        if os.path.exists(ip):
            df = pd.read_csv(ip, dtype=str).fillna("")
            idx = {r["run_id"]: dict(r) for _, r in df.iterrows()}
        runs = runner[list_name]
        ids = [runner["run_id"](c) for c in runs]
        for extra in sorted(set(idx) - set(ids)):
            rows.append({"batch": batch, "run_id": extra, "verdict": "NOT IN PROTOCOL", "reasons": "index row "
                         "for a configuration the protocol does not list", "notes": "", "out_dir": idx[extra].get("out_dir")})
        for cfg, rid in zip(runs, ids):
            rec = idx.get(rid, {})
            try:
                v, bad, notes = audit_run(cfg, rec, release, pipe_sha)
            except Exception as e:
                v, bad, notes = "FAIL", [f"audit error {type(e).__name__}: {e}"], []
            rows.append({"batch": batch, "run_id": rid, "verdict": v, "reasons": "; ".join(bad),
                         "notes": "; ".join(notes), "out_dir": rec.get("out_dir", "")})
    out = pd.DataFrame(rows)
    out_path = os.path.join(a.batch_dir, "confirmatory_v153_audit.csv")
    out.to_csv(out_path, index=False)
    print(f"Audit of release {proto['release']} (sha256 {proto['pipeline_sha256'][:16]}) protocol runs"
          + (f", and the dense-constellation arm, release {specs[-1][3]} (sha256 {specs[-1][4][:16]})"
             if len(specs) > len(BATCHES) else "") + ":")
    for _, _, batch, _, _ in specs:
        s = out[out.batch == batch].verdict.value_counts().to_dict()
        print(f"  {batch}: " + ", ".join(f"{k} {v}" for k, v in sorted(s.items())))
    for _, r in out[out.verdict != "PASS"].iterrows():
        print(f"  [{r.verdict}] {r.run_id}: {r.reasons}")
    n_notes = int((out.notes != "").sum())
    if n_notes:
        print(f"  {n_notes} runs with notes on shared input caches (not failures); see the table.")
    n_pass = int((out.verdict == "PASS").sum())
    print(f"{n_pass}/{sum(len(r[l]) for r, l, _, _, _ in specs)} protocol runs PASS. Table: {out_path}")
    return out


if __name__ == "__main__":
    main()
