# =============================================================================
# ARF RESUME CELL (Colab, 24 October 2026): finish the real-fleet batch, then
# resume the dense-constellation arm in an order that survives disconnections
# =============================================================================
# PART = "eo"     the 7 failed real-fleet runs, release 1.5.3, through
#                 arf_batch_runner.py itself (the 49 complete runs are skipped after
#                 their manifests are checked); the failure record is copied first;
#                 then the audit. About 40 min. Do NOT use the 1.5.3 launch cell any
#                 more: it would put the A2 protocol back over A3.
# PART = "dense"  the remaining runs of the arm, release 1.5.4, same 27
#                 configurations and same index (exploratory_v154_dense), run in
#                 this order: windows 2-8 at N = 100, 500, 1000 (21 runs, each under
#                 about 2 h), then window 1 at N = 4000 (about 7.5 h), then N = 8000
#                 (about 15 h in one piece) unless INCLUDE_8000 = False. Completed
#                 runs are skipped; a run cut off part-way starts again.
# One notebook only; keep the tab open. Secrets are tested for presence only.
# =============================================================================
PART = "eo"                                       # "eo", then "dense"
INCLUDE_8000 = True
ROOT = "/content/drive/MyDrive"
CODE = f"{ROOT}/ARF_code"
BATCH = f"{ROOT}/ARF_batch"
SHA_153 = "14426c89c916ffcaac032cb94ec186ef7a7f8011748044c99f2386585d917950"
SHA_154 = "2b618ec442e71cf896e1104c6da89c94f2d6b6785a04c4a9236b92c21115a1cb"

import os, sys, json, shutil, hashlib, runpy, subprocess

if ROOT.startswith("/content/drive") and not os.path.isdir("/content/drive/MyDrive"):
    from google.colab import drive
    drive.mount("/content/drive")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


assert sha(os.path.join(CODE, "algorithmic_reconstruction_framework.py")) == SHA_153, "ARF_code does not hold release 1.5.3"
assert sha(os.path.join(CODE, "algorithmic_reconstruction_framework_v154.py")) == SHA_154, "ARF_code does not hold release 1.5.4"
proto = json.load(open(os.path.join(CODE, "arf_confirmatory_protocol.json")))
assert [a["id"] for a in proto["amendments"]][-1] == "A3" and proto["pipeline_sha256"] == SHA_153, \
    "ARF_code/arf_confirmatory_protocol.json is not the A3 protocol (was the 1.5.3 launch cell run again?)"
try:
    from google.colab import userdata
    missing = [k for k in ("SPACETRACK_USERNAME", "SPACETRACK_PASSWORD") if not userdata.get(k)]
except Exception as e:
    missing = [f"Colab secrets unavailable ({type(e).__name__})"]
assert not missing, f"missing Space-Track secrets: {missing}"
print("releases 1.5.3 and 1.5.4 and the A3 protocol in place; secrets present")

if PART == "eo":
    idx = os.path.join(BATCH, "confirmatory_v153_exploratory_index.csv")
    keep = os.path.join(BATCH, "confirmatory_v153_exploratory_index_20261022_errors.csv")
    if os.path.exists(idx) and not os.path.exists(keep):
        shutil.copy2(idx, keep)
        print(f"failure record kept: {keep}")
    runpy.run_path(os.path.join(CODE, "arf_batch_runner.py"), run_name="__main__")
    r = subprocess.run([sys.executable, os.path.join(CODE, "arf_batch_audit.py"), "--code", CODE, "--batch-dir", BATCH],
                       capture_output=True, text=True)
    print(r.stdout[-6000:], r.stderr[-3000:])
elif PART == "dense":
    # read-only look at the element sets the operational filter counted as malformed (window 1)
    cache = os.path.join(BATCH, "shared_cache", "arf_gp_history_starlink_20260803T000000_7d.json")
    if os.path.exists(cache):
        bad = {}
        for name, l1, l2 in json.load(open(cache))["records"]:
            try:
                int(l1[2:7]); float(l1[20:32]); float("0." + l2[26:33].strip()); float(l2[52:63])
            except (ValueError, TypeError, IndexError):
                kind = ("catalogue number not numeric (Alpha-5)" if isinstance(l1, str) and not l1[2:7].strip().isdigit()
                        else "other")
                bad.setdefault(kind, []).append((name, str(l1)[:24]))
        print("malformed element sets in window 1:", {k: len(v) for k, v in bad.items()},
              "| examples:", {k: v[:3] for k, v in bad.items()})
    R = runpy.run_path(os.path.join(CODE, "arf_batch_runner_v154.py"), run_name="arf_batch_runner_v154")
    D, W = R["DENSE"], R["PROTOCOL_WINDOWS"]
    small = [c for w in W[1:] for c in D if c["analysis_epoch"] == w]
    w1 = [c for c in D if c["analysis_epoch"] == W[0]]
    order = small + [c for c in w1 if c["active_n"] < 8000] + ([c for c in w1 if c["active_n"] == 8000] if INCLUDE_8000 else [])
    ids = [R["run_id"](c) for c in order]
    assert len(set(ids)) == len(ids) and set(ids) <= {R["run_id"](c) for c in D}, "order is not the A3 run list"
    print(f"resuming the dense-constellation arm: {len(order)} of {len(D)} configurations in disconnection-safe "
          "order (completed runs are skipped)" + ("" if INCLUDE_8000 else "; N = 8000 left out (INCLUDE_8000 = False)"))
    R["run_batch"](order, "exploratory_v154_dense")
else:
    raise ValueError("PART must be 'eo' or 'dense'")
