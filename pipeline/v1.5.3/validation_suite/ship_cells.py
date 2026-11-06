"""Split a release file into its 14 cells as separate .py files (verbatim
slices), check that header + cells == the full file byte for byte, and write a
SHA-256 list. Usage: ship_cells.py PIPELINE OUT_DIR"""
import sys, os, re, hashlib

P, OUT = sys.argv[1], sys.argv[2]
NAMES = {"0": "arf_cell0_config.py", "1": "arf_cell1_env_setup.py", "2": "arf_cell2_tle_acquisition.py",
         "3": "arf_cell3_validation.py", "4": "arf_cell4_access_detection.py", "5": "arf_cell5_refinement.py",
         "6": "arf_cell6_selection.py", "6B": "arf_cell6b_elbow_rederivation.py", "7": "arf_cell7_figures.py",
         "7B": "arf_cell7b_coverage_vs_time.py", "7C": "arf_cell7c_onset_robustness.py",
         "7D": "arf_cell7d_acquisition_validation.py", "8": "arf_cell8_hypothesis_validation.py",
         "9": "arf_cell9_manifest.py"}
raw = open(P, "rb").read()
src = raw.decode("utf-8")
hdr = re.compile(r'^# =+\n# CELL ([0-9A-Z]+):', re.M)
ms = list(hdr.finditer(src))
assert [m.group(1) for m in ms] == list(NAMES), [m.group(1) for m in ms]
os.makedirs(OUT, exist_ok=True)
head = src[:ms[0].start()]
parts = [head]
for i, m in enumerate(ms):
    body = src[m.start(): ms[i + 1].start() if i + 1 < len(ms) else len(src)]
    parts.append(body)
    with open(os.path.join(OUT, NAMES[m.group(1)]), "wb") as f:
        f.write(body.encode("utf-8"))
assert "".join(parts).encode("utf-8") == raw, "header + cells differ from the full file"
full = os.path.join(OUT, "algorithmic_reconstruction_framework.py")
with open(full, "wb") as f:
    f.write(raw)
ver = re.search(r'^__version__ = "([^"]+)"', src, re.M).group(1)
for n in NAMES.values():
    s = open(os.path.join(OUT, n), encoding="utf-8").read()
    compile(s, n, "exec")                      # each cell parses on its own
    g = re.search(r'_ARF_CELL_VERSION = "([^"]+)"', s)
    assert (g and g.group(1) == ver) or n == "arf_cell0_config.py", (n, g and g.group(1))
print(f"release {ver}: {len(NAMES)} cells written; header + cells == full file; every cell compiles and carries "
      f"its version guard")
