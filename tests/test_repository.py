"""Repository checks that need no network, no Orekit and no licensed data.

Run from the repository root:  python -m pytest -q tests
  1. The release of record is the file the protocol was frozen against (SHA-256).
  2. The protocol carries amendments A1 to A13, in order.
  3. Every Python file in the repository compiles.
  4. The manuscript tables regenerate from results/ and equal the reference copies.
  5. The result figures and the A12/A13 figure regenerate (each script asserts its
     inputs against the stored summaries before drawing).
  6. No file holds a Space-Track element set or a WDPA/WDPCA boundary file.
"""
import hashlib, json, os, py_compile, re, subprocess, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "manuscript_tools"
REF = ROOT / "results" / "reference_tables"
RELEASE = ROOT / "pipeline" / "v1.5.3" / "algorithmic_reconstruction_framework.py"
PROTOCOL = ROOT / "protocol" / "arf_confirmatory_protocol.json"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(args, cwd):
    env = dict(os.environ, ARF_ROOT=str(ROOT), MPLBACKEND="Agg")
    r = subprocess.run([sys.executable, *map(str, args)], cwd=cwd, env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    return r


def test_release_of_record_matches_protocol():
    p = json.loads(PROTOCOL.read_text())
    assert sha256(RELEASE) == p["pipeline_sha256"]
    assert sha256(RELEASE).startswith("14426c89")


def test_protocol_amendments_in_order():
    ids = [a["id"] for a in json.loads(PROTOCOL.read_text())["amendments"]]
    assert ids == [f"A{i}" for i in range(1, 14)]


def test_every_python_file_compiles(tmp_path):
    files = [f for f in ROOT.rglob("*.py") if ".git" not in f.parts]
    assert len(files) > 50
    for f in files:
        py_compile.compile(str(f), cfile=str(tmp_path / "x.pyc"), doraise=True)


def test_tables_regenerate(tmp_path):
    (tmp_path / "tables").mkdir()
    run([TOOLS / "make_tables.py"], cwd=tmp_path)
    for name in ("tab_confirm.tex", "tab_benchmark.tex", "tab_a4_exact_economy.tex",
                 "tab_a9_budget.tex", "tab_a9_failsignal.tex", "tab_a8_pooled.tex",
                 "tab_a11_volume.tex"):
        assert (tmp_path / "tables" / name).read_text() == (REF / name).read_text(), name


def test_a12_a13_figure_and_table(tmp_path):
    summ = ROOT / "results" / "summary" / "a12_a13_summary.json"
    run([TOOLS / "figscripts" / "fig_a12_a13.py", summ, tmp_path / "fig.png", tmp_path / "tab_a13.tex"],
        cwd=tmp_path)
    assert (tmp_path / "tab_a13.tex").read_text() == (REF / "tab_a13.tex").read_text()
    assert (tmp_path / "fig.png").stat().st_size > 10_000


def test_result_figures(tmp_path):
    run([TOOLS / "figscripts" / "fig_results.py", str(tmp_path) + "/"], cwd=tmp_path)
    assert len(list(tmp_path.glob("*.png"))) >= 5


TLE_LINE = re.compile(r"^1 \d{5}[A-Z ] \d{5}[A-Z ]{0,3}\s+\d{5}\.\d+", re.M)
FORBIDDEN_EXT = {".tle", ".3le", ".shp", ".shx", ".dbf", ".gpkg", ".geojson", ".kml", ".kmz"}


def test_no_licensed_data():
    for f in ROOT.rglob("*"):
        if ".git" in f.parts or not f.is_file():
            continue
        assert f.suffix.lower() not in FORBIDDEN_EXT, f
        assert "wdpa" not in f.name.lower() and "wdpca" not in f.name.lower(), f
        if f.suffix.lower() in {".txt", ".csv", ".json", ".py", ".log", ".md"}:
            assert not TLE_LINE.search(f.read_text(errors="ignore")), f
