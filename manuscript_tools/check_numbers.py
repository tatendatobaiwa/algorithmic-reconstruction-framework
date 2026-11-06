"""Every numeric token in a section must be registered (with its source) in facts_methods.json, or be layout.
Usage: check_numbers.py SECTION.tex FACTS.json [FACTS2.json ...]"""
import sys, re, json
t = open(sys.argv[1]).read(); F = {}
for fp in sys.argv[2:]:
    F.update(json.load(open(fp)))
t = re.sub(r"(?m)%.*$", "", t)
t = re.sub(r"\\(cite[pt]?|ref|label|eqref|input)\{[^}]*\}", "", t)
t = re.sub(r"\\begin\{tabular\}\{[^\n]*\}", "", t)                    # column specs
t = re.sub(r"\\(arraystretch|setlength)[^\n]*", "", t)
t = re.sub(r"[0-9a-f]{8}\}", "}", t)                                   # digest prefix in \texttt
LAYOUT = {"1.2", "1.9", "2.6", "2.3", "3.6", "7.9", "256", "00", "1.5.3", "1.5.2", "1.5.0", "1.5.4"}
bad = []
for m in re.finditer(r"(?<![A-Za-z_\\])\d+(?:\{,\}\d+|,\d{3})*(?:\.\d+)*", t):
    n = m.group(0).replace("{,}", ",")
    if n in F or n in LAYOUT or re.fullmatch(r"A\d+", n):
        continue
    if re.fullmatch(r"1\.5\.\d", n):
        continue
    bad.append((n, t[max(0, m.start() - 40):m.end() + 20].replace("\n", " ")))
for n, ctx in bad:
    print(f"UNREGISTERED {n!r}: ...{ctx}...")
print(f"{len(bad)} unregistered numbers"); sys.exit(1 if bad else 0)
