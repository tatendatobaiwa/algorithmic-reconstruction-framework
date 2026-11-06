"""Repository-relative locations of the result files read by the manuscript scripts.
Set ARF_ROOT to override the repository root."""
import os

def _root():
    if os.environ.get("ARF_ROOT"):
        return os.environ["ARF_ROOT"]
    d = os.path.dirname(os.path.abspath(__file__))
    while not os.path.isdir(os.path.join(d, "protocol")):
        nd = os.path.dirname(d)
        if nd == d:
            raise RuntimeError("repository root not found; set ARF_ROOT")
        d = nd
    return d

ROOT = _root()
D = os.path.join(ROOT, "results", "summary") + os.sep        # aggregated summaries (was ms/data)
A = os.path.join(ROOT, "results", "amendments") + os.sep     # per-amendment results (was arf/)
FACTS = os.path.join(D, "facts.json")                        # number registry of the manuscript
