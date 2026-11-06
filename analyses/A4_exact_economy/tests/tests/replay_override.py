"""Replay Cells 6-8 of a release on saved real inputs (the parity procedure of
arf_colab_parity.replay), with Cell 0 settings overridden, to test sensitivity
settings and measure their cost. Usage: replay_override.py PIPELINE PARITY REFINED TILES ORBIT_MIN ROI FLEET CASE NSATS JSON_OVERRIDES"""
import sys, json, time, importlib.util, inspect, textwrap, resource
P, PAR, REF, TILES, ORB, ROI, FLEET, CASE, NS, OV = sys.argv[1:]
spec = importlib.util.spec_from_file_location("par", PAR); par = importlib.util.module_from_spec(spec); spec.loader.exec_module(par)
src = inspect.getsource(par.replay)
src = src.replace('    ns["h3"] = h3', '    ns["h3"] = h3\n    ns.update(OVERRIDES)').replace(
    "def replay(pipeline_path,", "def replay_ov(OVERRIDES, pipeline_path,")
g = dict(vars(par)); exec(src, g)
import pandas as pd
tiles = pd.read_csv(TILES).iloc[:, 0].tolist()
ov = json.loads(OV)
t0 = time.time()
ns = g["replay_ov"](ov, P, REF, tiles, float(ORB), ROI, FLEET, CASE, int(NS))
print(json.dumps({"overrides": ov, "seconds": round(time.time() - t0, 1), "out_dir": ns["_out_dir"],
                  "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024}))
