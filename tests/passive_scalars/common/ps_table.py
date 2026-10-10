"""Summary table of the passive-scalar test matrix: worst value of every metric per
(test, species set, profile) over the concentrations and slope limiters.

    python ps_table.py RUNS LABEL [LABEL2 ...]
"""
import collections
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ps_check  # noqa: E402

KEYS = ("top", "child_abs", "child", "ionmass", "slot")


def rows(runs, label):
    agg = collections.defaultdict(lambda: collections.defaultdict(float))
    verd = collections.defaultdict(lambda: [0, 0])
    global strict
    strict = collections.defaultdict(lambda: [0, 0])
    for f in glob.glob(os.path.join(runs, label, "*", "*", "check.json")):
        r = json.load(open(f))
        c = r["case"]
        prof = c["prof"] + "".join(f"_{k}{c[k]}" for k in ("rho_contrast", "w", "m") if k in c)
        key = (c["test"], c["set"], prof)
        for k in KEYS:
            if k in r["values"]:
                agg[key][k] = max(agg[key][k], r["values"][k])
        verd[key][0] += r["verdict"] == "PASS"; verd[key][1] += 1
        # strict: held to what the consistent scheme guarantees, at round-off
        g = ps_check.guaranteed("cma", c)
        bad = any((r["values"][k] < ps_check.ROUNDOFF[k]) if k == "neg" else (r["values"][k] > ps_check.ROUNDOFF[k])
                  for k in g if k in r["values"])
        strict[key][0] += bad; strict[key][1] += 1
    return agg, verd


def main():
    runs = sys.argv[1]
    for label in sys.argv[2:]:
        agg, verd = rows(runs, label)
        if not agg:
            continue
        print(f"\n[{label}] worst value per (test, set, profile); PASS/cases; strict = cases failing"
              " the round-off guarantees of the consistent scheme")
        print(f"{'test':7s} {'set':9s} {'profile':22s}" + "".join(f"{k:>11s}" for k in KEYS) + "   pass  strict")
        for key in sorted(agg):
            v = agg[key]
            print(f"{key[0]:7s} {key[1]:9s} {key[2]:22s}" + "".join(f"{v[k]:11.1e}" for k in KEYS)
                  + f"   {verd[key][0]}/{verd[key][1]}   {strict[key][0]}/{strict[key][1]} fail")


if __name__ == "__main__":
    main()
