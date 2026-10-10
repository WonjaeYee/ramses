"""Record the values measured with the legacy scheme as the bounds it is held to.

    python record_bounds.py RUNS/<label> [OUT=../bounds_legacy.json]

Reads every RUNS/<label>/<test>/<case>/check.json of a legacy (ion slot = x*rho) run and
stores its metric values per case id; ps_check then fails a legacy run whose error grows
beyond twice the recorded value (a regression guard for the current scheme).
"""
import glob
import json
import os
import sys

src = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                                          "bounds_legacy.json")
b = {}
for f in sorted(glob.glob(os.path.join(src, "*", "*", "check.json"))):
    r = json.load(open(f))
    if r.get("scheme") != "legacy":
        continue
    b[r["id"]] = {k: float(f"{v:.3e}") for k, v in r["values"].items() if k != "neg"}
json.dump(b, open(out, "w"), indent=0, sort_keys=True)
print(f"{len(b)} cases -> {out}")
