"""Run a matrix of advection-accuracy cases and collect the errors.

    acc_run.py --matrix full|nightly --out OUT --builds BUILDS [--tree TREE] [-j 8] [--keep] [--only S1,S2]

BUILDS/<solver>_<ndim>d/ramses_acc<ndim>d are the executables; with --tree, the missing ones are built from
TREE (acc_build.sh).  Every case runs in OUT/runs/<case>/ (snapshots deleted after the check unless --keep;
acc.json, profile.npz, run.log stay).  OUT/results.json holds every case's acc.json, plus the wall time.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import acc_cases as A  # noqa: E402
import acc_check  # noqa: E402


def exe_for(builds, solver, ndim, tree):
    exe = os.path.join(os.path.abspath(builds), f"{solver}_{ndim}d", f"ramses_acc{ndim}d")
    if not os.path.exists(exe):
        if not tree:
            raise SystemExit(f"missing {exe} (pass --tree to build it)")
        r = subprocess.run([os.path.join(HERE, "acc_build.sh"), tree, solver, str(ndim),
                            os.path.dirname(exe)], capture_output=True, text=True,
                           env=dict(os.environ, PS_PYTHON=sys.executable))
        if r.returncode != 0 or not os.path.exists(exe):
            raise SystemExit(f"build {solver} {ndim}d failed: {r.stderr.strip()}")
    return exe


def one(c, exe, out, keep, timeout):
    d = os.path.join(out, "runs", A.case_name(c))
    shutil.rmtree(d, ignore_errors=True)
    A.prepare(c, d)
    t0 = time.time()
    try:
        with open(os.path.join(d, "run.log"), "w") as log:
            subprocess.run([exe, "acc.nml"], cwd=d, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    except subprocess.TimeoutExpired:
        pass
    wall = time.time() - t0
    try:
        r = acc_check.check(d)
    except Exception as e:                       # noqa: BLE001 - a broken run is a result, not a crash
        r = dict(case=A.case_name(c), status="FAIL", error=f"check: {e!r}"[:300])
    r["wall_s"] = wall
    json.dump(r, open(os.path.join(d, "acc.json"), "w"), indent=1)
    if not keep:
        for o in os.listdir(d):
            if o.startswith("output_"):
                shutil.rmtree(os.path.join(d, o), ignore_errors=True)
        os.remove(os.path.join(d, "ic_table.dat"))
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", default="nightly")
    ap.add_argument("--out", required=True)
    ap.add_argument("--builds", required=True)
    ap.add_argument("--tree")
    ap.add_argument("-j", type=int, default=4)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--timeout", type=float, default=600)
    a = ap.parse_args()
    only = a.only.split(",")
    cases = [c for c in A.matrix(a.matrix) if any(o in A.case_name(c) for o in only)]
    exes = {}
    for c in cases:
        key = (c["solver"], c["ndim"])
        if key not in exes:
            exes[key] = exe_for(a.builds, *key, a.tree)
    a.out = os.path.abspath(a.out)
    os.makedirs(a.out, exist_ok=True)
    # longest first: low Mach and large N take the most steps
    cases.sort(key=lambda c: -(1 + 1 / c["mach"]) * c["N"] ** (c["ndim"] + 1) * c["periods"])
    t0 = time.time()
    with ProcessPoolExecutor(a.j) as ex:          # processes: the snapshot reader changes directory
        futs = [ex.submit(one, c, exes[(c["solver"], c["ndim"])], a.out, a.keep, a.timeout) for c in cases]
        res = [f.result() for f in futs]
    res.sort(key=lambda r: r["case"])
    json.dump(dict(matrix=a.matrix, wall_s=time.time() - t0, ncases=len(res), results=res),
              open(os.path.join(a.out, "results.json"), "w"), indent=1)
    bad = [r for r in res if r["status"] != "OK"]
    print(f"{len(res)} cases in {time.time() - t0:.0f} s, {len(bad)} failed")
    for r in bad:
        print(f"  FAIL {r['case']}: {r.get('error', '')}")


if __name__ == "__main__":
    main()
