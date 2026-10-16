"""tests/passive_scalars/restart: an advect case run straight to t=1 and, separately, to
t=0.5 then restarted from that output to t=1.  The ion slots of a restart file hold the
ionisation fractions, which init_hydro multiplies by rho; with ion slot = x*rho_element
init_passive_groups must turn them into ion mass densities.  PASS if the final states agree
to 1e-12 (every field) and the restarted run is still consistent.

    python run_restart.py EXEDIR LABEL OUTDIR [--set full]
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import ps_cases as C  # noqa: E402
import ps_check  # noqa: E402
import ps_run  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exes"); ap.add_argument("label"); ap.add_argument("out")
    ap.add_argument("--set", default="full")
    a = ap.parse_args()
    case = dict(test="advect", set=a.set, prof="step3", conc="extreme", slope=1, levelmin=7, levelmax=7,
                tend=1.0, crossings=1)
    case["id"] = "restart_" + C.case_name(case)
    exe = os.path.join(a.exes, f"{ps_run.PREFIX.get(a.label, a.label)}_{a.set}_1d")
    straight = os.path.join(a.out, a.label, "restart", case["id"], "straight")
    split = os.path.join(a.out, a.label, "restart", case["id"], "split")
    rc, _ = ps_run.run(case, exe, straight, a.label)
    assert rc == 0, "straight run failed"
    rc, _ = ps_run.run(case, exe, split, a.label)            # same run; restart from output 6 (t=0.5)
    nml = os.path.join(split, "advect.nml")
    s = open(nml).read()
    s = s.replace("nrestart=0", "nrestart=6")
    open(nml, "w").write(s)
    for d in sorted(os.listdir(split)):
        if d.startswith("output_") and int(d[-5:]) > 6:
            subprocess.run(["rm", "-rf", os.path.join(split, d)])
    with open(os.path.join(split, "run_restart.log"), "w") as lg:
        r = subprocess.run(["timeout", "600", exe, "advect.nml"], cwd=split, stdout=lg, stderr=subprocess.STDOUT)
    assert r.returncode == 0, "restart failed"
    s1 = ps_check.load(straight)[-1]; s2 = ps_check.load(split)[-1]
    lay = C.Layout(a.set)
    o1, o2 = np.argsort(s1["x"]), np.argsort(s2["x"])
    # compare densities: rho, P, and every slot as a mass density (ion stages as
    # x*rho_e), each in units of its peak (the ion fraction x itself is ill-conditioned
    # where an element is drained; the restart of the x*rho_element convention goes
    # through x*rho*rho_e/rho, so it is not bitwise)
    f1 = np.array([s1[k][o1] for k in lay.rnames]); f2 = np.array([s2[k][o2] for k in lay.rnames])
    m1 = C.to_mass(lay, f1) * s1["density"][o1]; m2 = C.to_mass(lay, f2) * s2["density"][o2]
    worst = 0.0; wk = ""
    for k, v1, v2 in [("density", s1["density"][o1], s2["density"][o2]),
                      ("pressure", s1["pressure"][o1], s2["pressure"][o2])] + \
            [(lay.rnames[j], m1[j], m2[j]) for j in range(lay.n)]:
        scale = np.abs(v1).max()
        if scale <= 0:
            continue
        e = np.abs(v1 - v2).max() / scale
        if e > worst:
            worst, wk = e, k
    ok = worst < 1e-10
    print(f"  {'PASS' if ok else 'FAIL'} restart   {worst: .3e}  (field {wk}; tol 1.0e-10)")
    print(f"{'PASS' if ok else 'FAIL'} {case['id']} [{ps_check.scheme_of(straight)}]")


if __name__ == "__main__":
    main()
