"""Run the cases of one test of tests/passive_scalars with one executable family.

    python ps_run.py --test advect --exes EXEDIR --label base --out OUTDIR
                     [--quick] [--match SUBSTR] [--bounds FILE] [--record-bounds FILE]

For every case of the test (ps_cases.TESTS[test]) it writes OUTDIR/<label>/<test>/<id>/
(case.json, ic_table.dat, <test>.nml), runs EXEDIR/<label>_<set>_<ndim>d there and checks
the outputs (ps_check.py).  Prints one PASS/FAIL line per case and writes
OUTDIR/<label>/<test>/summary.json.  --record-bounds stores the measured values as the
bounds the legacy scheme is held to (bounds_legacy.json).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ps_cases as C  # noqa: E402
import ps_check  # noqa: E402

TESTDIR = os.path.join(HERE, "..")


# a label is the prefix of the executables (EXEDIR/<prefix>_<set>_<ndim>d), except
# 'cmaoff': the fixed build with passive_cma=.false. (new storage convention, legacy
# reconstruction), to separate the effect of the convention from that of CMA
PREFIX = {"cmaoff": "cma"}
EXTRA = {"cmaoff": "passive_cma=.false."}


def subs_for(case, lay):
    test = case["test"]
    ncell = 2 ** case["levelmax"]
    sub = dict(levelmin=case["levelmin"], levelmax=case["levelmax"], tend=case["tend"],
               delta_tout=case["tend"] / 10.0, slope=case["slope"], riemann=case.get("riemann", "hllc"),
               crossings=case.get("crossings", 0), rtz_cooling=".true." if case.get("chem") else ".false.",
               hydro_extra="", refine="", boundary="nboundary=0", calima_extra="", rt_extra="")
    if case.get("chem"):
        sub["rt_extra"] = "rtz_include_dust=.true.\nrtz_include_collisional_ionization=.true."
    if case.get("bc") == "outflow":
        sub["boundary"] = "nboundary=2\nibound_min=-1,+1\nibound_max=-1,+1\nbound_type=2,2"
    if test == "amr":
        # refine on the density jump and on the H mass fraction (composition)
        # (passive index k+1 of slot k: the H mass fraction and the HI stage)
        kh1 = lay.ion_idx["H"][0] + 1
        sub["refine"] = ("err_grad_d=0.05\nerr_grad_var(1)=1d-3\nerr_grad_floor(1)=1d-3\n"
                         f"err_grad_var({kh1})=0.02\nerr_grad_floor({kh1})=1d-3\n"
                         "interpol_var=1\ninterpol_type=2")
    if case.get("calima_extra"):
        sub["calima_extra"] = case["calima_extra"]
    if case.get("drift"):
        sub["calima_extra"] = (f"dust_tva=.true.\nuse_w_drift_test=.true.\nw_drift_test={case['drift']},0.0,0.0\n")
    return sub, ncell


def prepare(case, rundir, label=""):
    lay = C.Layout(case["set"])
    os.makedirs(rundir, exist_ok=True)
    for f in os.listdir(rundir):
        if f.startswith("output_"):
            shutil.rmtree(os.path.join(rundir, f))
    os.makedirs(os.path.join(rundir, "SEDtables"), exist_ok=True)
    sub, ncell = subs_for(case, lay)
    sub["hydro_extra"] = EXTRA.get(label, "")
    x = (np.arange(ncell) + 0.5) / ncell
    kw = {k: case[k] for k in ("rho_contrast", "w", "m", "u0") if k in case}
    rho, u, p, f = C.profile(lay, case["prof"], x, conc=case["conc"], **kw)
    C.write_ic_table(os.path.join(rundir, "ic_table.dat"), lay, x, rho, u, p, f)
    nml = os.path.join(rundir, case["test"] + ".nml")
    C.write_namelist(nml, lay, os.path.join(TESTDIR, case["test"], case["test"] + ".nml"), **sub)
    json.dump(case, open(os.path.join(rundir, "case.json"), "w"), indent=1)
    return nml


def run(case, exe, rundir, label="", timeout=600):
    nml = prepare(case, rundir, label)
    env = dict(os.environ, GFORTRAN_UNBUFFERED_ALL="1")
    t0 = time.time()
    with open(os.path.join(rundir, "run.log"), "w") as lg:
        r = subprocess.run(["timeout", str(timeout), exe, os.path.basename(nml)], cwd=rundir, stdout=lg,
                           stderr=subprocess.STDOUT, env=env)
    return r.returncode, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True, choices=sorted(C.TESTS))
    ap.add_argument("--exes", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--match", default="")
    ap.add_argument("--bounds", default=os.path.join(TESTDIR, "bounds_legacy.json"))
    ap.add_argument("--record-bounds", default=None)
    ap.add_argument("--noplot", action="store_true")
    a = ap.parse_args()
    cases = [c for c in C.TESTS[a.test](quick=a.quick) if a.match in c["id"]]
    results = []
    for case in cases:
        exe = os.path.join(a.exes, f"{PREFIX.get(a.label, a.label)}_{case['set']}_{case['ndim']}d")
        rundir = os.path.join(a.out, a.label, a.test, case["id"])
        if not os.path.exists(exe):
            print(f"SKIP {case['id']}: no {exe}"); continue
        rc, dt = run(case, exe, rundir, a.label)
        if rc != 0:
            tail = open(os.path.join(rundir, "run.log"), errors="replace").read().splitlines()[-3:]
            print(f"FAIL {case['id']}: exit {rc} ({' | '.join(tail)})")
            results.append(dict(id=case["id"], verdict="FAIL", error=f"exit {rc}")); continue
        print(f"--- {case['id']} ({dt:.1f} s)")
        try:
            results.append(ps_check.check(rundir, a.bounds, do_plot=not a.noplot))
        except Exception as e:
            print(f"FAIL {case['id']}: check raised {type(e).__name__}: {e}")
            results.append(dict(id=case["id"], verdict="FAIL", error=str(e)))
    os.makedirs(os.path.join(a.out, a.label, a.test), exist_ok=True)
    json.dump(results, open(os.path.join(a.out, a.label, a.test, "summary.json"), "w"), indent=1)
    if a.record_bounds:
        b = json.load(open(a.record_bounds)) if os.path.exists(a.record_bounds) else {}
        for r in results:
            if "values" in r:
                b[r["id"]] = {k: v for k, v in r["values"].items() if k != "neg"}
        json.dump(b, open(a.record_bounds, "w"), indent=1, sort_keys=True)
    npass = sum(r["verdict"] == "PASS" for r in results)
    print(f"== {a.test} [{a.label}]: {npass}/{len(results)} PASS")


if __name__ == "__main__":
    main()
