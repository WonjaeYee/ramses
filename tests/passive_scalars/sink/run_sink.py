"""tests/passive_scalars/sink: a main-sequence individual star (sink) blowing a MIST wind
into a uniform, ionised, low-metallicity 3D box.

    python run_sink.py EXE RUNDIR [--tend T] [--levelmin L]

EXE is a 3D build of the 'full_nodust' set with INDI_STAR=1 (build.sh ... full_nodust 3
EXE INDI_STAR=1).  Writes the IC table, ic_sink_indi and the namelist, runs, then checks
(ps_check metrics plus 'neutral': the injected metals are neutral and nothing ionises or
recombines, so the box mass of the neutral stage of every element must change exactly by
the change of the element mass; with ion slot = x*rho the injection adds the TOTAL injected
mass to every element's neutral slot, which is consistent with sum(ions)=rho but gives the
wrong ion masses).
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import ps_cases as C  # noqa: E402
import ps_check  # noqa: E402

MIST = os.environ.get("PS_MIST", "/Users/currodri/Documents/RAMSES_dev/rtzdata_for_public_ramses/"
                      "260518_mist-sample_v2.5_mcut4.0_with-seds.unf")


def write_ic_sink(path, boxlen, mstar_code, mets):
    vals = [1.0, boxlen / 2 + 0.01, boxlen / 2 + 0.01, boxlen / 2 + 0.01] + [0.0] * 17
    # id, msink, x, y, z, v(3), l(3), tform, acc_rate, del_mass, rho_gas, cs2, etherm, vgas(3), mbh, dmfsink
    line = f"{1:10d}" + "".join(f",{v:20.10E}" for v in vals)
    line += f",{6:10d},{mstar_code:20.10E},{0:10d},{0.0:20.10E}" + "".join(f",{m:20.10E}" for m in mets)
    with open(path, "w") as o:
        o.write(" # id,msink,x,y,z,vx,vy,vz,lx,ly,lz,tform,acc_rate,del_mass,rho_gas,cs**2,etherm,"
                "vx_gas,vy_gas,vz_gas,mbh,dmfsink,level,final_mass,evolution_flag,tms,metal1..metal10\n")
        o.write(" # units\n")
        o.write(line + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exe"); ap.add_argument("rundir")
    ap.add_argument("--tend", type=float, default=0.005)
    ap.add_argument("--levelmin", type=int, default=5)
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    lay = C.Layout("full_nodust")
    case = dict(test="sink", set="full_nodust", prof="uniform", conc="sink", slope=1, ndim=3,
                levelmin=a.levelmin, levelmax=a.levelmin + 1, tend=a.tend, bc="periodic", source=True)
    case["id"] = "sink_" + C.case_name(case)
    os.makedirs(a.rundir, exist_ok=True)
    for f in os.listdir(a.rundir):
        if f.startswith("output_"):
            subprocess.run(["rm", "-rf", os.path.join(a.rundir, f)])
    # ambient gas: n ~ 10, T ~ 1e4 K, 0.1 solar, fully ionised (stage 2 of every element)
    N = 2 ** a.levelmin; boxlen = 32.0
    x = (np.arange(N) + 0.5) / N * boxlen
    f = np.tile(C.composition(lay, Z=0.1, ion="ionised", trace=1e-10)[:, None], (1, N))
    # P/rho = (T/mu)/scale_T2 in these units (scale_T2 = m_H/k_B * (l/t)^2)
    scale_T2 = 1.6726e-24 / 1.380649e-16 * (3.0856776e18 / 3003860471423912.0) ** 2
    rho = np.full(N, 10.0); u = np.zeros(N); p = rho * (1e4 / 0.6) / scale_T2
    C.write_ic_table(os.path.join(a.rundir, "ic_table.dat"), lay, x, rho, u, p, f)
    # a 40 Msun star, the metal column of the example ic (element masses, solar-like)
    scale_m = 1.660539e-24 * 3.0856776e18 ** 3 / 1.989e33
    mets = [8.5619993119E+03, 1.1723943295E+03, 1.1723943295E+03, 1.1723943295E+03, 1.8187403636E+01,
            2.2940296605E+01, 2.7629529682E+01, 3.1927709140E+01, 3.6450971786E+01, 1.5109377310E+01]
    write_ic_sink(os.path.join(a.rundir, "ic_sink_indi"), boxlen, 40.0 / scale_m, mets)
    hydro_extra = "passive_cma=.false." if a.label == "cmaoff" else ""
    C.write_namelist(os.path.join(a.rundir, "sink.nml"), lay, os.path.join(HERE, "sink.nml"),
                     levelmin=a.levelmin, levelmax=a.levelmin + 1, tend=a.tend, delta_tout=a.tend / 5,
                     slope=1, riemann="hllc", hydro_extra=hydro_extra, refine="", boundary="nboundary=0",
                     rt_extra="", mist=MIST)
    json.dump(case, open(os.path.join(a.rundir, "case.json"), "w"), indent=1)
    cmd = ["timeout", "3000", a.exe, "sink.nml"]
    if a.np > 1:
        cmd = ["timeout", "3000", "mpirun", "-np", str(a.np), a.exe, "sink.nml"]
    with open(os.path.join(a.rundir, "run.log"), "w") as lg:
        r = subprocess.run(cmd, cwd=a.rundir, stdout=lg, stderr=subprocess.STDOUT,
                           env=dict(os.environ, GFORTRAN_UNBUFFERED_ALL="1"))
    if r.returncode != 0:
        print(f"FAIL {case['id']}: exit {r.returncode}"); sys.exit(1)
    out = ps_check.check(a.rundir, None)
    # neutral-stage mass vs element mass
    snaps = ps_check.load(a.rundir)
    d0, d1 = snaps[0], snaps[-1]
    worst = 0.0; lines = []
    for e, k in lay.elem_idx.items():
        n1 = lay.ion_idx[e][0]
        def tot(d, which):
            rho_ = d["density"]; dv = d["dx"] ** 3
            xe = d[lay.rnames[k]]
            return (rho_ * xe * (d[lay.rnames[n1]] if which == "n" else 1.0) * dv).sum()
        de = tot(d1, "e") - tot(d0, "e"); dn = tot(d1, "n") - tot(d0, "n")
        if abs(de) > 0:
            err = abs(dn - de) / abs(de)
            worst = max(worst, err)
            lines.append(f"    {e:2s}: element mass +{de:.4e}, neutral mass +{dn:.4e}, rel. error {err:.2e}")
    print("\n".join(lines))
    ok = worst < 1e-8
    verdict = "PASS" if ok and out["verdict"] == "PASS" else "FAIL"
    print(f"  {'PASS' if ok else 'FAIL'} neutral   {worst: .3e}  (tol 1.0e-08)")
    out["values"]["neutral"] = worst; out["verdict"] = verdict
    json.dump(out, open(os.path.join(a.rundir, "check.json"), "w"), indent=1)
    print(f"{verdict} {case['id']} [{out['scheme']}] (with neutral)")


if __name__ == "__main__":
    main()
