"""Sweep the passive-scalar test matrix with the Python replica of the RAMSES hydro step.

    python run_replica.py [--quick] [--out DIR] [--schemes legacy-A,cma-B,...]

For every case and scheme it records the consistency and conservation metrics and writes
DIR/replica_results.json and a PASS/FAIL table (DIR/replica_summary.txt).  Schemes are
<scheme>-<convention>: legacy-A is RAMSES-RTZ as it is, cma-B the recommended scheme.
Tolerances (round-off) apply to the properties a scheme guarantees:
  slot masses      : every scheme (flux form)                       1e-12
  sum(top)=rho     : cma, cma_y, renorm                             1e-12
  sum(ions)=elem.  : cma, cma_y, renorm   1e-13 in units of the element's peak density
                     (child_abs; the error relative to the local element density, 'child',
                     grows where a cell is drained of the element and is reported only)
  ion mass y*rho_e : conv B (the slot is the ion mass)              1e-12
  positivity       : all (min fraction >= -1e-14)
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import ps1d as P  # noqa: E402
import ps_cases as C  # noqa: E402

TOL = dict(slot=1e-12, top=1e-12, child=1e-10, child_abs=1e-13, ionmass=1e-12, neg=-1e-14)


def guarantees(scheme, conv):
    g = {"slot", "neg"}
    if scheme in ("cma", "cma_y", "renorm"):
        g |= {"top", "child_abs"}
    if conv == "B" and scheme != "renorm":
        g |= {"ionmass"}
    return g


def metrics(lay, hist, conv, sim, f0=None):
    d0 = hist[0][1]
    top = max(h[1]["top_err"] for h in hist)
    child = max(h[1]["child_err"] for h in hist)
    child_abs = max(h[1]["child_abs"] for h in hist)
    m0 = d0["spec_mass"]; tot = abs(m0).sum()
    keep = np.abs(m0) > 1e-14 * tot
    slot = max(np.max(np.abs(h[1]["spec_mass"][keep] - m0[keep]) / np.abs(m0[keep])) for h in hist)
    im0 = d0["ion_mass"]
    sig = np.zeros(lay.n, bool)
    for p, ks in lay.children.items():
        sig[ks] = im0[ks] > 1e-6 * m0[p]
    allk = im0 > 0
    ion_sig = max(np.max(np.abs(h[1]["ion_mass"][sig] - im0[sig]) / im0[sig]) for h in hist) if sig.any() else 0.0
    ion_all = max(np.max(np.abs(h[1]["ion_mass"][allk] - im0[allk]) / im0[allk]) for h in hist) if allk.any() else 0.0
    neg = min(h[1]["min_frac"] for h in hist)
    out = dict(top=top, child=child, child_abs=child_abs, slot=slot, ionmass=ion_sig, ionmass_all=ion_all, neg=neg,
               nstep=sim.nstep, renorm_cells=sim.counters["renorm_cells"])
    if f0 is not None:   # periodic advection over one crossing: compare with the initial state
        f = hist[-1][1]["fnat"]
        tp = lay.top
        e_top = np.mean([np.abs(f[k] - f0[k]).sum() / np.abs(f0[k]).sum() for k in tp if f0[k].sum() > 0])
        # ion MASS densities (y * element fraction): the ion fraction of a mixed cell is
        # dominated by the metal-rich parcel, so y itself is not compared with the
        # unmixed exact solution
        mi, mi0 = C.to_mass(lay, f), C.to_mass(lay, f0)
        ions = [k for k in range(lay.n) if lay.kind[k] == "ion" and sig[k]]
        e_ion = np.mean([np.abs(mi[k] - mi0[k]).sum() / np.abs(mi0[k]).sum() for k in ions]) if ions else 0.0
        out.update(L1_top=float(e_top), L1_ion=float(e_ion))
    return out


def run_case(case, scheme, conv, N=128):
    lay = C.Layout(case["set"])
    x = (np.arange(N) + 0.5) / N
    kw = {k: case[k] for k in ("rho_contrast", "w", "m") if k in case}
    rho, u, p, f = C.profile(lay, case["prof"], x, conc=case["conc"], **kw)
    bc = "outflow" if case["prof"] == "sod3" else "periodic"
    drift = None
    if case.get("drift"):
        wv = case["drift"]
        drift = lambda xx: wv * (0.5 + 0.5 * np.sin(2 * np.pi * xx))  # noqa: E731
    sim = P.Sim(lay, x, rho, u, p, f.copy(), conv=conv, scheme=scheme, slope_type=case["slope"],
                riemann=case.get("riemann", "hllc"), bc=bc, drift=drift,
                counterflux=case.get("counterflux", False), gamma=1.4)
    tend = case.get("tend", 1.0)
    hist = sim.run(tend, nout=10)
    f0 = f if (bc == "periodic" and drift is None and abs(tend - 1.0) < 1e-12) else None
    return metrics(lay, hist, conv, sim, f0)


def all_cases(quick):
    cases = C.advect_cases(quick=quick)
    for s in ["HHeCO", "full"]:
        for conc in ["solar", "extreme"]:
            cases.append(dict(test="sod", set=s, prof="sod3", conc=conc, slope=1, tend=0.2, u0=0.0))
            cases.append(dict(test="sod", set=s, prof="sod3", conc=conc, slope=2, tend=0.2, u0=0.0))
            cases.append(dict(test="advect_llf", set=s, prof="step3", conc=conc, slope=1, riemann="llf"))
            for cf in (False, True):
                cases.append(dict(test="tva", set=s, prof="smooth", m=1, conc=conc, slope=1,
                                  drift=0.3, counterflux=cf, tend=0.5))
    return cases


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=os.environ.get("PS_OUT", "."))
    ap.add_argument("--schemes", default="legacy-A,legacy-B,cma-A,cma-B,cma_y-B,renorm-A")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    schemes = [s.split("-") for s in a.schemes.split(",")]
    res = []
    t0 = time.time()
    for case in all_cases(a.quick):
        for sch, conv in schemes:
            m = run_case(case, sch, conv)
            g = guarantees(sch, conv)
            fails = [k for k in g if (m[k] < TOL[k] if k == "neg" else m[k] > TOL[k])]
            res.append(dict(case=case, scheme=sch, conv=conv, metrics=m, fails=fails))
    json.dump(res, open(os.path.join(a.out, "replica_results.json"), "w"), indent=1, default=float)
    lines = []
    for r in res:
        m = r["metrics"]
        lines.append(f"{'PASS' if not r['fails'] else 'FAIL'} {r['scheme']+'-'+r['conv']:9s} "
                     f"{C.case_name({'test': 0, **r['case']}):70s} top {m['top']:.1e} child {m['child']:.1e} ({m['child_abs']:.1e}) "
                     f"slot {m['slot']:.1e} ionmass {m['ionmass']:.1e} (all {m['ionmass_all']:.1e}) "
                     f"min {m['neg']:.1e}" + (f" L1top {m['L1_top']:.2e} L1ion {m['L1_ion']:.2e}" if 'L1_top' in m else "")
                     + (f"  [{','.join(sorted(r['fails']))}]" if r["fails"] else ""))
    open(os.path.join(a.out, "replica_summary.txt"), "w").write("\n".join(lines) + "\n")
    npass = sum(1 for r in res if not r["fails"])
    print("\n".join(lines))
    print(f"replica: {npass}/{len(res)} PASS ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
