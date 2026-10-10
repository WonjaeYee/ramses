"""Consistency and conservation check of a RAMSES run (tests/passive_scalars).

    python ps_check.py RUNDIR            (RUNDIR/case.json written by ps_run.py)

Reads every output of RUNDIR, computes the metrics below over time and prints one line
per metric with PASS/FAIL against its tolerance, then an overall verdict; writes
RUNDIR/check.json and RUNDIR/check.png.

Metrics (outputs hold top species as fractions of rho and ion stages as fractions of
their element, in either storage convention):
  top        max |sum(elements, CO, dust, PAH)/rho - 1|
  child      max |sum(ion stages (+H2))/1 - 1| over cells holding the element
  child_abs  max |sum(ion stage densities) - rho_e| / max over the box of rho_e: the same
             error in units of the element's peak density, which is the scale of the
             round-off of a consistent flux-form scheme (where a cell is drained of an
             element, its content falls far below the fluxes that crossed it, and the
             relative error 'child' grows as the content falls)
  slot       max relative change of the box mass of every advected slot (in the x*rho
             convention the ion slot is x*rho, in the x*rho_element one the ion mass)
  ionmass    max relative change of the box mass of every ion stage, x*rho_e (what the
             chemistry sees), over stages holding > 1e-6 of their element
  mass       relative change of the box mass (rho)
  neg        min fraction (positivity)
  L1_top, L1_ion  (periodic advection over whole crossings) L1 error of the top fractions
             and of the ion mass densities with respect to the initial profile
Tolerances: round-off where the scheme guarantees the property; otherwise the measured
bound of the legacy scheme (bounds_legacy.json, if present), else reported only (INFO).
With PS_STRICT=1 every scheme is held to what the consistent scheme guarantees.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ps_cases as C  # noqa: E402

ROUNDOFF = dict(top=1e-12, child=1e-10, child_abs=1e-12, slot=1e-12, ionmass=1e-12,
                mass=1e-12, neg=-1e-14)


def load(rundir):
    vis = os.environ.get("PS_VISU", os.path.join(HERE, "..", "..", "visu"))
    sys.path.insert(0, vis)
    import contextlib
    import io
    import visu_ramses as V
    outs = sorted(glob.glob(os.path.join(rundir, "output_0*")))
    snaps = []
    cwd = os.getcwd()
    os.chdir(rundir)
    try:
        for o in outs:
            n = int(o.split("_")[-1])
            with contextlib.redirect_stdout(io.StringIO()):
                d = V.load_snapshot(n)
            dd = d["data"]
            if "HYDROGEN" not in dd and "metallicity" in dd:   # NMETALS=1 builds (H only)
                dd["HYDROGEN"] = dd["metallicity"]
            snaps.append(dd)
    finally:
        os.chdir(cwd)
    return snaps


def scheme_of(rundir):
    """'cma' / 'legacyB' (fix build with passive_cma off) / 'legacy' (x*rho convention)."""
    log = open(os.path.join(rundir, "run.log"), errors="replace").read()
    if "ion slot = x_ion * rho_element" in log:
        return "cma" if "passive_cma = T" in log else "legacyB"
    return "legacy"


def guaranteed(scheme, case):
    g = {"slot", "mass", "neg"}
    if scheme == "cma":
        g |= {"top", "child_abs", "ionmass"}
        if case.get("drift") and not case.get("counterflux", True):
            g -= {"top"}
    if scheme == "legacyB":
        g |= {"ionmass"}
    if case.get("bc", "periodic") != "periodic":
        g -= {"slot", "mass", "ionmass"}
    if case.get("chem"):
        g -= {"ionmass", "slot"}      # the chemistry changes the ion masses
    if case.get("source"):
        g -= {"slot", "mass", "ionmass"}
    return g


def metrics(lay, snaps, case, conv):
    ndim = case.get("ndim", 1)
    names = lay.rnames
    res = dict(t=[], top=[], child=[], child_abs=[], slot=[], ionmass=[], mass=[], neg=[])
    M0 = None
    for d in snaps:
        rho = d["density"]; dv = d["dx"] ** ndim
        f = np.array([d[nm] for nm in names])           # natural fractions (n, ncell)
        m = C.to_mass(lay, f)                            # mass fractions of every slot
        top = np.abs(f[lay.top].sum(axis=0) - 1).max()
        ch = 0.0; cha = 0.0
        for p, ks in lay.children.items():
            ok = f[p] > 0
            if ok.any():
                ch = max(ch, np.abs(f[ks][:, ok].sum(axis=0) - 1).max())
            ue = rho * f[p]; us = rho * m[ks].sum(axis=0)
            if ue.max() > 0:
                cha = max(cha, np.abs(us - ue).max() / ue.max())
        slotfrac = m.copy()
        if conv == "A":
            for p, ks in lay.children.items():
                slotfrac[ks] = f[ks]                     # x*rho slots
        slot_mass = (rho * slotfrac * dv).sum(axis=1)
        ion_mass = (rho * m * dv).sum(axis=1)
        mass = (rho * dv).sum()
        if M0 is None:
            M0 = (slot_mass, ion_mass, mass)
            el_mass = (rho * m * dv).sum(axis=1)
            sig = np.zeros(lay.n, bool)
            for p, ks in lay.children.items():
                sig[ks] = ion_mass[ks] > 1e-6 * el_mass[p]
            keep = np.abs(slot_mass) > 1e-14 * np.abs(slot_mass).sum()
        res["t"].append(d["time"])
        res["top"].append(top); res["child"].append(ch); res["child_abs"].append(cha)
        res["slot"].append(np.max(np.abs(slot_mass[keep] - M0[0][keep]) / np.abs(M0[0][keep])))
        res["ionmass"].append(np.max(np.abs(ion_mass[sig] - M0[1][sig]) / M0[1][sig]) if sig.any() else 0.0)
        res["mass"].append(abs(mass - M0[2]) / M0[2])
        res["neg"].append(min(f[lay.top].min(), f[lay.parent >= 0].min() if (lay.parent >= 0).any() else 0))
    # accuracy over whole crossings of a periodic box (uniform grid only)
    if case.get("bc", "periodic") == "periodic" and case.get("crossings") and not case.get("drift") \
            and not case.get("source") and not case.get("chem") and case.get("levelmin") == case.get("levelmax"):
        d0, d1 = snaps[0], snaps[-1]
        o0, o1 = np.argsort(d0["x"]), np.argsort(d1["x"])
        f0 = np.array([d0[nm] for nm in names])[:, o0]; f1 = np.array([d1[nm] for nm in names])[:, o1]
        m0, m1 = C.to_mass(lay, f0), C.to_mass(lay, f1)
        e_top = np.mean([np.abs(f1[k] - f0[k]).sum() / np.abs(f0[k]).sum() for k in lay.top if f0[k].sum() > 0])
        ions = [k for k in range(lay.n) if lay.kind[k] == "ion" and sig[k]]
        e_ion = np.mean([np.abs(m1[k] - m0[k]).sum() / np.abs(m0[k]).sum() for k in ions]) if ions else 0.0
        res["L1_top"] = [float(e_top)]; res["L1_ion"] = [float(e_ion)]
    return res


def plot(rundir, lay, snaps, res, case, verdicts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    d = snaps[-1]
    o = np.argsort(d["x"])
    f = np.array([d[nm] for nm in lay.rnames])[:, o]
    x = d["x"][o]
    ax[0].semilogy(x, np.abs(f[lay.top].sum(axis=0) - 1) + 1e-17, label="|sum(top)/rho - 1|")
    ce = np.zeros_like(x)
    for p, ks in lay.children.items():
        okp = f[p] > 0
        ce = np.maximum(ce, np.where(okp, np.abs(f[ks].sum(axis=0) - 1), 0))
    ax[0].semilogy(x, ce + 1e-17, label="max_e |sum(ions)/rho_e - 1|")
    ax[0].set_xlabel("x"); ax[0].set_title(f"consistency at t={d['time']:.3g}"); ax[0].legend(fontsize=8)
    ax[0].set_ylim(1e-17, 1)
    t = np.array(res["t"])
    for k in ("top", "child", "slot", "ionmass", "mass"):
        ax[1].semilogy(t, np.array(res[k]) + 1e-17, marker="o", ms=3, label=k)
    ax[1].set_xlabel("t"); ax[1].set_title("max error vs time"); ax[1].legend(fontsize=8)
    ax[1].set_ylim(1e-17, 10)
    d0 = snaps[0]; o0 = np.argsort(d0["x"])
    for k in list(lay.top[:3]) + [ks[min(1, len(ks) - 1)] for ks in list(lay.children.values())[:2]]:
        ax[2].semilogy(d0["x"][o0], d0[lay.rnames[k]][o0], ":", lw=1)
        ax[2].semilogy(x, f[k], "-", lw=1, label=lay.rnames[k])
    ax[2].set_xlabel("x"); ax[2].set_title("fractions: t=0 (dotted), final (solid)"); ax[2].legend(fontsize=7)
    v = "PASS" if all(s == "PASS" or s == "INFO" for s in verdicts.values()) else "FAIL"
    fig.suptitle(f"{case.get('test')} {C.case_name(case)} [{scheme_of(rundir)}] -> {v}", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(rundir, "check.png"), dpi=90)
    plt.close(fig)


def check(rundir, bounds_file=None, do_plot=True):
    case = json.load(open(os.path.join(rundir, "case.json")))
    lay = C.Layout(case["set"])
    snaps = load(rundir)
    if len(snaps) < 2:
        print(f"FAIL {rundir}: {len(snaps)} outputs"); return dict(verdict="FAIL", error="no outputs")
    sch = scheme_of(rundir)
    conv = "A" if sch == "legacy" else "B"
    res = metrics(lay, snaps, case, conv)
    g = guaranteed(sch, case)
    if os.environ.get("PS_STRICT"):     # hold any scheme to what the consistent scheme guarantees
        g = guaranteed("cma", case)
    bounds = {}
    if bounds_file and os.path.exists(bounds_file):
        bounds = json.load(open(bounds_file)).get(case["id"], {})
    verdicts = {}; vals = {}
    for k in ("top", "child", "child_abs", "slot", "ionmass", "mass", "neg", "L1_top", "L1_ion"):
        if k not in res:
            continue
        v = min(res[k]) if k == "neg" else max(res[k])
        vals[k] = float(v)
        if k in g:
            ok = v >= ROUNDOFF[k] if k == "neg" else v <= ROUNDOFF[k]
            verdicts[k] = "PASS" if ok else "FAIL"
            tol = ROUNDOFF[k]
        elif sch == "legacy" and k in bounds:
            tol = 2.0 * bounds[k] + 1e-15
            verdicts[k] = "PASS" if v <= tol else "FAIL"
        else:
            verdicts[k] = "INFO"; tol = None
        print(f"  {verdicts[k]:4s} {k:9s} {v: .3e}" + (f"  (tol {tol:.1e})" if tol is not None else ""))
    overall = "FAIL" if "FAIL" in verdicts.values() else "PASS"
    out = dict(id=case["id"], case=case, scheme=sch, values=vals, verdicts=verdicts, verdict=overall)
    json.dump(out, open(os.path.join(rundir, "check.json"), "w"), indent=1)
    if do_plot:
        try:
            plot(rundir, lay, snaps, res, case, verdicts)
        except Exception as e:  # plotting must never decide the verdict
            print("  (plot failed:", e, ")")
    print(f"{overall} {case['id']} [{sch}]")
    return out


if __name__ == "__main__":
    check(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
