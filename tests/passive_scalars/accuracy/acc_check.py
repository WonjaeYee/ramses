"""Errors of one advection-accuracy run against the exact solution: writes acc.json in the run directory.

    acc_check.py <rundir> [<rundir> ...]

Per scalar (top-hat, Gaussian, sine), with q the cell values and q_e the exact cell averages at the
snapshot's own time (RAMSES overshoots tout by up to one step):
  L1, L2, Linf          mean |q - q_e|, rms, max over the cells
  mass                  |sum rho q (end) - sum rho q (start)| / sum rho q (start)  (flux form: round-off)
  over, under           max(q) - max(q_e) and min(q_e) - min(q) where positive (new extrema), qmin
  tv_ratio              total variation of the profile along s against the exact one (> 1: oscillations)
  top-hat   width       10-90% rise of the two edges [cells], mean of the two
  Gaussian  peak        max(q)/max(q_e);  D_num = d(sigma^2)/(2t) / (v h), the numerical diffusivity in
                        units of v*dx (first-order upwind: (1 - nu)/2; first-order HLL/LLF, subsonic: (c_s/v - nu)/2
                        and (1 + c_s/v - nu)/2), meaningful while sigma (final width) < ~0.12;  shift =
                        centroid error [cells]
  sine      amp, shift  amplitude ratio and displacement [cells] of the k=1 Fourier mode;  amp_theory for
                        slope_type 0 in 1-D hydro: |G|^n of the first-order scheme with the Riemann solver's
                        own scalar flux (upwind for hllc/exact/acoustic; the acoustic speeds for llf/hll)
The profile along s = k.x is the cell values binned by their s (one bin per cell in 1-D and grid-aligned
2-D, the diagonals of the grid otherwise).  nu = v_x dt/dx (mean dt = t/nstep).
"""
import contextlib
import glob
import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import acc_cases as A  # noqa: E402


def _visu():
    vis = os.environ.get("ACC_VISU") or os.path.join(HERE, "..", "..", "visu")
    if vis not in sys.path:
        sys.path.insert(0, vis)
    import visu_ramses as V
    return V


def load(rundir, n):
    V = _visu()
    cwd = os.getcwd()
    os.chdir(rundir)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            d = V.load_snapshot(n)
    finally:
        os.chdir(cwd)
    return d


def profile_bins(s, h, q):
    """Bin values by s (multiples of h/2), return sorted s and the bin means."""
    key = np.round(s / (0.5 * h)).astype(np.int64)
    u, inv = np.unique(key, return_inverse=True)
    sums = np.bincount(inv, weights=q)
    cnt = np.bincount(inv)
    return u * 0.5 * h, sums / cnt


def edge_width(s, q):
    """Mean 10-90% width [in s] of the rising and falling edges of a periodic top-hat profile (nan if the
    profile no longer reaches 0.9)."""
    n = len(q)
    ds = 1.0 / n if n > 1 else 1.0
    widths = []
    for sign in (+1, -1):
        y = sign * q
        lo, hi = (0.1, 0.9) if sign > 0 else (-0.9, -0.1)
        i = int(np.argmax(np.roll(y, -1) - y))           # steepest rise between i and i+1
        y = np.roll(y, n // 2 - i)                        # edge in the middle, away from the wrap
        i = n // 2
        a = i
        while a > 0 and y[a] > lo:
            a -= 1
        b = i + 1
        while b < n - 1 and y[b] < hi:
            b += 1
        if y[a] > lo or y[b] < hi:
            widths.append(np.nan)
            continue
        xa = a + (lo - y[a]) / (y[a + 1] - y[a])
        xb = b - 1 + (hi - y[b - 1]) / (y[b] - y[b - 1])
        widths.append((xb - xa) * ds)
    return float(np.mean(widths))


def theory_upwind_amp(c, nu_steps, nstep, cs):
    """|G|^n of first-order (slope 0) 1-D advection with the scalar flux of riemann; None if not modelled."""
    if not (c["solver"] == "hydro" and c["ndim"] == 1 and c["slope"] == 0):
        return None
    v = c["mach"]
    r = c["riemann"]
    if r in ("hllc", "exact", "acoustic"):
        al, be = v, 0.0
    elif r == "llf":
        S = abs(v) + cs
        al, be = 0.5 * (v + S), 0.5 * (v - S)
    elif r == "hll":
        SL, SR = min(v - cs, 0.0), max(v + cs, 0.0)
        al, be = SR * (v - SL) / (SR - SL), SL * (SR - v) / (SR - SL)
    else:
        return None
    lam = nu_steps / v                       # dt/h
    th = 2 * math.pi / c["N"]
    G = 1 - lam * (al * (1 - np.exp(-1j * th)) + be * (np.exp(1j * th) - 1))
    return float(abs(G) ** nstep), float(np.angle(G) * nstep)


def check(rundir):
    c = json.load(open(os.path.join(rundir, "case.json")))
    c["k"] = tuple(c["k"])
    outs = sorted(int(p.rsplit("_", 1)[1]) for p in glob.glob(os.path.join(rundir, "output_[0-9]*")))
    res = dict(case=A.case_name(c), **{k: c[k] for k in ("solver", "ndim", "slope", "riemann", "mach", "N",
                                                         "courant", "periods")}, k=list(c["k"]))
    if len(outs) < 2:
        res["status"] = "FAIL"
        res["error"] = "no final snapshot"
        log = os.path.join(rundir, "run.log")
        if os.path.exists(log):
            tail = [l.strip() for l in open(log, errors="replace") if l.strip()][-3:]
            res["error"] += ": " + " | ".join(tail)[:300]
        json.dump(res, open(os.path.join(rundir, "acc.json"), "w"), indent=1)
        return res
    d0, d1 = load(rundir, outs[0]), load(rundir, outs[-1])
    D0, D1 = d0["data"], d1["data"]
    t = float(D1["time"])
    nstep = int(d1["info"]["nstep_coarse"])
    ndim, N = c["ndim"], c["N"]
    h = 1.0 / N
    vel = c["vel"]
    kx, ky = (c["k"][0], c["k"][1]) if ndim > 1 else (c["k"][0], 0)
    kn = math.hypot(kx, ky)
    shift = (kx * vel[0] + ky * vel[1]) * t
    xc = np.column_stack([D1["x"]] + ([D1["y"]] if ndim > 1 else []))
    order = np.lexsort(xc.T[::-1])
    xc = xc[order]
    s = (kx * xc[:, 0] + (ky * xc[:, 1] if ndim > 1 else 0.0)) % 1.0
    exact = A.exact_on_cells(xc, h, c["k"], ndim, shift)
    rho1 = D1["density"][order]
    rho0 = D0["density"]
    nu = vel[0] * (t / nstep) / h
    cs = math.sqrt(A.GAMMA * float(np.mean(D1["pressure"])) / float(np.mean(rho1)))
    res.update(t=t, nstep=nstep, nu=nu, rho_dev=float(np.abs(rho1 - 1).max()))
    prof = {}
    for ip, p in enumerate(A.PROFILES):
        name = f"scalar_{ip:02d}"
        q = D1[name][order]
        qe = exact[p]
        m0 = float(np.sum(rho0 * D0[name]))
        m1 = float(np.sum(rho1 * q))
        r = dict(L1=float(np.mean(np.abs(q - qe))), L2=float(np.sqrt(np.mean((q - qe) ** 2))),
                 Linf=float(np.abs(q - qe).max()), mass=abs(m1 - m0) / m0,
                 over=float(max(0.0, q.max() - qe.max())), under=float(max(0.0, qe.min() - q.min())),
                 qmin=float(q.min()), qmax=float(q.max()))
        sb, qb = profile_bins(s, h, q)   # uniform spacing along s (h, or h/2 on the diagonal)
        _, qeb = profile_bins(s, h, qe)
        tv = lambda y: float(np.abs(np.diff(np.concatenate([y, y[:1]]))).sum())   # noqa: E731
        r["tv_ratio"] = tv(qb) / tv(qeb)
        h_phys = h                                       # cell size; s-distances / |k| are physical
        if p == "tophat":
            r["width"] = edge_width(sb, qb) / kn / h_phys
            r["width_exact"] = edge_width(sb, qeb) / kn / h_phys
        if p == "gauss":
            cen = (A.G_C + shift) % 1.0
            ds = (sb - cen + 0.5) % 1.0 - 0.5
            def mom(y):
                w = np.clip(y, 0, None)
                m = float(np.sum(w * ds) / np.sum(w))
                return m, float(np.sum(w * (ds - m) ** 2) / np.sum(w))
            m_r, v_r = mom(qb)
            m_e, v_e = mom(qeb)
            r["peak"] = float(qb.max() / qeb.max())
            r["D_num"] = (v_r - v_e) / kn ** 2 / (2 * t) / (c["mach"] * h_phys)
            r["sigma"] = math.sqrt(v_r) / kn          # > ~0.12: the profile wraps the box, D_num saturates
            r["shift"] = (m_r - m_e) / kn / h_phys
        if p == "sine":
            c1 = np.sum(qb * np.exp(-2j * math.pi * sb))
            c1e = np.sum(qeb * np.exp(-2j * math.pi * sb))
            r["amp"] = float(abs(c1) / abs(c1e))
            dphi = float(np.angle(c1 / c1e))
            r["shift"] = -dphi / (2 * math.pi) / kn / h_phys
            th = theory_upwind_amp(c, nu, nstep, cs)
            if th is not None:
                # the theory's phase is relative to the start; the exact solution's phase is -2 pi shift
                r["amp_theory"] = th[0]
        res[p] = r
        prof[f"q_{p}"], prof[f"qe_{p}"] = qb, qeb
    np.savez_compressed(os.path.join(rundir, "profile.npz"), s=sb, **prof)
    res["status"] = "OK"
    json.dump(res, open(os.path.join(rundir, "acc.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    for d in sys.argv[1:]:
        r = check(d)
        if r["status"] != "OK":
            print(f"{r['case']}: {r['status']} {r.get('error', '')}")
            continue
        print(f"{r['case']}: nu={r['nu']:.3f} nstep={r['nstep']}  " + "  ".join(
            f"{p}: L1={r[p]['L1']:.3e} mass={r[p]['mass']:.1e}" for p in A.PROFILES)
              + f"  D_num={r['gauss']['D_num']:.4f} amp={r['sine']['amp']:.4f}"
              + (f" (theory {r['sine']['amp_theory']:.4f})" if "amp_theory" in r["sine"] else ""))
