"""Pure-Python 1D replica of the RAMSES MUSCL-Hancock hydro step (hydro/umuscl.f90 +
hydro/godunov_utils.f90), used to study how passive scalars (elements, ions, CO, dust,
PAH) are advected and how their sums stay consistent with rho.

What is replicated (1D, NENER=0):
  ctoprim   q = (rho, u, P, X_k = U_k/rho)            rho = max(U_1, smallr)
  uslope    slope_type 1 (minmod), 2 (MC), 7 (van Leer), 8 (generalised, slope_theta)
  trace1d   MUSCL-Hancock predictor; passive: qp = a - dax/2 (1+u dt/dx), qm = a + dax/2 (1-u dt/dx)
            rho clipped back to the cell value if the face value < smallr
  riemann   'hllc' (as riemann_hllc: passive flux = F_mass * X_upwind by sign(ustar))
            'llf'  (as riemann_llf:  passive flux = LLF of rho X)
  update    U^{n+1} = U^n - dt/dx (F_{i+1/2} - F_{i-1/2})

Passive-scalar layout (mirrors RAMSES-RTZ/CALIMA): "top" species (elements, CO, dust bins,
PAH bins) whose densities should sum to rho, and "children" (ion stages, H2) whose
densities should sum to their parent element.  Storage convention for the children:
  conv='A'  child slot = y * rho      (y = ion fraction of the element; RAMSES-RTZ now)
  conv='B'  child slot = y * rho_e    (ion mass density; origin/auto_debug)

Remedies (option `scheme`):
  'legacy'   as RAMSES: every scalar limited independently, no constraint
  'cma'      Consistent Multi-fluid Advection (Plewa & Mueller 1999) at the face states:
             top fractions renormalised to sum to 1, then children renormalised to sum to
             their parent's face fraction (B) or to 1 (A).  Nested.
  'cma_y'    (B only) reconstruct the children as y = U_child/U_parent instead of U/rho,
             normalise y to 1 at the face, face child fraction = X_parent,face * y_face
  'renorm'   legacy + post-step renormalisation of the cell values (what RTZ's chemistry
             and auto_debug's interpol cap effectively do)
Optional TVA-like drift of the dust/PAH bins relative to the mixture (`drift`), with or
without the opposite (counter) flux on the gas species (`counterflux`).
"""
import numpy as np

SMALLR = 1e-10
SMALLC = 1e-10


import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from ps_cases import Layout  # noqa: E402  (same layout and constraint tree as the RAMSES builds)


def natural_to_cons(lay, rho, fnat, conv):
    """Conservative passive densities from natural fractions. fnat: (n, ncell)."""
    U = np.zeros_like(fnat)
    U[lay.top] = rho * fnat[lay.top]
    for p, ks in lay.children.items():
        base = rho if conv == "A" else rho * fnat[p]
        U[ks] = base * fnat[ks]
    return U


def cons_to_natural(lay, rho, U, conv):
    f = np.zeros_like(U)
    f[lay.top] = U[lay.top] / rho
    for p, ks in lay.children.items():
        base = rho if conv == "A" else np.maximum(U[p], 1e-300)
        f[ks] = U[ks] / base
    return f


# ----------------------------------------------------------------------------------------------
# hydro kernel
# ----------------------------------------------------------------------------------------------
def slopes(q, slope_type, theta=1.5):
    """q: (nv, ncell incl. ghosts). Limited slopes in the interior, zero at the outermost cells."""
    dq = np.zeros_like(q)
    dl = q[:, 1:-1] - q[:, :-2]
    dr = q[:, 2:] - q[:, 1:-1]
    if slope_type in (1, 2):
        f = slope_type
        dlft, drgt = f * dl, f * dr
        dcen = 0.5 * (dlft + drgt) / f
        dsgn = np.sign(dcen); dsgn[dsgn == 0] = 1.0
        dlim = np.minimum(np.abs(dlft), np.abs(drgt))
        dlim = np.where(dlft * drgt <= 0, 0.0, dlim)
        dq[:, 1:-1] = dsgn * np.minimum(dlim, np.abs(dcen))
    elif slope_type == 7:
        with np.errstate(invalid="ignore", divide="ignore"):
            v = 2 * dl * dr / (dl + dr)
        dq[:, 1:-1] = np.where(dl * dr <= 0, 0.0, v)
    elif slope_type == 8:
        dcen = 0.5 * (dl + dr)
        dsgn = np.sign(dcen); dsgn[dsgn == 0] = 1.0
        dlim = np.minimum(theta * np.abs(dl), theta * np.abs(dr))
        dlim = np.where(dl * dr <= 0, 0.0, dlim)
        dq[:, 1:-1] = dsgn * np.minimum(dlim, np.abs(dcen))
    elif slope_type == 0:
        pass
    else:
        raise ValueError(slope_type)
    return dq


def hllc_flux(rl, ul, pl, rr, ur, pr, gamma):
    smallp = SMALLC ** 2 / gamma
    entho = 1.0 / (gamma - 1.0)
    rl = np.maximum(rl, SMALLR); pl = np.maximum(pl, rl * smallp)
    rr = np.maximum(rr, SMALLR); pr = np.maximum(pr, rr * smallp)
    el = pl * entho; er = pr * entho
    etotl = el + 0.5 * rl * ul * ul; etotr = er + 0.5 * rr * ur * ur
    cl = np.sqrt(np.maximum(gamma * pl / rl, SMALLC ** 2))
    cr = np.sqrt(np.maximum(gamma * pr / rr, SMALLC ** 2))
    SL = np.minimum(ul, ur) - np.maximum(cl, cr)
    SR = np.maximum(ul, ur) + np.maximum(cl, cr)
    rcl = rl * (ul - SL); rcr = rr * (SR - ur)
    ustar = (rcr * ur + rcl * ul + (pl - pr)) / (rcr + rcl)
    pstar = (rcr * pl + rcl * pr + rcl * rcr * (ul - ur)) / (rcr + rcl)
    rstarl = rl * (SL - ul) / (SL - ustar); rstarr = rr * (SR - ur) / (SR - ustar)
    etsl = ((SL - ul) * etotl - pl * ul + pstar * ustar) / (SL - ustar)
    etsr = ((SR - ur) * etotr - pr * ur + pstar * ustar) / (SR - ustar)
    c1 = SL > 0; c2 = (~c1) & (ustar > 0); c3 = (~c1) & (~c2) & (SR > 0)
    c4 = ~(c1 | c2 | c3)
    ro = np.select([c1, c2, c3, c4], [rl, rstarl, rstarr, rr])
    uo = np.select([c1, c2, c3, c4], [ul, ustar, ustar, ur])
    po = np.select([c1, c2, c3, c4], [pl, pstar, pstar, pr])
    eo = np.select([c1, c2, c3, c4], [etotl, etsl, etsr, etotr])
    fm = ro * uo
    return fm, ro * uo * uo + po, (eo + po) * uo, ustar


def llf_flux(rl, ul, pl, rr, ur, pr, gamma):
    rl = np.maximum(rl, SMALLR); rr = np.maximum(rr, SMALLR)
    cl = np.sqrt(gamma * pl / rl); cr = np.sqrt(gamma * pr / rr)
    cmax = np.maximum(np.abs(ul) + cl, np.abs(ur) + cr)
    El = pl / (gamma - 1) + 0.5 * rl * ul ** 2; Er = pr / (gamma - 1) + 0.5 * rr * ur ** 2
    f1 = 0.5 * (rl * ul + rr * ur - cmax * (rr - rl))
    f2 = 0.5 * (rl * ul ** 2 + pl + rr * ur ** 2 + pr - cmax * (rr * ur - rl * ul))
    f3 = 0.5 * (ul * (El + pl) + ur * (Er + pr) - cmax * (Er - El))
    return f1, f2, f3, cmax


class Sim:
    def __init__(self, lay, x, rho, u, p, fnat, conv="A", scheme="legacy", slope_type=1,
                 riemann="hllc", gamma=1.4, courant=0.8, bc="periodic", theta=1.5,
                 drift=None, counterflux=False, tva_face_clip=True):
        self.lay, self.conv, self.scheme = lay, conv, scheme
        self.slope_type, self.riemann, self.gamma = slope_type, riemann, gamma
        self.courant, self.bc, self.theta = courant, bc, theta
        self.x = x; self.dx = x[1] - x[0]; self.N = len(x)
        self.U = np.zeros((3 + lay.n, self.N))
        self.U[0] = rho; self.U[1] = rho * u
        self.U[2] = p / (gamma - 1) + 0.5 * rho * u * u
        self.U[3:] = natural_to_cons(lay, rho, fnat, conv)
        self.drift = drift              # None or callable(x) -> drift velocity of the TVA bins
        self.counterflux = counterflux
        self.t = 0.0
        self.counters = {"renorm_cells": 0, "neg_face": 0, "neg_cell": 0}

    # ghost cells (2 each side)
    def ghost(self, U):
        if self.bc == "periodic":
            return np.concatenate([U[:, -2:], U, U[:, :2]], axis=1)
        return np.concatenate([U[:, :1], U[:, :1], U, U[:, -1:], U[:, -1:]], axis=1)

    def dt_cfl(self):
        rho = np.maximum(self.U[0], SMALLR); u = self.U[1] / rho
        p = (self.gamma - 1) * (self.U[2] - 0.5 * rho * u * u)
        c = np.sqrt(self.gamma * np.maximum(p, 1e-30) / rho)
        vmax = np.max(np.abs(u) + c)
        if self.drift is not None:
            vmax = max(vmax, np.max(np.abs(u) + np.abs(self.drift(self.x))))
        return self.courant * self.dx / vmax

    def normalise_faces(self, X):
        """Nested CMA on face fractions X (nspec, nface)."""
        lay = self.lay
        X = np.maximum(X, 0.0)
        s = X[lay.top].sum(axis=0)
        X[lay.top] /= np.where(s > 0, s, 1.0)
        for p, ks in lay.children.items():
            s = X[ks].sum(axis=0)
            target = 1.0 if self.conv == "A" else X[p]
            X[ks] *= np.where(s > 0, target / np.where(s > 0, s, 1.0), 0.0)
        return X

    def step(self, dt):
        lay, g = self.lay, self.gamma
        Ug = self.ghost(self.U)
        rho = np.maximum(Ug[0], SMALLR)
        u = Ug[1] / rho
        p = np.maximum((g - 1) * (Ug[2] - 0.5 * rho * u * u), rho * SMALLC ** 2 / g / (g - 1) * (g - 1))
        X = Ug[3:] / rho                                    # ctoprim: X = U/rho
        if self.scheme == "cma_y":
            # children as fractions of the parent element
            for pp, ks in lay.children.items():
                X[ks] = Ug[3 + ks] / np.maximum(Ug[3 + pp], 1e-300)
        q = np.vstack([rho, u, p, X])
        dq = slopes(q, self.slope_type, self.theta)
        dtdx = dt / self.dx
        r, dr, du, dp_ = q[0], dq[0], dq[1], dq[2]
        sr0 = -u * dr - du * r
        sp0 = -u * dp_ - du * g * p
        su0 = -u * du - dp_ / r
        qp = np.empty_like(q); qm = np.empty_like(q)
        qp[0] = r - 0.5 * dr + sr0 * dtdx * 0.5; qm[0] = r + 0.5 * dr + sr0 * dtdx * 0.5
        qp[0] = np.where(qp[0] < SMALLR, r, qp[0]); qm[0] = np.where(qm[0] < SMALLR, r, qm[0])
        qp[1] = u - 0.5 * du + su0 * dtdx * 0.5; qm[1] = u + 0.5 * du + su0 * dtdx * 0.5
        qp[2] = p - 0.5 * dp_ + sp0 * dtdx * 0.5; qm[2] = p + 0.5 * dp_ + sp0 * dtdx * 0.5
        a, da = q[3:], dq[3:]
        sa0 = -u * da
        qp[3:] = a - 0.5 * da + sa0 * dtdx * 0.5
        qm[3:] = a + 0.5 * da + sa0 * dtdx * 0.5
        # faces i-1/2 for interior cells: left state = qm[cell i-1], right = qp[cell i]
        # ghost array index j = cell + 2; faces between j=1|2 ... N+1|N+2
        L = qm[:, 1:self.N + 2]; R = qp[:, 2:self.N + 3]
        XL, XR = L[3:].copy(), R[3:].copy()
        self.counters["neg_face"] += int(np.sum(XL < 0) + np.sum(XR < 0))
        if self.scheme == "cma_y":
            for pp, ks in lay.children.items():
                for XX in (XL, XR):
                    yy = np.maximum(XX[ks], 0.0); s = yy.sum(axis=0)
                    yy /= np.where(s > 0, s, 1.0)
                    XX[ks] = yy
            # top-level normalisation, then child fraction = parent * y
            for XX in (XL, XR):
                XX[lay.top] = np.maximum(XX[lay.top], 0)
                s = XX[lay.top].sum(axis=0); XX[lay.top] /= np.where(s > 0, s, 1.0)
                for pp, ks in lay.children.items():
                    XX[ks] = XX[ks] * XX[pp]
        elif self.scheme == "cma":
            XL = self.normalise_faces(XL); XR = self.normalise_faces(XR)
        if self.riemann == "hllc":
            fm, fmom, fen, ustar = hllc_flux(L[0], L[1], L[2], R[0], R[1], R[2], g)
            Xup = np.where(ustar > 0, XL, XR)
            if self.scheme == "cma_y" and self.conv == "A":
                raise ValueError("cma_y needs conv B")
            fX = fm * Xup
        else:
            fm, fmom, fen, cmax = llf_flux(L[0], L[1], L[2], R[0], R[1], R[2], g)
            rl = np.maximum(L[0], SMALLR); rr = np.maximum(R[0], SMALLR)
            fX = 0.5 * (rl * L[1] * XL + rr * R[1] * XR - cmax * (rr * XR - rl * XL))
        F = np.vstack([fm, fmom, fen, fX])                   # (nv, N+1)
        if self.drift is not None:
            F = F + self.drift_flux(dt)
        self.U = self.U - dtdx * (F[:, 1:] - F[:, :-1])
        self.t += dt
        if self.scheme == "renorm":
            self.renormalise()

    def drift_flux(self, dt):
        """First-order upwind drift flux of the TVA bins at the faces, plus the opposite
        (counter) flux on the gas species if requested.  Returns (nv, N+1)."""
        lay = self.lay
        Ug = self.ghost(self.U)
        xf = np.concatenate([[self.x[0] - 0.5 * self.dx], self.x + 0.5 * self.dx])
        w = self.drift(xf)
        F = np.zeros((3 + lay.n, self.N + 1))
        Ul, Ur = Ug[:, 1:self.N + 2], Ug[:, 2:self.N + 3]
        up = np.where(w > 0, Ul[3 + lay.tva], Ur[3 + lay.tva])
        Fd = w * up
        F[3 + lay.tva] = Fd
        if self.counterflux:
            Ftot = Fd.sum(axis=0)                         # mass flux the gas must carry back
            # gas moves with -Ftot: upwind cell with respect to the gas drift
            gl = Ul[3 + lay.gas_top]; gr = Ur[3 + lay.gas_top]
            gup = np.where(-Ftot > 0, gl, gr)
            frac = gup / np.maximum(gup.sum(axis=0), 1e-300)
            F[3 + lay.gas_top] = -Ftot * frac
            if self.conv == "B":   # ions follow their element
                for pp, ks in lay.children.items():
                    if pp in lay.gas_top:
                        cu = np.where(-Ftot > 0, Ul[3 + ks], Ur[3 + ks])
                        y = cu / np.maximum(cu.sum(axis=0), 1e-300)
                        F[3 + ks] = F[3 + pp] * y
        return F

    def renormalise(self):
        lay = self.lay
        P = self.U[3:]
        P[:] = np.maximum(P, 0.0)
        s = P[lay.top].sum(axis=0)
        bad = np.abs(s / self.U[0] - 1) > 1e-12
        self.counters["renorm_cells"] += int(bad.sum())
        P[lay.top] *= self.U[0] / s
        for p, ks in lay.children.items():
            s = P[ks].sum(axis=0)
            target = self.U[0] if self.conv == "A" else P[p]
            P[ks] *= np.where(s > 0, target / np.where(s > 0, s, 1), 0)

    # ------------------------------------------------------------------------------------------
    def diagnostics(self):
        lay = self.lay
        rho = self.U[0]; P = self.U[3:]
        d = {}
        d["mass"] = rho.sum() * self.dx
        d["spec_mass"] = P.sum(axis=1) * self.dx
        d["top_err"] = np.max(np.abs(P[lay.top].sum(axis=0) / rho - 1))
        ce = 0.0
        for p, ks in lay.children.items():
            target = rho if self.conv == "A" else P[p]
            m = target > 1e-300
            if m.any():
                ce = max(ce, np.max(np.abs(P[ks][:, m].sum(axis=0) / target[m] - 1)))
        d["child_err"] = ce
        ca = 0.0                          # same, in units of the parent's peak density
        for p, ks in lay.children.items():
            target = rho if self.conv == "A" else P[p]
            if target.max() > 0:
                ca = max(ca, np.max(np.abs(P[ks].sum(axis=0) - target)) / target.max())
        d["child_abs"] = ca
        d["min_frac"] = np.min(P / rho)
        f = cons_to_natural(lay, rho, P, self.conv)
        d["fnat"] = f
        # physical ion masses as the chemistry would use them: y * rho_e
        im = np.zeros(lay.n)
        for p, ks in lay.children.items():
            im[ks] = (f[ks] * P[p]).sum(axis=1) * self.dx
        d["ion_mass"] = im
        return d

    def run(self, tend, nout=20, maxstep=200000):
        hist = [(self.t, self.diagnostics())]
        tout = np.linspace(0, tend, nout + 1)[1:]
        io = 0; n = 0
        while self.t < tend * (1 - 1e-12) and n < maxstep:
            dt = min(self.dt_cfl(), tend - self.t)
            self.step(dt); n += 1
            if self.t >= tout[io] * (1 - 1e-12):
                hist.append((self.t, self.diagnostics())); io += 1
        self.nstep = n
        return hist
