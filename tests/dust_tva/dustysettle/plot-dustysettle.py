#!/usr/bin/env python3
"""DustySettle (Lebreuilly et al. 2019, sect. 4.5; Price & Laibe 2015): dust settling to the mid-plane of a disk,
checked against the exact solution by characteristics.

The vertical of a disk at r = 50 au around 1 Msun (H/r = 0.05), 1 mm grains of 3 g/cm^3, eps0 = 4.99e-3, 256 cells
over 20 au, 1, 5 and 10 orbits. In 1D, gravity_type=2 with softening r is exactly f = -GM y/(y^2 + r^2)^(3/2). The
gas stays in hydrostatic equilibrium, so the drift w(y) = (1-eps) t_s(y) f(y) does not depend on time, and the dust
equation is solved exactly along dy/dt = w: rho_d(y, t) = rho_d0(y0) w(y0)/w(y), with y0 the foot of the
characteristic. t_s is the code's Epstein time, sqrt(pi gamma/8) s_grain a / (rho_g c_s). The drift is not capped
(tva_wmax_cs = 0), as in L19.

Checks: eps at the mid-plane against the exact value at every output; the dust mass within one scale height
(|y| < H = 2.5 au), the settled dust, against the exact one (a robust measure: by 10 orbits the sharp edge of the
settling dust is inside that range, and a point-wise error there measures the front's sub-cell position); the total
dust mass to round-off.

Usage: ./plot-dustysettle.py [namelist]      (runs ./ramses_dust_test1d, or $RAMSES_EXE)
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

NML = sys.argv[1] if len(sys.argv) > 1 else "dustysettle.nml"
NAME = os.path.splitext(NML)[0]
EXE = os.environ.get("RAMSES_EXE", "./ramses_dust_test1d")
# dustysettle.nml, code units (1 au, 1/Omega at 50 au, 1e-13 g/cm^3)
GM, R, YMID, L = 1.25e5, 50.0, 10.0, 20.0
RHO0, P0, EPS0, GAMMA = 6.0, 37.3129, 4.99e-3, 1.6666667
SA = (3.0 / 1e-13) * (0.1 / 1.495978707e13)           # s_grain * a in code units
TOL_MID, TOL_INNER, TOL_MASS = 2e-2, 2e-2, 1e-12


def gas(y):
    cs2 = P0 / RHO0
    rho = RHO0 * np.exp(GM * (1.0 / np.sqrt(y ** 2 + R ** 2) - 1.0 / R) / cs2)
    return rho, cs2 * rho


def drift(y):
    rho, p = gas(y)
    rg = (1.0 - EPS0) * rho
    cs = np.sqrt(GAMMA * p / rg)
    tau = np.sqrt(np.pi * GAMMA / 8.0) * SA / (rg * cs)
    return (1.0 - EPS0) * tau * (-GM * y / (y ** 2 + R ** 2) ** 1.5)


YG = np.geomspace(1e-10, L / 2, 40001)
_f = YG / np.abs(drift(YG))
G = -np.concatenate([np.cumsum((np.diff(np.log(YG)) * 0.5 * (_f[1:] + _f[:-1]))[::-1])[::-1], [0.0]])   # time to the edge


def exact_rhod(y, t):
    """rho_d(y, t): the foot y0 of each characteristic from G(y0) = G(y) + t; no dust where it lies past the edge."""
    a = np.clip(np.abs(y), YG[0], L / 2)
    g0 = np.interp(np.log(a), np.log(YG), G) + t
    y0 = np.exp(np.interp(np.minimum(g0, 0.0), G, np.log(YG)))
    return np.where(g0 > 0.0, 0.0, EPS0 * gas(y0)[0] * np.abs(drift(y0)) / np.abs(drift(a)))


def exact(y, t):
    rd = exact_rhod(y, t)
    return rd / ((1.0 - EPS0) * gas(np.clip(np.abs(y), YG[0], L / 2))[0] + rd)


print(f"--- Running DustySettle ({NML}) ---")
for old in glob.glob("output_*"):
    shutil.rmtree(old)
subprocess.run([EXE, NML], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def load(n):
    d = visu_ramses.load_snapshot(n)["data"]
    o = np.argsort(d["x"])
    return {k: (d[k][o] if np.ndim(d[k]) else d[k]) for k in ("x", "dx", "density", "DustBin_01", "time")}

snaps, n = [], 1
while os.path.isdir(f"output_{n:05d}"):
    snaps.append(load(n)); n += 1
mass = lambda s: np.sum(s["density"] * s["DustBin_01"] * s["dx"])
failures = []
fig, ax = plt.subplots(figsize=(5.5, 4))
for s, c in zip(snaps[1:], ("C0", "C1", "C2")):
    y, e, t = s["x"] - YMID, s["DustBin_01"], s["time"]
    ex = exact(y, t)
    i0 = np.argmin(np.abs(y))
    dmid = e[i0] / ex[i0] - 1.0
    inner = np.abs(y) < 2.5                                       # |y| < H
    # the exact mass over the same cells' extent, by the trapezoid rule on a fine grid
    yf = np.linspace((y - 0.5 * s["dx"])[inner].min(), (y + 0.5 * s["dx"])[inner].max(), 20001)
    rf = exact_rhod(yf, t)
    m_in, m_in_ex = np.sum((s["density"] * e * s["dx"])[inner]), np.sum(0.5 * (rf[1:] + rf[:-1]) * np.diff(yf))
    din = m_in / m_in_ex - 1.0
    dm = (mass(s) - mass(snaps[0])) / mass(snaps[0])
    print(f"  {t / (2 * np.pi):5.2f} orbits: eps(0) {e[i0]:.5e} exact {ex[i0]:.5e} ({dmid:+.2%})  mass within H {din:+.2%}  dM/M {dm:.1e}")
    if abs(dmid) > TOL_MID:
        failures.append(f"t = {t / (2 * np.pi):.2f} orbits: mid-plane eps off by {dmid:+.2%} (tolerance {TOL_MID:.0%})")
    if abs(din) > TOL_INNER:
        failures.append(f"t = {t / (2 * np.pi):.2f} orbits: dust mass within H off by {din:+.2%} (tolerance {TOL_INNER:.0%})")
    if abs(dm) > TOL_MASS:
        failures.append(f"t = {t / (2 * np.pi):.2f} orbits: dust mass error {dm:.2e} > {TOL_MASS:.0e}")
    ax.semilogy(y, e, "o", color=c, ms=2.5, mfc="none", label=f"{t / (2 * np.pi):.0f} orbits")
    ax.semilogy(y, np.maximum(ex, 1e-12), "-", color=c, lw=1.0)
ax.plot([], [], "k-", lw=1.0, label="exact (characteristics)")
ax.set_xlim(-L / 2, L / 2); ax.set_ylim(1e-6, 2e-2); ax.set_xlabel("y [au]"); ax.set_ylabel("dust fraction")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{NAME}.pdf", bbox_inches="tight")
fig.savefig(f"{NAME}.png", dpi=200, bbox_inches="tight")
print(f"  wrote {NAME}.pdf / .png")

if failures:
    for f in failures:
        print(f"FAILED: {f}")
    print(f"{NAME:30s} failed")
    sys.exit(1)
print(f"{NAME:30s} passed")
