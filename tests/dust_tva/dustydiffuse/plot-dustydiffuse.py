#!/usr/bin/env python3
"""DustyDiffuse (Lebreuilly et al. 2019, sect. 4.2): pressure-driven dust diffusion, checked against the
converged solution of the equation the code solves.

The mixture is static and isothermal (c_s = 1) with tau_s = 0.1, so the drift w = -(1-eps) tau_s c_s^2 deps/dx
turns the dust update into a nonlinear diffusion,
    deps/dt = tau_s c_s^2 d/dx [ eps (1-eps) deps/dx ],
from a parabola of peak 0.1 and half-width 0.2. The reference is that equation solved by a method of lines on
4000 cells. The Barenblatt-Pattle solution (its eps << 1 limit) is drawn for comparison only: the factor 1-eps
puts it 1-1.5% off the solution.

Namelists:
  dustydiffuse.nml      tva_scheme='calima': L19's AMR set-up (levels 4-10), MC, the scheme's own time step
  dustydiffuse_l19.nml  tva_scheme='l19' (with its drift pre-pass): 2^8 cells, MinMod, a fixed dt = 1.5e-3, inside
                        the scheme's stability limit (at the eq. 26 time step alone it is unstable from 2^8 cells)

Checks: the relative L2 error against the reference at every output; the dust mass to round-off; eps >= 0.

Usage: ./plot-dustydiffuse.py [namelist]      (runs ./ramses_dust_test1d, or $RAMSES_EXE)
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
from scipy.integrate import solve_ivp
from scipy.sparse import diags

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

NML = sys.argv[1] if len(sys.argv) > 1 else "dustydiffuse.nml"
NAME = os.path.splitext(NML)[0]
EXE = os.environ.get("RAMSES_EXE", "./ramses_dust_test1d")
TS, EPS0, XC = 0.1, 0.1, 0.2            # tau_s c_s^2, peak and half-width of the initial parabola
TOL_L2 = {"dustydiffuse": 1.5e-2, "dustydiffuse_l19": 3.0e-3}.get(NAME, 1.5e-2)
TOL_MASS = 1.0e-12


def reference(times, n=4000):
    """deps/dt = d/dx[TS eps(1-eps) deps/dx], periodic, by a method of lines (face diffusivity from the face mean)."""
    dx = 1.0 / n
    x = (np.arange(n) + 0.5) * dx
    e0 = np.clip(EPS0 * (1.0 - ((x - 0.5) / XC) ** 2), 0.0, None)

    def rhs(t, e):
        ef = 0.5 * (e + np.roll(e, -1))
        flux = -TS * ef * (1.0 - ef) * (np.roll(e, -1) - e) / dx
        return -(flux - np.roll(flux, 1)) / dx

    jac = diags([1.0] * 5, [-(n - 1), -1, 0, 1, n - 1], shape=(n, n))
    sol = solve_ivp(rhs, (0.0, max(times)), e0, t_eval=sorted(times), method="BDF", rtol=1e-8, atol=1e-12,
                    jac_sparsity=jac)
    return x, {t: sol.y[:, k] for k, t in enumerate(sol.t)}


def barenblatt(x, t):
    c = (EPS0 * XC / np.sqrt(6.0)) ** (2.0 / 3.0)
    tt = TS * (t + c ** 3 / (TS * EPS0 ** 3))
    return np.clip(tt ** (-1.0 / 3.0) * (c - (x - 0.5) ** 2 / (6.0 * tt ** (2.0 / 3.0))), 0.0, None)


print(f"--- Running DustyDiffuse ({NML}) ---")
for old in glob.glob("output_*"):
    shutil.rmtree(old)
subprocess.run([EXE, NML], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

snaps = []
n = 1
while os.path.isdir(f"output_{n:05d}"):
    d = visu_ramses.load_snapshot(n)["data"]
    o = np.argsort(d["x"])
    snaps.append({k: (d[k][o] if np.ndim(d[k]) else d[k]) for k in ("x", "dx", "density", "DustBin_01", "time")})
    n += 1
times = [s["time"] for s in snaps[1:]]
xr, ref = reference(times)

mass0 = np.sum(snaps[0]["density"] * snaps[0]["DustBin_01"] * snaps[0]["dx"])
failures = []
fig, ax = plt.subplots(figsize=(5.5, 4))
cols = plt.get_cmap("viridis")(np.linspace(0.0, 0.8, len(times)))
print(f"{'t':>8} {'L2 vs reference':>16} {'L2 vs B-P':>10} {'dM/M':>10} {'min eps':>10}")
for s, c in zip(snaps[1:], cols):
    x, e, dx, t = s["x"], s["DustBin_01"], s["dx"], s["time"]
    r = np.interp(x, xr, ref[t])
    l2 = np.sqrt(np.sum(dx * (e - r) ** 2) / np.sum(dx * r ** 2))
    l2bp = np.sqrt(np.sum(dx * (e - barenblatt(x, t)) ** 2) / np.sum(dx * barenblatt(x, t) ** 2))
    dm = (np.sum(s["density"] * e * dx) - mass0) / mass0
    print(f"{t:8.4f} {l2:16.3e} {l2bp:10.3e} {dm:10.1e} {e.min():10.2e}")
    if l2 > TOL_L2:
        failures.append(f"t = {t:.3f}: L2 against the reference {l2:.3e} > {TOL_L2:.1e}")
    if abs(dm) > TOL_MASS:
        failures.append(f"t = {t:.3f}: dust mass error {dm:.2e} > {TOL_MASS:.0e}")
    if e.min() < 0.0:
        failures.append(f"t = {t:.3f}: negative dust fraction {e.min():.2e}")
    ax.plot(x, e, "o", color=c, ms=2.5, mfc="none", label=f"$t = {t:.1f}$")
    ax.plot(xr, ref[t], "-", color=c, lw=1.2)
    ax.plot(xr, barenblatt(xr, t), "--", color="0.4", lw=0.7)
ax.plot([], [], "k-", lw=1.2, label="reference (converged PDE)")
ax.plot([], [], "--", color="0.4", lw=0.7, label="Barenblatt-Pattle")
ax.set_xlim(0.2, 0.8); ax.set_xlabel("x"); ax.set_ylabel("dust fraction")
ax.legend(fontsize=8, loc="upper right")
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
