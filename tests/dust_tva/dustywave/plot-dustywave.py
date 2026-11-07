#!/usr/bin/env python3
"""DustyWave (Lebreuilly et al. 2019, sect. 4.4, their Figs. 5-6): a sound wave in a dusty mixture, checked
against the linearised terminal-velocity equations, which the code solves.

A uniform mixture (rho0 = 2, eps0 = 0.5, c_s = 1, isothermal up to gamma - 1 = 1.2e-6) carries one mode,
delta rho = rho0 A sin(kx), delta v = A c_s sin(kx), delta P = (1-eps0) c_s^2 delta rho, eps unperturbed, with
k = 2 pi and A = 1e-4. The drag coefficient is K = 50 (t_s = eps rho / K), on 64 cells at dt = 1e-4 to t = 4.5.
The reference is the exact solution of the linearised TVA equations. The linearised two-fluid solution (Laibe &
Price 2011) differs from it by 0.3% in the dust amplitude at this K, the error of the approximation, so it is not
used for the checks.

Namelists:
  dustywave.nml      tva_scheme='calima'
  dustywave_l19.nml  tva_scheme='l19' (with its drift pre-pass)
Both use MC slopes: with MinMod the gas solver alone takes 3-5% off the amplitude at 64 cells.

Checks: the L2 error of delta rho_d and of v against the reference, relative to the amplitude; the amplitude of
delta rho_d; the dust mass to round-off.

Usage: ./plot-dustywave.py [namelist]      (runs ./ramses_dust_test1d, or $RAMSES_EXE)
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
from scipy.linalg import expm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

NML = sys.argv[1] if len(sys.argv) > 1 else "dustywave.nml"
NAME = os.path.splitext(NML)[0]
EXE = os.environ.get("RAMSES_EXE", "./ramses_dust_test1d")
RHO0, EPS0, CS, KDRAG, AMP, KW = 2.0, 0.5, 1.0, 50.0, 1e-4, 2.0 * np.pi
TOL = {"dustywave": dict(l2_rhod=2e-2, l2_v=3e-2, amp=1e-2),
       "dustywave_l19": dict(l2_rhod=2e-2, l2_v=3e-2, amp=1e-2)}.get(NAME, dict(l2_rhod=2e-2, l2_v=3e-2, amp=1e-2))
TOL_MASS = 1e-12


def tva(x, t):
    """The linearised one-species TVA equations: drho, v and deps, exact in time (matrix exponential)."""
    tau = EPS0 * RHO0 / KDRAG
    ik, k2 = 1j * KW, KW ** 2
    dpdr, dpde = CS ** 2 * (1.0 - EPS0), -CS ** 2 * RHO0          # dP = dpdr drho + dpde deps
    m = np.zeros((3, 3), dtype=complex)
    m[0, 1] = -ik * RHO0
    m[1, 0], m[1, 2] = -ik * dpdr / RHO0, -ik * dpde / RHO0
    c = k2 * EPS0 * (tau - EPS0 * tau) / RHO0                     # d deps/dt = k^2 eps (t_s - eps t_s) dP / rho0
    m[2, 0], m[2, 2] = c * dpdr, c * dpde
    y = expm(m * t) @ np.array([RHO0 * AMP, AMP, 0.0], dtype=complex)
    dr, v, de = (np.imag(yk * np.exp(ik * x)) for yk in y)
    return dict(drho_d=EPS0 * dr + RHO0 * de, v=v)


def amplitude(f, x):
    return np.hypot(2.0 * np.mean(f * np.sin(KW * x)), 2.0 * np.mean(f * np.cos(KW * x)))


print(f"--- Running DustyWave ({NML}) ---")
for old in glob.glob("output_*"):
    shutil.rmtree(old)
subprocess.run([EXE, NML], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def load(n):
    d = visu_ramses.load_snapshot(n)["data"]
    o = np.argsort(d["x"])
    return {k: (d[k][o] if np.ndim(d[k]) else d[k]) for k in ("x", "dx", "density", "velocity_x", "DustBin_01", "time")}

d0, d = load(1), load(2)
x, t = d["x"], d["time"]
num = dict(drho_d=d["density"] * d["DustBin_01"] - RHO0 * EPS0, v=d["velocity_x"])
ref = tva(x, t)
rms = lambda a: np.sqrt(np.mean(a ** 2))
l2_rhod = rms(num["drho_d"] - ref["drho_d"]) / amplitude(ref["drho_d"], x)
l2_v = rms(num["v"] - ref["v"]) / amplitude(ref["v"], x)
damp = amplitude(num["drho_d"], x) / amplitude(ref["drho_d"], x) - 1.0
mass = lambda s: np.sum(s["density"] * s["DustBin_01"] * s["dx"])
dm = (mass(d) - mass(d0)) / mass(d0)
print(f"  t = {t:.3f}: L2(drho_d)/amp {l2_rhod:.3e}  L2(v)/amp {l2_v:.3e}  drho_d amplitude {damp:+.3%}  dM/M {dm:.1e}")

failures = []
if l2_rhod > TOL["l2_rhod"]:
    failures.append(f"L2(drho_d)/amplitude {l2_rhod:.3e} > {TOL['l2_rhod']:.0e}")
if l2_v > TOL["l2_v"]:
    failures.append(f"L2(v)/amplitude {l2_v:.3e} > {TOL['l2_v']:.0e}")
if abs(damp) > TOL["amp"]:
    failures.append(f"drho_d amplitude off by {damp:+.3%} (tolerance {TOL['amp']:.0%})")
if abs(dm) > TOL_MASS:
    failures.append(f"dust mass error {dm:.2e} > {TOL_MASS:.0e}")

xe = np.linspace(0.0, 1.0, 400)
fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
for ax, key, lab in zip(axs, ("v", "drho_d"), ("velocity perturbation", "dust density perturbation")):
    ax.plot(x, num[key] / AMP, "o", ms=3, mfc="none", label="RAMSES")
    ax.plot(xe, tva(xe, t)[key] / AMP, "k-", lw=1.2, label="linear TVA")
    ax.set_xlabel("x"); ax.set_ylabel(f"{lab} / A")
axs[0].legend(fontsize=8)
fig.suptitle(f"{NAME}, K = {KDRAG:g}, t = {t:.2f}")
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
