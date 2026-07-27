#!/usr/bin/env python3
"""DustyAdvect3D: 3D constant-drift advection of a dust cube, checked against the
analytic solution.

Unlike the six 1D TVA tests, this exercises all three flux directions at once, with a
deliberately anisotropic drift w = (1.0, 0.5, 0.25) so a transverse index mix-up cannot
cancel out -- that is the failure mode which let the c_s_face bug in
get_dust_courant_dt survive (it used the x-neighbour on the y and z faces).

The gas is static and all dust chemistry is off, so the dust bins are advected and
nothing else. The analytic solution is therefore the initial cube translated by w*t,
periodic in the box.

Checks (physical tolerances, not a bit-for-bit reference comparison -- the latter, as
used by dustyshell_mulbin at 1e-12, breaks on any compiler or decomposition change):
  1. total dust mass conserved to round-off
  2. dust density non-negative
  3. per-direction centroid displacement equals w*t

Produces dustyadvect3d.pdf: profiles along each axis against the analytic top hat.

Usage
-----
  ./plot-dustyadvect3d.py [nout] [--nodust-suffix N]

`nout` defaults to the last output. The dust field is summed over all DustBin_* columns
present, so the script works for any NDUST.
"""

import glob
import os
import re
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

# --- must match the namelist ---------------------------------------------------------
W_DRIFT = np.array([1.0, 0.5, 0.25])
BOXLEN = 1.0
CUBE_CTR = np.array([0.25, 0.25, 0.25])
CUBE_LEN = 0.25                       # full width (length_x); 'square' is a max-norm box
AXES = ("x", "y", "z")


def dust_total(d):
    """Sum every DustBin_* field, so this works for any NDUST."""
    keys = sorted(k for k in d if re.match(r"^DustBin_\d+$", str(k)))
    if not keys:
        raise SystemExit("no DustBin_* fields in the output -- was CALIMA compiled in?")
    return sum(np.asarray(d[k]) for k in keys), keys


def last_output():
    outs = sorted(glob.glob("output_[0-9]" * 1 + "*"))
    outs = [o for o in outs if re.match(r"^output_\d{5}$", o)]
    if not outs:
        raise SystemExit("no output_NNNNN directories found")
    return int(outs[-1].split("_")[1])


nout = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else last_output()
data = visu_ramses.load_snapshot(nout)
d = data["data"]

x, y, z = (np.asarray(d[k]) for k in ("x", "y", "z"))
dx = np.asarray(d["dx"])
rho = np.asarray(d["density"])
eps, binkeys = dust_total(d)
t = float(d["time"])

vol = dx ** 3
rho_dust = rho * eps
mass = float(np.sum(rho_dust * vol))

# Reference the t=0 snapshot rather than trying to infer the initial cube amplitude
# from the evolved field -- upwind advection diffuses the edges, so a percentile of the
# final state underestimates the peak and would fake a conservation error.
ref = visu_ramses.load_snapshot(1)["data"]
eps0, _ = dust_total(ref)
rho0 = np.asarray(ref["density"])
vol0 = np.asarray(ref["dx"]) ** 3
m_expect = float(np.sum(rho0 * eps0 * vol0))
eps_bg = float(np.min(eps0))
eps_hi = float(np.max(eps0))
err_mass = abs(mass - m_expect) / m_expect

eps_min = float(np.min(eps))

# Centroid of the dust in excess of the background: the uniform background carries no
# information and would otherwise dominate the first moment.
excess = np.clip(rho_dust - eps_bg, 0.0, None) * vol
tot = float(np.sum(excess))
com = np.array([float(np.sum(excess * c) / tot) for c in (x, y, z)])
expected = np.mod(CUBE_CTR + W_DRIFT * t, BOXLEN)
delta = np.mod(com - expected + 0.5 * BOXLEN, BOXLEN) - 0.5 * BOXLEN

print("DustyAdvect3D   t = %.6f   bins summed: %s" % (t, ",".join(binkeys)))
print("  dust mass       = %.10e  (analytic %.10e, rel err %.3e)" % (mass, m_expect, err_mass))
print("  min dust        = %.6e" % eps_min)
print("  centroid        = [%.5f %.5f %.5f]" % tuple(com))
print("  analytic        = [%.5f %.5f %.5f]" % tuple(expected))
print("  offset          = [%+.5f %+.5f %+.5f]" % tuple(delta))

# --- figure: profile along each axis vs the analytic top hat --------------------------
fig, axarr = plt.subplots(1, 3, figsize=(13, 4))
cell = BOXLEN / 2 ** int(round(np.log2(BOXLEN / np.min(dx))))
for ia, (ax, coord) in enumerate(zip(axarr, (x, y, z))):
    # collapse onto this axis: mass-weighted mean dust fraction per slab
    nb = int(round(BOXLEN / np.min(dx)))
    edges = np.linspace(0.0, BOXLEN, nb + 1)
    num, _ = np.histogram(coord, bins=edges, weights=rho_dust * vol)
    den, _ = np.histogram(coord, bins=edges, weights=vol)
    prof = num / np.maximum(den, 1e-300)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    ax.step(ctr, prof, where="mid", label="RAMSES", lw=1.6)

    # analytic: top hat of width CUBE_LEN centred on the drifted position, periodic
    c0 = np.mod(CUBE_CTR[ia] + W_DRIFT[ia] * t, BOXLEN)
    off = np.mod(ctr - c0 + 0.5 * BOXLEN, BOXLEN) - 0.5 * BOXLEN
    # the other two directions still average the cube over the box, so the slab mean of
    # the enhanced region is diluted by (CUBE_LEN/BOXLEN)^2
    dilut = (CUBE_LEN / BOXLEN) ** 2
    ana = eps_bg + (eps_hi - eps_bg) * dilut * (np.abs(off) <= 0.5 * CUBE_LEN)
    ax.plot(ctr, ana, "--", color="crimson", label="analytic", lw=1.4)

    ax.axvline(com[ia], color="C0", alpha=0.5, ls=":", label="centroid")
    ax.axvline(expected[ia], color="crimson", alpha=0.5, ls=":", label="analytic centroid")
    ax.set_xlabel(AXES[ia])
    ax.set_ylabel("mean dust fraction" if ia == 0 else "")
    ax.set_title(r"%s:  $w=%.2f$,  $wt=%.3f$" % (AXES[ia], W_DRIFT[ia], W_DRIFT[ia] * t))
    if ia == 0:
        ax.legend(fontsize=8)
fig.suptitle("DustyAdvect3D   t = %.4f   mass err %.2e   centroid offset [%+.4f %+.4f %+.4f]"
             % (t, err_mass, delta[0], delta[1], delta[2]))
fig.tight_layout()
fig.savefig("dustyadvect3d.pdf", bbox_inches="tight")
fig.savefig("dustyadvect3d.png", dpi=110, bbox_inches="tight")
print("  wrote dustyadvect3d.pdf / .png")

# Keep the harness's bookkeeping happy (it expects a .tex per test).
with open("dustyadvect3d.tex", "w") as fh:
    fh.write(" \n")

failures = []
if err_mass > 1.0e-10:
    failures.append("dust mass not conserved: relative error %.3e" % err_mass)
if eps_min < 0.0:
    failures.append("negative dust density: min = %.6e" % eps_min)
tol_pos = 1.5 * np.min(dx)      # upwind diffuses the edges; the centroid must not shift
if np.any(np.abs(delta) > tol_pos):
    failures.append("centroid off by [%+.4f %+.4f %+.4f], tolerance %.4f"
                    % (delta[0], delta[1], delta[2], tol_pos))

if failures:
    for f in failures:
        print("FAILED: %s" % f)
    print("%-30s %s" % ("dustyadvect3d", "failed"))
    sys.exit(1)

print("%-30s %s" % ("dustyadvect3d", "passed"))
