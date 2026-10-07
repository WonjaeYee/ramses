"""DustyIRtrap: verify the trapped-IR contribution to the TVA dust drift.

Setup (hydro/condinit.f90 :: dustyirtrap_condinit and dustyirtrap.nml): uniform
static gas, a LOG RAMP in dust-to-gas ratio spanning ~2 decades, and a LINEAR
trapped-IR radiation pressure ramp held fixed by TVA_TEST_IRTRAP.

With NDUST=1 and NPAH=0 the Rosseland opacity share is s_1 = 1 identically, so
none of the predictions below depend on the dust optical tables:

    grad(P_trap) = -Ptrap_0*Ptrap_dl/boxlen              (constant, negative)
    D_1(rho_d)   = -grad(P_trap) * (1/rho_d - 1/rho_mix)
    t_s(eps)     = sqrt(pi*gamma/8) * s_grain*a_grain / ((1-eps)*rho_mix*c_s)
    w_d          = t_s*D_1 - sum_j eps_j t_s,j D_j

c_s = sqrt(gamma*p_region/rho_mix) independently of eps, because the NENER
energy is removed before the thermal pressure (see the erad subtraction in
calculate_drag_rad_fluxes).

WHY THIS IS A FIRST-STEP CHECK. At step 1 the state is still exactly the
analytic IC, so every quantity has a closed form. Later steps do not, for two
reasons: (a) D_k propto 1/rho_k is a positive feedback (a dust-poor cell drifts
faster and evacuates further) which in a real run is regulated by
tau propto chi_R propto rho_d -> f_trap -> 0, but here P_trap is frozen so that
regulation is absent; (b) the trapped-IR flux F = rho_d*w_d propto t_s*(1-eps)^2 DECREASES with
rho_d, so its characteristic speed dF/drho_d < 0 while w_d > 0; the flux routine
upwinds on sign(w_face), i.e. against the characteristic, which is
anti-diffusive. The domain boundaries seed it: the zero-gradient ghost fill
halves grad(P_trap) in the two edge cells (-5e-6 instead of -1e-5), an O(1)
error in their drift. dustyirtrap_runaway.png shows a 2*dx mode eating inward
from both boundaries while the interior stays smooth at the 1e-4 level.
"""
import io
import os
import subprocess
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "axes.labelsize": 10, "axes.titlesize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.linewidth": 0.8, "legend.frameon": False, "figure.facecolor": "white",
})

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

EXEC, NML, LOG = "./ramses_irtrap_test1d", "dustyirtrap.nml", "run.log"

# ---------------------------------------------------------------- parameters
# must mirror dustyirtrap.nml and dustyirtrap_condinit
gamma, boxlen = 1.4, 10.0
d_region, p_region = 0.204537956, 0.14609854
eps_min, eps_max = 3.0e-3, 3.0e-1
Ptrap_0, Ptrap_dl = 2.0e-4, 0.5
scale_l, scale_d, scale_t = 3.08567758149137e18, 1.0e-24, 3.1556926e13
# read the grain properties from the namelist rather than duplicating them:
# asize must match the optical table shipped for DustBin_01 (0.01 micron), and
# a mismatch here would silently shift t_s.
def _nml(key, default):
    for line in io.open(NML):
        line = line.split("!")[0].strip()
        if "=" in line and line.split("=", 1)[0].strip().lower() == key.lower():
            v = line.split("=", 1)[1].strip().rstrip(",").split(",")[0].strip()
            try:
                return float(v.replace("d", "e").replace("D", "e"))
            except ValueError:
                return default
    return default

a_grain_cm, s_grain = _nml("asize", 0.01) * 1e-4, _nml("sgrain", 2.2)
levelmin = 7
TOL = 1.0e-5          # relative tolerance on the code-vs-analytic comparison

rho_mix = d_region
epstein_coef = np.sqrt(np.pi * gamma / 8.0)
sa_code = (s_grain / scale_d) * (a_grain_cm / scale_l)
grad_Ptrap = -Ptrap_0 * Ptrap_dl / boxlen

# condinit pre-divides the gas pressure by (1-eps), so the TVA solver sees a
# UNIFORM Pg = p_region and grad(P_gas) = 0. The sound speed is then
# c_s = sqrt(gamma*Pg/rho_gas) with rho_gas = (1-eps)*rho_mix, i.e. eps-dependent.
def cs_analytic(rho_d):
    eps = rho_d / rho_mix
    return np.sqrt(gamma * p_region / ((1.0 - eps) * rho_mix))
c_s = np.sqrt(gamma * p_region / rho_mix)   # reference value at eps -> 0

def D_analytic(rho_d):
    return -grad_Ptrap * (1.0 / rho_d - 1.0 / rho_mix)

def ts_analytic(rho_d):
    eps = rho_d / rho_mix
    return epstein_coef * sa_code / ((1.0 - eps) * rho_mix * cs_analytic(rho_d))

def Ptrap_analytic(x):
    return Ptrap_0 * (1.0 - Ptrap_dl * (x / boxlen - 0.5))

def eps_ic(x):
    return eps_min * (eps_max / eps_min) ** (x / boxlen)

print("--- DustyIRtrap analytic setup (code units) ---")
print(f"  c_s          = {c_s:.6f}")
print(f"  grad(P_trap) = {grad_Ptrap:.6e}")
print(f"  eps range    = {eps_min:.3e} .. {eps_max:.3e}")
print(f"  |w_d| range  ~ {ts_analytic(eps_min*rho_mix)*D_analytic(eps_min*rho_mix)/c_s:.3f}"
      f" .. {ts_analytic(eps_max*rho_mix)*D_analytic(eps_max*rho_mix)/c_s:.4f} c_s\n")

# ------------------------------------------------------------------ run
if not os.path.exists("output_00001") or not os.path.exists(LOG):
    print(f"--- Running {EXEC} {NML} (log -> {LOG}) ---")
    env = dict(os.environ, GFORTRAN_UNBUFFERED_ALL="1")
    with io.open(LOG, "w") as fh:
        subprocess.run([EXEC, NML], check=True, stdout=fh,
                       stderr=subprocess.STDOUT, env=env)

# --------------------------------------------- parse the code's own diagnostic
rows = []
for line in io.open(LOG, errors="replace"):
    if "IRTRAP_DIAG" not in line:
        continue
    f = line.replace("=", "= ").split()
    d = {}
    for i, tok in enumerate(f):
        if tok.endswith("=") and i + 1 < len(f):
            try:
                d[tok[:-1]] = float(f[i + 1])
            except ValueError:
                pass
    if d.get("step") == 1.0:
        rows.append(d)

if not rows:
    print("FAIL: no step-1 IRTRAP_DIAG lines in run.log "
          "(need condinit_kind='dustyirtrap' and rt_isIRtrap=.true.)")
    sys.exit(1)

rho_d = np.array([r["rho_d"] for r in rows])
D_code = np.array([r["D"] for r in rows])
ts_code = np.array([r["t_s"] for r in rows])
share = np.array([r["share"] for r in rows])
tsD_code = np.array([r["ts_D"] for r in rows])
gPt_code = np.array([r["grad_Ptrap"] for r in rows])
o = rho_d.argsort()
rho_d, D_code, ts_code, share, tsD_code, gPt_code = (
    a[o] for a in (rho_d, D_code, ts_code, share, tsD_code, gPt_code))

D_ref, ts_ref = D_analytic(rho_d), ts_analytic(rho_d)
err = lambda a, b: np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
e_D, e_ts, e_tsD = err(D_code, D_ref), err(ts_code, ts_ref), err(tsD_code, ts_ref * D_ref)
e_gPt, e_share = err(gPt_code, grad_Ptrap), np.abs(share - 1.0)

# --------------------------------------------------------------- snapshots
snaps = []
for idx in (1, 2, 3):
    try:
        s = visu_ramses.load_snapshot(idx)
    except Exception:
        continue
    q = s["data"]
    oo = q["x"].argsort()
    snaps.append(dict(t=float(np.atleast_1d(q["time"])[0]),
                      x=np.asarray(q["x"])[oo],
                      eps=np.asarray(q["DustBin_01"])[oo],
                      Pt=np.asarray(q["non_thermal_pressure_01"])[oo],
                      pg=np.asarray(q["pressure"])[oo]))

# ------------------------------------------------------------------- plots
# Two separable concerns, checked separately:
#   (1) the difference stencil that produces grad(P_trap). This is the SAME
#       operator the pre-existing gas-pressure term uses, not something the
#       trapped-IR work introduced. Away from the domain edges it must equal the
#       analytic constant; at the outermost cells a one-sided difference is used
#       by design, so those legitimately differ.
#   (2) the trapped-IR driver itself, D = -grad(P_trap)*(s_k/rho_d - 1/rho_mix).
#       Checking D against the code's OWN reported gradient isolates the part
#       this work added, independently of (1).
D_from_code_grad = -gPt_code * (share / rho_d - 1.0 / rho_mix)
e_D_struct = err(D_code, D_from_code_grad)          # concern (2)
interior = np.isclose(gPt_code, grad_Ptrap, rtol=1e-6)
n_edge = int((~interior).sum())

fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.4))
axA, axB, axC, axD = axes.ravel()

# A: frozen trapped-IR pressure profile, numerical vs analytic
xf = np.linspace(0, boxlen, 400)
if snaps:
    s0 = snaps[0]
    axA.plot(xf, Ptrap_analytic(xf), "-", color="#888888", lw=3, alpha=.6,
             label="analytic ramp")
    axA.plot(s0["x"], s0["Pt"], "--", color="#1565c0", lw=1.8, label="t=0 (code)")
    if len(snaps) > 1:
        axA.plot(snaps[-1]["x"], snaps[-1]["Pt"], ":", color="#b71c1c", lw=1.8,
                 label=f"t={snaps[-1]['t']:.3g} (code, frozen)")
    gi = np.gradient(s0["Pt"], s0["x"])[2:-2].mean()
    axA.set_title(f"A. frozen $P_{{\\rm trap}}(x)$\ninterior $\\nabla P$ = {gi:.4e}"
                  f" vs analytic {grad_Ptrap:.4e}")
    axA.set_xlabel("x [pc]"); axA.set_ylabel(r"$P_{\rm trap}$ [code]")
    axA.legend(loc="upper right", fontsize=8); axA.grid(True, ls=":", alpha=.3)

# B: initial dust ramp, numerical vs analytic
if snaps:
    axB.semilogy(xf, eps_ic(xf), "-", color="#888888", lw=3, alpha=.6, label="analytic IC")
    axB.semilogy(s0["x"], s0["eps"], "--", color="#1565c0", lw=1.8, label="t=0 (code)")
    axB.set_title("B. initial dust-to-gas ratio (2-decade ramp)")
    axB.set_xlabel("x [pc]"); axB.set_ylabel(r"$\epsilon_1$")
    axB.legend(loc="lower right", fontsize=8); axB.grid(True, ls=":", alpha=.3)

# C: THE functional check -- D(rho_d) code vs analytic over 2 decades
rr = np.logspace(np.log10(rho_d.min()*.8), np.log10(rho_d.max()*1.2), 300)
axC.loglog(rr, D_analytic(rr), "-", color="#888888", lw=3, alpha=.6,
           label=r"analytic $-\nabla P_{\rm trap}(1/\rho_d-1/\rho_{\rm mix})$")
axC.loglog(rr, -grad_Ptrap/rr, ":", color="#b71c1c", lw=1.4,
           label=r"leading term only ($1/\rho_d$)")
axC.loglog(rho_d[interior], D_code[interior], "o", ms=4.5, mfc="none",
           mec="#1565c0", mew=1.2, label=f"code, step 1 ({int(interior.sum())} cells)")
if n_edge:
    axC.loglog(rho_d[~interior], np.abs(D_code[~interior]), "x", ms=6,
               color="#b71c1c", label=f"domain-edge cells ({n_edge}, |D|)")
axC.set_title("C. drift driver $D_1$ vs dust density\n"
              "(gap to the dotted line = barycentric $-1/\\rho_{\\rm mix}$ term)")
axC.set_xlabel(r"$\rho_d$ [code]"); axC.set_ylabel(r"$D_1$ [code]")
axC.legend(loc="lower left", fontsize=8); axC.grid(True, which="both", ls=":", alpha=.3)

# D: relative errors
axD.semilogy(rho_d, np.maximum(e_D_struct, 1e-17), "o", ms=4, mfc="none", mec="#1565c0",
             label=r"$D_1$ vs code's own $\nabla P_{\rm trap}$")
axD.semilogy(rho_d, np.maximum(e_ts, 1e-17), "s", ms=4, mfc="none", mec="#1b5e20",
             label=r"$t_s$")
axD.semilogy(rho_d[interior], np.maximum(e_D[interior], 1e-17), "^", ms=4,
             mfc="none", mec="#6a1b9a", label=r"$D_1$ vs fully analytic (interior)")
axD.axhline(TOL, color="#b71c1c", ls="--", lw=1.2, label=f"tolerance {TOL:g}")
axD.set_xscale("log"); axD.set_ylim(1e-17, 1.0)
axD.set_title("D. relative error, code vs analytic (step 1)")
axD.set_xlabel(r"$\rho_d$ [code]"); axD.set_ylabel("relative error")
axD.legend(loc="upper left", fontsize=7.5); axD.grid(True, which="both", ls=":", alpha=.3)

fig.suptitle("DustyIRtrap: trapped-IR radiation pressure on TVA dust drift", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("dustyirtrap.png", dpi=170, bbox_inches="tight")
print("Saved dustyirtrap.png")

if len(snaps) > 1:
    f2, a2 = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for sn, col in zip(snaps, ["k", "#1565c0", "#b71c1c"]):
        a2[0].semilogy(sn["x"], sn["eps"], color=col, lw=1.6, label=f"t={sn['t']:.3g}")
        a2[1].plot(sn["x"], sn["pg"], color=col, lw=1.6, label=f"t={sn['t']:.3g}")
    a2[0].set_xlabel("x [pc]"); a2[0].set_ylabel(r"$\epsilon_1$")
    a2[0].set_title("dust evolution, frozen $P_{\\rm trap}$: a $2\\Delta x$ mode\n"
                    "eats inward from BOTH boundaries; interior stays smooth")
    a2[1].axhline(0, color="k", lw=.8)
    a2[1].set_xlabel("x [pc]"); a2[1].set_ylabel(r"$P_{\rm gas}$ [code]")
    a2[1].set_title("gas pressure (stays positive; note $q_{\\rm neul}$ is\n"
                    "pre-divided by $1-\\epsilon$, so the ramp is by construction)")
    for a in a2:
        a.legend(fontsize=8); a.grid(True, ls=":", alpha=.3)
    f2.tight_layout(); f2.savefig("dustyirtrap_runaway.png", dpi=170, bbox_inches="tight")
    print("Saved dustyirtrap_runaway.png")

# ------------------------------------------------------------------ verdict
print(f"\n--- step-1 check over {len(rho_d)} cells, "
      f"rho_d in [{rho_d.min():.3e}, {rho_d.max():.3e}] ---")
print(f"{'quantity':>34} {'max rel.err':>13}   {'verdict':>7}")
ok = True
checks = [
    ("share - 1", e_share, None),
    ("t_s", e_ts, None),
    ("D_1 vs code's own grad (the new term)", e_D_struct, None),
    ("D_1 vs fully analytic (interior only)", e_D[interior], None),
    ("grad(P_trap) (interior only)", e_gPt[interior], None),
]
for name, e, _ in checks:
    m = float(np.max(e)) if np.size(e) else float("nan")
    good = m <= TOL
    ok &= good
    print(f"{name:>34} {m:13.2e}   {'ok' if good else 'FAIL':>7}")
print(f"\ncells using a one-sided difference at the domain edge: {n_edge} "
      f"(expected 2-4 with bound_type=2)")
if n_edge > 4:
    print("FAIL: too many cells deviate from the constant analytic gradient")
    ok = False
if not np.all(D_code[interior] > 0):
    print("FAIL: interior D_1 must be > 0 (negative P_trap gradient drives dust to +x)")
    ok = False
bary = np.abs(D_analytic(rho_d)/(-grad_Ptrap/rho_d) - 1.0).max()
print(f"barycentric -1/rho_mix term contributes up to {100*bary:.1f}% over the sampled range")
if bary < 0.05:
    print("WARN: barycentric term barely exercised; widen eps_max")
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
