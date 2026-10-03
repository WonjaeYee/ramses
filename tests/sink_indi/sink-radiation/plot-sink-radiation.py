"""Photon injection by a main-sequence individual-star sink and the HII region it grows.

A 300 Msun MIST star (rt_sink, injected in the central cell) shines into static, pure-H gas
(n_H = 100 cm^-3, T fixed at 1e4 K by rt_Tconst, reduced speed of light c_r = 1e-3 c). The
sink's ionising rate Q(t) is read from the sink table in the log. Analytic references:
  - Stromgren radius R_S = (3 Q/(4 pi alpha n_H^2))^(1/3), cases A and B;
  - the I-front with a finite light speed (Shapiro et al. 2006), integrated with the logged Q(t):
        dR/dt = c_r (Q - 4/3 pi R^3 alpha_B n^2) / (Q + 4 pi R^2 c_r n - 4/3 pi R^3 alpha_B n^2);
  - the ionising-photon budget: int Q dt = N_HII + N_gamma + int alpha_B n_e n_p dV dt, i.e. the
    photons emitted are in flight, have ionised an atom, or were lost to a recombination
    (rt_otsa = .false.: ground-state recombinations are re-emitted, so the net loss is case B).
The front radius is the ionised volume, R_I = (3/(4 pi) sum x_HII dV)^(1/3).

While the front is within a few cells of the source it is not resolved: a cell's light-crossing
time (4.9 kyr) exceeds the recombination time (1.2 kyr), and the explicit M1 (GLF) step moves
photons up to a cell per RT step, so a numerical precursor runs ahead of c_r t at ~dx/dt. The
front is therefore compared in cells at every output, and relatively once c_r t >= R_S.
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import sink_indi_utils as U  # noqa: E402

ALPHA_A, ALPHA_B = 4.18e-13, 2.59e-13   # cm^3/s at 1e4 K (Osterbrock & Ferland 2006)
E_ION = 13.6
TOL_RS = 0.10        # late-time front radius vs the case-B Stromgren radius
TOL_FRONT = 0.15     # front radius vs the finite-light-speed solution, once c_r t >= R_S
TOL_FRONT_DX = 2.0   # ... and within this many cells at every output (see below)
TOL_BUDGET = 0.05    # ionising-photon budget closure
T_LATE_KYR = 100.0


def front_ode(t_eval, t_q, q, n, c_r, alpha, nstep=20000):
    """I-front radius [cm] at t_eval [s] for the rate q(t_q) [1/s], RK4 from R = 0."""
    ts = np.linspace(0.0, t_eval.max(), nstep + 1)

    def f(t, r):
        qt = np.interp(t, t_q, q)
        rec = 4/3*np.pi*r**3*alpha*n*n
        return c_r*max(qt - rec, 0.0)/(qt - rec + 4*np.pi*r*r*c_r*n)

    r, out, h = 0.0, [0.0], ts[1] - ts[0]
    for t in ts[:-1]:
        k1 = f(t, r)
        k2 = f(t + h/2, r + h/2*k1)
        k3 = f(t + h/2, r + h/2*k2)
        k4 = f(t + h, r + h*k3)
        r += h/6*(k1 + 2*k2 + 2*k3 + k4)
        out.append(r)
    return np.interp(t_eval, ts, out)


log = U.read_log(U.find_log("sink-radiation"))
if log["t_yr"].size == 0:
    print("No sink table found in run.log or ../../test_suite.log")
    print("Test sink-radiation: FAILED")
    sys.exit()
t_log = log["t_yr"]*U.YR
q_ion = log["Q_ion"]
dt_log = np.diff(np.concatenate([[0.0], t_log]))
q_cum = np.cumsum(q_ion*dt_log)    # the rate printed after a step was used during that step

nouts = [n for n in U.outputs() if n > 1]
with open(f"output_{nouts[0]:05d}/info_rt_{nouts[0]:05d}.txt") as f:
    l0 = np.array(re.search(r"groupL0\s+\[eV\]\s+=(.*)", f.read()).group(1).split(), float)
ion_groups = [g + 1 for g in range(len(l0)) if l0[g] >= E_ION - 1e-6]

rows, snaps = [], {}
for n in nouts:
    d = U.load(n)
    c, info = d["data"], d["info"]
    V, x, nH = c["vol_cm3"], c["H_02"], c["nH"]
    t = info["time"]*info["unit_t"]
    rows.append(dict(
        t=t, R=(3*np.sum(x*V)/(4*np.pi))**(1/3), n=np.sum(nH*V)/np.sum(V),
        N_HII=np.sum(x*nH*V), N_g=sum(np.sum(c[f"Np_{g}"]*V) for g in ion_groups),
        rec_B=np.sum(ALPHA_B*(x*nH)**2*V), Q=np.interp(t, t_log, q_ion), Q_cum=np.interp(t, t_log, q_cum)))
    snaps[n] = d
c_r = info["c_red"]
t = np.array([r["t"] for r in rows])
R = np.array([r["R"] for r in rows])
n_h = np.mean([r["n"] for r in rows])
Q = np.array([r["Q"] for r in rows])
rec_B = np.array([r["rec_B"] for r in rows])
# photons lost to recombinations: trapezoid over the outputs, from zero at t = 0
lost = np.concatenate([[0.0], np.cumsum(0.5*(rec_B[1:] + rec_B[:-1])*np.diff(t))]) + 0.5*rec_B[0]*t[0]
budget = (np.array([r["N_HII"] + r["N_g"] for r in rows]) + lost)/np.array([r["Q_cum"] for r in rows])

r_s = {k: (3*Q/(4*np.pi*a*n_h**2))**(1/3) for k, a in (("A", ALPHA_A), ("B", ALPHA_B))}
t_fine = np.linspace(0.0, t.max(), 400)
front_B = front_ode(t_fine, t_log, q_ion, n_h, c_r, ALPHA_B)
front_A = front_ode(t_fine, t_log, q_ion, n_h, c_r, ALPHA_A)
front_at = np.interp(t, t_fine, front_B)
late = t/U.KYR >= T_LATE_KYR
resolved = c_r*t >= r_s["B"]
dx = snaps[nouts[-1]]["data"]["dx"].min()*snaps[nouts[-1]]["info"]["unit_l"]

checks = [
    (f"front radius / R_S(case B), t >= {T_LATE_KYR:g} kyr, worst",
     f"{(R[late]/r_s['B'][late] - 1)[np.argmax(np.abs(R[late]/r_s['B'][late] - 1))]:+.3f} (tol {TOL_RS:g})",
     np.all(np.abs(R[late]/r_s["B"][late] - 1) < TOL_RS)),
    ("front radius / finite-c_r solution, c_r t >= R_S, worst",
     f"{(R/front_at - 1)[resolved][np.argmax(np.abs(R/front_at - 1)[resolved])]:+.3f} (tol {TOL_FRONT:g})",
     resolved.any() and np.all(np.abs(R/front_at - 1)[resolved] < TOL_FRONT)),
    ("front radius - finite-c_r solution, every output, worst [cells]",
     f"{((R - front_at)/dx)[np.argmax(np.abs(R - front_at))]:+.2f} (tol {TOL_FRONT_DX:g})",
     np.all(np.abs(R - front_at) < TOL_FRONT_DX*dx)),
    ("ionising-photon budget (N_HII + N_gamma + lost)/int Q dt, worst",
     f"{(budget - 1)[np.argmax(np.abs(budget - 1))]:+.3f} (tol {TOL_BUDGET:g})",
     np.all(np.abs(budget - 1) < TOL_BUDGET)),
]

tk = t/U.KYR
fig, ax = plt.subplots(2, 2, figsize=(12, 9.5))
# (a) neutral fraction in the mid-plane at the last output, with R_S
d = snaps[nouts[-1]]
c = d["data"]
xs = np.array([d["sinks"]["x"][0], d["sinks"]["y"][0], d["sinks"]["z"][0]])
mid = np.abs(c["z"] - xs[2]) < 0.5*c["dx"]
ul = d["info"]["unit_l"]/U.PC
ext = d["info"]["boxlen"]*ul
nx = int(round(d["info"]["boxlen"]/c["dx"][mid][0]))
img = np.full((nx, nx), np.nan)
ix = np.floor(c["x"][mid]/c["dx"][mid]).astype(int)
iy = np.floor(c["y"][mid]/c["dx"][mid]).astype(int)
img[iy, ix] = c["H_01"][mid]
im = ax[0, 0].imshow(np.log10(np.maximum(img, 1e-8)), origin="lower", extent=(0, ext, 0, ext), cmap="Greys",
                     vmin=-6, vmax=0)
fig.colorbar(im, ax=ax[0, 0], label=r"log$_{10}$ x$_{\rm HI}$")
th = np.linspace(0, 2*np.pi, 200)
for k, ls in (("B", "-"), ("A", "--")):
    rr = r_s[k][-1]/U.PC
    ax[0, 0].plot(xs[0]*ul + rr*np.cos(th), xs[1]*ul + rr*np.sin(th), ls, color="C3", lw=1.5,
                  label=f"R$_S$ case {k}")
ax[0, 0].set_xlabel("x [pc]")
ax[0, 0].set_ylabel("y [pc]")
ax[0, 0].set_title(f"mid-plane, t = {tk[-1]:.0f} kyr")
ax[0, 0].legend(frameon=False, fontsize=9, loc="upper right", labelcolor="w")
# (b) front radius
ax[0, 1].plot(t_fine/U.KYR, front_B/U.PC, "-", color="C0", lw=2, label=r"finite c$_r$, case B")
ax[0, 1].plot(t_fine/U.KYR, front_A/U.PC, "--", color="C0", lw=1.5, label=r"finite c$_r$, case A")
ax[0, 1].plot(t_fine/U.KYR, c_r*t_fine/U.PC, ":", color="0.5", lw=1.5, label=r"c$_r$ t")
ax[0, 1].plot(tk, r_s["B"]/U.PC, "-", color="C3", lw=1, label=r"R$_S$(Q(t)), case B")
ax[0, 1].errorbar(tk, R/U.PC, yerr=TOL_FRONT_DX*dx/U.PC, fmt="o", color="C1", ms=8, mfc="none", capsize=3,
                  label=f"RAMSES ($\\pm${TOL_FRONT_DX:g} cells)")
ax[0, 1].set_ylim(0, 1.3*r_s["B"].max()/U.PC)
ax[0, 1].set_xlabel("t [kyr]")
ax[0, 1].set_ylabel(r"R$_I$ [pc]")
ax[0, 1].set_title(f"Q$_{{\\rm ion}}$ = {q_ion.min():.2e}-{q_ion.max():.2e} s$^{{-1}}$, n$_H$ = {n_h:.4g} cm$^{{-3}}$",
                   fontsize=10)
ax[0, 1].legend(frameon=False, fontsize=9, loc="lower right")
# (c) ionised-fraction profiles
for k, n in enumerate(nouts[::2]):
    c = snaps[n]["data"]
    edges = np.arange(0, c["r_pc"].max() + 1.5, 1.5)
    idx = np.digitize(c["r_pc"], edges)
    prof = [c["H_02"][idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(edges))]
    tn = snaps[n]["info"]["t_kyr"]
    ax[1, 0].plot(0.5*(edges[1:] + edges[:-1]), prof, "-", color=plt.cm.viridis(k/max(1, len(nouts[::2]) - 1)),
                  lw=2, label=f"{tn:.0f} kyr")
    ax[1, 0].axvline(np.interp(tn*U.KYR, t_fine, front_B)/U.PC, color=plt.cm.viridis(k/max(1, len(nouts[::2]) - 1)),
                     lw=1, ls=":")
ax[1, 0].set_xlim(0, 24)
ax[1, 0].set_xlabel("r [pc]")
ax[1, 0].set_ylabel(r"x$_{\rm HII}$ (shell mean)")
ax[1, 0].set_title("profiles; dotted: finite-c$_r$ front")
ax[1, 0].legend(frameon=False, fontsize=9)
# (d) photon budget
ax[1, 1].axhspan(1 - TOL_BUDGET, 1 + TOL_BUDGET, color="0.9", label=f"tolerance $\\pm${TOL_BUDGET:g}")
ax[1, 1].plot(tk, budget, "o", ms=8, mfc="none", label=r"(N$_{\rm HII}$ + N$_\gamma$ + lost)/$\int$Q dt")
ax[1, 1].plot(tk, rec_B/Q, "s", ms=8, mfc="none", label=r"recombinations$_B$/Q")
ax[1, 1].plot(tk, [r["N_g"]/r["Q_cum"] for r in rows], "^", ms=8, mfc="none", label=r"N$_\gamma$/$\int$Q dt")
ax[1, 1].axhline(1, color="0.4", lw=1)
ax[1, 1].set_xlabel("t [kyr]")
ax[1, 1].set_ylabel("ratio")
ax[1, 1].set_title("ionising-photon budget")
ax[1, 1].legend(frameon=False, fontsize=9, loc="center right")
fig.tight_layout()
fig.savefig("sink-radiation.pdf", bbox_inches="tight")

U.report("sink-radiation", checks)
