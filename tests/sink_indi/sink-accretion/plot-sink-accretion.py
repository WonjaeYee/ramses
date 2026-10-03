"""Bondi accretion onto a pre-main-sequence individual-star sink, against the analytic rate.

The sink (seed 0.03 Msun, final mass 300 Msun) sits at rest in uniform, static-temperature gas
(n_H = 100 cm^-3, T = 1e4 K). RAMSES's Bondi rate (Krumholz+2004, pm/sink_particle.f90) is
    Mdot = 4 pi rho_inf (G M)^2 / c^3,   rho_inf = <rho>/alpha(x),   x = ir_cloud dx_min/(2 r_B),
with c^2 = (gamma-1) e_th/rho, so at constant rho and c the mass follows
    M(t) = M0/(1 - A M0 t),   A = 4 pi rho_inf G^2/c^3.
Checks, from the sink files of every output:
  - the sink's acc_rate equals the formula evaluated on the sink's own rho_gas, cs**2 and msink;
  - the mass gained, M(t) - M0, equals the analytic M(t) - M0;
  - the sink stays pre-main-sequence (msink < final_mass, evolution_flag = 1).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import sink_indi_utils as U  # noqa: E402

TOL_RATE = 1e-4   # acc_rate vs formula on the same state
TOL_MASS = 2e-2   # accreted mass vs analytic (the sink file prints 11 digits of msink)
IR_CLOUD = 2

ALPHA_TAB = np.array([
    820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427, 277.295, 237.845, 204.1,
    175.23, 150.524, 129.377, 111.27, 95.7613, 82.4745, 71.0869, 61.3237, 52.9498, 45.7644,
    39.5963, 34.2989, 29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714, 11.4717,
    10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812, 4.94699, 4.44279, 4.00497, 3.6246,
    3.29395, 3.00637, 2.75612, 2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,
    1.62439])


def bondi_alpha(x, ndim=3):
    """rho/rho_inf of the critical Bondi solution at x = r/r_B (bondi_alpha in sink_particle.f90)."""
    xmin, xmax, n = 0.01, 2.0, len(ALPHA_TAB)
    if x <= xmin:
        return np.exp(1.5)/4/np.sqrt(2*x**ndim)
    if x >= xmax:
        return np.exp(1/x)
    i = int(np.floor((n - 1)*np.log(x/xmin)/np.log(xmax/xmin)))
    x0 = np.exp(np.log(xmin) + i*np.log(xmax/xmin)/(n - 1))
    x1 = np.exp(np.log(xmin) + (i + 1)*np.log(xmax/xmin)/(n - 1))
    return ALPHA_TAB[i]*(ALPHA_TAB[i + 1]/ALPHA_TAB[i])**(np.log(x/x0)/np.log(x1/x0))


def bondi_rate(m, rho, c2, dx_min):
    """Mdot [g/s] for sink mass m [g], gas density rho [g/cm3], c2 [cm2/s2], dx_min [cm]."""
    r_b = U.G_CGS*m/c2
    return 4*np.pi*rho/bondi_alpha(IR_CLOUD*0.5*dx_min/r_b)*r_b**2*np.sqrt(c2)


nouts = U.outputs()
S = [U.read_sink_csv(n) for n in nouts]
ul, ud, ut = S[0]["unit_l"], S[0]["unit_d"], S[0]["unit_t"]
um, uv = ud*ul**3, ul/ut
with open(f"output_{nouts[0]:05d}/info_{nouts[0]:05d}.txt") as f:
    info = {k.strip(): v for k, v in (l.split("=", 1) for l in f if "=" in l)}
dx_min = float(info["boxlen"])*0.5**int(info["levelmax"])*ul

t = np.array([s["time"] for s in S])*ut                 # s
m = np.array([s["msink"] for s in S])*um                # g
acc = np.array([s["acc_rate"] for s in S])*um/ut        # g/s
live = np.array([s["cs**2"] > 0 for s in S])            # the first output precedes the first rate
pred = np.array([bondi_rate(s["msink"]*um, s["rho_gas"]*ud, s["cs**2"]*uv**2, dx_min) if ok else np.nan
                 for s, ok in zip(S, live)])

# analytic M(t) from the medium of the first live output (rho and c stay constant)
s1 = S[np.argmax(live)]
rho, c2 = s1["rho_gas"]*ud, s1["cs**2"]*uv**2
A = bondi_rate(1.0, rho, c2, dx_min)                    # Mdot = A M^2 while alpha(x) ~ 1
m0 = m[0]
m_an = m0/(1 - A*m0*(t - t[0]))
dm_sim, dm_an = m - m0, m_an - m0

r_rate = acc[live]/pred[live] - 1
r_mass = dm_sim[1:]/dm_an[1:] - 1
checks = [
    ("acc_rate vs Bondi formula, max |rel. diff|", f"{np.max(np.abs(r_rate)):.2e} (tol {TOL_RATE:g})",
     np.max(np.abs(r_rate)) < TOL_RATE),
    ("accreted mass vs analytic M(t), final", f"{r_mass[-1]:+.2e} (tol {TOL_MASS:g})", abs(r_mass[-1]) < TOL_MASS),
    ("pre-main-sequence throughout (msink < final_mass, flag 1)", "",
     all(s["msink"] < s["final_mass"] and s["evolution_flag"] == 1 for s in S)),
]

tm = (t - t[0])/(1e6*U.YR)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].plot(tm[live], acc[live]*U.YR/U.MSUN, "o", ms=8, mfc="none", label="RAMSES acc_rate")
ax[0].plot(tm, A*m_an**2*U.YR/U.MSUN, "-", lw=2, label=r"$4\pi\rho_\infty G^2M(t)^2/c^3$")
ax[0].set_ylim(0.95*A*m0**2*U.YR/U.MSUN, 1.05*A*m0**2*U.YR/U.MSUN)
ax[0].ticklabel_format(axis="y", useOffset=False)
ax[0].set_xlabel("t [Myr]")
ax[0].set_ylabel(r"$\dot M$ [M$_\odot$ yr$^{-1}$]")
ax[0].set_title("Bondi rate")
ax[1].plot(tm[1:], dm_sim[1:]/U.MSUN, "o", ms=8, mfc="none", label="RAMSES msink - M$_0$")
ax[1].plot(tm, dm_an/U.MSUN, "-", lw=2, label=r"analytic $M_0/(1-AM_0t) - M_0$")
ax[1].set_xlabel("t [Myr]")
ax[1].set_ylabel(r"accreted mass [M$_\odot$]")
ax[1].set_title(f"M$_0$ = {m0/U.MSUN:.4g} M$_\\odot$, n$_H$ = {rho/U.MH:.4g} cm$^{{-3}}$, "
                f"c = {np.sqrt(c2)/1e5:.3g} km/s", fontsize=10)
ax[2].axhspan(-TOL_MASS, TOL_MASS, color="0.9", label=f"mass tolerance $\\pm${TOL_MASS:g}")
ax[2].plot(tm[live], r_rate, "s", ms=8, mfc="none", label="acc_rate / formula - 1")
# the sink file prints msink with 11 digits: +-0.5 in the last one, for both masses in M(t) - M0
err = np.sqrt(2)*0.5*10**(np.floor(np.log10(np.abs(m[1:]/um))) - 10)*um/dm_an[1:]
ax[2].errorbar(tm[1:], r_mass, yerr=err, fmt="o", ms=8, mfc="none", capsize=3,
               label="accreted mass / analytic - 1")
ax[2].axhline(0, color="0.4", lw=1)
ax[2].set_xlabel("t [Myr]")
ax[2].set_ylabel("relative residual")
ax[2].set_title("residuals")
for a in ax:
    a.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("sink-accretion.pdf", bbox_inches="tight")

U.report("sink-accretion", checks)
