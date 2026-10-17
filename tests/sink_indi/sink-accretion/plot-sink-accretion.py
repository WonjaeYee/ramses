"""Bondi growth of a pre-main-sequence individual-star sink, and its switch to the main sequence.

A 2 Msun seed (final mass 20 Msun) sits at rest in uniform gas at fixed temperature (n_H = 1.8e6
cm^-3, T = 1e3 K, rt_Tconst). RAMSES's Bondi rate (Krumholz+2004, pm/sink_particle.f90) is
    Mdot = 4 pi rho_inf (G M)^2 / c^3,   rho_inf = <rho>/alpha(x),   x = ir_cloud dx_min/(2 r_B),
with c^2 = (gamma-1) e_th/rho. It is a subgrid rate: gravity enters only through r_B = G M/c^2, so
the test runs without Poisson, which is consistent while r_B << r_acc (here r_B < r_acc/70 and
alpha - 1 < 1.4%). In a uniform medium the mass follows t(M) = int_M0^M dM'/Mdot(M'), close to
M0/(1 - A M0 t). Checks:
  - the sink's acc_rate equals the formula on the sink's own rho_gas, cs**2 and msink (outputs);
  - the growth: the time at which the logged mass reaches M, against the analytic t(M);
  - the switch: the sink reaches the main sequence at the analytic t(M_F), with msink within 1% of
    the final mass, then stops accreting and, within 10 kyr, emits ionising photons (MIST).
The growth and the photon rates are read from the sink table in the log (run.log, or this test's
part of tests/test_suite.log).
"""
import os
import sys

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import sink_indi_utils as U  # noqa: E402

TOL_RATE = 1e-4     # acc_rate vs formula on the same state
TOL_TIME = 2e-2     # t(M) vs analytic, in units of the analytic time to the final mass
TOL_OVER = 1e-2     # msink overshoot of the final mass at the switch
T_ON_KYR = 10.0     # ionising photons must start within this time of the switch: the MIST rates
                    # are zero below the first age node and ramp up over ~60 kyr for 20 Msun
IR_CLOUD = 2
MSUN_CODE = 1.9891e33   # M_sun in amr/constants.f90, used for the log's Msun columns

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


L = U.read_log(U.find_log("sink-accretion"))
if L["t_yr"].size == 0:
    print("No sink table found in run.log or ../../test_suite.log")
    print("Test sink-accretion: FAILED")
    sys.exit()

# --- outputs: the rate check and the medium -------------------------------------------------------
nouts = U.outputs()
S = [U.read_sink_csv(n) for n in nouts]
ul, ud, ut = S[0]["unit_l"], S[0]["unit_d"], S[0]["unit_t"]
um, uv = ud*ul**3, ul/ut
with open(f"output_{nouts[0]:05d}/info_{nouts[0]:05d}.txt") as f:
    info = {k.strip(): v for k, v in (l.split("=", 1) for l in f if "=" in l)}
dx_min = float(info["boxlen"])*0.5**int(info["levelmax"])*ul

t_out = np.array([s["time"] for s in S])*ut
m_out = np.array([s["msink"] for s in S])*um
pre = np.array([s["cs**2"] > 0 and s["evolution_flag"] == 1 for s in S])   # rate computed, pre-MS
acc = np.array([s["acc_rate"] for s in S])*um/ut
pred = np.array([bondi_rate(s["msink"]*um, s["rho_gas"]*ud, s["cs**2"]*uv**2, dx_min) if ok else np.nan
                 for s, ok in zip(S, pre)])
s1 = S[np.argmax(pre)]
rho, c2 = s1["rho_gas"]*ud, s1["cs**2"]*uv**2
m0 = m_out[0]
m_f = S[0]["final_mass"]*um

# --- analytic growth: t(M) by quadrature, with alpha(x(M)) -----------------------------------------
m_grid = np.geomspace(m0, m_f, 4000)
dtdm = 1/np.array([bondi_rate(m, rho, c2, dx_min) for m in m_grid])
t_an = np.concatenate([[0.0], np.cumsum(0.5*(dtdm[1:] + dtdm[:-1])*np.diff(m_grid))])
t_switch_an = t_an[-1]

# --- the growth as logged every step ---------------------------------------------------------------
t_g, m_g = L["t_yr"]*U.YR, L["M"]*MSUN_CODE
mdot_g, q_g, tag = L["Mdot"]*MSUN_CODE/U.YR, L["Q_ion"], L["tag"]
grow = (tag == 1) & (m_g <= (1 - 1e-3)*m_f)
r_time = (t_g[grow] - np.interp(m_g[grow], m_grid, t_an))/t_switch_an
ms = tag == 0
t_switch = t_g[np.argmax(ms)] if ms.any() else np.nan
m_switch = m_g[np.argmax(ms)] if ms.any() else np.nan
lit = q_g > 0
t_on = (t_g[np.argmax(lit)] - t_switch)/U.KYR if lit.any() else np.nan
on = ms & (t_g - t_switch >= T_ON_KYR*U.KYR)

r_rate = acc[pre]/pred[pre] - 1
checks = [
    ("acc_rate vs Bondi formula, max |rel. diff|", f"{np.max(np.abs(r_rate)):.2e} (tol {TOL_RATE:g})",
     np.max(np.abs(r_rate)) < TOL_RATE),
    ("growth: t(M) vs analytic, worst, per t(M_F)", f"{r_time[np.argmax(np.abs(r_time))]:+.2e} (tol {TOL_TIME:g})",
     np.all(np.abs(r_time) < TOL_TIME)),
    ("reaches the main sequence", f"at t = {t_switch/U.YR/1e6:.4f} Myr" if ms.any() else "never", ms.any()),
    ("switch time vs analytic t(M_F)", f"{t_switch/t_switch_an - 1:+.2e} (tol {TOL_TIME:g})",
     ms.any() and abs(t_switch/t_switch_an - 1) < TOL_TIME),
    ("msink at the switch / final mass - 1", f"{m_switch/m_f - 1:+.2e} (in [0, {TOL_OVER:g}))",
     ms.any() and 0 <= m_switch/m_f - 1 < TOL_OVER),
    ("main sequence: no accretion", f"max Mdot {np.max(mdot_g[ms]) if ms.any() else np.nan:.1e} g/s",
     ms.any() and np.all(mdot_g[ms] == 0)),
    (f"ionising photons: none before the switch, on within {T_ON_KYR:g} kyr and after",
     f"on {t_on:.1f} kyr after the switch, Q_ion at the end {q_g[-1]:.2e} s^-1",
     on.any() and np.all(q_g[~ms] == 0) and np.all(q_g[on] > 0)),
]

# --- figure ---------------------------------------------------------------------------------------
myr = 1e6*U.YR
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
ax[0].plot(t_an/myr, m_grid/U.MSUN, "-", color="C1", lw=2.5, label=r"analytic $t(M)=\int dM/\dot M$")
ax[0].plot(t_g/myr, m_g/U.MSUN, "-", color="C0", lw=1, label="RAMSES (log, every step)")
ax[0].plot(t_out/myr, m_out/U.MSUN, "o", color="C0", ms=8, mfc="none", label="RAMSES outputs")
ax[0].axvline(t_switch_an/myr, color="0.4", ls=":", lw=1.5, label="analytic t(M$_F$)")
if ms.any():
    ax[0].axvspan(t_switch/myr, t_g.max()/myr, color="0.92", label="main sequence (RAMSES)")
ax[0].set_yscale("log")
ax[0].set_xlabel("t [Myr]")
ax[0].set_ylabel(r"M$_{\rm sink}$ [M$_\odot$]")
ax[0].set_title(f"n$_H$ = {rho/U.MH:.3g} cm$^{{-3}}$, c = {np.sqrt(c2)/1e5:.3g} km/s, "
                f"M$_F$ = {m_f/U.MSUN:.3g} M$_\\odot$", fontsize=10)
mm = np.geomspace(m0, m_f, 200)
ax[1].plot(mm/U.MSUN, [bondi_rate(m, rho, c2, dx_min)*U.YR/U.MSUN for m in mm], "-", color="C1", lw=2.5,
           label=r"$4\pi\rho_\infty G^2M^2/c^3$")
ax[1].plot(m_g[tag == 1]/U.MSUN, mdot_g[tag == 1]*U.YR/U.MSUN, "-", color="C0", lw=1, label="RAMSES (log)")
ax[1].plot(m_out[pre]/U.MSUN, acc[pre]*U.YR/U.MSUN, "o", color="C0", ms=8, mfc="none", label="RAMSES acc_rate")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_xlabel(r"M$_{\rm sink}$ [M$_\odot$]")
ax[1].set_ylabel(r"$\dot M$ [M$_\odot$ yr$^{-1}$]")
ax[1].set_title("Bondi rate vs mass (pre-main-sequence)")
ax[2].axhspan(-TOL_TIME, TOL_TIME, color="0.9", label=f"tolerance $\\pm${TOL_TIME:g}")
ax[2].plot(m_g[grow]/U.MSUN, r_time, "-", color="C0", lw=1.5, label=r"(t - t$_{an}$(M)) / t$_{an}$(M$_F$)")
ax[2].plot(m_out[pre]/U.MSUN, r_rate, "s", color="C2", ms=8, mfc="none", label="acc_rate / formula - 1")
ax[2].axhline(0, color="0.4", lw=1)
ax[2].set_xscale("log")
ax[2].set_xlabel(r"M$_{\rm sink}$ [M$_\odot$]")
ax[2].set_ylabel("relative residual")
ax[2].set_title("residuals")
ticks = [m for m in (2, 3, 5, 10, 20, 50, 100, 300) if m0/U.MSUN*0.99 <= m <= m_f/U.MSUN*1.01]
for a in ax[1:]:
    a.set_xticks(ticks, [f"{m:g}" for m in ticks])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].set_yticks(ticks, [f"{m:g}" for m in ticks])
ax[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for a in ax:
    a.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("sink-accretion.pdf", bbox_inches="tight")

U.report("sink-accretion", checks)
