"""DustyLev: dust levitation -- the absolute normalisation of the IR radiation
force against gravity, plus what an isolated layer can and cannot say about the
trapped-IR channel.

Setup (hydro/condinit.f90 :: dustylev_condinit + dustylev_*.nml): an isolated
Gaussian dusty layer (sigma = boxlen/10, uniform dust mass fraction eps = 1e-2)
in a uniform downward gravity field (gravity_type=1), illuminated from x=0 by a
directed IR beam. Density falls to a floor at BOTH boundaries, so there is no
wall pressure force and the box-integrated momentum budget closes:

    d<v>/dt = g*(f_E - 1),    f_E = a_rad/g,   a_rad = kappa_R*eps*F_E/c

Everything on the right is taken from the code's OWN state: F_E from the dumped
photon flux times the group energy in info_rt (which rt_init recomputes from
the band, so the namelist group_egy is not used), and kappa_R from the
Rosseland table the run itself writes to SEDtables/.

Usage:  python3 plot-dustylev.py [dustylev_bal.nml]
"""
import glob
import io
import os
import re
import shutil
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
    "axes.linewidth": 0.8, "legend.frameon": False, "figure.facecolor": "white",
})
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

_cand = sorted(glob.glob("./ramses_*1d"))
EXEC = _cand[0] if _cand else "./ramses_lev_test1d"
NML = sys.argv[1] if len(sys.argv) > 1 else "dustylev_bal.nml"
TAG = os.path.splitext(os.path.basename(NML))[0].replace("dustylev_", "")
LOG = f"run_{TAG}.log"
DUSTBIN = 1

c_cgs, eV, a_r = 2.99792458e10, 1.602176634e-12, 7.565767e-15


# ------------------------------------------------------- namelist parameters
def nml_val(key, default=None):
    for line in io.open(NML):
        line = line.split("!")[0].strip()
        if "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip().lower() != key.lower():
            continue
        v = rhs.strip().rstrip(",").split(",")[0].strip()
        try:
            return float(v.replace("d", "e").replace("D", "e"))
        except ValueError:
            return default
    return default


boxlen = nml_val("boxlen", 20.0)
levelmin = int(nml_val("levelmin", 7))
d_region = nml_val("d_region", 1.0)
g_code = -nml_val("gravity_params", -0.7142857)     # gravity_params(1) = -g
f_c = nml_val("rt_c_fraction", 1.0e-3)
n_source = nml_val("rt_n_source(1)", 0.0)
scale_l, scale_d, scale_t = nml_val("units_length"), nml_val("units_density"), nml_val("units_time")
scale_v = scale_l / scale_t
scale_E = scale_d * scale_v ** 2
a_grain = nml_val("asize", 0.01) * 1e-4
s_grain = nml_val("sgrain", 2.2)
m_grain = 4.0 / 3.0 * np.pi * a_grain ** 3 * s_grain
eps_lev = 1.0e-2                                     # mirrors condinit
ncell = 2 ** levelmin
dx = boxlen / ncell
g_cgs = g_code * scale_l / scale_t ** 2
# rt/rt_init.f90 sets gamma_rad(1) for the trapped-IR NENER slot. It is now
# 1/3 + 1 (RT15 eq. 46: P = E_t/3 with E_t physical, since Np2Ep already
# applies c_red/c). Before that was fixed it was rt_c_fraction/3 + 1, which
# applied c_red/c a second time. Override for a legacy run.
gamma_rad_m1 = float(os.environ.get("DUSTYLEV_GAMMA_RAD_M1", 1.0 / 3.0))

# ------------------------------------------------------------------ run
# All four variants write to the same output_0000N directories, so a marker
# records which one produced them; otherwise switching namelists silently
# analyses the previous run.
MARK = ".variant"
if os.path.exists(MARK) and io.open(MARK).read().strip() != TAG:
    for d in sorted(os.listdir(".")):
        if d.startswith("output_"):
            shutil.rmtree(d, ignore_errors=True)
if not os.path.isdir("output_00002") or not os.path.exists(LOG):
    print(f"--- Running {EXEC} {NML} (log -> {LOG}) ---")
    for d in sorted(os.listdir(".")):
        if d.startswith("output_"):
            shutil.rmtree(d, ignore_errors=True)
    os.makedirs("SEDtables", exist_ok=True)
    with io.open(LOG, "w") as fh:
        subprocess.run([EXEC, NML], check=True, stdout=fh,
                       stderr=subprocess.STDOUT,
                       env=dict(os.environ, GFORTRAN_UNBUFFERED_ALL="1"))
io.open(MARK, "w").write(TAG + "\n")

# ---- group photon energy: rt_init RECOMPUTES it from the band, so read it back
egy_eV = None
for f in sorted(glob.glob("output_*/info_rt_*.txt")):
    m = re.search(r"egy\s+\[eV\]\s*=\s*([0-9.Ee+-]+)", io.open(f).read())
    if m:
        egy_eV = float(m.group(1))
        break
if egy_eV is None:
    egy_eV = nml_val("group_egy", 0.01)
    print("(WARNING: could not read egy from info_rt; falling back to the namelist)")
egy = egy_eV * eV

# ---- Rosseland mean, from the table the run itself wrote
_ROSS = "SEDtables/rosseland_mean_DustBin_%02d.list" % DUSTBIN
_rT, _rS = np.loadtxt(_ROSS, unpack=True)          # log10(T) , sigma [cm^2/grain]


def kappa_R(T):
    """Rosseland mean cross section per gram of DUST [cm^2/g]."""
    return np.interp(np.log10(np.maximum(T, 1.0)), _rT, _rS) / m_grain


def T_rad_of(Np):
    """The code's radiation temperature: T = (e_gamma*Np*f_c/a_r)^(1/4)."""
    return np.maximum((np.maximum(Np, 0.0) * egy * f_c / a_r) ** 0.25, 1.0)


# -------------------------------------------- read hydro AND rt snapshots
def load(idx):
    h = visu_ramses.load_snapshot(idx)["data"]
    o = np.asarray(h["x"]).argsort()
    out = f"output_{idx:05d}"
    stage = os.path.abspath(f"_rtstage_{idx:05d}")
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(os.path.join(stage, out))
    for f in os.listdir(out):
        src, dst = os.path.join(out, f), os.path.join(stage, out)
        if f.startswith("rt_") and ".out" in f:
            shutil.copy(src, os.path.join(dst, f.replace("rt_", "hydro_", 1)))
        elif f == "rt_file_descriptor.txt":
            shutil.copy(src, os.path.join(dst, "hydro_file_descriptor.txt"))
        elif f.startswith("amr_"):
            shutil.copy(src, os.path.join(dst, f))
        elif f.endswith(".txt") and f != "hydro_file_descriptor.txt" \
                and not f.startswith("info_rt"):
            shutil.copy(src, os.path.join(dst, f))
    cwd = os.getcwd()
    try:
        os.chdir(stage)
        r = visu_ramses.load_snapshot(idx)["data"]
    finally:
        os.chdir(cwd)
        shutil.rmtree(stage, ignore_errors=True)
    ro = np.asarray(r["x"]).argsort()
    pick = lambda d, od: {k: (np.asarray(v)[od]
                              if np.ndim(v) and np.size(v) == np.size(od) else v)
                          for k, v in d.items()}
    return np.asarray(h["x"])[o], pick(h, o), pick(r, ro)


snaps = []
for i in range(1, 40):
    if not os.path.isdir(f"output_{i:05d}"):
        break
    xx, H, R = load(i)
    snaps.append((float(np.atleast_1d(H["time"])[0]), xx, H, R))
if len(snaps) < 4:
    print("FAIL: need at least four snapshots")
    sys.exit(1)

rtc_code = f_c * c_cgs / scale_v


def derived(x, H, R):
    """Per-cell quantities, in cgs unless noted."""
    rho_mix = np.asarray(H["density"]) * scale_d
    eps = np.asarray(H["DustBin_%02d" % DUSTBIN])          # already a fraction
    rho_d = eps * rho_mix
    Np = np.asarray(R["photon_flux_01"]) / rtc_code        # [1/cm3]
    Fp = np.asarray(R["photon_flux_01_x"]) * scale_v       # [1/cm2/s]
    F_E = Fp * egy                                         # [erg/cm2/s]
    T_rad = T_rad_of(Np)
    kR = kappa_R(T_rad)
    chi_R = kR * rho_d                                     # [1/cm]
    a_rad = chi_R * F_E / c_cgs / rho_mix                  # [cm/s2], RT15 eq.27
    P_trap = np.asarray(H["non_thermal_pressure_01"])      # code units
    # cooling_fine builds tau_code = 1.5*dx*chi_R and then f_trap = exp(-1/tau),
    # i.e. exactly RT15 eq. 50, exp(-2/(3*tau_c)) with tau_c = dx*chi_R.
    tau_code = 1.5 * dx * scale_l * chi_R
    f_trap_pred = np.where(tau_code > 0,
                           np.exp(-1.0 / np.maximum(tau_code, 1e-30)), 0.0)
    # what the run ACTUALLY partitioned: E_trap is uold(iIRtrapVar), recovered
    # from the dumped P_trap via gamma_rad(1)-1, and it is already a physical
    # energy density (cooling_fine applies Np2Ep = f_c). The streaming part is
    # egy*Np*f_c, the same conversion.
    E_trap = P_trap / gamma_rad_m1                         # code units
    E_stream = Np * egy * f_c / scale_E                    # code units
    f_trap_ach = E_trap / np.maximum(E_trap + E_stream, 1e-300)
    return dict(rho_mix=rho_mix, eps=eps, rho_d=rho_d, Np=Np, F_E=F_E,
                T_rad=T_rad, kR=kR, chi_R=chi_R, a_rad=a_rad,
                P_trap=P_trap, tau_code=tau_code, f_trap=f_trap_pred,
                f_trap_ach=f_trap_ach, E_trap=E_trap, E_stream=E_stream,
                v=np.asarray(H["velocity_x"]))


# --------------------------------------------------------- momentum budget
t = np.array([s[0] for s in snaps])
M = np.zeros_like(t)
vbar = np.zeros_like(t)
a_pred = np.zeros_like(t)          # <a_rad> = int(chi_R*F_E/c)dx / int(rho)dx
f_trap_max = np.zeros_like(t)
Ptrap_max = np.zeros_like(t)
for k, (tt, x, H, R) in enumerate(snaps):
    D = derived(x, H, R)
    dl = dx * scale_l
    M[k] = np.sum(D["rho_mix"]) * dl
    vbar[k] = np.sum(D["rho_mix"] * D["v"]) / np.sum(D["rho_mix"])
    a_pred[k] = np.sum(D["chi_R"] * D["F_E"] / c_cgs) * dl / M[k]
    f_trap_max[k] = D["f_trap"].max()
    Ptrap_max[k] = D["P_trap"].max()

# The first snapshot(s) contain the light-crossing transient: the beam needs
# boxlen/c_red to fill the box, during which the layer is in free fall. Fit the
# slope only after 3x that, and after the first output.
t_fill = 3.0 * boxlen / rtc_code
mask = (t > max(t_fill, t[1])) & (t > 0)
if mask.sum() < 3:
    mask = t > t[1]
A = np.vstack([t[mask], np.ones(mask.sum())]).T
slope, icept = np.linalg.lstsq(A, vbar[mask], rcond=None)[0]
slope_cgs = slope * scale_v / scale_t

fE_meas = 1.0 + slope / g_code
fE_pred = float(np.mean(a_pred[mask])) / g_cgs
a_pred_code = fE_pred * g_code

print(f"--- DustyLev setup ({TAG}: {NML}) ---")
print(f"  boxlen/ncell     = {boxlen:g} pc / {ncell}   dx = {dx:.5f} pc")
print(f"  g                = {g_code:.7f} code = {g_cgs:.4e} cm/s2")
print(f"  group egy        = {egy_eV:.5g} eV   (from info_rt; rt_init recomputes it)")
print(f"  rho_peak         = {d_region*scale_d:.3e} g/cm3  (n_H ~ {d_region*scale_d/(1.4*1.6726e-24):.2e})")
print(f"  Sigma_gas(t=0)   = {M[0]:.5e} g/cm2")
print(f"  m_grain          = {m_grain:.4e} g")
D0 = derived(*snaps[-1][1:])
_pk = int(np.argmax(D0["rho_mix"]))
print(f"  at the layer peak, t = {t[-1]:.3f}:")
print(f"     F_E           = {D0['F_E'][_pk]:.5e} erg/cm2/s")
print(f"     T_rad         = {D0['T_rad'][_pk]:.4f} K   kappa_R = {D0['kR'][_pk]:.4f} cm2/g_dust")
print(f"     tau_code      = {D0['tau_code'][_pk]:.4e} (= 1.5*tau_c)"
      f"   f_trap: pred {D0['f_trap'][_pk]:.5f}, achieved {D0['f_trap_ach'][_pk]:.5f}")
print(f"  tau_R(layer)     = {np.sum(D0['chi_R'])*dx*scale_l:.4e}")
print()
print(f"  measured d<v>/dt = {slope:+.6f} code/Myr = {slope_cgs:+.4e} cm/s2")
print(f"  predicted        = {a_pred_code - g_code:+.6f} code/Myr"
      f"   [g*(f_E-1), f_E from the code's own chi_R and F_E]")
print(f"  f_E measured     = {fE_meas:.5f}")
print(f"  f_E predicted    = {fE_pred:.5f}")
print(f"  free fall would be {-g_code:+.6f} code/Myr")
print()

# ----------------------------------------------------------------- checks
checks = []


def check(name, ok, detail):
    checks.append((name, bool(ok), detail))


TRAPPING = f_trap_max.max() > 1e-3
Dl = derived(*snaps[-1][1:])

if not TRAPPING:
    # The clean quantitative test: the ABSOLUTE normalisation of the IR
    # radiation force against gravity. Quoted as an error in the net
    # acceleration relative to g -- the meaningful scale, since at f_E = 1 the
    # net acceleration is zero and a relative error on it is undefined.
    err_a = abs(slope - (a_pred_code - g_code)) / g_code
    check("1. net accel = g*(f_E-1)", err_a < 0.02,
          f"|measured-predicted|/g = {err_a:.3e}  (f_E: {fE_meas:.4f} vs {fE_pred:.4f})")
    check("2. no trapping (consistent)", Ptrap_max.max() == 0.0,
          f"max f_trap(pred) = {f_trap_max.max():.3e}, max P_trap = {Ptrap_max.max():.3e}")
else:
    # With trapping on, a per-snapshot prediction of the net acceleration is
    # NOT well posed: the dumped Fp and Np have already been multiplied by
    # (1-f_trap) and de-/re-partitioned relative to the state the momentum
    # transfer inside cool_step actually saw. What IS well posed is a bracket.
    #
    # An isolated layer also gets NO net force from the trapped channel:
    # int(-dP_trap/dx)dx = P_trap(0) - P_trap(L) = 0 because P_trap vanishes in
    # the transparent ambient at both ends. So the net momentum can only come
    # from what still streams, and the measured f_E must sit strictly between
    # free fall (nothing delivered, f_E = 0) and full delivery (f_E = f_E,thin
    # limited by 1-exp(-tau_R), computed from the incident beam).
    F_inc = n_source * egy                                  # incident beam [erg/cm2/s]
    tau_R = float(np.sum(Dl["chi_R"]) * dx * scale_l)
    fE_full = (F_inc / c_cgs) * (1.0 - np.exp(-tau_R)) / (g_cgs * M[-1])
    frac = fE_meas / fE_full
    print(f"  incident beam    = {F_inc:.5e} erg/cm2/s")
    print(f"  f_E if ALL of the absorbed momentum reached the gas = {fE_full:.4f}")
    print(f"  -> the layer receives {100*frac:.1f}% of it "
          f"(f_E measured {fE_meas:.4f}); the rest is momentum the trapping "
          f"closure holds back.")
    print()
    check("1. net accel between free fall and full delivery",
          0.0 < fE_meas < fE_full,
          f"0 < {fE_meas:.4f} < {fE_full:.4f}  ({100*frac:.1f}% delivered)")

    # 2. RT15 eq. 50: the partition itself. f_trap achieved vs predicted, both
    #    from the code's own Rosseland extinction. Only where trapping is
    #    actually significant.
    sel = Dl["f_trap"] > 0.05
    if sel.sum() >= 3:
        r = Dl["f_trap_ach"][sel] / Dl["f_trap"][sel]
        med = float(np.median(r))
        check("2. f_trap follows RT15 eq.50", abs(med - 1.0) < 0.30,
              f"median achieved/predicted = {med:.4f} over {sel.sum()} cells "
              f"(max f_trap {Dl['f_trap'].max():.4f})")
    else:
        check("2. f_trap follows RT15 eq.50", False,
              f"only {sel.sum()} cells with f_trap > 0.05 -- variant mis-sized")

    # 3. the trapped force integral over an isolated layer
    net_trap = (Dl["P_trap"][0] - Dl["P_trap"][-1]) / max(Dl["P_trap"].max(), 1e-300)
    check("3. net trapped force ~ 0 for an isolated layer", abs(net_trap) < 1e-6,
          f"[P_trap(0)-P_trap(L)]/max(P_trap) = {net_trap:.3e} "
          f"(trapped channel only redistributes momentum internally)")

# common
rho_end = Dl["rho_mix"]
edge = max(rho_end[0], rho_end[-1]) / rho_end.max()
check("4. layer clear of the boundaries", edge < 1e-2,
      f"max(rho_edge)/rho_peak = {edge:.3e}")
dM = abs(M[-1] / M[0] - 1.0)
check("5. mass conserved", dM < 5e-3, f"|dM/M| = {dM:.3e}")

# --------------------------------------------- reported, not pass/fail
if TRAPPING:
    sel = Dl["P_trap"] > 1e-3 * Dl["P_trap"].max()
    # Rebuild E_trap INDEPENDENTLY of gamma_rad, from the partition itself:
    # E_trap = f_trap * (total physical IR energy density), with the total
    # recovered from the streaming remainder, Np_stream = (1-f_trap)*Np_tot.
    Np_tot = Dl["Np"][sel] / np.maximum(1.0 - Dl["f_trap_ach"][sel], 1e-30)
    E_ind = Dl["f_trap_ach"][sel] * Np_tot * egy * f_c / scale_E
    ratio = float(np.median(Dl["P_trap"][sel] / np.maximum(E_ind, 1e-300)))
    print("  --- trapped-pressure normalisation cross-check ---")
    print(f"  measured P_trap / E_trap = {ratio:.5e};  RT15 eq. 46 requires "
          f"1/3 = {1/3:.5e};  ratio = {ratio*3:.4f}")
    print(f"  (1 = satisfied. Was {f_c:g} = rt_c_fraction before rt/rt_init.f90's"
          f" gamma_rad(1) was corrected from rt_c_fraction/3+1 to 4/3.)")
    print()

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(2, 2, figsize=(9.6, 6.8))
(A, B), (C, Dp) = ax

# Projections at several times: the layer profile, and the local Eddington
# ratio, at 5 epochs spread over the run.
NT = min(5, len(snaps))
sel_t = np.linspace(0, len(snaps) - 1, NT).astype(int)
cmap = plt.cm.viridis(np.linspace(0.05, 0.9, NT))
for cc, j in zip(cmap, sel_t):
    A.plot(snaps[j][1], np.asarray(snaps[j][2]["density"]), "-", color=cc, lw=1.1,
           label=f"t = {t[j]:.2f}")
A.set_yscale("log")
A.set_xlabel("x [pc]")
A.set_ylabel(r"$\rho$ [code]")
A.set_title(f"A. the dusty layer at {NT} times")
A.legend(fontsize=7)

B.plot(t, vbar, "ko", ms=4, label=r"measured $\langle v\rangle$")
tf = np.linspace(t[mask].min(), t[-1], 50)
B.plot(tf, slope * tf + icept, "C0-", lw=1.2,
       label=rf"fit: $d\langle v\rangle/dt$ = {slope:+.4f}")
B.plot(tf, (a_pred_code - g_code) * (tf - t[mask].min()) + (slope * t[mask].min() + icept),
       "C3--", lw=1.4, label=rf"$g(f_E-1)$ = {a_pred_code-g_code:+.4f}")
B.plot(tf, -g_code * (tf - t[mask].min()) + (slope * t[mask].min() + icept),
       "0.6", ls=":", lw=1.2, label=rf"free fall ($-g$)")
B.axvspan(t[0], max(t_fill, t[1]), color="0.9", zorder=0)
B.set_xlabel("t [Myr]")
B.set_ylabel(r"$\langle v\rangle$ [code]")
B.set_title(f"B. momentum budget   ($f_E$ = {fE_meas:.4f} vs {fE_pred:.4f})")
B.legend(fontsize=7.5, loc="best")

for cc, j in zip(cmap, sel_t):
    Dj = derived(*snaps[j][1:])
    C.plot(snaps[j][1], Dj["a_rad"] / g_cgs, "-", color=cc, lw=1.1)
C.axhline(1.0, color="0.6", ls=":", lw=1.0)
C.axhline(fE_pred, color="C3", ls="--", lw=1.2,
          label=rf"mass-weighted = {fE_pred:.3f}")
C.set_yscale("log")
C.set_xlabel("x [pc]")
C.set_ylabel(r"$f_E$ = $a_{\rm rad}/g$")
C.set_title(rf"C. local Eddington ratio, same {NT} times")
C.legend(fontsize=8)

if TRAPPING:
    sel = Dl["P_trap"] > 0
    Dp.plot(snaps[-1][1], Dl["P_trap"], "C1-", lw=1.2, label=r"$P_{\rm trap}$ [code]")
    Dp.plot(snaps[-1][1], np.asarray(snaps[-1][2]["pressure"]), "0.6", lw=1.0,
            label=r"$P_{\rm gas}$ [code]")
    Dp2 = Dp.twinx()
    Dp2.plot(snaps[-1][1], Dl["f_trap"], "C2--", lw=1.2, label=r"$f_{\rm trap}$")
    Dp2.set_ylabel(r"$f_{\rm trap}$", color="C2")
    Dp2.tick_params(axis="y", colors="C2")
    Dp.set_yscale("log")
    Dp.set_xlabel("x [pc]")
    Dp.set_ylabel("pressure [code]")
    Dp.set_title(r"D. trapped IR: $P_{\rm trap}$ vs $P_{\rm gas}$ and $f_{\rm trap}$")
    Dp.legend(fontsize=8, loc="upper left")
else:
    # RT15 eqs. 85-86 for reference, so the isolated-layer test can be put on
    # the same axes as tests/dust_tva/dustylevatm.
    tauV, tauF = [], []
    for j in range(len(snaps)):
        Dj = derived(*snaps[j][1:])
        dl = dx * scale_l
        tauV.append(np.sum(Dj["chi_R"]) * dl)
        Fm = np.mean(np.asarray(snaps[j][3]["photon_flux_01_x"]))
        tauF.append(np.sum(Dj["chi_R"] * np.asarray(snaps[j][3]["photon_flux_01_x"]))
                    * dl / (len(Dj["chi_R"]) * Fm) if Fm != 0 else np.nan)
    Dp.plot(t, tauV, "C0-o", ms=3, lw=1.1, label=r"$\tau_V$ (RT15 eq. 85)")
    Dp2 = Dp.twinx()
    Dp2.plot(t, np.array(tauF) / np.array(tauV), "C1--s", ms=3, lw=1.1)
    Dp2.set_ylabel(r"$\tau_F/\tau_V$", color="C1")
    Dp2.tick_params(axis="y", colors="C1")
    Dp.set_yscale("log")
    Dp.set_xlabel("t [Myr]")
    Dp.set_ylabel(r"$\tau_V$")
    Dp.set_title(r"D. RT15 optical-depth diagnostics")
    Dp.legend(fontsize=7.5, loc="best")

fig.suptitle(f"DustyLev ({TAG}): IR radiation force vs gravity", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.965])
OUT = f"dustylev_{TAG}.png"
fig.savefig(OUT, dpi=150)
print(f"Saved {OUT}")

# ------------------------------------------------------------------ report
print()
npass = 0
for name, ok, detail in checks:
    print(f"  {'ok  ' if ok else 'FAIL'}  {name:48s} {detail}")
    npass += ok
print(f"\nRESULT: {'PASS' if npass == len(checks) else 'FAIL'}"
      f"  ({npass}/{len(checks)} checks)")
sys.exit(0 if npass == len(checks) else 1)
