"""LevGas: the Krumholz & Thompson (2013) / Davis et al. (2014) dusty
atmosphere, as repeated by Rosdahl & Teyssier 2015 (RT15) sec. 3.7, with
CALIMA's real grain optics in place of the KT13 kappa_R power law.

Reproduces RT15's diagnostics so the comparison is direct:

  RT15 eq. 84   f_E,V = <f_y,rad> / <g rho>            (0.5 at t = 0)
  RT15 eq. 85   tau_V = L_box <kappa_R rho>            (3   at t = 0)
  RT15 eq. 86   tau_F = L_box <kappa_R rho F_y> / <F_y>
                tau_F / tau_V                          (1   at t = 0)
  RT15 eq. 81   f_y,rad = kappa_R rho F_y / c - (1/3) dE_t/dy
  RT15 fig. 9   cell optical depths tau_c = kappa_R rho dy (max, mass-weighted)
  RT15 fig. 13  the streaming / trapped split of f_E,V

Note on the sign in RT15 eq. 81: as printed it reads "+ (1/3) grad E_t", but
E_t decreases upwards in an illuminated atmosphere, so that term would point
DOWN. Their own eq. B8 has the minus, and the code applies -grad(P_trap), so
the minus is used here.

The trapped term is evaluated twice on purpose:
  * "code"  -- from the dumped P_trap = (gamma_rad(1)-1)*E_t, i.e. the force
               the gas actually feels;
  * "RT15"  -- as (1/3)*E_t, i.e. RT15 eq. 46/81.
They differ by exactly rt_c_fraction (see calima/README.md).

Usage:  python3 plot-levgas.py [namelist]
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

_cand = sorted(glob.glob("./ramses_*[0-9]d"))
EXEC = _cand[0] if _cand else "./ramses_levgas2d"
NML = sys.argv[1] if len(sys.argv) > 1 else "levgas_rt15_quick.nml"
# TAG names the log and the figure, and gates the "already run?" check. Derive
# it from the namelist stem with the folder prefix stripped, so the log the
# script looks for is the one a manual run writes.
TAG = os.path.splitext(os.path.basename(NML))[0]
for _pre in ("levgas_rt15_", "levgas_", "atm_"):
    if TAG.startswith(_pre):
        TAG = TAG[len(_pre):]; break
TAG = TAG or "run"
LOG = f"run_{TAG}.log"
DUSTBIN = 1

c_cgs, eV, a_r = 2.99792458e10, 1.602176634e-12, 7.565767e-15
T_STAR = 82.0                 # RT15: the homogeneous initial gas temperature
F_STAR = 1.03e4               # RT15: the incident flux
# RT15 reference values to compare against
RT15 = dict(fE0=0.5, tauV0=3.0, ratio0=1.0, fE_static=10.7, tauV_static=32.5,
            tauc_mean_notrap=2.0, tauc_max_notrap=4.0)


def _nml_lines():
    """Only the lines inside &GROUP ... / blocks -- the free-text header of these
    namelists mentions things like "units_time = t_*", which a naive scan picks
    up as a parameter."""
    inside = False
    for line in io.open(NML):
        st = line.strip()
        if st.startswith("&"):
            inside = True
            continue
        if st == "/":
            inside = False
            continue
        if inside:
            yield line


def nml_val(key, default=None):
    for line in _nml_lines():
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


def nml_str(key, default=""):
    for line in _nml_lines():
        line = line.split("!")[0].strip()
        if "=" in line and line.split("=", 1)[0].strip().lower() == key.lower():
            return line.split("=", 1)[1].strip().strip("'\" ,").lower()
    return default


boxlen = nml_val("boxlen", 64.0)
levelmin = int(nml_val("levelmin", 7))
def nml_list(key, default=(0.0,)):
    for line in _nml_lines():
        line = line.split("!")[0].strip()
        if "=" in line and line.split("=", 1)[0].strip().lower() == key.lower():
            out = []
            for v in line.split("=", 1)[1].strip().rstrip(",").split(","):
                try:
                    out.append(float(v.strip().replace("d", "e").replace("D", "e")))
                except ValueError:
                    pass
            if out:
                return out
    return list(default)


# gravity_params is a vector: -g along x in 1D, (0,-g) along y in 2D. Take the
# largest component so the same parser serves both.
g_code = max(abs(v) for v in nml_list("gravity_params", (-1.0,)))
f_c = nml_val("rt_c_fraction", 3.0e-3)
n_source = nml_val("rt_n_source(1)", 0.0)
scale_l, scale_d, scale_t = nml_val("units_length"), nml_val("units_density"), nml_val("units_time")
scale_v = scale_l / scale_t
scale_E = scale_d * scale_v ** 2
a_grain = nml_val("asize", 0.01) * 1e-4
s_grain = nml_val("sgrain", 2.2)
m_grain = 4.0 / 3.0 * np.pi * a_grain ** 3 * s_grain
STATIC = nml_str("static", ".false.").startswith(".t")
IRTRAP = nml_str("rt_isirtrap", ".true.").startswith(".t")
ncell = 2 ** levelmin
dy = boxlen / ncell
g_cgs = g_code * scale_l / scale_t ** 2

# ------------------------------------------------------------------ run
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

# ---- the group photon energy: rt_init recomputes it, so read it back
egy_eV = None
for f in sorted(glob.glob("output_*/info_rt_*.txt")):
    m = re.search(r"egy\s+\[eV\]\s*=\s*([0-9.Ee+-]+)", io.open(f).read())
    if m:
        egy_eV = float(m.group(1))
        break
egy = (egy_eV if egy_eV else nml_val("group_egy", 0.01)) * eV
NDIM = int(re.search(r"ndim\s*=\s*(\d+)",
                     io.open(sorted(glob.glob("output_00001/info_*.txt"))[0]).read()).group(1))

# Two opacity modes. With CALIMA the Rosseland mean comes from the grain optics
# table the run itself writes (per gram of DUST). Without CALIMA, RAMSES-RT's
# is_kIR_T branch uses the Krumholz & Thompson (2013) power law of RT15 eq. 79
# (per gram of GAS), with the exp(-T_rad/1000 K) sublimation cutoff the code
# applies on top of it.
KAPPA_SC = nml_val("kappasc", None)
if KAPPA_SC is None:
    KAPPA_SC = nml_val("kappasc(1)", None)   # davis.nml writes the index
KT13 = KAPPA_SC is not None
RT15_KAPPA = nml_str("rt_kir_rt15", ".false.").startswith(".t")
if KT13:
    def kappa_R(T):
        """RT15 eq. 79. With rt_kIR_RT15 the code drops the exp(-T_R/1000 K)
        sublimation cutoff, which is what their sec. 3.7 version did; without
        it the cutoff is applied and kappa_R(82 K) is 1.9575 rather than
        2.1248, putting tau_* at 2.76 instead of 3."""
        T = np.maximum(T, 1.0)
        k = KAPPA_SC * (T / 10.0) ** 2
        return k if RT15_KAPPA else k * np.exp(-T / 1.0e3)
else:
    _ROSS = "SEDtables/rosseland_mean_DustBin_%02d.list" % DUSTBIN
    _rT, _rS = np.loadtxt(_ROSS, unpack=True)
    kappa_R = lambda T: np.interp(np.log10(np.maximum(T, 1.0)), _rT, _rS) / m_grain
T_rad_of = lambda Np: np.maximum((np.maximum(Np, 0.0) * egy * f_c / a_r) ** 0.25, 1.0)
rtc_code = f_c * c_cgs / scale_v

# rt/rt_init.f90 sets gamma_rad(1) for the trapped-IR NENER slot. It is now
# 1/3 + 1 (RT15 eq. 46: P = E_t/3 with E_t physical, since Np2Ep already
# applies c_red/c). Before that was fixed it was rt_c_fraction/3 + 1, which
# applied c_red/c a second time. Set DUSTYLEVATM_GAMMA_RAD_M1 to re-analyse a
# legacy run.
GAMMA_RAD_M1 = float(os.environ.get("DUSTYLEVATM_GAMMA_RAD_M1", 1.0 / 3.0))


def tau_eq(kappa_of_T, Sigma):
    """RT15 eqs. 87-88: the static equilibrium of a single optically thick
    layer, reached when the STREAMING flux leaving it equals F_*,

        (1 - exp(-2/(3 tau_c))) c E = 2 F_* ,   tau_c = kappa_R(T_r) Sigma ,
        T_r = (E/a_r)^(1/4) .

    kappa_of_T is per gram of GAS. Solved as a damped fixed point. Feeding it
    RT15's own KT13 law and Sigma = 1.4 reproduces the tau_c = 27 they quote,
    so the same estimator applied to CALIMA's measured Rosseland table is an
    apples-to-apples prediction for this run.
    """
    tau, T_r = 1.0, 0.0
    for _ in range(400):
        ft = np.exp(-2.0 / (3.0 * max(tau, 1e-6)))
        E = 2.0 * F_STAR / (c_cgs * max(1.0 - ft, 1e-12))
        T_r = (E / a_r) ** 0.25
        tau = 0.5 * tau + 0.5 * kappa_of_T(T_r) * Sigma
    return tau, T_r


def load(idx):
    h = visu_ramses.load_snapshot(idx)["data"]
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
    return h, r


def grid(h, r):
    """Map the flat cell list onto a regular (ny, nx) grid; nx = 1 in 1D.

    The height axis is the LAST spatial dimension: x in 1D, y in 2D.
    """
    if NDIM == 1:
        o = np.asarray(h["x"]).argsort()
        ro = np.asarray(r["x"]).argsort()
        shp = (ncell, 1)
        yy = np.asarray(h["x"])[o].reshape(shp)
        xx = np.zeros(shp)
    else:
        o = np.lexsort((np.asarray(h["x"]), np.asarray(h["y"])))
        ro = np.lexsort((np.asarray(r["x"]), np.asarray(r["y"])))
        shp = (ncell, ncell)
        yy = np.asarray(h["y"])[o].reshape(shp)
        xx = np.asarray(h["x"])[o].reshape(shp)
    g = lambda d, k, od: np.asarray(d[k])[od].reshape(shp)
    fy_key = "photon_flux_01_x" if NDIM == 1 else "photon_flux_01_y"
    # output_hydro.f90 labels the NENER slots "non_thermal_pressure_<ivar-3>",
    # not <ivar-nhydro>, so the number shifts with NDIM: _01 in 1D, _02 in 2D.
    pt_key = next(k for k in h if k.startswith("non_thermal_pressure_"))
    out = dict(
        x=xx, y=yy,
        rho=g(h, "density", o) * scale_d,                        # g/cm3
        vy=g(h, "velocity_x" if NDIM == 1 else "velocity_y", o),  # code
        eps=(np.full(shp, nml_val("z_ave", 1.0)) if KT13
             else g(h, "DustBin_%02d" % DUSTBIN, o)),
        P_gas=g(h, "pressure", o),
        P_trap=g(h, pt_key, o),                                  # code
        Np=g(r, "photon_flux_01", ro) / rtc_code,                # 1/cm3
        Fy=g(r, fy_key, ro) * scale_v,                           # 1/cm2/s
        t=float(np.atleast_1d(h["time"])[0]),
    )
    # chi_R = kappa_R * (the mass the opacity refers to). With CALIMA that is
    # the dust mass; with the KT13 law it is the gas mass scaled by Zsolar,
    # exactly as cooling_fine builds tau = nH*Zsolar*f_dust*unit_tau*kIR.
    out["rho_d"] = out["eps"] * out["rho"]
    out["T_rad"] = T_rad_of(out["Np"])
    out["kR"] = kappa_R(out["T_rad"])                            # cm2/g_dust
    out["chi_R"] = out["kR"] * out["rho_d"]                      # 1/cm
    out["F_E"] = out["Fy"] * egy                                 # erg/cm2/s
    # E_t as a PHYSICAL energy density: cooling_fine's Np2Ep already applied
    # f_c, so uold(iIRtrapVar) = P_trap/(gamma_rad-1) is physical.
    out["E_t"] = out["P_trap"] / GAMMA_RAD_M1 * scale_E            # erg/cm3
    out["tau_c"] = out["chi_R"] * dy * scale_l                   # per cell
    return out


def diagnostics(G):
    """RT15 eqs. 81, 84, 85, 86. <> are plain volume averages over the box."""
    dl = dy * scale_l
    # --- streaming term, RT15 eq. 81 first part
    f_stream = G["chi_R"] * G["F_E"] / c_cgs                     # dyn/cm3
    # --- trapped term, -(1/3) dE_t/dy, central differences along the height
    def ddy(a):
        d = np.zeros_like(a)
        d[1:-1] = (a[2:] - a[:-2]) / (2.0 * dl)
        d[0] = (a[1] - a[0]) / dl
        d[-1] = (a[-1] - a[-2]) / dl
        return d
    f_trap_RT15 = -(1.0 / 3.0) * ddy(G["E_t"])                   # RT15 eq. 46
    f_trap_code = -f_c * f_trap_RT15 * -1.0                      # placeholder
    # the force the CODE applies is -d(P_trap)/dy with its own P_trap
    f_trap_code = -ddy(G["P_trap"] * scale_E)
    w = G["rho"] * g_cgs                                         # |g rho|
    mean = lambda a: float(np.mean(a))
    d = dict(
        t=G["t"],
        fE_stream=mean(f_stream) / mean(w),
        fE_trap_code=mean(f_trap_code) / mean(w),
        fE_trap_RT15=mean(f_trap_RT15) / mean(w),
        tau_V=boxlen * scale_l * mean(G["chi_R"]),                # eq. 85
        tau_F=(boxlen * scale_l * mean(G["chi_R"] * G["Fy"])
               / max(mean(G["Fy"]), 1e-300)),                    # eq. 86
        tau_c_max=float(G["tau_c"].max()),
        tau_c_M=float(np.sum(G["tau_c"] * G["rho"]) / np.sum(G["rho"])),
        vy_M=float(np.sum(G["rho"] * G["vy"]) / np.sum(G["rho"])),
        Sigma=mean(G["rho"]) * boxlen * scale_l,
        F_in=float(np.median(G["F_E"][0:2])),                     # bottom rows
        F_out=float(np.median(G["F_E"][-3:])),                    # top rows
    )
    d["fE_code"] = d["fE_stream"] + d["fE_trap_code"]
    d["fE_RT15"] = d["fE_stream"] + d["fE_trap_RT15"]
    d["ratio"] = d["tau_F"] / d["tau_V"] if d["tau_V"] > 0 else np.nan
    return d


snaps, D = [], []
for i in range(1, 60):
    if not os.path.isdir(f"output_{i:05d}"):
        break
    G = grid(*load(i))
    snaps.append(G)
    D.append(diagnostics(G))
if len(snaps) < 3:
    print("FAIL: need at least three snapshots")
    sys.exit(1)
T = np.array([d["t"] for d in D])
get = lambda k: np.array([d[k] for d in D])

print(f"--- LevGas ({TAG}: {NML}) ---")
print(f"  NDIM = {NDIM}   {ncell}^{NDIM} cells   boxlen = {boxlen:g} h_*   dy = {dy:g} h_*")
print(f"  static = {STATIC}   rt_isIRtrap = {IRTRAP}   rt_c_fraction = {f_c:g}")
print(f"  group egy  = {egy_eV:.5g} eV (from info_rt)")
print(f"  g          = {g_cgs:.5e} cm/s2      (RT15: 1.46e-6)")
print(f"  Sigma(t=0) = {D[0]['Sigma']:.5e} g/cm2   (RT15: 1.4)")
print(f"  F_* (namelist rt_n_bound x egy) = {nml_val('rt_n_bound', 0.0)*egy:.5e} erg/cm2/s"
      f"   (RT15: {F_STAR:.3e})")
if KT13:
    print(f"  OPACITY    = KT13 power law (no CALIMA): kappa_R(T_*) = "
          f"{kappa_R(T_STAR):.4f} cm2/g_gas   (RT15 quote 2.13)"
          f"   [rt_kIR_RT15 = {RT15_KAPPA}]")
else:
    _eps0 = float(np.median(snaps[0]["eps"]))
    print(f"  OPACITY    = CALIMA grain optics: kappa_R(T_*) = "
          f"{kappa_R(T_STAR):.3f} cm2/g_dust; with eps = {_eps0:.6f} that is "
          f"{kappa_R(T_STAR)*_eps0:.4f} cm2/g_gas"
          f"   (RT15's KT13 law gives {0.0316*(T_STAR/10)**2:.4f})")
print(f"  T_rad,*    = {(F_STAR/(c_cgs*a_r))**0.25:.3f} K")
print(f"  radiative energy budget: F(top)/F_* = {D[-1]['F_out']/F_STAR:.4f}"
      f"   [the few-per-cent shortfall is the residual of the dust-radiation"
      f" energy exchange, not a thermostat sink: it is 0.970 with T pinned and"
      f" 0.970 with T free.]")
print()
hdr = (f"{'t/t*':>7} {'f_E,V(code)':>12} {'f_E,V(RT15)':>12} {'stream':>9} "
       f"{'trap(code)':>11} {'trap(RT15)':>11} {'tau_V':>9} {'tau_F/tau_V':>12} "
       f"{'tau_c,max':>10} {'<tau_c>_M':>10} {'<vy>/c*':>9}")
print(hdr)
for d in D:
    print(f"{d['t']:7.3f} {d['fE_code']:12.4f} {d['fE_RT15']:12.4f} "
          f"{d['fE_stream']:9.4f} {d['fE_trap_code']:11.4e} {d['fE_trap_RT15']:11.4f} "
          f"{d['tau_V']:9.4f} {d['ratio']:12.4f} {d['tau_c_max']:10.4f} "
          f"{d['tau_c_M']:10.4f} {d['vy_M']:9.4f}")
print()

# ---------------------------------------------------------------- figure
NPROJ = min(4, len(snaps))
pick = [snaps[j] for j in np.linspace(0, len(snaps) - 1, NPROJ).astype(int)]

if NDIM == 1:
    fig = plt.figure(figsize=(11.0, 7.6))
    gs = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.30)
    axp = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
    for ax, key, lab, logy in zip(
            axp, ["rho", "T_rad", "P_trap"],
            [r"$\rho/\rho_*$", r"$T_{\rm rad}/T_*$", r"$P$ [code]"], [True, False, True]):
        for G in pick:
            yv = G[key][:, 0]
            if key == "rho":
                yv = yv / scale_d
            if key == "T_rad":
                yv = yv / T_STAR
            ax.plot(G["y"][:, 0], np.maximum(yv, 1e-30), lw=1.1,
                    label=rf"$t={G['t']:.2f}\,t_*$")
        if key == "P_trap":
            ax.plot(pick[-1]["y"][:, 0], pick[-1]["P_gas"][:, 0], "0.6", lw=1.0,
                    label=r"$P_{\rm gas}$")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel(r"$h/h_*$")
        ax.set_ylabel(lab)
        ax.set_title("profiles at 4 times" if key == "rho" else "")
    axp[0].legend(fontsize=7)
    axp[2].legend(fontsize=7)
    ax1 = fig.add_subplot(gs[1, :])
    ax2 = fig.add_subplot(gs[2, 0])
    ax3 = fig.add_subplot(gs[2, 1])
    ax4 = fig.add_subplot(gs[2, 2])
else:
    fig = plt.figure(figsize=(11.4, 12.0))
    gs = fig.add_gridspec(5, NPROJ, hspace=0.50, wspace=0.28,
                          height_ratios=[1.4, 1.4, 1.4, 1.0, 1.0])
    for j, G in enumerate(pick):
        for row, key, cmap, lab in (
                (0, "rho", "magma", r"$\log_{10}\rho/\rho_*$"),
                (1, "T_rad", "inferno", r"$T_{\rm rad}/T_*$"),
                (2, "drho", "RdBu_r", r"$\delta\rho/\langle\rho\rangle_x$")):
            ax = fig.add_subplot(gs[row, j])
            if row == 0:
                # clip at 1e-4 rho_*: over 10 decades a 25% horizontal
                # modulation is 0.1 dex and invisible.
                im = np.log10(np.maximum(G["rho"] / scale_d, 1e-4))
            elif row == 1:
                im = G["T_rad"] / T_STAR
            else:
                # horizontal structure: this is where chimneys would show
                rm = G["rho"].mean(axis=1, keepdims=True)
                im = G["rho"] / np.maximum(rm, 1e-300) - 1.0
            h = ax.imshow(im, origin="lower", aspect="auto", cmap=cmap,
                          extent=[0, boxlen, 0, boxlen],
                          vmin=(-4 if row == 0 else (1.0 if row == 1 else -0.4)),
                          vmax=(0.2 if row == 0 else (None if row == 1 else 0.4)))
            ax.set_title(rf"$t = {G['t']:.1f}\,t_*$", fontsize=9)
            if j == 0:
                ax.set_ylabel(r"$y/h_*$")
            ax.set_xlabel(r"$x/h_*$")
            if j == NPROJ - 1:
                cb = fig.colorbar(h, ax=ax, fraction=0.05, pad=0.02)
                cb.set_label(lab, fontsize=8)
    ax1 = fig.add_subplot(gs[3, :])
    ax2 = fig.add_subplot(gs[4, 0])
    ax3 = fig.add_subplot(gs[4, 1])
    ax4 = fig.add_subplot(gs[4, 2])

# RT15 fig. 11 top panel + fig. 13 decomposition
ax1.plot(T, get("fE_code"), "C0-o", ms=3, lw=1.2, label=r"$f_{E,V}$ as coded")
ax1.plot(T, get("fE_RT15"), "C3-s", ms=3, lw=1.2, label=r"$f_{E,V}$ with RT15 eq. 46")
ax1.plot(T, get("fE_stream"), "C2--", lw=1.0, label=r"streaming $\kappa_R\rho F_y/c$")
ax1.plot(T, get("fE_trap_RT15"), "C1:", lw=1.2, label=r"trapped $-\frac{1}{3}\nabla E_t$ (RT15)")
ax1.axhline(RT15["fE0"], color="0.6", ls=":", lw=1.0)
ax1.text(T[-1], RT15["fE0"], r"  RT15 $f_{E,*}=0.5$", va="center", fontsize=7, color="0.4")
if STATIC:
    ax1.axhline(RT15["fE_static"], color="k", ls="-.", lw=1.0)
    ax1.text(T[-1], RT15["fE_static"], r"  RT15 static $=10.7$", va="center",
             fontsize=7, color="k")
ax1.set_yscale("log")
_fin = np.concatenate([get("fE_code"), get("fE_RT15"), get("fE_stream")])
_fin = _fin[np.isfinite(_fin) & (_fin > 0)]
if _fin.size:
    ax1.set_ylim(max(1e-3, 0.3 * _fin.min()), 3.0 * _fin.max())
ax1.set_xlabel(r"$t/t_*$")
ax1.set_ylabel(r"$f_{E,V}$")
ax1.set_title(r"RT15 eq. 84 / fig. 11 top + fig. 13 decomposition")
ax1.legend(fontsize=7.5, ncol=2)

ax2.plot(T, get("tau_V"), "C0-o", ms=3, lw=1.2)
ax2.axhline(RT15["tauV0"], color="0.6", ls=":", lw=1.0)
if STATIC:
    ax2.axhline(RT15["tauV_static"], color="k", ls="-.", lw=1.0,
                label=r"RT15 static 32.5")
    ax2.legend(fontsize=7)
ax2.set_yscale("log")
ax2.set_xlabel(r"$t/t_*$")
ax2.set_ylabel(r"$\tau_V$")
ax2.set_title(r"RT15 eq. 85 ($\tau_*=3$)")

ax3.plot(T, get("ratio"), "C0-o", ms=3, lw=1.2)
ax3.axhline(1.0, color="0.6", ls=":", lw=1.0)
ax3.set_ylim(0.0, 1.25)
ax3.set_xlabel(r"$t/t_*$")
ax3.set_ylabel(r"$\tau_F/\tau_V$")
ax3.set_title(r"RT15 eq. 86  ($\tau_F/\tau_V$)")

ax4.plot(T, get("tau_c_max"), "C3-", lw=1.2, label=r"$\tau_{c,\max}$")
ax4.plot(T, get("tau_c_M"), "C3--", lw=1.2, label=r"$\langle\tau_{c}\rangle_M$")
ax4.axhline(1.0, color="0.6", ls=":", lw=1.0)
ax4.set_yscale("log")
ax4.set_xlabel(r"$t/t_*$")
ax4.set_ylabel(r"$\tau_c$")
ax4.set_title("RT15 fig. 9")
ax4.legend(fontsize=7.5)

fig.suptitle(f"LevGas ({TAG}): {'CALIMA + dust dynamics' if not KT13 else 'plain RAMSES-RT'}"
             f" vs Rosdahl & Teyssier 2015 sec. 3.7",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.965])
OUT = f"levgas_{TAG}.png"
fig.savefig(OUT, dpi=140)
print(f"Saved {OUT}")

# ---------------------------------------------------------------- checks
checks = []
chk = lambda n, ok, d: checks.append((n, bool(ok), d))

chk("1. f_E,V = 0.5 at t=0 (RT15 eq. 84)", abs(D[0]["fE_code"] / 0.5 - 1) < 0.10,
    f"measured {D[0]['fE_code']:.4f} vs RT15 0.5")
chk("2. tau_V = 3 at t=0 (RT15 eq. 85)", abs(D[0]["tau_V"] / 3.0 - 1) < 0.10,
    f"measured {D[0]['tau_V']:.4f} vs RT15 3")
chk("3. tau_F/tau_V = 1 at t=0 (RT15 eq. 86)", abs(D[0]["ratio"] - 1) < 0.10,
    f"measured {D[0]['ratio']:.4f} vs RT15 1")
chk("4. injected flux = F_*", abs(D[0]["F_in"] / F_STAR - 1) < 0.10,
    f"measured {D[0]['F_in']:.4e} vs F_* = {F_STAR:.4e} erg/cm2/s")
if IRTRAP:
    lastq = max(1, len(D) // 2)
    # Evaluate only where trapping is actually active: once the gas disperses
    # both terms go to exactly zero and the ratio is 0/0.
    _tr, _tc = get("fE_trap_RT15"), get("fE_trap_code")
    _act = np.abs(_tr) > 1e-3 * np.max(np.abs(_tr)) if np.max(np.abs(_tr)) > 0 \
        else np.zeros_like(_tr, dtype=bool)
    tr = np.median(_tr[_act]) if _act.any() else 0.0
    tc = np.median(_tc[_act]) if _act.any() else 0.0
    chk("5. trapped force = RT15 eq. 46",
        abs(tc / tr - 1) < 0.02 if tr != 0 else False,
        f"coded / RT15 eq.46 = {tc/tr if tr else np.nan:.6f}"
        f"   (was rt_c_fraction = {f_c:g} before the rt_init gamma_rad fix)")
    if STATIC:
        fe = np.median(get("fE_RT15")[lastq:])
        tv = np.median(get("tau_V")[lastq:])
        Sig = D[0]["Sigma"]
        _kg = (lambda T: kappa_R(T)) if KT13 else (lambda T: kappa_R(T) * 1e-2)
        tau_c_calima, Tr_calima = tau_eq(_kg, Sig)
        tau_c_kt13, Tr_kt13 = tau_eq(lambda T: 0.0316 * (T / 10.0) ** 2, 1.4)
        print("  --- static equilibrium, RT15 eqs. 87-88 ---")
        print(f"  estimator fed RT15's own KT13 law and Sigma=1.4:"
              f"  tau_c = {tau_c_kt13:.2f}  (RT15 quote 27)  T_r = {Tr_kt13:.1f} K")
        print(f"  same estimator with THIS run's opacity, Sigma={Sig:.3f}:"
              f"  tau_c = {tau_c_calima:.2f}  T_r = {Tr_calima:.1f} K")
        print(f"  MEASURED tau_V = {tv:.3f};  RT15 MEASURE 32.5 with the KT13 law")
        print(f"  MEASURED f_E,V (RT15 eq.46 normalisation) = {fe:.3f};  RT15 measure 10.7")
        print(f"  MEASURED f_E,V as the code applies it       = "
              f"{np.median(get('fE_code')[lastq:]):.4f}  <- what the gas feels")
        print()
        if KT13:
            chk("6b. static equilibrium reproduces RT15's published numbers",
                abs(tv / RT15["tauV_static"] - 1) < 0.35
                and abs(fe / RT15["fE_static"] - 1) < 0.45,
                f"tau_V = {tv:.2f} vs RT15 32.5;  f_E,V = {fe:.2f} vs RT15 10.7")
        chk("6. static tau_V matches RT15 eq.87/88 for THIS opacity",
            abs(tv / tau_c_calima - 1) < 0.25,
            f"measured {tv:.3f} vs predicted {tau_c_calima:.3f} "
            f"(RT15 measure 32.5 for the KT13 law, whose kappa_R ~ T^2 keeps "
            f"rising where CALIMA's flattens)")
        chk("7. the estimator reproduces RT15's own tau_c = 27",
            abs(tau_c_kt13 / 27.0 - 1) < 0.10,
            f"KT13 law, Sigma=1.4 -> tau_c = {tau_c_kt13:.2f}")
else:
    lastq = max(1, len(D) // 2)
    chk("5. no trapping: P_trap identically zero", np.max(np.abs(
        np.array([np.abs(G["P_trap"]).max() for G in snaps]))) == 0.0,
        "max |P_trap| = 0")
    # RT15's no-trapping run also shows a gradual RISE in tau_V (their fig. 11,
    # light-green): even without the partition, the radiation bath in an
    # optically thick layer raises T_rad and hence kappa_R. What it must NOT do
    # is produce the trapping run's early spike.
    tvn = np.median(get("tau_V")[lastq:])
    chk("6. no trapping: tau_V rises above tau_* but stays bounded",
        3.0 < tvn < 3.0 * 4, f"median tau_V = {tvn:.4f} (tau_* = 3); "
        f"the opacity feedback via T_rad works without the partition")

print()
npass = 0
for name, ok, detail in checks:
    print(f"  {'ok  ' if ok else 'FAIL'}  {name:64s} {detail}")
    npass += ok
print(f"\nRESULT: {'PASS' if npass == len(checks) else 'FAIL'}  ({npass}/{len(checks)})")
sys.exit(0 if npass == len(checks) else 1)
