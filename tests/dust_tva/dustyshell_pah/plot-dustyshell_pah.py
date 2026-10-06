"""DustyShell_PAH: the drift of the dust AND PAH bins under a uniform flux, against the TVA.

Every species starts as the same Gaussian shell and is pushed by the uniform flux of
tva_test_mode = TVA_TEST_SPRESS (dust_radpressure.f90; no force on the gas, static gas). Each
then drifts at
    w_k = t_s,k D_k - sum_l eps_l t_s,l D_l,   D_k = grad P_g / rho + a_k - sum_l eps_l a_l,
    a_k = sigma_rp,k F / (m_k c),   P_g = (1 - eps) P,
with Epstein t_s,k = sqrt(pi gamma / 8) rho_gr,k a_k / (rho_g c_s). The pressure term matters:
once the 0.1 um shells run ahead, the gradient of eps across the slower shells moves those by up
to ~10% of their radiation drift. The prediction is w_k in every cell of every output (centred
gradient, as the code), weighted by the shell's excess over the floor and integrated over the
outputs (trapezoid). Nothing is calibrated:
- sigma_rp,k is the group average of the CALIMA tables (averaged_cross_section_*Bin_NN.txt) with
  the code's weighting (initialize_cross_sections_from_blackbody_dust_pah: int f lam sigma dlam /
  int f lam dlam, f a 1e5 K Planck function, on 1000 points linear in lambda, with the table
  interpolated in log sigma - log lambda as DustTable%interpolate does);
- the PAH radius and mass follow dust_utils (Nc_to_a, Nc_to_mass), with rho_gr = spah;
- F = rt_n_source group_egy (the SPRESS flux), with the group energy RTZ computes from the source
  spectrum (info_rt; 7.644 eV for 5-10 eV, not the 8 eV of the namelist: the 0.9556 factor of
  plot-dustyshell_mulbin.py); rho_g and c_s come from the output.
In 5-10 eV the neutral and ionised PAH cross sections are equal, so no charge state enters.
The measured drift is the displacement of each shell's centroid (excess over the 1e-8 floor). The
radiation-only drift (no pressure term, eps -> 0) is printed for reference.

Usage: python3 plot-dustyshell_pah.py [dustyshell_pah.nml]   (runs ./ramses_*1d if there is no output)
"""
import glob
import io
import os
import subprocess
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses  # noqa: E402

NML = sys.argv[1] if len(sys.argv) > 1 else "dustyshell_pah.nml"
TOL = 0.02          # |measured / predicted - 1| per species


def nml_arr(key, default=None):
    for line in io.open(NML):
        line = line.split("!")[0].strip()
        if "=" in line and line.split("=", 1)[0].strip().lower() == key.lower():
            out = []
            for v in line.split("=", 1)[1].strip().rstrip(",").split(","):
                v = v.strip().strip("'\"")
                try:
                    out.append(float(v.replace("d", "e").replace("D", "e")))
                except ValueError:
                    out.append(v)
            return out
    return default


if not os.path.isdir("output_00002"):
    exe = sorted(glob.glob("./ramses_*1d"))
    print(f"--- running {exe[0]} {NML} ---")
    with io.open("run.log", "w") as fh:
        subprocess.run([exe[0], NML], check=True, stdout=fh, stderr=subprocess.STDOUT)

c_cgs, h_pl, k_B, eV, amu = 2.99792458e10, 6.62607015e-27, 1.380649e-16, 1.602176634e-12, 1.66053906660e-24
scale_l, scale_d, scale_t = (nml_arr(k)[0] for k in ("units_length", "units_density", "units_time"))
scale_v = scale_l / scale_t
gamma = nml_arr("gamma", [1.4])[0]
E0, E1 = nml_arr("groupL0")[0], nml_arr("groupL1")[0]
Eg = [float(ln.split("=")[1]) for ln in io.open("output_00002/info_rt_00002.txt") if ln.strip().startswith("egy")][0]
F = nml_arr("rt_n_source(1)")[0] * Eg * eV                       # erg cm^-2 s^-1
tables = nml_arr("dust_tables_dir")[0]


def band_cs(path, cols):
    """Group-averaged cross section of the table columns (lambda, C_rp), as the code does."""
    rows = []
    for ln in io.open(path):
        p = ln.replace("|", " ").split()
        if not ln.startswith("#") and len(p) >= max(cols) + 1:
            try:
                rows.append([float(p[i]) for i in cols])
            except ValueError:
                pass
    lam_t, cs_t = np.array(rows).T
    lam = np.linspace(12398.4 / E1, 12398.4 / E0, 1000) * 1e-8
    x = h_pl * c_cgs / (lam * k_B * 1e5)
    f = 2 * h_pl * c_cgs ** 2 / lam ** 5 / np.expm1(np.minimum(x, 700))
    sig = 10 ** np.interp(np.log10(lam * 1e8), np.log10(lam_t), np.log10(cs_t))
    return np.trapezoid(f * lam * sig, lam) / np.trapezoid(f * lam, lam)


def nc_to_nh(nc):
    return int(0.5 * nc + 0.5) if nc <= 25 else int(2.5 * np.sqrt(nc) + 0.5) if nc <= 100 else int(0.25 * nc + 0.5)


species = []                                                      # (name, a [cm], rho_gr, m [g], sigma_rp [cm^2])
for k, nc in enumerate(nml_arr("pah_nc", [])):
    nc = int(nc)
    a = 1e-7 * (nc / 468.0) ** (1 / 3)
    m = (nc * 12.0107 + nc_to_nh(nc) * 1.0080) * amu   # dust_utils Nc_to_mass, constants.f90
    path = os.path.join(tables, f"averaged_cross_section_PAHBin_{k + 1:02d}.txt")
    s_n, s_i = band_cs(path, (0, 3)), band_cs(path, (4, 7))
    if abs(s_i / s_n - 1) > 1e-6:
        print(f"note: PAHBin_{k + 1:02d} neutral and ionised cross sections differ in this group "
              f"({s_n:.4e}, {s_i:.4e}); the prediction takes the neutral one")
    species.append((f"PAHBin_{k + 1:02d}", a, nml_arr("spah")[k], m, s_n))
for k, (asz, sg) in enumerate(zip(nml_arr("asize"), nml_arr("sgrain"))):
    a = asz * 1e-4
    path = os.path.join(tables, f"averaged_cross_section_DustBin_{k + 1:02d}.txt")
    species.append((f"DustBin_{k + 1:02d}", a, sg, 4 / 3 * np.pi * a ** 3 * sg, band_cs(path, (0, 3))))


def snap(i):
    d = visu_ramses.load_snapshot(i)["data"]
    o = np.argsort(d["x"])
    return float(np.atleast_1d(d["time"])[0]), {k: np.asarray(v)[o] for k, v in d.items()
                                                  if np.ndim(v) and np.size(v) == np.size(o)}


nsnap = len([d for d in os.listdir(".") if d.startswith("output_")])
snaps = [snap(i) for i in range(1, nsnap + 1)]
t0, S0 = snaps[0]
t1, S1 = snaps[-1]
x = S1["x"]
eps_floor = 1e-8
names = [n for n, *_ in species]
acc = np.array([s * F / (m * c_cgs) for _, _, _, m, s in species])          # cm s^-2, uniform
geo = np.array([np.sqrt(np.pi * gamma / 8) * rg * a for _, a, rg, _, _ in species])


def centroid(S, name):
    ex = np.maximum(S[name] - eps_floor, 0.0)
    return (ex * S["x"]).sum() / ex.sum()


def drift(S):
    """w_k in every cell [cm/s], (species, cell), as dust_dynamics computes it."""
    eps = np.array([S[n] for n in names])
    et = eps.sum(axis=0)
    rho = S["density"] * scale_d
    P = S["pressure"] * scale_d * scale_v ** 2
    Pg = (1 - et) * P
    c_s = np.sqrt(gamma * Pg / ((1 - et) * rho))
    ts = geo[:, None] / ((1 - et) * rho * c_s)
    D = np.gradient(Pg, S["x"] * scale_l) / rho + acc[:, None] - (eps * acc[:, None]).sum(axis=0)
    return ts * D - (eps * ts * D).sum(axis=0), ts, c_s


def shell_mean(S, w):
    ex = np.maximum(np.array([S[n] for n in names]) - eps_floor, 0.0)
    return (ex * w).sum(axis=1) / ex.sum(axis=1)


tt = np.array([t for t, _ in snaps]) * scale_t
wbar = np.array([shell_mean(S, drift(S)[0]) for _, S in snaps])               # (output, species)
dx_pred = np.trapezoid(wbar, tt, axis=0)                                       # cm
dx_meas = np.array([centroid(S1, n) - centroid(S0, n) for n in names]) * scale_l
w1, ts1, cs1 = drift(S1)
ipk = np.argmax(S1[names[0]])
w_rad = acc * geo / (S1["density"][ipk] * scale_d * cs1[ipk])                # eps -> 0, no pressure
w_pred, w_meas = dx_pred / (tt[-1] - tt[0]), dx_meas / (tt[-1] - tt[0])

print(f"group energy {Eg:.4f} eV, F = {F:.4e} erg/cm2/s, c_s = {cs1[ipk] / 1e5:.4f} km/s, "
      f"{nsnap} outputs, t = {t0:g}-{t1:g}")
print(f"{'species':11s} {'a [um]':>9s} {'sigma_rp [cm2]':>15s} {'t_s [yr]':>10s} {'w_rad [km/s]':>13s} "
      f"{'w_TVA [km/s]':>13s} {'w_meas [km/s]':>14s} {'meas/TVA':>9s}")
ok = True
for k, (n, a, _, _, s) in enumerate(species):
    r = w_meas[k] / w_pred[k]
    ok &= abs(r - 1) < TOL
    print(f"{n:11s} {a * 1e4:9.5f} {s:15.4e} {ts1[k, ipk] / 3.15576e7:10.4g} {w_rad[k] / 1e5:13.5e} "
          f"{w_pred[k] / 1e5:13.5e} {w_meas[k] / 1e5:14.5e} {r:9.4f}")

fig, ax = plt.subplots(figsize=(7.0, 4.4))
cols = ["#0f766e", "#5eead4", "#1565c0", "#b71c1c", "#1b5e20", "#6a1b9a"]
ax.plot(S0["x"], S0[species[0][0]], ":", color="0.4", lw=1.5, label="t = 0 (every species)")
x0 = np.array([centroid(S0, n) for n, *_ in species])
for n, dp, x_0, col in zip(names, dx_pred, x0, cols):
    ax.plot(x, S1[n], "-", color=col, lw=1.8, label=f"{n} (code)")
    xs = x_0 + dp / scale_l
    ax.axvline(xs, color=col, lw=1.0, ls="--")
ax.set_xlabel("x [pc]")
ax.set_ylabel(r"mass fraction $\epsilon_k$")
ax.set_title(f"DustyShell_PAH, t = {t1:g}: shells (lines) and predicted centroids (dashed)", fontsize=10)
ax.legend(fontsize=7, frameon=False, ncol=2)
fig.tight_layout()
fig.savefig("dustyshell_pah.png", dpi=150)
print("saved dustyshell_pah.png")
print(f"{'PASS' if ok else 'FAIL'}: every species drifts at its TVA prediction to within {TOL:.0%}")
sys.exit(0 if ok else 1)
