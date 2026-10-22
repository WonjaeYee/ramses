#!/usr/bin/env python3
"""Write heI_recombination_cloudy.dat, the He I recombination data RTZ reads from its data directory
(data_dir, next to rad_rec.dat). Everything in it comes from Cloudy's He-like atom:

1. Case-B branching, on (T, n_e). One-zone Cloudy models at constant T,
   'case B 10 no photoionization no pdest'. That makes n^1P -> 1^1S optically thick (tau = 1e10), turns off
   excited-state photoionization and turns off the destruction of line photons by background opacity.
   Collisional transfer between the n=2 levels is included. Per case-B He+ -> He0 recombination
   (R_B = alpha_B n_e n(He+)):
     fM1  = n(2^3S) A(2^3S -> 1^1S) / R_B       the 19.82 eV decays
     f2ph = n(2^1S) A_2ph(2^1S) / R_B           the two-photon decays, with every 2^1P arrival converted
                                                to 2^1S (584 A fully trapped)
     f1P  = n(2^1P) A(2^1P -> 2^1S) / R_B       the arrivals in 2^1P
   RTZ decides per cell which 2^1P arrivals escape as 584 A photons (heI_recrad_cell).
2. The radiative recombination coefficient to He0, total and to 1^1S, from Cloudy's data/he_iso_recomb.dat.
3. The He I two-photon spectrum (2^1S -> 1^1S, 20.616 eV), from the He I row of Cloudy's
   source/atmdat_2photon.cpp: photons above 13.598 eV per decay, their mean energy, and the mean energy of
   the rest.

The values are written with the precision RTZ was validated with (4 decimals for the branching, the 6 of
he_iso_recomb.dat for alpha).

    python make_heI_recombination_table.py --cloudy /path/to/cloudy/source/cloudy.exe [--out FILE] [-j 6]

The Cloudy tree (data/ and source/) is found from the executable's path, or given with --cloudy-root.
Made with Cloudy master e6555102.
"""
import argparse
import os
import re
import subprocess
import tempfile
from multiprocessing import Pool

import numpy as np
from scipy.interpolate import CubicSpline

trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2: trapezoid

T_GRID = [4000.0, 5000.0, 7000.0, 10000.0, 14000.0, 20000.0, 28000.0]
LOGNH_GRID = [0.5 * i for i in range(13)]                # log n_H = 0 ... 6
A_M1 = 1.272e-4      # 2^3S -> 1^1S [s^-1], Cloudy he_transprob.dat
A_2PH = 51.02        # 2^1S -> 1^1S two-photon [s^-1]
A_2058 = 1.976e6     # 2^1P -> 2^1S [s^-1]
LGT_ALPHA = (1.0, 8.0, 0.25)                              # log T range and step of the alpha table
E_2PH, E_H = 20.616, 13.598                               # 2^1S - 1^1S, and the H0 edge [eV]


def he_alpha_rows(root):
    """log10 alpha [cm^3 s^-1] to He0, total and to 1^1S, on Cloudy's log T = 0, 0.25, ..., 10 grid."""
    lines = open(os.path.join(root, "data", "he_iso_recomb.dat")).read().split("\n")
    nlev, ntemp = int(lines[0].split()[1]), int(lines[0].split()[3])
    rows = {}
    for line in lines[1:]:
        p = line.split()
        if len(p) > 3 and p[0] == "1":                    # element 1 = He
            rows[int(p[1])] = np.array([float(x) for x in p[2:2 + ntemp]])
    lgT = 0.25 * np.arange(ntemp)
    return lgT, rows[nlev], rows[0]


def hlike_rows(root):
    """log10 alpha [cm^3 s^-1] of the H-like sequence (data/h_iso_recomb.dat): (H I total, H I to 1s, He II total,
    He II to 1s) on Cloudy's log T = 0, 0.25, ..., 10 grid."""
    lines = open(os.path.join(root, "data", "h_iso_recomb.dat")).read().split("\n")
    nlev, ntemp = int(lines[0].split()[1]), int(lines[0].split()[3])
    rows = {}
    for line in lines[1:]:
        p = line.split()
        if len(p) > 3 and p[0] in ("0", "1"):              # element 0 = H, 1 = He (He II)
            rows[int(p[0]), int(p[1])] = np.array([float(x) for x in p[2:2 + ntemp]])
    return rows[0, nlev], rows[0, 0], rows[1, nlev], rows[1, 0]


def write_hlike(root, out):
    """hlike_recombination_cloudy.dat: the H I and He II recombination of Cloudy's H-like atoms, total and to 1s,
    which split RTZ's recombination emission (ground state / case B) and give its OTSA case-B rates."""
    tot_H, gnd_H, tot_He, gnd_He = hlike_rows(root)
    i0 = int(round(LGT_ALPHA[0] / 0.25))
    n_a = int(round((LGT_ALPHA[1] - LGT_ALPHA[0]) / LGT_ALPHA[2])) + 1
    f = lambda v: " ".join("%.6f" % x for x in v[i0:i0 + n_a])  # noqa: E731
    with open(out, "w") as o:
        o.write("# H-like recombination data for RAMSES-RTZ, from Cloudy's H-like atoms (data/h_iso_recomb.dat,\n"
                "# Cloudy master e6555102). Made by utils/py/rtz/make_heI_recombination_table.py; read by\n"
                "# recombination_module (data_dir). n, log10 T of the first point, step; then log10 alpha [cm^3 s^-1]:\n"
                "# H I total, H I to 1s, He II total, He II to 1s.\n")
        o.write(f"{n_a} {LGT_ALPHA[0]:.2f} {LGT_ALPHA[2]:.2f}\n")
        for v in (tot_H, gnd_H, tot_He, gnd_He):
            o.write(f(v) + "\n")


def two_photon(root):
    """He I two-photon split: photons above E_H per decay, their mean energy, mean energy of the others."""
    src = open(os.path.join(root, "source", "atmdat_2photon.cpp")).read()
    num = lambda s: np.array([float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)])  # noqa: E731
    y = num(re.search(r"He2NuEnergies\[NCRS_HE\]\s*=\s*\{(.*?)\}", src, re.S).group(1))
    table = re.search(r"He2NuShapeFunc\[29\]\[NCRS_HE\]\s*=\s*\{\s*\{(.*?)\}", src, re.S).group(1)
    phi = num(table)[:len(y)]                             # the He I row: photons per unit y, 2 per decay
    spline = CubicSpline(y, phi, bc_type="natural")       # Cloudy's spline (natural end conditions)
    yy = np.linspace(0.0, 1.0, 200001)
    ff = spline(np.minimum(yy, 1.0 - yy))                 # symmetric about y = 1/2
    above = yy >= E_H / E_2PH
    n_tot = trapz(ff, yy)
    n_ion = trapz(ff[above], yy[above])
    e_ion = E_2PH * trapz((yy * ff)[above], yy[above]) / n_ion
    e_low = (E_2PH * trapz(yy * ff, yy) - n_ion * e_ion) / (n_tot - n_ion)
    return n_tot, n_ion, e_ion, e_low


def run_zone(args):
    exe, workdir, lgT_a, lga_tot, lga_gnd, T, lognh = args
    tag = f"T{int(T)}_n{lognh:g}"
    open(os.path.join(workdir, tag + ".in"), "w").write(
        f"title He I case-B branching T={T} log n={lognh}\n"
        "blackbody 5e4\nionization parameter -2\n"
        f"hden {lognh}\nconstant temperature {np.log10(T):.6f}\n"
        "case B 10 no photoionization no pdest\nstop zone 1\nset dr -2\n"
        f'save species densities last "{tag}.pop" "He[1:7]"\n'
        f'save element helium last "{tag}.he"\n'
        f'save physical conditions last "{tag}.phy"\n')
    subprocess.run([exe, "-r", tag], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    pop = np.loadtxt(os.path.join(workdir, tag + ".pop"), comments="#")[1:]
    he = np.loadtxt(os.path.join(workdir, tag + ".he"), comments="#")[1:]
    ne = np.loadtxt(os.path.join(workdir, tag + ".phy"), comments="#")[3]
    nHep = he[1] * pop[0] / he[0]                         # n(He+) = x(He+) n(He0) / x(He0)
    a_tot = 10 ** np.interp(np.log10(T), lgT_a, lga_tot)
    a_gnd = 10 ** np.interp(np.log10(T), lgT_a, lga_gnd)
    RB = (a_tot - a_gnd) * ne * nHep
    return ne, pop[1] * A_M1 / RB, pop[2] * A_2PH / RB, pop[6] * A_2058 / RB


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cloudy", required=True, help="the Cloudy executable")
    ap.add_argument("--cloudy-root", help="the Cloudy tree (data/, source/); default: from the executable's path")
    ap.add_argument("--out", default="heI_recombination_cloudy.dat")
    ap.add_argument("-j", type=int, default=6, help="Cloudy models run at a time")
    ap.add_argument("--hlike-only", metavar="FILE", help="write only hlike_recombination_cloudy.dat (no Cloudy runs)")
    a = ap.parse_args()
    exe = os.path.abspath(a.cloudy)
    root = a.cloudy_root or os.path.dirname(os.path.dirname(os.path.dirname(exe)))   # <root>/source/<sys>/cloudy.exe
    if not os.path.isdir(os.path.join(root, "data")):
        root = os.path.dirname(os.path.dirname(exe))                                   # <root>/source/cloudy.exe
    if a.hlike_only:
        write_hlike(root, a.hlike_only)
        print(f"wrote {a.hlike_only}")
        return
    lgT_a, lga_tot, lga_gnd = he_alpha_rows(root)
    n_tot, n_ion, e_ion, e_low = two_photon(root)

    with tempfile.TemporaryDirectory() as work:
        jobs = [(exe, work, lgT_a, lga_tot, lga_gnd, T, n) for n in LOGNH_GRID for T in T_GRID]
        with Pool(a.j) as pool:
            res = np.array(pool.map(run_zone, jobs)).reshape(len(LOGNH_GRID), len(T_GRID), 4)
    lgne = np.log10(res[:, :, 0]).mean(axis=1)            # one n_e per density (H and He fully ionized)

    i0 = int(round(LGT_ALPHA[0] / 0.25))
    n_a = int(round((LGT_ALPHA[1] - LGT_ALPHA[0]) / LGT_ALPHA[2])) + 1

    f = lambda v, fmt: " ".join(fmt % x for x in v)      # noqa: E731
    with open(a.out, "w") as o:
        o.write("# He I recombination data for RAMSES-RTZ, from Cloudy's He-like atom (Cloudy master e6555102).\n"
                "# Made by utils/py/rtz/make_heI_recombination_table.py; read by recombination_module (data_dir).\n"
                "# 1. Case-B branching per He+ -> He0 case-B recombination, one-zone Cloudy models at constant T,\n"
                "#    'case B 10 no photoionization no pdest': nT, nne; log10 T [K]; log10 n_e [cm^-3]; then fM1\n"
                "#    (2^3S -> 1^1S, 19.82 eV), f2ph (2^1S two-photon, every 2^1P arrival converted) and f1P\n"
                "#    (arrivals in 2^1P), each nne rows of nT values.\n"
                "# 2. Radiative recombination to He0 (data/he_iso_recomb.dat): n, log10 T of the first point, step;\n"
                "#    log10 alpha_total [cm^3 s^-1]; log10 alpha(1^1S).\n"
                "# 3. Two-photon decay 2^1S -> 1^1S (source/atmdat_2photon.cpp, He I): photons above 13.598 eV per\n"
                "#    decay, their mean energy [eV], the mean energy of the others [eV].\n")
        o.write(f"{len(T_GRID)} {len(LOGNH_GRID)}\n")
        o.write(f([np.log10(T) for T in T_GRID], "%.5f") + "\n")
        o.write(f(lgne, "%.4f") + "\n")
        for k in (1, 2, 3):
            for j in range(len(LOGNH_GRID)):
                o.write(f(res[j, :, k], "%.4f") + "\n")
        o.write(f"{n_a} {LGT_ALPHA[0]:.2f} {LGT_ALPHA[2]:.2f}\n")
        o.write(f(lga_tot[i0:i0 + n_a], "%.6f") + "\n")
        o.write(f(lga_gnd[i0:i0 + n_a], "%.6f") + "\n")
        o.write(f"{n_ion:.3f} {e_ion:.2f} {e_low:.2f}\n")
    print(f"wrote {a.out}: two-photon {n_tot:.4f} photons per decay, {n_ion:.4f} above 13.6 eV "
          f"(mean {e_ion:.3f} eV), the others {e_low:.3f} eV")


if __name__ == "__main__":
    main()
