#!/usr/bin/env python3
"""Read and check the grain-optics tables CALIMA writes at start-up (init_dust_mean_cross_sections).

RAMSES writes them on rank 1 into sed_dir (the namelist's SED directory, or $RAMSES_SED_DIR; SEDtables/ in
the tests), two columns per line:

  dust bin NN (DustBin_NN)
    rosseland_mean_DustBin_NN.list   log10 T [K]        Rosseland-mean absorption cross section [cm^2/grain]
    planck_mean_DustBin_NN.list      log10 T [K]        Planck-mean absorption cross section    [cm^2/grain]
    planck_power_DustBin_NN.list     log10 T [K]        log10 P_emit = log10(4 sigma_SB <C_abs>_P T^4) [erg/s/grain]
    planck_deriv_DustBin_NN.list     log10 T [K]        log10 dP_emit/dT                       [erg/s/K/grain]
    tdust_tab_DustBin_NN.list        log10 P_emit       log10 T [K]  (the inverse table: T_dust of an absorbed power)
  PAH bin NN (PAHBin_NN), each charge state, _n neutral and _i ionised
    rosseland_mean_PAHBin_NN_{n,i}.list, planck_mean_PAHBin_NN_{n,i}.list   log10 T [K]  cross section [cm^2/PAH]

Usage:
    python calima_tables.py SEDTABLES_DIR [--png tables.png]
prints every table found and checks it (finite values, an increasing grid, positive cross sections, the
inverse T_dust table against the Planck power), exit code 1 if a check fails; --png also plots them.
From python:
    import calima_tables as ct
    logT, sig = ct.read(d, 'rosseland_mean', 1)              # dust bin 1
    logT, sig = ct.read(d, 'planck_mean', 2, pah='i')        # PAH bin 2, ionised
    ct.cross_section(d, 'rosseland_mean', 1, T=30.)          # interpolated in log T, as RAMSES does
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

DUST_KINDS = ("rosseland_mean", "planck_mean", "planck_power", "planck_deriv", "tdust_tab")
PAH_KINDS = ("rosseland_mean", "planck_mean")


def path(d, kind, ibin, pah=None):
    """File name of a table: kind of DUST_KINDS (dust) or PAH_KINDS (pah='n' or 'i')."""
    if pah is None:
        return os.path.join(d, f"{kind}_DustBin_{ibin:02d}.list")
    return os.path.join(d, f"{kind}_PAHBin_{ibin:02d}_{pah}.list")


def read(d, kind, ibin, pah=None):
    """The two columns of a table (see the module docstring for what they are)."""
    x, y = np.loadtxt(path(d, kind, ibin, pah), unpack=True)
    return x, y


def cross_section(d, kind, ibin, T, pah=None):
    """A mean cross section [cm^2 per grain or PAH] at temperature(s) T [K], linear in log10 T and clamped
    to the table, as RAMSES interpolates it."""
    logT, sig = read(d, kind, ibin, pah)
    return np.interp(np.log10(np.asarray(T, dtype=float)), logT, sig)


def bins(d):
    """The dust and PAH bin numbers that have tables in d."""
    dust = sorted({int(m.group(1)) for f in glob.glob(os.path.join(d, "*_DustBin_*.list"))
                   if (m := re.search(r"_DustBin_(\d+)\.list$", f))})
    pah = sorted({int(m.group(1)) for f in glob.glob(os.path.join(d, "*_PAHBin_*.list"))
                  if (m := re.search(r"_PAHBin_(\d+)_[ni]\.list$", f))})
    return dust, pah


def check(d):
    """Print and check every table in d; returns the list of problems."""
    bad = []
    dust, pah = bins(d)
    if not dust and not pah:
        return [f"no CALIMA tables in {d}"]
    entries = [(k, i, None) for i in dust for k in DUST_KINDS] + [(k, i, c) for i in pah for k in PAH_KINDS
                                                                  for c in "ni"]
    for kind, ibin, c in entries:
        f = path(d, kind, ibin, c)
        name = os.path.basename(f)
        if not os.path.exists(f):
            bad.append(f"{name}: missing")
            continue
        x, y = read(d, kind, ibin, c)
        msg = f"{name:40s} {len(x):5d} rows, x {x.min():9.3f}..{x.max():9.3f}, y {y.min():11.4e}..{y.max():11.4e}"
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
            bad.append(f"{name}: non-finite values")
        if np.any(np.diff(x) <= 0):
            bad.append(f"{name}: the first column is not increasing")
        if kind in ("rosseland_mean", "planck_mean") and np.any(y <= 0):
            bad.append(f"{name}: non-positive cross sections")
        print(msg)
    for ibin in dust:
        # the inverse table must give back the temperature of the Planck power it was resampled from
        try:
            logT, logP = read(d, "planck_power", ibin)
            logPi, logTi = read(d, "tdust_tab", ibin)
        except OSError:
            continue
        err = np.max(np.abs(np.interp(logTi, logT, logP) - logPi))
        print(f"DustBin_{ibin:02d}: T_dust table against the Planck power: max |d log P| = {err:.2e}")
        if err > 1e-2:
            bad.append(f"DustBin_{ibin:02d}: tdust_tab is not the inverse of planck_power ({err:.2e} dex)")
    return bad


def plot(d, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dust, pah = bins(d)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    for i in dust:
        for kind, ls in (("rosseland_mean", "-"), ("planck_mean", "--")):
            x, y = read(d, kind, i)
            ax[0].loglog(10**x, y, ls, label=f"DustBin_{i:02d} {kind.split('_')[0]}")
        x, y = read(d, "planck_power", i)
        ax[2].loglog(10**x, 10**y, label=f"DustBin_{i:02d}")
    for i in pah:
        for c, ls in (("n", "-"), ("i", ":")):
            x, y = read(d, "rosseland_mean", i, c)
            ax[1].loglog(10**x, y, ls, label=f"PAHBin_{i:02d} {c} Rosseland")
    ax[0].set(xlabel="T [K]", ylabel="cross section [cm$^2$/grain]", title="dust: Rosseland (-), Planck (--)")
    ax[1].set(xlabel="T [K]", ylabel="cross section [cm$^2$/PAH]", title="PAHs: neutral (-), ionised (:)")
    ax[2].set(xlabel="T [K]", ylabel="P$_{emit}$ [erg/s/grain]", title="dust emitted power")
    for a in ax:
        a.legend(fontsize=7, frameon=False)
        a.grid(alpha=0.3)
    fig.savefig(out, dpi=110)
    print(f"saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("sed_dir")
    ap.add_argument("--png", default=None)
    a = ap.parse_args()
    bad = check(a.sed_dir)
    if a.png:
        plot(a.sed_dir, a.png)
    for b in bad:
        print("FAIL:", b)
    print("RESULT:", "FAIL" if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
