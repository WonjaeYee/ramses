#!/usr/bin/env python3
"""RTZ's per-ion metal-line cooling tables, from Cloudy, on a (T, n_e) grid:
<data_dir>/high_T_cooling/<ELEMENT>/all_cool_ne.dat, read by rtz_coolrates_module.

Each entry is log10 Lambda_ion(T, n_e) = log10[cooling / (n_e n_ion)] [erg cm^3 s^-1] of one ion stage alone (its lines
and its free-free; not Cloudy's total, see run_one), from a
one-zone Cloudy model with only that ion present (the setup of the earlier high_T_cooling/make_crates.py, which made
the n_e = 1 cm^-3 tables all_cool.dat): fixed T, electrons at n_e, no photo- or collisional ionization, no charge
transfer, hydrogen at 1e-15 cm^-3.

Why n_e: the earlier tables are the low-density limit (n_e = 1). Lines whose critical density is near the gas density
are then overestimated, because their upper levels are collisionally de-excited instead of radiating: at n_e = 300
cm^-3 and 7000 K, Cloudy's N III cools 0.62x the low-density value (57 um, n_crit ~ 1e3 cm^-3), S III 0.91x, O II,
S II, Fe II/III ~0.95x. RTZ's nitrogen cooling was 1.6x Cloudy's on the matched 1-D H II slab at Z_sun.

Why n_e only: Lambda_ion is the ion's level populations, set by collisions with electrons (and radiative decay). Varying
everything else at T = 7000 K, n_e = 300 cm^-3 changes it by < 0.1%: the ion's own abundance (so metallicity and
abundance ratios, which enter through n_ion, as RTZ multiplies), n(H+) = 300 (proton collisions), n(H0) = 300, ISM
grains; a 42900 K blackbody at the ladder's flux (pumping) by <= 2%. Neutral-H collisions matter only in mostly neutral
gas, for ground-term fine-structure lines, which RTZ's fine-structure solver already treats with H0, H2 and He.
(RAMSES_dev/studies/rt_groups_resolution/ladder_final/diag/crates_density.)

What is computed: only the ion stages whose cooling depends on n_e. A screening pass runs every stage at the top of
the n_e grid (1e6 cm^-3) over all T; a stage whose log Lambda moves by more than 0.004 dex (1%) from n_e = 1, at any T
where it cools within 1e-3 of its peak, is computed at every n_e of the grid. The others (most of the high stages,
whose lines have critical densities far above 1e6 cm^-3) keep the n_e = 1 values at every n_e, and the n_e = 1 block
is the earlier all_cool.dat, which these models reproduce exactly. The flagged stages are listed in the file header.

File format (all_cool_ne.dat): '#' comment lines; then "n_ne lg_ne_min dlg_ne"; then n_ne blocks, one per n_e
(increasing), of Z+1 rows (ion stages I..Z+1) of the 121 log10 Lambda on log T = 3, 3.05, ..., 9 (temperatures.dat).

    python make_metal_cooling_tables.py <cloudy.exe> <data_dir>/high_T_cooling [-j 10] [--lgne 0 6 0.5]
Needs Cloudy only to regenerate the tables, not to run RAMSES (~1 h on 10 cores). Each model
grid runs on one core (Cloudy's 'grid ... sequential'; its default parallel grid oversubscribes the machine when the
jobs run side by side); finished jobs are kept in <out>/.cache, so an interrupted run resumes.
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

ELEMENTS = [("carbon", 6), ("nitrogen", 7), ("oxygen", 8), ("neon", 10), ("magnesium", 12), ("silicon", 14),
            ("sulphur", 16), ("iron", 26)]
ALL = ["HELIUM", "LITHIUM", "MAGNESIUM", "BERYLLIUM", "CARBON", "NITROGEN", "SULPHUR", "OXYGEN", "NEON", "SILICON",
       "IRON", "BORON", "FLUORINE", "SODIUM", "ALUMINIUM", "PHOSPHORUS", "CHLORINE", "ARGON", "POTASSIUM", "CALCIUM",
       "SCANDIUM", "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "COBALT", "NICKEL", "COPPER", "ZINC"]
NT = 121                                   # log T = 3 .. 9, step 0.05 (temperatures.dat)
SYM = {"carbon": "C", "nitrogen": "N", "oxygen": "O", "neon": "Ne", "magnesium": "Mg", "silicon": "Si", "sulphur": "S",
       "iron": "Fe"}


def deck(el, an, istate, lgne, lgtmax=9.0):
    z = ["0.0"] * (an + 1)
    z[istate - 1] = "1.0"
    lines = ["constant temperature 4.2 vary", f"grid 3 {lgtmax:g} 0.05 sequential", "stop zone 1", "set dr 0", "hden -15.0",
             f"set eden {lgne} log", 'save cooling each "cooling.dat" no hash', "element hydrogen ionization 1.0 0.0",
             "no ionization reevaluation", "no collisional ionization", "no photoionization", "no charge transfer",
             "no auger effect", "no radiation pressure", "no molecules", f"element {el} abundance 15.0",
             f"element {el} ionization {' '.join(z)}"]
    lines += [f"element {e} off" for e in ALL if e.lower() != el]
    return "\n".join(lines) + "\n"


def run_one(exe, el, an, istate, lgne, lgtmax=9.0, cache=None):
    """log10 Lambda (per n_e n_ion) on log T = 3 .. lgtmax (step 0.05) for one ion stage at one n_e"""
    key = os.path.join(cache, f"{el}_{istate}_{lgne:g}_{lgtmax:g}.npy") if cache else None
    if key and os.path.exists(key):
        return np.load(key)
    nt = int(round((lgtmax - 3.0) / 0.05)) + 1
    d = tempfile.mkdtemp(prefix=f"mct_{el}_{istate}_")
    try:
        open(os.path.join(d, "m.in"), "w").write(deck(el, an, istate, lgne, lgtmax))
        subprocess.run(["nice", exe, "-r", "m"], cwd=d, capture_output=True, check=False)
        dat = np.loadtxt(os.path.join(d, "cooling.dat"))
        if dat.shape[0] != nt:
            raise RuntimeError(f"{el} {istate} lg n_e {lgne}: {dat.shape[0]} temperatures, see {d}")
        # the ion's own cooling [erg cm^-3 s^-1], n_ion = 10^(-15 + 15) = 1 cm^-3: its element's column (its lines) plus
        # the metals' free-free (this ion's), not the total, whose electron-electron bremsstrahlung (eeff, ~ n_e^2) is
        # nil at n_e = 1 but dominates above 1e6 K at n_e = 1e6 (at n_e = 1 the two agree to < 1e-4)
        hdr = [h.strip() for h in open(os.path.join(d, "cooling.dat")).readline().lstrip("#").rstrip("\n").split("\t")]
        own = dat[:, hdr.index(SYM[el])] + dat[:, hdr.index("FF_M")]
        out = np.log10(np.maximum(own, 1e-300) / 10.0 ** lgne)
        if key:
            np.save(key, out)
        return out
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cloudy")
    ap.add_argument("out", help="the high_T_cooling directory of the RTZ data (one subdirectory per element, with the "
                                "n_e = 1 tables all_cool.dat)")
    ap.add_argument("-j", type=int, default=10, help="Cloudy models run at a time (one core each)")
    ap.add_argument("--lgne", type=float, nargs=3, default=(0.0, 6.0, 0.5), metavar=("MIN", "MAX", "STEP"))
    a = ap.parse_args()
    exe = os.path.abspath(a.cloudy)
    lo, hi, st = a.lgne
    grid = np.round(np.arange(lo, hi + 0.5 * st, st), 6)
    cache = os.path.join(a.out, ".cache")
    os.makedirs(cache, exist_ok=True)
    ions = [(el, an, i) for el, an in ELEMENTS for i in range(1, an + 2)]
    base = {el: np.loadtxt(os.path.join(a.out, el.upper(), "all_cool.dat")) for el, _ in ELEMENTS}   # n_e = 1
    print(f"screening {len(ions)} ion stages at lg n_e = {grid[-1]:g}", flush=True)
    with ThreadPoolExecutor(a.j) as ex:
        top = dict(zip(ions, ex.map(lambda q: run_one(exe, *q, grid[-1], 9.0, cache), ions)))
    flagged = []
    for el, an, i in ions:
        b = base[el][i - 1]
        live = b > b.max() - 3
        dev = np.max(np.abs(top[(el, an, i)][live] - b[live]))
        if dev > 0.004:
            flagged.append((el, an, i))
        print(f"  {el} {i}: max |d log Lambda| = {dev:.4f} dex{'  -> n_e grid' if dev > 0.004 else ''}", flush=True)
    jobs = [(el, an, i, g, 9.0, cache) for el, an, i in flagged for g in grid if 0.0 < g < grid[-1]]
    print(f"{len(flagged)} stages depend on n_e: {len(jobs)} more one-core Cloudy grids", flush=True)
    with ThreadPoolExecutor(a.j) as ex:
        res = dict(zip([j[:4] for j in jobs], ex.map(lambda q: run_one(exe, *q), jobs)))
    for el, an, i in flagged:
        res[(el, an, i, grid[-1])] = top[(el, an, i)]
    for el, an in ELEMENTS:
        d = os.path.join(a.out, el.upper())
        fl = [i for e, _, i in flagged if e == el]
        with open(os.path.join(d, "all_cool_ne.dat"), "w") as f:
            f.write(f"# log10 Lambda [erg cm^3 s^-1] per ion and electron of {el}, ion stages I..{an + 1} (rows) on log T = "
                    "3..9 step 0.05 (columns),\n# one block per log10 n_e [cm^-3]. Made by "
                    "utils/py/rtz/make_metal_cooling_tables.py from one-zone Cloudy models; read by rtz_coolrates_module.\n"
                    f"# n_e dependence computed for stages {fl} (> 1% at n_e = 1e6); the others: the n_e = 1 values.\n")
            f.write(f"{len(grid)} {grid[0]:.2f} {st:.2f}\n")
            for g in grid:
                for i in range(1, an + 2):
                    row = base[el][i - 1].copy()
                    if (el, an, i, g) in res:
                        row = res[(el, an, i, g)]
                    f.write(" ".join(f"{v:.6f}" for v in row) + "\n")
        print(f"wrote {d}/all_cool_ne.dat ({an + 1} ions x {len(grid)} n_e x {NT} T)")


if __name__ == "__main__":
    main()
