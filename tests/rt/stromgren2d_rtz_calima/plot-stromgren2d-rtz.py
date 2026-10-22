#!/usr/bin/env python3
"""2D Stromgren test (tests/rt/stromgren2d_rtz_calima: RTZ + CALIMA, 8 groups,
4e4 K blackbody point source in the corner) with several grain-charging models,
one run directory each, named after its charging_model: WDB06isrf (uniform-ISRF
tables at the local G0) and WDB06rt (per-group charge balance).

Per snapshot and run: radial profiles about the source of x_HII, T and n_e.
At one snapshot, on a random subsample of cells, the grain potential of every
dust bin and the total PE heating and recombination cooling per H as each
run's model gives them on its own cell state: WDB06rt = pyCALIMA
rt_group_charging.solve_cell (RAMSES computes the same; compare_ramses_rtgroups.py
checks that), WDB06isrf =
charging_isrf_tables.isrf_lookup at G0 of the groups with 5.6 < egy < 13.6 eV.

The outputs are read directly (read_ramses: leaf cells of the uniform grid);
yt cannot parse the RTZ info_rt file nor 2D runs with boundary regions.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pycalima.models.dust_charge import charging_isrf_tables as cit
from pycalima.models.dust_charge import rt_group_charging as rg


HERE = Path.cwd()          # holds the run directories, one per charging_model
ISRF_TABLES = RTG_TABLES = None   # set from the command line
SOURCE_PC = np.array([0.0, 0.0])
M_U, KB, C = 1.66054e-24, 1.380649e-16, 2.99792458e10
M_HION, M_HEION = 1.67262192e-24, 6.6446573e-24
# element, mass [amu], number of ion stages (hydro_file_descriptor order)
ELEMENTS = (("HYDROGEN", "H", 1.00794, 2), ("HELIUM", "He", 4.002602, 3), ("CARBON", "C", 12.0107, 7),
            ("NITROGEN", "N", 14.0067, 6), ("OXYGEN", "O", 15.9994, 6), ("NEON", "Ne", 20.1797, 6),
            ("MAGNESIUM", "Mg", 24.305, 6), ("SILICON", "Si", 28.0855, 6), ("SULFUR", "S", 32.065, 6),
            ("IRON", "Fe", 55.845, 6))
BINS = ("DustBin_01", "DustBin_02", "DustBin_03", "DustBin_04")
GRAIN_RHO = (2.2, 2.2, 3.3, 3.3)   # sgrain [g/cm^3]
BIN_LABEL = ("graphite 0.01 um", "graphite 0.1 um", "silicate 0.005 um", "silicate 0.1 um")
MODELS = {"WDB06isrf": ("uniform-ISRF tables (WDB06isrf)", "#eb6834", "--"), "WDB06rt": ("per-group WDB06rt", "#2a78d6", "-")}
RUNS = ()                  # (run directory = model, label, color, linestyle), from the command line


def _records(path):
    with open(path, "rb") as f:
        buf = f.read()
    out, i = [], 0
    while i < len(buf):
        n = int(np.frombuffer(buf, "<i4", 1, i)[0])
        out.append(buf[i + 4:i + 4 + n])
        i += n + 8
    return out


def _info(path):
    d = {}
    for ln in path.read_text().splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            try:
                d[k.strip()] = float(v.split()[0])
            except ValueError:
                pass
    return d


def read_ramses(out):
    """Leaf cells of a RAMSES output: positions [box units] and the hydro and
    rt variables by descriptor name (code units)."""
    iout = out.name.split("_")[1]
    info = _info(out / f"info_{iout}.txt")
    ncpu, ndim = int(info["ncpu"]), int(info["ndim"])
    names = {}
    for kind in ("hydro", "rt"):
        names[kind] = [ln.split(",")[1].strip() for ln in (out / f"{kind}_file_descriptor.txt").read_text().splitlines()
                       if not ln.startswith("#")]
    cols = {k: [] for k in names["hydro"] + names["rt"]}
    pos = []
    for icpu in range(1, ncpu + 1):
        a = _records(out / f"amr_{iout}.out{icpu:05d}")
        nx = np.frombuffer(a[2], "<i4")
        nlev, nbound = int(np.frombuffer(a[3], "<i4")[0]), int(np.frombuffer(a[5], "<i4")[0])
        numbl = np.frombuffer(a[21], "<i4").reshape(nlev, ncpu)
        k = 23
        numbb = np.frombuffer(a[k + 2], "<i4").reshape(nlev, nbound) if nbound else None
        k = k + 3 if nbound else k
        k += 1 + 1 + 1 + 3          # headf.., ordering, bound_key, coarse son/flag1/cpu_map
        grids = {}
        for il in range(nlev):
            for ib in range(ncpu + nbound):
                nc = numbl[il, ib] if ib < ncpu else numbb[il, ib - ncpu]
                if nc == 0:
                    continue
                xg = np.column_stack([np.frombuffer(a[k + 3 + d], "<f8") for d in range(ndim)])
                son = np.column_stack([np.frombuffer(a[k + 4 + 3 * ndim + c], "<i4") for c in range(2 ** ndim)])
                if ib == icpu - 1:
                    grids[il] = (xg, son)
                k += 4 + ndim + 2 * ndim + 3 * 2 ** ndim
        for kind in ("hydro", "rt"):
            h = _records(out / f"{kind}_{iout}.out{icpu:05d}")
            nvar = int(np.frombuffer(h[1], "<i4")[0])
            j = 6
            for il in range(nlev):
                for ib in range(ncpu + nbound):
                    nc = int(np.frombuffer(h[j + 1], "<i4")[0])
                    j += 2
                    if nc == 0:
                        continue
                    if ib == icpu - 1:
                        xg, son = grids[il]
                        dx = 0.5 ** (il + 1)
                        for c in range(2 ** ndim):
                            leaf = son[:, c] == 0
                            for v in range(nvar):
                                cols[names[kind][v]].append(np.frombuffer(h[j + c * nvar + v], "<f8")[leaf])
                            if kind == "hydro":
                                off = np.array([((c >> d) & 1) - 0.5 for d in range(ndim)]) * dx
                                pos.append(xg[leaf] + off - nx[:ndim] // 2)
                    j += 2 ** ndim * nvar
    cells = {k: np.concatenate(v) for k, v in cols.items()}
    return np.concatenate(pos) * info["boxlen"], cells, info


def group_energies(out):
    txt = next(out.glob("info_rt_*.txt")).read_text().splitlines()
    return np.array([float(ln.split("=")[1]) for ln in txt if ln.strip().startswith("egy")])


def cell_state(out):
    pos, c, info = read_ramses(out)
    rt = _info(next(out.glob("info_rt_*.txt")))
    r = lambda k: c[k]  # noqa: E731
    ul, ud, ut = info["unit_l"], info["unit_d"], info["unit_t"]
    rho = c["density"] * ud
    s = dict(t=info["time"] * ut / 3.15576e13, r=np.linalg.norm(pos * ul / 3.0856776e18 - SOURCE_PC, axis=1),
             rho=rho)
    ne = np.zeros_like(rho)
    ntot = np.zeros_like(rho)
    for name, sym, A, nst in ELEMENTS:
        n = rho * r(name) / (A * M_U)
        x = np.column_stack([r(f"{sym}_{j + 1:02d}") for j in range(nst)])
        ne += n * (x * np.arange(nst)).sum(axis=1)
        ntot += n
        s["n_" + sym] = n
        s["x_" + sym] = x
    s["ne"] = ne
    s["T"] = c["pressure"] * ud * (ul / ut) ** 2 / (KB * (ntot + ne))
    s["Np"] = np.column_stack([c[f"photon_flux_{g + 1:02d}"] * rt["unit_np"] for g in range(8)])
    s["rho_dust"] = np.column_stack([rho * r(b) for b in BINS])
    s["egy"] = group_energies(out)
    s["c_red"] = 0.001 * C
    return s


def radial_median(r, y, edges):
    idx = np.digitize(r, edges) - 1
    return np.array([np.median(y[idx == i]) if np.any(idx == i) else np.nan for i in range(edges.size - 1)])


def potentials(s, model, tables, isrftabs, cells):
    """Grain potential [V] per bin, and PE heating and recombination cooling per H per bin
    [erg/s], on the given cells, as the run's charging model gives them."""
    phi = np.zeros((len(BINS), cells.size))
    heat = np.full_like(phi, np.nan)
    cool = np.full_like(phi, np.nan)
    gts = [rg.GroupChargingTable(tab, egy=s["egy"], kernels="pchip") for tab in tables]
    G0_bg = 0.0
    for j, c in enumerate(cells):
        T, ne = s["T"][c], s["ne"][c]
        n_gr = [s["rho_dust"][c, b] / (4.0 / 3.0 * np.pi * tables[b]["a"] ** 3 * GRAIN_RHO[b]) for b in range(len(BINS))]
        if model == "WDB06rt":
            ions = [(s["n_H"][c] * s["x_H"][c, 1], 1.0, M_HION), (s["n_He"][c] * s["x_He"][c, 1], 1.0, M_HEION),
                    (s["n_He"][c] * s["x_He"][c, 2], 2.0, M_HEION)]
            res = rg.solve_cell(gts, n_gr, s["Np"][c], s["c_red"], T, ne, ions, G0_bg=G0_bg)
            Z = [r["Zmean"] for r in res]
            G = [r["Gamma"] for r in res]
            L = [r["Lambda"] + r["Lambda_auto"] for r in res]
        else:
            # RTZ passes CALIMA the G0 of these groups plus the (unattenuated) background
            sel = (s["egy"] > 5.6) & (s["egy"] < 13.6)
            G0 = float(np.sum(s["Np"][c][sel] * s["c_red"] * s["egy"][sel]) * 1.602176634e-12 / 1.6e-3) + G0_bg
            res = [cit.isrf_lookup(t, G0, T, ne) for t in isrftabs]
            Z = [r["Zmean"] for r in res]
            G = [r["Gamma"] for r in res]
            L = [r["Lambda"] for r in res]
        for b in range(len(BINS)):
            phi[b, j] = Z[b] * rg.E2_EV_CM / tables[b]["a"]
            heat[b, j] = n_gr[b] * G[b] / s["n_H"][c]
            cool[b, j] = n_gr[b] * L[b] / s["n_H"][c]
    return phi, heat, cool


def main(n_cells=600, phi_snap=3):
    outs = {m: sorted((HERE / m).glob("output_*")) for m, _, _, _ in RUNS}
    n = min(len(v) for v in outs.values())
    if n < 2:
        raise SystemExit("need at least two snapshots in each run")
    snaps = sorted({1, max(2, n // 3), max(2, 2 * n // 3), n, phi_snap})
    tables = [rg.read_bin_table(RTG_TABLES / f"dust_charging_rtgroups_{b}.dat") for b in BINS]
    isrftabs = [cit.read_isrf_table(ISRF_TABLES / f"dust_charging_isrf_{b}.dat") for b in BINS] if ISRF_TABLES else None
    edges = np.linspace(0.0, 1000.0, 41)
    rc = 0.5 * (edges[1:] + edges[:-1])

    fig = plt.figure(figsize=(16, 14))
    gs = fig.add_gridspec(4, 4)
    axx, axT, axe = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[1, 0:2])
    axN = fig.add_subplot(gs[1, 2:4])
    axp = [fig.add_subplot(gs[2, k]) for k in range(4)]
    axh, axc = fig.add_subplot(gs[3, 0:2]), fig.add_subplot(gs[3, 2:4])
    last = {}
    rng = np.random.default_rng(1)
    for m, label, color, ls in RUNS:
        for k, i in enumerate(snaps):
            s = cell_state(outs[m][i - 1])
            a = 0.3 + 0.7 * (k + 1) / len(snaps)
            lab = f"{label}, t = {s['t']:.2f} Myr" if k == len(snaps) - 1 else None
            axx.plot(rc, radial_median(s["r"], s["x_H"][:, 1], edges), color=color, alpha=a, ls=ls, label=lab)
            axT.semilogy(rc, radial_median(s["r"], s["T"], edges), color=color, alpha=a, ls=ls)
            axe.semilogy(rc, radial_median(s["r"], s["ne"], edges), color=color, alpha=a, ls=ls)
            if i == phi_snap:
                last[m] = s
                for g in (2, 3, 4):
                    axN.semilogy(rc, radial_median(s["r"], s["Np"][:, g], edges), color=color, ls=ls,
                                 lw=1 + 0.5 * (g - 2), label=f"{m}: group {g + 1} ({s['egy'][g]:.1f} eV)")
    m0 = RUNS[0][0]
    cells = rng.choice(last[m0]["r"].size, size=min(n_cells, last[m0]["r"].size), replace=False)
    for m, label, color, ls in RUNS:
        s = last[m]
        phi, heat, cool = potentials(s, m, tables, isrftabs, cells)
        e2 = np.linspace(0, s["r"][cells].max(), 25)
        r2 = 0.5 * (e2[1:] + e2[:-1])
        for b in range(len(BINS)):
            axp[b].scatter(s["r"][cells], phi[b], s=4, color=color, alpha=0.3, lw=0)
            axp[b].plot(r2, radial_median(s["r"][cells], phi[b], e2), color=color, lw=2, ls=ls, label=label)
        if np.isfinite(heat).any():
            axh.semilogy(r2, radial_median(s["r"][cells], heat.sum(axis=0), e2), color=color, lw=2, ls=ls, label=label)
            axc.semilogy(r2, radial_median(s["r"][cells], cool.sum(axis=0), e2), color=color, lw=2, ls=ls, label=label)
        np.savez(HERE / f"stromgren_phi_{m}.npz", r=s["r"][cells], phi=phi, T=s["T"][cells], ne=s["ne"][cells],
                 Np=s["Np"][cells], xHII=s["x_H"][cells, 1], egy=s["egy"], t=s["t"], heat=heat, cool=cool)
    axx.set_ylabel(r"$x_{\rm HII}$")
    axT.set_ylabel("T (K)")
    axe.set_ylabel(r"$n_e$ (cm$^{-3}$)")
    axN.set_ylabel(r"$N_g$ (cm$^{-3}$)")
    axh.set_ylabel(r"PE heating per H, all dust bins (erg s$^{-1}$)")
    axc.set_ylabel(r"recombination cooling per H (erg s$^{-1}$)")
    for ax in (axx, axT, axe, axN, axh, axc):
        ax.set_xlabel("r from source (pc)")
        ax.grid(color="0.92", lw=0.6)
    axx.legend(frameon=False, fontsize=7)
    axN.legend(frameon=False, fontsize=6)
    axh.legend(frameon=False, fontsize=8)
    for b, ax in enumerate(axp):
        ax.set_title(BIN_LABEL[b], fontsize=10)
        ax.set_xlabel("r from source (pc)")
        ax.axhline(0.0, color="0.6", lw=0.8)
        ax.grid(color="0.92", lw=0.6)
    axp[0].set_ylabel(r"grain potential $\phi$ (V), t = %.2f Myr" % last[m0]["t"])
    axp[0].legend(frameon=False, fontsize=8)
    fig.suptitle("2D Stromgren, RTZ + CALIMA, by grain-charging model (line style per model)", fontsize=12)
    fig.tight_layout()
    out = HERE / "stromgren_charging_models.png"
    fig.savefig(out, dpi=130)
    print("Saved", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rtgroups-tables", required=True, help="dir with dust_charging_rtgroups_DustBin_XX.dat")
    ap.add_argument("--isrf-tables", default=None, help="dir with dust_charging_isrf_DustBin_XX.dat (WDB06isrf)")
    ap.add_argument("--runs", nargs="+", default=["WDB06isrf", "WDB06rt"],
                    help="run directories, named after their charging_model")
    ap.add_argument("--phi-snapshot", type=int, default=3, help="snapshot for the grain potentials (3: 2 Myr)")
    args = ap.parse_args()
    RTG_TABLES = Path(args.rtgroups_tables)
    ISRF_TABLES = Path(args.isrf_tables) if args.isrf_tables else None
    RUNS = tuple((m,) + MODELS[m] for m in args.runs)
    main(phi_snap=args.phi_snapshot)
