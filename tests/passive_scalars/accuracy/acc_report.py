"""Tables and figures from an advection-accuracy run (acc_run.py).

    acc_report.py <OUT> [--label TEXT]

Reads OUT/results.json (and OUT/runs/*/profile.npz), writes OUT/summary.md and OUT/fig_*.png.
"""
import argparse
import json
import math
import os
import sys

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import acc_cases as A  # noqa: E402

SLOPE_NAME = {0: "0 first order", 1: "1 minmod", 2: "2 MC", 3: "3 multi-D MC", 4: "4 superbee",
              5: "5 ultrabee (= 0: density only)", 6: "6 unlimited (= 0: density only)", 7: "7 van Leer", 8: "8 gen. MC (theta 1.5)"}
# colour follows the slope type in every figure
SLOPE_COLOR = {0: "#7f7f7f", 1: "#1f77b4", 2: "#d62728", 3: "#ff7f0e", 4: "#2ca02c", 5: "#9467bd",
               6: "#8c564b", 7: "#17becf", 8: "#e377c2"}
PROF_TITLE = {"tophat": "top-hat", "gauss": "Gaussian", "sine": "sine"}


def load(out):
    R = json.load(open(os.path.join(out, "results.json")))
    return [r for r in R["results"] if r["status"] == "OK"], R


def pick(res, **kw):
    out = []
    for r in res:
        if all((tuple(r[k]) if k == "k" else r[k]) == v for k, v in kw.items()):
            out.append(r)
    return sorted(out, key=lambda r: r["N"])


def order(rs, prof, metric="L1"):
    """Least-squares convergence order over the two finest resolutions."""
    if len(rs) < 2:
        return float("nan")
    a, b = rs[-2], rs[-1]
    return math.log(a[prof][metric] / b[prof][metric]) / math.log(b["N"] / a["N"])


def base(**kw):
    d = dict(courant=0.8, periods=1)
    d.update(kw)
    return d


def fig_convergence(res, out, solver="hydro", ndim=1, k=(1, 0), riemann="hllc", mach=1.0, slopes=None,
                    fname="fig_convergence_hydro1d.png", title=None):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    slopes = slopes or sorted({r["slope"] for r in pick(res, **base(solver=solver, ndim=ndim, k=k))})
    for ip, p in enumerate(A.PROFILES):
        for s in slopes:
            rs = pick(res, **base(solver=solver, ndim=ndim, k=k, riemann=riemann, mach=mach, slope=s))
            if not rs:
                continue
            N = [r["N"] for r in rs]
            ax[ip].loglog(N, [r[p]["L1"] for r in rs], "o-", color=SLOPE_COLOR[s], lw=2, ms=6,
                          label=f"{SLOPE_NAME[s]}  (p={order(rs, p):.2f})")
        Nr = np.array([32, 512])
        ref = {"tophat": ((0.5, "N$^{-1/2}$"), (2 / 3, "N$^{-2/3}$")), "gauss": ((1, "N$^{-1}$"), (2, "N$^{-2}$")),
               "sine": ((1, "N$^{-1}$"), (2, "N$^{-2}$"))}[p]
        y0 = ax[ip].get_ylim()[1] * 0.5
        for (e, lab), ls in zip(ref, (":", "--")):
            ax[ip].loglog(Nr, y0 * (Nr / Nr[0]) ** (-e), ls, color="0.4", lw=1, label=lab)
        ax[ip].set(xlabel="N cells", ylabel="L1 error", title=PROF_TITLE[p])
        nticks(ax[ip], A.NS_1D)
        ax[ip].legend(frameon=False, fontsize=7.5)
        ax[ip].grid(alpha=0.25, which="both")
    fig.suptitle(title or f"{solver} {ndim}-D, {riemann}, Mach {mach:g}: passive scalar L1 error after one period "
                 "(p: order between the two finest N)")
    fig.savefig(os.path.join(out, fname), dpi=110)
    plt.close(fig)


def fig_profiles(res, out, N=64):
    sl = [0, 1, 2, 4, 5, 6, 7, 8]
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.4), constrained_layout=True)
    for ip, p in enumerate(A.PROFILES):
        first = True
        for s in sl:
            rs = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann="hllc", mach=1.0, slope=s, N=N))
            if not rs:
                continue
            z = np.load(os.path.join(out, "runs", rs[0]["case"], "profile.npz"))
            if first:
                ax[0, ip].step(z["s"], z[f"qe_{p}"], where="mid", color="k", lw=2.5, label="exact cell averages")
                first = False
            ax[0, ip].plot(z["s"], z[f"q_{p}"], "-", color=SLOPE_COLOR[s], lw=1.5, label=SLOPE_NAME[s])
        ax[0, ip].set(title=f"{PROF_TITLE[p]}: hllc, Mach 1, N={N}", xlabel="x", ylabel="q")
        first = True
        for r_, ls in (("hllc", "-"), ("hll", "--"), ("llf", ":")):
            rs = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann=r_, mach=0.1, slope=2, N=N))
            if not rs:
                continue
            z = np.load(os.path.join(out, "runs", rs[0]["case"], "profile.npz"))
            if first:
                ax[1, ip].step(z["s"], z[f"qe_{p}"], where="mid", color="k", lw=2.5, label="exact cell averages")
                first = False
            ax[1, ip].plot(z["s"], z[f"q_{p}"], ls, color=SLOPE_COLOR[2], lw=1.8,
                           label=f"MC, {r_}: min {rs[0][p]['qmin']:+.3f}")
        ax[1, ip].set(title=f"{PROF_TITLE[p]}: MC, Mach 0.1, N={N}, by Riemann solver", xlabel="x", ylabel="q")
        for a in ax[:, ip]:
            a.grid(alpha=0.25)
        ax[0, ip].legend(frameon=False, fontsize=7)
        ax[1, ip].legend(frameon=False, fontsize=8)
    if "tophat" in A.PROFILES:
        for a in ax[:, 0]:
            a.set_xlim(0.2, 0.8)
    fig.suptitle("Final profiles after one period (top: slope limiters; bottom: Riemann solvers at low Mach)")
    fig.savefig(os.path.join(out, "fig_profiles.png"), dpi=110)
    plt.close(fig)


def nticks(ax, Ns):
    ax.set_xticks(Ns, [str(n) for n in Ns])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def heat(ax, M, rows, cols, fmt, title, cmap="viridis_r", log=True):
    data = np.array(M, float)
    show = np.log10(np.clip(data, 1e-30, None)) if log else data
    im = ax.imshow(np.ma.masked_invalid(show), cmap=cmap, aspect="auto")
    fin = show[np.isfinite(show)]
    lo, span = (fin.min(), max(np.ptp(fin), 1e-30)) if fin.size else (0.0, 1.0)
    for i in range(len(rows)):
        for j in range(len(cols)):
            if np.isfinite(data[i, j]):
                ax.text(j, i, fmt(data[i, j]), ha="center", va="center", fontsize=7.5,
                        color="w" if (show[i, j] - lo) > 0.5 * span else "k")   # both maps: high = dark
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)), rows, fontsize=8)
    ax.set_title(title, fontsize=10)
    return im


def fig_riemann(res, out, N=128):
    slopes = list(A.HYDRO_SLOPES_1D)
    cols = [(r, m) for m in (0.1, 1.0, 10.0) for r in A.HYDRO_RIEMANN]
    lab = [f"{r} M{m:g}" for r, m in cols]
    mats = {k: [] for k in ("sine", "gauss", "ext")}
    for s in slopes:
        row = {k: [] for k in mats}
        for r_, m in cols:
            rs = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann=r_, mach=m, slope=s, N=N))
            if not rs:
                for k in row:
                    row[k].append(np.nan)
                continue
            x = rs[0]
            row["sine"].append(x["sine"]["L1"])
            row["gauss"].append(x["gauss"]["L1"])
            row["ext"].append(max(max(x[p]["over"], x[p]["under"]) for p in A.PROFILES))
        for k in mats:
            mats[k].append(row[k])
    fig, ax = plt.subplots(3, 1, figsize=(13, 13), constrained_layout=True)
    rl = [SLOPE_NAME[s] for s in slopes]
    heat(ax[0], mats["sine"], rl, lab, lambda v: f"{v:.1e}", f"sine L1 error, hydro 1-D, N={N}")
    heat(ax[1], mats["gauss"], rl, lab, lambda v: f"{v:.1e}", f"Gaussian L1 error, hydro 1-D, N={N}")
    e = np.array(mats["ext"])
    heat(ax[2], np.where(e > 1e-12, e, np.nan), rl, lab, lambda v: f"{v:.0e}",
         f"new extrema max(overshoot, undershoot) over the 3 profiles, N={N} (blank: none above 1e-12)",
         cmap="magma_r")
    fig.suptitle("Slope limiter x Riemann solver x Mach number (Courant factor 0.8)")
    fig.savefig(os.path.join(out, "fig_riemann_mach.png"), dpi=110)
    plt.close(fig)


def fig_extrema_courant(res, out):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for j, r_ in enumerate(("hll", "llf")):
        for s in (1, 2, 4, 7, 8):
            for N, ls in ((64, "--"), (256, "-")):
                pts = []
                for cf in (0.2, 0.4, 0.8):
                    rs = pick(res, solver="hydro", ndim=1, k=(1, 0), riemann=r_, mach=0.1, slope=s, N=N,
                              courant=cf, periods=1)
                    if rs:
                        pts.append((cf, max(max(rs[0][p]["over"], rs[0][p]["under"]) for p in A.PROFILES)))
                if pts:
                    x, y = zip(*pts)
                    ax[j].semilogy(x, np.maximum(y, 1e-17), "o" + ls, color=SLOPE_COLOR[s], lw=1.8,
                                   label=f"{SLOPE_NAME[s]}, N={N}")
        ax[j].set(xlabel="courant_factor", ylabel="max new extremum", ylim=(1e-17, 1),
                  title=f"{r_}, Mach 0.1: over/undershoots vs Courant factor (1e-17: none)")
        ax[j].legend(frameon=False, fontsize=7, ncol=2)
        ax[j].grid(alpha=0.25, which="both")
    fig.savefig(os.path.join(out, "fig_extrema_courant.png"), dpi=110)
    plt.close(fig)


def fig_2d(res, out):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for ip, p in enumerate(A.PROFILES):
        for s in A.HYDRO_SLOPES_2D:
            for k, ls, tag in (((1, 1), "-", "diagonal"), ((1, 0), ":", "aligned")):
                rs = pick(res, **base(solver="hydro", ndim=2, k=k, riemann="hllc", mach=1.0, slope=s))
                if rs:
                    ax[ip].loglog([r["N"] for r in rs], [r[p]["L1"] for r in rs], "o" + ls, color=SLOPE_COLOR[s],
                                  lw=1.8, ms=5, label=f"{SLOPE_NAME[s]} {tag}")
        for s in A.MHD_SLOPES_2D:
            rs = pick(res, **base(solver="mhd", ndim=2, k=(1, 1), riemann="hlld", mach=1.0, slope=s))
            if rs:
                ax[ip].loglog([r["N"] for r in rs], [r[p]["L1"] for r in rs], "s--", color=SLOPE_COLOR[s], lw=1,
                              ms=4, mfc="none", label=f"MHD hlld {s} diagonal")
        ax[ip].set(xlabel="N cells per side", ylabel="L1 error", title=PROF_TITLE[p])
        nticks(ax[ip], A.NS_2D)
        ax[ip].grid(alpha=0.25, which="both")
    ax[0].legend(frameon=False, fontsize=6.5, ncol=2)
    fig.suptitle("2-D: advection along the diagonal (solid) and along x (dotted), hllc, Mach 1; MHD hlld (squares)")
    fig.savefig(os.path.join(out, "fig_2d.png"), dpi=110)
    plt.close(fig)


def fig_mhd(res, out, N=128):
    slopes = list(A.MHD_SLOPES_1D)
    cols = [(r, m) for m in (0.1, 1.0) for r in A.MHD_RIEMANN]
    lab = [f"{r} M{m:g}" for r, m in cols]
    M1, M2 = [], []
    for s in slopes:
        r1, r2 = [], []
        for r_, m in cols:
            rs = pick(res, **base(solver="mhd", ndim=1, k=(1, 0), riemann=r_, mach=m, slope=s, N=N))
            r1.append(rs[0]["sine"]["L1"] if rs else np.nan)
            r2.append(max(max(rs[0][p]["over"], rs[0][p]["under"]) for p in A.PROFILES) if rs else np.nan)
        M1.append(r1)
        M2.append(r2)
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), constrained_layout=True)
    rl = [SLOPE_NAME[s] for s in slopes]
    heat(ax[0], M1, rl, lab, lambda v: f"{v:.1e}", f"MHD 1-D sine L1 error, N={N} (beta=10 transverse field)")
    e = np.array(M2)
    heat(ax[1], np.where(e > 1e-12, e, np.nan), rl, lab, lambda v: f"{v:.0e}", "new extrema (blank: none)",
         cmap="magma_r")
    fig.savefig(os.path.join(out, "fig_mhd.png"), dpi=110)
    plt.close(fig)


def fmt(v, f=".2e"):
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, f)


def summary(res, R, out, label):
    L = [f"# Passive-scalar advection accuracy {label}", "",
         f"{R['ncases']} cases ({len(res)} completed), wall {R['wall_s']:.0f} s.  Uniform rho=1, c_s=1, periodic "
         "unit box, one period unless noted; errors against exact cell averages.  nu = v dt/dx.", ""]
    mass = max(max(r[p]["mass"] for p in A.PROFILES) for r in res)
    L += [f"**Mass conservation:** worst |dM/M| over all runs and scalars = {mass:.1e}.", ""]
    th = [(r["case"], r["sine"]["amp"], r["sine"]["amp_theory"]) for r in res if "amp_theory" in r["sine"]]
    if th:
        worst = max(th, key=lambda x: abs(x[1] - x[2]))
        L += [f"**First-order check (slope_type 0):** the sine amplitude matches |G|^n of the first-order scheme with "
              f"each solver's scalar flux in {len(th)} runs; worst |amp - theory| = {abs(worst[1] - worst[2]):.1e} "
              f"({worst[0]}).", ""]
    # hydro 1-D table
    L += ["## Hydro 1-D, hllc, Mach 1, N = 128", "",
          "| slope | nu | L1 top-hat | L1 Gauss | L1 sine | order top/Gauss/sine | top-hat width [cells] | "
          "D_num [v dx] | sine amp | sine shift [cells] | new extrema | TV top-hat |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for s in A.HYDRO_SLOPES_1D:
        rs = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann="hllc", mach=1.0, slope=s))
        x = [r for r in rs if r["N"] == 128]
        if not x:
            continue
        x = x[0]
        ords = "/".join(f"{order(rs, p):.2f}" for p in A.PROFILES)
        ext = max(max(x[p]["over"], x[p]["under"]) for p in A.PROFILES)
        L.append(f"| {SLOPE_NAME[s]} | {x['nu']:.2f} | {fmt(x['tophat']['L1'])} | {fmt(x['gauss']['L1'])} | "
                 f"{fmt(x['sine']['L1'])} | {ords} | {fmt(x['tophat'].get('width'), '.1f')} | "
                 f"{fmt(x['gauss']['D_num'], '.4f')} | {x['sine']['amp']:.5f} | {x['sine']['shift']:+.3f} | "
                 f"{fmt(ext, '.0e') if ext > 1e-12 else 'none'} | {x['tophat']['tv_ratio']:.3f} |")
    # Riemann x Mach
    L += ["", "## Riemann solver and Mach number (hydro 1-D, N = 128): Gaussian D_num [v dx] / new extrema", "",
          "| slope | " + " | ".join(f"{r_} M{m:g}" for m in (0.1, 1.0, 10.0) for r_ in A.HYDRO_RIEMANN) + " |",
          "|---|" + "---|" * (3 * len(A.HYDRO_RIEMANN))]
    for s in A.HYDRO_SLOPES_1D:
        cells = []
        for m in (0.1, 1.0, 10.0):
            for r_ in A.HYDRO_RIEMANN:
                x = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann=r_, mach=m, slope=s, N=128))
                if not x:
                    cells.append("—")
                    continue
                x = x[0]
                ext = max(max(x[p]["over"], x[p]["under"]) for p in A.PROFILES)
                d = x["gauss"]["D_num"] if x["gauss"].get("sigma", 0) < 0.12 else float("nan")
                cells.append(f"{fmt(d, '.3f')}" + (f" / **{ext:.0e}**" if ext > 1e-12 else ""))
        L.append(f"| {SLOPE_NAME[s]} | " + " | ".join(cells) + " |")
    L += ["", "D_num is left out (—) where the Gaussian has spread over the box (sigma > 0.12).", ""]
    # extrema vs Courant
    L += ["## New extrema vs Courant factor (hydro 1-D, Mach 0.1, N = 256): max over/undershoot", "",
          "| slope | hll C=0.2 | hll 0.4 | hll 0.8 | llf 0.2 | llf 0.4 | llf 0.8 |", "|---|---|---|---|---|---|---|"]
    for s in (1, 2, 4, 7, 8):
        cells = []
        for r_ in ("hll", "llf"):
            for cf in (0.2, 0.4, 0.8):
                x = pick(res, solver="hydro", ndim=1, k=(1, 0), riemann=r_, mach=0.1, slope=s, N=256, courant=cf,
                         periods=1)
                ext = max(max(x[0][p]["over"], x[0][p]["under"]) for p in A.PROFILES) if x else None
                cells.append("—" if ext is None else ("none" if ext <= 1e-12 else f"{ext:.0e}"))
        L.append(f"| {SLOPE_NAME[s]} | " + " | ".join(cells) + " |")
    # Courant and periods with hllc
    L += ["", "## Courant factor and 10 periods (hydro 1-D, hllc, Mach 1, N = 256): L1 sine / L1 top-hat", "",
          "| slope | C=0.8 | C=0.4 | C=0.8, 10 periods |", "|---|---|---|---|"]
    for s in (1, 2, 7, 8):
        cells = []
        for kw in (dict(courant=0.8, periods=1), dict(courant=0.4, periods=1), dict(courant=0.8, periods=10)):
            x = pick(res, solver="hydro", ndim=1, k=(1, 0), riemann="hllc", mach=1.0, slope=s, N=256, **kw)
            cells.append(f"{fmt(x[0]['sine']['L1'])} / {fmt(x[0]['tophat']['L1'])}" if x else "—")
        L.append(f"| {SLOPE_NAME[s]} | " + " | ".join(cells) + " |")
    # 2-D
    L += ["", "## 2-D (hllc, Mach 1, N = 128): diagonal vs grid-aligned, L1 Gauss / L1 sine / new extrema", "",
          "| slope | diagonal hllc | diagonal hll | aligned hllc | diagonal/aligned (Gauss) |", "|---|---|---|---|---|"]
    for s in A.HYDRO_SLOPES_2D:
        cells, g = [], []
        for r_, k in (("hllc", (1, 1)), ("hll", (1, 1)), ("hllc", (1, 0))):
            x = pick(res, **base(solver="hydro", ndim=2, k=k, riemann=r_, mach=1.0, slope=s, N=128))
            if not x:
                cells.append("—")
                g.append(np.nan)
                continue
            x = x[0]
            ext = max(max(x[p]["over"], x[p]["under"]) for p in A.PROFILES)
            cells.append(f"{fmt(x['gauss']['L1'])} / {fmt(x['sine']['L1'])} / "
                         + ("none" if ext <= 1e-12 else f"{ext:.0e}"))
            g.append(x["gauss"]["L1"])
        L.append(f"| {SLOPE_NAME[s]} | " + " | ".join(cells) + f" | {g[0] / g[2]:.2f} |")
    # MHD
    L += ["", "## MHD 1-D (beta = 10 transverse field, N = 128): L1 sine / new extrema", "",
          "| slope | " + " | ".join(f"{r_} M{m:g}" for m in (0.1, 1.0) for r_ in A.MHD_RIEMANN) + " |",
          "|---|" + "---|" * (2 * len(A.MHD_RIEMANN))]
    for s in A.MHD_SLOPES_1D:
        cells = []
        for m in (0.1, 1.0):
            for r_ in A.MHD_RIEMANN:
                x = pick(res, **base(solver="mhd", ndim=1, k=(1, 0), riemann=r_, mach=m, slope=s, N=128))
                if not x:
                    cells.append("—")
                    continue
                x = x[0]
                ext = max(max(x[p]["over"], x[p]["under"]) for p in A.PROFILES)
                cells.append(f"{fmt(x['sine']['L1'])}" + (f" / **{ext:.0e}**" if ext > 1e-12 else ""))
        L.append(f"| {SLOPE_NAME[s]} | " + " | ".join(cells) + " |")
    hy = pick(res, **base(solver="hydro", ndim=1, k=(1, 0), riemann="hllc", mach=1.0, slope=2, N=128))
    mh = pick(res, **base(solver="mhd", ndim=1, k=(1, 0), riemann="hlld", mach=1.0, slope=2, N=128))
    if hy and mh:
        L += ["", f"MC: hydro hllc sine L1 {hy[0]['sine']['L1']:.3e} (nu {hy[0]['nu']:.3f}) vs MHD hlld "
              f"{mh[0]['sine']['L1']:.3e} (nu {mh[0]['nu']:.3f})."]
    L += ["", "Figures: fig_convergence_hydro1d.png, fig_convergence_mhd1d.png, fig_profiles.png, "
          "fig_riemann_mach.png, fig_extrema_courant.png, fig_2d.png, fig_mhd.png", ""]
    open(os.path.join(out, "summary.md"), "w").write("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    res, R = load(a.out)
    fig_convergence(res, a.out)
    fig_convergence(res, a.out, solver="mhd", riemann="hlld", fname="fig_convergence_mhd1d.png",
                    title="MHD 1-D, hlld, Mach 1: passive scalar L1 error after one period")
    fig_profiles(res, a.out)
    fig_riemann(res, a.out)
    fig_extrema_courant(res, a.out)
    fig_2d(res, a.out)
    fig_mhd(res, a.out)
    summary(res, R, a.out, a.label)
    print(open(os.path.join(a.out, "summary.md")).read())


if __name__ == "__main__":
    main()
