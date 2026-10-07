#!/usr/bin/env python3
"""
Projections and radial profiles for the 2D Stromgren test, comparing plain
RAMSES-RT against the RT+CALIMA build with dust.

  plain RT : ../stromgren2d_plainrt   (RT=1 RTZ=0 CALIMA=0)
  RT+CALIMA: .                        (RT=1 RTZ=0 CALIMA=1, NDUST=1)

Both runs use the same namelist physics and the same snapshot times, so any
difference is the dust. Run the two simulations first; this script only reads.

Writes: stromgren2d_maps_rt.png, stromgren2d_maps_calima.png,
        stromgren2d_profiles.png
"""
import os, re, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

RT_DIR     = os.path.abspath("../stromgren2d_plainrt")
CAL_DIR    = os.path.abspath(".")
RT_NML     = "stromgren2d_rt.nml"
CAL_NML    = "stromgren2d_dust.nml"
SNAPS      = [2, 3, 4, 5]           # output_00001 is t=0 and uniform

# physical constants, matching amr/constants.f90 so T agrees with the code
mH, kB = 1.6738233e-24, 1.3806490e-16

def nml_val(path, key, default):
    txt = io = open(path).read()
    m = re.search(r"^\s*%s\s*=\s*([0-9eEdD.+-]+)" % key, txt, re.M | re.I)
    return float(m.group(1).lower().replace("d", "e")) if m else default

def load(d, n):
    cwd = os.getcwd(); os.chdir(d)
    try:   return visu_ramses.load_snapshot(n)["data"]
    finally: os.chdir(cwd)

def derived(data, X, has_dust):
    """Physical fields from a snapshot. Ion slots shift by one when CALIMA adds
    DustBin_01 at ivar 5, so pick the names by presence rather than position."""
    d   = np.asarray(data["density"],  float)
    p   = np.asarray(data["pressure"], float)
    if "scalar_00" in data:                 # plain RT: slots 5,6
        xHI, xHII = np.asarray(data["scalar_00"], float), np.asarray(data["scalar_01"], float)
    else:                                   # CALIMA: slot 5 is dust, ions are 6,7
        xHI, xHII = np.asarray(data["scalar_01"], float), np.asarray(data["scalar_02"], float)
    ul, ut, ud = float(np.ravel(data["unit_l"])[0]), float(np.ravel(data["unit_t"])[0]), float(np.ravel(data["unit_d"])[0])
    scale_v  = ul / ut
    scale_T2 = mH / kB * scale_v**2         # amr/units.f90:32
    # getMu with Y=0: mu = 1/(X*(0.5+0.5 xHI+1.5 xHII)); ->1 neutral, 0.5 ionised, 2 molecular
    mu  = 1.0 / (X * (0.5 + 0.5 * xHI + 1.5 * xHII))
    out = dict(nH   = d * ud / mH * X,
               T    = p / d * scale_T2 * mu,
               xHII = xHII, xHI = xHI,
               x = np.asarray(data["x"], float), y = np.asarray(data["y"], float),
               t = float(np.ravel(data["time"])[0]))
    vx, vy = np.asarray(data["velocity_x"], float), np.asarray(data["velocity_y"], float)
    r = np.hypot(out["x"], out["y"])
    out["vr"] = np.where(r > 0, (vx * out["x"] + vy * out["y"]) / np.maximum(r, 1e-30), 0.0) * scale_v / 1e5
    out["r"]  = r
    out["eps"] = np.asarray(data["DustBin_01"], float) if has_dust else None
    return out

def togrid(f, x, y):
    """The run is levelmin=levelmax=7, so the mesh is a uniform 128^2 grid."""
    n  = int(round(np.sqrt(f.size)))
    dx = 1.0 / n
    ix = np.round(x / dx - 0.5).astype(int)
    iy = np.round(y / dx - 0.5).astype(int)
    img = np.full((n, n), np.nan)
    img[iy, ix] = f
    return img

def radial(f, r, w=None, nbin=64, rmax=1.0):
    edges = np.linspace(0, rmax, nbin + 1)
    idx   = np.digitize(r, edges) - 1
    ok    = (idx >= 0) & (idx < nbin)
    ww    = np.ones_like(f) if w is None else w
    num   = np.bincount(idx[ok], weights=(f * ww)[ok], minlength=nbin)
    den   = np.bincount(idx[ok], weights=ww[ok],       minlength=nbin)
    prof  = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)
    return 0.5 * (edges[1:] + edges[:-1]), prof

# ------------------------------------------------------------------ load
X_rt  = nml_val(os.path.join(RT_DIR,  RT_NML),  "X", 1.0)
X_cal = nml_val(os.path.join(CAL_DIR, CAL_NML), "X", 1.0)
RT  = [derived(load(RT_DIR,  n), X_rt,  False) for n in SNAPS]
CAL = [derived(load(CAL_DIR, n), X_cal, True)  for n in SNAPS]
for a, b in zip(RT, CAL):
    assert abs(a["t"] - b["t"]) < 1e-9, "snapshot times differ: %g vs %g" % (a["t"], b["t"])
times = [a["t"] for a in RT]
print("snapshot times [Myr]:", ", ".join("%.2f" % t for t in times))

# ------------------------------------------------------------------ maps
def maps(runs, cols, fname, title):
    nr, nc = len(runs), len(cols)
    fig, ax = plt.subplots(nr, nc, figsize=(3.05 * nc, 2.95 * nr),
                           squeeze=False, constrained_layout=True)
    for i, R in enumerate(runs):
        for j, (key, lab, cmap, norm, lim) in enumerate(cols):
            a = ax[i][j]
            f = key(R) if callable(key) else R[key]
            img = togrid(f, R["x"], R["y"])
            kw = dict(cmap=cmap, origin="lower", extent=[0, 1, 0, 1])
            if norm == "div":
                # depletion and enhancement each get half the colour range, so
                # a few-per-cent shell is visible next to a total evacuation.
                kw["norm"] = TwoSlopeNorm(vmin=lim[0], vcenter=lim[1], vmax=lim[2])
            elif norm == "log":
                pos = img[np.isfinite(img) & (img > 0)]
                vmin = lim[0] if lim else (np.percentile(pos, 1) if pos.size else 1e-30)
                vmax = lim[1] if lim else (np.percentile(pos, 99.9) if pos.size else 1.0)
                kw["norm"] = LogNorm(vmin=max(vmin, 1e-30), vmax=max(vmax, vmin * 10))
            elif lim:
                kw.update(vmin=lim[0], vmax=lim[1])
            im = a.imshow(img, **kw)
            cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
            if norm == "div":
                # value-proportional bar, so label the narrow enhancement side
                cb.set_ticks([lim[0], 0.5 * (lim[0] + lim[1]), lim[1],
                              lim[1] + 0.5 * (lim[2] - lim[1]), lim[2]])
                cb.ax.tick_params(labelsize=8)
            a.set_xticks([]); a.set_yticks([])
            if i == 0: a.set_title(lab, fontsize=10)
            if j == 0: a.set_ylabel("t = %.2f Myr" % R["t"], fontsize=10)
    fig.suptitle(title, fontsize=12)
    fig.savefig(fname, dpi=115)
    plt.close(fig)
    print("  saved", fname)

nH_lim = (3e-2, 4e-1)
T_lim  = (1e2, 3e4)
maps(RT, [("nH",   r"$n_{\rm H}$ [cm$^{-3}$]", "viridis", "log", nH_lim),
          ("xHII", r"$x_{\rm HII}$",           "inferno", None,  (0, 1)),
          ("T",    r"$T$ [K]",                 "magma",   "log", T_lim)],
     "stromgren2d_maps_rt.png", "2D Stromgren -- plain RAMSES-RT (RT=1, CALIMA=0)")

maps(CAL, [("nH",   r"$n_{\rm H}$ [cm$^{-3}$]", "viridis", "log", nH_lim),
           ("xHII", r"$x_{\rm HII}$",           "inferno", None,  (0, 1)),
           ("T",    r"$T$ [K]",                 "magma",   "log", T_lim),
           (lambda R: R["eps"] / 1e-2, r"$\epsilon/\epsilon_0$  (dust)", "RdBu_r", "div", (0.0, 1.0, 1.15))],
     "stromgren2d_maps_calima.png",
     "2D Stromgren -- RT + CALIMA, 1% dust, drift + radiation pressure")

# difference map: where does the dust actually change the ionisation
fig, ax = plt.subplots(1, len(SNAPS), figsize=(3.4 * len(SNAPS), 3.2), constrained_layout=True)
for j, (A, B) in enumerate(zip(RT, CAL)):
    dimg = togrid(B["xHII"] - A["xHII"], A["x"], A["y"])
    m = np.nanmax(np.abs(dimg)) or 1e-12
    im = ax[j].imshow(dimg, cmap="RdBu_r", vmin=-m, vmax=m, origin="lower", extent=[0, 1, 0, 1])
    fig.colorbar(im, ax=ax[j], fraction=0.046, pad=0.02)
    ax[j].set_title("t = %.2f Myr\nmax|$\\Delta$| = %.1e" % (A["t"], m), fontsize=9)
    ax[j].set_xticks([]); ax[j].set_yticks([])
fig.suptitle(r"$x_{\rm HII}$(RT+CALIMA) $-$ $x_{\rm HII}$(plain RT): the dust effect", fontsize=12)
fig.savefig("stromgren2d_maps_diff.png", dpi=115); plt.close(fig)
print("  saved stromgren2d_maps_diff.png")

# ------------------------------------------------------------------ profiles
panels = [("xHII", r"$x_{\rm HII}$",              False),
          ("xHI",  r"$x_{\rm HI}$",               True ),
          ("T",    r"$T$ [K]",                    True ),
          ("nH",   r"$n_{\rm H}$ [cm$^{-3}$]",    True ),
          ("vr",   r"$v_r$ [km s$^{-1}$]",        False),
          ("eps",  r"$\epsilon/\epsilon_0$  (dust), full range", True ),
          ("epszoom", r"$\epsilon/\epsilon_0$  -- shell zoom",  False),
          ("DxHII", r"$\Delta x_{\rm HII}$  (CALIMA $-$ RT)", False),
          ("DT",    r"$\Delta T$ [K]  (CALIMA $-$ RT)",        False)]
fig, ax = plt.subplots(3, 3, figsize=(14.4, 11.0), constrained_layout=True)
cmap = plt.get_cmap("viridis")
cols = [cmap(v) for v in np.linspace(0.05, 0.85, len(SNAPS))]
for k, (key, lab, logy) in enumerate(panels):
    a = ax.flat[k]
    for i, (A, B) in enumerate(zip(RT, CAL)):
        if key.startswith("D"):                      # residual panels
            base = key[1:]
            r, pa = radial(A[base], A["r"]); _, pb = radial(B[base], B["r"])
            a.plot(r, pb - pa, "-", color=cols[i], lw=1.5, label="t = %.2f Myr" % B["t"])
            a.axhline(0.0, color="0.6", lw=0.8)
            continue
        if key.startswith("eps"):
            r, p = radial(B["eps"], B["r"]); p = p / 1e-2
        else:
            r, pa = radial(A[key], A["r"]); a.plot(r, pa, "--", color=cols[i], lw=1.3)
            r, p  = radial(B[key], B["r"])
        a.plot(r, p, "-", color=cols[i], lw=1.7, label="t = %.2f Myr" % B["t"])
        if key == "epszoom":
            k = np.nanargmax(np.where(r < 0.95, p, np.nan))
            a.plot(r[k], p[k], "o", color=cols[i], ms=4)
    a.set_xlabel("r from source [box units]"); a.set_ylabel(lab)
    if logy: a.set_yscale("log")
    a.set_xlim(0, 1.0)
    if k == 0:
        a.legend(fontsize=8, loc="lower left")
        a.text(0.97, 0.95, "solid: RT+CALIMA\ndashed: plain RT", transform=a.transAxes,
               ha="right", va="top", fontsize=8)
    if key == "xHI":
        a.set_ylim(1e-12, 2.0)
    if key == "eps":
        a.set_ylim(2e-3, 1.5)
        a.axhline(1.0, color="0.5", lw=0.9, ls=":")
        a.text(0.03, 1.05, "initial 1% dust", fontsize=8, color="0.4")
    if key == "epszoom":
        # the pile-up is only a few per cent, so it needs its own scale
        a.set_ylim(0.975, 1.12)
        a.axhline(1.0, color="0.5", lw=0.9, ls=":")
        a.text(0.03, 1.103, "dots mark the shell peak", fontsize=8, color="0.4")
    if key == "DxHII":
        a.text(0.5, 0.95, "plain-RT dashed curves in the other panels lie\n"
                          "under the solid ones; this is the actual difference",
               transform=a.transAxes, ha="center", va="top", fontsize=7.5, color="0.35")
fig.suptitle("2D Stromgren, radial profiles from the source at (0,0) -- "
             "solid: RT+CALIMA (1% dust), dashed: plain RT", fontsize=12)
fig.savefig("stromgren2d_profiles.png", dpi=115); plt.close(fig)
print("  saved stromgren2d_profiles.png")

# ------------------------------------------------------------------ numbers
print("\n  t[Myr]   R_ion(RT)  R_ion(CAL)   <T>_ion RT   <T>_ion CAL   max|dxHII|   min eps      max eps")
for A, B in zip(RT, CAL):
    def rion(R):
        r, p = radial(R["xHII"], R["r"])
        ok = np.isfinite(p) & (p > 0.5)
        return r[ok].max() if ok.any() else np.nan
    ia, ib = A["xHII"] > 0.5, B["xHII"] > 0.5
    print("  %6.2f   %9.4f  %9.4f   %10.1f   %11.1f   %10.3e   %9.3e  %9.3e"
          % (A["t"], rion(A), rion(B),
             A["T"][ia].mean() if ia.any() else np.nan,
             B["T"][ib].mean() if ib.any() else np.nan,
             np.abs(B["xHII"] - A["xHII"]).max(), B["eps"].min(), B["eps"].max()))
