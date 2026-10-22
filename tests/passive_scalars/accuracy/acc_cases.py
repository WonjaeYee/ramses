"""Passive-scalar advection accuracy: profiles, exact cell averages, the case matrix, run directories.

Three scalars ride on a uniform state (rho = 1, c_s = 1, uniform velocity, periodic unit box): a top-hat
(|s - 0.5| < 0.15), a Gaussian (sigma = 0.05, centred at 0.5) and a sine (0.5 + 0.5 sin 2 pi s).  The profiles
are functions of s = kx x + ky y (k = (1, 0) aligned with the grid, (1, 1) the diagonal); the velocity is
along k, and the run ends when the profile has moved by an integer number of periods in s, so the exact
solution is the initial one.  Both the initial table and the exact solution are exact cell averages
(antiderivatives of the profiles, no quadrature).
"""
import math
import os

import numpy as np
from scipy.special import erf

PROFILES = ("tophat", "gauss", "sine")
TOP_LO, TOP_HI = 0.35, 0.65
G_C, G_SIG = 0.5, 0.05
GAMMA = 1.4
BETA_MHD = 10.0          # plasma beta of the uniform MHD field (transverse to the flow)


# ---------------------------------------------------------------- antiderivatives on the unwrapped line, s >= 0
def _F1(name, s):
    s = np.asarray(s, float)
    if name == "tophat":
        m, r = np.floor(s), s - np.floor(s)
        return (TOP_HI - TOP_LO) * m + np.clip(r - TOP_LO, 0.0, TOP_HI - TOP_LO)
    if name == "gauss":
        a = G_SIG * math.sqrt(math.pi / 2)
        out = np.zeros_like(s)
        for n in range(-2, int(np.ceil(s.max())) + 3):
            out += a * (erf((s - G_C - n) / (math.sqrt(2) * G_SIG)) + 1.0)
        return out
    if name == "sine":
        return 0.5 * s - 0.5 * np.cos(2 * math.pi * s) / (2 * math.pi)
    raise KeyError(name)


def _F2(name, s):
    s = np.asarray(s, float)
    if name == "tophat":
        w = TOP_HI - TOP_LO
        m, r = np.floor(s), s - np.floor(s)
        G = np.where(r < TOP_LO, 0.0, np.where(r < TOP_HI, 0.5 * (r - TOP_LO) ** 2,
                                                0.5 * w * w + w * (r - TOP_HI)))
        return 0.5 * w * m * m + w * m * r + G
    if name == "gauss":
        a = G_SIG * math.sqrt(math.pi / 2)
        b = math.sqrt(2) * G_SIG
        out = np.zeros_like(s)
        for n in range(-2, int(np.ceil(s.max())) + 3):
            z = (s - G_C - n) / b
            out += a * b * (z * erf(z) + np.exp(-z * z) / math.sqrt(math.pi) + z)
        return out
    if name == "sine":
        return 0.25 * s * s - 0.5 * np.sin(2 * math.pi * s) / (4 * math.pi ** 2)
    raise KeyError(name)


def cell_average(name, x0, h, k, ndim, shift=0.0):
    """Average of f(k.x - shift) over the cells [x0, x0 + h] (x0: (n, ndim)), k integer, k[0] > 0."""
    x0 = np.atleast_2d(np.asarray(x0, float))
    kx = k[0]
    ky = k[1] if ndim > 1 else 0
    off = math.ceil(abs(shift)) + 2 - shift          # integer shift + offset keeps s >= 0, periodicity exact
    if ndim == 1 or ky == 0:
        a = kx * x0[:, 0] + off
        return (_F1(name, a + kx * h) - _F1(name, a)) / (kx * h)
    a = kx * x0[:, 0] + ky * x0[:, 1] + off
    return (_F2(name, a + kx * h + ky * h) - _F2(name, a + ky * h) - _F2(name, a + kx * h) + _F2(name, a)) \
        / (kx * ky * h * h)


def exact_on_cells(xc, h, k, ndim, shift=0.0):
    """Exact averages of the three profiles on cells of centres xc (n, ndim) and size h: dict name -> (n,)."""
    x0 = np.atleast_2d(xc) - 0.5 * h
    out = {p: cell_average(p, x0, h, k, ndim, shift) for p in PROFILES}
    for p in ("tophat", "gauss"):                    # non-negative: drop the round-off of the differences
        out[p] = np.maximum(out[p], 0.0)
    return out


def ic_rows(N, k, ndim):
    """s values a cell centre takes, and the exact cell averages there (one table row each)."""
    h = 1.0 / N
    ksum = k[0] + (k[1] if ndim > 1 else 0)
    if ndim == 1 or ksum == k[0]:
        s = (np.arange(N) + 0.5) * h * k[0] % 1.0
        s = np.unique(np.round(s / (0.5 * h)) * 0.5 * h)
    else:
        s = np.arange(2 * N + 1) * 0.5 * h                      # s_c = h (kx(i+1/2) + ky(j+1/2)) mod 1
    xc = np.zeros((s.size, ndim))
    xc[:, 0] = s / ksum
    if ndim > 1 and ksum != k[0]:
        xc[:, 1] = s / ksum
    return s, exact_on_cells(xc, h, k, ndim)


# ---------------------------------------------------------------- cases
def case_name(c):
    g = "" if c["k"] == (1, 0) else f"_k{c['k'][0]}{c['k'][1]}"
    cf = "" if c["courant"] == 0.8 else f"_C{c['courant']:g}"
    pe = "" if c["periods"] == 1 else f"_P{c['periods']}"
    return (f"{c['solver'][0]}{c['ndim']}d_s{c['slope']}_{c['riemann']}_M{c['mach']:g}_N{c['N']}"
            f"{g}{cf}{pe}")


def _case(solver, ndim, slope, riemann, mach, N, k=(1, 0), courant=0.8, periods=1):
    # slope_type 8 (generalised MC) uses slope_theta = 1.5: it is not in the hydro_params namelist
    return dict(solver=solver, ndim=ndim, slope=slope, riemann=riemann, mach=mach, N=N, k=tuple(k),
                courant=courant, periods=periods)


NS_1D = (32, 64, 128, 256, 512)
NS_2D = (32, 64, 128)
HYDRO_SLOPES_1D = (0, 1, 2, 4, 5, 6, 7, 8)       # 3 = 2 in 1-D; 6 is the unlimited (unstable) slope
HYDRO_RIEMANN = ("llf", "hll", "hllc", "exact", "acoustic")
HYDRO_SLOPES_2D = (0, 1, 2, 3, 7, 8)             # 4-6 exist in 1-D only
MHD_SLOPES_1D = (0, 1, 2)                         # the MHD 1-D slopes
MHD_SLOPES_2D = (0, 1, 2, 3)
MHD_RIEMANN = ("llf", "hll", "roe", "hlld", "upwind")


def matrix(name):
    out = []
    if name == "full":
        for s in HYDRO_SLOPES_1D:
            for r in HYDRO_RIEMANN:
                for m in (0.1, 1.0, 10.0):
                    out += [_case("hydro", 1, s, r, m, N) for N in NS_1D]
        for s in (1, 2, 7, 8):                                   # Courant number and many periods
            out += [_case("hydro", 1, s, "hllc", 1.0, N, courant=0.4) for N in NS_1D]
            out += [_case("hydro", 1, s, "hllc", 1.0, N, periods=10) for N in (64, 256)]
        for s in (1, 2, 4, 7, 8):                                # new extrema with the non-upwind fluxes
            for r in ("hll", "llf"):
                for cf in (0.2, 0.4):
                    out += [_case("hydro", 1, s, r, 0.1, N, courant=cf) for N in (64, 256)]
        for s in (1, 2):
            out += [_case("mhd", 1, s, "llf", 1.0, N, courant=cf) for cf in (0.2, 0.4) for N in (64, 256)]
        for s in HYDRO_SLOPES_2D:
            for r in ("hllc", "hll"):
                out += [_case("hydro", 2, s, r, 1.0, N, k=(1, 1)) for N in NS_2D]
            out += [_case("hydro", 2, s, "hllc", 1.0, N, k=(1, 0)) for N in NS_2D]
        for s in MHD_SLOPES_1D:
            for r in MHD_RIEMANN:
                for m in (0.1, 1.0):
                    out += [_case("mhd", 1, s, r, m, N) for N in NS_1D]
        for s in MHD_SLOPES_2D:
            for r in ("hlld", "llf"):
                out += [_case("mhd", 2, s, r, 1.0, N, k=(1, 1)) for N in NS_2D]
    elif name == "nightly":
        for s in (0, 1, 2, 7, 8):
            for r in ("hllc", "hll"):
                for m in (0.1, 1.0):
                    out += [_case("hydro", 1, s, r, m, N) for N in (64, 256)]
        for s in (1, 2, 3):
            out += [_case("hydro", 2, s, "hllc", 1.0, 64, k=(1, 1))]
        for s in (1, 2):
            for r in ("hlld", "llf"):
                out += [_case("mhd", 1, s, r, 1.0, N) for N in (64, 256)]
        out += [_case("mhd", 2, 2, "hlld", 1.0, 64, k=(1, 1))]
    else:
        raise KeyError(name)
    return out


# ---------------------------------------------------------------- run directory
NML = """! tests/passive_scalars/accuracy: {name}
! three passive scalars (top-hat, Gaussian, sine) carried by a uniform flow through a periodic box
&RUN_PARAMS
hydro=.true.
poisson=.false.
nrestart=0
ncontrol=100000
nstepmax=100000000
verbose=.false.
/
&AMR_PARAMS
levelmin={level}
levelmax={level}
ngridmax={ngridmax}
nexpand=1
boxlen=1.0
/
&INIT_PARAMS
condinit_kind='acc_table'
/
&OUTPUT_PARAMS
foutput=100000000
noutput=1
tout={tend:.17g}
/
&HYDRO_PARAMS
gamma={gamma}
courant_factor={courant}
scheme='muscl'
slope_type={slope}
riemann='{riemann}'
{riemann2d}pressure_fix=.false.
/
&COOLING_PARAMS
cooling=.false.
/
"""


def velocity(c):
    """Flow speed (Mach x c_s) along k, and the time for an integer number of periods in s."""
    kx, ky = c["k"] if c["ndim"] > 1 else (c["k"][0], 0)
    kn = math.hypot(kx, ky)
    v = c["mach"] * 1.0
    vel = (v * kx / kn, v * ky / kn, 0.0)
    return vel, c["periods"] / (v * kn)


def prepare(c, rundir):
    os.makedirs(rundir, exist_ok=True)
    N, ndim = c["N"], c["ndim"]
    level = int(round(math.log2(N)))
    assert 2 ** level == N
    s, avg = ic_rows(N, c["k"], ndim)
    rho, P = 1.0, 1.0 / GAMMA
    vel, tend = velocity(c)
    cols = [s, np.full_like(s, rho), np.full_like(s, vel[0]), np.full_like(s, P)]
    if c["solver"] == "mhd":
        B = math.sqrt(2 * P / BETA_MHD)        # P_mag = B^2/2
        Bxyz = (0.0, B, 0.0) if ndim == 1 else (0.0, 0.0, B)   # transverse to the flow (out of plane in 2-D)
        cols += [np.full_like(s, b) for b in Bxyz]
    cols += [avg[p] for p in PROFILES]
    tab = np.column_stack(cols)
    with open(os.path.join(rundir, "ic_table.dat"), "w") as f:
        f.write(f"{tab.shape[0]} {tab.shape[1]}\n")
        np.savetxt(f, tab, fmt="%.17e")
    kx, ky = c["k"] if ndim > 1 else (c["k"][0], 0)
    with open(os.path.join(rundir, "ic_geom.dat"), "w") as f:
        f.write(f"{kx} {ky} 0 {vel[0]:.17e} {vel[1]:.17e} {vel[2]:.17e}\n")
    ngridmax = 2 * (N ** ndim) // (2 ** ndim) + 2000
    r2d = f"riemann2d='{c['riemann']}'\n" if (c["solver"] == "mhd" and ndim > 1) else ""
    if c["solver"] == "mhd":     # slope_mag_type defaults to slope_type and has no 3; B is uniform anyway
        r2d += f"slope_mag_type={min(c['slope'], 2)}\n"
    with open(os.path.join(rundir, "acc.nml"), "w") as f:
        f.write(NML.format(name=case_name(c), level=level, ngridmax=ngridmax, tend=tend, gamma=GAMMA,
                           courant=c["courant"], slope=c["slope"], riemann=c["riemann"],
                           riemann2d=r2d))
    import json
    json.dump(dict(c, k=list(c["k"]), tend=tend, vel=list(vel)), open(os.path.join(rundir, "case.json"), "w"),
              indent=1)
    return tend
