"""Readers shared by the sink_indi plot scripts (tests/sink_indi).

load(nout)      cells of output_<nout>: hydro and RT variables, cgs helpers, the sinks
read_log(path)  the sink table printed every ncontrol steps (INDIVIDUAL_SINK_STARS format)
"""
import os
import re

import numpy as np

import visu_ramses

MSUN = 1.98847e33
YR = 3.15576e7
KYR = 1e3*YR
PC = 3.0856776e18
MH = 1.6726e-24
KB = 1.380649e-16
G_CGS = 6.674e-8
C_CGS = 2.99792458e10


def _load_rt(nout):
    """The RT file has the hydro file's record layout: read it with visu_ramses by swapping
    the file prefix and the descriptor."""
    fname, desc = visu_ramses.generate_fname, visu_ramses.read_descriptor
    visu_ramses.generate_fname = lambda n, ftype="", cpuid=1: fname(n, "rt" if ftype == "hydro" else ftype, cpuid)
    visu_ramses.read_descriptor = lambda f: desc(f.replace("hydro_file_descriptor", "rt_file_descriptor"))
    try:
        return visu_ramses.load_snapshot(nout)
    finally:
        visu_ramses.generate_fname, visu_ramses.read_descriptor = fname, desc


def load(nout, rt=True):
    """Cells of output_<nout>. Adds r_pc (distance to sink 1), vol_cm3, nH (pure H: rho/m_H),
    and, with rt, Np_<g> [photons cm^-3] for every group and info["c_red"] [cm/s]."""
    d = visu_ramses.load_snapshot(nout)
    c, info = d["data"], d["info"]
    if rt:
        r = _load_rt(nout)
        for k, v in r["data"].items():
            c.setdefault(k, v)
        with open(f"output_{nout:05d}/info_rt_{nout:05d}.txt") as f:
            txt = f.read()
        unit_pf = float(re.search(r"unit_pf\s*=\s*(\S+)", txt).group(1))
        info["c_red"] = float(re.search(r"rt_c_frac\s*=\s*(\S+)", txt).group(1))*C_CGS
        # photon_flux_<g> is the photon density in flux units, rt_c*Np (rt_output_hydro)
        g = 1
        while f"photon_flux_{g:02d}" in c:
            c[f"Np_{g}"] = c[f"photon_flux_{g:02d}"]*unit_pf/info["c_red"]
            g += 1
        info["ngroups"] = g - 1
    ul, ud, ut = info["unit_l"], info["unit_d"], info["unit_t"]
    s = d["sinks"]
    xs = np.array([s["x"][0], s["y"][0], s["z"][0]]) if s["nsinks"] > 0 else 0.5*info["boxlen"]*np.ones(3)
    c["r_pc"] = np.sqrt((c["x"] - xs[0])**2 + (c["y"] - xs[1])**2 + (c["z"] - xs[2])**2)*ul/PC
    c["vol_cm3"] = (c["dx"]*ul)**3
    c["nH"] = c["density"]*ud/MH
    info["t_kyr"] = info["time"]*ut/KYR
    return d


def outputs():
    """Output numbers present in the working directory, in order."""
    return sorted(int(d[7:]) for d in os.listdir(".") if re.fullmatch(r"output_\d{5}", d))


def read_sink_csv(nout):
    """Sink 1 of output_<nout>/sink_<nout>.csv as a dict of floats (code units), plus t (code)."""
    with open(f"output_{nout:05d}/sink_{nout:05d}.csv") as f:
        lines = f.readlines()
    keys = lines[0].replace("#", "").replace(" ", "").strip().split(",")
    s = dict(zip(keys, (float(v) for v in lines[2].split(","))))
    with open(f"output_{nout:05d}/info_{nout:05d}.txt") as f:
        info = dict(l.split("=", 1) for l in f if "=" in l)
    s.update({k.strip(): float(info[k]) for k in info if k.strip() in ("time", "unit_l", "unit_d", "unit_t")})
    return s


def find_log(test):
    """Text of this test's RAMSES log: run.log if present, else this test's section of the
    test suite log (run_test_suite.sh appends every run to tests/test_suite.log)."""
    if os.path.exists("run.log"):
        with open("run.log", errors="replace") as f:
            return f.read()
    suite = os.path.join("..", "..", "test_suite.log")
    with open(suite, errors="replace") as f:
        txt = f.read()
    start = txt.rfind(f": sink_indi/{test}\n")
    return txt[start:] if start >= 0 else ""


_ROW = re.compile(r"^\s+\d+\s+\d+(\s+[-+0-9.E]+){14,}\s*$")


def read_log(txt):
    """Sink table rows of a RAMSES log text, one per printout, with the simulation time printed
    just before it. Returns a dict of arrays: t_yr, M, M_F [Msun], Mdot [Msun/yr], age_Myr,
    rho [g/cm3], Q_sub, Q_ion [s^-1] (sink 1 only)."""
    keys = ("M", "M_F", "x", "y", "z", "vx", "vy", "vz", "spin", "Mdot", "age_Myr", "rho", "Q_sub", "Q_ion")
    out = {k: [] for k in ("t_yr",) + keys}
    t = np.nan
    for line in txt.splitlines():
        if "simulation time [yr]" in line:
            t = float(line.split("=")[1])
        elif _ROW.match(line) and line.split()[0] == "1":
            v = [float(x) for x in line.split()[2:2 + len(keys)]]
            out["t_yr"].append(t)
            for k, x in zip(keys, v):
                out[k].append(x)
    return {k: np.array(v) for k, v in out.items()}


def report(name, checks):
    """Print each check and the suite's PASSED/FAILED line. checks: (label, value, ok)."""
    for label, value, ok in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {label}: {value}")
    print(f"Test {name}: {'PASSED' if all(ok for _, _, ok in checks) else 'FAILED'}")
