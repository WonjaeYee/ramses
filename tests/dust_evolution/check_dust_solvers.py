#!/usr/bin/env python3
"""Accuracy of a CALIMA dust solver test (calima/dust_test.f90) against the exact solution of its rate laws.

    check_dust_solvers.py <dust_*_test.log> [--ref <log>] [--tol 1e-2]

The test holds the gas state fixed, so each process has a rate law with constant coefficients. The log's header
gives the rates of every process at t = 0 ('# Rate ...'), from which the coefficients follow:
  linear loss  (dust sputtering, sublimation, PAH sputtering, photolysis)   y' = -k y,  k = loss / y(0)
  cluster evaporation   y2' = -k y2, y1' = +k y2                             (bin 2 into bin 1)
  coalescence           y1' = -a y1^2, y2' = +a y1^2,  a = loss / y1(0)^2
  both together         a Riccati equation in y1 (y1 + y2 constant), with its equilibrium
  freezing              y_p' = -y_p sum_d c_pd y_d,  y_d' = y_d sum_p c_pd y_p   (c_pd from each grain's gain)
  accretion             y_i' = K_i L(t) y_i, L = min_e n_e / (f_e sqrt(m_e)) of the bin's elements, which the
                        growth depletes (assumes a constant sticking factor: T_dust >= 100 K, or no ice change)
The first four in closed form, the last two integrated to 1e-12 (scipy). Coagulation and shattering have no
closed form: give a converged run of the same test as --ref (e.g. RK54 with errmax = 1e-8 at a small dt).
Reports, per bin holding more than 1e-6 of the dust + PAH mass, the largest relative error over the outputs, the
relative change of the total (gas + grains) mass, and the cost; exit 1 if an error exceeds --tol.
Species: 1..npah are the PAH bins, npah+1..npah+ndust the dust bins (the solver's order).
"""
import argparse
import re
import sys

import numpy as np

NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"
LINEAR = {"sputtering", "sublimation", "photolysis"}


def parse(path):
    h = {"rates": [], "bin_el": {}, "gas_mass": {}, "pah": {}}
    rows = []
    for line in open(path):
        if not line.strip():
            continue
        if line.startswith("#"):
            m = re.match(rf"# Rate (dust|pah) (\S+) species (\d+) gain =\s+({NUM}) loss =\s+({NUM})", line)
            if m:
                h["rates"].append((m.group(1), m.group(2), int(m.group(3)), float(m.group(4)), float(m.group(5))))
                continue
            m = re.match(rf"# GasElement (\d+) atomic_number =\s+(\d+) mass_g =\s+({NUM})", line)
            if m:
                h["gas_mass"][int(m.group(1))] = float(m.group(3))
                continue
            m = re.match(rf"# Bin (\d+) element (\d+) index =\s+(\d+) mfraction =\s+({NUM}) mass_g =\s+({NUM})", line)
            if m:
                h["bin_el"].setdefault(int(m.group(1)), []).append((int(m.group(3)), float(m.group(4)), float(m.group(5))))
                continue
            m = re.match(rf"# PAH (\d+) mpah =\s+({NUM}) apah_cm =\s+({NUM}) nc =\s+(\d+) is_cluster =\s+([TF])", line)
            if m:
                h["pah"][int(m.group(1))] = dict(mpah=float(m.group(2)), cluster=m.group(5) == "T")
                continue
            m = re.match(r"# ([A-Za-z_]+) =\s+(\S+)", line)
            if m and m.group(1) not in ("Time",):
                v = m.group(2)
                try:
                    v = float(v) if v not in ("T", "F") else v == "T"
                except ValueError:
                    pass
                h[m.group(1)] = v
            m = re.match(rf"# ode_accepted =\s+(\d+) rejected =\s+(\d+) reduced =\s+(\d+)", line)
            if m:
                h["ode"] = tuple(int(x) for x in m.groups())
            continue
        rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    nd, npah = int(h["ndust"]), int(h["npah"])
    nel = a.shape[1] - 2 - nd - npah
    h["t"] = a[:, 0] * 3.15576e7
    h["gas"] = a[:, 2:2 + nel]
    h["y"] = np.hstack([a[:, 2 + nel + nd:], a[:, 2 + nel:2 + nel + nd]])   # solver order: PAH, then dust
    h["nel"], h["nd"], h["npah"] = nel, nd, npah
    return h


def total_mass(h):
    m = np.array([h["gas_mass"].get(i + 1, 0.0) for i in range(h["nel"])])
    return (h["gas"] * m).sum(axis=1) + h["y"].sum(axis=1)


def coeffs(h):
    """per process name: arrays of gain and loss per species at t = 0"""
    n = h["nd"] + h["npah"]
    c = {}
    for kind, name, s, g, l in h["rates"]:
        key = (kind, name)
        c.setdefault(key, (np.zeros(n), np.zeros(n)))
        c[key][0][s - 1] += g
        c[key][1][s - 1] += l
    return c


def exact(h):
    """the exact solution at the output times, or None if the active processes have no exact form here"""
    t, y0 = h["t"], h["y"][0]
    npah = h["npah"]
    c = coeffs(h)
    names = {k[1] if k[0] == "pah" else "dust_" + k[1] for k in c}
    y = np.tile(y0, (len(t), 1))
    lin = {"dust_sputtering", "dust_sublimation", "sputtering", "photolysis"}
    if names <= lin:
        for key, (g, l) in c.items():
            for s in np.nonzero(l > 0)[0]:
                y[:, s] *= np.exp(-l[s] / y0[s] * t)
        return y
    if names == {"cluster_evaporation"}:
        g, l = c[("pah", "cluster_evaporation")]
        for s in np.nonzero(l > 0)[0]:
            k = l[s] / y0[s]
            y[:, s] = y0[s] * np.exp(-k * t)
            y[:, s - 1] = y0[s - 1] + y0[s] * -np.expm1(-k * t)
        return y
    if names == {"coalescence"}:
        g, l = c[("pah", "coalescence")]
        for s in np.nonzero(l > 0)[0]:
            a = l[s] / y0[s] ** 2
            y[:, s] = y0[s] / (1 + a * y0[s] * t)
            y[:, s + 1] = y0[s + 1] + (y0[s] - y[:, s])
        return y
    if names == {"coalescence", "cluster_evaporation"} and npah == 2:
        a = c[("pah", "coalescence")][1][0] / y0[0] ** 2
        k = c[("pah", "cluster_evaporation")][1][1] / y0[1]
        # a (y1 - y+)(y1 - y-) with y+ the equilibrium; written without cancellation for k >> a S (stiff
        # evaporation), where bin 2's equilibrium y2* = 4 a k S^2 / (k + disc)^2 is ~ a S / k of the total
        S = y0[0] + y0[1]
        disc = np.sqrt(k * k + 4 * a * k * S)
        D = disc / a                                   # y+ - y-
        y2s = 4 * a * k * S * S / (k + disc) ** 2       # S - y+
        w0 = y0[1] - y2s                               # bin 2 above its equilibrium
        e = np.exp(-disc * t)
        w = D * w0 * e / (D - w0 * -np.expm1(-disc * t))
        y[:, 1] = y2s + w
        y[:, 0] = S - y[:, 1]
        return y
    if names in ({"freezing"}, {"dust_accretion"}):
        from scipy.integrate import solve_ivp
        if names == {"freezing"}:
            g, l = c[("pah", "freezing")]
            P = [p for p in range(npah) if l[p] > 0]
            D = [d for d in range(npah, len(y0)) if g[d] > 0]
            # c_pd: each grain's gain is sum_p c_pd y_p y_d; with one freezing PAH bin it is exact
            if len(P) != 1:
                return None
            p = P[0]
            cpd = {d: g[d] / (y0[p] * y0[d]) for d in D}

            def f(_, z):
                r = np.zeros_like(z)
                for d in D:
                    r[d] = cpd[d] * z[p] * z[d]
                r[p] = -sum(r[d] for d in D)
                return r
        else:
            g, l = c[("dust", "accretion")]
            nd = h["nd"]
            mass = h["gas_mass"]
            rho0 = h["gas"][0] * np.array([mass.get(i + 1, 0.0) for i in range(h["nel"])])  # g cm-3
            # bins of one composition share the gas: group them by their element list
            groups = {}
            for b in range(1, nd + 1):
                groups.setdefault(tuple(e for e, _, _ in h["bin_el"][b]), []).append(b)

            def L(rho, b):
                return min(rho[e - 1] / (f * np.sqrt(m)) for e, f, m in h["bin_el"][b])

            K = {b: g[npah + b - 1] / (y0[npah + b - 1] * L(rho0, b)) if y0[npah + b - 1] > 0 else 0.0
                 for b in range(1, nd + 1)}

            def f(_, z):
                r = np.zeros_like(z)
                for els, bins in groups.items():
                    grow = sum(z[npah + b - 1] - y0[npah + b - 1] for b in bins)
                    rho = rho0.copy()
                    for e, fr, _ in h["bin_el"][bins[0]]:
                        rho[e - 1] = rho0[e - 1] - fr * grow
                    for b in bins:
                        r[npah + b - 1] = K[b] * max(L(rho, b), 0.0) * z[npah + b - 1]
                return r
        sol = solve_ivp(f, (0, t[-1]), y0, method="Radau", t_eval=t, rtol=1e-12, atol=1e-40 * max(y0.max(), 1e-300))
        return sol.y.T if sol.success else None
    return None


def compare(h, ref):
    y = h["y"]
    keep = y[0] > 1e-6 * y[0].sum() if y[0].sum() > 0 else np.ones(y.shape[1], bool)
    keep |= ref.max(axis=0) > 1e-6 * ref[0].sum()
    err = np.abs(y - ref) / np.maximum(np.abs(ref), 1e-6 * ref[0].sum())
    return {int(s): float(err[1:, s].max()) for s in np.nonzero(keep)[0]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("log")
    ap.add_argument("--ref", help="a converged log of the same test (coagulation, shattering)")
    ap.add_argument("--tol", type=float, default=1e-2)
    a = ap.parse_args()
    h = parse(a.log)
    if a.ref:
        r = parse(a.ref)
        ref = np.array([np.interp(h["t"], r["t"], r["y"][:, s]) for s in range(r["y"].shape[1])]).T
        kind = "reference " + a.ref
    else:
        ref = exact(h)
        kind = "exact rate-law solution"
        if ref is None:
            print(f"{a.log}: no exact solution for these processes ({sorted({k[1] for k in coeffs(h)})}); give --ref")
            return 2
    e = compare(h, ref)
    M = total_mass(h)
    lab = lambda s: f"PAH {s + 1}" if s < h["npah"] else f"dust {s - h['npah'] + 1}"  # noqa: E731
    print(f"{a.log}: solver {int(h['dust_solver_type'])}, dt = {h['dt_yr']:.3g} yr, {len(h['t']) - 1} steps; against the {kind}")
    for s, v in e.items():
        print(f"   {lab(s):8s} max relative error {v:.3e}")
    print(f"   total mass change {abs(M[-1] - M[0]) / M[0]:.2e}; cpu {h.get('cpu_s', float('nan')):.3g} s"
          + (f"; ODE substeps {h['ode'][0]} accepted, {h['ode'][1]} rejected" if "ode" in h else ""))
    bad = [lab(s) for s, v in e.items() if v > a.tol]
    print("   PASS" if not bad else f"   FAIL (> {a.tol:g}): {', '.join(bad)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
