"""Shared definitions for the passive-scalar conservation tests (tests/passive_scalars).

* SPECIES_SETS: element/ion/molecule/dust selections, each with the Makefile variables of
  the RAMSES-RTZ(/CALIMA) build that carries exactly that layout.
* Layout: the passive-scalar index layout of such a build (RAMSES order: elements, CO,
  PAH bins, dust bins, ion blocks per element, H2), with the constraint tree
  (top species sum to rho, ions/H2 sum to their element).
* compositions / profiles: initial conditions, built as PHYSICAL mixtures of a few gas
  parcels (so that the same physical state is stored in either ion convention).
* CASES: the parameterised test matrix (one definition per test, many cases).
* write_ic_table / write_namelist: RAMSES inputs for a case (patch/condinit.f90 reads
  ic_table.dat when condinit_kind='table').

Natural fractions ("fnat") are used everywhere outside the solver: top species as mass
fractions of rho, ion stages (and H2, as the fraction of H nuclei in H2) as fractions of
their element.  That is also what RAMSES writes to its outputs for the ions in both
conventions (slot/rho in the x*rho convention, slot/rho_element in the x*rho_element one).
"""
import itertools
import os
import numpy as np

ELEMENTS = [  # symbol, Z, RAMSES name, solar mass fraction (approx.)
    ("H", 1, "HYDROGEN", 0.7381), ("He", 2, "HELIUM", 0.2485), ("C", 6, "CARBON", 2.37e-3),
    ("N", 7, "NITROGEN", 6.93e-4), ("O", 8, "OXYGEN", 5.73e-3), ("Ne", 10, "NEON", 1.26e-3),
    ("Mg", 12, "MAGNESIUM", 7.08e-4), ("Si", 14, "SILICON", 6.65e-4), ("S", 16, "SULFUR", 3.10e-4),
    ("Fe", 26, "IRON", 1.29e-3)]
MAKE_ION_VAR = {"H": "N_HYDROGEN_IONS", "He": "N_HELIUM_IONS", "C": "N_CARBON_IONS",
                "N": "N_NITROGEN_IONS", "O": "N_OXYGEN_IONS", "Ne": "N_NEON_IONS",
                "Mg": "N_MAGNESIUM_IONS", "Si": "N_SILICON_IONS", "S": "N_SULFUR_IONS",
                "Fe": "N_IRON_IONS"}
ALL = [e[0] for e in ELEMENTS]

SPECIES_SETS = {
    # name: elements, ion stages per element (default Z+1), H2, CO, CALIMA bins
    "H":        dict(elements=["H"], h2=False, co=False, ndust=0, npah=0),
    "HHe_H2":   dict(elements=["H", "He"], h2=True, co=False, ndust=0, npah=0),
    "HHeCO":    dict(elements=["H", "He", "C", "O"], nions={"C": 3, "O": 3}, h2=True, co=True,
                     ndust=4, npah=2, ndchemtype=1),
    "full":     dict(elements=ALL, h2=True, co=True, ndust=4, npah=2, ndchemtype=2),
    "full_np0": dict(elements=ALL, h2=False, co=False, ndust=4, npah=0, ndchemtype=2),
    # 3D sink test (individual stars need all ten elements; no CALIMA, to keep it light)
    "full_nodust": dict(elements=ALL, h2=True, co=True, ndust=0, npah=0),
}


def make_vars(setname, ndim=1):
    """Makefile variable overrides for the RAMSES build of a species set."""
    s = SPECIES_SETS[setname]
    z = {e[0]: e[1] for e in ELEMENTS}
    v = dict(COMPILER="GNU", MPIF90="mpif90", SOLVER="hydro", USE_TURB="0", INDI_STAR="0",
             NDIM=str(ndim), RT="0", RTZ="1", NENER="0", NGROUPS="1", NPSCAL="0",
             NMETALS=str(len(s["elements"])), N_H2="1" if s["h2"] else "0",
             CO="1" if s["co"] else "0")
    for e, var in MAKE_ION_VAR.items():
        v[var] = str((s.get("nions") or {}).get(e, z[e] + 1)) if e in s["elements"] else "0"
    if s["ndust"] + s["npah"] > 0:
        v.update(CALIMA="1", NDUST=str(s["ndust"]), NPAH=str(s["npah"]),
                 NDCHEMTYPE=str(s.get("ndchemtype", 1)))
    else:
        v.update(CALIMA="0")
    return v


class Layout:
    def __init__(self, setname):
        s = SPECIES_SETS[setname]
        z = {e[0]: e[1] for e in ELEMENTS}
        self.setname = setname
        self.names, self.parent, self.kind, self.rnames = [], [], [], []
        rname = {e[0]: e[2] for e in ELEMENTS}
        self.elem_idx = {}
        for e in s["elements"]:
            self.elem_idx[e] = len(self.names)
            self._add(e, -1, "elem", rname[e])
        if s["co"]:
            self._add("CO", -1, "CO", "CO")
        for j in range(s["npah"]):
            self._add(f"PAH{j+1}", -1, "pah", f"PAHBin_{j+1:02d}")
        for j in range(s["ndust"]):
            self._add(f"dust{j+1}", -1, "dust", f"DustBin_{j+1:02d}")
        self.ion_idx = {}
        for e in s["elements"]:
            n = (s.get("nions") or {}).get(e, z[e] + 1)
            self.ion_idx[e] = []
            for j in range(n):
                self.ion_idx[e].append(len(self.names))
                self._add(f"{e}{j+1}", self.elem_idx[e], "ion", f"{e}_{j+1:02d}")
        if s["h2"]:
            self.ion_idx["H"].append(len(self.names))
            self._add("H2", self.elem_idx["H"], "ion", "H2")
        self.parent = np.array(self.parent)
        self.n = len(self.names)
        self.top = np.where(self.parent < 0)[0]
        self.children = {int(p): np.where(self.parent == p)[0] for p in self.top
                         if np.any(self.parent == p)}
        self.tva = np.array([k for k in range(self.n) if self.kind[k] in ("dust", "pah")], int)
        self.gas_top = np.array([k for k in self.top if self.kind[k] in ("elem", "CO")], int)
        self.ndust, self.npah = s["ndust"], s["npah"]
        self.h2, self.co = s["h2"], s["co"]

    def _add(self, name, parent, kind, rname):
        self.names.append(name); self.parent.append(parent); self.kind.append(kind)
        self.rnames.append(rname)


# ----------------------------------------------------------------------------------------------
# compositions (natural fractions) and physical mixtures
# ----------------------------------------------------------------------------------------------
def composition(lay, Z=1.0, dtg=0.0, ion="neutral", trace=1e-10, co=0.0, pah=0.0, seed=0):
    rng = np.random.default_rng(seed)
    f = np.zeros(lay.n)
    sol = {e[0]: e[3] for e in ELEMENTS}
    dust = dtg if lay.ndust else 0.0
    pahs = pah if lay.npah else 0.0
    co = co if lay.co else 0.0
    gas = 1.0 - dust - pahs - co
    xs = {e: sol[e] * (Z if e not in ("H", "He") else 1.0) for e in lay.elem_idx}
    tot = sum(xs.values())
    for e, k in lay.elem_idx.items():
        f[k] = gas * xs[e] / tot
    if lay.co:
        f[lay.names.index("CO")] = co
    for kind, total in (("dust", dust), ("pah", pahs)):
        ks = [k for k in lay.tva if lay.kind[k] == kind]
        if ks:
            w = rng.uniform(0.2, 1.0, len(ks))
            f[ks] = total * w / w.sum()
    for e, ks in lay.ion_idx.items():
        n = len(ks)
        if ion == "neutral":
            y = np.full(n, trace); y[0] = 1.0
        elif ion == "ionised":
            y = np.full(n, trace); y[min(1, n - 1)] = 1.0
        elif ion == "hot":
            c = (n - 1) * 0.6
            y = np.exp(-0.5 * ((np.arange(n) - c) / max(1.0, n / 6.0)) ** 2) + trace
        elif ion == "random":
            y = 10 ** rng.uniform(np.log10(trace), 0, n)
        else:
            raise ValueError(ion)
        f[ks] = y / y.sum()
    return f


def to_mass(lay, fnat):
    """Natural fractions -> physical mass fractions of every slot (ions as ion mass / rho)."""
    m = fnat.copy()
    for p, ks in lay.children.items():
        m[ks] = fnat[ks] * fnat[p]
    return m


def from_mass(lay, m):
    f = m.copy()
    for p, ks in lay.children.items():
        f[ks] = np.where(m[p] > 0, m[ks] / np.where(m[p] > 0, m[p], 1.0), 0.0)
    return f


def mix(lay, comps, weights):
    """Physical mixture of parcels with natural fractions comps[j] and MASS weights
    weights[j] (arrays over cells, summing to 1). Returns natural fractions (n, ncell)."""
    m = sum(np.outer(to_mass(lay, c), w) for c, w in zip(comps, weights))
    return from_mass(lay, m)


CONCENTRATIONS = {
    # three parcels A, B, C (natural fractions): metallicity, dust-to-gas, ionisation
    "solar": [dict(Z=1.0, dtg=1e-2, ion="neutral", co=1e-4, pah=1e-3, seed=1),
              dict(Z=1.0, dtg=5e-3, ion="ionised", co=1e-6, pah=1e-4, seed=2),
              dict(Z=1.0, dtg=1e-3, ion="hot", co=0.0, pah=1e-5, seed=3)],
    "extreme": [dict(Z=1e-4, dtg=1e-6, ion="neutral", trace=1e-10, co=1e-9, pah=1e-8, seed=4),
                dict(Z=3.0, dtg=0.1, ion="hot", trace=1e-10, co=1e-3, pah=1e-2, seed=5),
                dict(Z=0.1, dtg=1e-3, ion="random", trace=1e-10, co=0.0, pah=1e-5, seed=6)],
}


# ----------------------------------------------------------------------------------------------
# profiles: x in [0,1] (cell centres), returns rho, u, p, fnat (n, N)
# ----------------------------------------------------------------------------------------------
def _tanh(x, x0, w):
    return 0.5 * (1 + np.tanh((x - x0) / max(w, 1e-12)))


def profile(lay, name, x, conc="solar", rho_contrast=1.0, w=0.0, m=1, u0=1.0, **_):
    A, B, C = (composition(lay, **c) for c in CONCENTRATIONS[conc])
    N = len(x); dx = x[1] - x[0]
    rho = np.ones(N); u = np.full(N, u0); p = np.ones(N)
    one = np.ones(N)
    if name == "step2":          # two parcels only: collinear, independent limiting is exact
        wb = ((x > 0.25) & (x < 0.75)).astype(float)
        f = mix(lay, [A, B], [1 - wb, wb])
        rho = np.where(wb > 0, rho_contrast, 1.0)
    elif name == "step3":        # A | B (narrow, 2 cells) | C | A
        wb = ((x > 0.25) & (x < 0.25 + 2 * dx)).astype(float)
        wc = ((x >= 0.25 + 2 * dx) & (x < 0.75)).astype(float)
        f = mix(lay, [A, B, C], [1 - wb - wc, wb, wc])
        rho = 1.0 + (rho_contrast - 1.0) * wc
    elif name == "fronts":       # two smooth fronts of width w (cells), offset by 3 cells
        wb = _tanh(x, 0.3, w * dx) * (1 - _tanh(x, 0.7, w * dx))
        wc = _tanh(x, 0.3 + 3 * dx, w * dx) * (1 - _tanh(x, 0.7 - 3 * dx, w * dx))
        wc = 0.5 * wc; wb = 0.5 * wb
        f = mix(lay, [A, B, C], [1 - wb - wc, wb, wc])
    elif name == "smooth":       # periodic smooth mixture of A, B, C, wavenumber m
        wb = 0.4 * (0.5 + 0.5 * np.sin(2 * np.pi * m * x))
        wc = 0.4 * (0.5 + 0.5 * np.cos(2 * np.pi * m * x + 1.0))
        f = mix(lay, [A, B, C], [1 - wb - wc, wb, wc])
    elif name == "spike":        # smooth A/B background + a 1-cell C delta + a 3-cell B spike
        wb = 0.3 * (0.5 + 0.5 * np.sin(2 * np.pi * x))
        i0 = int(0.3 * N); i1 = int(0.6 * N)
        wc = np.zeros(N); wc[i0] = 1.0
        wb = np.where(wc > 0, 0.0, wb)
        wb[i1:i1 + 3] = 1.0
        f = mix(lay, [A, B, C], [1 - wb - wc, wb, wc])
    elif name == "sod3":         # Sod tube; A | B | C | A, contacts swept by the waves
        rho = np.where(x < 0.5, 1.0, 0.125); p = np.where(x < 0.5, 1.0, 0.1); u = 0 * x
        wb = ((x > 0.3) & (x < 0.5)).astype(float)
        wc = ((x >= 0.5) & (x < 0.62)).astype(float)
        f = mix(lay, [A, B, C], [1 - wb - wc, wb, wc])
    elif name == "uniform":      # one parcel at rest (source terms only); rho, p from rho_contrast, w
        f = mix(lay, [A], [one])
        rho = rho_contrast * one; p = w * one; u = 0 * x
    else:
        raise ValueError(name)
    return rho, u, p, f


# ----------------------------------------------------------------------------------------------
# test matrix
# ----------------------------------------------------------------------------------------------
SKIP_IN_NAME = ("test", "id", "ndim", "levelmin", "levelmax", "tend", "crossings", "bc", "u0",
                "calima_extra", "chem")


def case_name(c):
    return "_".join(f"{k}-{v}" for k, v in c.items() if k not in SKIP_IN_NAME)


ALL_SETS = ["H", "HHe_H2", "HHeCO", "full", "full_np0"]
CALIMA_SETS = ["HHeCO", "full", "full_np0"]


def _finish(cases):
    for c in cases:
        c.setdefault("ndim", 1)
        c["id"] = c["test"] + "_" + case_name(c)
    return cases


def advect_cases(quick=False):
    """advect: composition profile through a periodic box, uniform grid, 1 crossing."""
    profs = [("step2", {}), ("step3", {}), ("step3", {"rho_contrast": 100.0}),
             ("fronts", {"w": 0.5}), ("fronts", {"w": 4.0}),
             ("smooth", {"m": 1}), ("smooth", {"m": 8}), ("spike", {})]
    out = []
    for s, (p, kw), conc, st in itertools.product(ALL_SETS, profs, ["solar", "extreme"], [1, 2]):
        if quick and (st == 2 or conc == "solar" or p in ("fronts", "smooth") and kw.get("m", 0) != 8
                      and kw.get("w", 0) != 0.5):
            continue
        c = dict(test="advect", set=s, prof=p, conc=conc, slope=st)
        c.update(kw)
        c.update(levelmin=7, levelmax=7, tend=1.0, crossings=1)
        out.append(c)
    return _finish(out)


def sod_cases(quick=False):
    """sod: Sod tube with three compositions; the rarefaction, contact and shock sweep the
    composition contacts (outflow boundaries)."""
    out = []
    for s, conc, st in itertools.product(ALL_SETS, ["solar", "extreme"], [1, 2]):
        if quick and st == 2:
            continue
        out.append(dict(test="sod", set=s, prof="sod3", conc=conc, slope=st, levelmin=8, levelmax=8,
                        tend=0.2, bc="outflow", u0=0.0))
    return _finish(out)


def amr_cases(quick=False):
    """amr: as advect, on an AMR grid refined on the density and composition gradients, so
    that refined patches follow the contacts (prolongation of new cells, coarse-fine fluxes,
    restriction)."""
    out = []
    for s, (p, kw), conc in itertools.product(ALL_SETS, [("step3", {"rho_contrast": 100.0}), ("smooth", {"m": 4}),
                                                       ("spike", {})], ["solar", "extreme"]):
        if quick and conc == "solar":
            continue
        c = dict(test="amr", set=s, prof=p, conc=conc, slope=1, levelmin=6, levelmax=8, tend=1.0)
        c.update(kw)
        out.append(c)
    return _finish(out)


def tva_cases(quick=False):
    """tva: CALIMA dust and PAH bins drifting relative to the mixture (TVA with an imposed
    drift velocity, use_w_drift_test) through uniform gas with a smooth composition."""
    out = []
    for s, conc in itertools.product(CALIMA_SETS, ["solar", "extreme"]):
        if quick and conc == "solar":
            continue
        out.append(dict(test="tva", set=s, prof="smooth", m=1, conc=conc, slope=1, levelmin=7, levelmax=7,
                        tend=0.5, drift=0.3, u0=0.0))
    return _finish(out)


def source_cases(quick=False):
    """source: a uniform box at rest with the RTZ chemistry and the CALIMA dust processes on:
    cold dense gas (CO and H2 formation) and hot gas
    (sputtering returns metals). No transport: isolates the source terms."""
    out = []
    # (dust_accretion forces cmp_sigma_turb, which is 3D-only: the cold regime has the
    # molecule chemistry only in 1D)
    regimes = [("cold", dict(rho_contrast=1e3, w=3.3), ""),
               ("hot", dict(rho_contrast=1.0, w=139.0), "dust_sputtering=.true.\n")]
    for s, (reg, kw, extra) in itertools.product(CALIMA_SETS + ["HHe_H2"], regimes):
        if quick and s != "full":
            continue
        c = dict(test="source", set=s, prof="uniform", conc="solar", slope=1, regime=reg,
                 levelmin=3, levelmax=3, tend=20.0 if reg == "cold" else 2.0, chem=True)
        c.update(kw); c["calima_extra"] = extra
        out.append(c)
    return _finish(out)


TESTS = {"advect": advect_cases, "sod": sod_cases, "amr": amr_cases, "tva": tva_cases,
         "source": source_cases}


def write_ic_table(path, lay, x, rho, u, p, fnat):
    """ic_table.dat for patch/condinit.f90 (columns x rho vx P q_passive...)."""
    N = len(x)
    with open(path, "w") as o:
        o.write(f"{N} {4 + lay.n}\n")
        for i in range(N):
            row = [x[i], rho[i], u[i], p[i]] + list(fnat[:, i])
            o.write(" ".join(f"{v:.17e}" for v in row) + "\n")


def calima_block(lay, extra=""):
    if lay.ndust + lay.npah == 0:
        return ""
    s = SPECIES_SETS[lay.setname]
    if s.get("ndchemtype", 1) == 1:
        comp = "dust_composition(1,6) = 1d0\n"
        per = f"dustbins_per_chemtype={lay.ndust}\n"
        asz = ",".join(["0.01d0", "0.1d0", "0.005d0", "0.1d0"][:lay.ndust])
        sgr = ",".join(["2.2d0"] * lay.ndust)
        amn = ",".join(["4d-3", "3d-2", "1d-3", "2d-2"][:lay.ndust])
        amx = ",".join(["3d-2", "1d0", "2d-2", "1d0"][:lay.ndust])
    else:
        comp = ("dust_composition(1,6) = 1d0\ndust_composition(2,8)  = 4d0\n"
                "dust_composition(2,12) = 1d0\ndust_composition(2,14) = 1d0\n"
                "dust_composition(2,26) = 1d0\n")
        per = "dustbins_per_chemtype=2,2\n"
        asz, sgr = "0.01d0,0.1d0,0.005d0,0.1d0", "2.2d0,2.2d0,3.3d0,3.3d0"
        amn, amx = "4d-3,3d-2,1d-3,2d-2", "3d-2,1d0,2d-2,1d0"
    pah = ""
    if lay.npah:
        pah = ("pah_nc = 54,418\npah_nc_min = 12,101\npah_nc_max = 101,29952\n"
               "spah = 2d0,2d0\npah_ncharge_states = 4, 4\n")
    data = os.environ.get("RTZ_DATA", "/Users/currodri/Documents/RAMSES_dev/rtzdata_for_public_ramses/")
    return ("&CALIMA_PARAMS\n" + comp + per + f"asize={asz}\nsgrain={sgr}\namin={amn}\namax={amx}\n"
            + pah + "charging_model='WDB06isrf'\n" + extra
            + f"dust_tables_dir='{data}calima_data/'\n/\n")


def write_namelist(path, lay, template, **subs):
    """Fill a namelist template ({key} placeholders) and append the RT/CALIMA blocks."""
    data = os.environ.get("RTZ_DATA", "/Users/currodri/Documents/RAMSES_dev/rtzdata_for_public_ramses/")
    s = SPECIES_SETS[lay.setname]
    d = dict(data_dir=data, isHe=".true." if "He" in s["elements"] else ".false.",
             isH2=".true." if s["h2"] else ".false.", isCO=".true." if s["co"] else ".false.",
             calima="", calima_extra="")
    d.update(subs)
    txt = open(template).read()
    for k, v in d.items():
        txt = txt.replace("{" + k + "}", str(v))
    txt += calima_block(lay, d["calima_extra"])
    open(path, "w").write(txt)
