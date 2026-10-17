#!/usr/bin/env python3
"""One-zone tests of CALIMA grain charging / PE heating on the local RT groups.

Each run is a periodic 3D box of 2^3 identical cells (levelmin = levelmax = 1)
at fixed density and temperature (rt_TConst), lit by one 'square' source per
photon group covering the whole box. A square source sets Np = rt_n_source /
c_red in every cell each step, so rt_n_source(g) is the photon flux of group g
[photons cm^-2 s^-1]; the fluxes follow a chosen SED normalised to G0 Habing in
6-13.6 eV. Exactly one dust bin has a non-zero abundance per run
(NDUST = 4, NPAH = 2, PAHs at zero).

    python make_runs.py --outdir RUNS --exe-new ../../bin_onezone/ramses3d \
        --exe-ref ../../../staging-ref/bin_onezone/ramses3d \
        --rtgroups-tables /path/to/pyCALIMA/model_data/dust_charging_rtgroups_data

writes RUNS/<model>/<source>_bin<k>/{namelist.nml,run.sh} for model = new
(charging_model='WDB06rt', debug dump on) and ref (charging_model='WDB06isrf':
the uniform-ISRF tables at the G0 of the groups; --exe-ref needs a build that has
it), plus a tables directory that links the standard CALIMA tables and the
RT-group ones. Compare with pycalima diagnostics/dust_charge/compare_ramses_rtgroups.py.
"""
import argparse
import os
from pathlib import Path

import numpy as np

C = 2.99792458e10
EV = 1.602176634e-12
KB_EV = 8.617333262e-5
U_HAB = 5.33e-14
M_H = 1.6735575e-24
UNITS_D, UNITS_T, UNITS_L = 1.66e-24, 3.1557e13, 3.0857e18

GROUP_L0 = [0.1, 1.0, 5.6, 11.2, 13.6, 15.2, 24.59, 54.42]
GROUP_L1 = [1.0, 5.6, 11.2, 13.6, 15.2, 24.59, 54.42, 500.0]

# mass fractions, as the CALIMA slab runs
X_H, Y_HE = 0.73769344, 2.49279831e-01
METALS = [X_H, Y_HE, 1.0671119504e-03, 6.9294363934e-04, 4.3716649936e-03, 1.2567811419e-03,
          1.8925730137e-04, 6.6595143634e-05, 3.0975659757e-04, 9.9721662758e-05]
N_IONS = [2, 3, 7, 6, 6, 6, 6, 6, 6, 6]           # H He C N O Ne Mg Si S Fe
DUST_MASS_FRACTION = 1e-3


def sed(name, E):
    if name.startswith('BB'):
        T = float(name[2:])
        return E ** 3 / np.expm1(np.minimum(E / (KB_EV * T), 700.0))
    if name == 'quasar':  # Sazonov et al. (2004), WDB06 eq. 29
        return np.where(E < 10, 0.0798 * E ** -0.6,
                        np.where(E < 2e3, E ** -1.7 * np.exp(E / 2e3),
                                 2.94e-3 * E ** -0.8 * np.exp(-E / 2e5)))
    raise ValueError(name)


def group_fluxes(name, G0):
    """Photon flux per group [cm^-2 s^-1] of SED `name` at G0 Habing (6-13.6 eV)."""
    E = np.geomspace(0.1, 1000.0, 20000)
    uE = sed(name, E)
    b = (E >= 6) & (E <= 13.6)
    uE = uE * G0 * U_HAB / np.trapezoid(uE[b], E[b])
    F = []
    for l0, l1 in zip(GROUP_L0, GROUP_L1):
        m = (E >= l0) & (E < l1)
        F.append(C * np.trapezoid(uE[m] / (E[m] * EV), E[m]))
    return np.array(F)


def namelist(source, G0, dust_bin, nH, T, model, tables_dir, data_dir, nstep, free_T=False):
    F = group_fluxes(source, G0)
    d = nH * M_H / (X_H * UNITS_D)
    p = 2.3 * nH * 1.380649e-16 * T / (UNITS_D * (UNITS_L / UNITS_T) ** 2)
    L = []
    L.append(f"! CALIMA one-zone RT-group charging test: {source} at G0 = {G0:g}, "
             f"n_H = {nH:g}, T = {T:g} K, dust bin {dust_bin} only, model = {model}")
    L.append("&RUN_PARAMS\nhydro=.true.\nrt=.true.\npoisson=.false.\npic=.false.\nnrestart=0\n"
             f"ncontrol=1\nnstepmax={nstep}\nnsubcycle=10*1\nverbose=.false.\n/")
    L.append("&AMR_PARAMS\nlevelmin=1\nlevelmax=1\nngridmax=100\nnexpand=1\nboxlen=1.0\n/")
    reg = ["&INIT_PARAMS", "nregion=1", "region_type(1)='square'", "x_center(1)=0.5", "y_center(1)=0.5",
           "z_center(1)=0.5", "length_x(1)=10", "length_y(1)=10", "length_z(1)=10", "exp_region(1)=10.0",
           f"d_region(1)={d:.10e}", f"p_region(1)={p:.10e}"]
    iv = 1
    for m in METALS:
        reg.append(f"var_region(1,{iv})={m:.10e}")
        iv += 1
    reg.append(f"var_region(1,{iv})=0.0          ! CO")
    iv += 1
    for k in range(2):
        reg.append(f"var_region(1,{iv})=0.0          ! PAHBin_0{k+1}")
        iv += 1
    for k in range(1, 5):
        reg.append(f"var_region(1,{iv})={DUST_MASS_FRACTION if k == dust_bin else 0.0:.10e}   ! DustBin_0{k}")
        iv += 1
    for n in N_IONS:
        for s in range(n):
            reg.append(f"var_region(1,{iv})={1.0 if s == 0 else 1e-6:.6e}")
            iv += 1
    reg.append(f"var_region(1,{iv})=0.0          ! H2")
    reg.append("/")
    L.append("\n".join(reg))
    L.append("&OUTPUT_PARAMS\nnoutput=1\ntout=1d10\n/")
    L.append("&HYDRO_PARAMS\ngamma=1.666667\ncourant_factor=0.8\nscheme='muscl'\nslope_type=1\n"
             f"riemann='hllc'\npressure_fix=.true.\nbeta_fix=0.5\ndata_dir='{data_dir}'\n/")
    L.append("&COOLING_PARAMS\ncooling=.true.\nmetal=.true.\nneq_chem=.true.\nz_ave=1\n/")
    L.append(f"&UNITS_PARAMS\nunits_density={UNITS_D}\nunits_time={UNITS_T}\nunits_length={UNITS_L}\n/")
    rt = ["&RT_PARAMS", f"X={X_H}", f"Y={Y_HE}", "rt_flux_scheme='glf'", "rt_courant_factor=0.8",
          "rt_c_fraction=0.001", "rt_smooth=.true.", "rt_nsubcycle=10", "rt_otsa=.true.",
          "rt_is_init_xion=.false.", "rt_output_coolstats=.true.", "rtz_cooling=.true.",
          f"rt_TConst={0.0 if free_T else T}", "isHe=.true.", "isH2_rtz=.false.", "isCO_rtz=.false.",
          "rtz_include_collisional_ionization=.true.", "rtz_include_photoionization=.true.",
          "rtz_include_charge_exchange=.true.", "rtz_include_dust_recombination=.true.",
          "rtz_include_cosmic_ray_ionization=.false.", "rtz_primary_cosmic_ray_ionization_rate=0.d0",
          "rtz_include_HM12_UVB=.false.", "rtz_UV_background_G0=0.d0",
          "rt_nsource=8", "rt_source_type=8*'square'",
          "rt_src_x_center=8*0.5", "rt_src_y_center=8*0.5", "rt_src_z_center=8*0.5",
          "rt_src_length_x=8*10.0", "rt_src_length_y=8*10.0", "rt_src_length_z=8*10.0",
          "rt_exp_source=8*10.0", "rt_src_group=1,2,3,4,5,6,7,8", "rt_u_source=8*0.0",
          "rt_n_source=" + ",".join(f"{f:.6e}" for f in F), "/"]
    L.append("\n".join(rt))
    L.append("&RT_GROUPS\ngroupL0 = " + ", ".join(f"{x:g}" for x in GROUP_L0)
             + "\ngroupL1 = " + ", ".join(f"{x:g}" for x in GROUP_L1) + "\n/")
    cal = ["&CALIMA_PARAMS", "dust_pe_heating=.true.", "dust_pe_heating_isrf=.false.",
           "pah_pe_heating=.true.",
           f"charging_model='{'WDB06rt' if model == 'new' else 'WDB06isrf'}'",
           "dust_composition(1,6) = 1d0", "dust_composition(2,8)  = 4d0", "dust_composition(2,12) = 1d0",
           "dust_composition(2,14) = 1d0", "dust_composition(2,26) = 1d0", "dustbins_per_chemtype=2,2",
           "asize=0.01d0,0.1d0,0.005d0,0.1d0", "sgrain=2.2d0,2.2d0,3.3d0,3.3d0",
           "amin = 4d-3,3d-2,1d-3,2d-2", "amax = 3d-2,1d0,2d-2,1d0",
           "pah_nc = 54,418", "pah_nc_min = 12,101", "pah_nc_max = 101,29952", "spah = 2d0,2d0",
           "pah_ncharge_states = 4, 4", f"dust_tables_dir='{tables_dir}/'"]
    if model == 'new':
        cal += ["dust_rtgroups_debug=.true.", "dust_rtgroups_debug_max=4000"]
    cal.append("/")
    L.append("\n".join(cal))
    return "\n\n".join(L) + "\n", F


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--exe-new', required=True)
    ap.add_argument('--exe-ref', required=True)
    ap.add_argument('--rtgroups-tables', required=True,
                    help='directory with dust_charging_rtgroups_DustBin_XX.dat (pyCALIMA export)')
    ap.add_argument('--calima-tables',
                    default='/Users/currodri/Documents/RAMSES_dev/rtzdata_for_public_ramses/calima_data')
    ap.add_argument('--data-dir', default='/Users/currodri/Documents/RAMSES_dev/rtzdata_for_public_ramses/')
    ap.add_argument('--sources', nargs='+', default=['BB4e4:1e3', 'BB1e5:1e4', 'quasar:1e3', 'BB1e4:10'],
                    help='SED:G0 pairs (SED = BB<T> or quasar)')
    ap.add_argument('--nH', type=float, default=100.0)
    ap.add_argument('--T', type=float, default=8000.0)
    ap.add_argument('--nstep', type=int, default=20)
    ap.add_argument('--free-T', action='store_true', help='let T evolve from --T (rt_TConst = 0)')
    args = ap.parse_args()

    out = Path(args.outdir).resolve()
    tables = out / 'tables'
    tables.mkdir(parents=True, exist_ok=True)
    for src in list(Path(args.calima_tables).iterdir()) + sorted(Path(args.rtgroups_tables).glob('*.dat')):
        dst = tables / src.name
        if not dst.exists():
            os.symlink(src.resolve(), dst)
    for model, exe in (('new', args.exe_new), ('ref', args.exe_ref)):
        for spec in args.sources:
            source, G0 = spec.split(':')
            for k in range(1, 5):
                run = out / model / f'{source}_G{G0}_bin{k}{"_freeT" if args.free_T else ""}'
                run.mkdir(parents=True, exist_ok=True)
                # relative: CALIMA builds some table file names in character(len=128)
                text, F = namelist(source, float(G0), k, args.nH, args.T, model, '../../tables', args.data_dir,
                                   args.nstep, args.free_T)
                (run / 'namelist.nml').write_text(text)
                (run / 'SEDtables').mkdir(exist_ok=True)   # RAMSES writes its per-bin SED tables here
                (run / 'run.sh').write_text(f"#!/bin/sh\ncd {run}\n{Path(exe).resolve()} namelist.nml > run.log 2>&1\n")
                os.chmod(run / 'run.sh', 0o755)
    print(f'wrote runs under {out}')


if __name__ == '__main__':
    main()
