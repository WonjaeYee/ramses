import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'axes.linewidth': 0.8,
    'axes.facecolor': 'none',
    'figure.facecolor': 'none',
    'legend.frameon': False,
    'lines.linewidth': 2.5,
})

pi = 3.14159265

def parse_log(filename):
    params = {}
    bin_props = {}
    gas_elements = {}
    
    data_lines = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # Parse parameters from comments
                match_nh = re.match(r'# nH\s*=\s*([\d.E+-]+)', line)
                if match_nh:
                    params['nH'] = float(match_nh.group(1))
                match_tk = re.match(r'# Tk\s*=\s*([\d.E+-]+)', line)
                if match_tk:
                    params['Tk'] = float(match_tk.group(1))
                match_ne = re.match(r'# ne\s*=\s*([\d.E+-]+)', line)
                if match_ne:
                    params['ne'] = float(match_ne.group(1))
                match_mu = re.match(r'# (?:test_mu|local_mu)\s*=\s*([\d.E+-]+)', line)
                if match_mu:
                    params['test_mu'] = float(match_mu.group(1))
                match_ndust = re.match(r'# ndust\s*=\s*(\d+)', line)
                if match_ndust:
                    params['ndust'] = int(match_ndust.group(1))
                match_npah = re.match(r'# npah\s*=\s*(\d+)', line)
                if match_npah:
                    params['npah'] = int(match_npah.group(1))
                match_dt = re.match(r'# dt_yr\s*=\s*([\d.E+-]+)', line)
                if match_dt:
                    params['dt_yr'] = float(match_dt.group(1))
                match_acc = re.match(r'# dust_accretion\s*=\s*([TFtf])', line)
                if match_acc:
                    params['dust_accretion'] = (match_acc.group(1).lower() == 't')
                match_sput = re.match(r'# dust_sputtering\s*=\s*([TFtf])', line)
                if match_sput:
                    params['dust_sputtering'] = (match_sput.group(1).lower() == 't')
                match_coag = re.match(r'# dust_coagulation\s*=\s*([TFtf])', line)
                if match_coag:
                    params['dust_coagulation'] = (match_coag.group(1).lower() == 't')
                match_shat = re.match(r'# dust_shattering\s*=\s*([TFtf])', line)
                if match_shat:
                    params['dust_shattering'] = (match_shat.group(1).lower() == 't')
                match_tdir = re.match(r'# dust_tables_dir\s*=\s*(.*)', line)
                if match_tdir:
                    params['dust_tables_dir'] = match_tdir.group(1).strip()
                
                # Gas elements
                match_gelem = re.match(r'# GasElement\s+(\d+)\s+atomic_number\s*=\s*([-\d]+)\s+mass_g\s*=\s*([\d.E+-]+)', line)
                if match_gelem:
                    idx = int(match_gelem.group(1)) - 1
                    Z = int(match_gelem.group(2))
                    mass_g = float(match_gelem.group(3))
                    gas_elements[idx] = {'Z': Z, 'mass_g': mass_g}
                
                # Parse bin properties
                match_k0 = re.match(r'# Bin\s+(\d+)\s+k0_acc\s*=\s*([\d.E+-]+)\s+nhmax_acc\s*=\s*([\d.E+-]+)', line)
                if match_k0:
                    bin_idx = int(match_k0.group(1)) - 1
                    if bin_idx not in bin_props:
                        bin_props[bin_idx] = {'elements': []}
                    bin_props[bin_idx]['k0_acc'] = float(match_k0.group(2))
                    bin_props[bin_idx]['nhmax_acc'] = float(match_k0.group(3))
                    
                match_extra1 = re.match(
                    r'# Bin\s+(\d+)\s+asize\s*=\s*([\d.E+-]+)\s+sgrain\s*=\s*([\d.E+-]+)\s+mgrain\s*=\s*([\d.E+-]+)'
                    r'\s+amin\s*=\s*([\d.E+-]+)\s+amax\s*=\s*([\d.E+-]+)\s+surf_energy\s*=\s*([\d.E+-]+)'
                    r'\s+interact_group\s*=\s*(\d+)', line)
                if match_extra1:
                    bin_idx = int(match_extra1.group(1)) - 1
                    if bin_idx not in bin_props:
                        bin_props[bin_idx] = {'elements': []}
                    bin_props[bin_idx]['asize'] = float(match_extra1.group(2))
                    bin_props[bin_idx]['sgrain'] = float(match_extra1.group(3))
                    bin_props[bin_idx]['mgrain'] = float(match_extra1.group(4))
                    bin_props[bin_idx]['amin'] = float(match_extra1.group(5))
                    bin_props[bin_idx]['amax'] = float(match_extra1.group(6))
                    bin_props[bin_idx]['surf_energy'] = float(match_extra1.group(7))
                    bin_props[bin_idx]['interact_group'] = int(match_extra1.group(8))
                    
                match_extra2 = re.match(
                    r'# Bin\s+(\d+)\s+tensile_strength\s*=\s*([\d.E+-]+)\s+Youngs_modulus\s*=\s*([\d.E+-]+)\s+Poisson_ratio\s*=\s*([\d.E+-]+)', line)
                if match_extra2:
                    bin_idx = int(match_extra2.group(1)) - 1
                    if bin_idx not in bin_props:
                        bin_props[bin_idx] = {'elements': []}
                    bin_props[bin_idx]['tensile_strength'] = float(match_extra2.group(2))
                    bin_props[bin_idx]['Youngs_modulus'] = float(match_extra2.group(3))
                    bin_props[bin_idx]['Poisson_ratio'] = float(match_extra2.group(4))
                    
                match_elem = re.match(r'# Bin\s+(\d+)\s+element\s+(\d+)\s+index\s*=\s*(\d+)\s+mfraction\s*=\s*([\d.E+-]+)\s+mass_g\s*=\s*([\d.E+-]+)', line)
                if match_elem:
                    bin_idx = int(match_elem.group(1)) - 1
                    el_idx = int(match_elem.group(3))
                    mfraction = float(match_elem.group(4))
                    mass_g = float(match_elem.group(5))
                    if bin_idx not in bin_props:
                        bin_props[bin_idx] = {'elements': []}
                    bin_props[bin_idx]['elements'].append({
                        'index': el_idx,
                        'mfraction': mfraction,
                        'mass_g': mass_g
                    })
            else:
                # Data lines
                data_lines.append([float(x) for x in line.split()])
                
    params.setdefault('dust_accretion', True)
    params.setdefault('dust_sputtering', False)
    params.setdefault('dust_coagulation', False)
    params.setdefault('dust_shattering', False)
    params.setdefault('test_mu', 1.0)
    
    # Precompute mgrain_min and mgrain_max from amin, amax, sgrain
    for bin_idx, props in bin_props.items():
        if 'amin' in props and 'sgrain' in props:
            amin_cm = props['amin'] * 1.0e-4
            props['mgrain_min'] = (4.0/3.0) * pi * (amin_cm**3) * props['sgrain']
        if 'amax' in props and 'sgrain' in props:
            amax_cm = props['amax'] * 1.0e-4
            props['mgrain_max'] = (4.0/3.0) * pi * (amax_cm**3) * props['sgrain']
            
    data = np.array(data_lines)
    return params, bin_props, gas_elements, data

def sigmoid_function(k, x0, x):
    normx = x / x0
    return 1.0 / (1.0 + np.exp(-k * (normx - 1.0)))

def load_sputtering_table(filepath):
    with open(filepath, 'r') as f:
        # Skip 6 comment lines
        for _ in range(6):
            f.readline()
        # Read nT, nphi
        nT, nphi = map(int, f.readline().split())
        # Read phi_grid
        phi_grid = np.array([float(x) for x in f.readline().split()])
        # Find index closest to 0
        iphi0 = np.argmin(np.abs(phi_grid))
        
        # Read T_grid and rates
        T_grid = []
        rates = []
        for _ in range(nT):
            line_parts = [float(x) for x in f.readline().split()]
            T_grid.append(line_parts[0])
            rates.append(line_parts[1 + iphi0])
            
        return np.array(T_grid), np.array(rates)

def interpolate_sputtering_rate(lT, T_grid, rates):
    if lT <= T_grid[0]:
        return rates[0]
    elif lT >= T_grid[-1]:
        return rates[-1]
    
    idx = np.searchsorted(T_grid, lT) - 1
    t = (lT - T_grid[idx]) / (T_grid[idx+1] - T_grid[idx])
    return (1.0 - t) * rates[idx] + t * rates[idx+1]

def OC07_function(x):
    return 3.2 - (1.0 + x) + 2.0 / (1.0 + x) * (1.0 / 2.6 + (x**3) / (1.6 + x))

def grain_relative_velocity(model, T, rho_gas, nH, v_turb, local_mu, inject_L,
                            target_a, target_s, target_m,
                            projectile_a, projectile_s, projectile_m):
    kB = 1.3806490e-16  # erg/K
    mH = 1.6738233e-24  # g
    
    cs_gas = np.sqrt(5.0 / 3.0 * kB * T / (mH * local_mu))
    v_th = np.sqrt(8.0 / pi) * cs_gas
    
    if model == 'Ormel2007':
        dV_thermal = np.sqrt(8.0 * kB * T * (target_m + projectile_m) / (target_m * projectile_m))
        
        e2instatC = 2.3070775e-17
        rc = e2instatC / (kB * T)
        mfp = 1.0 / (nH * rc**2)
        tau_L = inject_L / v_turb
        Re = 3.0 * v_turb * inject_L / (cs_gas * mfp)
        tau_eta = tau_L / np.sqrt(Re)
        
        ts_target = target_s * target_a / (rho_gas * v_th)
        ts_projectile = projectile_s * projectile_a / (rho_gas * v_th)
        
        # Target must have ts_target >= ts_projectile
        if ts_target < ts_projectile:
            ts_target, ts_projectile = ts_projectile, ts_target
            target_m, projectile_m = projectile_m, target_m
            target_s, projectile_s = projectile_s, target_s
            target_a, projectile_a = projectile_a, target_a
            
        St_target = ts_target / tau_L
        St_projectile = ts_projectile / tau_L
        
        Stmin = tau_eta / tau_L
        if ts_target < tau_eta:
            if St_target + St_projectile > 0.0:
                dV_turb = np.sqrt(1.5) * v_turb * np.sqrt((St_target - St_projectile) / (St_target + St_projectile)) * \
                          np.sqrt((St_target**2 / (St_target + Stmin)) - (St_projectile**2 / (St_projectile + Stmin)))
            else:
                dV_turb = 0.0
        elif (tau_eta <= ts_target) and (ts_target < tau_L):
            dV_turb = np.sqrt(1.5) * v_turb * np.sqrt(OC07_function(St_projectile / St_target) * St_target)
        else:
            dV_turb = np.sqrt(1.5) * v_turb * np.sqrt(1.0 / (1.0 + St_target) + 1.0 / (1.0 + St_projectile))
            
        return np.sqrt(dV_thermal**2 + dV_turb**2)
    else:
        return v_turb

def compute_shattered_fragments_direct(ii, kk, v_rel, bin_props, npah, ndust, ii1, ii2, slope_frag_func=0.4333):
    mgrain1 = bin_props[ii]['mgrain']
    mgrain2 = bin_props[kk]['mgrain']
    
    # Calculate shear modulus and catastrophic specific energy
    # shear_modulus = Youngs_modulus / (2 * (1 + Poisson_ratio))
    # catastrophic_spec_energy = shear_modulus / (2 * sgrain)
    sm1 = bin_props[ii]['Youngs_modulus'] / (2.0 * (1.0 + bin_props[ii]['Poisson_ratio']))
    cat_egy1 = sm1 / (2.0 * bin_props[ii]['sgrain'])
    
    E_imp = 0.5 * (mgrain1 * mgrain2) / (mgrain1 + mgrain2) * v_rel**2
    phi = E_imp / (mgrain1 * cat_egy1)
    m_ej = phi / (1.0 + phi) * mgrain1
    
    m_remnant = mgrain1 - m_ej
    m_max = 2e-2 * m_ej
    m_min = 1e-6 * m_max
    
    prefactor = 0.0
    if m_ej > 1e-99:
        m_max_pow = m_max**slope_frag_func
        m_min_pow = m_min**slope_frag_func
        denom = m_max_pow - m_min_pow
        if abs(denom) > 1e-99:
            prefactor = m_ej / denom
            
    if ii == 2 and kk == 2:
        print(f"DEBUG BIN 3 self-shattering:")
        print(f"  v_rel = {v_rel:.15e}")
        print(f"  E_imp = {E_imp:.15e}")
        print(f"  phi = {phi:.15e}")
        print(f"  m_ej = {m_ej:.15e}")
        print(f"  m_remnant = {m_remnant:.15e}")
        print(f"  m_max = {m_max:.15e}")
        print(f"  m_min = {m_min:.15e}")
        print(f"  prefactor = {prefactor:.15e}")
            
    chi_frag_dest = 0.0
    chi_frag_dust = np.zeros(ndust)
    
    if prefactor > 0.0:
        if m_min < bin_props[ii1]['mgrain_min']:
            chi_frag_dest = prefactor * (min(bin_props[ii1]['mgrain_min'], m_max)**slope_frag_func - m_min**slope_frag_func)
            
        for ll in range(ii1, ii2 + 1):
            if m_min >= bin_props[ll]['mgrain_max'] or m_max < bin_props[ll]['mgrain_min']:
                chi_frag_dust[ll] = 0.0
            else:
                chi_frag_dust[ll] = prefactor * (min(bin_props[ll]['mgrain_max'], m_max)**slope_frag_func - \
                                                 max(bin_props[ll]['mgrain_min'], m_min)**slope_frag_func)
                
    remnant_assigned = False
    if m_remnant > 0.0:
        if ii1 <= ii <= ii2:
            if bin_props[ii]['mgrain_min'] <= m_remnant <= bin_props[ii]['mgrain_max']:
                chi_frag_dust[ii] += m_remnant
                remnant_assigned = True
                
        if not remnant_assigned:
            for ll in range(ii1, ii2 + 1):
                if bin_props[ll]['mgrain_min'] <= m_remnant < bin_props[ll]['mgrain_max']:
                    chi_frag_dust[ll] += m_remnant
                    remnant_assigned = True
                    break
                    
        if not remnant_assigned:
            nearest_idx = ii1
            min_logdist = abs(np.log(max(m_remnant, 1e-99) / bin_props[ii1]['mgrain']))
            for ll in range(ii1 + 1, ii2 + 1):
                logdist = abs(np.log(max(m_remnant, 1e-99) / bin_props[ll]['mgrain']))
                if logdist < min_logdist:
                    min_logdist = logdist
                    nearest_idx = ll
            chi_frag_dust[nearest_idx] += m_remnant
            remnant_assigned = True
            
    if ii == 2 and kk == 2:
        print(f"  BEFORE NORM: chi_frag_dest={chi_frag_dest:.15e}, chi_frag_dust={chi_frag_dust}")
    m_tot = chi_frag_dest + np.sum(chi_frag_dust)
    if m_tot > 1e-99:
        chi_frag_dest /= m_tot
        chi_frag_dust /= m_tot
    if ii == 2 and kk == 2:
        print(f"  AFTER NORM: chi_frag_dest={chi_frag_dest:.15e}, chi_frag_dust={chi_frag_dust}")
        
    return chi_frag_dest, chi_frag_dust

def integrate_rk54(t_eval, y0, n_elements, ndust, bin_props, gas_elements, params, prefactor, nH, errmax=0.1):
    """
    Adaptive Cash-Karp RK5(4) integration in Python matching the exact physics and strategy in RAMSES.
    """
    # Group bins by their constituent elements to identify chemical types
    chemtypes = {}
    for b_idx, props in bin_props.items():
        el_tuple = tuple(sorted([el['index'] for el in props['elements']]))
        if el_tuple not in chemtypes:
            chemtypes[el_tuple] = []
        chemtypes[el_tuple].append(b_idx)
        
    yr2sec = 3.15576000e7
    Myr2sec = 3.15576000e13
    
    # Cash-Karp Coefficients
    c2, c3, c4, c5, c6 = 0.2, 0.3, 0.6, 1.0, 0.875
    a21 = 0.2
    a31, a32 = 3.0/40.0, 9.0/40.0
    a41, a42, a43 = 0.3, -0.9, 1.2
    a51, a52, a53, a54 = -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0
    a61, a62, a63, a64, a65 = 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0
    
    b1 = 37.0/378.0
    b3 = 250.0/621.0
    b4 = 125.0/594.0
    b6 = 512.0/1771.0
    
    e1 = 37.0/378.0 - 2825.0/27648.0
    e3 = 250.0/621.0 - 18575.0/48384.0
    e4 = 125.0/594.0 - 13525.0/55296.0
    e5 = -277.0/14336.0
    e6 = 512.0/1771.0 - 0.25
    
    Tk = params['Tk']
    npah = params.get('npah', 0)
    
    # Pre-load sputtering tables
    sputtering_tables = {}
    if params.get('dust_sputtering'):
        tables_dir = params['dust_tables_dir']
        for ii in range(ndust):
            sputtering_tables[ii] = {}
            for iel_idx, elem_info in gas_elements.items():
                Zi = elem_info['Z']
                if Zi <= 0:
                    continue
                filename = f"sputtering_DustBin_{ii+1:02d}_Z_{Zi}"
                filepath = os.path.join(tables_dir, filename)
                if os.path.exists(filepath):
                    T_grid, rates = load_sputtering_table(filepath)
                    sputtering_tables[ii][iel_idx] = (T_grid, rates)
                    
    # Precompute coagulation coefficients (Aoyama2017)
    if params.get('dust_coagulation'):
        for ii in range(ndust):
            asize_cm = bin_props[ii]['asize'] * 1e-4
            mgrain = bin_props[ii]['mgrain']
            bin_props[ii]['k0_coa'] = 0.5 * 4.0 * pi * (asize_cm**2) * 1.0e7 / mgrain
            
    # Precompute relative velocities and shattering fragments
    cached_v_rel_dust_dust = np.zeros((ndust, ndust))
    cached_chi_frag_dest = np.zeros((ndust, ndust))
    cached_chi_frag_dust = {}
    
    if params.get('dust_shattering'):
        amu2g = 1.6605390e-24
        local_rho = nH * (1.4 * amu2g)
        local_sigma = 1.0e5
        local_dx = 1.0e18
        local_mu = params.get('test_mu', 1.0)
        
        groups = {}
        for idx in range(ndust):
            grp = bin_props[idx]['interact_group']
            if grp not in groups:
                groups[grp] = []
            groups[grp].append(idx)
            
        for grp, bin_indices in groups.items():
            ii1 = min(bin_indices)
            ii2 = max(bin_indices)
            for ii in bin_indices:
                for kk in bin_indices:
                    v_rel = grain_relative_velocity(
                        'Ormel2007', Tk, local_rho, nH, local_sigma, local_mu, local_dx,
                        bin_props[ii]['asize']*1e-4, bin_props[ii]['sgrain'], bin_props[ii]['mgrain'],
                        bin_props[kk]['asize']*1e-4, bin_props[kk]['sgrain'], bin_props[kk]['mgrain']
                    )
                    cached_v_rel_dust_dust[ii, kk] = v_rel
                    
                    dest, frag_dust = compute_shattered_fragments_direct(
                        ii, kk, v_rel, bin_props, npah, ndust, ii1, ii2,
                        slope_frag_func=params.get('slope_frag_func', 0.4333)
                    )
                    cached_chi_frag_dest[ii, kk] = dest
                    cached_chi_frag_dust[(ii, kk)] = frag_dust
        
    def derivs(y):
        dydt = np.zeros_like(y)
        y_gas = y[:n_elements]
        y_dust = y[n_elements:n_elements+ndust]
        kmax = 0.0
        
        # 1. Accretion rate computation
        if params.get('dust_accretion'):
            for el_tuple, bin_indices in chemtypes.items():
                first_bin = bin_indices[0]
                # Find the limiting element in this chemical type
                limit_rate = np.inf
                for el_info in bin_props[first_bin]['elements']:
                    e_idx = el_info['index'] - 1
                    f_e = el_info['mfraction']
                    m_e = el_info['mass_g']
                    pseudo_rate = y_gas[e_idx] / (f_e * np.sqrt(m_e))
                    if pseudo_rate < limit_rate:
                        limit_rate = pseudo_rate
                
                nhmax_acc = bin_props[first_bin]['nhmax_acc']
                tacc_max = 5.0
                log_nH = np.log10(max(nH, 1.0e-10))
                log_nhmax_acc = np.log10(nhmax_acc)
                sfunc = sigmoid_function(tacc_max, log_nhmax_acc, log_nH)
                        
                total_rate_type = 0.0
                for b_idx in bin_indices:
                    k0_acc = bin_props[b_idx]['k0_acc']
                    rate = limit_rate * k0_acc * prefactor
                    
                    # Apply accretion rate smoothing
                    if rate > 0.0 and sfunc > 0.0:
                        tacc_log = np.log10(1.0 / (rate * Myr2sec))
                        tacc_log = (1.0 - sfunc) * tacc_log + sfunc * tacc_max
                        rate = 1.0 / ((10.0**tacc_log) * Myr2sec)
                    
                    kmax = max(kmax, abs(rate))
                        
                    bin_rate = rate * y_dust[b_idx]
                    dydt[n_elements + b_idx] += bin_rate
                    total_rate_type += bin_rate
                    
                # Deplete gas elements
                for el_info in bin_props[first_bin]['elements']:
                    e_idx = el_info['index'] - 1
                    f_e = el_info['mfraction']
                    dydt[e_idx] -= total_rate_type * f_e

        # 2. Thermal Sputtering rate computation
        if params.get('dust_sputtering'):
            lT = np.log10(Tk)
            for ii in range(ndust):
                rate_total = 0.0
                for iel_idx, (T_grid, rates) in sputtering_tables[ii].items():
                    if y_gas[iel_idx] < 1e-40:
                        continue
                    irate = interpolate_sputtering_rate(lT, T_grid, rates)
                    n_gas = y_gas[iel_idx] / gas_elements[iel_idx]['mass_g']
                    rate_total += (10.0**irate) * n_gas
                    
                asize = bin_props[ii]['asize']
                rate1 = 3.0 * rate_total / asize / yr2sec
                kmax = max(kmax, abs(rate1))
                
                mass_loss = rate1 * y_dust[ii]
                dydt[n_elements + ii] -= mass_loss
                
                for el in bin_props[ii]['elements']:
                    e_idx = el['index'] - 1
                    mf = el['mfraction']
                    dydt[e_idx] += mass_loss * mf

        # 3. Coagulation rate computation (Aoyama2017)
        if params.get('dust_coagulation'):
            local_Jeans = 4.81973044e19 * np.sqrt(Tk / nH)
            local_dx = 1.0e18
            # Check conditions matching Aoyama2017_coagulation_rate
            if not (Tk > 1.0e4 or nH < 0.1 or local_Jeans > 4.0 * local_dx):
                groups = {}
                for ii in range(ndust):
                    grp = bin_props[ii]['interact_group']
                    if grp not in groups:
                        groups[grp] = []
                    groups[grp].append(ii)
                    
                for grp, bin_indices in groups.items():
                    bin_indices = sorted(bin_indices)
                    for ii in bin_indices[:-1]:
                        rate1 = bin_props[ii]['k0_coa'] * y_dust[ii] / nH
                        kmax = max(kmax, rate1)
                        rate2 = rate1 * y_dust[ii]
                        dydt[n_elements + ii] -= rate2
                        dydt[n_elements + ii + 1] += rate2

        # 4. Shattering rate computation (Hirashita2015)
        if params.get('dust_shattering'):
            groups = {}
            for idx in range(ndust):
                grp = bin_props[idx]['interact_group']
                if grp not in groups:
                    groups[grp] = []
                groups[grp].append(idx)
                
            for grp, bin_indices in groups.items():
                ii1 = min(bin_indices)
                ii2 = max(bin_indices)
                for jj in range(ii2, ii1 - 1, -1):
                    v_rel = cached_v_rel_dust_dust[jj, jj]
                    chi_frag_dest = cached_chi_frag_dest[jj, jj]
                    chi_frag = cached_chi_frag_dust[(jj, jj)]
                    
                    asize_cm = bin_props[jj]['asize'] * 1e-4
                    coll_factor = np.sqrt(8.0 / (3.0 * pi)) * 4.0 * pi * (asize_cm**2) * v_rel
                    
                    rate1 = coll_factor * y_dust[jj] / bin_props[jj]['mgrain']
                    kmax = max(kmax, abs(rate1))
                    
                    loss = rate1 * y_dust[jj]
                    dydt[n_elements + jj] -= loss
                    
                    for ll in range(ii1, ii2 + 1):
                        dydt[n_elements + ll] += rate1 * chi_frag[ll] * y_dust[jj]
                        
                    if chi_frag_dest > 1e-10:
                        rate_dest = rate1 * chi_frag_dest * y_dust[jj]
                        for el in bin_props[jj]['elements']:
                            e_idx = el['index'] - 1
                            mf = el['mfraction']
                            dydt[e_idx] += rate_dest * mf
                            
        return dydt, kmax

    # Integrate step-by-step
    y = np.copy(y0)
    ys = [np.copy(y0)]
    
    total_steps = 0
    total_rejected = 0
    
    # Integrate between evaluation points (test_dt intervals)
    for i in range(len(t_eval) - 1):
        t_start = t_eval[i] * yr2sec
        t_end = t_eval[i+1] * yr2sec
        dt_step = t_end - t_start
        
        tau = 0.0
        h = dt_step
        first_call = True
        
        while tau < dt_step:
            h = min(h, dt_step - tau)
            
            # Compute derivatives at current state to get k1 and kmax
            k1, kmax = derivs(y)
            
            if kmax == 0.0:
                break
                
            if first_call:
                h_local = min(1.0 / kmax, h)
            else:
                h_local = h
                
            # Stage 2
            y2 = y + h_local * a21 * k1
            y2[:n_elements] = np.maximum(y2[:n_elements], 0.0)
            k2, _ = derivs(y2)
            
            # Stage 3
            y3 = y + h_local * (a31 * k1 + a32 * k2)
            y3[:n_elements] = np.maximum(y3[:n_elements], 0.0)
            k3, _ = derivs(y3)
            
            # Stage 4
            y4 = y + h_local * (a41 * k1 + a42 * k2 + a43 * k3)
            y4[:n_elements] = np.maximum(y4[:n_elements], 0.0)
            k4, _ = derivs(y4)
            
            # Stage 5
            y5 = y + h_local * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
            y5[:n_elements] = np.maximum(y5[:n_elements], 0.0)
            k5, _ = derivs(y5)
            
            # Stage 6
            y6 = y + h_local * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
            y6[:n_elements] = np.maximum(y6[:n_elements], 0.0)
            k6, _ = derivs(y6)
            
            # Compute 5th-order solutions
            y_new = y + h_local * (b1 * k1 + b3 * k3 + b4 * k4 + b6 * k6)
            
            # Error estimation
            error_gas = np.abs(h_local * (e1 * k1[:n_elements] + e3 * k3[:n_elements] + e4 * k4[:n_elements] + e5 * k5[:n_elements] + e6 * k6[:n_elements])) / np.maximum(np.abs(y[:n_elements]), 1.0e-40)
            error_dust = np.abs(h_local * (e1 * k1[n_elements:] + e3 * k3[n_elements:] + e4 * k4[n_elements:] + e5 * k5[n_elements:] + e6 * k6[n_elements:])) / np.maximum(np.abs(y[n_elements:]), 1.0e-40)
            
            max_error = max(np.max(error_gas), np.max(error_dust))
            
            accepted = (max_error <= errmax)
            
            scale = 0.9 * (errmax / max(max_error, 1.0e-10))**0.2
            h_new = h_local * min(2.0, max(0.1, scale))
            
            if accepted:
                y = y_new
                y[:n_elements] = np.maximum(y[:n_elements], 0.0)
                y[n_elements:] = np.maximum(y[n_elements:], 0.0)
                tau += h_local
                first_call = False
                h = h_new
                total_steps += 1
            else:
                h = h_new
                total_rejected += 1
                
        ys.append(np.copy(y))
        
    print(f"Python RK54 stats: accepted steps = {total_steps}, rejected steps = {total_rejected}")
    return np.array(ys)

def main():
    if len(sys.argv) > 1:
        log_file = sys.argv[1]
    else:
        log_file = 'dust_accretion_test.log'
        
    if not os.path.exists(log_file):
        print(f"Error: {log_file} not found. Please specify a log file or run the compiled test binary first.")
        print("Usage: python plot_dust_test.py <log_file>")
        return
        
    print(f"Reading {log_file}...")
    params, bin_props, gas_elements, data = parse_log(log_file)
    
    t_yr = data[:, 0]
    ne = data[:, 1]
    
    ndust = params['ndust']
    npah = params['npah']
    n_elements = data.shape[1] - 2 - ndust - npah
    
    gas_idx_start = 2
    dust_idx_start = 2 + n_elements
    
    gas_data = data[:, gas_idx_start:dust_idx_start]
    dust_data = data[:, dust_idx_start:dust_idx_start+ndust]
    
    # Map back to old variables to preserve compatibility with existing plotting code
    bin_elements = {ii: props['elements'] for ii, props in bin_props.items()}
    bin_k0_acc = [bin_props[ii].get('k0_acc', 0.0) for ii in range(ndust)]
    bin_nhmax_acc = [bin_props[ii].get('nhmax_acc', 0.0) for ii in range(ndust)]
    
    el_atomic_masses_g = np.zeros(n_elements)
    for iel_idx, elem_info in gas_elements.items():
        el_atomic_masses_g[iel_idx] = elem_info['mass_g']
        
    # Initial state vector for python integrator (using mass densities)
    y0 = np.zeros(n_elements + ndust)
    y0[:n_elements] = gas_data[0, :] * el_atomic_masses_g
    y0[n_elements:] = dust_data[0, :]
    
    Tk = params['Tk']
    prefactor = np.sqrt(Tk) / (1.0 + 1.0e-4 * Tk**1.5)
    nH = params['nH']
    
    print("Performing reference Python integration...")
    ys_python = integrate_rk54(t_yr, y0, n_elements, ndust, bin_props, gas_elements, params, prefactor, nH)
    
    # Plot results
    fig, axes = plt.subplots(2, 1, figsize=(5.5, 6.5), sharex=True)
    
    # 1. Plot dust density evolution
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, ndust))
    y_analytic_final = None
    
    # Only do analytic check if it is pure accretion (no other processes)
    do_analytic = params.get('dust_accretion') and not (
        params.get('dust_sputtering') or 
        params.get('dust_coagulation') or 
        params.get('dust_shattering')
    )
    
    for i in range(ndust):
        # RAMSES simulation data
        axes[0].plot(t_yr, dust_data[:, i] / dust_data[0, i], 'o', color=colors[i], label=f'Ramses Bin {i+1}')
        # Python integration data
        if ys_python is not None:
            axes[0].plot(t_yr, ys_python[:, n_elements + i] / y0[n_elements + i], '--', color=colors[i], alpha=0.8, label=f'Python Bin {i+1}')
        
            # Analytic solution for Bin 1 (Carbonaceous dust growth)
            if do_analytic and i == 0 and 0 in bin_elements:
                elems = bin_elements[0]
                if len(elems) == 1 and elems[0]['index'] == 6:
                    f_c = elems[0]['mfraction']
                    m_c = elems[0]['mass_g']
                    y_dust_0 = y0[n_elements]
                    y_gas_c_0 = y0[5] # Carbon is index 6 (0-based 5)
                    y_tot = y_gas_c_0 / f_c + y_dust_0
                    
                    # Check sfunc smoothing for analytical solution at t=0
                    nhmax_acc = bin_nhmax_acc[0]
                    tacc_max = 5.0
                    sfunc = sigmoid_function(tacc_max, np.log10(nhmax_acc), np.log10(max(nH, 1.0e-10)))
                    
                    rate_0 = (y_gas_c_0 / (f_c * np.sqrt(m_c))) * bin_k0_acc[0] * prefactor
                    if rate_0 > 0.0 and sfunc > 0.0:
                        Myr2sec = 3.15576000e13
                        tacc_log = np.log10(1.0 / (rate_0 * Myr2sec))
                        tacc_log = (1.0 - sfunc) * tacc_log + sfunc * tacc_max
                        rate_0_smooth = 1.0 / ((10.0**tacc_log) * Myr2sec)
                    else:
                        rate_0_smooth = rate_0
                        
                    C_coeff = rate_0_smooth / y_gas_c_0
                    
                    # Logistic equation using initial smoothed rate
                    yr2sec = 3.15576000e7
                    t_sec = t_yr * yr2sec
                    y_analytic = y_tot / (1.0 + ((y_tot - y_dust_0) / y_dust_0) * np.exp(-C_coeff * y_tot * t_sec))
                    
                    axes[0].plot(t_yr, y_analytic / y_dust_0, ':', color='#b71c1c', linewidth=2.5, label='Analytic Solution (Logistic)')
                    y_analytic_final = y_analytic[-1]
                
    axes[0].set_ylabel(r'Dust Mass Density ($\rho(t)/\rho(0)$)')
    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[0].legend(loc='lower right',frameon=False)
    axes[0].grid(True, linestyle=':', alpha=0.3, color='#aaaaaa')
    
    # 2. Plot gas metal abundance evolution
    # Dynamically find all elements present in the dust compositions
    track_elements = set()
    for b_idx, props in bin_props.items():
        for el in props['elements']:
            track_elements.add(el['index'])
    track_elements = sorted(list(track_elements))
    
    element_symbols = {
        1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 12: 'Mg', 14: 'Si', 16: 'S', 26: 'Fe'
    }
    element_colors = {
        6: '#212121',
        8: '#b71c1c',
        12: '#1b5e20',
        14: '#1565c0',
        26: '#e65100',
        1: '#616161',
        2: '#4a148c',
        7: '#006064',
        10: '#880e4f',
        16: '#827717',
    }
    
    for idx, el_idx in enumerate(track_elements):
        if el_idx <= n_elements:
            e_idx = el_idx - 1
            symbol = element_symbols.get(el_idx, f'El {el_idx}')
            color = element_colors.get(el_idx, plt.cm.tab10(idx % 10))
            # Number density from RAMSES
            axes[1].plot(t_yr, gas_data[:, e_idx] / gas_data[0, e_idx], '.', color=color, label=f'Ramses {symbol}')
            # Number density from Python
            if ys_python is not None:
                gas_num_py = ys_python[:, e_idx] / el_atomic_masses_g[e_idx]
                axes[1].plot(t_yr, gas_num_py / gas_num_py[0], '-', color=color, alpha=0.8, label=f'Python {symbol}')
            
    axes[1].set_xlabel('Time (yr)')
    axes[1].set_xscale('log')
    axes[1].set_ylabel('Gas-Phase Abundance ($n(t)/n(0)$)')
    axes[1].legend(loc='lower left',frameon=False,ncol=2)
    axes[1].grid(True, linestyle=':', alpha=0.3, color='#aaaaaa')
    
    plt.tight_layout()
    
    suffix = 'accretion'
    if params.get('dust_sputtering'):
        suffix = 'sputtering'
    elif params.get('dust_coagulation'):
        suffix = 'coagulation'
    elif params.get('dust_shattering'):
        suffix = 'shattering'
    plot_file = f'dust_{suffix}_test.png'
    
    plt.savefig(plot_file, dpi=300, transparent=True)
    print(f"Plot saved to {plot_file}")
    
    # Print comparison report
    print("\n" + "="*80)
    print("                      CALIMA DUST SOLVER TEST REPORT")
    print("="*80)
    print(f"Final simulation time: {t_yr[-1]:.3e} yr")
    print(f"Gas density (nH): {nH:.3e} cm-3, Temperature (Tk): {Tk:.3e} K, ne: {ne[0]:.3e} cm-3, mu: {params.get('test_mu'):.3f}")
    
    active_procs = []
    if params.get('dust_accretion'): active_procs.append('accretion')
    if params.get('dust_sputtering'): active_procs.append('sputtering')
    if params.get('dust_coagulation'): active_procs.append('coagulation')
    if params.get('dust_shattering'): active_procs.append('shattering')
    print(f"Active processes: {', '.join(active_procs)}")
    
    print("-"*80)
    if ys_python is not None:
        print(f"{'Species/Bin':<15} | {'Ramses Final (norm)':<20} | {'Python Final (norm)':<20} | {'Rel. Difference':<15}")
        print("-"*80)
        for i in range(ndust):
            ramses_val = dust_data[-1, i] / dust_data[0, i]
            python_val = ys_python[-1, n_elements + i] / y0[n_elements + i]
            rel_diff = abs(ramses_val - python_val) / max(ramses_val, 1e-40)
            print(f"Dust Bin {i+1:<6} | {ramses_val:<20.6e} | {python_val:<20.6e} | {rel_diff:<15.6e}")
            
        if y_analytic_final is not None:
            y_dust_0 = y0[n_elements]
            analytic_norm = y_analytic_final / y_dust_0
            ramses_norm = dust_data[-1, 0] / dust_data[0, 0]
            rel_diff_analytic = abs(ramses_norm - analytic_norm) / max(ramses_norm, 1e-40)
            print("-"*80)
            print(f"Analytic Bin 1  | {analytic_norm:<20.6e} | (relative to Ramses)       | {rel_diff_analytic:<15.6e}")
    print("="*80 + "\n")
    
if __name__ == '__main__':
    main()
