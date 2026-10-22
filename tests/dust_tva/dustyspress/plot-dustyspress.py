import sys
import os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
mpl.rcParams.update({
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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

def run_simulation():
    """Run the simulation once using the parameter file."""
    print("--- Running DustySpress simulation ---")
    subprocess.run(
        ["./ramses_dust_test1d", "dustyspress.nml"],
        check=True
    )

# Run the simulation
run_simulation()

# =====================================================================
# 1. PHYSICAL CONSTANTS (CGS) & ENVIRONMENTAL ASSUMPTIONS
# =====================================================================
pc_to_cm    = 3.08567758149137e18  # 1 pc in cm
myr_to_sec  = 3.1556926e13         # 1 Myr in s
c_cgs       = 2.99792458e10        # Speed of light (cm/s)

# Parameters tailored to hit the ~3 Myr timeline objective:
rho_0       = 1.0e-24              # Background mixture density (g/cm3)
c_s         = 2.0e4                # Sound speed (2e4 cm/s = 0.2 km/s)
# Note: we set units_density = 1e-24, d_region = 1.0 in nml, so physical rho_0 = 1e-24 g/cm3.
# We also set p_region = 0.02987 in nml, so physical c_s = 2e4 cm/s.
c_fraction  = 0.01
c_sim       = c_fraction * c_cgs   # Reduced speed of light in simulation
group_egy   = 8.0                  # Photon group energy (eV)
eV2erg      = 1.60217663e-12       # eV to erg conversion
E_photon    = group_egy * eV2erg   # Energy per photon (erg)
F_photon    = 1e9                  # Fixed photon flux (photons/cm2/s)
F_0         = F_photon * E_photon  # Injected energy flux (erg/s/cm2)

# =====================================================================
# 2. GRAIN PROPERTIES & KINEMATICS (0.1 micron Graphite Grain)
# =====================================================================
a_1         = 1.0e-5               # Grain radius (0.1 micron = 1e-5 cm)
rho_s       = 2.2                  # Material density of graphite (g/cm3)
sigma_pr_1  = 1.5577448922252549e-12  # Given radiation cross section (cm2)

# Mass of a single spherical grain:
m_grain     = (4.0 / 3.0) * np.pi * (a_1**3) * rho_s

# Mass-specific dust opacity:
kappa_dust  = sigma_pr_1 / m_grain  # ~169.03 cm2/g

# Epstein aerodynamic drag stopping time (t_s):
t_s         = (rho_s * a_1) / (rho_0 * c_s)  # ~1.1e15 seconds (~34.8 Myr)

# Exact terminal drift velocity (u_drift) from the TVA approximation loop:
u_drift_pc_myr = 0.77692  # For use_w_drift_test = .true. with w_drift_test = 0.77692
u_drift_cgs = (u_drift_pc_myr * pc_to_cm) / myr_to_sec

print(f"--- KINEMATIC REPORT FOR BENCHMARK ---")
print(f"Grain mass (m_grain):      {m_grain:.5e} g")
print(f"Specific Opacity (kappa):  {kappa_dust:.2f} cm2/g")
print(f"Stopping time (t_s):       {t_s / myr_to_sec:.2f} Myr")
print(f"Terminal drift velocity:   {u_drift_cgs * 1e-5:.5f} km/s ({u_drift_pc_myr:.5f} pc/Myr)\n")

# =====================================================================
# 3. INITIAL CONDITIONS & SIMULATION REGIME
# =====================================================================
boxlen      = 10.0                 # 10 pc domain
x0          = 2.0                  # Initial shell center (pc)
sigma_shell = 0.4                  # Initial characteristic width (pc)
eps_base    = 1.0e-8               # Ambient dust background
delta_eps   = 1.0e-4                 # Shell injection amplitude

# Setting the 4 desired snapshot targets matching our tend calculation
output_times = np.array([1.0,2.0,4.0,6.0])  # in Myr
x_eval = np.linspace(0.0, boxlen, 300)

# Initial Gaussian Profile function
def initial_profile(x):
    return eps_base + delta_eps * np.exp(-((x - x0)**2) / (2.0 * sigma_shell**2))

# =====================================================================
# 4. PLOT OVERLAY
# =====================================================================
fig, ax = plt.subplots(figsize=(5.5, 4))

# Plot analytical initial condition
ax.plot(x_eval, initial_profile(x_eval), 'k-', alpha=0.4, linewidth=2.5, label='Initial (Analytic)')

# Load and plot numerical initial condition (Snapshot 1)
try:
    snap_init = visu_ramses.load_snapshot(1)
    order = snap_init["data"]["x"].argsort()
    x_num = snap_init["data"]["x"][order]
    f_dust_num = snap_init["data"]["DustBin_01"][order]
    ax.plot(x_num, f_dust_num, 'k--', linewidth=2.5, label='Initial (Numerical)')
except Exception as e:
    print(f"Error loading initial snapshot: {e}")

# Loop and calculate spatial shifts for each specified output time
colors = ['#1565c0', '#b71c1c', '#1b5e20', '#6a1b9a']
for i, t in enumerate(output_times):
    # Pure hyperbolic advection shift: x_shifted = x0 + u_drift * t
    x_shifted = x0 + (u_drift_pc_myr * t)
    
    # Evaluate shifted analytical distribution profile
    eps_analytic = eps_base + delta_eps * np.exp(-((x_eval - x_shifted)**2) / (2.0 * sigma_shell**2))
    
    ax.plot(x_eval, eps_analytic, '-', color=colors[i], linewidth=2.5, alpha=0.5,
             label=f't = {t:.2f} Myr (Analytic)')
             
    # Load and plot matching numerical snapshot (Snapshot 2 is t=0.74, 3 is t=1.48, etc.)
    snap_idx = i + 2
    try:
        snap = visu_ramses.load_snapshot(snap_idx)
        order = snap["data"]["x"].argsort()
        x_num = snap["data"]["x"][order]
        f_dust_num = snap["data"]["DustBin_01"][order]
        ax.plot(x_num, f_dust_num, '--', color=colors[i], linewidth=2.5,
                 label=f't = {t:.2f} Myr (Numerical)')
    except Exception as e:
        print(f"Error loading snapshot {snap_idx}: {e}")

# Grid presentation properties
ax.set_xlim(0.0, boxlen)
ax.set_ylim(0.0, delta_eps * 1.2)
ax.set_xlabel('Position x [pc]', fontsize=12)
ax.set_ylabel('Dust Mass Fraction $\epsilon_1$', fontsize=12)
ax.legend(loc='upper right', frameon=False)
ax.grid(True, linestyle=':', alpha=0.3, color='#aaaaaa')
fig.tight_layout()

# Save figure
fig.savefig("dustyspress.png", dpi=300, bbox_inches="tight", format='png', transparent=True)
print("Saved convergence plot to dustyspress.png")