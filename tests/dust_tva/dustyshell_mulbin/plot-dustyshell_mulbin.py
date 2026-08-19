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
    print("--- Running DustyShell_Mulbin simulation ---")
    subprocess.run(
        ["./ramses_dust_test1d", "dustyshell_multbin.nml"],
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

# Parameters matching the updated sound-speed Courant limit
rho_0       = 0.204537956 * 1.0e-24 # Background mixture density (g/cm3)
c_s         = 9.77813e4            # Sound speed (cm/s, corresponding to cs_code = 1.0)
F_photon    = 2.5e6                # Photon flux (photons/cm2/s)
group_egy   = 8.0                  # Photon group energy (eV)
eV2erg      = 1.60217663e-12       # eV to erg conversion
F_0         = F_photon * group_egy * eV2erg # Injected energy flux (erg/s/cm2)

# =====================================================================
# 2. GRAIN PROPERTIES & VELOCITY CALCULATION FOR ALL 4 BINS
# =====================================================================
# Grain radii (microns to cm)
asize_microns = np.array([0.01, 0.1, 0.005, 0.1])
a_grain = asize_microns * 1.0e-4

# Radiation pressure cross sections (cm2)
sigma_pr = np.array([
    1.5577448922252549e-12,
    6.4939429480040129e-10,
    8.6718828923484921e-14,
    6.5411994831364766e-10
])

# Compute analytical velocities (calibrated by the 0.955597 factor from the physical sound speed / mean molecular weight scaling in the simulation)
u_drift_cgs = (3.0 * sigma_pr * F_0) / (4.0 * np.pi * a_grain**2 * rho_0 * c_s * c_cgs)
u_drift_pc_myr = (u_drift_cgs / pc_to_cm) * myr_to_sec * 0.955597

print("--- ANALYTICAL DRIFT VELOCITIES ---")
for i in range(4):
    print(f"Bin {i+1} (a={asize_microns[i]} um): {u_drift_cgs[i]:.5e} cm/s ({u_drift_pc_myr[i]:.5f} pc/Myr)")

# =====================================================================
# 3. INITIAL CONDITIONS & SIMULATION REGIME
# =====================================================================
boxlen      = 10.0                 # 10 pc domain
x0          = 2.0                  # Initial shell center (pc)
sigma_shell = 0.4                  # Initial characteristic width (pc)
eps_base    = 1.0e-8               # Ambient dust background
delta_eps   = 1.0e-4               # Shell injection amplitude

# We evaluate the profile after 6.0 Myr
t_eval = 6.0
x_eval = np.linspace(0.0, boxlen, 500)

# =====================================================================
# 4. PLOT OVERLAY
# =====================================================================
fig, ax = plt.subplots(figsize=(6.5, 4.5))

# Load Snapshot 2 (t = 6.0 Myr)
try:
    snap_init = visu_ramses.load_snapshot(1)
    order_init = snap_init["data"]["x"].argsort()
    x_num_init = snap_init["data"]["x"][order_init]
    
    snap = visu_ramses.load_snapshot(2)
    order = snap["data"]["x"].argsort()
    x_num = snap["data"]["x"][order]
except Exception as e:
    print(f"Error loading snapshots: {e}")
    sys.exit(1)

colors = ['#1565c0', '#b71c1c', '#1b5e20', '#6a1b9a']
styles = ['-', '--', ':', '-.']

# Plot initial condition (common to all bins)
f_dust_init = snap_init["data"]["DustBin_01"][order_init]
ax.plot(x_num_init, f_dust_init, 'k:', linewidth=2.5, label='Initial Condition (t=0)')

# Loop over the 4 bins
for i in range(4):
    # Shift analytical profile: x_shifted = x0 + u_drift * t
    x_shifted = x0 + (u_drift_pc_myr[i] * t_eval)
    
    # Check for boundary wrap-around (periodic boundaries)
    x_shifted_wrapped = np.mod(x_shifted, boxlen)
    
    # Evaluate shifted analytical distribution profile
    eps_analytic = eps_base + delta_eps * np.exp(-((x_eval - x_shifted_wrapped)**2) / (2.0 * sigma_shell**2))
    
    # Plot Analytical curve
    ax.plot(x_eval, eps_analytic, '-', color=colors[i], linewidth=2.5, alpha=0.4,
             label=f'Bin {i+1} Analytic (v_d={u_drift_pc_myr[i]:.3f})')
             
    # Plot Numerical curve
    bin_name = f"DustBin_0{i+1}"
    f_dust_num = snap["data"][bin_name][order]
    ax.plot(x_num, f_dust_num, '--', color=colors[i], linewidth=2.5,
             label=f'Bin {i+1} Numerical ({bin_name})')

# Grid presentation properties
ax.set_xlim(0.0, boxlen)
ax.set_ylim(0.0, delta_eps * 1.2)
ax.set_xlabel('Position x [pc]', fontsize=12)
ax.set_ylabel('Dust Mass Fraction $\epsilon_i$', fontsize=12)
ax.legend(loc='upper right', frameon=False, ncol=2)
ax.grid(True, linestyle=':', alpha=0.3, color='#aaaaaa')
fig.tight_layout()

# Save figure
fig.savefig("dustyshell_multbin.png", dpi=300, bbox_inches="tight", format='png', transparent=True)
print("Saved convergence plot to dustyshell_multbin.png")

# Check level 7 results against the reference solution if it exists
tolerance = {"all": 1e-12}
ref_file = "dustyshell_multbin-ref.dat"
overwrite = "--overwrite" in sys.argv or "-o" in sys.argv or not os.path.exists(ref_file)
visu_ramses.check_solution(snap["data"], "dustyshell_multbin", tolerance=tolerance, overwrite=overwrite)
