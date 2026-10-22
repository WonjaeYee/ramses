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
    print("--- Running DustyGauss simulation ---")
    subprocess.run(
        ["mpirun", "-np", "8", "./ramses_dust_test1d", "dustygauss.nml"],
        check=True
    )

# Run the simulation
run_simulation()

# =====================================================================
# 1. PHYSICAL CONSTANTS & ENVIRONMENTAL ASSUMPTIONS
# =====================================================================
u_drift_pc_myr = 0.77692  # Advection velocity in pc/Myr (matching the namelist u_region)
boxlen      = 10.0                 # 10 pc domain
x0          = 2.0                  # Initial shell center (pc)
sigma_shell = 0.4                  # Initial characteristic width (pc)
eps_base    = 1.0e-8               # Ambient dust background
delta_eps   = 1.0e-4                 # Shell injection amplitude

# Setting the 4 desired snapshot targets
x_eval = np.linspace(0.0, boxlen, 300)

# Initial Gaussian Profile function
def initial_profile(x):
    return eps_base + delta_eps * np.exp(-((x - x0)**2) / (2.0 * sigma_shell**2))

# =====================================================================
# 2. PLOT OVERLAY
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

# Loop and calculate spatial shifts for each snapshot
colors = ['#1565c0', '#b71c1c', '#1b5e20', '#6a1b9a']
for i in range(4):
    snap_idx = i + 2
    try:
        snap = visu_ramses.load_snapshot(snap_idx)
        t = snap["info"]["time"]  # Read actual time from snapshot info
        
        # Pure hyperbolic advection shift: x_shifted = x0 + u_drift * t
        x_shifted = x0 + (u_drift_pc_myr * t)
        
        # Evaluate shifted analytical distribution profile
        eps_analytic = eps_base + delta_eps * np.exp(-((x_eval - x_shifted)**2) / (2.0 * sigma_shell**2))
        
        ax.plot(x_eval, eps_analytic, '-', color=colors[i], linewidth=2.5, alpha=0.5,
                 label=f't = {t:.2f} Myr (Analytic)')
                 
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
ax.set_ylabel('Dust Mass Fraction $\\epsilon_1$', fontsize=12)
ax.legend(loc='upper right', frameon=False)
ax.grid(True, linestyle=':', alpha=0.3, color='#aaaaaa')
fig.tight_layout()

# Save figure
fig.savefig("dustygauss.pdf", bbox_inches="tight", transparent=True)
fig.savefig("dustygauss.png", bbox_inches="tight", dpi=300, transparent=True)
print("Saved convergence plot to dustygauss.pdf and dustygauss.png")
