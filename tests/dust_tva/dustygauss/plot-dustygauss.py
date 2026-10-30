import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../utils/py/calima")))
from calima_style import use_style, tidy, save  # noqa: E402
use_style()

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
fig, ax = plt.subplots(figsize=(10, 5))

# Plot analytical initial condition
ax.plot(x_eval, initial_profile(x_eval), 'k-', alpha=0.4, label='Initial (Analytic)')

# Load and plot numerical initial condition (Snapshot 1)
try:
    snap_init = visu_ramses.load_snapshot(1)
    order = snap_init["data"]["x"].argsort()
    x_num = snap_init["data"]["x"][order]
    f_dust_num = snap_init["data"]["DustBin_01"][order]
    ax.plot(x_num, f_dust_num, 'k--', label='Initial (Numerical)')
except Exception as e:
    print(f"Error loading initial snapshot: {e}")

# Loop and calculate spatial shifts for each snapshot
colors = ['C0', 'C1', 'C2', 'C3']
for i in range(4):
    snap_idx = i + 2
    try:
        snap = visu_ramses.load_snapshot(snap_idx)
        t = snap["info"]["time"]  # Read actual time from snapshot info
        
        # Pure hyperbolic advection shift: x_shifted = x0 + u_drift * t
        x_shifted = x0 + (u_drift_pc_myr * t)
        
        # Evaluate shifted analytical distribution profile
        eps_analytic = eps_base + delta_eps * np.exp(-((x_eval - x_shifted)**2) / (2.0 * sigma_shell**2))
        
        ax.plot(x_eval, eps_analytic, '-', color=colors[i], alpha=0.5,
                 label=f't = {t:.2f} Myr (Analytic)')
                 
        order = snap["data"]["x"].argsort()
        x_num = snap["data"]["x"][order]
        f_dust_num = snap["data"]["DustBin_01"][order]
        ax.plot(x_num, f_dust_num, '--', color=colors[i],
                 label=f't = {t:.2f} Myr (Numerical)')
    except Exception as e:
        print(f"Error loading snapshot {snap_idx}: {e}")

# Grid presentation properties
ax.set_xlim(0.0, boxlen)
ax.set_ylim(0.0, delta_eps * 1.2)
ax.set_xlabel('Position x [pc]')
ax.set_ylabel('Dust Mass Fraction $\\epsilon_1$')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
tidy(ax)
fig.tight_layout()

# Save figure
save(fig, "dustygauss")
print("Saved convergence plot to dustygauss.pdf and dustygauss.png")
