import sys
import os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

def run_simulation():
    """Run the simulation once using the parameter file."""
    print("--- Running DustyDiffuse simulation ---")
    subprocess.run(
        ["./ramses_dust_test1d", "dustydiffuse.nml"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

# Run the simulation
run_simulation()

# Create plot
fig, ax = plt.subplots(figsize=(9, 6))

# Define colors for the outputs (Snapshot 1: initial condition, Snapshots 2, 3, 4, 5: evolution)
output_indices = [1, 2, 3, 4, 5]
colors = {
    1: "#777777",   # Dark Gray (initial condition)
    2: "#5c3c92",   # Purple
    3: "#077b8a",   # Teal
    4: "#2e8b57",   # Sea Green
    5: "#d9534f"   # Soft Red
}

# Physical parameters for the analytical Barenblatt-Pattle solution
L = 1.0  # Box size
x_center = 0.5 * L
eps_0 = 0.1
x_c = 0.5 * L
gamma = 1.4
c_s = 1.0  # Enforced isothermal sound speed

eta_0 = 0.1  # Enforced diffusion coefficient (t_s * c_s**2)
t_0 = (x_c**2) / (6.0 * eps_0 * eta_0)
C = eps_0 * (eta_0 * t_0)**(1.0/3.0)

# Generate a fine grid for plotting a smooth analytic line
x_fine = np.linspace(0.0, L, 2000)

for idx in output_indices:
    # Load snapshot data
    snapshot = visu_ramses.load_snapshot(idx)
    t = snapshot["data"]["time"]
    
    # Sort and extract numerical solution
    order = snapshot["data"]["x"].argsort()
    x_num = snapshot["data"]["x"][order]
    rho_num = snapshot["data"]["density"][order]
    f_dust_num = snapshot["data"]["DustBin_01"][order]
    rho_dust_num = rho_num * f_dust_num
    
    # Evaluate analytical solution at time t
    S = eta_0 * (t + t_0)
    dist_periodic = (x_fine - x_center + L/2.0) % L - L/2.0
    term = C - dist_periodic**2 / (6.0 * S**(2.0/3.0))
    rho_dust_analytic = (1.0 / S**(1.0/3.0)) * np.maximum(term, 0.0)
    
    # Plot analytical solution (solid line)
    label_analytic = "Analytic Initial" if idx == 1 else f"t={t:.1f} (Analytic)"
    ax.plot(
        x_fine, 
        rho_dust_analytic, 
        "-", 
        color=colors[idx], 
        linewidth=2.0, 
        alpha=0.4,
        label=label_analytic
    )
    
    # Plot numerical solution (dashed line)
    label_numerical = "Numerical Initial" if idx == 1 else f"Numerical t={t:.1f}"
    ax.plot(
        x_num,
        rho_dust_num,
        "--",
        color=colors[idx],
        linewidth=1.5,
        label=label_numerical
    )

# Formatting
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("Dust Density", fontsize=12)
ax.set_title("DustyDiffuse: Non-linear Diffusion Evolution", fontsize=14, fontweight="bold")
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(0.0, L)
ax.set_ylim(0.0, 0.12)

fig.tight_layout()
fig.savefig("dustydiffuse.pdf", bbox_inches="tight")
print("Saved convergence plot to dustydiffuse.pdf")
