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
fig, ax = plt.subplots(figsize=(5.5, 4))

# Define colors for the outputs (Snapshot 1: initial condition, Snapshots 2, 3, 4, 5: evolution)
output_indices = [1, 2, 3, 4, 5]
colors = {
    1: "#424242",   # Near Black (initial condition)
    2: "#4a148c",   # Dark Purple
    3: "#006064",   # Dark Teal
    4: "#1b5e20",   # Dark Green
    5: "#b71c1c",   # Dark Red
}

# Physical parameters for the analytical Barenblatt-Pattle solution
L = 1.0  # Box size
x_center = 0.5 * L
eps_0 = 0.1
x_c = 0.5 * L
gamma = 1.4
c_s = 1.0  # Enforced isothermal sound speed

xc = 0.2
eps_0 = 0.1
ts = 0.1
cs = 1
C = (eps_0*xc/np.sqrt(6))**(2./3.)
t0 = C**3. / (ts*cs**2.*eps_0**3.)

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
    
    if t != 0:
        # Evaluate analytical solution at time t
        t_eff = t + t0
        eps_analytic = (ts*cs**2.*t_eff)**(-1./3.) * (C - 1./6.*(x_num-0.5)**2./(ts*cs**2*t_eff)**(2./3.))
        eps_analytic = np.maximum(0.0, eps_analytic)
        rho_dust_analytic = rho_num * eps_analytic
        # Plot analytical solution (solid line)
        label_analytic = f"t={t:.1f} (Analytic)"
        ax.plot(
            x_num, 
            rho_dust_analytic, 
            "-", 
            color=colors[idx], 
            linewidth=2.5,
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
        linewidth=2.5,
        label=label_numerical
    )

# Formatting
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("Dust Density", fontsize=12)
ax.grid(True, linestyle=":", alpha=0.3, color="#aaaaaa")
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(0.0, L)
ax.set_ylim(0.0, 0.12)

fig.tight_layout()
fig.savefig("dustydiffuse.pdf", bbox_inches="tight", transparent=True)
fig.savefig("dustydiffuse.png", dpi=300, bbox_inches="tight", transparent=True)
print("Saved convergence plot to dustydiffuse.pdf / dustydiffuse.png")
