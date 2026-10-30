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
fig, ax = plt.subplots(figsize=(7, 5))

# Define colors for the outputs (Snapshot 1: initial condition, Snapshots 2, 3, 4, 5: evolution)
output_indices = [1, 2, 3, 4, 5]
colors = {
    1: "k",    # initial condition
    2: "C0",
    3: "C1",
    4: "C2",
    5: "C3",
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
        label=label_numerical
    )

# Formatting
ax.set_xlabel("x")
ax.set_ylabel("Dust Density")
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
ax.set_xlim(0.0, L)
ax.set_ylim(0.0, 0.12)
tidy(ax)

fig.tight_layout()
save(fig, "dustydiffuse")
print("Saved convergence plot to dustydiffuse.pdf / dustydiffuse.png")
