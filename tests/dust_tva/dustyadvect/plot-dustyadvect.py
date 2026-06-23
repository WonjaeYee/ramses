import sys
import os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../visu")))
import visu_ramses

def run_simulation(level):
    """Modify namelist and run the simulation for a given level."""
    # Read template namelist
    with open("dustyadvect.nml", "r") as f:
        lines = f.readlines()
    
    # Generate modified namelist
    new_lines = []
    for line in lines:
        if "levelmin=" in line:
            new_lines.append(f"levelmin={level}\n")
        elif "levelmax=" in line:
            new_lines.append(f"levelmax={level}\n")
        elif "ngridmax=" in line:
            # Set a safe ngridmax to prevent grid overflow
            new_lines.append(f"ngridmax={max(200, 2**(level + 1))}\n")
        else:
            new_lines.append(line)
            
    temp_nml = f"temp_lvl{level}.nml"
    with open(temp_nml, "w") as f:
        f.writelines(new_lines)
        
    print(f"--- Running simulation for level {level} ---")
    # Execute RAMSES test binary
    subprocess.run(
        ["./ramses_dust_test1d", temp_nml],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    
    # Remove temporary namelist
    if os.path.exists(temp_nml):
        os.remove(temp_nml)

# Levels to run
levels = [7, 9, 11, 13]
results = {}

# Run simulations and load snapshots
for lvl in levels:
    run_simulation(lvl)
    
    data_init = visu_ramses.load_snapshot(1)
    data_final = visu_ramses.load_snapshot(2)
    
    # Extract and sort initial conditions
    order_init = data_init["data"]["x"].argsort()
    x_init = data_init["data"]["x"][order_init]
    rho_init = data_init["data"]["density"][order_init]
    f_dust_init = data_init["data"]["DustBin_01"][order_init]
    rho_dust_init = rho_init * f_dust_init
    
    # Extract and sort final solution
    order_final = data_final["data"]["x"].argsort()
    x_final = data_final["data"]["x"][order_final]
    rho_final = data_final["data"]["density"][order_final]
    f_dust_final = data_final["data"]["DustBin_01"][order_final]
    rho_dust_final = rho_final * f_dust_final
    
    results[lvl] = {
        "x_init": x_init,
        "rho_dust_init": rho_dust_init,
        "x_final": x_final,
        "rho_dust_final": rho_dust_final,
        "data_init": data_init,
        "data_final": data_final
    }

# Create plot comparing all levels
fig, ax = plt.subplots(figsize=(9, 6))

# Define premium color palette
colors = {
    7: "#5c3c92",   # Purple
    9: "#077b8a",   # Teal
    11: "#2e8b57",  # Sea Green
    13: "#d9534f"   # Soft Red
}

# Plot analytic solution using the highest resolution initial condition (level 13)
L = results[13]["data_init"]["data"]["boxlen"]
t = results[13]["data_final"]["data"]["time"]
w_x = 1.0

x_init_ref = results[13]["x_init"]
rho_dust_init_ref = results[13]["rho_dust_init"]

# Generate fine grid for plotting a smooth analytic line
x_fine = np.linspace(0.0, L, 2000)
x_shifted = (x_fine - w_x * t) % L

# Extend periodically to avoid boundary interpolation errors
x_pad = np.concatenate([x_init_ref - L, x_init_ref, x_init_ref + L])
rho_dust_pad = np.concatenate([rho_dust_init_ref, rho_dust_init_ref, rho_dust_init_ref])

rho_dust_analytic = np.interp(x_shifted, x_pad, rho_dust_pad)

# Plot analytic solution
ax.plot(x_fine, rho_dust_analytic, "-", color="black", linewidth=2.0, label="Analytic")

# Plot numerical solutions for each level
for lvl in levels:
    ax.plot(
        results[lvl]["x_final"],
        results[lvl]["rho_dust_final"],
        label=f"Level {lvl}",
        color=colors[lvl],
        linewidth=1.5
    )

# Formatting
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("Dust Density", fontsize=12)
ax.set_title("DustyAdvect: Grid Resolution Convergence", fontsize=14, fontweight="bold")
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(0.0, L)
ax.set_ylim(0.0, 0.12)

fig.tight_layout()
fig.savefig("dustyadvect.pdf", bbox_inches="tight")

# Check level 7 results against the reference solution
tolerance = {"all": 1e-12}
visu_ramses.check_solution(results[7]["data_final"]["data"], "dustyadvect", tolerance=tolerance)
