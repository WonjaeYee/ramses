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

# ====================================================================
# TOP-HAT ADVECTION CONVERGENCE TEST
# Dust top-hat advected at w_x = 1.0 for one full period (t = 1.0).
# Analytic solution at t=1.0 is identical to t=0 (periodic, L=1).
# ====================================================================

def run_simulation(level):
    """Generate a level-specific namelist and run the simulation."""
    with open("dustyadvect.nml", "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if "levelmin=" in line:
            new_lines.append(f"levelmin={level}\n")
        elif "levelmax=" in line:
            new_lines.append(f"levelmax={level}\n")
        elif "ngridmax=" in line:
            new_lines.append(f"ngridmax={max(200, 2**(level + 1))}\n")
        else:
            new_lines.append(line)

    temp_nml = f"temp_lvl{level}.nml"
    with open(temp_nml, "w") as f:
        f.writelines(new_lines)

    print(f"--- Running simulation for level {level} ({2**level} cells) ---")
    cmd = ["mpirun", "-np", "4", "./ramses_dust_test1d", temp_nml]
    if "--no-mpi" in sys.argv:
        cmd = ["./ramses_dust_test1d", temp_nml]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if os.path.exists(temp_nml):
        os.remove(temp_nml)


levels = [7, 9, 11, 13]
results = {}

for lvl in levels:
    run_simulation(lvl)

    data_init  = visu_ramses.load_snapshot(1)
    data_final = visu_ramses.load_snapshot(2)

    order_init  = data_init["data"]["x"].argsort()
    order_final = data_final["data"]["x"].argsort()

    x_init      = data_init["data"]["x"][order_init]
    rho_init    = data_init["data"]["density"][order_init]
    fd_init     = data_init["data"]["DustBin_01"][order_init]

    x_final     = data_final["data"]["x"][order_final]
    rho_final   = data_final["data"]["density"][order_final]
    fd_final    = data_final["data"]["DustBin_01"][order_final]

    results[lvl] = {
        "x_init":      x_init,
        "rho_dust_init":  rho_init * fd_init,
        "x_final":     x_final,
        "rho_dust_final": rho_final * fd_final,
        "data_init":   data_init,
        "data_final":  data_final,
    }

# ---- analytic solution (top-hat shifted back to t=0 position) ----------
# Use the highest-resolution initial condition as the reference shape.
L    = results[13]["data_init"]["data"]["boxlen"]
t    = results[13]["data_final"]["data"]["time"]
w_x  = 1.0

x_ref      = results[13]["x_init"]
rho_ref    = results[13]["rho_dust_init"]

x_fine    = np.linspace(0.0, L, 4000)
x_shifted = (x_fine - w_x * t) % L

x_pad     = np.concatenate([x_ref - L, x_ref, x_ref + L])
rho_pad   = np.tile(rho_ref, 3)
rho_analytic = np.interp(x_shifted, x_pad, rho_pad)

# ---- plot ---------------------------------------------------------------
colors = {7: "C0", 9: "C1", 11: "C2", 13: "C3"}
linestyles = {7: "-", 9: "--", 11: ":", 13: "-."}

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(x_fine, rho_analytic, "-", color="k", label="Analytic", zorder=5)

for lvl in levels:
    ax.plot(
        results[lvl]["x_final"],
        results[lvl]["rho_dust_final"],
        label=f"Level {lvl}  ({2**lvl} cells)",
        color=colors[lvl],
        linestyle=linestyles[lvl],
        drawstyle="steps-mid",
    )

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Dust density $\rho_\mathrm{d}$")
ax.legend(loc="center")
ax.set_xlim(0.0, L)
tidy(ax)

fig.tight_layout()
save(fig, "dustyadvect")
print("--- Top-hat convergence plot saved as dustyadvect.pdf / dustyadvect.png ---")

# Reference solution check (level 7)
tolerance = {"all": 1e-12}
overwrite = "--overwrite" in sys.argv or "-o" in sys.argv or not os.path.exists("dustyadvect-ref.dat")
visu_ramses.check_solution(results[7]["data_final"]["data"], "dustyadvect", tolerance=tolerance, overwrite=overwrite)
