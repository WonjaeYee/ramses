import sys
import os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

# def run_simulation(level):
#     """Modify namelist and run the simulation for a given level."""
#     # Read template namelist
#     with open("dustyadvect.nml", "r") as f:
#         lines = f.readlines()
    
#     # Generate modified namelist
#     new_lines = []
#     for line in lines:
#         if "levelmin=" in line:
#             new_lines.append(f"levelmin={level}\n")
#         elif "levelmax=" in line:
#             new_lines.append(f"levelmax={level}\n")
#         elif "ngridmax=" in line:
#             # Set a safe ngridmax to prevent grid overflow
#             new_lines.append(f"ngridmax={max(200, 2**(level + 1))}\n")
#         else:
#             new_lines.append(line)
            
#     temp_nml = f"temp_lvl{level}.nml"
#     with open(temp_nml, "w") as f:
#         f.writelines(new_lines)
        
#     print(f"--- Running simulation for level {level} ---")
#     cmd = ["./ramses_dust_test1d", temp_nml]
#     if "--no-mpi" not in sys.argv:
#         cmd = ["mpirun", "-np", "8"] + cmd
#     # Execute RAMSES test binary
#     subprocess.run(
#         cmd,
#         check=True,
#         stdout=subprocess.DEVNULL,
#         stderr=subprocess.DEVNULL
#     )
    
#     # Remove temporary namelist
#     if os.path.exists(temp_nml):
#         os.remove(temp_nml)

# # Levels to run
# levels = [7, 9, 11, 13]
# results = {}

# # Run simulations and load snapshots
# for lvl in levels:
#     run_simulation(lvl)
    
#     data_init = visu_ramses.load_snapshot(1)
#     data_final = visu_ramses.load_snapshot(2)
    
#     # Extract and sort initial conditions
#     order_init = data_init["data"]["x"].argsort()
#     x_init = data_init["data"]["x"][order_init]
#     rho_init = data_init["data"]["density"][order_init]
#     f_dust_init = data_init["data"]["DustBin_01"][order_init]
#     rho_dust_init = rho_init * f_dust_init
    
#     # Extract and sort final solution
#     order_final = data_final["data"]["x"].argsort()
#     x_final = data_final["data"]["x"][order_final]
#     rho_final = data_final["data"]["density"][order_final]
#     f_dust_final = data_final["data"]["DustBin_01"][order_final]
#     rho_dust_final = rho_final * f_dust_final
    
#     results[lvl] = {
#         "x_init": x_init,
#         "rho_dust_init": rho_dust_init,
#         "x_final": x_final,
#         "rho_dust_final": rho_dust_final,
#         "data_init": data_init,
#         "data_final": data_final
#     }

# # Create plot comparing all levels
# fig, ax = plt.subplots(figsize=(9, 6))

# # Define premium color palette
# colors = {
#     7: "#5c3c92",   # Purple
#     9: "#077b8a",   # Teal
#     11: "#2e8b57",  # Sea Green
#     13: "#d9534f"   # Soft Red
# }

# # Plot analytic solution using the highest resolution initial condition (level 13)
# L = results[13]["data_init"]["data"]["boxlen"]
# t = results[13]["data_final"]["data"]["time"]
# w_x = 1.0

# x_init_ref = results[13]["x_init"]
# rho_dust_init_ref = results[13]["rho_dust_init"]

# # Generate fine grid for plotting a smooth analytic line
# x_fine = np.linspace(0.0, L, 2000)
# x_shifted = (x_fine - w_x * t) % L

# # Extend periodically to avoid boundary interpolation errors
# x_pad = np.concatenate([x_init_ref - L, x_init_ref, x_init_ref + L])
# rho_dust_pad = np.concatenate([rho_dust_init_ref, rho_dust_init_ref, rho_dust_init_ref])

# rho_dust_analytic = np.interp(x_shifted, x_pad, rho_dust_pad)

# # Plot analytic solution
# ax.plot(x_fine, rho_dust_analytic, "-", color="black", linewidth=2.0, label="Analytic")

# # Plot numerical solutions for each level
# for lvl in levels:
#     ax.plot(
#         results[lvl]["x_final"],
#         results[lvl]["rho_dust_final"],
#         label=f"Level {lvl}",
#         color=colors[lvl],
#         linewidth=1.5
#     )

# # Formatting
# ax.set_xlabel("x", fontsize=12)
# ax.set_ylabel("Dust Density", fontsize=12)
# ax.set_title("DustyAdvect: Grid Resolution Convergence", fontsize=14, fontweight="bold")
# ax.grid(True, linestyle=":", alpha=0.6)
# ax.legend(fontsize=10, loc="upper right")
# ax.set_xlim(0.0, L)
# ax.set_ylim(0.0, 0.12)

# fig.tight_layout()
# fig.savefig("dustyadvect.pdf", bbox_inches="tight")

# # Check level 7 results against the reference solution
# tolerance = {"all": 1e-12}
# overwrite = "--overwrite" in sys.argv or "-o" in sys.argv or not os.path.exists("dustyadvect-ref.dat")
# visu_ramses.check_solution(results[7]["data_final"]["data"], "dustyadvect", tolerance=tolerance, overwrite=overwrite)

# ====================================================================
# L2 NORM CONVERGENCE TEST (Lebreuilly et al. 2019)
# ====================================================================
def run_l2_simulation(level, slope):
    """Modify namelist for smooth Gaussian IC and run simulation."""
    with open("dustyadvect.nml", "r") as f:
        lines = f.readlines()
    
    new_lines = []
    in_init_params = False
    
    for line in lines:
        if "&INIT_PARAMS" in line:
            in_init_params = True
            new_lines.append(line)
            new_lines.append("condinit_kind='dustygauss'\n")
            new_lines.append("nregion=1\n")
            new_lines.append("region_type(1)='square'\n")
            new_lines.append("x_center=0.5\n")
            new_lines.append("length_x=1.0\n")
            new_lines.append("d_region=1.0\n")
            new_lines.append("u_region=0.0\n")
            new_lines.append("p_region=1.0\n")
            continue
        
        if in_init_params:
            if "/" in line:
                in_init_params = False
                new_lines.append(line)
            continue
            
        if "levelmin=" in line:
            new_lines.append(f"levelmin={level}\n")
        elif "levelmax=" in line:
            new_lines.append(f"levelmax={level}\n")
        elif "ngridmax=" in line:
            new_lines.append(f"ngridmax={max(200, 2**(level + 1))}\n")
        elif "slope_type=" in line:
            new_lines.append(f"slope_type={slope}\n")
        elif "tout=" in line:
            new_lines.append("tout=0.01\n")
        else:
            new_lines.append(line)
            
    temp_nml = f"temp_l2_lvl{level}_slope{slope}.nml"
    with open(temp_nml, "w") as f:
        f.writelines(new_lines)
        
    print(f"--- Running L2 simulation for level {level}, slope_type={slope} ---")
    print(f"Namelist saved in: {temp_nml}")
    cmd = ["./ramses_dust_test1d", temp_nml]
    if "--no-mpi" not in sys.argv:
        cmd = ["mpirun", "-np", "8"] + cmd
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    
    # if os.path.exists(temp_nml):
    #     os.remove(temp_nml)

l2_levels = [4]
l2_errors = {0: [], 1: []}
l2_profiles = {}

for slope in [0,1]:
    for level in l2_levels:
        run_l2_simulation(level, slope)
        data_final = visu_ramses.load_snapshot(2)
        
        x_final = data_final["data"]["x"]
        rho_final = data_final["data"]["density"]
        f_dust_final = data_final["data"]["DustBin_01"]
        Ncell = 2**level
        rho_dust_final = rho_final * f_dust_final
        
        # Analytical solution shifted by w_x * t (where w_x = 1.0)
        t = data_final["data"]["time"]
        L = data_final["data"]["boxlen"]
        dist = x_final - 0.5 * L - 1.0 * t
        dist = (dist + 0.5 * L) % L - 0.5 * L
        rho_dust_anal = 0.01 + 0.1 * np.exp(- (dist / (0.25 * L))**2)
        
        L2_err = np.sqrt(np.sum((rho_dust_final - rho_dust_anal)**2)/Ncell)
        l2_errors[slope].append(L2_err)
        print(f"  Level {level}, L2 error = {L2_err:.4e}")
        
        # Save profiles for selected levels to compare to analytic
        if slope == 1 and level in l2_levels:
            order = x_final.argsort()
            l2_profiles[level] = {
                "x": x_final[order],
                "rho_dust": rho_dust_final[order]
            }

# Save the L2 convergence plot
fig2, ax2 = plt.subplots(figsize=(8, 6))
dx_values = [1.0 / 2**lvl for lvl in l2_levels]

ax2.loglog(dx_values, l2_errors[0], "o--", color="#d9534f", linewidth=1.5, label="First Order (slope_type=0)")
ax2.loglog(dx_values, l2_errors[1], "s-", color="#077b8a", linewidth=1.5, label="MinMod (slope_type=1)")

# Reference slopes
ax2.loglog(dx_values, l2_errors[0][0] * (np.array(dx_values)/dx_values[0])**1, ":", color="gray", label="1st Order Slope")
ax2.loglog(dx_values, l2_errors[1][0] * (np.array(dx_values)/dx_values[0])**2, "-.", color="gray", label="2nd Order Slope")

ax2.set_xlabel("dx", fontsize=12)
ax2.set_ylabel("L2 Error", fontsize=12)
ax2.set_title("L2 Error Convergence for Gaussian Advection", fontsize=14, fontweight="bold")
ax2.grid(True, which="both", linestyle=":", alpha=0.6)
ax2.legend(fontsize=10, loc="lower right")

fig2.tight_layout()
fig2.savefig("dustyadvect_l2.pdf", bbox_inches="tight")
print("--- L2 convergence plot saved as dustyadvect_l2.pdf ---")

# Save the final solution plot for the Gaussian profile convergence
fig3, ax3 = plt.subplots(figsize=(9, 6))

# Define premium color palette for levels 4 to 9
colors = {
    4: "#d9534f",   # Soft Red
    5: "#f0ad4e",   # Soft Orange
    6: "#5cb85c",   # Soft Green
    7: "#5bc0de",   # Soft Blue
    8: "#0275d8",   # Strong Blue
    9: "#5c3c92"    # Soft Purple
}

# Generate fine grid for plotting a smooth analytic line for the Gaussian
L = 1.0  # Box size
t = 0.01  # L2 simulation end time
x_fine = np.linspace(0.0, L, 2000)
dist_fine = x_fine - 0.5 * L - 1.0 * t
dist_fine = (dist_fine + 0.5 * L) % L - 0.5 * L
rho_dust_analytic_gauss = 0.01 + 0.1 * np.exp(- (dist_fine / (0.25 * L))**2)

# Plot analytic solution
ax3.plot(x_fine, rho_dust_analytic_gauss, "-", color="black", linewidth=2.0, label="Analytic")

# Plot numerical profiles
for lvl in l2_levels:
    if lvl in l2_profiles:
        ax3.plot(
            l2_profiles[lvl]["x"],
            l2_profiles[lvl]["rho_dust"],
            label=f"Level {lvl}",
            color=colors[lvl],
            linewidth=1.5
        )

# Formatting
ax3.set_xlabel("x", fontsize=12)
ax3.set_ylabel("Dust Density", fontsize=12)
ax3.set_title("DustyAdvect: Gaussian Profile Convergence (slope_type=1)", fontsize=14, fontweight="bold")
ax3.grid(True, linestyle=":", alpha=0.6)
ax3.legend(fontsize=10, loc="upper right")
ax3.set_xlim(0.0, L)
ax3.set_ylim(0.0, 0.12)

fig3.tight_layout()
fig3.savefig("dustyadvect_gauss_profile.pdf", bbox_inches="tight")
print("--- Gaussian profile convergence plot saved as dustyadvect_gauss_profile.pdf ---")


