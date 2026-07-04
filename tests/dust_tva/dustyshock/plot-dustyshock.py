import sys
import os
import subprocess
import numpy as np
from scipy.optimize import brentq
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../visu")))
import visu_ramses

# ---------------------------------------------------------------------------
# Physical / coupling constants for this test  (from dustyshock.nml)
# ---------------------------------------------------------------------------
K     = 1000.0  # Drag coefficient
GAMMA = 1.4     # Adiabatic index
EPS   = 0.5     # Uniform initial dust fraction

# Initial left / right states  (d_region, u_region, p_region in nml)
RHO_MIX_L, U_L, P_L = 1.0,   0.0, 1.0
RHO_MIX_R, U_R, P_R = 0.125, 0.0, 0.1
X0 = 0.5   # diaphragm position

# ---------------------------------------------------------------------------
# Analytic Riemann solver for the dusty Sod shock tube
# ---------------------------------------------------------------------------
# The mixture behaves like an effective single fluid whose sound speed is
#   c_eff = sqrt(gamma * P / rho_mix) = c_s * sqrt(1 - eps)
# because P is the gas pressure while rho_mix = rho_gas / (1 - eps).
# The effective gamma is the same (gamma = 1.4); only the mixture density
# appears in the denominator of c².
# In the perfectly-coupled analytic limit K → ∞:  v_drift → 0, so
#   v_gas = v_dust = v_mix.
# ---------------------------------------------------------------------------

def _pressure_function(P_star, P_k, rho_k, c_k, gamma):
    """
    Contribution from wave k to the pressure equation.
    Rarefaction  when P_star <= P_k, shock otherwise.
    """
    if P_star <= P_k:
        # Isentropic rarefaction
        return (2.0 * c_k / (gamma - 1.0)) * (
            (P_star / P_k) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0
        )
    else:
        # Rankine-Hugoniot shock
        A_k = 2.0 / ((gamma + 1.0) * rho_k)
        B_k = (gamma - 1.0) / (gamma + 1.0) * P_k
        return (P_star - P_k) * np.sqrt(A_k / (P_star + B_k))


def sod_analytic(rho_L, u_L, P_L, rho_R, u_R, P_R, gamma, x_arr, x0, t):
    """
    Exact solution of the Riemann problem at time t on grid x_arr.

    Uses the *mixture* density (rho_L, rho_R) so that
        c_L = sqrt(gamma * P_L / rho_L)
    already incorporates the c_eff = cs * sqrt(1-eps) correction.

    Returns
    -------
    rho_mix : mixture density
    u_mix   : mixture (= gas = dust) velocity  [perfectly coupled limit]
    P_gas   : gas pressure
    """
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)

    # ---- find P* by Brent's root-finding -----------------------------------
    def pressure_eq(P_star):
        return (
            _pressure_function(P_star, P_L, rho_L, c_L, gamma)
            + _pressure_function(P_star, P_R, rho_R, c_R, gamma)
            + (u_R - u_L)
        )

    P_star = brentq(pressure_eq, 1.0e-12, max(P_L, P_R) * 1.0e4, xtol=1.0e-12)

    # ---- u* ----------------------------------------------------------------
    u_star = 0.5 * (u_L + u_R) + 0.5 * (
        _pressure_function(P_star, P_R, rho_R, c_R, gamma)
        - _pressure_function(P_star, P_L, rho_L, c_L, gamma)
    )

    # ---- left star state (rarefaction → isentropic) ------------------------
    rho_star_L = rho_L * (P_star / P_L) ** (1.0 / gamma)
    c_star_L   = np.sqrt(gamma * P_star / rho_star_L)

    # ---- right star state (shock → Rankine-Hugoniot) -----------------------
    rho_star_R = rho_R * (
        (P_star / P_R + (gamma - 1.0) / (gamma + 1.0))
        / ((gamma - 1.0) / (gamma + 1.0) * P_star / P_R + 1.0)
    )

    # ---- wave speeds -------------------------------------------------------
    xi_head_L = u_L - c_L          # head of left rarefaction fan
    xi_tail_L = u_star - c_star_L  # tail of left rarefaction fan
    xi_contact = u_star             # contact discontinuity
    S_R = u_R + c_R * np.sqrt(
        (gamma + 1.0) / (2.0 * gamma) * P_star / P_R
        + (gamma - 1.0) / (2.0 * gamma)
    )                               # right shock speed

    # ---- map onto x --------------------------------------------------------
    xi = (x_arr - x0) / t   # similarity variable

    rho_out = np.empty_like(xi)
    u_out   = np.empty_like(xi)
    P_out   = np.empty_like(xi)

    # Masks for each region
    m1 = xi <= xi_head_L
    m2 = (xi > xi_head_L) & (xi <= xi_tail_L)
    m3 = (xi > xi_tail_L) & (xi <= xi_contact)
    m4 = (xi > xi_contact) & (xi <= S_R)
    m5 = xi > S_R

    # Region 1: undisturbed left state
    rho_out[m1] = rho_L;   u_out[m1] = u_L;      P_out[m1] = P_L

    # Region 2: left rarefaction fan (isentropic, self-similar)
    if np.any(m2):
        s2       = xi[m2]
        u_fan    = 2.0 / (gamma + 1.0) * (
            (gamma - 1.0) / 2.0 * u_L + c_L + s2
        )
        c_fan    = u_fan - s2                          # c = u - xi in fan
        rho_out[m2] = rho_L * (c_fan / c_L) ** (2.0 / (gamma - 1.0))
        u_out[m2]   = u_fan
        P_out[m2]   = P_L   * (c_fan / c_L) ** (2.0 * gamma / (gamma - 1.0))

    # Region 3: left star region
    rho_out[m3] = rho_star_L;  u_out[m3] = u_star;  P_out[m3] = P_star

    # Region 4: right star region
    rho_out[m4] = rho_star_R;  u_out[m4] = u_star;  P_out[m4] = P_star

    # Region 5: undisturbed right state
    rho_out[m5] = rho_R;   u_out[m5] = u_R;      P_out[m5] = P_R

    return rho_out, u_out, P_out


# ---------------------------------------------------------------------------
# TVA velocity decomposition  (numerical solution)
# ---------------------------------------------------------------------------

def derive_velocities(x, rho, eps, pressure, v_mix):
    """
    Derive gas and dust velocities from mixture variables (TVA).

        t_s     = eps * (1-eps) * rho / K
        v_drift = t_s * nabla(P_gas) / rho_gas
        v_gas   = v_mix - eps * v_drift
        v_dust  = v_mix + (1-eps) * v_drift
    """
    rho_gas = (1.0 - eps) * rho
    t_s     = eps * (1.0 - eps) * rho / K
    nabla_P = np.gradient(pressure, x)
    v_drift = t_s * nabla_P / rho_gas
    return v_mix - eps * v_drift, v_mix + (1.0 - eps) * v_drift


# ---------------------------------------------------------------------------
# Run the simulation
# ---------------------------------------------------------------------------
def run_simulation():
    """Run the dustyshock simulation using the parameter file."""
    print("--- Running DustyShock simulation ---")
    subprocess.run(
        ["./ramses_dust_test1d", "dustyshock.nml"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

run_simulation()

# ---------------------------------------------------------------------------
# Load the final output (output_00002, t ≈ 0.245)
# ---------------------------------------------------------------------------
# RAMSES writes output_00001 as the initial condition and output_00002 as
# the scheduled output at tout = 0.245.
snapshot = visu_ramses.load_snapshot(2)
t = snapshot["data"]["time"]
print(f"Plotting snapshot at t = {t:.4f}")

order    = snapshot["data"]["x"].argsort()
x        = snapshot["data"]["x"][order]
rho      = snapshot["data"]["density"][order]
pressure = snapshot["data"]["pressure"][order]
eps      = snapshot["data"]["DustBin_01"][order]
v_mix    = snapshot["data"]["velocity_x"][order]

# Derived numerical quantities
rho_gas  = (1.0 - eps) * rho
rho_dust = eps * rho
v_gas, v_dust = derive_velocities(x, rho, eps, pressure, v_mix)

# ---------------------------------------------------------------------------
# Analytic solution on a fine grid
# ---------------------------------------------------------------------------
x_fine = np.linspace(x.min(), x.max(), 4000)

rho_mix_a, u_mix_a, P_gas_a = sod_analytic(
    RHO_MIX_L, U_L, P_L,
    RHO_MIX_R, U_R, P_R,
    GAMMA, x_fine, X0, t,
)

# In the perfectly-coupled limit eps is passively advected and stays at EPS=0.5
rho_gas_a  = (1.0 - EPS) * rho_mix_a
rho_dust_a = EPS          * rho_mix_a
# v_drift → 0 as K → ∞, so v_gas = v_dust = v_mix
v_mix_a    = u_mix_a

# ---------------------------------------------------------------------------
# Plotting — 2×2 panel figure
# ---------------------------------------------------------------------------
color_gas    = "#1f77b4"   # muted blue  (numerical gas)
color_dust   = "#d62728"   # brick red   (numerical dust)
color_prs    = "#ff7f0e"   # orange      (numerical pressure)
color_analy  = "#444444"   # dark grey   (analytic)

lw_num   = 1.8
lw_ana   = 1.5
ls_ana   = "--"

fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
ax_rho_gas  = axes[0, 0]
ax_prs      = axes[1, 0]
ax_rho_dust = axes[0, 1]
ax_vel      = axes[1, 1]

# --- Top-left : gas density -------------------------------------------------
ax_rho_gas.plot(x,      rho_gas,   color=color_gas,   lw=lw_num, label="Numerical")
ax_rho_gas.plot(x_fine, rho_gas_a, color=color_analy, lw=lw_ana, ls=ls_ana, label="Analytic")
ax_rho_gas.set_ylabel(r"$\rho_\mathrm{gas}$", fontsize=13)
ax_rho_gas.set_title("Gas Density", fontsize=12, fontweight="bold")
ax_rho_gas.grid(True, linestyle=":", alpha=0.5)
ax_rho_gas.legend(fontsize=10, loc="upper right")

# --- Bottom-left : gas pressure ---------------------------------------------
ax_prs.plot(x,      pressure, color=color_prs,   lw=lw_num, label="Numerical")
ax_prs.plot(x_fine, P_gas_a,  color=color_analy, lw=lw_ana, ls=ls_ana, label="Analytic")
ax_prs.set_xlabel("$x$", fontsize=13)
ax_prs.set_ylabel(r"$P_\mathrm{gas}$", fontsize=13)
ax_prs.set_title("Gas Pressure", fontsize=12, fontweight="bold")
ax_prs.grid(True, linestyle=":", alpha=0.5)
ax_prs.legend(fontsize=10, loc="upper right")

# --- Top-right : dust density -----------------------------------------------
ax_rho_dust.plot(x,      rho_dust,   color=color_dust,  lw=lw_num, label="Numerical")
ax_rho_dust.plot(x_fine, rho_dust_a, color=color_analy, lw=lw_ana, ls=ls_ana, label="Analytic")
ax_rho_dust.set_ylabel(r"$\rho_\mathrm{dust}$", fontsize=13)
ax_rho_dust.set_title("Dust Density", fontsize=12, fontweight="bold")
ax_rho_dust.grid(True, linestyle=":", alpha=0.5)
ax_rho_dust.legend(fontsize=10, loc="upper right")

# --- Bottom-right : gas & dust velocities -----------------------------------
# Numerical: gas (solid blue) and dust (dashed red) differ via v_drift
ax_vel.plot(x, v_gas,  color=color_gas,   lw=lw_num, ls="-",   label=r"$v_\mathrm{gas}$ (num)")
ax_vel.plot(x, v_dust, color=color_dust,  lw=lw_num, ls="--",  label=r"$v_\mathrm{dust}$ (num)")
# Analytic: perfectly coupled → single mixture velocity
ax_vel.plot(x_fine, v_mix_a, color=color_analy, lw=lw_ana, ls=ls_ana,
            label=r"$v_\mathrm{mix}$ (analytic, $K\!\to\!\infty$)")
ax_vel.set_xlabel("$x$", fontsize=13)
ax_vel.set_ylabel("Velocity", fontsize=13)
ax_vel.set_title("Gas & Dust Velocities", fontsize=12, fontweight="bold")
ax_vel.grid(True, linestyle=":", alpha=0.5)
ax_vel.legend(fontsize=10, loc="upper right")

# Global formatting
fig.suptitle(
    rf"DustyShock — $t = {t:.4f}$  "
    r"($c_\mathrm{eff} = c_s\sqrt{1-\varepsilon}$, TVA)",
    fontsize=14,
    fontweight="bold",
    y=1.01,
)
for ax in axes.flat:
    ax.set_xlim(x.min(), x.max())

fig.tight_layout()
outfile = "dustyshock.pdf"
fig.savefig(outfile, bbox_inches="tight")
print(f"Saved plot to {outfile}")
