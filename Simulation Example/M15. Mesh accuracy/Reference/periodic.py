__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device.poisson_linear import Solver as PoissonSolver
from qtcad.device.poisson_linear import SolverParams as PoissonSolverParams

# Set up paths
script_dir = Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"
path_out.mkdir(exist_ok=True)

# Load the mesh
scaling = 1e-9
file = path_mesh / "periodic_fdsoi.msh"
mesh = Mesh(scaling, file)

# Define gate biases
back_gate_bias = -0.5
N = 3  # Number of gates
B = 1e-1  # 100 mV amplitude
tau = 1e-9  # 1 ns period
Omega = 2 * np.pi / tau


# define the gate bias function:
def phi_j(t, j):
    """
    Calculates the time-dependent gate bias of the clavier gate j.

    Args:
        t (float): Time in seconds.
        j (int): Index of the clavier gate.
    Returns:
        float: Time-dependent gate bias.

    """
    return B * np.cos(Omega * t - ((2 * np.pi * j) / N))


# Define the device object
dvc = Device(mesh, conf_carriers="e")
dvc.set_temperature(0.1)

# Create the regions
dvc.new_region("oxide", mt.SiO2)
dvc.new_region("gate_oxide", mt.HfO2)
dvc.new_region("buried_oxide", mt.SiO2)
dvc.new_region("channel", mt.Si)

# Periodic boundary conditions:
dvc.new_periodic_bnd("left_bnd", "right_bnd")

# Define the gate biases for the first clavier gate activation
phi_A0 = 1
phi_A1 = 0
phi_A2 = 0
phi_back_gate = 0

# Apply Dirichlet boundary conditions
dvc.new_dirichlet_bnd("A0", phi_A0)
dvc.new_dirichlet_bnd("A1", phi_A1)
dvc.new_dirichlet_bnd("A2", phi_A2)
dvc.new_dirichlet_bnd("back_gate_bnd", phi_back_gate)

# Define parameters for frozen boundary to be applied later to the back gate
frozen_values = {
    "material": mt.Si,
    "doping": 1e15 * 1e6,  # Doping in cm⁻³ converted to m⁻³
    "doping_type": "n",
    "binding_energy": 46e-3 * ct.e,  # Binding energy in joules
}

# Define typical metallic work function (mid-gap) in joules
Ew = mt.Si.Eg / 2 + mt.Si.chi

# Configure the linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence tolerance for the error
params_poisson.maxiter = 200

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson)

# Clavier gate biases for each activation
activate_bias = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
Na = len(activate_bias)

# Define the linecut
x, y, z = dvc.mesh.glob_nodes.T
ymin = np.min(y)
ymax = np.max(y)
height = -1e-9  # 1 nm below the silicon channel and the gate oxide

# Initialize a list for saving the linecuts
potential_linecuts = []

for idn, bias in enumerate(activate_bias):
    print(f"Calculating potential for activation {idn + 1} of {Na}")
    # Update gate biases
    dvc.set_applied_potential("A0", bias[0])
    dvc.set_applied_potential("A1", bias[1])
    dvc.set_applied_potential("A2", bias[2])
    dvc.set_applied_potential("back_gate_bnd", bias[3])

    # Solve
    poisson_slv.solve()

    # Calculate the confinement potential
    dvc.set_V_from_phi()

    # Save potential and confinement potential
    distance, linecut_pot = an.linecut(
        dvc.mesh, dvc.phi, (0, ymin, height), (0, ymax, height), method="pyvista"
    )
    potential_linecuts.append([distance, linecut_pot])

# Save the linecuts
output_file = path_out / "potential_linecuts.npy"
np.save(output_file, np.array(potential_linecuts))

# Load the linecuts
potential_linecuts = np.load(output_file)


def confinement_potential(
    t, potential_linecuts=potential_linecuts, Ew=Ew, frozen_values=frozen_values
):
    """
    Calculates the potential along the linecut.

    Args:
        t (float): Time in seconds.
        potential_linecuts (2d array): Pre-computed potential linecuts for
            each gate activation.
        Ew (float): Metallic work function in joules.
        frozen_values (dict): Values of the frozen boundary conditions.

    Returns:
        distance (1d array): Distance along the linecut in meters.
        potential (1d array): Electrostatic potential linecut at time t (in eV).
        confinement_pot (1d array): Confinement potential linecut at time t (in eV).
    """

    # Initialize potentials
    potential = 0
    # reference potential
    phi_ref = dvc.phi_F[0, 0]  # Reference potential from the device
    for j in range(N):
        # Apply the gate boundary condition for each clavier gate
        potential += (phi_j(t, j) + phi_ref - Ew / ct.e) * potential_linecuts[j, 1, :]

    # Frozen boundary condition at the back gate
    T = 0.1  # Temperature in K
    gc = frozen_values["material"].gc
    mc = frozen_values["material"].mc
    chi = frozen_values["material"].chi

    doping = frozen_values["doping"]
    binding_energy = frozen_values["binding_energy"]

    x = gc * np.sqrt(mc * 3) * ((T) ** 1.5 * ct.C0)
    EC = binding_energy / 2.0 + ct.kB * T / 2.0 * np.log(2 * x / doping)
    phi_bi = -(EC + chi) / ct.e

    # Apply back gate contribution
    potential += (back_gate_bias + phi_ref + phi_bi) * potential_linecuts[3, 1, :]
    confinement_pot = -ct.e * (potential - phi_ref) - chi

    # Return distance, electrostatic potential, and confinement potential
    return potential_linecuts[0, 0, :], potential, confinement_pot / ct.e


# calculate the confinement potential at t=tau/3
dist, phi_t, V_t = confinement_potential(tau / 3)
# Plot the potential linecut
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Distance along linecut [nm]")
ax.set_ylabel("Electrostatic potential ($\\varphi$) [mV]")
ax.plot(dist / 1e-9, phi_t / 1e-3, "-r")
ax.axhline(phi_t[0] / 1e-3, color="black", lw=0.5, ls="--", alpha=0.5)

# plot the confinement potential
ax2 = ax.twinx()
ax2.set_ylabel("Confinement potential ($V$) [meV]", color="blue")
ax2.plot(dist / 1e-9, V_t / 1e-3, "-b", label="Confinement potential")
ax2.tick_params(axis="y", labelcolor="blue")
ax2.axhline(V_t[0] / 1e-3, color="black", lw=0.5, ls="--", alpha=0.5)

plt.show()
fig.savefig(path_out / "confinement_potential_t0.png", dpi=150)

# Varying the gate voltages with respect to time
time_steps = 10
time = np.linspace(tau / 3, 2 * tau / 3, time_steps)

# plot the confinement potential linecuts
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Distance along linecut [nm]")
ax.set_ylabel("confinement potential ($V$) [meV]")
for idt, t in enumerate(time):
    dist, phi_t, V_t = confinement_potential(t)
    ax.plot(dist / 1e-9, V_t / 1e-3, label=f"t = {t * 1e9:.2f} ns")
ax.annotate(
    "", (13, 30), (12, 30), arrowprops=dict(facecolor="black", arrowstyle="-|>")
)
ax.text(12.5, 40, "Time Progression", fontsize=10, ha="center", va="center")
ax.legend()
plt.show()
fig.savefig(path_out / "confinement_potential_time_sweep.png", dpi=150)

# Fit a second-order polynomial near the potential's minimum.
fit_range = 5e-9
# Lists to store the fitting parameters.
coef_list = []
x_min_list = []
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Distance along linecut [nm]")
ax.set_ylabel("confinement potential ($V$) [meV]")
for idt, t in enumerate(time):
    dist, _, V_t = confinement_potential(t)
    # Find the minimum of the linecut.
    min_idx = np.argmin(V_t)
    min_dist = dist[min_idx]

    # Select the range around the minimum for fitting
    fit_mask = (dist >= min_dist - fit_range) & (dist <= min_dist + fit_range)
    fit_dist = dist[fit_mask]
    fit_potential = V_t[fit_mask]

    # Fit a second-order polynomial
    coeffs = np.polyfit(fit_dist, fit_potential, 2)
    coef_list.append(coeffs)

    # create a polynomial function from the coefficients
    poly = np.poly1d(coeffs)

    # Calculate the approximated polynomial's minimum
    x_min = -coeffs[1] / (2 * coeffs[0])
    x_min_list.append(x_min)  # save the minimum of the polynomial
    y_min = poly(x_min)

    # Plot the results
    line = ax.plot(dist / 1e-9, V_t / 1e-3, label=f"t = {t * 1e9:.2f} ns")
    ax.plot(fit_dist / 1e-9, poly(fit_dist) / 1e-3, "--", color=line[0].get_color())
    ax.plot(x_min / 1e-9, y_min / 1e-3, "ko")
ax.annotate(
    "", (13, 30), (12, 30), arrowprops=dict(facecolor="black", arrowstyle="-|>")
)
ax.text(12.5, 40, "Time Progression", fontsize=10, ha="center", va="center")
ax.legend()
plt.show()
fig.savefig(path_out / "confinement_potential_parabolic_fits.png", dpi=150)

# Extract the position, speed and size of the quantum dot.
coef_list = np.array(coef_list)
x_min_list = np.array(x_min_list)
# Shift the position to start from 0
x_min_list = x_min_list - x_min_list[0]
w = np.sqrt(2 * ct.e * coef_list[:, 0] / (0.19 * ct.me))
size = np.sqrt(ct.hbar / (0.19 * ct.me * w))  # Size of the quantum dot
speed = np.gradient(x_min_list, time)

# Plot the position of the quantum dot over time
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Time [ns]")
ax.set_ylabel("Position [nm]")
ax.fill_between(
    time * 1e9,
    (x_min_list - size / 2) / 1e-9,
    (x_min_list + size / 2) / 1e-9,
    color="gray",
    alpha=0.3,
    label="Quantum Dot Size",
)
ax.plot(time * 1e9, x_min_list / 1e-9, "o-", label="Position")
ax.legend()
plt.show()
fig.savefig(path_out / "dot_position_vs_time.png", dpi=150)

# Plot the speed of the quantum dot over time
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Time [ns]")
ax.set_ylabel("Speed [m/s]")
ax.plot(time * 1e9, speed, "o-")
plt.show()
fig.savefig(path_out / "dot_speed_vs_time.png", dpi=150)

# Plot the Size of the quantum dot over time
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("Time [ns]")
ax.set_ylabel("Size [nm]")
ax.plot(time * 1e9, size / 1e-9, "o-")
plt.show()
fig.savefig(path_out / "dot_size_vs_time.png", dpi=150)

# Save the shuttling summary (position, speed, size vs. time) for the
# analysis notebook
np.savetxt(
    path_out / "shuttling_summary.csv",
    np.column_stack([time, x_min_list, speed, size]),
    header="time (s), position (m), speed (m/s), size (m)",
)
