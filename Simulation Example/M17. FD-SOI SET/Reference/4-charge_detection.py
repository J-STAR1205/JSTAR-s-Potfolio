__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from qtcad.device import io
from qtcad.device import constants as ct
from qtcad.device.mesh3d import SubMesh
from qtcad.device import SubDevice
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from double_dot_fdsoi import (
    get_double_dot_fdsoi,
    V_set,
    save_slice,
    state_probability_density,
)
from chemical_potential import compute_chemical_potential, schrod_params

# Files -----------------------------------------------------------------------
script_dir = pathlib.Path(__file__).parent.resolve()

path_out = script_dir / "output"
path_local_out = path_out / "charge_sensing"
phi_in_file = path_out / "coulomb_peaks" / f"phi_vset{V_set:.2f}_N13.hdf5"
la_in_file = path_out / "coulomb_peaks" / "lever_arm.txt"

# Visualization
wf_qubit_slice_file = path_local_out / "wf_qubit_slice.png"
rho_qubit_slice_file = path_local_out / "rho_qubit_slice.png"
lever_arm_plot_file = path_local_out / "lever_arm_plot.svg"

# Setup ----------------------------------------------------------------------

# Create the device
d, set_region, qd_region = get_double_dot_fdsoi(mesh_file_name="refined_dqdfdsoi.msh")
# load the potential
d.set_potential(io.load(phi_in_file, var_name="phi"))
d.set_V_from_phi()

# Qubit quantum dot -----------------------------------------------------------
# QD subdevice
qd_mesh = SubMesh(d.mesh, qd_region)
qd_device = SubDevice(d, qd_mesh)
qd_device.set_insulator_boundaries()
# Solve the Schrödinger equation
schrod_solver = SchrodingerSolver(qd_device, solver_params=schrod_params)
schrod_solver.solve()
qd_device.print_energies()

# Lowest valley-resolved state:
# axis 1 -> state index, axis 2 -> valley component index
gs_state = 0
gs_valley = 0
gs_wavefunction = qd_device.eigenfunctions[:, gs_state, gs_valley]
gs_density = state_probability_density(qd_device, gs_state)

# Charge density from the QD ground state
rho0 = np.zeros_like(d.phi)
rho0[qd_device.mesh.nodes_in_parent] = -ct.e * gs_density
d.set_vol_charge_density(rho0)

# Visualize the QD ground-state wavefunction component and charge density
# (Ground-state wavefunction component)
save_slice(
    qd_device,
    gs_wavefunction,
    wf_qubit_slice_file,
    label="Ground-state component, $\\psi_{0,0}$ [1 / m$^{3/2}$]",
)
# (Charge density)
save_slice(d, d.rho_0 / ct.e, rho_qubit_slice_file, label="$\\rho_0$ [$e$ / m$^3$]")

# Compute the chemical potential ----------------------------------------------
set_mesh = SubMesh(d.mesh, set_region)
set_device = SubDevice(d, set_mesh)
set_device.set_insulator_boundaries()

# V_set values
V_SET = np.array([V_set, V_set - 0.05])
# Number of particles
Nmin = 12
Nmax = 13
# Generate data
Mu = compute_chemical_potential(
    V_SET, Nmin, Nmax, d, set_device, out_path=path_local_out
)

# Coulomb peak position -------------------------------------------------------

# Fit a line: y = alpha*x + b
alpha, b = np.polyfit(V_SET, Mu / ct.e, 1)
print("--------------------------------------------------")
print(f"Charge Sensing: Lever Arm [eV/V] = {alpha}, Intercept [eV] = {b}")
peak_pos = -b / alpha
print(f"Charge Sensing: Coulomb Peak Position [V]: {peak_pos}")

# Load lever arm data from previous Coulomb peak calculation
charge_neutral_la = np.loadtxt(la_in_file)
lever_arm = charge_neutral_la[:, 0]
intercept = charge_neutral_la[:, 1]
V_peaks = charge_neutral_la[:, 2]

# Compare peak positions
d_pos = np.abs(peak_pos[0] - V_peaks[0])
print(f"Shift in peak position [V]: {d_pos}")

# Plot
plt.figure(figsize=(8, 5))
# (Add data to plot)
V_line = np.linspace(
    np.min([V_peaks[0], peak_pos[0]]) - d_pos,
    np.max([V_peaks[0], peak_pos[0]]) + d_pos,
    200,
)
plt.axhline(y=0, color="black", linewidth=3)  # E_F = 0
plt.plot(
    V_line,
    alpha * V_line + b,
    color="blue",
    linewidth=2.0,
    linestyle="--",
    label="$Q_{\\mathrm{qubit}}=1$",
)  # Plot data for charged qubit dot

alpha_q1 = lever_arm[0]
b_q1 = intercept[0]
mu_line = alpha_q1 * V_line + b_q1
plt.plot(
    V_line,
    mu_line,
    color="red",
    linewidth=2.0,
    linestyle="--",
    label="$Q_{\\mathrm{qubit}}=0$",
)  # Plot data for empty qubit dot

# (Labels, legend, grid)
plt.xlabel("$V_{SET}$ [V]", fontsize=16)
plt.ylabel("$\\mu$ [eV]", fontsize=16)
plt.legend(fontsize=14)
plt.grid(True)
# (Save)
plt.tight_layout()
plt.savefig(lever_arm_plot_file, dpi=300)
