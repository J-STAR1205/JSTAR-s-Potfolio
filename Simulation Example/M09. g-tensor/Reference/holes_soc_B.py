__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from progress.bar import ChargingBar as Bar
import pathlib
from qtcad.device import materials as mt
from qtcad.device import constants as ct
from qtcad.device.pauli import J_x, J_y
from qtcad.device import Device
from qtcad.device.schrodinger import Solver, SolverParams
from qtcad.device.mesh3d import Mesh

# -------------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------------
# Constants
nm = 1e-9
A = 1e-10
meV = ct.e / 1000

# Names of files ---------------------------------------------------------------
path = pathlib.Path(__file__).parent.resolve()

# Output
n_no_soc = "hole_energies_vs_B"
n_soc = "hole_energies_withSOC_vs_B"
n_kz2 = "kz2_holes"
n_plot = "holes_energies_vs_B"

# Input
n_mesh = "box.msh"

# Initialize the meshes --------------------------------------------------------
path_mesh = str(path / "meshes" / f"{n_mesh}")
# Load the mesh
scaling = nm
mesh = Mesh(scaling, path_mesh)

# Create device ------------------------------------------------------------
# Create device from mesh
d = Device(mesh, conf_carriers="h", hole_kp_model="luttinger_kohn_foreman")

# Create regions
d.new_region("domain", mt.GaAs)

# Then create insulator boundaries
d.new_insulator("bnd")

# -------------------------------------------------------------------------------
# Solving Schrodinger Equation with a magnetic field
# -------------------------------------------------------------------------------
print("-------------------------------------------------------")
print("Solving Schrodinger Equation including B field")
print("-------------------------------------------------------")
# Set up solver parameters -------------------------------------------------
# Configure and create the Schrodinger solver
params_schrod = SolverParams()
params_schrod.num_states = 10
params_schrod.verbose = False


# Set up B field ---------------------------------------------------------------
B = 10
B_vec = np.linspace(0, 1, 6) * B
# use a progress bar
progress_bar = Bar("B loop", max=len(B_vec))

data = np.zeros((len(B_vec), params_schrod.num_states + 1))

for i, B0 in enumerate(B_vec):
    # include B in the calculation
    d.set_Bfield(np.array([0, 0, B0]), gauge="landau")

    # Solve --------------------------------------------------------------------
    # Create Schrodinger solver
    ss = Solver(d, solver_params=params_schrod)
    ss.solve()

    # Print results
    print("\n-----------------------------")
    print(f"B = {B0}")
    d.print_energies()
    # Update progress bar
    progress_bar.next()
    print("\n-----------------------------\n")

    # Save results to file
    data[i, 0] = B0
    data[i, 1:] = d.energies / ct.e
    with open(str(path / "output" / f"{n_no_soc}.txt"), "w") as f:
        np.savetxt(f, data)

# -------------------------------------------------------------------------------
# Solving Schrodinger Equation with a magnetic field and SOC
# -------------------------------------------------------------------------------
print("-------------------------------------------------------")
print("Solving Schrodinger Equation including B field and SOC")
print("-------------------------------------------------------")

# Set up SOC -------------------------------------------------------------------
# Compute <k_z^2>
kz2 = np.array([d.get_k_squared(i) for i in range(params_schrod.num_states)])

with open(str(path / "output" / f"{n_kz2}.txt"), "w") as f:
    np.savetxt(f, kz2)

# Coupling coefficient
beta = 82 * ct.e * A**3
# Add SOC to the device
d.set_soc([beta * kz2[0] * J_x, -beta * kz2[0] * J_y, None])

# use a progress bar
progress_bar = Bar("B loop", max=len(B_vec))

data = np.zeros((len(B_vec), params_schrod.num_states + 1))

for i, B0 in enumerate(B_vec):
    # include B in the calculation
    d.set_Bfield(np.array([0, 0, B0]), gauge="landau")

    # Solve --------------------------------------------------------------------
    # Create Schrodinger solver
    ss = Solver(d, solver_params=params_schrod)
    ss.solve()

    # Print results
    print("\n-----------------------------")
    print(f"B = {B0}")
    d.print_energies()
    # Update progress bar
    progress_bar.next()
    print("\n-----------------------------\n")

    # Save results to file
    data[i, 0] = B0
    data[i, 1:] = d.energies / ct.e
    with open(str(path / "output" / f"{n_soc}.txt"), "w") as f:
        np.savetxt(f, data)

# -------------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------------
# <k_z^2> ----------------------------------------------------------------------
kz2_data = np.loadtxt(str(path / "output" / f"{n_kz2}.txt"))

# Analytic value of <k_z^2>
Lz = 3 * nm
kz2_analytic = (np.pi / Lz) ** 2

# Plot
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(kz2_data * nm**2, "o", label="Data")
ax1.axhline(np.average(kz2_data) * nm**2, ls="--", label="Average")
ax1.axhline(kz2_analytic * nm**2, ls="-.", color="orange", label="Analytic")

ax1.set_xlabel(r"$state$", fontsize=16)
ax1.set_ylabel(r"$\left<k_z^2\right>$ ($\mathrm{nm}^{-2}$)", fontsize=16)
ax1.legend()

# Save plot
plt.tight_layout()
plt.savefig(str(path / "output" / f"{n_kz2}.png"))

# Energies ---------------------------------------------------------------------

# Load data
SOC_data = np.loadtxt(str(path / "output" / f"{n_soc}.txt"))
no_SOC_data = np.loadtxt(str(path / "output" / f"{n_no_soc}.txt"))

# Plot
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(SOC_data[:, 0], SOC_data[:, 1:] * ct.e / meV, marker="o", linestyle="-")
ax1.plot(no_SOC_data[:, 0], no_SOC_data[:, 1:] * ct.e / meV, "--")

ax1.set_xlabel(r"$B_0 (T)$", fontsize=16)
ax1.set_ylabel(r"$E$ ($meV$)", fontsize=16)

# Creating legend
dashed = mlines.Line2D([], [], color="black", ls="--", label="Without SOC")
circle = mlines.Line2D(
    [], [], color="black", marker="o", linestyle="None", markersize=10, label="With SOC"
)
ax1.legend(handles=[dashed, circle])

# Save plot
plt.tight_layout()
plt.savefig(str(path / "output" / f"{n_plot}.png"))
