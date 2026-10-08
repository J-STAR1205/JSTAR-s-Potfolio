__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import pathlib
from progress.bar import ChargingBar as Bar
from qtcad.device import materials as mt
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device.schrodinger import Solver, SolverParams

# -------------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------------
# Constants --------------------------------------------------------------------
nm = 1e-9
meV = ct.e / 1000

# Names of files ---------------------------------------------------------------
path = pathlib.Path(__file__).parent.resolve()

# Output
n_bw = "holes_strain_bw"
n_energy = "holes_strain_energies"
n_plot = "holes_strain"

# Input
n_mesh = "box_order2.msh"

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
# Solving Schrodinger Equation under biaxial tensile strain
# -------------------------------------------------------------------------------
# Set up solver parameters -----------------------------------------------------
params = SolverParams()
params.num_states = 10
params.verbose = False

# Set up strain ----------------------------------------------------------------
# Biaxial tensile strain
R = np.linspace(0, 0.007, 5)
strain = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
# use a progress bar
progress_bar = Bar("Strain loop", max=len(R))

nrg = []
bw = []

for r in R:
    # Set strain
    d.set_strain(r * strain)

    # Solve --------------------------------------------------------------------
    # Create Schrodinger solver
    s = Solver(d, solver_params=params)
    s.solve()

    # Print results
    print("------------------------")
    print(f"e_xx = e_yy = {r}")
    d.print_energies()
    progress_bar.next()  # Update progress bar
    print("\n------------------------")

    # Store energies
    nrg.append(d.energies)
    # Store band weights
    bw.append(d.band_weight())

# Finish progress bar
progress_bar.finish()

# Save data --------------------------------------------------------------------
# Energy
nrg = np.array(nrg)
E = np.zeros((nrg.shape[0], nrg.shape[1] + 1))
E[:, 0] = R
E[:, 1:] = nrg
np.savetxt(str(path / "data" / f"{n_energy}.txt"), E)
# Band weights
bw = np.array(bw)
BW = np.zeros((bw.shape[0], bw.shape[1], 2))
BW[..., 0] = bw[..., 0] + bw[..., 1]  # Heavy hole weight
BW[..., 1] = bw[..., 2] + bw[..., 3]  # Light hole weight
np.savetxt(str(path / "data" / f"{n_bw}.txt"), BW[:, 0, :])

# -------------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------------
# Create plot
fig = plt.figure()
# Plot energies
ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(E[:, 0], E[:, 1:] / meV)
ax1.set_xlabel(r"biaxial strain, $\varepsilon_{xx} = \varepsilon_{yy}$", fontsize=14)
ax1.set_ylabel(r"$E$ ($meV$)", fontsize=14)
# Plot band weights
ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(E[:, 0], BW[:, 0, 0], label="HH weight")
ax2.plot(E[:, 0], BW[:, 0, 1], label="LH weight")
ax2.set_xlabel(r"biaxial strain, $\varepsilon_{xx} = \varepsilon_{yy}$", fontsize=14)
ax2.set_ylabel(r"ground-state weights", fontsize=14)
ax2.legend()

# Save plot
plt.tight_layout()
plt.savefig(str(path / "output" / f"{n_plot}.png"))
