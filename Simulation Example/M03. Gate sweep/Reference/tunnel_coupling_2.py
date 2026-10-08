__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from qtcad.device import io
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

# Set up paths
script_dir = Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# Load the mesh
scaling = 1e-9
file = path_mesh / "refined_dqdfdsoi.msh"
mesh = Mesh(scaling, file)

# Load the detuning and coupling
path = path_out / "detuning_and_coupling.txt"
detuning, coupling = np.loadtxt(path)

# Define gate biases including the detuning found previously
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias = 0.51
plunger_gate_2_bias = 0.59 + detuning
barrier_gate_3_bias = 0.5

# Define the device
dvc = get_double_dot_fdsoi(
    mesh,
    back_gate_bias,
    barrier_gate_1_bias,
    plunger_gate_1_bias,
    barrier_gate_2_bias,
    plunger_gate_2_bias,
    barrier_gate_3_bias,
)

# List of regions forming the double quantum dot region
dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]

# Configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence threshold (tolerance) for the error
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson)

# Instantiate Schrodinger solver parameters
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-6  # Tolerance on energies in eV

# Load the electric potential from the previous run
phi = io.load(path_out / "tunnel_coupling_final_potential.hdf5", var_name="phi")
dvc.set_potential(phi)

# Create the submesh object for the dot region
submesh = SubMesh(dvc.mesh, dot_region_list)

# Define barrier gate sweep
V_barrier_vec = np.linspace(0.57, 0.51, 7)

# Initialize the energy level array
energy_arr = np.zeros((len(V_barrier_vec), 11))
energy_arr[:, 0] = V_barrier_vec

for idx, V_barrier_2 in enumerate(V_barrier_vec):
    print("=" * 80)
    print("Solving at barrier gate 2 bias: {:.3f} V".format(V_barrier_2))

    # Set barrier height
    dvc.set_applied_potential("barrier_gate_2_bnd", V_barrier_2)

    # Self-consistent solution
    poisson_slv.solve(initialize=False)

    # Calculate the confinement potential from the electric potential
    # The confinement potential is stored in dvc.V
    dvc.set_V_from_phi()

    # Create subdevice for the quantum dots
    subdvc = SubDevice(dvc, submesh)

    # Then instantiate the Schrödinger solver
    schrod_slv = SchrodingerSolver(subdvc, solver_params=params_schrod)

    # Solve Schrödinger's equation
    schrod_slv.solve()

    # Store the solution as an initial guess for the next run
    params_schrod.guess = (subdvc.eigenfunctions, subdvc.energies)

    # Save energy levels
    energy_arr[idx, 1:] = subdvc.energies / ct.e

    np.savetxt(path_out / "tunnel_coupling_energies.txt", energy_arr)

# Calculate the splitting between ground and first excited states
splittings = energy_arr[:, 2] - energy_arr[:, 1]

# Fit a straight line in the log of the splitting
fit = np.polyfit(V_barrier_vec, np.log(splittings), 1)
coeff = np.exp(fit[1])
pref = fit[0]

# Plot the tunnel splittings
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_title("Barrier-gate control of the tunnel coupling")
ax.set_xlabel("Barrier gate 2 bias (V)")
ax.set_ylabel("Tunnel splitting (eV)")
ax.set_yscale("log")
ax.plot(V_barrier_vec, splittings, "ob")
ax.plot(
    V_barrier_vec, coeff * np.exp(pref * V_barrier_vec), "--k", label="Exponential fit"
)
plt.show()

# Save the largest tunnel splitting in a text file
np.savetxt(path_out / "tunnel_coupling_low_barrier.txt", np.array([splittings[0]]))
