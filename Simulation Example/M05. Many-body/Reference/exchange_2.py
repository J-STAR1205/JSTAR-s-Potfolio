__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import pathlib
from timeit import default_timer as timer
from qtcad.device import constants as ct
from qtcad.device import io
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

# Set up paths
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# Load the perturbative value of the exchange coupling from a previous
# calculation
path_exchange_perturbation = path_out / "exchange_perturbation.txt"
exchange_perturbation = np.loadtxt(path_exchange_perturbation)

# Load the mesh
scaling = 1e-9
path_file = str(path_mesh / "refined_dqdfdsoi.msh")
mesh = Mesh(scaling, path_file)

# Define the gate biases
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias = 0.57
plunger_gate_2_bias = 0.59
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

# Create the submesh object for the dot region
submesh = SubMesh(dvc.mesh, dot_region_list)

# Load the electric potential from the previous run to speed up Poisson
# calculation
phi = io.load(path_out / "tunnel_coupling_final_potential.hdf5", var_name="phi")
dvc.set_potential(phi)

# Configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence threshold (tolerance) for the error

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson)

# Solve Poisson's equation
poisson_slv.solve(initialize=False)

# Configure the Schrödinger solver
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-6  # Tolerance on energies in eV

# Get the potential energy from the band edge for usage in the Schrodinger
# solver
dvc.set_V_from_phi()

# Create a subdevice for the dot region
subdvc = SubDevice(dvc, submesh)

# Create a Schrodinger solver, and run it
schrod_solver = SchrodingerSolver(subdvc, solver_params=params_schrod)
schrod_solver.solve()
subdvc.print_energies()

# Solving the many-body hamiltonian
# Instantiate many-body solver
manybody_solver_params = ManyBodySolverParams()
manybody_solver_params.n_degen = 2
manybody_solver_params.num_states = 8
manybody_solver_params.num_particles = [2]
manybody_solver_params.overlap = True
manybody_solver = ManyBodySolver(subdvc, solver_params=manybody_solver_params)

# Solve the many-body hamiltonian
t0 = timer()
manybody_solver.solve()
print("Two-electron solution found in {} s".format(timer() - t0))

# Get the two-electron eigenenergies
two_electron_energies = subdvc.get_N_particle_subspace(2).eigval / ct.e
print("Two-electron energies")
print(two_electron_energies[0:8])
np.savetxt(path_out / "exchange_two_electron_energies.txt", two_electron_energies)

# Calculate the exchange interaction strength
exchange = two_electron_energies[1] - two_electron_energies[0]

# Print the output, and compare it with perturbation theory
print(
    "Exchange coupling (exact diagonalization): {:.3f} MHz".format(
        exchange * ct.e / ct.h / 1e6
    )
)
print(
    "Exchange coupling (perturbation theory): {:.3f} MHz".format(
        exchange_perturbation / ct.h / 1e6
    )
)
