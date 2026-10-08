__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import pathlib
from qtcad.device import constants as ct
from qtcad.device import io
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import SubDevice
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

# Set up paths
script_dir = pathlib.Path(__file__).parent.resolve()
# mesh path
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# Load the tunnel coupling from strong-tunneling simulation
path_tunnel = path_out / "tunnel_coupling_low_barrier.txt"
tunnel_coupling = np.loadtxt(path_tunnel)

# Load the mesh
scaling = 1e-9
path_file = str(path_mesh / "refined_dqdfdsoi.msh")
mesh = Mesh(scaling, path_file)

# Define the gate biases
detuning = 5e-3
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias = 0.57
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

# Create the submesh object for the dot region
submesh = SubMesh(dvc.mesh, dot_region_list)

# Load the electric potential from the previous run to speed up Poisson
# calculation
phi = io.load(path_out / "tunnel_coupling_final_potential.hdf5", var_name="phi")
dvc.set_potential(phi)

# Configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence threshold (tolerance) for the error

# Configure the Schrödinger solver
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-6  # Tolerance on energies in eV

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson)

# Solve Poisson's equation
poisson_slv.solve(initialize=False)

# Save the potential
path_potential = path_out / "detuned_potential.hdf5"
io.save(path_potential, dvc.phi, var_name="phi")

# Get the potential energy from the band edge for usage in the Schrodinger
# solver
dvc.set_V_from_phi()

# Create a subdevice for the dot region
subdvc = SubDevice(dvc, submesh)

# Create a Schrodinger solver, and run it
schrod_solver = SchrodingerSolver(subdvc, solver_params=params_schrod)
schrod_solver.solve()
subdvc.print_energies()

# Save the ground and first excited state wavefunctions into .vtu files
path_states = path_out / "exchange_localized_wavefunctions.vtu"
out_dict = {
    "Ground state": subdvc.eigenfunctions[:, 0],
    "First excited state": subdvc.eigenfunctions[:, 1],
}
io.save(path_states, out_dict, submesh)

# Instantiate a many-body solver
manybody_solver_params = ManyBodySolverParams()
manybody_solver_params.num_states = 2
manybody_solver_params.overlap = False
manybody_solver = ManyBodySolver(subdvc, solver_params=manybody_solver_params)

# Compute the matrix of Coulomb integrals
coulomb_mat = manybody_solver.get_coulomb_matrix(overlap=False)
print("Coulomb interaction matrix in a nearly localized single-electron basis set")
print(coulomb_mat[0:2, 0:2] / ct.e)

# Export the coulomb matrix
np.savetxt(path_out / "exchange_coulomb_localized.txt", coulomb_mat)

# Calculate the theoretical exchange coupling
exchange = (tunnel_coupling * ct.e) ** 2 / coulomb_mat[0, 0]
print("Theoretical exchange interaction: {:.3f} MHz".format(exchange / ct.h / 1e6))

# Export the theoretical exchange coupling
np.savetxt(path_out / "exchange_perturbation.txt", np.array([exchange]))
