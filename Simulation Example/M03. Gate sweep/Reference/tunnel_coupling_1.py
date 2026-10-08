__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.optimize import minimize_scalar
from timeit import default_timer as timer
from qtcad.device import constants as ct
from qtcad.device import io
from qtcad.device import analysis as an
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
# NOTE: path_geo is intentionally pointed at an ASCII-only path outside the
# project directory. The project directory contains non-ASCII (Korean)
# characters, and gmsh's XAO/XML reader (invoked internally by QTCAD's
# adaptive mesh refinement) fails to open .xao files when the path contains
# such characters, raising "Could not load XML file". The .xao content
# itself is identical; a plain-ASCII copy is used purely to work around
# this gmsh limitation.
path_geo = Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo\dqdfdsoi.xao"
)
path_out = script_dir / "output"

# Load the mesh
scaling = 1e-9
file = path_mesh / "dqdfdsoi.msh"
mesh = Mesh(scaling, file)

# Define gate biases
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias_low = 0.51
barrier_gate_2_bias_high = 0.57
plunger_gate_2_bias = 0.59
barrier_gate_3_bias = 0.5

# Define the device
dvc = get_double_dot_fdsoi(
    mesh,
    back_gate_bias,
    barrier_gate_1_bias,
    plunger_gate_1_bias,
    barrier_gate_2_bias_high,
    plunger_gate_2_bias,
    barrier_gate_3_bias,
)

# List of regions forming the double quantum dot region
dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]

# Configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence threshold (tolerance) for the error
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 50000
params_poisson.max_nodes = 1e5
params_poisson.maxiter_adapt = 30
params_poisson.maxiter = 200
params_poisson.refined_region = dot_region_list
params_poisson.h_refined = 0.8
params_poisson.refined_mesh_filename = path_mesh / "refined_dqdfdsoi.msh"

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=path_geo)

# Solve in the original gate bias configuration
poisson_slv.solve()

# Produce a linecut of the conduction band edge along the channel
x, y, z = dvc.mesh.glob_nodes.T
ymin = np.min(y)
ymax = np.max(y)
distance, linecut_high_bias = an.linecut(
    dvc.mesh, dvc.cond_band_edge(), (0, ymin, -1e-9), (0, ymax, -1e-9)
)

# Solve Poisson again at a lower central barrier gate bias
dvc.set_applied_potential("barrier_gate_2_bnd", barrier_gate_2_bias_low)
poisson_slv.solve(initialize=False)

# Produce another linecut of the conduction band edge along the channel
distance, linecut_low_bias = an.linecut(
    dvc.mesh, dvc.cond_band_edge(), (0, ymin, -1e-9), (0, ymax, -1e-9)
)

# Plot the linecuts
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("Distance along linecut (nm)")
ax.set_ylabel("Conduction band edge, $E_C$ (eV)")
ax.plot(
    distance / 1e-9,
    linecut_high_bias / ct.e,
    "-r",
    label=f"Barrier gate 2 bias: {barrier_gate_2_bias_high} V",
)
ax.plot(
    distance / 1e-9,
    linecut_low_bias / ct.e,
    "-b",
    label=f"Barrier gate 2 bias: {barrier_gate_2_bias_low} V",
)
ax.legend()
plt.show()

# Create a submesh including only the dot region
submesh = SubMesh(dvc.mesh, dot_region_list)

# Instantiate Schrodinger solver parameters
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-6  # Tolerance on energies in eV


# defining the function for optimization.
def function(x, gate_bias, path_states=None, path_phi=None):
    """Tunnel splitting to be optimized.

    Args:
        x (float): Gate bias detuning with respect to reference configuration.
        gate_bias (float): Gate bias in the reference configuration.
        path_states (str): Path in which to save the ground and first
            excited state wavefunctions in .vtu format. If None, do not save.
        path_phi (str): Path in which to save the electric potential in .hdf5
            format.

    Returns:
        float: The tunnel splitting (in eV).

    """

    print("=" * 80)
    print("Evaluating double-dot splitting for detuning: {} V".format(x))
    print("=" * 80)

    # Set the plunger gate bias
    dvc.set_applied_potential("plunger_gate_2_bnd", gate_bias + x)

    # Solve the non-linear Poisson equation
    poisson_slv.solve(initialize=False)

    # Get the potential energy from the band edge for usage in the Schrodinger
    # solver
    dvc.set_V_from_phi()

    # Create a subdevice for the dot region
    subdvc = SubDevice(dvc, submesh)

    # Instantiate and run a Schrödinger solver
    schrodinger_slv = SchrodingerSolver(subdvc, solver_params=params_schrod)
    schrodinger_slv.solve()

    # Store the solution as an initial guess for the next run
    params_schrod.guess = (subdvc.eigenfunctions, subdvc.energies)

    energies = subdvc.energies

    if path_states is not None:
        out_dict = {
            "Ground state": subdvc.eigenfunctions[:, 0],
            "First excited state": subdvc.eigenfunctions[:, 1],
        }
        io.save(path_states, out_dict, submesh)

    if path_phi is not None:
        out_dict = {"phi": dvc.phi}
        io.save(path_phi, out_dict)

    return (energies[1] - energies[0]) / ct.e


# Tune the bias of plunger gate 2 to find a symmetric double dot
# configuration
t0 = timer()
result = minimize_scalar(
    function,
    bounds=(-200e-6, 200e-6),
    args=(np.array([plunger_gate_2_bias])),
    method="bounded",
    options={"xatol": 1e-7},
)

# print results:
print(f"Optimum found at {result.x} V with the value: {result.fun} eV.")
print(f"Computation time: {timer() - t0} s")

# Save the optimum and detuning in a file
path_result = path_out / "detuning_and_coupling.txt"
tab = np.array([result.x, result.fun])
np.savetxt(path_result, tab)

# Save the ground and first excited state wave functions in .vtu format
# Save the electric potential in .hdf5 format
print("Solving Poisson and Schrödinger for the final detuning value.")
print("Saving eigenfunctions in .vtu format.")
path_states = str(path_out / "tunnel_coupling_wavefunctions_high_barrier.vtu")
path_phi = str(path_out / "tunnel_coupling_final_potential.hdf5")
function(
    result.x, np.array(plunger_gate_2_bias), path_states=path_states, path_phi=path_phi
)

# Save the ground and first excited state wave functions in .vtu format
# at high central barrier gate bias (strong tunneling)
dvc.set_applied_potential("barrier_gate_2_bnd", barrier_gate_2_bias_high)
path_states = str(path_out / "tunnel_coupling_wavefunctions_low_barrier.vtu")
function(result.x, np.array(plunger_gate_2_bias), path_states=path_states)
