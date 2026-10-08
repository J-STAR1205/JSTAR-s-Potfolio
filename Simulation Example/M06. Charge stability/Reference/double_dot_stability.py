__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import pathlib
from matplotlib import pyplot as plt
from progress.bar import ChargingBar as Bar
from qtcad.device import constants as ct
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from qtcad.device.leverarm_matrix import Solver as LeverArmSolver
from qtcad.device.leverarm_matrix import SolverParams as LeverArmSolverParams
from qtcad.transport.junction import Junction
from qtcad.transport.mastereq import add_spectrum
from helper.double_dot_fdsoi import get_double_dot_fdsoi

# Set up paths
script_dir = pathlib.Path(__file__).parent.resolve()

# Mesh path
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh / "dqdfdsoi.msh"))

# Define the gate bias parameters
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.6
barrier_gate_2_bias = 0.5
plunger_gate_2_bias = 0.6 + 1e-3
barrier_gate_3_bias = 0.5

# Define the device object from the function defined in the FD-SOI tutorial
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
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 50000
params_poisson.max_nodes = 1e5
params_poisson.maxiter_adapt = 30
params_poisson.maxiter = 200
params_poisson.refined_region = dot_region_list
params_poisson.h_refined = 0.8
# NOTE: these files are also read back in by gmsh during adaptive refinement
# (same Korean-path XAO/gmsh bug as geo_file above), so route them through
# the ASCII-only temp folder as well.
params_poisson.size_map_filename = r"C:\temp\qtcad_dqdfdsoi\refined_dqdfdsoi.pos"
params_poisson.refined_mesh_filename = r"C:\temp\qtcad_dqdfdsoi\refined_dqdfdsoi.msh"

# Number of single-electron orbital states to be considered
num_states = 4

# Instantiate Schrodinger solver's parameters
params_schrod = SchrodingerSolverParams()
params_schrod.num_states = num_states  # Number of states to consider

# Instantiate Poisson solver
# NOTE: geo_file is loaded internally by gmsh's Python API (gmsh.merge) during
# adaptive mesh refinement. gmsh's XAO/CAD reader fails with
# "Could not load XML file" when the path contains non-ASCII (Korean)
# characters, which this project's folder path does. As a workaround, a copy
# of dqdfdsoi.xao was placed at an ASCII-only path and is referenced here
# instead of the original (Korean-path) location.
geo_file_ascii = r"C:\temp\qtcad_dqdfdsoi\dqdfdsoi.xao"
poisson_slv = PoissonSolver(
    dvc, solver_params=params_poisson, geo_file=geo_file_ascii
)

# Solve Poisson's equation
poisson_slv.solve()

# Get the potential energy from the band edge for usage in the Schrodinger
# solver
dvc.set_V_from_phi()

# Create a submesh including only the dot region
submesh = SubMesh(dvc.mesh, dot_region_list)

# Create a subdevice for the dot region
subdvc = SubDevice(dvc, submesh)

# Create a Schrodinger solver
schrod_solver = SchrodingerSolver(subdvc, solver_params=params_schrod)

# Solve Schrodinger's equation
schrod_solver.solve()

# Output and save single-particle energy levels
subdvc.print_energies()
energies = subdvc.energies
np.save(path_out / "dqdfdsoi_energies.npy", energies)

# Plot single-electron eigenstates
x_slice = 0
z_slice = -2e-9
an.plot_slices(
    subdvc.mesh, subdvc.eigenfunctions[:, 0], title="Ground state", x=x_slice, z=z_slice
)
an.plot_slices(
    subdvc.mesh,
    subdvc.eigenfunctions[:, 1],
    title="First excited state",
    x=x_slice,
    z=z_slice,
)
an.plot_slices(
    subdvc.mesh,
    subdvc.eigenfunctions[:, 2],
    title="Second excited state",
    x=x_slice,
    z=z_slice,
)
an.plot_slices(
    subdvc.mesh,
    subdvc.eigenfunctions[:, 3],
    title="Third excited state",
    x=x_slice,
    z=z_slice,
)

# Instantiate the voltage bias vector
bias_vector = np.array([plunger_gate_1_bias, plunger_gate_2_bias])

# Bias labels
gate_labels = ["plunger_gate_1_bnd", "plunger_gate_2_bnd"]

# Solver params for the LeverArmSolver
lam_params = LeverArmSolverParams()
lam_params.pot_solver_params = params_poisson
lam_params.schrod_solver_params = params_schrod

# Instantiate lever arm matrix solver
slv = LeverArmSolver(
    dvc, gate_labels, bias_vector, dot_region=dot_region_list, solver_params=lam_params
)

# Calculate the lever arm matrix
bias_increment = 1e-3
lever_arm_matrix = slv.solve(bias_increment=bias_increment)

# Print and save the lever arm matrix
print("Lever arm matrix")
print(lever_arm_matrix)
np.save(path_out / "lever_arm_matrix.npy", lever_arm_matrix)

# Calculate the Coulomb interaction matrix
# Instantiate many-body solver
many_body_solver_params = ManyBodySolverParams()
many_body_solver_params.n_degen = 2
many_body_solver_params.num_states = num_states
slv = ManyBodySolver(subdvc, solver_params=many_body_solver_params)

# Compute, save, and print Coulomb interaction matrix without overlap terms
coulomb_no_overlap = slv.get_coulomb_matrix(overlap=False, verbose=True)
np.save(path_out / "coulomb_mat_no_overlap.npy", coulomb_no_overlap)
print("Coulomb interaction matrix (eV)")
print(coulomb_no_overlap / ct.e)

# Compute and save Coulomb interaction matrix with overlap terms
# coulomb_overlap = slv.get_coulomb_matrix(verbose=True, overlap=True)
# np.save(path_out/"full_coulomb_mat.npy", coulomb_overlap)

# Set the junction temperature to 10 K (instead of 100 mK) to make lines
# in the charge stability diagram thicker and decrease the resolution
# required to observe them
temperature_spec = 10

# Contact labels
gate_labels = ["source_bnd", "plunger_gate_1_bnd", "plunger_gate_2_bnd", "drain_bnd"]

# Instantiate a junction
many_body_solver_params.energies = energies
many_body_solver_params.overlap = False
many_body_solver_params.coulomb_mat = coulomb_no_overlap
many_body_solver_params.alpha = lever_arm_matrix
jc = Junction(
    many_body_solver_params=many_body_solver_params,
    temperature=temperature_spec,
    contact_labels=gate_labels,
)


def get_add_spectrum(junc, gate_labels, gate_biases, temperature, verbose=True):
    """Computes the response function of a double quantum dot for a given
    set of gate biases.

    Args:
       junc (transport junction object) : Junction object to consider.
       gate_labels (list): Labels of the gates at which biases are applied.
       gate_biases (ndarray) : Applied bias on each gate.
       temperature (float): Temperature to use in the particle addition
          spectrum.
       verbose (bool, optional): Verbosity level.

    Returns:
       float: The addition spectrum at the chosen gate bias configuration.
    """
    # Set list of biases with drain potential + increment
    junc.set_biases(gate_labels, gate_biases, verbose=verbose)

    # Calculate response function
    out = add_spectrum(junc, temperature=temperature)
    return out


# Base source and drain potentials
source_potential = 0
drain_potential = 0

# Extremal values of the gate biases
min_gate_1_bias = 30e-3
max_gate_1_bias = 140e-3
min_gate_2_bias = 30e-3
max_gate_2_bias = 140e-3

gate_1_biases = np.linspace(min_gate_1_bias, max_gate_1_bias, 90)
gate_2_biases = np.linspace(min_gate_2_bias, max_gate_2_bias, 90)
add_spectrum_mat = np.zeros((len(gate_1_biases), len(gate_2_biases)))

bartitle = "Calculating charge stability diagram."
progress_bar = Bar(bartitle, max=len(gate_1_biases) * len(gate_1_biases))
for idx_1, gate_1_bias in enumerate(gate_1_biases):
    for idx_2, gate_2_bias in enumerate(gate_2_biases):
        add_spectrum_mat[idx_1, idx_2] = get_add_spectrum(
            jc,
            gate_labels,
            np.array([source_potential, gate_1_bias, gate_2_bias, drain_potential]),
            temperature_spec,
            verbose=False,
        )
        np.savetxt(path_out / "addition_spectrum.txt", add_spectrum_mat)
        progress_bar.next()

progress_bar.finish()

# Plot the particle addition spectrum
# Transpose and flip the matrix to turn into an image with bias_1 along x
# and bias_2 along y
add_spec_to_plot = np.flip(np.transpose(add_spectrum_mat), axis=0)
fig, axs = plt.subplots()
axs.set_xlabel("$V_{g1}$ (V)", fontsize=16)
axs.set_ylabel("$V_{g2}$ (V)", fontsize=16)
vg1_ref = plunger_gate_1_bias
vg2_ref = plunger_gate_2_bias
diff_conds_map = axs.imshow(
    add_spec_to_plot / np.max(add_spec_to_plot),
    cmap="jet",
    interpolation="bilinear",
    extent=[
        min_gate_1_bias + vg1_ref,
        max_gate_1_bias + vg1_ref,
        min_gate_2_bias + vg2_ref,
        max_gate_2_bias + vg2_ref,
    ],
    aspect="auto",
)
fig.colorbar(diff_conds_map, ax=axs, label="Response (arb. units)")
fig.tight_layout()
plt.show()
