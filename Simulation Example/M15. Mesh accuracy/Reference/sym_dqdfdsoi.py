__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import io
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device
from qtcad.device import SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams

# Set up paths
script_dir = Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"
path_out.mkdir(exist_ok=True)

# NOTE: gmsh's XAO reader/writer fails on non-ASCII (Korean) paths with
# "Could not load XML file" (same bug hit in M12/M13/M14). The project
# directory is under a Korean-named folder, so geo_file reads/writes for
# symmetric_mesh() and the adaptive Poisson solver are staged through an
# ASCII-only path in C:\temp instead, then copied back to the project folder.
path_mesh_ascii = Path(r"C:\temp\m15_sym_dqdfdsoi")
path_mesh_ascii.mkdir(parents=True, exist_ok=True)
path_geo = path_mesh_ascii / "qdfdsoi.xao"

scaling = 1e-9
file = path_mesh_ascii / "qdfdsoi.msh"
hmesh = Mesh(scaling, file)
# NOTE: Mesh.show()/mesh.show() launch an interactive QTCAD/pyvista browser
# viewer (open_browser=True by default), which blocks headless script
# execution indefinitely. Disabled here for the same reason Builder.view()
# calls were disabled in M12/M13; the device geometry is instead
# reconstructed independently from the .geo parameters in the analysis
# notebook.
# hmesh.show()

mesh = hmesh.symmetric_mesh(
    0,
    1,
    0,
    -32.5e-9,
    mesh_file=str(path_mesh_ascii / "qdfdsoi_mirrored.msh"),
    geo_file=str(path_geo),
    verbose=True,
)
# mesh.show()

# Define gate biases
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias = 0.51

dvc = Device(hmesh, conf_carriers="e")
dvc.set_temperature(0.1)

# Create the regions
dvc.new_region("oxide", mt.SiO2)
dvc.new_region("oxide_dot", mt.SiO2)
dvc.new_region("gate_oxide", mt.HfO2)
dvc.new_region("gate_oxide_dot", mt.HfO2)
dvc.new_region("buried_oxide", mt.SiO2)
dvc.new_region("buried_oxide_dot", mt.SiO2)
dvc.new_region("channel", mt.Si)
dvc.new_region("channel_dot", mt.Si)
dvc.new_region("source", mt.Si, ndoping=1e20 * 1e6)
# (mirrored regions are omitted)

# Boundary conditions
Ew = mt.Si.Eg / 2 + mt.Si.chi

dvc.new_gate_bnd("barrier_gate_1_bnd", barrier_gate_1_bias, Ew)
dvc.new_gate_bnd("plunger_gate_1_bnd", plunger_gate_1_bias, Ew)
dvc.new_gate_bnd("barrier_gate_2_bnd", barrier_gate_2_bias, Ew)
dvc.new_ohmic_bnd("source_bnd")
dvc.new_frozen_bnd(
    "back_gate_bnd", back_gate_bias, mt.Si, 1e15 * 1e6, "n", 46 * 1e-3 * ct.e
)

# Mark quantum dot region
dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]
dvc.set_dot_region(dot_region_list)

# Configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-4
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 60000
params_poisson.max_nodes = 1e5
params_poisson.maxiter_adapt = 30
params_poisson.refined_region = dot_region_list
params_poisson.h_refined = 0.7
params_poisson.mirror_plane_coefs = (0, 1, 0, -32.5e-9)
params_poisson.refined_mesh_filename = path_mesh_ascii / "sym_refined_dqdfdsoi.msh"

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=path_geo)

# Solve Poisson equation
poisson_slv.solve()

# Produce a linecut of the conduction band edge along the channel
x, y, z = dvc.mesh.glob_nodes.T
ymin = np.min(y)
ymax = np.max(y)
distance, linecut_qdot = an.linecut(
    dvc.mesh, dvc.cond_band_edge(), (0, ymin, -1e-9), (0, ymax, -1e-9)
)

# Plot the linecuts
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("Distance along linecut (nm)")
ax.set_ylabel("Conduction band edge, $E_C$ (eV)")
ax.plot(distance / 1e-9, linecut_qdot / ct.e)
fig.savefig(path_out / "sym_cbe.png")
plt.show()

# Define the mirrored dot regions
dot_region_mirrored_list = [
    "oxide_dot",
    "gate_oxide_dot",
    "buried_oxide_dot",
    "channel_dot",
    "oxide_dot_mirrored",
    "gate_oxide_dot_mirrored",
    "buried_oxide_dot_mirrored",
    "channel_dot_mirrored",
]

# Create a submesh including only the dot region
submesh = SubMesh(dvc.mesh, dot_region_mirrored_list)

# Instantiate Schrodinger solver parameters
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-7  # Tolerance on energies in eV

# Calculating the confinement potential
dvc.set_V_from_phi()

# Create a subdevice for the dot region
subdvc = SubDevice(dvc, submesh)

# Instantiate and run a Schrödinger solver
schrodinger_slv = SchrodingerSolver(subdvc, solver_params=params_schrod)
schrodinger_slv.solve()

# Plot the ground state wavefunction
an.plot_slices(
    submesh,
    subdvc.eigenfunctions[:, 0],
    z=-2e-9,
    title="Ground state wavefunction",
    show_figure=False,
    path=path_out / "sym_ground_state.png",
)

# Linecut for ground and excited state
distance, linecut_ground = an.linecut(
    submesh, subdvc.eigenfunctions[:, 0], (0, ymin, -2e-9), (0, ymax, -2e-9)
)
distance, linecut_excited = an.linecut(
    submesh, subdvc.eigenfunctions[:, 1], (0, ymin, -2e-9), (0, ymax, -2e-9)
)

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("Distance along linecut (nm)")
ax.set_ylabel("Eigenfunctions (arb. units)")
ax.plot(distance / 1e-9, linecut_ground, "ro-", label=f"Ground state wavefunction")
ax.plot(
    distance / 1e-9, linecut_excited, "bx-", label=f"First excited state wavefunction"
)
ax.legend()
fig.savefig(path_out / "sym_eigenfunctions.png")
plt.show()

# Tunnel coupling calculation
energies = subdvc.energies
tunnel_coupling = np.abs(energies[1] - energies[0]) / ct.e
print(f"Tunnel coupling: {tunnel_coupling} eV")
