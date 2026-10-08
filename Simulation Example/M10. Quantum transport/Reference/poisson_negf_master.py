__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from pathlib import Path
from progress.bar import ChargingBar as Bar
from matplotlib import pyplot as plt
from copy import deepcopy
from scipy.interpolate import griddata
from scipy.signal import find_peaks
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device, SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.transport.negf_poisson import Solver as NEGFSolver
from qtcad.transport.negf_poisson import SolverParams as NEGFSolverParams
from qtcad.device.leverarm_matrix import Solver as LeverArmMatrixSolver
from qtcad.device.leverarm_matrix import SolverParams as LeverArmMatrixSolverParams
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from qtcad.transport.junction import Junction
from qtcad.transport.mastereq import seq_tunnel_curr

# Set up paths
script_dir = Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "fdsoi_poisson_negf_master.msh"
path_mesh_refined = script_dir / "meshes" / "refined_fdsoi.msh"
# NOTE: path_geo is intentionally pointed at an ASCII-only path outside the
# project directory. The project directory contains non-ASCII (Korean)
# characters, and gmsh's XAO/XML reader (invoked internally by QTCAD's
# adaptive mesh refinement) fails to open .xao files when the path contains
# such characters, raising "Could not load XML file". The .xao content
# itself is identical; a plain-ASCII copy is used purely to work around
# this gmsh limitation.
path_geo = Path(
    r"C:\temp\qtcad_fdsoi_poisson_negf_master\fdsoi_poisson_negf_master.xao"
)
path_out = script_dir / "output"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, path_mesh)

# Define the gate parameters
back_gate_bias = -0.2
doping_back_gate = 1e15 * 1e6
binding_energy = 46 * 1e-3 * ct.e
barrier_gate_1 = 0.6
plunger_gate = 0.707
barrier_gate_2 = 0.6
ndoping = 2e20 * 1e6
Ew = mt.Si.Eg / 2 + mt.Si.chi  # Midgap

# Define the device object
dvc = Device(mesh, conf_carriers="e")
dvc.set_temperature(0.1)

# Create the regions
dvc.new_region("top_ox", mt.SiO2_ideal)
dvc.new_region("source", mt.Si, ndoping=ndoping)
dvc.new_region("channel", mt.Si)
dvc.new_region("dot", mt.Si)
dvc.new_region("drain", mt.Si, ndoping=ndoping)
dvc.new_region("side_bot_ox", mt.SiO2_ideal)

# defining the region forming the quantum dot and the channel.
dvc.set_dot_region("dot")
channel_region_list = ["source", "channel", "dot", "drain"]

# define the boundaries
dvc.new_gate_bnd("barrier_gate_1_bnd", barrier_gate_1, Ew)
dvc.new_gate_bnd("plunger_gate_bnd", plunger_gate, Ew)
dvc.new_gate_bnd("barrier_gate_2_bnd", barrier_gate_2, Ew)
dvc.new_frozen_bnd(
    "back_gate_bnd", back_gate_bias, mt.Si, doping_back_gate, "n", binding_energy
)

# recreate the device with the extruded mesh for the channel.
dvc_negf = deepcopy(dvc)

# setup the Poisson parameters
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3  # Convergence threshold (tolerance) for the error
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 5e4
params_poisson.max_nodes = 1e5
params_poisson.maxiter_adapt = 100
params_poisson.maxiter = 2000
params_poisson.refined_region = channel_region_list
params_poisson.h_refined = 0.15
params_poisson.mesh_from_file = True
params_poisson.refined_mesh_filename = path_mesh_refined

# Configure the Schrodinger solver
params_schrodinger = SchrodingerSolverParams()
params_schrodinger.tol = 1e-9
params_schrodinger.num_states = 4

# Configure the NEGF-Poisson solver
solver_params = NEGFSolverParams()

# Define parameters pertaining to linecuts over which various quantities will
# be plotted. This line will cross the dot inside the channel.
t_top_ox = 2 * scaling  # top oxide thickness
t_si = 6 * scaling  # silicon film thickness
begin = (0, np.amin(mesh.glob_nodes[:, 1]), -t_top_ox - t_si / 5)
end = (0, np.amax(mesh.glob_nodes[:, 1]), -t_top_ox - t_si / 5)
show_figure = True  # whether plots are shown

# Solver params for the LeverArmSolver
lam_params = LeverArmMatrixSolverParams()
lam_params.pot_solver_params = params_poisson
lam_params.schrod_solver_params = params_schrodinger

# Instantiate Poisson solver
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=path_geo)

# Instantiate Poisson solver
poisson_slv.solve()

# Schrodinger equation calculations.
dvc.set_V_from_phi()

# Creating the dot region subdevice.
submesh = SubMesh(dvc.mesh, ["dot"])
subdevice = SubDevice(dvc, submesh)

# Instantiate Schrodinger solver
schrod_solver = SchrodingerSolver(subdevice, solver_params=params_schrodinger)

# Solve Schrodinger's equation
schrod_solver.solve()

# define the quantum lead boundary for the NEGF solver
# in the NEGF device.
dvc_negf.new_quantum_lead_bnd("source_bnd", essential=False)
dvc_negf.new_quantum_lead_bnd("drain_bnd", essential=False)

# Creating a submesh for the channel, source and drain.
submesh = SubMesh(mesh, channel_region_list)
submesh_source = SubMesh(mesh, ["source"])
submesh_drain = SubMesh(mesh, ["drain"])

# Restrict calculation of charge density on a subdevice
d_negf = SubDevice(dvc_negf, submesh)
d_source = SubDevice(dvc_negf, submesh_source)
d_drain = SubDevice(dvc_negf, submesh_drain)

# Interpolate the potential for the NEGF device and channel subdevices.
phi_interp = griddata(dvc.mesh.glob_nodes, dvc.phi, mesh.glob_nodes, method="nearest")
phi_interp_n = griddata(
    dvc.mesh.glob_nodes, dvc.phi, submesh.glob_nodes, method="nearest"
)
phi_interp_s = griddata(
    dvc.mesh.glob_nodes, dvc.phi, submesh_source.glob_nodes, method="nearest"
)
phi_interp_d = griddata(
    dvc.mesh.glob_nodes, dvc.phi, submesh_drain.glob_nodes, method="nearest"
)

# Set up the NEGF solver
Vds = 0.001  # drain-to-source voltage
NEGF_solver = NEGFSolver(
    d=dvc_negf,
    d_negf=d_negf,
    source_subd=d_source,
    drain_subd=d_drain,
    Vds=Vds,
    solver_params=solver_params,
)

# Set the potential from the adaptive mesh and update the Hamiltonian used in NEGF
NEGF_solver.d.set_potential(phi_interp)
NEGF_solver.d_negf.set_potential(phi_interp_n)
NEGF_solver.source_subd.set_potential(phi_interp_s)
NEGF_solver.drain_subd.set_potential(phi_interp_d)

# Update the solver and the Hamiltonian
NEGF_solver.update_hamiltonian()

# Compute the current
path_current = str(path_out / ("fdsoi_poisson_negf_master_current.png"))
current = NEGF_solver.get_current_adaptive(
    log_scale=False,
    show_figure=show_figure,
    path=path_current,
    title="Spectral current near equilibrium",
)
print("Current: " + str(current) + " [A]")

# compute LDOS
path_ldos = str(path_out / ("fdsoi_poisson_negf_master_ldos.png"))
NEGF_solver.plot_ldos(
    begin=begin,
    end=end,
    title="Local density of states",
    show_figure=show_figure,
    path=path_ldos,
    energies=np.linspace(-0.1 * ct.e, 0.03 * ct.e, 50),
)

# compute LDOS (zoomed-in)
path_ldos = str(path_out / ("fdsoi_poisson_negf_master_ldos_zoom.png"))
NEGF_solver.plot_ldos(
    begin=begin,
    end=end,
    title="Local density of states",
    show_figure=show_figure,
    path=path_ldos,
    energies=np.linspace(-0.02 * ct.e, 0.01 * ct.e, 50),
)

# Instantiate lever arm matrix solver
slv = LeverArmMatrixSolver(
    dvc,
    ["plunger_gate_bnd"],
    [plunger_gate],
    dot_region="dot",
    solver_params=lam_params,
)

# Calculate the lever arm matrix
bias_increment = 1e-3
lever_arm_matrix = slv.solve(bias_increment=bias_increment)

# Configure the many-body solver for the junction class.
params_mb = ManyBodySolverParams()
params_mb.n_degen = 2
params_mb.num_states = 2  # A larger value can be used for more accurate results
params_mb.alpha = lever_arm_matrix[0, 0]

# Initializing the Junction object used to calculate transport
jc = Junction(subdevice, many_body_solver_params=params_mb)

# Set a near-zero source-drain bias (same as the Vds used in NEGF solver)
jc.setVs(Vds)
jc.setVd(0)
# printing the Coulomb peak positions
print("Coulomb peak positions:")
print(jc.coulomb_peak_pos)

gate_bias_sweep = np.linspace(-0.1, 0.3, 400)
vec_Il = []
bartitle = "Calculating Coulomb peaks."
progress_bar = Bar(bartitle, max=len(gate_bias_sweep))
for idx, val in enumerate(gate_bias_sweep):
    jc.setVg(val)
    Il, Ir, prob = seq_tunnel_curr(jc)  # transport calculation
    vec_Il += [Il]
    progress_bar.next()
progress_bar.finish()

vec_Il = np.array(vec_Il)
np.savetxt(path_out / "PNEGF_current_scale.txt", vec_Il)

peaks, _ = find_peaks(vec_Il, 2e-20)
fig, ax = plt.subplots()
ax.plot(gate_bias_sweep, vec_Il, "-", label="Master equation")
ax.plot(gate_bias_sweep[peaks], vec_Il[peaks], "x", label="Coulomb peak maxima")
ax.set_ylabel("Current (arbitrary units)")
ax.set_xlabel("Gate bias (V)")
ax.legend()
fig.savefig(path_out / "PNEGF_current_featureless.png")
plt.show()

# scaling the values
scaling_value = current / vec_Il[peaks[0]]
vec_Il = scaling_value * vec_Il

fig, ax = plt.subplots()
ax.plot(gate_bias_sweep, vec_Il, "-")
ax.plot(gate_bias_sweep[peaks], vec_Il[peaks], "x")
ax.set_ylabel("Current (A)")
ax.set_xlabel("Gate bias (V)")
fig.savefig(path_out / "PNEGF_Current_scale.png")
plt.show()

# Calculating the charge stability diagram
source_drain_sweep = np.linspace(-0.1, 0.1, 150)
gate_bias_sweep = np.linspace(-0.1, 0.3, 150)
grid_Il = np.zeros((150, 150))
bartitle = "Calculating charge-stability diagram."
progress_bar = Bar(bartitle, max=len(gate_bias_sweep) * len(source_drain_sweep))
for idx, val_p in enumerate(gate_bias_sweep):
    jc.setVg(val_p)
    for idy, val in enumerate(source_drain_sweep):
        jc.setVs(val)
        Il, Ir, prob = seq_tunnel_curr(jc)  # transport calculation
        grid_Il[idx, idy] = Il
        progress_bar.next()
progress_bar.finish()

# Saving the data
np.savetxt(path_out / "PNEGF_CSD_scale.txt", grid_Il)

# Plotting the charge stability diagram
fig, ax = plt.subplots()
grid_Il = scaling_value * grid_Il
im = ax.imshow(
    np.transpose(grid_Il),
    interpolation="bilinear",
    cmap="bwr",
    extent=[-0.1, 0.3, -0.1, 0.1],
    aspect=2,
)
ax.set_ylabel(r"Drain--source bias $V_{DS}$ (V)")
ax.set_xlabel("Gate bias (V)")
cbar = fig.colorbar(im)
cbar.ax.set_ylabel("Current (A)")
fig.savefig(path_out / "PNEGF_Charge_diagram.png")
plt.show()
