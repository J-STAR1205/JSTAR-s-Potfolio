__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import pathlib
from progress.bar import ChargingBar as Bar
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.schrodinger import Solver
from qtcad.device.schrodinger import SolverParams as SingleBodySolverParams
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from qtcad.transport.junction import Junction
from qtcad.transport.lead import WKBLead
from qtcad.transport.mastereq import seq_tunnel_curr

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes/" / "lead_dot_lead.msh"
scalingFactor = 1e-9
mesh = Mesh(scalingFactor, str(path_mesh))

# Create the device from the mesh
d = Device(mesh, conf_carriers="e")

d.new_region("dot", mt.GaAs)
d.new_region("lead", mt.GaAs)
d.new_region("lead2", mt.GaAs)
d.new_region("barrier", mt.GaAs)
d.new_region("barrier2", mt.GaAs)

# Global nodes
x = mesh.glob_nodes[:, 0]
y = mesh.glob_nodes[:, 1]
z = mesh.glob_nodes[:, 2]

# Potential along y and z
# Position of the harmonic oscillator minimum
y0 = 0
z0 = 0
Ky = 2e-7  # Parabolic in y
Kz = 1e-4  # Parabolic in z
Vy = Ky / 2.0 * (y - y0) ** 2
Vz = Kz / 2.0 * (z - z0) ** 2
# Set the potential energy
d.set_V(Vy + Vz)

# Add the potential energy along x
d.add_to_V(12e-3 * ct.e, region="barrier")
d.add_to_V(12e-3 * ct.e, region="barrier2")
d.add_to_V(-10e-3 * ct.e, region="lead")
d.add_to_V(-10e-3 * ct.e, region="lead2")
d.add_to_V(-1e-2 * ct.e)

# Create submesh
dotmesh = SubMesh(d.mesh, ["dot", "barrier", "barrier2"])
# Create subdevice
dot = SubDevice(d, dotmesh)
# Solve the Schrodinger equation.
# Configure the Schrödinger solver
params_schrod = SingleBodySolverParams()
params_schrod.num_states = 10  # Number of energy levels to consider
params_schrod.tol = 1e-9  # Set the tolerance for convergence

# Create solver
s = Solver(dot, solver_params=params_schrod)
s.solve()  # solve

# Take linecuts to define lead potential
# Length of dot region
L_dot = 200 * scalingFactor
# Length of left lead
L_lead1 = 50 * scalingFactor
# Length of right lead
L_lead2 = 50 * scalingFactor

x0 = L_lead1 - 10 * scalingFactor
x0R = L_lead1 + L_dot + 10 * scalingFactor

L_intersection = (x0, y0, z0)
R_intersection = (x0R, y0, z0)

# Define the leads.
mass = 1 / mt.GaAs.Me_inv[0, 0]
R_lead = WKBLead(
    d,
    ["dot", "barrier", "barrier2"],
    R_intersection,
    "x",
    ["lead2"],
    mass=mass,
    n1=1,
    n2=1,
    extra_region=["lead"],
)
L_lead = WKBLead(
    d,
    ["dot", "barrier", "barrier2"],
    L_intersection,
    "x",
    ["lead"],
    mass=mass,
    n1=1,
    n2=1,
    extra_region=["lead2"],
)

# Plot lead wavefunctions for a given energy
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)
Etest = 5e-3 * ct.e
R_lead.plot_lead_wf(Etest, 0, 0, "x", "y", "z", path=str(out_dir / "R_lead_wf.png"))
L_lead.plot_lead_wf(Etest, 0, 0, "x", "y", "z", path=str(out_dir / "L_lead_wf.png"))

# Define the junction (left lead - dot - right lead)
many_body_solver_params = ManyBodySolverParams()
many_body_solver_params.num_states = 1
many_body_solver_params.n_degen = 2
junc = Junction(dot, L_lead, R_lead, many_body_solver_params=many_body_solver_params)

junc.set_temperature(1)  # set temperature in unit of K
# Tunnelling computed using the WKB approximation
# Set source and drain voltages to low-bias regime
junc.setVs(2.5e-3)  # set source voltage in volts
junc.setVd(0.0)  # set drain voltage in volts

# Compute Coulomb peaks
v_gate_rng = np.linspace(-0.001, 0.008, num=100)
currentl = []

# Progress bar for the calculation of the Coulomb peaks
number_points = v_gate_rng.size  # number of sampled voltages
progress_bar_peaks = Bar("Computing Coulomb peaks", max=number_points)

for V in v_gate_rng:  # loop over gate voltage
    junc.setVg(V)
    Il, Ir, prob = seq_tunnel_curr(junc)  # transport calculation
    currentl.append(-Il)
    progress_bar_peaks.next()

progress_bar_peaks.finish()

I = np.array(currentl)
np.savetxt(str(out_dir / "coulomb_peaks.txt"), np.column_stack([v_gate_rng, I]))

# Plot
fig, ax = plt.subplots(1)
ax.set_xlabel("Gate voltage (V)")
ax.set_ylabel("Current (nA)")
ax.plot(v_gate_rng, -I / 1e-9)

fig.savefig(str(out_dir / "coulomb_peaks.png"))

# For visualization purposes only: uniform rates assumed for all transitions.
junc.WKB = False  # 'Turn off' WKB approximation

# Compute Coulomb diamonds.
v_gate_rng = np.linspace(-0.001, 0.008, num=150)
v_left_rng = np.linspace(-0.01, 0.01, num=150)
diam = np.zeros((v_gate_rng.size, v_left_rng.size), dtype=float)

# Progress bar for the calculation of the Coulomb diamonds
number_grid_points = v_gate_rng.size * v_left_rng.size  # number of sampled grid points
progress_bar_diamond = Bar("Computing charge stability diagram", max=number_grid_points)

for i in range(v_gate_rng.size):  # loop over gate voltage
    v_gate = v_gate_rng[i]
    junc.setVg(v_gate)
    for j in range(v_left_rng.size):  # loop over source/drain voltage
        v_left = v_left_rng[j]
        junc.setVs(v_left)
        # drain voltage is assumed to be opposite to source
        junc.setVd(-v_left)
        Il, Ir, prob = seq_tunnel_curr(junc)  # transport calculation
        diam[i, j] = Il
        progress_bar_diamond.next()

progress_bar_diamond.finish()

diam = np.abs(diam[:, 1:] - diam[:, :-1])  # differential conductance

# plot results
fig, axs = plt.subplots()
axs.set_ylabel("$V_s(V)$", fontsize=16)
axs.set_xlabel("$V_g(V)$", fontsize=16)
diff_cond_map = axs.imshow(
    np.transpose(diam),
    interpolation="bilinear",
    extent=[v_gate_rng[0], v_gate_rng[-1], v_left_rng[-2], v_left_rng[0]],
    aspect="auto",
)
fig.colorbar(diff_cond_map, ax=axs, label="Differential conductance (arb. units)")
fig.tight_layout()
fig.savefig(str(out_dir / "coulomb_diamonds.png"))
np.save(str(out_dir / "coulomb_diamonds.npy"), diam)
