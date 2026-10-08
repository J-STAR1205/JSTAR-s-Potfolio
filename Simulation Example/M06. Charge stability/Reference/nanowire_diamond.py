__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from matplotlib import pyplot as plt
import pathlib
from progress.bar import ChargingBar as Bar
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh3d import Mesh
from qtcad.device import io
from qtcad.device import Device
from qtcad.device.schrodinger import Solver
from qtcad.device.many_body import SolverParams
from qtcad.transport.mastereq import seq_tunnel_curr
from qtcad.transport.junction import Junction

# Load the mesh
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes/" / "nanowire.msh"
path_in = script_dir / "output/" / "nanowire.hdf5"

mesh = Mesh(1e-9, str(path_mesh))
d = Device(mesh=mesh, conf_carriers="e")
d.new_region("channel_ox", mt.SiO2)
d.new_region("source_ox", mt.SiO2)
d.new_region("drain_ox", mt.SiO2)
d.new_region("channel", mt.Si)
d.new_region("source", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
d.new_region("drain", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
d.V = io.load(path_in, var_name="EC") * ct.e

# Solve Schrodinger's equation for single electrons
s = Solver(d)
s.solve()

# Many-body solver parameters
solver_params = SolverParams()
solver_params.num_states = 2  # Number of basis states to keep
solver_params.n_degen = 2  # spin degenerate system

# Initialize junction
jc = Junction(d, many_body_solver_params=solver_params)

jc.set_temperature(1)  # set temperature in unit of K

# Get Coulomb diamonds
v_gate_rng = np.linspace(0.45, 0.85, num=100)
v_source_rng = np.linspace(-0.08, 0.08, num=100)
diam = np.zeros((v_gate_rng.size, v_source_rng.size), dtype=float)

number_grid_points = (
    v_gate_rng.size * v_source_rng.size
)  # number of sampled grid points
progress_bar = Bar(
    "Computing charge stability diagram", max=number_grid_points
)  # initialize progress bar

# loop over gate voltage
for i, v_gate in enumerate(v_gate_rng):
    # Set gate voltage
    jc.setVg(v_gate)
    # loop over source/drain voltage
    for j, v_source in enumerate(v_source_rng):
        # Set source and drain voltages
        jc.setVs(v_source)
        jc.setVd(-v_source)  # Drain voltage opposite in sign to source

        # Compute currents
        Il, Ir, prob = seq_tunnel_curr(jc)  # transport calculation
        diam[i, j] = Il

        progress_bar.next()

progress_bar.finish()

# differential conductance
diam_diff_cond = np.abs(diam[:, 1:] - diam[:, :-1])

# plot results
fig, axs = plt.subplots()
axs.set_ylabel("$V_s(V)$", fontsize=16)
axs.set_xlabel("$V_g(V)$", fontsize=16)
diff_cond_map = axs.imshow(
    np.transpose(diam_diff_cond),
    interpolation="bilinear",
    extent=[v_gate_rng[0], v_gate_rng[-1], v_source_rng[-2], v_source_rng[0]],
    aspect="auto",
)
fig.colorbar(diff_cond_map, ax=axs, label="Differential conductance (arb. units)")
fig.tight_layout()
plt.show()
