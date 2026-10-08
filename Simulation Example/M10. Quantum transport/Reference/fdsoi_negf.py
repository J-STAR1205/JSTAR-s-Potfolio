__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
from os import cpu_count
import numpy as np
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.transport.negf_poisson import Solver, SolverParams

# Set up paths
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh / "fdsoi_negf.msh"))

# Define the device object
d = Device(mesh, conf_carriers="e")
d.set_temperature(0.1)

# Create the regions
ndoping = 2e20 * 1e6
d.new_region("top_ox", mt.SiO2_ideal)
d.new_region("source", mt.Si, ndoping=ndoping)
d.new_region("channel", mt.Si)
d.new_region("drain", mt.Si, ndoping=ndoping)
d.new_region("side_bot_ox", mt.SiO2_ideal)

# Set up boundary condition for back gate
back_gate_bias = -0.2
doping_back_gate = 1e15 * 1e6
binding_energy = 46 * 1e-3 * ct.e
d.new_frozen_bnd(
    "back_gate_bnd", back_gate_bias, mt.Si, doping_back_gate, "n", binding_energy
)

# Set boundary conditions for source and drain of FET for NEGF simulation
d.new_quantum_lead_bnd("source_bnd", essential=False)
d.new_quantum_lead_bnd("drain_bnd", essential=False)
Vds = 0.05  # drain-to-source voltage

# Restrict calculation of charge density on a subdevice
submesh = SubMesh(mesh, ["source", "channel", "drain"])
d_negf = SubDevice(d, submesh)

# NEGF-Poisson solver parameters
solver_params = SolverParams()
# Compute NEGF at different energy values in parallel
if cpu_count() >= 4:
    n_jobs = 4
else:
    n_jobs = 1
solver_params.n_jobs = n_jobs
mixing_param = [0.05, 0.05, 0.02]

# Define parameters pertaining to linecuts over which various quantities will
# be plotted
t_top_ox = 2 * scaling  # top oxide thickness
t_si = 3 * scaling  # silicon film thickness
begin = (0, np.amin(mesh.glob_nodes[:, 1]), -t_top_ox - t_si / 5)
end = (0, np.amax(mesh.glob_nodes[:, 1]), -t_top_ox - t_si / 5)
energies = np.linspace(-0.15 * ct.e, 0.05 * ct.e, 50)
show_figure = False  # whether plots are shown

# Define the voltages corresponding to off state, sequential tunneling regime,
# and on state, respectively
Ew = mt.Si.Eg / 2 + mt.Si.chi  # gate metal workfunction
barrier_gate_biases = np.array([0.7, 0.7, 1.7])
plunger_gate_biases = np.array([0.7, 0.76, 1.7])
current = np.zeros_like(barrier_gate_biases)  # drain-to-source current
regime_str = ["Off state", "Sequential tunneling", "On state"]

for ii, barrier_gate_bias in enumerate(barrier_gate_biases):
    # Set the boundary conditions for the top gates
    d.new_gate_bnd("barrier_gate_1_bnd", barrier_gate_bias, Ew)
    d.new_gate_bnd("plunger_gate_bnd", plunger_gate_biases[ii], Ew)
    d.new_gate_bnd("barrier_gate_2_bnd", barrier_gate_bias, Ew)

    # Set up NEGF-Poisson solver, and solve
    solver_params.mixing_param = mixing_param[ii]
    solver = Solver(d=d, d_negf=d_negf, Vds=Vds, solver_params=solver_params)
    solver.solve()

    # Plot and save various quantities of interest
    ii_str = "_" + str(ii)
    title = regime_str[ii]  # title of the plots
    # Conduction band edge
    path_band = str(path_out / ("band" + ii_str + ".png"))
    solver.plot_band_edge(
        begin=begin, end=end, title=title, show_figure=show_figure, path=path_band
    )
    # Local density of states
    path_ldos = str(path_out / ("ldos" + ii_str + ".png"))
    solver.plot_ldos(
        begin=begin,
        end=end,
        energies=energies,
        title=title,
        show_figure=show_figure,
        path=path_ldos,
    )
    # Spectral current (and current)
    path_spectral_current = str(path_out / ("spectral_current" + ii_str + ".png"))
    current[ii] = solver.get_current_adaptive(
        log_scale=False,
        title=title,
        show_figure=show_figure,
        path=path_spectral_current,
    )

# Print values of current in three investigated regimes
for ii, regime in enumerate(regime_str):
    print("-".join(regime.split()) + " regime current: " + str(current[ii]) + " A.")
