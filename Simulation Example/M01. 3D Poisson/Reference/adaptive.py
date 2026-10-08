__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import io
from qtcad.device import analysis
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.poisson import Solver, SolverParams
import pathlib

# Define some device parameters
Vgs = 2.0  # Gate-source bias
Ew = mt.Si.chi + mt.Si.Eg / 2  # Metal work function
Ltot = 20e-9  # Device height in m
radius = 2.5e-9  # Device radius in m

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "nanowire_adaptive.msh4"
path_geo = script_dir / "meshes" / "nanowire_adaptive.geo_unrolled"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, path_mesh)

# Create device from mesh and set temperature to 10 mK
d = Device(mesh, conf_carriers="e")
d.set_temperature(10e-3)

# Create device regions
d.new_region("channel_ox", mt.SiO2)
d.new_region("source_ox", mt.SiO2)
d.new_region("drain_ox", mt.SiO2)
d.new_region("channel", mt.Si)
d.new_region("source", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
d.new_region("drain", mt.Si, pdoping=5e20 * 1e6, ndoping=0)

# Create device boundaries
d.new_ohmic_bnd("source_bnd")
d.new_ohmic_bnd("drain_bnd")
d.new_gate_bnd("gate_bnd", Vgs, Ew)

# Configure the non-linear Poisson solver
params_poisson = SolverParams()
params_poisson.tol = 1e-5  # Convergence threshold (tolerance) for the error
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 0
params_poisson.refined_mesh_filename = (
    script_dir / "meshes" / "nanowire_adaptive_refined.msh4"
)

# Create an adaptive-mesh non-linear Poisson solver
s = Solver(d, solver_params=params_poisson, geo_file=path_geo)

# Self-consistent solution
s.solve()

# Plot band diagram
analysis.plot_bands(d, (0, 0, 0), (0, 0, Ltot))
