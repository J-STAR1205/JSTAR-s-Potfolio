__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import io
from qtcad.device import analysis
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
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

# Define the dot region as a list of region labels
dot_region = ["channel", "channel_ox"]

# Set up the dot region in which no classical charge is allowed
d.set_dot_region(dot_region)

# Configure the non-linear Poisson solver
params_poisson = SolverParams()
params_poisson.tol = 1e-5  # Convergence threshold (tolerance) for the error
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 10000
params_poisson.refined_region = dot_region
params_poisson.h_refined = 0.25

# Create an adaptive-mesh non-linear Poisson solver
s = PoissonSolver(d, solver_params=params_poisson, geo_file=path_geo)

# Self-consistent solution
s.solve()

# Plot band diagram
analysis.plot_bands(d, (0, 0, 0), (0, 0, Ltot))

# Get the potential energy from the band edge for usage in the Schrodinger
# solver
d.set_V_from_phi()

# Create a submesh including only the dot region
submesh = SubMesh(d.mesh, dot_region)
submesh.show()

# Create a subdevice for the dot region
subdevice = SubDevice(d, submesh)

# Create a Schrodinger solver
schrod_solver = SchrodingerSolver(subdevice)

# Solve Schrodinger's equation
schrod_solver.solve()

# Print energies
subdevice.print_energies()

# Plot the first three eigenfunctions
analysis.plot_slices(submesh, subdevice.eigenfunctions[:, 0], title="Ground state")
analysis.plot_slices(
    submesh, subdevice.eigenfunctions[:, 1], title="First excited state"
)
analysis.plot_slices(
    submesh, subdevice.eigenfunctions[:, 2], title="Second excited state"
)

# Define the xy plane for the slice
normal = (0, 0, 1)
origin = (0, 0, 1e-8)
# Plot the first eigenfunction in z=10 nm plane.
analysis.plot_slice(
    submesh,
    subdevice.eigenfunctions[:, 0],
    normal=normal,
    origin=origin,
    title="Ground state at $z$=10 nm",
)
