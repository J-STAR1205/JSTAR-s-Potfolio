__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import io
from qtcad.device import analysis
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
import pathlib

# Define some device parameters
Vgs = 0.5  # Gate-source bias
Ew = mt.Si.chi + mt.Si.Eg / 2  # Metal work function
Ltot = 20e-9  # Device height in m
radius = 2.5e-9  # Device radius in m

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "nanowire.msh"
path_hdf5 = script_dir / "output" / "nanowire.hdf5"
path_vtu = script_dir / "output" / "nanowire.vtu"
path_energies = script_dir / "output" / "nanowire_energies.txt"
path_linecuts_z = script_dir / "output" / "nanowire_linecuts_z.txt"
path_linecuts_x = script_dir / "output" / "nanowire_linecuts_x.txt"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# Plot the mesh
# mesh.show()  # disabled for headless/automated run (prerequisite script only)

# Create device from mesh and set statistics.
# It is important to set statistics before creating boundaries since
# they determine the boundary values. Failing to do so will produce
# inconsistent results.
d = Device(mesh)
d.statistics = "FD_approx"  # Analytic approximation to Fermi-Dirac statistics (default)

# Create regions first.
# Make sure each element is assigned to a region, otherwise default
# (silicon) parameters are used
d.new_region("channel_ox", mt.SiO2)
d.new_region("source_ox", mt.SiO2)
d.new_region("drain_ox", mt.SiO2)
d.new_region("channel", mt.Si)
d.new_region("source", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
d.new_region("drain", mt.Si, pdoping=5e20 * 1e6, ndoping=0)

# Then create boundaries
d.new_ohmic_bnd("source_bnd")
d.new_ohmic_bnd("drain_bnd")
d.new_gate_bnd("gate_bnd", Vgs, Ew)

# Visualize the device
# d.show()  # disabled for headless/automated run (prerequisite script only)

# Set specific camera position and target
# d.show(camera_position=[10, 10, 30], camera_target=[0, 0, 10], camera_offset=[0, 0, 0])

# Create and configure the non-linear Poisson solver
params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-5  # Convergence threshold (tolerance) for the error
s = PoissonSolver(d, solver_params=params_poisson)

# Self-consistent solution
s.solve()

# To save throughout the device: n, p, phi, EC, EV
arrays_dict = {
    "n": d.n / 1e6,
    "p": d.p / 1e6,
    "phi": d.phi,
    "EC": d.cond_band_edge() / ct.e,
    "EV": d.vlnce_band_edge() / ct.e,
}
io.save(path_hdf5, arrays_dict)

# A given variable may then be loaded with
phi = io.load(path_hdf5, "phi")

# These variables may also be saved in vtu format for external visualization
io.save(path_vtu, arrays_dict, mesh)

# Linecuts along z with x=y=0 for the same quantities
# Save n, p and phi throughout the device and along z axis
arrays = [d.n, d.p, d.phi]
array_names = ["n", "p", "phi"]
analysis.save_linecut(
    str(path_linecuts_z), mesh, arrays, (0, 0, 0), (0, 0, Ltot), array_names=array_names
)

# Linecuts along x with z=L/2 and y=0 for the same quantities
analysis.save_linecut(
    str(path_linecuts_x),
    mesh,
    arrays,
    (-radius, 0, Ltot / 2.0),
    (radius, 0, Ltot / 2.0),
    array_names=array_names,
)

# Plot linecut of the potential along the z axis
analysis.plot_linecut(
    mesh, d.phi, (0, 0, 0), (0, 0, Ltot), title="Electric potential (V)"
)

# Plot band diagram
analysis.plot_bands(d, (0, 0, 0), (0, 0, Ltot))

# Calculate the confinement potential from the electric potential
# found above, and set the confined carriers to be electrons
# The confinement potential is stored in d.V
d.conf_carriers = "e"
d.set_V_from_phi()

# Then instantiate the Schrodinger solver and its parameters
params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-5  # Tolerance on energies in eV
params_schrod.num_states = 5  # Number of states to consider
s = SchrodingerSolver(d, solver_params=params_schrod)

# Solve Schrodinger's equation
s.solve()

# Output and save energy levels
d.print_energies()
d.save_energies(str(path_energies))

# Get and print the quantum-dot geometric properties
dot_props = analysis.analyze_dot(d, verbose=True)

# Plot p and phi (with n and p converted to cm^-3)
# on orthogonal planar slices inside the device
analysis.plot_slices(mesh, d.p / 1e6, title="p (cm^-3)")
analysis.plot_slices(mesh, d.phi, title="phi (V)")

# Plot ground and first few excited state wavefunctions along the same slices
analysis.plot_slices(mesh, d.eigenfunctions[:, 0], title="Ground state")
analysis.plot_slices(mesh, d.eigenfunctions[:, 1], title="First excited state")
analysis.plot_slices(mesh, d.eigenfunctions[:, 2], title="Second excited state")
