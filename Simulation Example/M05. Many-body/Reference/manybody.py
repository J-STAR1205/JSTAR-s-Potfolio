__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import io
from qtcad.device import analysis as an
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.schrodinger import Solver as SingleBodySolver
import pathlib
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams

# Define some device parameters
Vgs = 0.5  # Gate-source bias
Ew = mt.Si.chi + mt.Si.Eg / 2  # Metal work function
Ltot = 20e-9  # Device height in m
radius = 2.5e-9  # Device radius in m

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes/" / "nanowire.msh"
path_hdf5 = script_dir / "output/" / "nanowire.hdf5"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# Create device from mesh and set statistics.
# It is important to set statistics before creating boundaries since
# they determine the boundary values, Failing to do so will produce
# inconsistent results.
dvc = Device(mesh, conf_carriers="e")
dvc.statistics = "FD_approx"  # Analytic approximation to Fermi-Dirac statistics

# Create regions first.
# Make sure each element is assigned to a region, otherwise default
# (silicon) parameters are used
dvc.new_region("channel_ox", mt.SiO2)
dvc.new_region("source_ox", mt.SiO2)
dvc.new_region("drain_ox", mt.SiO2)
dvc.new_region("channel", mt.Si)
dvc.new_region("source", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
dvc.new_region("drain", mt.Si, pdoping=5e20 * 1e6, ndoping=0)

# Then create boundaries
dvc.new_ohmic_bnd("source_bnd")
dvc.new_ohmic_bnd("drain_bnd")
dvc.new_gate_bnd("gate_bnd", Vgs, Ew)

# Load the electic potential
dvc.set_potential(io.load(path_hdf5, "phi"))

# Solve the single-electron problem
slv = SingleBodySolver(dvc)
slv.solve()
dvc.print_energies()

# Solve the many-body problem
solver_params = SolverParams()
solver_params.num_states = 3
solver_params.n_degen = 2
solver_params.alpha = 0.5
solver_params.num_particles = [0, 1, 2, 3]
slv = ManyBodySolver(dvc, solver_params=solver_params)
slv.solve()

# Many-body subspace properties
print("Many-body subspaces:", dvc.many_body_subspaces)
print("Number of subspaces:", len(dvc.many_body_subspaces))
print("Electron number in subspace 2:", dvc.many_body_subspaces[2].N)

print("Two-electron many-body basis set (integers): ")
print(dvc.get_N_particle_subspace(2).get_bas_set())

print("Two-electron many-body basis set (binary strings): ")
print(dvc.get_N_particle_subspace(2).get_bas_set(dtype="str"))

print("Two-electron many-body basis set (array): ")
print(dvc.get_N_particle_subspace(2).get_bas_set(dtype="array"))

print("Two-electron many-body eigenvalues (eV): ")
print(dvc.get_N_particle_subspace(2).eigval / ct.e)

print("Two-electron many-body ground state:")
print(dvc.get_N_particle_subspace(2).eigvec[:, 0])
print("Two-electron many-body first excited state:")
print(dvc.get_N_particle_subspace(2).eigvec[:, 1])

print("One-electron many-body basis set:")
print(dvc.get_N_particle_subspace(1).get_bas_set(dtype="array"))
print("One-electron many-body eigenstates:")
print(dvc.get_N_particle_subspace(1).eigvec)

# Many-body state objects
print("Two-electron ground-state coefficients and energy (eV)")
print(dvc.get_N_particle_subspace(2).get_eig_state(0).coeff)
print(dvc.get_N_particle_subspace(2).get_eig_state(0).energy / ct.e)

# Chemical potentials and Coulomb peak positions
print("Chemical potentials (eV)")
print(dvc.chem_potentials / ct.e)
print("Coulomb peak positions (V)")
print(dvc.coulomb_peak_pos)

# Plot electron density
an.plot_slices(
    mesh,
    dvc.get_many_body_density(1, 0),
    title="One-electron ground state particle density (1/m^3)",
)
an.plot_slices(
    mesh,
    dvc.get_many_body_density(1, 2),
    title="One-electron second-excited state particle density (1/m^3)",
)
an.plot_slices(
    mesh,
    dvc.get_many_body_density(2, 0),
    title="Two-electron ground state particle density (1/m^3)",
)
an.plot_slices(
    mesh,
    dvc.get_many_body_density(2, 2),
    title="Two-electron second-excited state particle density (1/m^3)",
)
an.plot_slices(
    mesh,
    dvc.get_many_body_density(3, 0),
    title="Three-electron ground state particle density (1/m^3)",
)
an.plot_slices(
    mesh,
    dvc.get_many_body_density(3, 5),
    title="Three-electron fifth-excited state particle density (1/m^3)",
)
